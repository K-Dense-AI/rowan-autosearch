"""Surrogate-assisted active-learning advisor.

The agent is still the optimizer; this tool is a quantitative co-pilot. It
learns a cheap surrogate from a run's own scored history
(``runs/<run_id>/iterations/*.json``) and uses it to:

  * **dedup** candidate SMILES against molecules already evaluated (so Rowan
    credits are never re-spent),
  * **rank** the agent's proposals by Expected Improvement (direction-aware),
  * **annotate** novelty (Tanimoto distance to the evaluated set), Murcko
    scaffold, predicted value + uncertainty, and RDKit constraint pass/fail,
  * **select** a chemically diverse top-k via MaxMin so the search does not
    collapse into a single chemotype.

Everything here is pure-local: no Rowan API calls. The command emits
paste-ready ``SMILES|PARENT|NOTE`` lines that pipe straight into
``rowan-score``, or structured JSON with ``--json``.
"""
from __future__ import annotations

import json
import math
from typing import Any

import click
import numpy as np

from . import state

# Below this many usable training points the surrogate is untrustworthy, so we
# fall back to pure novelty ranking and say so plainly.
MIN_TRAIN = 6
FP_RADIUS = 2
FP_BITS = 2048


# ---------------------------------------------------------------------------
# Featurization (Morgan / ECFP4) + canonicalization
# ---------------------------------------------------------------------------

_FP_GEN = None


def _fp_generator():
    global _FP_GEN
    if _FP_GEN is None:
        from rdkit.Chem import rdFingerprintGenerator

        _FP_GEN = rdFingerprintGenerator.GetMorganGenerator(
            radius=FP_RADIUS, fpSize=FP_BITS
        )
    return _FP_GEN


def canonical(smiles: str) -> str | None:
    """Canonical SMILES for exact dedup keys, or None if unparseable."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def morgan_fp(smiles: str):
    """RDKit ExplicitBitVect (for Tanimoto), or None for invalid SMILES."""
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return _fp_generator().GetFingerprint(mol)


def fp_to_numpy(fp) -> np.ndarray:
    from rdkit import DataStructs

    arr = np.zeros((FP_BITS,), dtype=np.float64)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr


def featurize(smiles: str) -> np.ndarray | None:
    """ECFP4 bit vector as a float numpy array, or None for invalid SMILES."""
    fp = morgan_fp(smiles)
    if fp is None:
        return None
    return fp_to_numpy(fp)


def murcko_scaffold(smiles: str) -> str:
    """Bemis-Murcko scaffold SMILES (empty string if unavailable)."""
    try:
        from rdkit.Chem.Scaffolds import MurckoScaffold

        return MurckoScaffold.MurckoScaffoldSmiles(smiles) or ""
    except Exception:
        return ""


# ---------------------------------------------------------------------------
# Run history -> training set
# ---------------------------------------------------------------------------

def evaluated_smiles(run_id: str) -> dict[str, dict[str, Any]]:
    """Canonical SMILES -> latest candidate record that actually hit Rowan.

    Used for dedup: any candidate with a numeric ``score`` consumed credits,
    regardless of whether it passed constraints.
    """
    out: dict[str, dict[str, Any]] = {}
    for c in state.all_candidates(run_id):
        if c.get("score") is None:
            continue
        key = canonical(c.get("smiles", "")) or c.get("smiles")
        if key:
            out[key] = c  # later iterations overwrite earlier duplicates
    return out


def training_set(run_id: str):
    """Return (X, y, fps, canon_list) for scored, constraint-passing molecules.

    De-duplicated by canonical SMILES (latest scoring wins). ``X`` is an
    (n, FP_BITS) array, ``y`` an (n,) array of scores, ``fps`` the matching
    RDKit bit vectors (for Tanimoto), and ``canon_list`` their canonical SMILES.
    """
    rows: dict[str, dict[str, Any]] = {}
    for c in state.all_candidates(run_id):
        if c.get("score") is None or not c.get("satisfies_constraints", True):
            continue
        key = canonical(c.get("smiles", "")) or c.get("smiles")
        if key:
            rows[key] = c

    X, y, fps, canon = [], [], [], []
    for key, c in rows.items():
        fp = morgan_fp(c["smiles"])
        if fp is None:
            continue
        fps.append(fp)
        X.append(fp_to_numpy(fp))
        y.append(float(c["score"]))
        canon.append(key)

    if not X:
        return np.empty((0, FP_BITS)), np.empty((0,)), [], []
    return np.vstack(X), np.asarray(y, dtype=np.float64), fps, canon


# ---------------------------------------------------------------------------
# Acquisition: Expected Improvement (no scipy dependency)
# ---------------------------------------------------------------------------

def _norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def _norm_pdf(x: float) -> float:
    return math.exp(-0.5 * x * x) / math.sqrt(2.0 * math.pi)


def expected_improvement(mu: float, sigma: float, best: float,
                         direction: str, xi: float = 0.01) -> float:
    """Standard Expected Improvement, sign-flipped for minimization.

    With ``sigma`` ~ 0 EI collapses to the (clipped) raw improvement, so a
    confident prediction better than the incumbent still scores positively.
    """
    if direction == "maximize":
        improvement = mu - best - xi
    else:
        improvement = best - mu - xi

    if sigma <= 1e-9:
        return max(improvement, 0.0)

    z = improvement / sigma
    return improvement * _norm_cdf(z) + sigma * _norm_pdf(z)


# ---------------------------------------------------------------------------
# Surrogate models
# ---------------------------------------------------------------------------

def tanimoto_kernel(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Pairwise Tanimoto similarity between two stacks of binary fingerprints."""
    inter = A @ B.T
    a = A.sum(axis=1)
    b = B.sum(axis=1)
    denom = a[:, None] + b[None, :] - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        sim = np.where(denom > 0, inter / denom, 0.0)
    return sim


class _RFPredictor:
    """Random forest; uncertainty from the spread across trees."""

    def __init__(self, model) -> None:
        self.model = model

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        per_tree = np.stack([est.predict(X) for est in self.model.estimators_], axis=0)
        return per_tree.mean(axis=0), per_tree.std(axis=0)


class _GPPredictor:
    """Gaussian process with a Tanimoto kernel, implemented in numpy.

    Gives better-calibrated uncertainty than the forest on fingerprints, at
    the cost of an O(n^3) Cholesky (fine for the dozens-to-hundreds of points
    a run produces).
    """

    def __init__(self, X: np.ndarray, y: np.ndarray,
                 noise_frac: float = 1e-3) -> None:
        self.X = X
        self.y_mean = float(y.mean())
        self.amplitude = float(y.var()) or 1.0
        residual = y - self.y_mean

        K = self.amplitude * tanimoto_kernel(X, X)
        noise = max(self.amplitude * noise_frac, 1e-8)
        jitter = noise
        for _ in range(6):
            try:
                self.L = np.linalg.cholesky(K + jitter * np.eye(len(X)))
                break
            except np.linalg.LinAlgError:
                jitter *= 10
        else:
            self.L = np.linalg.cholesky(K + 1e-3 * np.eye(len(X)))
        self.alpha = np.linalg.solve(self.L.T, np.linalg.solve(self.L, residual))

    def predict(self, X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        Ks = self.amplitude * tanimoto_kernel(X, self.X)
        mu = self.y_mean + Ks @ self.alpha
        v = np.linalg.solve(self.L, Ks.T)
        var = self.amplitude - np.sum(v ** 2, axis=0)
        sigma = np.sqrt(np.clip(var, 1e-9, None))
        return mu, sigma


def fit_surrogate(X: np.ndarray, y: np.ndarray, model: str = "rf"):
    """Fit a surrogate predictor exposing ``predict(X) -> (mu, sigma)``."""
    if model == "gp":
        return _GPPredictor(X, y)

    from sklearn.ensemble import RandomForestRegressor

    n_estimators = max(100, min(400, 20 * len(X)))
    rf = RandomForestRegressor(
        n_estimators=n_estimators,
        max_features="sqrt",
        bootstrap=True,
        random_state=42,
        n_jobs=-1,
    )
    rf.fit(X, y)
    return _RFPredictor(rf)


def cv_predictions(X: np.ndarray, y: np.ndarray, model: str = "rf",
                   max_folds: int = 5) -> tuple[np.ndarray, np.ndarray] | None:
    """Out-of-fold (y_true, y_pred) via k-fold (LOO for tiny n) CV."""
    n = len(X)
    if n < 3:
        return None

    k = min(max_folds, n)
    rng = np.random.RandomState(0)
    order = rng.permutation(n)
    folds = np.array_split(order, k)

    preds = np.full(n, np.nan)
    for fold in folds:
        test_idx = np.asarray(fold)
        train_idx = np.setdiff1d(order, test_idx, assume_unique=False)
        if len(train_idx) < 2 or len(test_idx) == 0:
            continue
        try:
            predictor = fit_surrogate(X[train_idx], y[train_idx], model=model)
            mu, _ = predictor.predict(X[test_idx])
        except Exception:
            continue
        preds[test_idx] = mu

    mask = ~np.isnan(preds)
    if mask.sum() < 2:
        return None
    return y[mask], preds[mask]


def cv_metrics(X: np.ndarray, y: np.ndarray, model: str = "rf",
               max_folds: int = 5) -> dict[str, float] | None:
    """k-fold (LOO for tiny n) CV R^2 and MAE so the agent can gauge trust."""
    cv = cv_predictions(X, y, model=model, max_folds=max_folds)
    if cv is None:
        return None
    yt, yp = cv
    mae = float(np.mean(np.abs(yt - yp)))
    ss_res = float(np.sum((yt - yp) ** 2))
    ss_tot = float(np.sum((yt - yt.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else float("nan")
    return {"r2": r2, "mae": mae, "n": int(len(yt))}


# ---------------------------------------------------------------------------
# Novelty + diversity
# ---------------------------------------------------------------------------

def novelty(fp, train_fps: list) -> float:
    """1 - max Tanimoto similarity to the evaluated set (1.0 if set empty)."""
    if not train_fps:
        return 1.0
    from rdkit import DataStructs

    sims = DataStructs.BulkTanimotoSimilarity(fp, train_fps)
    return float(1.0 - max(sims)) if sims else 1.0


def maxmin_select(fps: list, k: int) -> list[int]:
    """Greedy MaxMin diverse subset, seeded by the first (best-ranked) item.

    ``fps`` must already be in priority order; returns selected indices in
    selection order.
    """
    from rdkit import DataStructs

    n = len(fps)
    if n == 0:
        return []
    if k >= n:
        return list(range(n))

    selected = [0]
    min_dist = [1.0 - s for s in DataStructs.BulkTanimotoSimilarity(fps[0], fps)]
    min_dist[0] = -1.0  # never reselect

    while len(selected) < k:
        nxt = int(np.argmax(min_dist))
        if min_dist[nxt] < 0:
            break
        selected.append(nxt)
        sims = DataStructs.BulkTanimotoSimilarity(fps[nxt], fps)
        for i in range(n):
            d = 1.0 - sims[i]
            if d < min_dist[i]:
                min_dist[i] = d
        min_dist[nxt] = -1.0
    return selected


# ---------------------------------------------------------------------------
# Candidate parsing
# ---------------------------------------------------------------------------

def _parse_candidate(raw: str) -> dict[str, str]:
    parts = raw.split("|", 2)
    return {
        "smiles": parts[0].strip(),
        "parent_smiles": parts[1].strip() if len(parts) > 1 else "",
        "design_note": parts[2].strip() if len(parts) > 2 else "",
    }


def _candidates_from_parent(parent: str, strategy: str, n: int) -> list[dict[str, str]]:
    from . import propose

    fns = {
        "bioisostere": propose.bioisostere_variants,
        "scan-subst": propose.aromatic_substitution_variants,
        "brics": propose.brics_variants,
    }
    keys = list(fns) if strategy == "all" else [strategy]
    out: list[dict[str, str]] = []
    seen: set[str] = set()
    for key in keys:
        for smiles, note in fns[key](parent, n=n):
            if smiles in seen:
                continue
            seen.add(smiles)
            out.append({"smiles": smiles, "parent_smiles": parent, "design_note": note})
    return out


# ---------------------------------------------------------------------------
# Core ranking
# ---------------------------------------------------------------------------

def rank_candidates(
    run_id: str,
    candidates: list[dict[str, str]],
    *,
    cfg: dict[str, Any],
    xi: float = 0.01,
    model: str = "rf",
    respect_constraints: bool = True,
    include_seen: bool = False,
) -> dict[str, Any]:
    """Score candidates with the surrogate; return a structured result dict."""
    from . import score as score_mod

    direction = cfg.get("objective_direction", "maximize")
    constraints = cfg.get("constraints", {}) if respect_constraints else {}

    X_train, y_train, train_fps, _ = training_set(run_id)
    seen = evaluated_smiles(run_id)
    best = state.best_so_far(run_id)
    incumbent = float(best["score"]) if best else None

    modeled = len(X_train) >= MIN_TRAIN
    predictor = fit_surrogate(X_train, y_train, model=model) if modeled else None
    cv = cv_metrics(X_train, y_train, model=model) if modeled else None

    records: list[dict[str, Any]] = []
    for cand in candidates:
        smiles = cand["smiles"]
        rec: dict[str, Any] = {
            "smiles": smiles,
            "parent_smiles": cand.get("parent_smiles", ""),
            "design_note": cand.get("design_note", ""),
        }
        key = canonical(smiles)
        if key is None:
            rec.update(status="invalid", error="RDKit failed to parse SMILES")
            records.append(rec)
            continue

        rec["canonical_smiles"] = key
        prior = seen.get(key)
        rec["already_evaluated"] = prior is not None
        if prior is not None:
            rec["prior_score"] = prior.get("score")
            rec["prior_iter"] = prior.get("iter")

        desc = score_mod.local_descriptors(smiles)
        ok, failed = score_mod.check_constraints(desc, constraints)
        rec["satisfies_constraints"] = ok
        rec["constraint_failures"] = failed
        rec["scaffold"] = murcko_scaffold(smiles)

        fp = morgan_fp(smiles)
        rec["novelty"] = novelty(fp, train_fps)
        rec["_fp"] = fp

        if predictor is not None:
            mu, sigma = predictor.predict(fp_to_numpy(fp)[None, :])
            mu_f, sigma_f = float(mu[0]), float(sigma[0])
            rec["pred"] = mu_f
            rec["pred_sigma"] = sigma_f
            if incumbent is not None:
                rec["ei"] = expected_improvement(mu_f, sigma_f, incumbent, direction, xi)
            else:
                rec["ei"] = sigma_f  # nothing to improve on yet: chase information
        records.append(rec)

    # Decide eligibility for recommendation.
    def eligible(rec: dict[str, Any]) -> bool:
        if rec.get("status") == "invalid":
            return False
        if respect_constraints and not rec.get("satisfies_constraints", True):
            return False
        if not include_seen and rec.get("already_evaluated"):
            return False
        return True

    pool = [r for r in records if eligible(r)]
    rank_key = "ei" if modeled else "novelty"
    pool.sort(key=lambda r: r.get(rank_key, float("-inf")), reverse=True)

    return {
        "run_id": run_id,
        "direction": direction,
        "incumbent": incumbent,
        "model": model,
        "modeled": modeled,
        "n_train": int(len(X_train)),
        "min_train": MIN_TRAIN,
        "rank_key": rank_key,
        "cv": cv,
        "records": records,
        "ranked_pool": pool,
    }


def select_recommendations(result: dict[str, Any], top_k: int,
                           diverse: bool) -> list[dict[str, Any]]:
    pool = result["ranked_pool"]
    if not pool:
        return []
    if diverse and len(pool) > top_k:
        cap = max(top_k * 3, top_k)
        head = pool[:cap]
        idx = maxmin_select([r["_fp"] for r in head], top_k)
        return [head[i] for i in idx]
    return pool[:top_k]


def _note(rec: dict[str, Any], modeled: bool) -> str:
    base = rec.get("design_note") or ""
    bits = []
    if modeled and rec.get("ei") is not None:
        bits.append(f"EI={rec['ei']:.3f}")
        bits.append(f"pred={rec['pred']:.2f}+/-{rec['pred_sigma']:.2f}")
    bits.append(f"novelty={rec['novelty']:.2f}")
    if rec.get("scaffold"):
        bits.append(f"scaffold={rec['scaffold']}")
    annot = " ".join(bits)
    return f"{base} [{annot}]" if base else f"[{annot}]"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

@click.command()
@click.option("--run", "run_id", required=True)
@click.option("--candidate", "candidates", multiple=True,
              help="Candidate as 'SMILES|PARENT|NOTE', repeatable.")
@click.option("--from-parent", default=None,
              help="Auto-generate candidates from this parent via rowan-propose.")
@click.option("--strategy",
              type=click.Choice(["bioisostere", "scan-subst", "brics", "all"]),
              default="all", help="Strategy for --from-parent generation.")
@click.option("--n", default=20, type=int, help="Max variants per strategy for --from-parent.")
@click.option("--top-k", default=4, type=int, help="How many candidates to recommend.")
@click.option("--xi", default=0.01, type=float, help="EI exploration parameter.")
@click.option("--model", type=click.Choice(["rf", "gp"]), default="rf",
              help="Surrogate model: random forest (default) or Tanimoto-kernel GP.")
@click.option("--diverse/--no-diverse", default=True,
              help="Apply MaxMin diversity selection over the top pool.")
@click.option("--respect-constraints/--no-respect-constraints", default=True,
              help="Filter out candidates that fail the run's RDKit constraints.")
@click.option("--include-seen", is_flag=True,
              help="Include candidates already evaluated in this run.")
@click.option("--json", "as_json", is_flag=True, help="Emit structured JSON.")
def cli(run_id, candidates, from_parent, strategy, n, top_k, xi, model,
        diverse, respect_constraints, include_seen, as_json):
    """Pre-screen and rank candidates before spending Rowan credits.

    Recommends a diverse, non-duplicate, constraint-passing subset ranked by
    Expected Improvement (or novelty when there is too little data to model).
    Output is paste-ready ``SMILES|PARENT|NOTE`` for ``rowan-score``.
    """
    cfg = state.load_config(run_id)

    parsed = [_parse_candidate(c) for c in candidates if c.strip()]
    if from_parent:
        parsed.extend(_candidates_from_parent(from_parent, strategy, n))
    if not parsed:
        raise click.ClickException(
            "No candidates. Pass --candidate ... and/or --from-parent <smiles>."
        )

    result = rank_candidates(
        run_id, parsed, cfg=cfg, xi=xi, model=model,
        respect_constraints=respect_constraints, include_seen=include_seen,
    )
    recs = select_recommendations(result, top_k, diverse)

    if as_json:
        payload = {
            k: v for k, v in result.items() if k not in ("records", "ranked_pool")
        }
        strip = lambda r: {k: v for k, v in r.items() if k != "_fp"}
        payload["records"] = [strip(r) for r in result["records"]]
        payload["recommended"] = [strip(r) for r in recs]
        click.echo(json.dumps(payload, indent=2))
        return

    n_seen = sum(1 for r in result["records"] if r.get("already_evaluated"))
    n_invalid = sum(1 for r in result["records"] if r.get("status") == "invalid")
    n_violate = sum(1 for r in result["records"]
                    if r.get("status") != "invalid" and not r.get("satisfies_constraints", True))

    if result["modeled"]:
        cv = result["cv"]
        cv_str = (f"CV R^2={cv['r2']:.2f} MAE={cv['mae']:.2f} (n={cv['n']})"
                  if cv else "CV n/a")
        click.echo(
            f"[{run_id}] surrogate={result['model']} on n={result['n_train']} | "
            f"ranked by Expected Improvement | {cv_str}"
        )
        inc = result["incumbent"]
        click.echo(f"  incumbent {cfg['primary_metric']}="
                   f"{inc:.3f}" if inc is not None else "  no incumbent yet")
    else:
        click.echo(
            f"[{run_id}] cold start: only {result['n_train']}/{result['min_train']} "
            f"training points -> ranking by NOVELTY (no surrogate yet)."
        )
    click.echo(
        f"  inputs={len(result['records'])} | already-evaluated={n_seen} | "
        f"invalid={n_invalid} | constraint-violations={n_violate}"
    )

    if not recs:
        click.echo("\nNo eligible candidates to recommend "
                   "(all duplicates, invalid, or constraint-failing).")
        return

    click.echo(f"\nRecommended {len(recs)} candidate(s) to score next:\n")
    for rec in recs:
        click.echo(f"{rec['smiles']}|{rec.get('parent_smiles', '')}|{_note(rec, result['modeled'])}")

    click.echo("\nPaste the lines above into rowan-score (one --candidate each).")


if __name__ == "__main__":
    cli()
