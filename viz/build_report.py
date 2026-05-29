"""Generate plots and an HTML report for a run.

Produces under runs/<run_id>/:
  plots/progress.png    — score vs iteration, best-so-far overlay
  plots/pareto.png      — score vs molecular weight (constraint-aware)
  plots/grid.png        — 2D structures of top-K candidates with scores
  plots/genealogy.png   — parent→child edges across iterations
  report.html           — embeds all of the above, plus tables
"""
from __future__ import annotations

import base64
import io
import json
import shutil
from pathlib import Path
from typing import Any

import click
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from jinja2 import Template

from rowan_tools import state, surrogate

ROOT = Path(__file__).resolve().parent.parent
TEMPLATE_PATH = ROOT / "templates" / "report.html.j2"
BRAND_ASSETS_DIR = ROOT / "docs" / "brand"


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def _setup_style() -> None:
    plt.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "#fafafa",
        "axes.edgecolor": "#333",
        "axes.grid": True,
        "grid.color": "#e0e0e0",
        "grid.linestyle": "-",
        "grid.linewidth": 0.5,
        "font.family": "sans-serif",
        "font.size": 11,
    })


def plot_progress(cfg: dict[str, Any], cands: list[dict[str, Any]],
                  out: Path) -> None:
    _setup_style()
    fig, ax = plt.subplots(figsize=(9, 5.5))

    direction = cfg["objective_direction"]
    metric = cfg["primary_metric"]

    valid = [c for c in cands if c.get("score") is not None]
    if not valid:
        ax.text(0.5, 0.5, "No scored candidates yet", ha="center", va="center",
                transform=ax.transAxes, color="#666", fontsize=14)
        ax.set_axis_off()
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        return

    iters = np.array([c["iter"] for c in valid])
    scores = np.array([c["score"] for c in valid])
    sat = np.array([c.get("satisfies_constraints", True) for c in valid])

    # Jitter on x so overlapping points are visible.
    jx = iters + np.random.RandomState(42).uniform(-0.18, 0.18, size=len(iters))
    ax.scatter(jx[sat], scores[sat], s=70, c="#2563eb", edgecolor="white",
               linewidth=1.2, alpha=0.85, label="passes constraints", zorder=3)
    if (~sat).any():
        ax.scatter(jx[~sat], scores[~sat], s=50, c="#9ca3af", edgecolor="white",
                   linewidth=1.0, alpha=0.6, marker="x", label="violates constraints", zorder=2)

    # Best-so-far line (constraint-aware).
    by_iter = sorted(set(iters.tolist()))
    best = []
    running = None
    for n in by_iter:
        sub = [c for c in valid if c["iter"] <= n and c.get("satisfies_constraints", True)]
        if not sub:
            best.append(np.nan)
            continue
        cur = (max if direction == "maximize" else min)(sub, key=lambda c: c["score"])["score"]
        running = cur
        best.append(cur)
    ax.plot(by_iter, best, color="#dc2626", linewidth=2.5, marker="o",
            markersize=6, label="best-so-far", zorder=4)

    ax.set_xlabel("iteration")
    ax.set_ylabel(f"{metric}  ({'higher is better' if direction == 'maximize' else 'lower is better'})")
    ax.set_title(f"Optimization progress — {cfg['objective']}", fontsize=12)
    ax.legend(loc="best", framealpha=0.95)
    ax.set_xticks(by_iter)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_pareto(cfg: dict[str, Any], cands: list[dict[str, Any]],
                out: Path, x_key: str = "mw") -> None:
    _setup_style()
    fig, ax = plt.subplots(figsize=(8, 5.5))

    valid = [c for c in cands if c.get("score") is not None
             and c.get("metrics", {}).get(x_key) is not None]
    if not valid:
        ax.text(0.5, 0.5, "Pareto plot needs scored candidates",
                ha="center", va="center", transform=ax.transAxes, color="#666")
        ax.set_axis_off()
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        return

    xs = np.array([c["metrics"][x_key] for c in valid])
    ys = np.array([c["score"] for c in valid])
    iters = np.array([c["iter"] for c in valid])
    sat = np.array([c.get("satisfies_constraints", True) for c in valid])

    sc = ax.scatter(xs, ys, c=iters, cmap="viridis", s=80,
                    edgecolor=np.where(sat, "white", "red"), linewidth=1.4,
                    alpha=0.9)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("iteration")

    # Pareto front for minimizing x_key and optimizing score.
    direction = cfg["objective_direction"]
    pts = list(zip(xs.tolist(), ys.tolist()))
    front = []
    for i, (x, y) in enumerate(pts):
        dominated = False
        for j, (x2, y2) in enumerate(pts):
            if i == j:
                continue
            score_dominates = y2 >= y if direction == "maximize" else y2 <= y
            strictly_better = x2 < x or y2 != y
            if x2 <= x and score_dominates and strictly_better:
                dominated = True
                break
        if not dominated:
            front.append((x, y))
    front.sort(key=lambda p: p[0])
    if len(front) >= 2:
        fx, fy = zip(*front)
        ax.plot(fx, fy, color="#dc2626", linestyle="--", linewidth=1.5,
                alpha=0.75, label="pareto front")
        ax.legend()

    ax.set_xlabel(x_key)
    ax.set_ylabel(cfg["primary_metric"])
    ax.set_title(f"{cfg['primary_metric']} vs {x_key}")
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_grid(cfg: dict[str, Any], cands: list[dict[str, Any]],
              out: Path, k: int = 9) -> None:
    """2D structures of top-K candidates."""
    from rdkit import Chem
    from rdkit.Chem import Draw

    valid = [c for c in cands if c.get("score") is not None
             and c.get("satisfies_constraints", True)]
    if not valid:
        # Empty placeholder.
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No scored candidates yet", ha="center", va="center",
                transform=ax.transAxes, color="#666", fontsize=14)
        ax.set_axis_off()
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        return

    direction = cfg["objective_direction"]
    valid.sort(key=lambda c: c["score"], reverse=(direction == "maximize"))
    top = valid[:k]
    mols = []
    legends = []
    for c in top:
        m = Chem.MolFromSmiles(c["smiles"])
        if m is None:
            continue
        mols.append(m)
        legends.append(f"iter {c['iter']}  {cfg['primary_metric']}={c['score']:.2f}")

    if not mols:
        fig, ax = plt.subplots(figsize=(8, 4))
        ax.text(0.5, 0.5, "No drawable molecules", ha="center", va="center",
                transform=ax.transAxes, color="#666", fontsize=14)
        ax.set_axis_off()
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        return

    img = Draw.MolsToGridImage(mols, molsPerRow=3, subImgSize=(320, 260),
                               legends=legends, useSVG=False)
    img.save(out)


def plot_genealogy(cfg: dict[str, Any], cands: list[dict[str, Any]],
                   out: Path) -> None:
    _setup_style()
    fig, ax = plt.subplots(figsize=(10, 6))

    if not cands:
        ax.text(0.5, 0.5, "No candidates yet", ha="center", va="center",
                transform=ax.transAxes, color="#666")
        ax.set_axis_off()
        fig.savefig(out, dpi=140, bbox_inches="tight")
        plt.close(fig)
        return

    # Lay out: x = iteration, y = sorted within iteration by score.
    records = [
        {"node_id": f"{c.get('iter', 0)}:{i}:{c.get('smiles', '')}", "candidate": c}
        for i, c in enumerate(cands)
    ]
    by_iter: dict[int, list[dict[str, Any]]] = {}
    for record in records:
        by_iter.setdefault(record["candidate"]["iter"], []).append(record)

    direction = cfg["objective_direction"]
    pos: dict[str, tuple[float, float]] = {}
    node_for_candidate: dict[int, str] = {}
    nodes_by_smiles: dict[str, list[tuple[int, str]]] = {}
    for n, group in by_iter.items():
        group.sort(key=lambda r: (r["candidate"].get("score") if r["candidate"].get("score") is not None
                                  else (-1e9 if direction == "maximize" else 1e9)),
                   reverse=(direction == "maximize"))
        for i, record in enumerate(group):
            c = record["candidate"]
            node_id = record["node_id"]
            pos[node_id] = (n, -i)
            node_for_candidate[id(c)] = node_id
            nodes_by_smiles.setdefault(c["smiles"], []).append((n, node_id))

    for nodes in nodes_by_smiles.values():
        nodes.sort(key=lambda item: item[0])

    # Draw edges parent→child.
    for c in cands:
        p = c.get("parent_smiles")
        child_id = node_for_candidate.get(id(c))
        parent_nodes = [
            (n, node_id) for n, node_id in nodes_by_smiles.get(p, [])
            if n <= c["iter"]
        ]
        if p and parent_nodes and child_id in pos:
            _, parent_id = parent_nodes[-1]
            x1, y1 = pos[parent_id]
            x2, y2 = pos[child_id]
            ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                        arrowprops=dict(arrowstyle="->", color="#94a3b8",
                                        alpha=0.6, lw=1.0))

    # Draw nodes.
    xs, ys, colors = [], [], []
    for c in cands:
        node_id = node_for_candidate.get(id(c))
        if node_id not in pos:
            continue
        x, y = pos[node_id]
        xs.append(x); ys.append(y)
        if c.get("score") is None:
            colors.append("#9ca3af")
        elif not c.get("satisfies_constraints", True):
            colors.append("#fbbf24")
        else:
            colors.append("#2563eb")
    ax.scatter(xs, ys, c=colors, s=140, edgecolor="white", linewidth=1.4, zorder=3)

    # Highlight best-so-far.
    best = (max if direction == "maximize" else min)(
        [c for c in cands if c.get("score") is not None and c.get("satisfies_constraints", True)],
        key=lambda c: c["score"], default=None,
    )
    best_id = node_for_candidate.get(id(best)) if best else None
    if best_id in pos:
        x, y = pos[best_id]
        ax.scatter([x], [y], s=300, facecolors="none", edgecolors="#dc2626",
                   linewidth=2.5, zorder=4)

    ax.set_xlabel("iteration")
    ax.set_yticks([])
    ax.set_title("Candidate genealogy (parent → child)")
    ax.set_xticks(sorted(by_iter.keys()))
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def _placeholder(out: Path, msg: str) -> None:
    _setup_style()
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.text(0.5, 0.5, msg, ha="center", va="center", transform=ax.transAxes,
            color="#666", fontsize=13)
    ax.set_axis_off()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_chemical_space(cfg: dict[str, Any], cands: list[dict[str, Any]],
                        out: Path) -> None:
    """2D PCA of Morgan fingerprints, colored by score — chemotype coverage."""
    valid = []
    fps = []
    for c in cands:
        if c.get("score") is None:
            continue
        fp = surrogate.featurize(c.get("smiles", ""))
        if fp is None:
            continue
        valid.append(c)
        fps.append(fp)

    if len(valid) < 3:
        _placeholder(out, "Chemical-space map needs >= 3 scored candidates")
        return

    X = np.vstack(fps)
    try:
        from sklearn.decomposition import PCA

        coords = PCA(n_components=2, random_state=0).fit_transform(X)
    except Exception:
        _placeholder(out, "Could not compute PCA embedding")
        return

    _setup_style()
    fig, ax = plt.subplots(figsize=(8, 5.5))
    scores = np.array([c["score"] for c in valid])
    sat = np.array([c.get("satisfies_constraints", True) for c in valid])
    sc = ax.scatter(coords[:, 0], coords[:, 1], c=scores, cmap="viridis", s=90,
                    edgecolor=np.where(sat, "white", "red"), linewidth=1.4, alpha=0.9)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(cfg["primary_metric"])

    direction = cfg["objective_direction"]
    best = (max if direction == "maximize" else min)(valid, key=lambda c: c["score"])
    bi = valid.index(best)
    ax.scatter([coords[bi, 0]], [coords[bi, 1]], s=320, facecolors="none",
               edgecolors="#dc2626", linewidth=2.5, label="best", zorder=4)

    ax.set_xlabel("PC 1")
    ax.set_ylabel("PC 2")
    ax.set_title("Chemical space (Morgan-fingerprint PCA)")
    ax.legend(loc="best", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


def plot_parity(cfg: dict[str, Any], cands: list[dict[str, Any]],
                out: Path, model: str = "rf") -> None:
    """Surrogate cross-validated predicted-vs-actual — do we trust the model?"""
    rows: dict[str, dict[str, Any]] = {}
    for c in cands:
        if c.get("score") is None or not c.get("satisfies_constraints", True):
            continue
        key = surrogate.canonical(c.get("smiles", "")) or c.get("smiles")
        if key:
            rows[key] = c

    X, y = [], []
    for c in rows.values():
        fp = surrogate.featurize(c["smiles"])
        if fp is None:
            continue
        X.append(fp)
        y.append(float(c["score"]))

    if len(X) < surrogate.MIN_TRAIN:
        _placeholder(out, f"Parity plot needs >= {surrogate.MIN_TRAIN} "
                          f"constraint-passing candidates ({len(X)} so far)")
        return

    cv = surrogate.cv_predictions(np.vstack(X), np.asarray(y, dtype=float), model=model)
    if cv is None:
        _placeholder(out, "Not enough data for cross-validation")
        return

    y_true, y_pred = cv
    mae = float(np.mean(np.abs(y_true - y_pred)))
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - y_true.mean()) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 1e-12 else float("nan")

    _setup_style()
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(y_true, y_pred, s=80, c="#2563eb", edgecolor="white",
               linewidth=1.2, alpha=0.85)
    lo = float(min(y_true.min(), y_pred.min()))
    hi = float(max(y_true.max(), y_pred.max()))
    pad = 0.05 * (hi - lo + 1e-9)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color="#dc2626",
            linestyle="--", linewidth=1.5, alpha=0.8, label="ideal")
    ax.set_xlabel(f"measured {cfg['primary_metric']}")
    ax.set_ylabel(f"surrogate-predicted {cfg['primary_metric']}")
    ax.set_title(f"Surrogate parity ({model.upper()}, leave-out CV)  "
                 f"R^2={r2:.2f}  MAE={mae:.2f}  n={len(y_true)}")
    ax.legend(loc="best", framealpha=0.95)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------------------
# 3D coords for the best candidate (embedded as SDF for 3Dmol.js)
# ---------------------------------------------------------------------------

def best_sdf(best: dict[str, Any] | None) -> str:
    if not best:
        return ""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(best["smiles"])
    if m is None:
        return ""
    m = Chem.AddHs(m)
    try:
        AllChem.EmbedMolecule(m, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        return ""
    return Chem.MolToMolBlock(m)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_report(run_id: str, top_k: int = 9, live_refresh_seconds: int = 10) -> Path:
    """Build plots and report.html for a run, returning the HTML path."""
    cfg = state.load_config(run_id)
    cands = state.all_candidates(run_id)
    rd = state.run_dir(run_id)
    plots_dir = rd / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)
    brand_dir = rd / "brand"
    brand_dir.mkdir(parents=True, exist_ok=True)
    for logo_name in ("rowan-logo.svg", "rowan-mark.svg", "kdense-logo.png"):
        src = BRAND_ASSETS_DIR / logo_name
        if src.exists():
            shutil.copyfile(src, brand_dir / logo_name)

    plot_progress(cfg, cands, plots_dir / "progress.png")
    plot_pareto(cfg, cands, plots_dir / "pareto.png")
    plot_grid(cfg, cands, plots_dir / "grid.png", k=top_k)
    plot_genealogy(cfg, cands, plots_dir / "genealogy.png")
    plot_chemical_space(cfg, cands, plots_dir / "chemical_space.png")
    plot_parity(cfg, cands, plots_dir / "parity.png")

    best = state.best_so_far(run_id)
    sdf = best_sdf(best)

    iters = []
    for n in state.list_iterations(run_id):
        iters.append(state.load_iter(run_id, n))

    tmpl = Template(TEMPLATE_PATH.read_text())
    html = tmpl.render(
        cfg=cfg,
        iterations=iters,
        candidates=cands,
        best=best,
        best_sdf=sdf,
        n_iter=len(iters),
        n_candidates=len(cands),
        n_valid=sum(1 for c in cands if c.get("satisfies_constraints", True) and c.get("score") is not None),
        live_refresh_seconds=live_refresh_seconds,
    )
    report_path = rd / "report.html"
    report_path.write_text(html)
    return report_path


@click.command()
@click.option("--run", "run_id", required=True)
@click.option("--top-k", default=9, type=int, help="How many top molecules to show in the grid.")
def cli(run_id, top_k):
    click.echo(f"[{run_id}] {len(state.all_candidates(run_id))} candidates across {len(state.list_iterations(run_id))} iterations")
    report_path = build_report(run_id, top_k=top_k)
    click.echo(f"Wrote {report_path}")
    click.echo(f"Open it: file://{report_path.resolve()}")


if __name__ == "__main__":
    cli()
