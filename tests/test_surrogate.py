from __future__ import annotations

import json

import numpy as np
import pytest
from click.testing import CliRunner

from rowan_tools import state, surrogate


# ---------------------------------------------------------------------------
# Featurization + canonicalization
# ---------------------------------------------------------------------------

def test_featurize_is_deterministic_and_fixed_length():
    fp1 = surrogate.featurize("CCO")
    fp2 = surrogate.featurize("CCO")

    assert fp1 is not None
    assert fp1.shape == (surrogate.FP_BITS,)
    assert np.array_equal(fp1, fp2)


def test_featurize_returns_none_for_invalid_smiles():
    assert surrogate.featurize("not a smiles") is None


def test_canonical_collapses_equivalent_smiles_and_rejects_garbage():
    assert surrogate.canonical("OCC") == surrogate.canonical("CCO")
    assert surrogate.canonical("c1ccccc1") == surrogate.canonical("C1=CC=CC=C1")
    assert surrogate.canonical("not a smiles") is None


def test_morgan_fp_returns_none_for_invalid_smiles():
    assert surrogate.morgan_fp("not a smiles") is None
    assert surrogate.morgan_fp("CCO") is not None


def test_murcko_scaffold_strips_substituents_and_handles_acyclic():
    # Toluene and benzene share the benzene scaffold.
    assert surrogate.murcko_scaffold("Cc1ccccc1") == surrogate.murcko_scaffold("c1ccccc1")
    # Acyclic molecules have an empty Murcko scaffold.
    assert surrogate.murcko_scaffold("CCO") == ""


# ---------------------------------------------------------------------------
# Tanimoto kernel
# ---------------------------------------------------------------------------

def test_tanimoto_kernel_is_symmetric_with_unit_self_similarity():
    X = np.vstack([surrogate.featurize(s) for s in ["CCO", "CCN", "c1ccccc1"]])

    K = surrogate.tanimoto_kernel(X, X)

    assert K.shape == (3, 3)
    assert np.allclose(np.diag(K), 1.0)
    assert np.allclose(K, K.T)
    assert np.all((K >= 0) & (K <= 1.0 + 1e-9))


def test_tanimoto_kernel_handles_all_zero_rows_without_dividing_by_zero():
    zero = np.zeros((1, surrogate.FP_BITS))
    K = surrogate.tanimoto_kernel(zero, zero)

    assert K.shape == (1, 1)
    assert K[0, 0] == 0.0  # 0/0 is defined as 0, not NaN


# ---------------------------------------------------------------------------
# Expected Improvement
# ---------------------------------------------------------------------------

def test_expected_improvement_increases_with_mu_for_maximize():
    lo = surrogate.expected_improvement(1.0, 0.5, 0.0, "maximize")
    hi = surrogate.expected_improvement(2.0, 0.5, 0.0, "maximize")

    assert hi > lo > 0


def test_expected_improvement_mirrors_for_minimize():
    # For minimize, a *lower* predicted mu should improve more.
    worse = surrogate.expected_improvement(1.0, 0.5, 2.0, "minimize")
    better = surrogate.expected_improvement(0.0, 0.5, 2.0, "minimize")

    assert better > worse > 0


def test_expected_improvement_collapses_to_clipped_improvement_at_zero_sigma():
    # Confident and better than incumbent -> positive; confident and worse -> 0.
    assert surrogate.expected_improvement(2.0, 0.0, 1.0, "maximize") == pytest.approx(0.99)
    assert surrogate.expected_improvement(0.5, 0.0, 1.0, "maximize") == 0.0
    # Same logic mirrored for minimize: confident-and-lower beats the incumbent.
    assert surrogate.expected_improvement(0.0, 0.0, 1.0, "minimize") == pytest.approx(0.99)
    assert surrogate.expected_improvement(1.5, 0.0, 1.0, "minimize") == 0.0


def test_expected_improvement_grows_with_uncertainty_at_the_incumbent():
    # Sitting exactly on the incumbent, more uncertainty means more upside.
    low_sigma = surrogate.expected_improvement(1.0, 0.2, 1.0, "maximize", xi=0.0)
    high_sigma = surrogate.expected_improvement(1.0, 1.0, 1.0, "maximize", xi=0.0)

    assert high_sigma > low_sigma > 0


# ---------------------------------------------------------------------------
# Surrogate models + cross-validation
# ---------------------------------------------------------------------------

_TRAIN_SMILES = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccccc1O",
                 "CC(=O)O", "CCN", "c1ccncc1"]


def _toy_xy():
    X = np.vstack([surrogate.featurize(s) for s in _TRAIN_SMILES])
    y = np.linspace(-2.0, 1.0, len(_TRAIN_SMILES))
    return X, y


@pytest.mark.parametrize("model", ["rf", "gp"])
def test_fit_surrogate_predicts_mean_and_uncertainty(model):
    X, y = _toy_xy()
    predictor = surrogate.fit_surrogate(X, y, model=model)

    mu, sigma = predictor.predict(X)

    assert mu.shape == (len(X),)
    assert sigma.shape == (len(X),)
    assert np.all(sigma >= 0)


@pytest.mark.parametrize("model", ["rf", "gp"])
def test_cv_metrics_returns_r2_and_mae(model):
    X, y = _toy_xy()

    metrics = surrogate.cv_metrics(X, y, model=model)

    assert metrics is not None
    assert {"r2", "mae", "n"} <= set(metrics)
    assert np.isfinite(metrics["mae"])
    assert metrics["n"] >= 2


def test_cv_metrics_returns_none_for_tiny_data():
    X = np.vstack([surrogate.featurize("CCO"), surrogate.featurize("CCN")])
    assert surrogate.cv_metrics(X, np.array([1.0, 2.0])) is None


# ---------------------------------------------------------------------------
# Novelty + diversity
# ---------------------------------------------------------------------------

def test_novelty_is_zero_for_member_and_high_for_outsider():
    train_fps = [surrogate.morgan_fp(s) for s in _TRAIN_SMILES]

    assert surrogate.novelty(surrogate.morgan_fp("CCO"), train_fps) == pytest.approx(0.0)
    assert surrogate.novelty(surrogate.morgan_fp("FC(F)(F)c1ccc(Br)cc1"), train_fps) > 0.5


def test_novelty_is_one_against_empty_training_set():
    assert surrogate.novelty(surrogate.morgan_fp("CCO"), []) == 1.0


def test_maxmin_select_returns_distinct_indices():
    fps = [surrogate.morgan_fp(s) for s in _TRAIN_SMILES]

    picked = surrogate.maxmin_select(fps, 3)

    assert len(picked) == 3
    assert len(set(picked)) == 3
    assert picked[0] == 0  # seeded by the best-ranked item


def test_maxmin_select_returns_all_indices_when_k_exceeds_pool():
    fps = [surrogate.morgan_fp(s) for s in _TRAIN_SMILES[:3]]

    assert surrogate.maxmin_select(fps, 10) == [0, 1, 2]
    assert surrogate.maxmin_select([], 3) == []


# ---------------------------------------------------------------------------
# Ranking: cold start vs modeled
# ---------------------------------------------------------------------------

def _seed_run(run_id, smiles_scores, *, direction="maximize", constraints=None):
    state.save_config(run_id, {
        "run_id": run_id,
        "objective": "Maximize logS",
        "objective_direction": direction,
        "primary_metric": "logS",
        "workflow_type": "solubility",
        "constraints": constraints or {},
        "max_iterations": 20,
    })
    candidates = [
        {
            "smiles": s,
            "parent_smiles": "",
            "design_note": "seed",
            "score": float(score),
            "metrics": {"logS": float(score)},
            "satisfies_constraints": True,
            "constraint_failures": [],
        }
        for s, score in smiles_scores
    ]
    state.save_iter(run_id, 1, {
        "iter": 1,
        "rationale": "seed",
        "candidates": candidates,
        "iteration_best": None,
        "decision": "continue",
    })


# ---------------------------------------------------------------------------
# History -> training set
# ---------------------------------------------------------------------------

def test_evaluated_smiles_tracks_every_scored_candidate_by_canonical_key(isolated_runs):
    state.save_config("demo", {
        "run_id": "demo", "objective_direction": "maximize",
        "primary_metric": "logS", "workflow_type": "solubility",
        "constraints": {}, "max_iterations": 20,
    })
    state.save_iter("demo", 1, {"candidates": [
        {"smiles": "OCC", "score": -1.0, "satisfies_constraints": True},   # CCO, passes
        {"smiles": "c1ccccc1", "score": 0.5, "satisfies_constraints": False},  # scored but fails
        {"smiles": "CCN", "score": None, "satisfies_constraints": True},   # never hit Rowan
    ]})

    seen = surrogate.evaluated_smiles("demo")

    # Dedup is by canonical SMILES, and constraint-failing-but-scored still counts.
    assert surrogate.canonical("CCO") in seen
    assert surrogate.canonical("c1ccccc1") in seen
    # Unscored candidates never consumed credits, so they are not "seen".
    assert surrogate.canonical("CCN") not in seen


def test_training_set_dedups_and_drops_constraint_failures(isolated_runs):
    state.save_config("demo", {
        "run_id": "demo", "objective_direction": "maximize",
        "primary_metric": "logS", "workflow_type": "solubility",
        "constraints": {}, "max_iterations": 20,
    })
    state.save_iter("demo", 1, {"candidates": [
        {"smiles": "CCO", "score": -2.0, "satisfies_constraints": True},
        {"smiles": "c1ccccc1", "score": 9.0, "satisfies_constraints": False},
    ]})
    # Later iteration re-scores CCO (via a non-canonical spelling): newest wins.
    state.save_iter("demo", 2, {"candidates": [
        {"smiles": "OCC", "score": -1.5, "satisfies_constraints": True},
    ]})

    X, y, fps, canon = surrogate.training_set("demo")

    assert X.shape == (1, surrogate.FP_BITS)
    assert canon == [surrogate.canonical("CCO")]
    assert y.tolist() == [-1.5]  # latest score, not the original -2.0
    assert len(fps) == 1


def test_training_set_is_empty_for_a_fresh_run(isolated_runs):
    state.save_config("demo", {
        "run_id": "demo", "objective_direction": "maximize",
        "primary_metric": "logS", "workflow_type": "solubility",
        "constraints": {}, "max_iterations": 20,
    })

    X, y, fps, canon = surrogate.training_set("demo")

    assert X.shape == (0, surrogate.FP_BITS)
    assert y.shape == (0,)
    assert fps == [] and canon == []


def test_rank_candidates_cold_start_falls_back_to_novelty(isolated_runs):
    _seed_run("demo", [("CCO", -2.0), ("CCN", -1.0)])  # below MIN_TRAIN
    cfg = state.load_config("demo")

    result = surrogate.rank_candidates(
        "demo",
        [
            {"smiles": "CCS", "parent_smiles": "", "design_note": "thiol"},
            {"smiles": "FC(F)(F)c1ccc(Br)cc1", "parent_smiles": "", "design_note": "very different"},
        ],
        cfg=cfg,
    )

    assert result["modeled"] is False
    assert result["rank_key"] == "novelty"
    # The structurally most different candidate should rank first by novelty.
    assert result["ranked_pool"][0]["smiles"] == "FC(F)(F)c1ccc(Br)cc1"


def test_rank_candidates_modeled_uses_ei_and_flags_duplicates(isolated_runs):
    _seed_run("demo", list(zip(_TRAIN_SMILES, np.linspace(-2.0, 1.0, len(_TRAIN_SMILES)))))
    cfg = state.load_config("demo")

    result = surrogate.rank_candidates(
        "demo",
        [
            {"smiles": "OCC", "parent_smiles": "", "design_note": "dup of CCO"},
            {"smiles": "CC#N", "parent_smiles": "", "design_note": "novel"},
        ],
        cfg=cfg,
    )

    assert result["modeled"] is True
    assert result["rank_key"] == "ei"
    by_smiles = {r["smiles"]: r for r in result["records"]}
    assert by_smiles["OCC"]["already_evaluated"] is True
    assert by_smiles["CC#N"]["already_evaluated"] is False
    # Modeled records carry a surrogate prediction, uncertainty, and EI.
    assert {"pred", "pred_sigma", "ei"} <= set(by_smiles["CC#N"])
    # Duplicate is excluded from the recommendable pool.
    assert "OCC" not in {r["smiles"] for r in result["ranked_pool"]}


def test_rank_candidates_include_seen_keeps_duplicates_in_pool(isolated_runs):
    _seed_run("demo", list(zip(_TRAIN_SMILES, np.linspace(-2.0, 1.0, len(_TRAIN_SMILES)))))
    cfg = state.load_config("demo")

    result = surrogate.rank_candidates(
        "demo",
        [{"smiles": "OCC", "parent_smiles": "", "design_note": "dup of CCO"}],
        cfg=cfg,
        include_seen=True,
    )

    pool = {r["smiles"] for r in result["ranked_pool"]}
    assert "OCC" in pool  # not filtered out when include_seen=True


def test_rank_candidates_respect_constraints_toggle_admits_violators(isolated_runs):
    _seed_run(
        "demo",
        list(zip(_TRAIN_SMILES, np.linspace(-2.0, 1.0, len(_TRAIN_SMILES)))),
        constraints={"mw_max": 50},
    )
    cfg = state.load_config("demo")
    heavy = {"smiles": "CCCCCCCCCCCCO", "parent_smiles": "", "design_note": "too heavy"}

    filtered = surrogate.rank_candidates("demo", [heavy], cfg=cfg)
    unfiltered = surrogate.rank_candidates("demo", [heavy], cfg=cfg, respect_constraints=False)

    assert filtered["ranked_pool"] == []  # busts mw_max
    assert {r["smiles"] for r in unfiltered["ranked_pool"]} == {"CCCCCCCCCCCCO"}


# ---------------------------------------------------------------------------
# CLI smoke test (fully local, no Rowan)
# ---------------------------------------------------------------------------

def test_suggest_cli_dedups_filters_and_recommends(isolated_runs):
    _seed_run(
        "demo",
        list(zip(_TRAIN_SMILES, np.linspace(-2.0, 1.0, len(_TRAIN_SMILES)))),
        constraints={"mw_max": 100},
    )

    result = CliRunner().invoke(
        surrogate.cli,
        [
            "--run", "demo",
            "--candidate", "CCO||already tried",
            "--candidate", "not a smiles||garbage",
            "--candidate", "CCCCCCCCCCCCCCCCCCO||too heavy",  # MW ~270 > 100
            "--candidate", "CC#N||novel small nitrile",
            "--top-k", "3",
            "--json",
        ],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)

    assert payload["modeled"] is True
    assert payload["rank_key"] == "ei"

    records = {r.get("canonical_smiles", r["smiles"]): r for r in payload["records"]}
    assert records[surrogate.canonical("CCO")]["already_evaluated"] is True
    assert any(r.get("status") == "invalid" for r in payload["records"])
    heavy = records[surrogate.canonical("CCCCCCCCCCCCCCCCCCO")]
    assert heavy["satisfies_constraints"] is False

    recommended = {r["smiles"] for r in payload["recommended"]}
    assert surrogate.canonical("CC#N") in {
        surrogate.canonical(s) for s in recommended
    }
    # Duplicate, invalid, and constraint-busting candidates are never recommended.
    assert surrogate.canonical("CCO") not in {surrogate.canonical(s) for s in recommended}
    assert surrogate.canonical("CCCCCCCCCCCCCCCCCCO") not in {
        surrogate.canonical(s) for s in recommended
    }


def test_suggest_cli_errors_without_candidates(isolated_runs):
    _seed_run("demo", [("CCO", -2.0)])

    result = CliRunner().invoke(surrogate.cli, ["--run", "demo"])

    assert result.exit_code != 0
    assert "No candidates" in result.output


def test_suggest_cli_generates_candidates_from_parent(isolated_runs):
    _seed_run("demo", list(zip(_TRAIN_SMILES, np.linspace(-2.0, 1.0, len(_TRAIN_SMILES)))))

    result = CliRunner().invoke(
        surrogate.cli,
        ["--run", "demo", "--from-parent", "c1ccccc1",
         "--strategy", "scan-subst", "--n", "6", "--top-k", "3", "--json"],
    )

    assert result.exit_code == 0, result.output
    payload = json.loads(result.output)
    assert payload["modeled"] is True
    assert len(payload["records"]) > 0
    assert 1 <= len(payload["recommended"]) <= 3
    # Generated candidates all trace back to the requested parent.
    assert all(r.get("parent_smiles") == "c1ccccc1" for r in payload["recommended"])


def test_suggest_cli_cold_start_ranks_by_novelty_in_plain_text(isolated_runs):
    _seed_run("demo", [("CCO", -2.0), ("CCN", -1.0)])  # below MIN_TRAIN

    result = CliRunner().invoke(
        surrogate.cli,
        ["--run", "demo",
         "--candidate", "CCS||thiol",
         "--candidate", "FC(F)(F)c1ccc(Br)cc1||very different"],
    )

    assert result.exit_code == 0, result.output
    assert "cold start" in result.output
    assert "NOVELTY" in result.output
    # The structurally distant candidate is recommended ahead of the near one.
    lines = result.output.splitlines()
    rec_lines = [ln for ln in lines if ln.startswith(("CCS", "FC(F)"))]
    assert rec_lines and rec_lines[0].startswith("FC(F)")
