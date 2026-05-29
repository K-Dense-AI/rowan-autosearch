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
    # Duplicate is excluded from the recommendable pool.
    assert "OCC" not in {r["smiles"] for r in result["ranked_pool"]}


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
