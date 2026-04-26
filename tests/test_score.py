from __future__ import annotations

import os
import sys
import types

from click.testing import CliRunner
import pytest

from rowan_tools import score, state


def base_config(**overrides):
    cfg = {
        "run_id": "demo",
        "objective": "Maximize logS",
        "objective_direction": "maximize",
        "primary_metric": "logS",
        "workflow_type": "solubility",
        "metric_path": None,
        "constraints": {},
        "max_iterations": 5,
    }
    cfg.update(overrides)
    return cfg


def test_walk_supports_nested_dicts_and_list_indexes():
    payload = {"a": {"b": [{"c": 1.5}]}}

    assert score.walk(payload, "a.b.0.c") == 1.5
    assert score.walk(payload, "a.b.2.c") is None
    assert score.walk(payload, "a.b.nope.c") is None
    assert score.walk(payload, "a.b.0.c.extra") is None


def test_extract_metric_prefers_pinned_path_then_defaults_and_top_level():
    payload = {
        "custom": {"value": -2.5},
        "solubilities": {"O": {"solubilities": [-1.2]}},
        "logS": -0.3,
    }

    assert score.extract_metric("solubility", "logS", "custom.value", payload) == (
        -2.5,
        "custom.value",
    )
    assert score.extract_metric("solubility", "logS", "missing.path", payload) == (
        -1.2,
        "solubilities.O.solubilities.0",
    )
    assert score.extract_metric("unknown", "logS", None, payload) == (-0.3, "logS")
    assert score.extract_metric("unknown", "missing", None, payload) == (None, None)


def test_evaluate_optimization_objective_combines_rowan_local_and_target_terms():
    cfg = base_config(
        optimization_objective={
            "name": "composite",
            "terms": [
                {"name": "logS", "path": "solubility.logS", "goal": "maximize", "weight": 2},
                {"name": "mw", "source": "local", "metric": "mw", "goal": "minimize", "weight": 0.01},
                {"name": "pka", "path": "pka", "target": 4.5, "weight": 3},
            ],
        }
    )

    value, details, error = score.evaluate_optimization_objective(
        cfg,
        {"mw": 180.0},
        {"solubility": {"logS": -1.0}, "pka": 5.0},
    )

    assert error is None
    assert value == pytest.approx(-5.3)
    assert details["logS"]["contribution"] == pytest.approx(-2.0)
    assert details["mw"]["contribution"] == pytest.approx(-1.8)
    assert details["pka"]["contribution"] == pytest.approx(-1.5)


def test_evaluate_optimization_objective_reports_unresolved_terms():
    value, details, error = score.evaluate_optimization_objective(
        base_config(optimization_objective={"terms": [{"name": "missing", "path": "nope"}]}),
        {},
        {},
    )

    assert value is None
    assert details == {}
    assert "could not be resolved" in error


def test_evaluate_optimization_objective_reports_missing_target_value():
    value, details, error = score.evaluate_optimization_objective(
        base_config(optimization_objective={"terms": [{"name": "pka", "path": "pka", "goal": "target"}]}),
        {},
        {"pka": 5.0},
    )

    assert value is None
    assert details == {}
    assert "goal 'target' but no target value" in error


def test_check_constraints_ignores_unknown_rules_and_reports_failures():
    ok, failed = score.check_constraints(
        {"mw": 301.0, "logp_crippen": 4.2, "hbd": 1},
        {"mw_max": 300, "logp_max": 5, "unknown": 0},
    )

    assert ok is False
    assert failed == ["mw_max: mw=301.00 > 300"]


def test_local_descriptors_returns_empty_dict_for_invalid_smiles():
    assert score.local_descriptors("not a smiles") == {}


def test_score_candidate_handles_invalid_smiles_without_calling_rowan(monkeypatch):
    def fail_submit(*args, **kwargs):
        raise AssertionError("submit_one should not be called")

    monkeypatch.setattr(score, "submit_one", fail_submit)

    result = score.score_candidate(
        {"smiles": "not a smiles", "parent_smiles": "", "design_note": "bad"},
        base_config(),
        1,
    )

    assert result["score"] is None
    assert result["satisfies_constraints"] is False
    assert result["constraint_failures"] == ["invalid SMILES"]
    assert "RDKit failed" in result["error"]


def test_score_candidate_records_rowan_payload_metric_and_objective(monkeypatch):
    monkeypatch.setattr(score.state, "now_iso", lambda: "2026-01-01T00:00:00+00:00")
    monkeypatch.setattr(
        score,
        "local_descriptors",
        lambda smiles: {"mw": 100.0, "logp_crippen": 1.0},
    )
    monkeypatch.setattr(
        score,
        "submit_one",
        lambda smiles, cfg, name: {
            "uuid": "wf-123",
            "object_data": {"solubility": {"logS": -1.5}},
            "workflow_type": "solubility",
        },
    )

    result = score.score_candidate(
        {"smiles": "CCO", "parent_smiles": "CC", "design_note": "add OH"},
        base_config(
            metric_path="solubility.logS",
            constraints={"mw_max": 120},
            optimization_objective={
                "name": "objective",
                "terms": [
                    {"name": "logS", "path": "solubility.logS", "weight": 1},
                    {"name": "mw", "source": "local", "metric": "mw", "goal": "minimize", "weight": 0.01},
                ],
            },
        ),
        3,
    )

    assert result["workflow_uuid"] == "wf-123"
    assert result["metric_path_used"] == "solubility.logS"
    assert result["metrics"]["logS"] == -1.5
    assert result["score"] == pytest.approx(-2.5)
    assert result["optimization_objective"]["mw"]["contribution"] == pytest.approx(-1.0)
    assert result["satisfies_constraints"] is True


def test_score_candidate_preserves_submit_errors(monkeypatch):
    monkeypatch.setattr(score, "local_descriptors", lambda smiles: {"mw": 50.0})

    def fail_submit(*args, **kwargs):
        raise RuntimeError("network unavailable")

    monkeypatch.setattr(score, "submit_one", fail_submit)

    result = score.score_candidate({"smiles": "CCO"}, base_config(), 1)

    assert result["score"] is None
    assert "Rowan submit failed: network unavailable" in result["error"]
    assert "traceback" in result


def test_score_batch_docking_maps_scores_and_invalid_smiles(monkeypatch):
    class FakeResult:
        data = {"scores": {"CCO": -8.1}}
        scores = {"CCO": -8.1}

    class FakeWorkflow:
        uuid = "batch-123"

        def result(self, wait, poll_interval):
            assert wait is True
            assert poll_interval == 5
            return FakeResult()

    fake_rowan = types.SimpleNamespace(
        submit_batch_docking_workflow=lambda smiles_list, name, **params: FakeWorkflow()
    )
    monkeypatch.setitem(sys.modules, "rowan", fake_rowan)
    monkeypatch.setattr(score, "local_descriptors", lambda smiles: {} if smiles == "bad" else {"mw": 46.0})

    results = score.score_batch_docking(
        [{"smiles": "bad"}, {"smiles": "CCO", "parent_smiles": "CC", "design_note": "dock"}],
        base_config(workflow_type="batch_docking", primary_metric="docking_score"),
        1,
    )

    assert results[0]["error"] == "RDKit failed to parse SMILES"
    assert results[1]["workflow_uuid"] == "batch-123"
    assert results[1]["score"] == -8.1
    assert results[1]["metric_path_used"] == "scores.CCO"


def test_score_cli_dry_run_validates_candidates_without_api_key(isolated_runs, monkeypatch):
    monkeypatch.delenv("ROWAN_API_KEY", raising=False)
    state.save_config("demo", base_config())

    result = CliRunner().invoke(
        score.cli,
        [
            "--run",
            "demo",
            "--rationale",
            "Smoke test",
            "--candidate",
            "CCO|CC|add oxygen",
            "--dry-run",
        ],
    )

    assert result.exit_code == 0
    assert "(dry-run; not calling Rowan)" in result.output
    assert state.list_iterations("demo") == []


def test_score_cli_updates_latest_decision_without_candidates(isolated_runs, monkeypatch):
    monkeypatch.delenv("ROWAN_API_KEY", raising=False)
    state.save_config("demo", base_config())
    state.save_iter(
        "demo",
        1,
        {
            "iter": 1,
            "rationale": "Initial pass",
            "candidates": [{"smiles": "CCO", "score": 1.0, "satisfies_constraints": True}],
            "iteration_best": {"smiles": "CCO", "score": 1.0},
            "decision": "continue",
        },
    )
    monkeypatch.setattr(score, "rebuild_report", lambda run_id: state.run_dir(run_id) / "report.html")

    result = CliRunner().invoke(
        score.cli,
        [
            "--run",
            "demo",
            "--rationale",
            "No further variants are worth testing.",
            "--decision",
            "converged",
        ],
    )

    assert result.exit_code == 0, result.output
    assert state.list_iterations("demo") == [1]
    iteration = state.load_iter("demo", 1)
    assert iteration["decision"] == "converged"
    assert iteration["decision_rationale"] == "No further variants are worth testing."
    assert "Updated decision for iteration 1 to converged" in result.output


def test_score_cli_appends_wrapup_iteration_without_candidates(isolated_runs, monkeypatch):
    monkeypatch.delenv("ROWAN_API_KEY", raising=False)
    state.save_config("demo", base_config())
    state.save_iter(
        "demo",
        1,
        {
            "iter": 1,
            "rationale": "Initial pass",
            "candidates": [{"smiles": "CCO", "score": 1.0, "satisfies_constraints": True}],
            "iteration_best": {"smiles": "CCO", "score": 1.0},
            "decision": "converged",
        },
    )
    monkeypatch.setattr(score, "rebuild_report", lambda run_id: state.run_dir(run_id) / "report.html")

    result = CliRunner().invoke(
        score.cli,
        [
            "--run",
            "demo",
            "--rationale",
            "Final summary of the optimization path.",
        ],
    )

    assert result.exit_code == 0, result.output
    assert state.list_iterations("demo") == [1, 2]
    wrapup = state.load_iter("demo", 2)
    assert wrapup["rationale"] == "Final summary of the optimization path."
    assert wrapup["candidates"] == []
    assert wrapup["iteration_best"] is None
    assert wrapup["decision"] == "converged"
    assert "Saved wrap-up iteration" in result.output


def test_score_cli_writes_iteration_and_best_so_far(isolated_runs, monkeypatch):
    state.save_config("demo", base_config(constraints={"mw_max": 200}))
    monkeypatch.setenv("ROWAN_API_KEY", "test-key")
    monkeypatch.setattr(
        score,
        "score_candidate",
        lambda cand, cfg, iter_n: {
            **cand,
            "score": 2.0 if cand["smiles"] == "CCO" else 1.0,
            "metrics": {"logS": 2.0 if cand["smiles"] == "CCO" else 1.0},
            "satisfies_constraints": True,
            "constraint_failures": [],
        },
    )

    result = CliRunner().invoke(
        score.cli,
        [
            "--run",
            "demo",
            "--rationale",
            "Compare alcohols",
            "--candidate",
            "CCO|CC|add oxygen",
            "--candidate",
            "CO|C|smaller alcohol",
            "--max-workers",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output
    iteration = state.load_iter("demo", 1)
    assert iteration["iteration_best"] == {"smiles": "CCO", "score": 2.0}
    assert iteration["decision"] == "continue"
    assert "Best so far: CCO" in result.output


def test_score_cli_requires_api_key_outside_dry_run(isolated_runs, monkeypatch):
    state.save_config("demo", base_config())
    monkeypatch.delenv("ROWAN_API_KEY", raising=False)

    result = CliRunner().invoke(
        score.cli,
        ["--run", "demo", "--rationale", "Try one", "--candidate", "CCO||note"],
    )

    assert result.exit_code != 0
    assert "ROWAN_API_KEY not set" in result.output
