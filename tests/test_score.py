from __future__ import annotations

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


def test_extract_metric_ignores_pinned_path_that_resolves_to_non_number():
    payload = {"solubility": {"logS": "not-a-number"}, "logS": -0.7}

    # Pinned path resolves to a string, so it falls through to the top-level key.
    assert score.extract_metric("solubility", "logS", "solubility.logS", payload) == (-0.7, "logS")


def test_term_value_resolves_local_then_path_then_name_fallbacks():
    # source=local reads from the locally computed descriptor metrics.
    assert score._term_value({"source": "local", "metric": "mw"}, {"mw": 12.0}, {}) == (12.0, "mw")
    # An explicit path walks into Rowan object_data.
    assert score._term_value({"path": "a.b"}, {}, {"a": {"b": 3.0}}) == (3.0, "a.b")
    # With only a name, object_data wins over local metrics.
    assert score._term_value({"name": "x"}, {"x": 7.0}, {"x": 2.0}) == (2.0, "x")
    assert score._term_value({"name": "x"}, {"x": 7.0}, {}) == (7.0, "x")
    # Nothing resolvable -> (None, key).
    assert score._term_value({"name": "x"}, {}, {}) == (None, "x")


def test_workflow_params_accepts_either_config_key():
    assert score._workflow_params({"workflow_params": {"method": "fast"}}) == {"method": "fast"}
    assert score._workflow_params({"rowan_workflow_kwargs": {"solvents": ["water"]}}) == {
        "solvents": ["water"]
    }
    assert score._workflow_params({}) == {}


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


def test_check_constraints_enforces_min_bounds_and_passes_when_satisfied():
    failing_ok, failing = score.check_constraints({"mw": 80.0, "heavy_atoms": 3}, {"mw_min": 100, "ha_min": 5})
    assert failing_ok is False
    assert failing == ["mw_min: mw=80.00 < 100", "ha_min: heavy_atoms=3.00 < 5"]

    passing_ok, passing = score.check_constraints(
        {"mw": 150.0, "tpsa": 20.0, "hbd": 1}, {"mw_max": 200, "mw_min": 100, "tpsa_max": 40}
    )
    assert passing_ok is True
    assert passing == []


def test_check_constraints_skips_rules_with_missing_descriptors():
    # logp_crippen absent from the descriptor dict -> the rule is skipped, not failed.
    ok, failed = score.check_constraints({"mw": 120.0}, {"logp_max": 3})

    assert ok is True
    assert failed == []


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


def test_submit_one_applies_solubility_defaults_and_signature_dispatch(monkeypatch):
    captured = {}

    class FakeWorkflow:
        uuid = "sol-1"

        def result(self, wait, poll_interval):
            assert wait is True and poll_interval == 5
            return types.SimpleNamespace(data={"solubility": {"logS": -1.0}})

    def submit_solubility_workflow(initial_smiles=None, name=None, method=None,
                                   solvents=None, temperatures=None):
        captured.update(
            initial_smiles=initial_smiles, name=name, method=method,
            solvents=solvents, temperatures=temperatures,
        )
        return FakeWorkflow()

    monkeypatch.setitem(
        sys.modules, "rowan",
        types.SimpleNamespace(submit_solubility_workflow=submit_solubility_workflow),
    )

    result = score.submit_one(
        "CCO",
        base_config(workflow_params={"method": "fastsolv"}),
        "job-name",
    )

    assert result == {
        "uuid": "sol-1",
        "object_data": {"solubility": {"logS": -1.0}},
        "workflow_type": "solubility",
    }
    # User-supplied params override the kingfisher defaults; the rest stay.
    assert captured["initial_smiles"] == "CCO"
    assert captured["name"] == "job-name"
    assert captured["method"] == "fastsolv"
    assert captured["solvents"] == ["water"]
    assert captured["temperatures"] == [298.15]


def test_submit_one_dispatches_to_named_workflow_for_non_solubility(monkeypatch):
    captured = {}

    class FakeWorkflow:
        uuid = "pka-1"

        def result(self, wait, poll_interval):
            return types.SimpleNamespace(data={"pka": 4.2})

    def submit_pka_workflow(initial_smiles=None, name=None):
        captured.update(initial_smiles=initial_smiles, name=name)
        return FakeWorkflow()

    monkeypatch.setitem(
        sys.modules, "rowan",
        types.SimpleNamespace(submit_pka_workflow=submit_pka_workflow),
    )

    result = score.submit_one("CCO", base_config(workflow_type="pka", primary_metric="pka"), "j")

    assert result["uuid"] == "pka-1"
    assert result["object_data"] == {"pka": 4.2}
    assert captured == {"initial_smiles": "CCO", "name": "j"}


def test_submit_one_raises_for_workflow_without_a_submitter(monkeypatch):
    monkeypatch.setitem(sys.modules, "rowan", types.SimpleNamespace())

    with pytest.raises(AttributeError, match="no submit function for workflow 'mystery'"):
        score.submit_one("CCO", base_config(workflow_type="mystery"), "j")


def test_score_candidate_flags_unextractable_metric(monkeypatch):
    monkeypatch.setattr(score, "local_descriptors", lambda smiles: {"mw": 46.0})
    monkeypatch.setattr(
        score, "submit_one",
        lambda smiles, cfg, name: {"uuid": "wf-1", "object_data": {"unrelated": 1}, "workflow_type": "solubility"},
    )

    result = score.score_candidate({"smiles": "CCO"}, base_config(), 1)

    assert result["score"] is None
    assert result["metric_path_used"] is None
    assert "Could not extract metric 'logS'" in result["error"]
    # The raw payload is retained so the agent can fix metric_path.
    assert result["object_data"] == {"unrelated": 1}


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


def test_score_batch_docking_short_circuits_when_all_smiles_invalid(monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("Rowan should not be called when nothing is valid")

    monkeypatch.setitem(sys.modules, "rowan", types.SimpleNamespace(submit_batch_docking_workflow=fail))
    monkeypatch.setattr(score, "local_descriptors", lambda smiles: {})

    results = score.score_batch_docking(
        [{"smiles": "bad"}], base_config(workflow_type="batch_docking"), 1
    )

    assert len(results) == 1
    assert results[0]["score"] is None
    assert results[0]["error"] == "RDKit failed to parse SMILES"


def test_score_batch_docking_propagates_submit_failure_to_each_candidate(monkeypatch):
    def boom(smiles_list, name, **params):
        raise RuntimeError("docking service down")

    monkeypatch.setitem(
        sys.modules, "rowan",
        types.SimpleNamespace(submit_batch_docking_workflow=boom),
    )
    monkeypatch.setattr(score, "local_descriptors", lambda smiles: {"mw": 46.0})

    results = score.score_batch_docking(
        [{"smiles": "CCO"}, {"smiles": "CCN"}],
        base_config(workflow_type="batch_docking", primary_metric="docking_score"),
        2,
    )

    assert [r["score"] for r in results] == [None, None]
    assert all("docking service down" in r["error"] for r in results)
    assert all("traceback" in r for r in results)


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


def test_score_cli_rejects_empty_smiles(isolated_runs):
    state.save_config("demo", base_config())

    result = CliRunner().invoke(
        score.cli,
        ["--run", "demo", "--rationale", "oops", "--candidate", "|CC|empty smiles"],
    )

    assert result.exit_code != 0
    assert "Empty SMILES" in result.output


def test_score_cli_refuses_decision_before_any_iteration(isolated_runs):
    state.save_config("demo", base_config())

    result = CliRunner().invoke(
        score.cli,
        ["--run", "demo", "--rationale", "converged already?", "--decision", "converged"],
    )

    assert result.exit_code != 0
    assert "Cannot record a decision before any iterations exist" in result.output
    assert state.list_iterations("demo") == []


def test_score_cli_dry_run_without_candidates_is_a_noop(isolated_runs):
    state.save_config("demo", base_config())

    result = CliRunner().invoke(
        score.cli,
        ["--run", "demo", "--rationale", "nothing yet", "--dry-run"],
    )

    assert result.exit_code == 0, result.output
    assert "no candidates to score" in result.output
    assert state.list_iterations("demo") == []


def test_score_cli_routes_batch_docking_workflow(isolated_runs, monkeypatch):
    state.save_config("demo", base_config(workflow_type="batch_docking", primary_metric="docking_score"))
    monkeypatch.setenv("ROWAN_API_KEY", "test-key")
    monkeypatch.setattr(score, "rebuild_report", lambda run_id: state.run_dir(run_id) / "report.html")

    def fake_batch(candidates, cfg, iter_n):
        return [
            {**candidates[0], "score": -9.3, "metrics": {"docking_score": -9.3},
             "satisfies_constraints": True, "constraint_failures": []},
        ]

    monkeypatch.setattr(score, "score_batch_docking", fake_batch)

    result = CliRunner().invoke(
        score.cli,
        ["--run", "demo", "--rationale", "Dock one ligand", "--candidate", "CCO||dock"],
    )

    assert result.exit_code == 0, result.output
    iteration = state.load_iter("demo", 1)
    assert iteration["iteration_best"] == {"smiles": "CCO", "score": -9.3}
    assert "docking_score=-9.300" in result.output
