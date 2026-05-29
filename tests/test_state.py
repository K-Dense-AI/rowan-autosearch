from __future__ import annotations

import json
from datetime import datetime

import click
from click.testing import CliRunner
import pytest

from rowan_tools import state


# ---------------------------------------------------------------------------
# Pure helpers
# ---------------------------------------------------------------------------

def test_now_iso_is_timezone_aware_and_second_resolution():
    stamp = state.now_iso()
    parsed = datetime.fromisoformat(stamp)

    assert parsed.tzinfo is not None  # always UTC-aware
    assert parsed.microsecond == 0  # truncated to seconds


def test_parse_jsonish_handles_json_numbers_and_plain_strings():
    assert state.parse_jsonish('["water", "ethanol"]') == ["water", "ethanol"]
    assert state.parse_jsonish('{"method": "fast"}') == {"method": "fast"}
    assert state.parse_jsonish("298.15") == 298.15
    assert state.parse_jsonish("kingfisher") == "kingfisher"


def test_parse_key_value_options_rejects_malformed_values():
    assert state.parse_key_value_options(("mw_max=250", "label=acid"), "--constraint") == {
        "mw_max": 250.0,
        "label": "acid",
    }

    with pytest.raises(click.ClickException, match="expected key=value"):
        state.parse_key_value_options(("mw_max",), "--constraint")


def test_init_cli_creates_config_and_directories(isolated_runs):
    result = CliRunner().invoke(
        state.cli,
        [
            "init",
            "--run",
            "demo",
            "--objective",
            "Maximize logS",
            "--direction",
            "maximize",
            "--metric",
            "logS",
            "--workflow",
            "solubility",
            "--start-smiles",
            "CCO",
            "--metric-path",
            "solubility.logS",
            "--constraint",
            "mw_max=250",
            "--workflow-param",
            'solvents=["water"]',
            "--objective-term",
            '{"name":"mw","source":"local","metric":"mw","goal":"minimize","weight":0.01}',
            "--max-iter",
            "7",
            "--candidates-per-iter",
            "3",
        ],
    )

    assert result.exit_code == 0, result.output
    cfg = json.loads((isolated_runs / "demo" / "config.json").read_text())
    assert cfg["run_id"] == "demo"
    assert cfg["constraints"] == {"mw_max": 250.0}
    assert cfg["workflow_params"] == {"solvents": ["water"]}
    assert cfg["optimization_objective"]["terms"][0]["name"] == "mw"
    assert (isolated_runs / "demo" / "iterations").is_dir()
    assert (isolated_runs / "demo" / "plots").is_dir()


def test_init_cli_requires_force_for_existing_config(isolated_runs):
    args = [
        "init",
        "--run",
        "demo",
        "--objective",
        "Minimize score",
        "--direction",
        "minimize",
        "--metric",
        "score",
        "--workflow",
        "descriptors",
        "--start-smiles",
        "CCO",
    ]

    runner = CliRunner()
    assert runner.invoke(state.cli, args).exit_code == 0
    cfg = json.loads((isolated_runs / "demo" / "config.json").read_text())
    assert cfg["max_iterations"] == 50
    result = runner.invoke(state.cli, args)

    assert result.exit_code != 0
    assert "Use --force to overwrite" in result.output


def test_iteration_listing_and_best_so_far_are_constraint_and_direction_aware(isolated_runs):
    cfg = {
        "run_id": "demo",
        "objective": "Optimize",
        "objective_direction": "maximize",
        "primary_metric": "logS",
        "workflow_type": "solubility",
        "max_iterations": 5,
    }
    state.save_config("demo", cfg)
    state.save_iter(
        "demo",
        2,
        {
            "candidates": [
                {
                    "smiles": "CCO",
                    "score": 1.0,
                    "metrics": {"logS": 1.0},
                    "satisfies_constraints": True,
                },
                {
                    "smiles": "CCCCCCCC",
                    "score": 10.0,
                    "metrics": {"logS": 10.0},
                    "satisfies_constraints": False,
                },
            ]
        },
    )
    state.save_iter(
        "demo",
        1,
        {
            "candidates": [
                {
                    "smiles": "CO",
                    "score": 2.0,
                    "metrics": {"logS": 2.0},
                    "satisfies_constraints": True,
                }
            ]
        },
    )
    (isolated_runs / "demo" / "iterations" / "notes.json").write_text("{}")

    assert state.list_iterations("demo") == [1, 2]
    best = state.best_so_far("demo")
    assert best["smiles"] == "CO"
    assert best["iter"] == 1

    cfg["objective_direction"] = "minimize"
    state.save_config("demo", cfg)
    assert state.best_so_far("demo")["smiles"] == "CCO"


def test_status_and_list_cli_render_run_summary(isolated_runs):
    state.save_config(
        "demo",
        {
            "run_id": "demo",
            "objective": "Maximize logS",
            "objective_direction": "maximize",
            "primary_metric": "logS",
            "workflow_type": "solubility",
            "max_iterations": 3,
        },
    )
    state.save_iter(
        "demo",
        1,
        {
            "candidates": [
                {
                    "smiles": "CCO",
                    "score": 1.25,
                    "metrics": {"logS": 1.25},
                    "satisfies_constraints": True,
                }
            ]
        },
    )

    status = CliRunner().invoke(state.cli, ["status", "--run", "demo"])
    listing = CliRunner().invoke(state.cli, ["list"])

    assert status.exit_code == 0
    assert "Best so far:   CCO" in status.output
    assert "Found at iter: 1" in status.output
    assert listing.exit_code == 0
    assert "demo" in listing.output


def test_status_uses_score_when_primary_metric_key_is_absent(isolated_runs):
    state.save_config(
        "demo",
        {
            "run_id": "demo",
            "objective": "Optimize composite",
            "objective_direction": "maximize",
            "primary_metric": "composite",
            "workflow_type": "solubility",
            "max_iterations": 3,
        },
    )
    state.save_iter(
        "demo",
        1,
        {
            "candidates": [
                {
                    "smiles": "CCO",
                    "score": 1.25,
                    "metrics": {"logS": -1.0, "mw": 46.0},
                    "satisfies_constraints": True,
                }
            ]
        },
    )

    status = CliRunner().invoke(state.cli, ["status", "--run", "demo"])

    assert status.exit_code == 0
    assert "composite=1.25" in status.output


# ---------------------------------------------------------------------------
# Config + candidate aggregation internals
# ---------------------------------------------------------------------------

def test_load_config_raises_clear_error_when_missing(isolated_runs):
    with pytest.raises(FileNotFoundError, match="rowan-state init"):
        state.load_config("does-not-exist")


def _seed_two_iterations(run_id="demo", direction="maximize"):
    state.save_config(run_id, {
        "run_id": run_id,
        "objective": "Optimize",
        "objective_direction": direction,
        "primary_metric": "logS",
        "workflow_type": "solubility",
        "max_iterations": 5,
    })
    state.save_iter(run_id, 1, {"candidates": [
        {"smiles": "CCO", "score": 1.0, "satisfies_constraints": True},
    ]})
    state.save_iter(run_id, 2, {"candidates": [
        {"smiles": "CO", "score": 2.0, "satisfies_constraints": True},
        {"smiles": "CN", "score": None, "satisfies_constraints": True},
    ]})


def test_all_candidates_flattens_every_iteration_and_tags_origin(isolated_runs):
    _seed_two_iterations()

    cands = state.all_candidates("demo")

    assert [c["smiles"] for c in cands] == ["CCO", "CO", "CN"]
    assert [c["iter"] for c in cands] == [1, 2, 2]


def test_best_so_far_returns_none_when_no_candidate_qualifies(isolated_runs):
    state.save_config("demo", {
        "run_id": "demo",
        "objective": "Optimize",
        "objective_direction": "maximize",
        "primary_metric": "logS",
        "workflow_type": "solubility",
        "max_iterations": 5,
    })
    # One unscored, one scored-but-constraint-failing: neither is eligible.
    state.save_iter("demo", 1, {"candidates": [
        {"smiles": "CCO", "score": None, "satisfies_constraints": True},
        {"smiles": "CCCCCCCC", "score": 9.0, "satisfies_constraints": False},
    ]})

    assert state.best_so_far("demo") is None


# ---------------------------------------------------------------------------
# init CLI edge cases
# ---------------------------------------------------------------------------

def _init_args(**overrides):
    args = {
        "--run": "demo",
        "--objective": "Maximize logS",
        "--direction": "maximize",
        "--metric": "logS",
        "--workflow": "solubility",
        "--start-smiles": "CCO",
    }
    args.update(overrides)
    flat = ["init"]
    for k, v in args.items():
        flat.extend([k, str(v)])
    return flat


def test_init_force_overwrites_existing_config(isolated_runs):
    runner = CliRunner()
    assert runner.invoke(state.cli, _init_args()).exit_code == 0

    result = runner.invoke(
        state.cli,
        _init_args(**{"--objective": "Minimize logP", "--max-iter": 99}) + ["--force"],
    )

    assert result.exit_code == 0, result.output
    cfg = state.load_config("demo")
    assert cfg["objective"] == "Minimize logP"
    assert cfg["max_iterations"] == 99


def test_init_rejects_objective_term_that_is_not_a_json_object(isolated_runs):
    result = CliRunner().invoke(
        state.cli,
        _init_args() + ["--objective-term", "[1, 2, 3]"],
    )

    assert result.exit_code != 0
    assert "expected a JSON object" in result.output


def test_init_omits_optional_sections_when_not_supplied(isolated_runs):
    assert CliRunner().invoke(state.cli, _init_args()).exit_code == 0

    cfg = state.load_config("demo")
    assert "workflow_params" not in cfg
    assert "optimization_objective" not in cfg
    assert cfg["candidates_per_iter"] == 4


# ---------------------------------------------------------------------------
# status / list CLI fallbacks
# ---------------------------------------------------------------------------

def test_status_reports_no_best_before_any_eligible_candidate(isolated_runs):
    state.save_config("demo", {
        "run_id": "demo",
        "objective": "Maximize logS",
        "objective_direction": "maximize",
        "primary_metric": "logS",
        "workflow_type": "solubility",
        "max_iterations": 3,
    })

    status = CliRunner().invoke(state.cli, ["status", "--run", "demo"])

    assert status.exit_code == 0
    assert "Best so far:   (none yet)" in status.output
    assert "Iterations:    0 / 3" in status.output


def test_list_runs_reports_empty_when_no_runs_directory(isolated_runs):
    result = CliRunner().invoke(state.cli, ["list"])

    assert result.exit_code == 0
    assert "No runs yet." in result.output
