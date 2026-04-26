from __future__ import annotations

import json

import click
from click.testing import CliRunner
import pytest

from rowan_tools import state


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
