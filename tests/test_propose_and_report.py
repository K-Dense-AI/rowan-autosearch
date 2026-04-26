from __future__ import annotations

from click.testing import CliRunner

from rowan_tools import propose, state
from viz import build_report


def test_proposal_helpers_return_empty_for_invalid_smiles():
    assert propose.bioisostere_variants("not a smiles") == []
    assert propose.aromatic_substitution_variants("not a smiles") == []
    assert propose.brics_variants("not a smiles") == []


def test_aromatic_substitution_variants_are_unique_and_limited():
    variants = propose.aromatic_substitution_variants("c1ccccc1", n=3)

    assert len(variants) <= 3
    assert len({smiles for smiles, _ in variants}) == len(variants)
    assert all("add -" in note for _, note in variants)


def test_propose_cli_prints_pipe_delimited_candidates():
    result = CliRunner().invoke(
        propose.cli,
        ["--smiles", "c1ccccc1", "--strategy", "scan-subst", "--n", "2", "--seed", "1"],
    )

    assert result.exit_code == 0, result.output
    lines = [line for line in result.output.splitlines() if line.strip()]
    assert 1 <= len(lines) <= 2
    for line in lines:
        parts = line.split("|")
        assert len(parts) == 3
        assert parts[1] == "c1ccccc1"
        assert parts[2].startswith("add -")


def test_plot_functions_create_placeholder_images_for_empty_data(tmp_path):
    cfg = {
        "objective": "Maximize logS",
        "objective_direction": "maximize",
        "primary_metric": "logS",
    }
    outputs = [
        tmp_path / "progress.png",
        tmp_path / "pareto.png",
        tmp_path / "grid.png",
        tmp_path / "genealogy.png",
    ]

    build_report.plot_progress(cfg, [], outputs[0])
    build_report.plot_pareto(cfg, [], outputs[1])
    build_report.plot_grid(cfg, [], outputs[2])
    build_report.plot_genealogy(cfg, [], outputs[3])

    assert all(path.exists() and path.stat().st_size > 0 for path in outputs)


def test_plot_grid_handles_scored_but_undrawable_molecules(tmp_path):
    cfg = {
        "objective": "Maximize logS",
        "objective_direction": "maximize",
        "primary_metric": "logS",
    }

    out = tmp_path / "grid.png"
    build_report.plot_grid(
        cfg,
        [{"iter": 1, "smiles": "not a smiles", "score": 1.0, "satisfies_constraints": True}],
        out,
    )

    assert out.exists() and out.stat().st_size > 0


def test_plot_genealogy_handles_duplicate_smiles(tmp_path):
    cfg = {
        "objective": "Maximize logS",
        "objective_direction": "maximize",
        "primary_metric": "logS",
    }
    out = tmp_path / "genealogy.png"

    build_report.plot_genealogy(
        cfg,
        [
            {"iter": 1, "smiles": "CCO", "score": 1.0, "satisfies_constraints": True},
            {"iter": 2, "smiles": "CCO", "parent_smiles": "CCO", "score": 1.5, "satisfies_constraints": True},
            {"iter": 2, "smiles": "CO", "parent_smiles": "CCO", "score": 0.5, "satisfies_constraints": True},
        ],
        out,
    )

    assert out.exists() and out.stat().st_size > 0


def test_best_sdf_returns_empty_for_missing_or_invalid_best():
    assert build_report.best_sdf(None) == ""
    assert build_report.best_sdf({"smiles": "not a smiles"}) == ""


def test_report_cli_renders_html_from_isolated_run(isolated_runs, monkeypatch):
    state.save_config(
        "demo",
        {
            "run_id": "demo",
            "objective": "Maximize logS",
            "objective_direction": "maximize",
            "primary_metric": "logS",
            "workflow_type": "solubility",
            "max_iterations": 4,
        },
    )
    state.save_iter(
        "demo",
        1,
        {
            "iter": 1,
            "rationale": "Initial alcohol check",
            "decision": "continue",
            "iteration_best": {"smiles": "CCO", "score": 1.2},
            "candidates": [
                {
                    "smiles": "CCO",
                    "score": 1.2,
                    "metrics": {"logS": 1.2, "mw": 46.0},
                    "satisfies_constraints": True,
                    "constraint_failures": [],
                    "design_note": "baseline",
                }
            ],
        },
    )
    monkeypatch.setattr(build_report, "best_sdf", lambda best: "")

    result = CliRunner().invoke(build_report.cli, ["--run", "demo"])

    assert result.exit_code == 0, result.output
    report = isolated_runs / "demo" / "report.html"
    assert report.exists()
    html = report.read_text()
    assert "demo" in html
    assert "Rowan" in html
    assert "K-Dense" in html
    assert "Initial alcohol check" in html
    assert "CCO" in html
    for name in ["progress.png", "pareto.png", "grid.png", "genealogy.png"]:
        assert (isolated_runs / "demo" / "plots" / name).exists()
    for name in ["rowan-logo.svg", "rowan-mark.svg", "kdense-logo.png"]:
        assert (isolated_runs / "demo" / "brand" / name).exists()
