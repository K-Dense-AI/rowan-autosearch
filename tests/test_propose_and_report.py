from __future__ import annotations

from click.testing import CliRunner

from rowan_tools import propose, state, surrogate
from viz import build_report


def _base_cfg(**overrides):
    cfg = {
        "objective": "Maximize logS",
        "objective_direction": "maximize",
        "primary_metric": "logS",
    }
    cfg.update(overrides)
    return cfg


def _scored(smiles, score, iter_n, *, mw=None, parent=None, ok=True):
    metrics = {"logS": score}
    if mw is not None:
        metrics["mw"] = mw
    rec = {
        "iter": iter_n,
        "smiles": smiles,
        "score": score,
        "metrics": metrics,
        "satisfies_constraints": ok,
    }
    if parent is not None:
        rec["parent_smiles"] = parent
    return rec


def _is_valid_smiles(smiles):
    from rdkit import Chem

    return Chem.MolFromSmiles(smiles) is not None


# ---------------------------------------------------------------------------
# rowan-propose: mutation strategies
# ---------------------------------------------------------------------------

def test_proposal_helpers_return_empty_for_invalid_smiles():
    assert propose.bioisostere_variants("not a smiles") == []
    assert propose.aromatic_substitution_variants("not a smiles") == []
    assert propose.brics_variants("not a smiles") == []


def test_aromatic_substitution_variants_are_unique_and_limited():
    variants = propose.aromatic_substitution_variants("c1ccccc1", n=3)

    assert len(variants) <= 3
    assert len({smiles for smiles, _ in variants}) == len(variants)
    assert all("add -" in note for _, note in variants)


def test_bioisostere_variants_apply_known_swaps_and_stay_valid():
    variants = propose.bioisostere_variants("CC(=O)O", n=8)
    parent_canon = surrogate.canonical("CC(=O)O")

    assert variants  # carboxylic acid has several registered bioisosteres
    smis = {smiles for smiles, _ in variants}
    assert parent_canon not in smis  # never re-proposes the parent
    assert all(_is_valid_smiles(smiles) for smiles, _ in variants)
    assert all("." not in smiles for smiles, _ in variants)  # no fragmented products
    assert all(note.startswith("bioisostere:") for _, note in variants)


def test_brics_variants_recombine_fragments_into_new_molecules():
    parent = "CC(C)Cc1ccc(C(C)C(=O)O)cc1"  # decomposes into multiple BRICS fragments
    variants = propose.brics_variants(parent, n=4)

    assert variants
    assert len(variants) <= 4
    assert len({smiles for smiles, _ in variants}) == len(variants)
    assert all(_is_valid_smiles(smiles) for smiles, _ in variants)
    assert all(note == "BRICS recombination" for _, note in variants)


def test_propose_cli_all_strategy_emits_unique_lines():
    result = CliRunner().invoke(
        propose.cli,
        ["--smiles", "CC(=O)Oc1ccccc1", "--strategy", "all", "--n", "3", "--seed", "7"],
    )

    assert result.exit_code == 0, result.output
    lines = [ln for ln in result.output.splitlines() if ln.strip()]
    assert lines
    smis = [ln.split("|")[0] for ln in lines]
    assert len(smis) == len(set(smis))  # deduplicated across strategies
    for ln in lines:
        parts = ln.split("|")
        assert len(parts) == 3
        assert parts[1] == "CC(=O)Oc1ccccc1"  # parent column preserved


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


def test_plot_progress_and_pareto_render_with_scored_candidates(tmp_path):
    cfg = _base_cfg()
    cands = [
        _scored("CCO", 1.0, 1, mw=46.0),
        _scored("CCCCCCCC", 2.5, 1, mw=114.0, ok=False),  # violates constraints
        _scored("CO", 1.8, 2, mw=32.0),
    ]

    progress = tmp_path / "progress.png"
    pareto = tmp_path / "pareto.png"
    build_report.plot_progress(cfg, cands, progress)
    build_report.plot_pareto(cfg, cands, pareto)

    assert progress.exists() and progress.stat().st_size > 0
    assert pareto.exists() and pareto.stat().st_size > 0


def test_plot_grid_renders_top_k_structures(tmp_path):
    cfg = _base_cfg()
    cands = [
        _scored("CCO", 1.0, 1),
        _scored("c1ccccc1", 2.0, 1),
        _scored("CC(=O)O", 0.5, 2),
    ]

    out = tmp_path / "grid.png"
    build_report.plot_grid(cfg, cands, out, k=2)

    assert out.exists() and out.stat().st_size > 0


def test_plot_chemical_space_needs_three_points(tmp_path):
    cfg = _base_cfg()
    too_few = tmp_path / "few.png"
    enough = tmp_path / "enough.png"

    build_report.plot_chemical_space(cfg, [_scored("CCO", 1.0, 1)], too_few)
    build_report.plot_chemical_space(
        cfg,
        [_scored("CCO", 1.0, 1), _scored("c1ccccc1", 2.0, 1), _scored("CC(=O)O", 0.5, 2)],
        enough,
    )

    # Both render (the first as a placeholder), so just confirm files exist.
    assert too_few.exists() and too_few.stat().st_size > 0
    assert enough.exists() and enough.stat().st_size > 0


def test_plot_parity_renders_once_enough_data_exists(tmp_path):
    cfg = _base_cfg()
    smis = ["CCO", "CCCO", "CCCCO", "c1ccccc1", "c1ccccc1O", "CC(=O)O", "CCN"]
    cands = [_scored(s, float(i), 1) for i, s in enumerate(smis)]

    out = tmp_path / "parity.png"
    build_report.plot_parity(cfg, cands, out)

    assert len(cands) >= surrogate.MIN_TRAIN
    assert out.exists() and out.stat().st_size > 0


def test_best_sdf_returns_empty_for_missing_or_invalid_best():
    assert build_report.best_sdf(None) == ""
    assert build_report.best_sdf({"smiles": "not a smiles"}) == ""


def test_best_sdf_emits_a_3d_molblock_for_a_valid_molecule():
    block = build_report.best_sdf({"smiles": "CCO"})

    assert block  # non-empty molblock
    assert "V2000" in block  # standard SDF/MolBlock marker
    assert "3D" in block  # RDKit tags embedded conformers as 3D
    # Ethanol gains explicit hydrogens before embedding: C2H6O = 9 atoms.
    counts_line = block.splitlines()[3]
    assert int(counts_line[:3]) == 9


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
