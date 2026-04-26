# Examples

The `examples/` directory contains copy-pasteable run initializers. Use them as
templates for new searches and as regression examples for expected CLI shape.

## Available Examples

- `examples/aspirin_solubility.md`: maximize aqueous logS for aspirin-like
  acids.
- `examples/phenol_pka.md`: lower phenol pKa with electron-withdrawing
  substituents.
- `examples/quinone_redox.md`: tune benzoquinone toward a more positive
  reduction potential.
- `examples/ibuprofen_logp.md`: lower computed logP while preserving a simple
  arylpropionic-acid pharmacophore.
- `examples/aspirin_balanced_objective.md`: optimize a weighted objective that
  rewards solubility and lightly penalizes molecular weight.
- `examples/caffeine_tpsa.md`: push caffeine-like xanthines toward higher TPSA
  as a lower-CNS-penetration thought experiment.
- `examples/lidocaine_logd.md`: lower predicted logD for lidocaine-like local
  anesthetics.
- `examples/coumarin_dye_solubility.md`: water-solubilize a compact coumarin
  dye.
- `examples/anthraquinone_battery_redox.md`: tune anthraquinone-like organic
  battery scaffolds toward higher reduction potential.

## Choosing A Template

Use a solubility example when validating installation because the scorer has
explicit solubility defaults and built-in metric extraction paths.

Use a pKa or redox example when testing whether a chemical hypothesis changes an
electronic property in the intended direction. These workflows are more
sensitive to metric-path details, so score the baseline first.

Use the balanced objective example when the objective cannot be represented by
one Rowan property. Composite objectives are useful for trading off a Rowan
metric against local descriptors such as molecular weight, TPSA, or logP.

Use docking-family examples only when you already have the required Rowan
protein and pocket inputs. Docking runs should usually minimize
`docking_score`.

## Adapting An Example

1. Change `--run` to a new slug.
2. Rewrite `--objective` so the agent has a clear chemical goal.
3. Set `--start-smiles` to the scaffold or molecule you want to optimize.
4. Choose `--direction` based on whether larger or smaller values are better.
5. Confirm `--workflow`, `--metric`, and `--metric-path` match the Rowan output.
6. Add local constraints that reflect the project goal.
7. Score the starting molecule once before proposing variants.

## Example Agent Prompt

```text
Drive run phenol_pka per AGENTS.md. First run rowan-state status and inspect the
latest iteration. Propose four candidates close to the current best, each with a
clear parent SMILES and design note. Score them with rowan-score, rebuild the
report, and explain whether the run should continue.
```

## Baseline Check Pattern

For every adapted example, run a baseline check:

```bash
uv run rowan-score --run <run_id> \
  --rationale "Baseline score to confirm metric extraction before proposing variants." \
  --candidate "<START_SMILES>||baseline starting molecule"
```

If the score is missing, fix `metric_path` before continuing.
