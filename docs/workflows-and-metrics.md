# Workflows And Metrics

`rowan-score` submits Rowan workflows using the configured `workflow_type`, then
extracts a scalar metric from the returned `object_data`. That scalar becomes
the candidate `score` unless a composite objective is configured.

## Workflow Names

For most workflows, `rowan-score` calls:

```text
rowan.submit_<workflow_type>_workflow(...)
```

The workflow type is the suffix, for example `solubility` maps to
`rowan.submit_solubility_workflow`.

Common SMILES-first workflow types:

```text
solubility
pka
macropka
redox_potential
admet
bde
descriptors
conformer_search
tautomer_search
electronic_properties
fukui
hydrogen_bond_basicity
hydrogen_bond_donor_acceptor_strength
ion_mobility
irc
membrane_permeability
multistage_optimization
nmr
scan
solvent_dependent_conformers
strain
```

Some Rowan workflows require extra non-SMILES inputs such as proteins, pockets,
ligand sets, reaction endpoints, or spin states. Provide those through
`workflow_params` when the Rowan SDK submitter requires them.

## Built-In Metric Paths

When `metric_path` is not set, `rowan-score` tries workflow-specific paths in
order and uses the first numeric value.

| Workflow | Metric | Paths Tried |
| --- | --- | --- |
| `solubility` | `logS` | `solubilities.O.solubilities.0`, `solubility.logS`, `logS` |
| `solubility` | `logS_25C` | `solubility.logS_25C` |
| `pka` | `pka` | `strongest_acid`, `pka_values.0`, `pka` |
| `redox_potential` | `redox_potential` | `reduction_potential`, `oxidation_potential`, `redox_potential` |
| `admet` | `logP` | `properties.logP`, `logP` |
| `admet` | `logD` | `properties.logD`, `logD` |
| `admet` | `solubility` | `properties.solubility`, `solubility` |
| `bde` | `bde` | `bde`, `weakest_bond.bde` |
| `descriptors` | `logP` | `descriptors.logP` |
| `descriptors` | `tpsa` | `descriptors.tpsa` |
| `docking` | `docking_score` | `scores.0.score`, `scores.0` |
| `batch_docking` | `docking_score` | `scores` |

After workflow-specific paths, the scorer also checks a top-level key matching
`primary_metric`.

## Setting `metric_path`

Use `--metric-path` when the intended metric is not covered by the built-in
extractor:

```bash
uv run rowan-state init \
  --run custom_solubility \
  --objective "Maximize aqueous logS" \
  --direction maximize \
  --metric logS \
  --workflow solubility \
  --metric-path solubilities.O.solubilities.0 \
  --start-smiles "CC(=O)Oc1ccccc1C(=O)O"
```

Dot-path segments can traverse dictionaries and lists. For example,
`pka_values.0` means:

```text
object_data["pka_values"][0]
```

## Discovering A New Metric Path

For a new workflow, score one baseline molecule:

```bash
uv run rowan-score --run <run_id> \
  --rationale "Baseline metric-path check." \
  --candidate "<START_SMILES>||baseline"
```

Then inspect the latest `runs/<run_id>/iterations/NNNN.json`:

```text
candidates[*].object_data
```

Find the value you want to optimize and update `runs/<run_id>/config.json`:

```json
{
  "metric_path": "path.to.numeric.value"
}
```

Rescore candidates after fixing the path so their `score` fields are comparable.

## Batch Docking

`batch_docking` is special: the scorer submits one Rowan batch workflow for all
candidates in the iteration.

Example config excerpt:

```json
{
  "workflow_type": "batch_docking",
  "primary_metric": "docking_score",
  "objective_direction": "minimize",
  "workflow_params": {
    "protein": "protein-uuid-or-object-id",
    "pocket": [[0.0, 0.0, 0.0], [12.0, 12.0, 12.0]],
    "executable": "qvina2",
    "scoring_function": "vina",
    "exhaustiveness": 8
  }
}
```

Because docking scores are usually better when lower, initialize docking runs
with `--direction minimize`.
