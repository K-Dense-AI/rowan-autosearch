# Run Configuration

Every run is controlled by `runs/<run_id>/config.json`. The file is created by
`rowan-state init` and then read by `rowan-score`, `rowan-state status`, and
`rowan-report`.

## Minimal Shape

```json
{
  "run_id": "aspirin_solubility",
  "objective": "Maximize aqueous solubility while keeping MW < 250",
  "objective_direction": "maximize",
  "primary_metric": "logS",
  "workflow_type": "solubility",
  "start_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "metric_path": null,
  "constraints": {
    "mw_max": 250
  },
  "max_iterations": 12,
  "candidates_per_iter": 4,
  "created_at": "2026-04-26T00:00:00+00:00"
}
```

## Core Fields

- `run_id`: slug used for `runs/<run_id>/`.
- `objective`: plain-English goal. The agent uses this when proposing
  molecules.
- `objective_direction`: `maximize` or `minimize`.
- `primary_metric`: metric name stored in candidate `metrics` and optimized as
  `score` unless a composite objective is configured.
- `workflow_type`: Rowan workflow suffix, such as `solubility`, `pka`,
  `redox_potential`, `descriptors`, or `batch_docking`.
- `start_smiles`: root molecule for the search.
- `metric_path`: optional dot-path into Rowan `object_data`. Leave as `null`
  only when the built-in metric extractor knows the workflow and metric.
- `max_iterations`: search budget for the agent.
- `candidates_per_iter`: target batch size for the agent.

## Workflow Parameters

`workflow_params` are passed to Rowan submitter functions. Values passed through
the CLI are parsed as JSON when possible:

```bash
uv run rowan-state init \
  --run aspirin_solubility \
  --objective "Maximize aqueous logS" \
  --direction maximize \
  --metric logS \
  --workflow solubility \
  --start-smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --workflow-param method=kingfisher \
  --workflow-param solvents='["water"]' \
  --workflow-param temperatures='[298.15]'
```

This produces:

```json
{
  "workflow_params": {
    "method": "kingfisher",
    "solvents": ["water"],
    "temperatures": [298.15]
  }
}
```

For `solubility`, the scorer defaults to `method=kingfisher`,
`solvents=["water"]`, and `temperatures=[298.15]` unless overridden.

## Constraints

Constraints are checked locally with RDKit before Rowan results are interpreted.
Candidates that fail constraints are still scored and recorded, but they do not
count toward constraint-aware best-so-far.

Supported keys:

```text
mw_max      molecular weight maximum
mw_min      molecular weight minimum
logp_max    Crippen logP maximum
logp_min    Crippen logP minimum
tpsa_max    topological polar surface area maximum
tpsa_min    topological polar surface area minimum
hbd_max     H-bond donor maximum
hba_max     H-bond acceptor maximum
rotb_max    rotatable bond maximum
ha_max      heavy atom maximum
ha_min      heavy atom minimum
rings_max   ring count maximum
```

CLI example:

```bash
--constraint mw_max=350 --constraint logp_max=4 --constraint tpsa_max=140
```

## Composite Objectives

A composite objective replaces the primary extracted metric as `score`. Each
term can come from Rowan `object_data` or local RDKit descriptors.

```bash
uv run rowan-state init \
  --run aspirin_balanced \
  --objective "Maximize aqueous logS while mildly penalizing molecular weight" \
  --direction maximize \
  --metric objective_score \
  --workflow solubility \
  --start-smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --workflow-param method=kingfisher \
  --workflow-param solvents='["water"]' \
  --workflow-param temperatures='[298.15]' \
  --objective-term '{"name":"logS","path":"solubilities.O.solubilities.0","goal":"maximize","weight":1.0}' \
  --objective-term '{"name":"mw","source":"local","metric":"mw","goal":"minimize","weight":0.01}'
```

Equivalent config excerpt:

```json
{
  "optimization_objective": {
    "name": "objective_score",
    "terms": [
      {
        "name": "logS",
        "path": "solubilities.O.solubilities.0",
        "goal": "maximize",
        "weight": 1.0
      },
      {
        "name": "mw",
        "source": "local",
        "metric": "mw",
        "goal": "minimize",
        "weight": 0.01
      }
    ]
  }
}
```

Supported term goals:

- `maximize`: contribution is `weight * value`.
- `minimize`: contribution is `-weight * value`.
- `target`: contribution is `-weight * abs(value - target)`.

Local metric names include `mw`, `logp_crippen`, `tpsa`, `hbd`, `hba`, `rotb`,
`heavy_atoms`, and `rings`.
