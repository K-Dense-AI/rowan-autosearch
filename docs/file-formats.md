# File Formats

This page documents the files written under `runs/<run_id>/`.

## Run Directory

```text
runs/<run_id>/
├── config.json
├── iterations/
│   ├── 0001.json
│   ├── 0002.json
│   └── ...
├── plots/
│   ├── progress.png
│   ├── pareto.png
│   ├── grid.png
│   └── genealogy.png
└── report.html
```

The intended pattern is append-only. Add new iterations instead of editing old
ones, except for narrow cleanup such as removing an accidental duplicate or
fixing a decision field.

## `config.json`

`config.json` is the run contract. It tells the scorer what Rowan workflow to
call, how to extract or compute the score, what constraints to enforce, and how
many candidates the agent should propose.

Important fields:

```json
{
  "run_id": "aspirin_solubility",
  "objective": "Maximize aqueous solubility",
  "objective_direction": "maximize",
  "primary_metric": "logS",
  "workflow_type": "solubility",
  "start_smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "metric_path": null,
  "constraints": {
    "mw_max": 250
  },
  "workflow_params": {
    "method": "kingfisher",
    "solvents": ["water"],
    "temperatures": [298.15]
  },
  "max_iterations": 12,
  "candidates_per_iter": 4
}
```

See [Run configuration](run-configuration.md) for the full field guide.

## Iteration JSON

Each scoring call writes the next numbered file under `iterations/`.

Top-level fields:

```json
{
  "iter": 1,
  "created_at": "2026-04-26T00:00:00+00:00",
  "rationale": "Baseline score to confirm metric extraction.",
  "candidates": [],
  "iteration_best": null,
  "decision": "continue"
}
```

`iteration_best` is constraint-aware. If all candidates fail constraints or fail
to score, it is `null`.

## Candidate Records

Candidate records contain both local RDKit descriptors and Rowan-derived data:

```json
{
  "smiles": "CC(=O)Oc1ccccc1C(=O)O",
  "parent_smiles": null,
  "design_note": "baseline starting molecule",
  "scored_at": "2026-04-26T00:00:00+00:00",
  "metrics": {
    "mw": 180.159,
    "logp_crippen": 1.31,
    "tpsa": 63.6,
    "hbd": 1,
    "hba": 4,
    "rotb": 3,
    "heavy_atoms": 13,
    "rings": 1,
    "logS": -2.1
  },
  "satisfies_constraints": true,
  "constraint_failures": [],
  "score": -2.1,
  "workflow_uuid": "rowan-workflow-uuid",
  "metric_path_used": "solubilities.O.solubilities.0",
  "object_data": {}
}
```

When scoring fails, `score` is `null` and `error` explains why:

```json
{
  "smiles": "bad smiles",
  "error": "RDKit failed to parse SMILES",
  "score": null,
  "metrics": {},
  "satisfies_constraints": false,
  "constraint_failures": ["invalid SMILES"]
}
```

## Local Descriptor Keys

The scorer stores these RDKit descriptor keys in `metrics`:

```text
mw
logp_crippen
tpsa
hbd
hba
rotb
heavy_atoms
rings
```

The Rowan-derived primary metric is also stored in `metrics` when extraction
succeeds. Composite-objective scores add the objective name, usually
`objective_score`, and term details under `optimization_objective`.

## Report Files

`rowan-report` regenerates all report artifacts from `config.json` and
iteration JSON files:

- `progress.png`: all scored candidates by iteration plus best-so-far.
- `pareto.png`: score versus molecular weight.
- `grid.png`: 2D structures for top constraint-satisfying candidates.
- `genealogy.png`: parent-child edges across iterations.
- `report.html`: embedded summary, plots, candidate tables, rationales, and 3D
  view of the current winner when available.
