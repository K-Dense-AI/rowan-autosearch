# Scoring And Reporting

`rowan-score` appends scored candidates to a run. `rowan-report` turns the run
history into plots and an HTML report.

## Candidate Format

Each candidate is passed as one string:

```text
SMILES|PARENT_SMILES|DESIGN_NOTE
```

Example:

```bash
uv run rowan-score --run aspirin_solubility \
  --rationale "Test smaller polar aromatic substitutions for solubility gain while preserving acid functionality." \
  --candidate "COc1ccc(C(=O)O)c(OC)c1O|CC(=O)Oc1ccccc1C(=O)O|replace acetoxy with phenolic and add methoxy polarity" \
  --candidate "O=C(O)c1ccncc1O|CC(=O)Oc1ccccc1C(=O)O|pyridine ring swap to lower logP and add H-bond acceptor"
```

`PARENT_SMILES` may be empty for the baseline:

```bash
--candidate "CC(=O)Oc1ccccc1C(=O)O||baseline starting molecule"
```

## What Happens During Scoring

For each candidate, `rowan-score`:

1. Parses the SMILES with RDKit.
2. Computes local descriptors.
3. Checks local constraints.
4. Submits the configured Rowan workflow.
5. Waits for the Rowan result.
6. Saves the full Rowan `object_data`.
7. Extracts the primary metric or evaluates the composite objective.
8. Writes a new `iterations/NNNN.json` file.

For non-batch workflows, candidates are submitted concurrently with a
`ThreadPoolExecutor`. Use `--max-workers` to control concurrency:

```bash
uv run rowan-score --run <run_id> \
  --max-workers 2 \
  --rationale "..." \
  --candidate "..."
```

For `batch_docking`, the whole candidate list is submitted as one Rowan batch
workflow.

## Dry Runs

Use `--dry-run` to validate candidate formatting without calling Rowan:

```bash
uv run rowan-score --run <run_id> \
  --dry-run \
  --rationale "format check" \
  --candidate "<SMILES>|<PARENT>|<NOTE>"
```

Dry runs do not create iteration JSON files.

## Decisions

`--decision` records the agent's judgment for the iteration:

```bash
--decision continue
--decision converged
--decision stuck
```

The default is `continue`. Use `converged` or `stuck` on the scoring iteration
where the run outcome becomes clear.

You can also update the latest iteration's decision without rescoring or using
Rowan credits:

```bash
uv run rowan-score --run <run_id> \
  --rationale "The best score has plateaued after three strategies." \
  --decision stuck
```

If you omit `--candidate` and leave the decision as `continue`, `rowan-score`
adds a note-only wrap-up iteration with zero candidates and inherits the latest
iteration's decision.

## Iteration Output

A saved iteration looks like:

```json
{
  "iter": 2,
  "created_at": "2026-04-26T00:00:00+00:00",
  "rationale": "This iteration tests...",
  "candidates": [
    {
      "smiles": "O=C(O)c1ccccc1O",
      "parent_smiles": "CC(=O)Oc1ccccc1C(=O)O",
      "design_note": "deacetylate to phenolic salicylate",
      "metrics": {
        "mw": 138.122,
        "logp_crippen": 1.09,
        "tpsa": 57.53,
        "hbd": 2,
        "hba": 3,
        "rotb": 1,
        "heavy_atoms": 10,
        "rings": 1,
        "logS": -1.8
      },
      "satisfies_constraints": true,
      "constraint_failures": [],
      "score": -1.8,
      "workflow_uuid": "..."
    }
  ],
  "iteration_best": {
    "smiles": "O=C(O)c1ccccc1O",
    "score": -1.8
  },
  "decision": "continue"
}
```

The real file also includes `object_data`, which can be large. Keep it because
it is the audit trail for metric extraction and future analysis.

## Reports

Regenerate reports after scoring:

```bash
uv run rowan-report --run <run_id>
```

Use `--top-k` to control how many molecules appear in the 2D grid:

```bash
uv run rowan-report --run <run_id> --top-k 12
```

Report outputs:

```text
runs/<run_id>/plots/progress.png
runs/<run_id>/plots/pareto.png
runs/<run_id>/plots/grid.png
runs/<run_id>/plots/genealogy.png
runs/<run_id>/report.html
```

The report is constraint-aware: candidates that violate constraints remain
visible, but the best-so-far result only considers candidates with numeric
scores that satisfy constraints.
