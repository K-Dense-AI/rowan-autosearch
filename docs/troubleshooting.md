# Troubleshooting

This page covers the failures most likely to occur during a Rowan autosearch
run.

## `ROWAN_API_KEY not set`

`rowan-score` requires a Rowan API key.

Fix:

```bash
cp .env.example .env
```

Then set `ROWAN_API_KEY` in `.env`. The scorer loads `.env` from the repository
root.

## Invalid SMILES

Symptom in iteration JSON:

```json
{
  "score": null,
  "error": "RDKit failed to parse SMILES",
  "constraint_failures": ["invalid SMILES"]
}
```

Fix:

- Check parentheses, ring closure numbers, charges, and aromatic atom casing.
- Prefer canonical SMILES from RDKit or a trusted drawing tool.
- Use `rowan-score --dry-run` for formatting checks before spending Rowan
  credits.

## Metric Extraction Fails

Symptom:

```json
{
  "score": null,
  "error": "Could not extract metric ..."
}
```

Fix:

1. Open the latest `runs/<run_id>/iterations/NNNN.json`.
2. Inspect `candidates[*].object_data`.
3. Find the numeric value that matches the objective.
4. Add its dot-path to `runs/<run_id>/config.json` as `metric_path`.
5. Rescore the candidate or score a fresh baseline iteration.

Example:

```json
{
  "metric_path": "solubilities.O.solubilities.0"
}
```

## All Candidates Fail Constraints

Symptoms:

- Candidate `score` values may be numeric.
- `satisfies_constraints` is `false`.
- `constraint_failures` lists the violated local descriptor limits.
- `rowan-state status` does not treat them as best-so-far.

Fix:

- Propose smaller or less lipophilic modifications.
- Check whether the configured constraints reflect the project goal.
- Consider easing constraints only if they were arbitrary or exploratory.
- Do not chase a high score that is outside hard project constraints.

## Rowan Submission Fails

Symptom:

```json
{
  "score": null,
  "error": "Rowan submit failed: ..."
}
```

Fix:

- Confirm `ROWAN_API_KEY` is valid.
- Confirm `workflow_type` maps to an available Rowan SDK submitter.
- Confirm required `workflow_params` are present.
- Reduce `--max-workers` if failures look like rate limiting.
- Retry only failed candidates when possible to avoid duplicate successful
  records.

## Wrong Workflow Parameters

`--workflow-param` values are parsed as JSON when possible. Lists and objects
must be quoted for the shell:

```bash
--workflow-param solvents='["water"]'
--workflow-param temperatures='[298.15]'
--workflow-param pocket='[[0,0,0],[12,12,12]]'
```

Unquoted values are usually parsed as strings or numbers:

```bash
--workflow-param method=kingfisher
--workflow-param exhaustiveness=8
```

If Rowan rejects a parameter, compare `workflow_params` in `config.json` with
the current Rowan SDK submitter signature.

## Report Does Not Show A Winner

The report uses constraint-aware best-so-far logic. It will not show a winner if
no candidate has both:

- a numeric `score`
- `satisfies_constraints: true`

Fix scoring errors or constraints first, then rerun:

```bash
uv run rowan-report --run <run_id>
```

## `rowan-propose` Output Is Chemically Odd

`rowan-propose` is a brainstorming helper. It uses simple RDKit transformations
and may emit candidates that are syntactically valid but strategically weak.

Fix:

- Treat output as suggestions, not an auto-submit list.
- Filter for chemistry that matches the objective.
- Check constraints locally through `rowan-score --dry-run` or by estimating
  descriptor changes before scoring.

## Duplicate Or Bad Iteration Files

The run history is intended to be append-only. Avoid rewriting old iterations.

Acceptable cleanup:

- Remove an accidental duplicate candidate before further analysis.
- Fix a mistaken `decision` value.
- Drop failed duplicates created by a retry, if the successful result already
  exists elsewhere.

When in doubt, append a new iteration and explain the correction in the
`--rationale`.
