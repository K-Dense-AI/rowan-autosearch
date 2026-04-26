# Optimization Loop

The optimizer is external to this package. In normal use, an AI coding agent or
chemist reads the latest run state, proposes molecular variants, scores them
with Rowan, and decides what to try next.

## One Iteration

1. Read the run state.

   ```bash
   uv run rowan-state status --run <run_id>
   ```

   Also inspect `runs/<run_id>/config.json` and the latest
   `runs/<run_id>/iterations/NNNN.json`.

2. Form a chemical hypothesis.

   Use the current best molecule, objective direction, constraints, and recent
   failures. Examples:

   - Add polar functionality to improve aqueous solubility.
   - Add an electron-withdrawing group near an acidic site to lower pKa.
   - Reduce lipophilic substituent size to lower logP.
   - Extend conjugation or tune substituents to shift redox potential.

3. Propose candidates.

   Each candidate is a small record serialized as:

   ```text
   <SMILES>|<PARENT_SMILES>|<DESIGN_NOTE>
   ```

   Prefer small moves near the current best unless the run is explicitly
   exploring a new scaffold. Diversify the candidates in one iteration so a
   single Rowan batch tests several plausible ideas.

4. Score candidates.

   ```bash
   uv run rowan-score --run <run_id> \
     --rationale "What was observed, what is hypothesized, and what these candidates test." \
     --candidate "<SMILES>|<PARENT>|<NOTE>" \
     --candidate "<SMILES>|<PARENT>|<NOTE>"
   ```

5. Rebuild the report.

   ```bash
   uv run rowan-report --run <run_id>
   ```

6. Decide whether to continue.

   Continue when at least one candidate improved the score or the result points
   to a useful next hypothesis. Stop when the run has converged, hit its
   iteration budget, or stopped improving after multiple distinct strategies.

## Candidate Design Guidance

Good iterations are hypothesis-driven, not random enumeration.

- Keep the parent close to the current best for most iterations.
- Make one interpretable chemical change per candidate where possible.
- Include at least one conservative candidate and one more exploratory
  candidate when the budget allows.
- Respect local constraints before spending Rowan credits.
- Use `rowan-propose` only as a brainstorming helper, then filter its output.

Example:

```bash
uv run rowan-propose --smiles "<parent_smiles>" --strategy bioisostere --n 6
```

## Rationale Style

The `--rationale` is saved in the iteration JSON and displayed in the report.
Write it as a lab notebook paragraph:

```text
Iter 2 showed that para-methoxy improved logS but consumed MW headroom. This
iteration tests whether smaller polar substituents retain the solubility gain
while staying below MW 250, plus one pyridine ring swap to reduce logP.
```

Avoid generic rationales such as "trying variants".

## Stopping Criteria

Stop the run when any of these is true:

- Best score has not improved over the last five iterations and at least two
  distinct strategies were tested.
- The run has reached `max_iterations`.
- The current winner is clearly best and nearby variants only trade away
  constraint headroom.
- Rowan failures or metric-path ambiguity prevent meaningful comparison until
  configuration is fixed.

Record `--decision converged` or `--decision stuck` when the outcome is clear:

```bash
uv run rowan-score --run <run_id> \
  --rationale "Why the latest result is final." \
  --decision converged
```

With no `--candidate`, `rowan-score` updates the latest iteration's decision
without calling Rowan. A no-candidate run with the default `continue` decision
appends a note-only wrap-up iteration for final summary text and inherits the
latest iteration's decision.
