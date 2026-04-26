# Rowan autosearch — agent loop

You are driving an iterative property-optimization search over chemical space.
A starting molecule and an objective are pinned in `runs/<run_id>/config.json`.
Your job is to propose better molecules, score them with Rowan, and converge
on the best candidate that satisfies the constraints. **You are the optimizer.**

Do **not** inspect `docs/` unless the user specifically asks you to work on or
consult documentation.

## The loop

Each iteration:

1. **Read state.** Look at `runs/<run_id>/config.json` and the most recent
   `runs/<run_id>/iterations/*.json`. Run:
   ```
   uv run rowan-state status --run <run_id>
   ```
   Note the best-so-far score, what's been tried, and what constraints are
   active.
2. **Form a hypothesis.** What chemical change should move the metric in the
   right direction *and* respect constraints? Examples: "add polar group to
   raise solubility", "remove H-bond donor to lower TPSA", "introduce
   electron-withdrawing group on the aryl ring to shift redox potential".
3. **Propose N candidates.** Default `candidates_per_iter` is in config
   (typically 4). Each candidate must have a SMILES, a parent SMILES (the
   molecule it's derived from — usually the current best or last iteration's
   pick), and a one-line *design note* explaining the change.
4. **Score them.**
   ```
   uv run rowan-score --run <run_id> \
     --rationale "What ties these candidates together — your hypothesis." \
     --candidate "<SMILES>|<PARENT>|<NOTE>" \
     --candidate "<SMILES>|<PARENT>|<NOTE>" \
     ...
   ```
   This submits to Rowan in parallel, waits for results, parses the metric,
   checks constraints locally via RDKit, and writes the iteration JSON.
5. **Rebuild the report.**
   ```
   uv run rowan-report --run <run_id>
   ```
   Open `runs/<run_id>/report.html` to see progress visually.
6. **Decide.** Did the metric improve? Are you stuck? If converged or stuck,
   re-run `rowan-score` with `--decision converged` or `--decision stuck` (no
   new candidates needed — just amend the latest iteration's decision).

## Inspect before you commit

If you've never run this workflow before, **always score one candidate first**
to confirm `metric_path` is correct. The candidate's full `object_data` is
saved in the iteration JSON — read it and, if the metric extractor missed the
right path, edit `runs/<run_id>/config.json` to set `metric_path` to the
correct dot-path (e.g. `"solubility.logS_25C"`).

A failed metric extraction shows up as `score: null` plus an `error` field
pointing to the raw payload. Fix the path and rescore — don't keep guessing.

## Proposing good candidates

You're a medicinal chemist now. Some heuristics:

- **Stay near the current best** unless you're explicitly exploring. Most
  iterations should pick the current best as parent and make one small change.
- **Diversify within an iteration.** If you're submitting 4 candidates, don't
  make all 4 the same kind of change. Pair a polar swap with a steric one,
  for example.
- **Think about the metric mechanistically.** Aqueous solubility favors
  H-bond donors/acceptors and lower logP. pKa shifts with electron-withdrawing
  vs. donating groups near the ionizable site. Redox potentials shift with
  conjugation extent.
- **Respect constraints early.** RDKit-computed MW/logP/TPSA are in the
  iteration JSON; don't bother proposing a candidate that already busts MW.
- **Use `rowan-propose` only when stuck.** It's a brainstorming tool that
  applies bioisostere swaps, aromatic substitutions, or BRICS recombination
  to a parent. Treat its output as inspiration — pick the chemically sensible
  ones, drop the rest.
  ```
  uv run rowan-propose --smiles "<parent>" --strategy bioisostere --n 6
  ```

## Writing rationales

The `--rationale` is saved per iteration and shown in the report. Write it
like a lab notebook entry: **what you observed last iteration → what you
hypothesize → what these candidates test**. One paragraph. Future-you (or a
reviewer) should be able to read the report top-to-bottom and follow your
reasoning without seeing the conversation.

Bad: "Trying some variants."
Good: "Iter 3 showed that adding -OMe at the para position raised logS by
0.4. This iter probes whether two -OMe groups (para + meta) compound that
gain or just add MW without payoff, plus a sulfonamide swap as a more
aggressive polarity move."

## Convergence

Stop the loop when **any** of:

- Best score has not improved over the last 10 iterations and you've tried at
  least 5 distinct chemical strategies in that window.
- You've hit `max_iterations` from config.
- You've converged on a clearly winning candidate and further variants only
  trade away constraint headroom.

Mark the final iteration with `--decision converged`. Then write a short
summary at the end of the report by adding a final "wrap-up" iteration with
zero candidates: just `--rationale` describing what you found.

## Common pitfalls

- **Invalid SMILES** show up as `score: null` with `error: RDKit failed to
  parse SMILES`. Check parens and ring closures.
- **All candidates fail constraints.** Loosen the chemical move or check that
  the constraint values in config.json match the project goal.
- **Rowan rate limits / failures.** Re-run `rowan-score` with the same
  candidates; failed ones will resubmit. (Successful ones will create
  duplicates in a new iteration — better to manually drop the failed
  candidates from the prior iteration JSON if you re-score.)
- **Don't mutate prior iteration JSONs** except to fix `decision` or to drop
  a duplicate. Append-only keeps the genealogy plot honest.
