# Active-Learning Advisor (`rowan-suggest`)

`rowan-suggest` is a quantitative co-pilot for the proposer. It does **not**
take over the search — the agent (or chemist) is still the optimizer. What it
does is stop you from spending Rowan credits blind: it learns a cheap surrogate
from the run's *own* scored history and uses it to dedup, rank, and diversify
your proposals before any Rowan workflow is submitted.

Everything it does is local. It never calls Rowan.

## Why it exists

Each scored candidate triggers a full Rowan quantum/ML workflow, so the
dominant cost of a run is **credits per unit of real improvement**. Three
common ways credits get wasted:

1. **Re-scoring a molecule already evaluated** (easy to do across many
   iterations).
2. **Scoring proposals that will not beat the current best** when a model
   trained on what you already know could have told you so cheaply.
3. **Collapsing into one chemotype**, so each new candidate is barely
   different from the last.

`rowan-suggest` addresses all three.

## What it computes

For each candidate SMILES you pass in:

- **Canonicalization + dedup.** Canonical SMILES are matched against every
  molecule already scored in the run. Duplicates are flagged (with their prior
  score and iteration) and excluded from recommendations unless you pass
  `--include-seen`.
- **Surrogate prediction.** A model trained on the run history predicts the
  objective value and an uncertainty (`pred=mu +/- sigma`).
- **Expected Improvement (EI).** The acquisition score that balances "predicted
  to be good" against "uncertain, worth learning about", relative to the
  current best and the run's optimization direction.
- **Novelty.** `1 - max Tanimoto similarity` to the evaluated set (1.0 = a
  brand-new region of chemical space).
- **Scaffold.** The Bemis-Murcko scaffold, so you can see chemotype at a glance.
- **Constraint check.** The same RDKit gates used by `rowan-score`; failing
  candidates are dropped from recommendations (unless
  `--no-respect-constraints`).

It then returns a **diverse top-k** (greedy MaxMin over the highest-ranked
pool) as paste-ready `SMILES|PARENT|NOTE` lines.

## The surrogate model

- **`rf` (default)** — a random forest on Morgan/ECFP4 fingerprints.
  Uncertainty is the spread of predictions across trees. Robust and fast.
- **`gp`** — a Gaussian process with a Tanimoto kernel over the fingerprints.
  Better-calibrated uncertainty, slightly slower; a good choice once you have a
  few dozen points.

### Trust it only as far as cross-validation allows

`rowan-suggest` prints a cross-validated `R^2` and MAE for the surrogate. A
high `R^2` means the EI ranking is meaningful; a low or negative `R^2` means
the model cannot yet predict your objective from structure, so treat the
ranking as a weak prior and lean on your own chemical reasoning. The same
information is plotted as the **parity** panel in `report.html`.

### Cold start

Below ~6 usable training points the surrogate is untrustworthy, so the tool
**skips modeling entirely and ranks by novelty**, saying so plainly in its
output. This pushes early iterations toward broad chemical-space coverage,
which is exactly what you want before a model can be fit.

## Usage

Screen the agent's own proposals:

```bash
uv run rowan-suggest --run <run_id> \
  --candidate "<SMILES>|<PARENT>|<NOTE>" \
  --candidate "<SMILES>|<PARENT>|<NOTE>" \
  --top-k 4
```

Generate *and* rank variants from a parent in one shot (uses `rowan-propose`
under the hood):

```bash
uv run rowan-suggest --run <run_id> --from-parent "<parent>" \
  --strategy all --n 20 --top-k 4
```

Structured output for programmatic use:

```bash
uv run rowan-suggest --run <run_id> --candidate "..." --json
```

### Options

| Option | Default | Notes |
| --- | --- | --- |
| `--candidate` | — | Repeatable `SMILES|PARENT|NOTE` proposals to rank. |
| `--from-parent` | — | Auto-generate variants from this parent via `rowan-propose`. |
| `--strategy` | `all` | `bioisostere`, `scan-subst`, `brics`, or `all` (for `--from-parent`). |
| `--n` | `20` | Max variants per strategy for `--from-parent`. |
| `--top-k` | `4` | Number of candidates to recommend. |
| `--xi` | `0.01` | EI exploration parameter (higher = more exploratory). |
| `--model` | `rf` | `rf` (random forest) or `gp` (Tanimoto-kernel GP). |
| `--diverse / --no-diverse` | on | MaxMin diversity selection over the top pool. |
| `--respect-constraints / --no-respect-constraints` | on | Drop candidates failing the run's RDKit gates. |
| `--include-seen` | off | Include molecules already evaluated in this run. |
| `--json` | off | Emit structured JSON instead of paste-ready lines. |

## Where it fits in the loop

```text
propose ideas  ->  rowan-suggest (dedup + rank + diversify, local)  ->  rowan-score (Rowan)  ->  rowan-report
```

It slots between proposing and scoring. See [`../AGENTS.md`](../AGENTS.md) for
the full driver protocol.

## Important caveats

- **Advisory only.** It ranks; it never submits to Rowan. Override it whenever
  your chemical reasoning disagrees — a surrogate trained on a handful of
  points is a hint, not an oracle.
- **It optimizes the run's recorded `score`.** For composite objectives, that
  is the combined scalar, so the surrogate learns the same objective the agent
  is optimizing.
- **Fingerprint-based.** Molecules with the same Morgan fingerprint look
  identical to the model; very subtle electronic effects may be invisible to
  the surrogate even when Rowan resolves them.
