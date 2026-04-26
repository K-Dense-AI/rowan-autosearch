# Getting Started

This guide walks from a clean checkout to a first scored Rowan iteration.

## Requirements

- Python 3.12 or newer.
- `uv` is recommended for environment management.
- A Rowan API key from `labs.rowansci.com/account/api-keys`.

## Install

Create an environment and install the project in editable mode:

```bash
uv venv
uv pip install -e .
```

The package exposes these console scripts:

```text
rowan-state
rowan-score
rowan-propose
rowan-report
```

## Configure Rowan Credentials

Copy the example environment file and add your API key:

```bash
cp .env.example .env
```

Set:

```bash
ROWAN_API_KEY=<your Rowan API key>
```

Do not commit `.env`.

## Initialize A Run

This example starts an aspirin-like aqueous solubility search:

```bash
uv run rowan-state init \
  --run aspirin_solubility \
  --objective "Maximize aqueous solubility (logS) while keeping MW < 250 and preserving the carboxylic acid" \
  --direction maximize \
  --metric logS \
  --workflow solubility \
  --start-smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --workflow-param method=kingfisher \
  --workflow-param solvents='["water"]' \
  --workflow-param temperatures='[298.15]' \
  --constraint mw_max=250 \
  --candidates-per-iter 4 \
  --max-iter 12
```

This creates:

```text
runs/aspirin_solubility/
├── config.json
├── iterations/
└── plots/
```

## Check State

Before scoring, inspect the run:

```bash
uv run rowan-state status --run aspirin_solubility
```

Use this before every iteration. It shows the objective, workflow, iteration
count, and constraint-aware best-so-far molecule.

## Score The Baseline First

For a new workflow or metric, score one known molecule first. This confirms that
`rowan-score` can extract the intended metric from Rowan `object_data`.

```bash
uv run rowan-score --run aspirin_solubility \
  --rationale "Baseline aspirin score to confirm the solubility metric path before proposing variants." \
  --candidate "CC(=O)Oc1ccccc1C(=O)O||baseline starting molecule"
```

If the resulting candidate has `score: null`, inspect the latest
`runs/aspirin_solubility/iterations/NNNN.json`, find the correct value inside
`candidates[*].object_data`, and set `metric_path` in `config.json`. See
[Workflows and metrics](workflows-and-metrics.md).

## Generate The Report

After every scoring iteration, rebuild the report:

```bash
uv run rowan-report --run aspirin_solubility
open runs/aspirin_solubility/report.html
```

The report contains progress plots, a Pareto plot, a top-candidate structure
grid, a genealogy diagram, scored-candidate tables, and a 3D view of the current
winner when RDKit can embed it.

## Hand Off To An Agent

Once the baseline is verified, ask an agent to drive the loop:

```text
Drive run aspirin_solubility per AGENTS.md. Use rowan-state status, propose
four chemically sensible candidates per iteration, score them, rebuild the
report, and stop when the run has converged or is stuck.
```

The agent should keep each iteration rationale specific enough that the report
reads like a lab notebook.
