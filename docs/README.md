# Rowan Autosearch Documentation

This directory is the working handbook for `rowan-autosearch`: how to set up a
run, configure Rowan workflows, drive the agent loop, read outputs, and recover
from common failures.

## Start Here

- [Getting started](getting-started.md): install the package, configure
  credentials, initialize a first run, and score the baseline molecule.
- [Optimization loop](optimization-loop.md): the iteration cycle an AI agent or
  chemist follows while proposing and scoring molecules.
- [Run configuration](run-configuration.md): every important field in
  `runs/<run_id>/config.json`, including constraints and composite objectives.
- [Workflows and metrics](workflows-and-metrics.md): supported Rowan workflow
  names, common metric aliases, and how metric extraction works.
- [Scoring and reporting](scoring-and-reporting.md): candidate format, Rowan
  submission behavior, iteration JSON, plots, and `report.html`.
- [File formats](file-formats.md): concrete schemas for run directories,
  config files, and iteration records.
- [Examples](examples.md): recommended examples and when to use each one.
- [Troubleshooting](troubleshooting.md): fixes for API keys, invalid SMILES,
  metric-path misses, constraint failures, and Rowan submission errors.

## Mental Model

`rowan-autosearch` is not an autonomous optimizer by itself. It is a scoring and
bookkeeping harness for an external optimizer: an AI coding agent, chemist, or
future search algorithm proposes molecules; this project scores them with Rowan,
checks local RDKit constraints, and produces an auditable report.

Each run lives under `runs/<run_id>/` and follows an append-only pattern:

```text
runs/<run_id>/
├── config.json
├── iterations/
│   ├── 0001.json
│   ├── 0002.json
│   └── ...
├── plots/
└── report.html
```

The three most important commands are:

```bash
uv run rowan-state status --run <run_id>
uv run rowan-score --run <run_id> --rationale "..." --candidate "SMILES|PARENT|NOTE"
uv run rowan-report --run <run_id>
```

See [Getting started](getting-started.md) for a complete first run.
