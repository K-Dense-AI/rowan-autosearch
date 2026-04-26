# Example: maximize aspirin's aqueous solubility

```bash
uv run rowan-state init \
  --run aspirin_solubility \
  --objective "Maximize aqueous solubility (logS) starting from aspirin while keeping MW < 250 and preserving a carboxylic acid for activity" \
  --direction maximize \
  --metric logS \
  --workflow solubility \
  --start-smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --constraint mw_max=250 \
  --candidates-per-iter 4 \
  --max-iter 12
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `aspirin_solubility` per AGENTS.md. Score one baseline candidate
> first to confirm metric_path, then proceed with the loop.

The agent should:

1. Score the start SMILES alone first, confirm `logS` extracts cleanly.
2. Iter 2: try 4 polarity-increasing variants (replace acetyl with -OH,
   add hydroxyl groups, swap to a sulfonate ester, etc.).
3. Iter 3+: pick the best, branch from it, narrow toward the optimum.

Expected: the agent should land on something more soluble than aspirin
(logS > -2.4) while keeping the carboxylic acid intact.

## More examples

See `examples/README.md` for additional pKa, redox, descriptor, and composite
objective walkthroughs.
