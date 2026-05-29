# Example: maximize aspirin's aqueous solubility

```bash
uv run rowan-state init \
  --run aspirin_solubility \
  --objective "Maximize aqueous solubility (logS) starting from aspirin while keeping MW < 250 and preserving a carboxylic acid for activity" \
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

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `aspirin_solubility` per AGENTS.md. Score one baseline candidate
> first to confirm metric_path, then proceed with the loop.

The agent should:

1. Score the start SMILES alone first, confirm `logS` extracts cleanly.
2. Iter 2: try 4 polarity-increasing variants (hydrolyze the acetyl ester to
   the free phenol, add hydroxyl or carboxamide groups, or introduce a small
   polar substituent on the ring).
3. Iter 3+: pick the best, branch from it, narrow toward the optimum.

Expected: the agent should land on a candidate with higher predicted `logS`
than the aspirin baseline scored in step 1, while keeping the carboxylic acid
intact. (Aspirin's measured aqueous solubility is modest — logS ≈ -1.7 to -2.0
— so even a small, polarity-increasing change should register.)

## More examples

See `examples/README.md` for additional pKa, redox, descriptor, and composite
objective walkthroughs.
