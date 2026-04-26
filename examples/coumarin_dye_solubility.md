# Example: water-solubilize a coumarin dye

This example starts from coumarin and asks the agent to improve aqueous
solubility while keeping the compact fluorescent dye scaffold recognizable.

```bash
uv run rowan-state init \
  --run coumarin_dye_solubility \
  --objective "Maximize aqueous solubility of a coumarin-like dye while keeping MW < 300 and preserving the lactone chromophore" \
  --direction maximize \
  --metric logS \
  --workflow solubility \
  --start-smiles "O=c1oc2ccccc2cc1" \
  --workflow-param method=kingfisher \
  --workflow-param solvents='["water"]' \
  --workflow-param temperatures='[298.15]' \
  --constraint mw_max=300 \
  --constraint logp_max=3 \
  --constraint tpsa_max=120 \
  --candidates-per-iter 4 \
  --max-iter 12
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `coumarin_dye_solubility` per AGENTS.md. Score coumarin first to
> confirm aqueous logS extraction, then optimize solubility without replacing
> the dye core.

The agent should:

1. Treat the coumarin lactone and fused aromatic system as the parent scaffold.
2. Try hydroxyl, methoxy, amino, and small sulfonamide-like substitutions before
   larger side chains.
3. Balance polarity against the TPSA ceiling so candidates remain plausible
   compact dye analogs.

Expected: phenolic or amino coumarins should improve logS; very polar variants
may score well but should be checked against MW, logP, and TPSA.
