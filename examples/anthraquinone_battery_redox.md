# Example: tune anthraquinone for organic batteries

This redox example uses anthraquinone as an organic battery-like scaffold and
searches for analogs with a more positive reduction potential.

```bash
uv run rowan-state init \
  --run anthraquinone_battery_redox \
  --objective "Maximize reduction potential for anthraquinone-like organic battery candidates while keeping MW < 360 and logP < 4" \
  --direction maximize \
  --metric redox_potential \
  --workflow redox_potential \
  --start-smiles "O=C1c2ccccc2C(=O)c2ccccc21" \
  --constraint mw_max=360 \
  --constraint logp_max=4 \
  --constraint rings_max=4 \
  --candidates-per-iter 4 \
  --max-iter 12
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `anthraquinone_battery_redox` per AGENTS.md. Score anthraquinone
> first to confirm the redox metric path, then search for more positive
> reduction potential analogs.

The agent should:

1. Keep the anthraquinone redox core intact in early iterations.
2. Probe electron-withdrawing substituents such as fluoro, cyano, carbonyl, and
   sulfone-like groups on the outer rings.
3. If the search stalls, test solubilizing polar handles that might also shift
   electronics, while respecting the MW and ring-count constraints.

Expected: electron-withdrawing anthraquinones should move the reduction
potential positive, while polar handles can trade off redox gain against
descriptor headroom.
