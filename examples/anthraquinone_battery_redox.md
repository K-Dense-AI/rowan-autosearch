# Example: tune anthraquinone reduction potential

This redox example uses anthraquinone, the canonical scaffold for aqueous
organic redox-flow batteries, and searches for analogs with a more positive
reduction potential by adding electron-withdrawing groups.

**Mind the sign convention.** Anthraquinones are most often used as the
*negative* electrolyte (negolyte) in a flow cell, where a *more negative*
reduction potential raises cell voltage — that regime is reached with
electron-*donating* groups (–OH, –NH₂, –OMe), as in the 2,6-dihydroxy-
anthraquinone chemistries reported for high-voltage alkaline flow batteries.
This example deliberately tunes the *other* direction (more positive,
electron-withdrawing groups), which is the regime relevant to oxidant /
positive-electrode behavior, so the redox direction stays explicit. To design
a higher-voltage negolyte instead, flip `--direction` to `minimize` and probe
electron-donating substituents.

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

Expected: electron-withdrawing anthraquinones (e.g. fluoro, cyano, sulfone)
should move the reduction potential positive, while polar handles can trade off
redox gain against descriptor headroom. If you instead run the negolyte design
(`--direction minimize`), electron-donating hydroxy/amino analogs should move
the potential negative — the opposite, equally valid, optimization.
