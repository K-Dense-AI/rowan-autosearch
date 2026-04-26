# Example: maximize quinone reduction potential

This example asks the agent to tune a benzoquinone scaffold toward a more
positive reduction potential.

```bash
uv run rowan-state init \
  --run quinone_redox \
  --objective "Maximize the reduction potential of a benzoquinone-like scaffold while keeping MW < 300 and logP < 4" \
  --direction maximize \
  --metric redox_potential \
  --workflow redox_potential \
  --start-smiles "O=C1C=CC(=O)C=C1" \
  --constraint mw_max=300 \
  --constraint logp_max=4 \
  --candidates-per-iter 4 \
  --max-iter 12
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `quinone_redox` per AGENTS.md. Score benzoquinone first to confirm
> the redox metric path, then search for higher reduction potential variants.

The agent should:

1. Score the parent quinone before proposing variants.
2. Probe electron-withdrawing groups on the ring, such as fluoro, chloro,
   cyano, and trifluoromethyl substituents.
3. Keep changes close to the quinone core so the search tests electronics
   rather than scaffold replacement.

Expected: electron-withdrawing substitution should make reduction easier and
push the reduction potential more positive, subject to the local constraints.
