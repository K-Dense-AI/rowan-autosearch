# Example: lower phenol pKa

Use this when you want the agent to explore substituent effects around an acidic
phenol.

```bash
uv run rowan-state init \
  --run phenol_pka \
  --objective "Lower the pKa of phenol below 7 by adding electron-withdrawing substituents while keeping MW < 250" \
  --direction minimize \
  --metric pka \
  --workflow pka \
  --start-smiles "Oc1ccccc1" \
  --constraint mw_max=250 \
  --constraint logp_max=4 \
  --candidates-per-iter 4 \
  --max-iter 10
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `phenol_pka` per AGENTS.md. Score the unsubstituted phenol first to
> confirm pKa extraction, then optimize toward lower pKa.

The agent should:

1. Score phenol first and verify that `pka` maps to the acidic phenol value.
2. Test para and meta electron-withdrawing groups such as nitro, cyano,
   trifluoromethyl, and fluoro substituents.
3. If one substituent class works, branch into combinations or positional
   variants while monitoring MW and logP.

Expected: nitro or cyano substitution should lower the phenolic pKa relative to
phenol while staying under the MW and logP limits.
