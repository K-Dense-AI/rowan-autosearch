# Example: make caffeine less CNS-like by raising TPSA

This descriptor example starts from caffeine and asks the agent to make close
analogs with higher topological polar surface area. It is a deliberately
creative prompt for exploring how far a familiar stimulant can be pushed toward
lower passive CNS penetration while staying small.

```bash
uv run rowan-state init \
  --run caffeine_tpsa \
  --objective "Maximize TPSA for caffeine-like xanthines while keeping MW < 260, logP < 2, and the xanthine core recognizable" \
  --direction maximize \
  --metric tpsa \
  --workflow descriptors \
  --start-smiles "Cn1c(=O)c2c(ncn2C)n(C)c1=O" \
  --constraint mw_max=260 \
  --constraint logp_max=2 \
  --constraint hbd_max=3 \
  --constraint hba_max=8 \
  --candidates-per-iter 4 \
  --max-iter 10
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `caffeine_tpsa` per AGENTS.md. Score caffeine first to confirm the
> descriptor `tpsa` path, then search for higher-TPSA xanthine analogs.

The agent should:

1. Keep the xanthine ring system intact for the first few iterations.
2. Replace one N-methyl group at a time with polar side chains or hydroxyalkyl
   groups.
3. Watch MW and H-bonding constraints so the search does not simply append a
   bulky polar tail.

Expected: hydroxyalkyl or small amide-bearing analogs should raise TPSA while
preserving the core and staying under the local descriptor limits.
