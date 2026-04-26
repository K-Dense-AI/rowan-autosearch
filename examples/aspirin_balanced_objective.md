# Example: balance aspirin solubility and molecular weight

This example demonstrates a weighted objective. The primary score rewards
aqueous logS and lightly penalizes local RDKit molecular weight.

```bash
uv run rowan-state init \
  --run aspirin_balanced \
  --objective "Maximize aqueous solubility for aspirin-like acids while mildly penalizing molecular weight" \
  --direction maximize \
  --metric objective_score \
  --workflow solubility \
  --start-smiles "CC(=O)Oc1ccccc1C(=O)O" \
  --workflow-param method=kingfisher \
  --workflow-param solvents='["water"]' \
  --workflow-param temperatures='[298.15]' \
  --constraint mw_max=275 \
  --constraint logp_max=3 \
  --objective-term '{"name":"logS","path":"solubilities.O.solubilities.0","goal":"maximize","weight":1.0}' \
  --objective-term '{"name":"mw","source":"local","metric":"mw","goal":"minimize","weight":0.01}' \
  --candidates-per-iter 4 \
  --max-iter 12
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `aspirin_balanced` per AGENTS.md. Score aspirin first, confirm the
> logS path in the composite objective, then optimize the weighted score.

The agent should:

1. Score aspirin once and confirm that the `logS` term resolves from
   `solubilities.O.solubilities.0`.
2. Prefer compact polarity-increasing changes over heavy substituent additions.
3. Compare candidates by `objective_score`, not logS alone, because extra
   molecular weight carries a small penalty.

Expected: the best candidate may have slightly lower raw logS than a heavier
analog, but a better weighted score because it preserves MW headroom.
