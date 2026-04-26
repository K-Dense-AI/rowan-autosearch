# Example: lower ibuprofen-like logP

This descriptor example uses Rowan's descriptor workflow to reduce lipophilicity
while preserving an arylpropionic-acid motif.

```bash
uv run rowan-state init \
  --run ibuprofen_logp \
  --objective "Minimize logP for an ibuprofen-like arylpropionic acid while preserving the carboxylic acid and keeping MW between 180 and 300" \
  --direction minimize \
  --metric logP \
  --workflow descriptors \
  --start-smiles "CC(C)Cc1ccc(C(C)C(=O)O)cc1" \
  --constraint mw_min=180 \
  --constraint mw_max=300 \
  --constraint hbd_max=2 \
  --constraint hba_max=5 \
  --candidates-per-iter 4 \
  --max-iter 10
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `ibuprofen_logp` per AGENTS.md. Score the parent first, verify the
> descriptor `logP` path, and search for lower-logP analogs.

The agent should:

1. Keep the carboxylic acid and arylpropionic-acid relationship intact.
2. Replace lipophilic alkyl bulk with smaller or more polar substituents.
3. Avoid adding too many donors or acceptors, since the objective is a modest
   lipophilicity reduction rather than a solubility-only search.

Expected: shortening the isobutyl group or introducing modest polarity should
lower logP without exceeding the descriptor constraints.
