# Example: make lidocaine less lipophilic

This ADMET example treats lidocaine as a starting point and asks the agent to
lower predicted logD without destroying the recognizable anilide anesthetic
motif.

```bash
uv run rowan-state init \
  --run lidocaine_logd \
  --objective "Minimize logD for lidocaine-like local anesthetics while preserving the anilide and tertiary amine motifs" \
  --direction minimize \
  --metric logD \
  --workflow admet \
  --start-smiles "CCN(CC)CC(=O)Nc1c(C)cccc1C" \
  --constraint mw_max=325 \
  --constraint logp_max=4 \
  --constraint hbd_max=2 \
  --constraint hba_max=6 \
  --constraint rotb_max=8 \
  --candidates-per-iter 4 \
  --max-iter 12
```

Then in Claude Code, Cursor, Codex, OpenCode, or another agentic coding system:

> Drive run `lidocaine_logd` per AGENTS.md. Score lidocaine first, verify `logD`
> extraction from the ADMET payload, then search for lower-logD analogs.

The agent should:

1. Preserve the amide linker and a basic tertiary amine in early iterations.
2. Try smaller amine substituents, heteroatom-containing side chains, and
   reduced aryl methyl substitution.
3. Avoid over-polarizing into permanently charged or obviously non-drug-like
   structures unless the run is stuck.

Expected: trimming hydrophobic aryl or dialkylamine bulk should reduce logD
while keeping enough motif similarity for a meaningful analog search.
