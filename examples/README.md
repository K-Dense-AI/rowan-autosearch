# Examples

Each example initializes a run under `runs/<run_id>/`, then hands the loop to an
agentic coding system such as Claude Code, Cursor, Codex, OpenCode, or another
agent that can follow `AGENTS.md`.

- `aspirin_solubility.md`: maximize aqueous logS for aspirin-like acids.
- `phenol_pka.md`: lower phenol pKa with electron-withdrawing substituents.
- `quinone_redox.md`: tune benzoquinone toward a more positive reduction
  potential.
- `ibuprofen_logp.md`: lower computed logP while preserving a simple
  arylpropionic-acid pharmacophore.
- `aspirin_balanced_objective.md`: optimize a weighted objective that rewards
  solubility and lightly penalizes molecular weight.
- `caffeine_tpsa.md`: push caffeine-like xanthines toward higher TPSA as a
  lower-CNS-penetration thought experiment.
- `lidocaine_logd.md`: lower predicted logD for lidocaine-like local
  anesthetics.
- `coumarin_dye_solubility.md`: water-solubilize a compact coumarin dye.
- `anthraquinone_battery_redox.md`: tune anthraquinone-like organic battery
  scaffolds toward higher reduction potential.

For a new workflow or metric path, score the starting molecule once, inspect the
latest iteration JSON, and set `metric_path` in the run config if extraction
misses the intended Rowan result.
