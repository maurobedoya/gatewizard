# Study cases

End-to-end walkthroughs live here. API method pages still include the numbered scripts under `tests/<module>_examples/`.

- [STP10 in the desktop app](7aaq-gui.md) — PDB 7AAQ from Visualize through Analysis: MemPro orientation, a titration figure on the oriented PDB, pH 6, HID histidines, a pre-oriented POPC build, 0.15 M KCl, and 100 ns of NPT production.
- [STP10 with the API](7aaq-api.md) — the same system from `tests/study_cases/7aaq/`: MemPro, preparation, the pre-oriented membrane build, and the OpenMM inputs.
- [Prepare systems](prepare-systems.md) — POPC membranes for the systems in `prepare_systems.py`: what each PDB keeps, what the structure checker found, how long the finished builds took, and what to change before a new run.
