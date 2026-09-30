# Structure Manager examples

API-only scripts for `StructureManager` (no GUI). Included in
`docs/api/structure_manager.md` and run by
`pytest tests/test_structure_manager.py::TestStructureManagerExamples`.

```bash
python tests/structure_manager_examples/structure_manager_example_01.py
```

| # | Description |
|---|-------------|
| 01 | Create a StructureManager and inspect defaults |
| 02 | Load a PDB and print summary info |
| 03 | Query chains, residues, secondary structure |
| 04 | Select atoms by criteria |
| 05 | Select by chain and residue range |
| 06 | Auto-detect molecules |
| 07 | Rename chains and residues |
| 08 | Renumber residues |
| 09 | Delete atoms and save PDB |
| 10 | Full workflow: load → select → edit → save |
| 11 | Reassign secondary structure |
| 12 | Rotate around an axis |
| 13 | Translate and center |
| 14 | Align to an axis |
| 15 | Primary and secondary axes |
