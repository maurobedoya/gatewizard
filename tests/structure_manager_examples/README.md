# Viewer Examples

These examples demonstrate the `StructureManager` API for loading,
inspecting, selecting, editing, and saving protein structures.

All examples use the API only (no GUI) and can be run with:

```bash
python tests/viewer_examples/structure_manager_example_01.py
```

## Examples

| # | Description |
|---|-------------|
| 01 | Create a StructureManager and inspect defaults |
| 02 | Load a PDB structure and print summary info |
| 03 | Query chains, residues, secondary structure |
| 04 | Select atoms by criteria (protein, backbone, ligand, etc.) |
| 05 | Select by chain and residue range |
| 06 | Auto-detect molecules (protein, water, ligands) |
| 07 | Rename chains and residues |
| 08 | Renumber residues |
| 09 | Delete atoms and save modified PDB |
| 10 | Full workflow: load → select → edit → save |
| 11 | Reassign secondary structure (psique, heuristic, pdb_records) |
| 12 | Rotate atoms around an axis |
| 13 | Translate and center structure |
| 14 | Align structure to an axis |
| 15 | Align with primary and secondary axes |
