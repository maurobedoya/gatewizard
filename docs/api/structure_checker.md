# Structure checker

Look at a structure and report what is missing, broken, or ambiguous. The check does not build missing loops and does not score geometry.

Alternate conformations are the PDB altLoc letter. The likelihood is the occupancy, the column just before the B-factor. The most likely letter is the highest occupancy. A tie keeps the earlier letter (A before B).

## Import

```python
from gatewizard.core.structure_checker import (
    diagnose_structure,
    apply_altloc_choices,
    keep_most_likely_altlocs,
)
```

### diagnose_structure

```python
diagnose_structure(path) -> dict
```

Read a PDB, mmCIF, or GRO file with MDAnalysis. Only the first model is inspected. GRO files have no altLoc, so that list is empty.

The result has `summary` counts and these lists:

- `missing_atoms` — heavy atoms absent from a standard amino acid, including Amber names such as HID, ASH, and CYX. Each entry has the chain, residue number, atom name, and a `selection`.
- `missing_residues` — numbering jumps inside a chain, and SEQRES names that have no coordinates.
- `gaps` — backbone breaks where the previous C and the next N of consecutively numbered residues are farther than about 2 Å.
- `altlocs` — residues that have two or more altLoc letters. A ligand, water, lipid, or other hetero group is its own residue: the identity is chain, residue number, insertion code, and residue name. `default` is the most likely letter for that residue. Each option includes `occupancy`, how many atoms use that letter, and `likely`. `selection` is the residue, so one conformation can be shown at a time. Atoms that overlap and have no altLoc letter are listed too, including two copies of one ligand in the same site. Those groups are marked `inferred`.
- `chain_pairs` — protein chains that are at least 80% identical or nearly the same length. A pair says whether the sequences match, which residues are substituted, and whether one chain is longer.
- `issues` — more than one model, duplicate atoms, unknown residue names, occupancy 0, and a note when the protein has no hydrogens. Each entry has a `selection` for the atoms it describes.

### apply_altloc_choices

```python
apply_altloc_choices(path, choices=None, output_path=None) -> dict
```

Write a PDB that keeps one altLoc for each disordered residue. `choices` is a list of `{chain, resid, icode, resname, altloc}`. `resname` keeps a ligand separate from a protein residue that uses the same number. Without `resname`, the letter applies only when one residue has that chain, number, and insertion code. Add `atom` only to override a single atom. A residue with no choice keeps the highest occupancy. Kept atoms have a blank altLoc and occupancy 1. `ATOM` and `HETATM` stay in their original records. SEQRES, LINK, and CONECT stay in the file. The original file is not changed.

When `output_path` is omitted, the file is written next to the source as `{name}_checked.pdb`.

The result has `output_file`, `atoms_kept`, and `atoms_removed`.

### keep_most_likely_altlocs

```python
keep_most_likely_altlocs(path, output_path) -> dict
```

Write a PDB that keeps the higher-occupancy letter of every alternate conformation, including ligands, waters, and other molecules. A tie keeps the earlier letter. A file with no altLoc and no overlapping copies is copied unchanged. This is the same choice as Preparation → **Keep the most likely alternate conformations**.
