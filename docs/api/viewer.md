# Viewer Module

Module for loading, inspecting, editing, and visualizing molecular structures.
Uses MDAnalysis for PDB parsing and VTK for 3D rendering.

- Load PDB files from disk or by PDB ID
- Query chains, residues, secondary structure
- Select atoms by criteria (protein, backbone, ligand, chain, range, etc.)
- Auto-detect molecules (protein, water, ligands)
- Edit: rename chains/residues, renumber residues, delete atoms
- Save modified structures as PDB
- Full VTK 3D visualization with multiple representations

## Import

```python
from gatewizard.core.viewer import MolecularViewer, Selection
from gatewizard import MolecularViewer  # Also available at top level
```

## Class: MolecularViewer

Main API class for programmatic structure viewing and editing.

### Constructor

```python
MolecularViewer()
```

Creates a new viewer instance. No structure is loaded initially.

### Example 1: Create a MolecularViewer and inspect defaults
```python
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()
print(f"MolecularViewer created: {viewer}")
print(f"Structure loaded: {viewer.structure is not None}")
```

---

## Loading Methods

### Method: load_structure()

```python
viewer.load_structure(filepath: str) -> dict
```

Load a PDB file from disk.

**Parameters:**
- `filepath` (str): Path to the PDB file.

**Returns:** Dictionary with keys `n_atoms`, `n_residues`, `n_chains`, `n_bonds`, `title`.

**Raises:** `ViewerError` if file not found or parse fails.

### Example 2: Load a PDB structure and print summary info
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

# Create a minimal PDB for testing
pdb_content = """\
HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    info = viewer.load_structure(tmp_path)
    print(f"Loaded: {info['n_atoms']} atoms, {info['n_residues']} residues")
    print(f"Chains: {info['n_chains']}, Bonds: {info['n_bonds']}")
    print(f"Title: {info.get('title', 'N/A')}")
finally:
    os.unlink(tmp_path)
```

---

### Method: load_from_pdb_id()

```python
viewer.load_from_pdb_id(pdb_id: str, output_dir: str = ".") -> dict
```

Download a PDB from RCSB and load it.

**Parameters:**
- `pdb_id` (str): 4-character PDB identifier (e.g. `"1CRN"`).
- `output_dir` (str): Directory to save the downloaded file.

**Returns:** Same dict as `load_structure`.

---

## Query Methods

### Method: get_structure_info()

```python
viewer.get_structure_info() -> dict
```

Return summary of the loaded structure.

**Returns:** Dictionary with `n_atoms`, `n_residues`, `n_chains`, `n_bonds`, `title`.

---

### Method: get_chains()

```python
viewer.get_chains() -> dict
```

Return chain IDs and their residue counts.

**Returns:** `{"A": 150, "B": 120, ...}`

---

### Method: get_residues()

```python
viewer.get_residues(chain_id: str = None) -> list
```

List residues, optionally filtered by chain.

**Parameters:**
- `chain_id` (str, optional): Filter by chain. If `None`, returns all residues.

**Returns:** List of dicts with keys `name`, `seq_id`, `chain_id`, `n_atoms`, `ss`.

---

### Method: get_secondary_structure_summary()

```python
viewer.get_secondary_structure_summary() -> dict
```

Count residues by secondary structure type.

**Returns:** `{"H": 45, "E": 30, "C": 75, ...}`

### Example 3: Query chains, residues, secondary structure
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00  0.00           O
ATOM      9  N   ALA B   1      11.000   2.000   3.000  1.00  0.00           N
ATOM     10  CA  ALA B   1      12.000   2.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    chains = viewer.get_chains()
    print(f"Chains: {chains}")

    residues = viewer.get_residues(chain_id='A')
    print(f"Chain A residues: {len(residues)}")
    for r in residues:
        print(f"  {r['name']} {r['seq_id']} ({r['n_atoms']} atoms, SS: {r['ss']})")

    ss = viewer.get_secondary_structure_summary()
    print(f"SS summary: {ss}")
finally:
    os.unlink(tmp_path)
```

---

## Selection Methods

### Method: select_by_criteria()

```python
viewer.select_by_criteria(criteria: str, extra: str = "") -> list
```

Select atoms using predefined criteria.

**Parameters:**
- `criteria` (str): One of `"All"`, `"Protein"`, `"Backbone"`, `"Sidechain"`, `"Water"`, `"Ligand"`, `"Chain..."`, `"Residue range..."`.
- `extra` (str): Required for `"Chain..."` (chain ID) and `"Residue range..."` (e.g. `"A:10-50"`).

**Returns:** List of atom indices.

### Example 4: Select atoms by criteria (protein, backbone, ligand, etc.)
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
HETATM    6  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
HETATM    7  C1  LIG A 200      30.000  30.000  30.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    protein = viewer.select_by_criteria('Protein')
    print(f"Protein atoms: {len(protein)}")

    backbone = viewer.select_by_criteria('Backbone')
    print(f"Backbone atoms: {len(backbone)}")

    sidechain = viewer.select_by_criteria('Sidechain')
    print(f"Sidechain atoms: {len(sidechain)}")

    water = viewer.select_by_criteria('Water')
    print(f"Water atoms: {len(water)}")

    ligand = viewer.select_by_criteria('Ligand')
    print(f"Ligand atoms: {len(ligand)}")
finally:
    os.unlink(tmp_path)
```

### Example 5: Select by chain and residue range
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      5  N   ALA B   1      11.000   2.000   3.000  1.00  0.00           N
ATOM      6  CA  ALA B   1      12.000   2.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    chain_a = viewer.select_by_criteria('Chain...', 'A')
    print(f"Chain A atoms: {len(chain_a)}")

    rng = viewer.select_by_criteria('Residue range...', 'A:1-2')
    print(f"A:1-2 atoms: {len(rng)}")

    all_atoms = viewer.select_by_criteria('All')
    print(f"All atoms: {len(all_atoms)}")
finally:
    os.unlink(tmp_path)
```

---

### Method: auto_detect_molecules()

```python
viewer.auto_detect_molecules() -> list[Selection]
```

Automatically group atoms into protein, water, and individual ligand selections.

**Returns:** List of `Selection` objects with sensible defaults (representation, color scheme).

### Example 6: Auto-detect molecules (protein, water, ligands)
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
HETATM    4  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
HETATM    5  C1  LIG A 200      30.000  30.000  30.000  1.00  0.00           C
HETATM    6  C2  LIG A 200      31.000  30.000  30.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    selections = viewer.auto_detect_molecules()
    for sel in selections:
        print(f"Selection '{sel.name}': {len(sel.atom_indices)} atoms, "
              f"rep={sel.representation}, cs={sel.color_scheme}")
finally:
    os.unlink(tmp_path)
```

---

## Edit Methods

### Method: rename_chain()

```python
viewer.rename_chain(old_chain: str, new_chain: str) -> int
```

Rename all atoms in a chain.

**Parameters:**
- `old_chain` (str): Current chain ID.
- `new_chain` (str): New chain ID (1 character).

**Returns:** Number of atoms renamed.

---

### Method: rename_residues()

```python
viewer.rename_residues(chain_id: str, start: int, end: int, new_name: str) -> int
```

Rename residues in a range.

**Parameters:**
- `chain_id` (str): Chain to modify.
- `start`, `end` (int): Residue number range (inclusive).
- `new_name` (str): New residue name.

**Returns:** Number of atoms renamed.

### Example 7: Rename chains and residues
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Rename chain A -> X
    count = viewer.rename_chain('A', 'X')
    print(f"Renamed {count} atoms from chain A to X")
    chains = viewer.get_chains()
    print(f"Chains after rename: {chains}")

    # Rename residue
    count = viewer.rename_residues('X', 1, 1, 'MET')
    print(f"Renamed {count} atoms to MET")
    residues = viewer.get_residues('X')
    for r in residues:
        print(f"  {r['name']} {r['seq_id']}")
finally:
    os.unlink(tmp_path)
```

---

### Method: renumber_residues()

```python
viewer.renumber_residues(chain_id: str, start: int, end: int, new_start: int) -> int
```

Renumber residues sequentially from `new_start`.

**Returns:** Number of atoms renumbered.

### Example 8: Renumber residues
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A  10       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A  10       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A  11       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A  11       5.000   1.000   3.000  1.00  0.00           C
ATOM      5  N   ALA A  12       7.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  ALA A  12       8.000   1.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Renumber residues 10-12 to start at 1
    count = viewer.renumber_residues('A', 10, 12, new_start=1)
    print(f"Renumbered {count} atoms")
    residues = viewer.get_residues('A')
    for r in residues:
        print(f"  {r['name']} {r['seq_id']}")
finally:
    os.unlink(tmp_path)
```

---

### Method: delete_atoms()

```python
viewer.delete_atoms(indices: list) -> int
```

Remove atoms by index. Rebuilds residues, chains, and bonds.

**Returns:** Number of atoms removed.

---

### Method: save_pdb()

```python
viewer.save_pdb(filepath: str) -> str
```

Write the current structure to a PDB file.

**Returns:** Absolute path of the saved file.

**Raises:** `ViewerError` if no structure is loaded.

### Example 9: Delete atoms and save modified PDB
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
HETATM    5  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

out_path = tmp_path + '_out.pdb'

try:
    viewer.load_structure(tmp_path)

    # Delete water atoms
    water_idx = viewer.select_by_criteria('Water')
    print(f"Deleting {len(water_idx)} water atoms")
    removed = viewer.delete_atoms(water_idx)
    print(f"Removed {removed} atoms")

    info = viewer.get_structure_info()
    print(f"After deletion: {info['n_atoms']} atoms")

    # Save modified structure
    saved = viewer.save_pdb(out_path)
    print(f"Saved to: {os.path.basename(saved)}")
finally:
    os.unlink(tmp_path)
    if os.path.exists(out_path):
        os.unlink(out_path)
```

---

## Secondary Structure Assignment

### Method: assign_secondary_structure()

```python
viewer.assign_secondary_structure(method: str = 'auto') -> dict
```

Reassign secondary structure using a specific method.

**Parameters:**
- `method` (str): Assignment method. One of:
    - `'auto'` – PDB HELIX/SHEET records → psique → heuristic (default, same as initial load).
    - `'psique'` – Use the psique tool (raises `ViewerError` if psique is not available).
    - `'heuristic'` – CA-angle heuristic (always available).
    - `'pdb_records'` – Only read HELIX/SHEET from the PDB file (raises `ViewerError` if none found).

**Returns:** Updated secondary structure summary `{"H": n, "E": n, ...}`.

**Raises:** `ViewerError` if the requested method is not available or fails.

### Example 11: Reassign secondary structure (psique, heuristic, pdb_records)

Uses backbone atoms from residues 1–20 of PDB 2MVJ (an alpha-helical region)
so that `psique` can compute secondary structure from geometry. No HELIX/SHEET
header records are included — this forces the `auto` method to fall through to
`psique` instead of just reading PDB records.

```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer, ViewerError

viewer = MolecularViewer()

# Backbone atoms (N, CA, C, O) from residues 1-20 of 2MVJ.
# Contains an alpha-helix (residues 5-19) that psique can detect.
pdb_content = """\
HEADER    BACKBONE OF 2MVJ RESIDUES 1-20
ATOM      1  N   MET A   1      12.463  -9.010  13.613  1.00  1.00           N
ATOM      2  CA  MET A   1      12.443  -7.552  13.317  1.00  1.00           C
ATOM      3  C   MET A   1      11.302  -6.897  14.099  1.00  1.00           C
ATOM      4  O   MET A   1      10.705  -7.511  14.982  1.00  1.00           O
ATOM      9  N   LYS A   2      11.011  -5.641  13.765  1.00  1.00           N
ATOM     10  CA  LYS A   2       9.947  -4.887  14.428  1.00  1.00           C
ATOM     11  C   LYS A   2       8.585  -5.480  14.086  1.00  1.00           C
ATOM     12  O   LYS A   2       7.578  -5.120  14.696  1.00  1.00           O
ATOM     18  N   PHE A   3       8.562  -6.358  13.080  1.00  1.00           N
ATOM     19  CA  PHE A   3       7.325  -6.995  12.609  1.00  1.00           C
ATOM     20  C   PHE A   3       6.359  -5.955  12.027  1.00  1.00           C
ATOM     21  O   PHE A   3       5.390  -6.304  11.351  1.00  1.00           O
ATOM     29  N   TYR A   4       6.650  -4.674  12.270  1.00  1.00           N
ATOM     30  CA  TYR A   4       5.852  -3.563  11.761  1.00  1.00           C
ATOM     31  C   TYR A   4       6.211  -3.358  10.285  1.00  1.00           C
ATOM     32  O   TYR A   4       5.510  -2.668   9.545  1.00  1.00           O
ATOM     41  N   THR A   5       7.302  -3.996   9.874  1.00  1.00           N
ATOM     42  CA  THR A   5       7.765  -3.930   8.497  1.00  1.00           C
ATOM     43  C   THR A   5       6.719  -4.546   7.560  1.00  1.00           C
ATOM     44  O   THR A   5       6.375  -3.963   6.532  1.00  1.00           O
ATOM     48  N   ILE A   6       6.190  -5.718   7.941  1.00  1.00           N
ATOM     49  CA  ILE A   6       5.155  -6.387   7.142  1.00  1.00           C
ATOM     50  C   ILE A   6       3.900  -5.510   7.126  1.00  1.00           C
ATOM     51  O   ILE A   6       3.257  -5.341   6.090  1.00  1.00           O
ATOM     56  N   LYS A   7       3.556  -4.976   8.296  1.00  1.00           N
ATOM     57  CA  LYS A   7       2.368  -4.139   8.440  1.00  1.00           C
ATOM     58  C   LYS A   7       2.417  -2.939   7.502  1.00  1.00           C
ATOM     59  O   LYS A   7       1.446  -2.649   6.802  1.00  1.00           O
ATOM     65  N   LEU A   8       3.537  -2.227   7.513  1.00  1.00           N
ATOM     66  CA  LEU A   8       3.683  -1.040   6.683  1.00  1.00           C
ATOM     67  C   LEU A   8       3.602  -1.416   5.212  1.00  1.00           C
ATOM     68  O   LEU A   8       2.948  -0.751   4.419  1.00  1.00           O
ATOM     73  N   ALA A   9       4.262  -2.505   4.854  1.00  1.00           N
ATOM     74  CA  ALA A   9       4.249  -2.968   3.470  1.00  1.00           C
ATOM     75  C   ALA A   9       2.803  -3.133   2.979  1.00  1.00           C
ATOM     76  O   ALA A   9       2.452  -2.708   1.880  1.00  1.00           O
ATOM     78  N   LYS A  10       1.968  -3.727   3.829  1.00  1.00           N
ATOM     79  CA  LYS A  10       0.551  -3.916   3.507  1.00  1.00           C
ATOM     80  C   LYS A  10      -0.150  -2.562   3.402  1.00  1.00           C
ATOM     81  O   LYS A  10      -0.936  -2.325   2.485  1.00  1.00           O
ATOM     87  N   PHE A  11       0.147  -1.681   4.350  1.00  1.00           N
ATOM     88  CA  PHE A  11      -0.448  -0.340   4.375  1.00  1.00           C
ATOM     89  C   PHE A  11      -0.163   0.373   3.058  1.00  1.00           C
ATOM     90  O   PHE A  11      -1.068   0.881   2.404  1.00  1.00           O
ATOM     98  N   LEU A  12       1.106   0.371   2.669  1.00  1.00           N
ATOM     99  CA  LEU A  12       1.514   0.988   1.405  1.00  1.00           C
ATOM    100  C   LEU A  12       0.887   0.229   0.237  1.00  1.00           C
ATOM    101  O   LEU A  12       0.425   0.834  -0.728  1.00  1.00           O
ATOM    106  N   GLY A  13       0.891  -1.100   0.324  1.00  1.00           N
ATOM    107  CA  GLY A  13       0.333  -1.932  -0.737  1.00  1.00           C
ATOM    108  C   GLY A  13      -1.169  -1.735  -0.857  1.00  1.00           C
ATOM    109  O   GLY A  13      -1.790  -2.148  -1.836  1.00  1.00           O
ATOM    110  N   GLY A  14      -1.757  -1.091   0.145  1.00  1.00           N
ATOM    111  CA  GLY A  14      -3.194  -0.806   0.137  1.00  1.00           C
ATOM    112  C   GLY A  14      -3.420   0.524  -0.569  1.00  1.00           C
ATOM    113  O   GLY A  14      -4.313   0.667  -1.404  1.00  1.00           O
ATOM    114  N   ILE A  15      -2.571   1.493  -0.221  1.00  1.00           N
ATOM    115  CA  ILE A  15      -2.618   2.835  -0.802  1.00  1.00           C
ATOM    116  C   ILE A  15      -2.152   2.813  -2.253  1.00  1.00           C
ATOM    117  O   ILE A  15      -2.687   3.540  -3.088  1.00  1.00           O
ATOM    122  N   VAL A  16      -1.114   2.020  -2.551  1.00  1.00           N
ATOM    123  CA  VAL A  16      -0.570   1.990  -3.912  1.00  1.00           C
ATOM    124  C   VAL A  16      -1.698   1.740  -4.925  1.00  1.00           C
ATOM    125  O   VAL A  16      -1.785   2.408  -5.956  1.00  1.00           O
ATOM    129  N   ARG A  17      -2.576   0.798  -4.587  1.00  1.00           N
ATOM    130  CA  ARG A  17      -3.730   0.477  -5.416  1.00  1.00           C
ATOM    131  C   ARG A  17      -4.695   1.668  -5.437  1.00  1.00           C
ATOM    132  O   ARG A  17      -5.274   1.987  -6.472  1.00  1.00           O
ATOM    140  N   ALA A  18      -4.880   2.304  -4.276  1.00  1.00           N
ATOM    141  CA  ALA A  18      -5.800   3.442  -4.169  1.00  1.00           C
ATOM    142  C   ALA A  18      -5.396   4.577  -5.114  1.00  1.00           C
ATOM    143  O   ALA A  18      -6.249   5.205  -5.738  1.00  1.00           O
ATOM    145  N   MET A  19      -4.092   4.828  -5.217  1.00  1.00           N
ATOM    146  CA  MET A  19      -3.578   5.878  -6.081  1.00  1.00           C
ATOM    147  C   MET A  19      -3.933   5.591  -7.539  1.00  1.00           C
ATOM    148  O   MET A  19      -4.546   6.409  -8.217  1.00  1.00           O
ATOM    153  N   LEU A  20      -3.542   4.417  -8.002  1.00  1.00           N
ATOM    154  CA  LEU A  20      -3.823   3.999  -9.374  1.00  1.00           C
ATOM    155  C   LEU A  20      -5.327   3.888  -9.559  1.00  1.00           C
ATOM    156  O   LEU A  20      -5.879   4.265 -10.593  1.00  1.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Default SS (assigned automatically at load time)
    # No HELIX records -> auto falls through to psique -> detects helix
    print("SS after load (auto):")
    print(f"  {viewer.get_secondary_structure_summary()}")

    # Reassign using the heuristic method
    ss = viewer.assign_secondary_structure('heuristic')
    print(f"SS after heuristic: {ss}")

    # Reassign from PDB HELIX/SHEET records (none in this file)
    try:
        ss = viewer.assign_secondary_structure('pdb_records')
        print(f"SS after pdb_records: {ss}")
    except ViewerError as e:
        print(f"pdb_records: {e}")

    # psique: runs the executable, writes temp PDB with HELIX records, parses SS
    try:
        ss = viewer.assign_secondary_structure('psique')
        print(f"SS after psique: {ss}")
    except ViewerError as e:
        print(f"psique error: {e}")

    # Auto method (same priority as load: pdb_records -> psique -> heuristic)
    ss = viewer.assign_secondary_structure('auto')
    print(f"SS after auto: {ss}")
finally:
    os.unlink(tmp_path)
```

**Expected output:**
```
SS after load (auto):
  {'C': 5, 'H': 15}
SS after heuristic: {'C': 4, 'H': 16}
pdb_records: No HELIX/SHEET records found in PDB file
SS after psique: {'C': 5, 'H': 15}
SS after auto: {'C': 5, 'H': 15}
```

---

## Full Workflow

### Example 10: Full workflow: load → select → edit → save
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A  50       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A  50       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A  50       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A  50       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A  51       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A  51       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A  51       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A  51       6.500   2.000   3.000  1.00  0.00           O
HETATM    9  O   HOH A 300      20.000  20.000  20.000  1.00  0.00           O
HETATM   10  C1  LIG B   1      30.000  30.000  30.000  1.00  0.00           C
HETATM   11  C2  LIG B   1      31.000  30.000  30.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name
out_path = tmp_path + '_edited.pdb'

try:
    # 1. Load
    info = viewer.load_structure(tmp_path)
    print(f"Loaded: {info['n_atoms']} atoms, {info['n_chains']} chains")

    # 2. Inspect
    print(f"Chains: {viewer.get_chains()}")
    sels = viewer.auto_detect_molecules()
    for s in sels:
        print(f"  Detected: {s.name} ({len(s.atom_indices)} atoms)")

    # 3. Edit: rename chain A -> X, renumber residues
    viewer.rename_chain('A', 'X')
    viewer.renumber_residues('X', 50, 51, new_start=1)

    # 4. Delete water
    water = viewer.select_by_criteria('Water')
    viewer.delete_atoms(water)

    # 5. Save
    viewer.save_pdb(out_path)
    print(f"Saved edited structure: {os.path.basename(out_path)}")

    # 6. Verify
    viewer2 = MolecularViewer()
    info2 = viewer2.load_structure(out_path)
    print(f"Verified: {info2['n_atoms']} atoms, chains: {viewer2.get_chains()}")
    residues = viewer2.get_residues('X')
    for r in residues:
        print(f"  {r['name']} {r['seq_id']}")
finally:
    os.unlink(tmp_path)
    if os.path.exists(out_path):
        os.unlink(out_path)
```

---

## Class: Selection

A named subset of atoms with display properties.

```python
Selection(name, atom_indices, *, representation='ball_stick',
          color_scheme='element', uniform_color=None, visible=True, ...)
```

### Attributes

| Attribute | Type | Default | Description |
|-----------|------|---------|-------------|
| `name` | str | — | Display name |
| `atom_indices` | list[int] | — | Indices into structure.atoms |
| `representation` | str | `'ball_stick'` | `'vdw'`, `'ball_stick'`, `'sticks'`, `'cartoon'`, `'tube_ss'`, `'backbone'`, `'surface'` |
| `color_scheme` | str | `'element'` | `'element'`, `'chain'`, `'ss'`, `'uniform'` |
| `uniform_color` | tuple | `None` | RGB tuple `(r, g, b)` when color_scheme is `'uniform'` |
| `visible` | bool | `True` | Show/hide |
| `quality` | int | `3` | 1–5, controls mesh resolution |
| `opacity` | float | `0.5` | Surface opacity |

---

## Representations

The viewer supports seven molecular representations:

| Key | Name | Description |
|-----|------|-------------|
| `vdw` | VDW (Spacefill) | Atoms as spheres at van der Waals radii |
| `ball_stick` | Ball & Stick | Small spheres + bond sticks |
| `sticks` | Sticks | Bond sticks only |
| `cartoon` | Cartoon | Ribbon diagram with helix/sheet/coil |
| `tube_ss` | Tube SS | Colored tubes by secondary structure |
| `backbone` | Backbone | CA trace as tube |
| `surface` | Surface | Molecular surface |

---

## GUI: VisualizeFrame

The `VisualizeFrame` in `gatewizard.gui.frames.visualize` provides the full
interactive VTK-based 3D viewer with:

- Load/download PDB files
- Multiple selections with independent representations and colors
- Drag-reorder selections
- Per-selection quality, size, material, and SS color settings
- Edit operations (rename chain, rename/renumber residues, delete atoms)
- SSAO ambient occlusion, shadows, depth cueing
- Save high-resolution images (PNG, JPEG, TIFF, BMP) with configurable scale
- Save/load viewpoints as JSON
- Save modified structures as PDB
