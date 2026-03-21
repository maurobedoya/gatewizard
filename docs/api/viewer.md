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
```python
import os
import tempfile
from gatewizard.core.viewer import MolecularViewer, ViewerError

viewer = MolecularViewer()

# A small structure with HELIX/SHEET records
pdb_content = """\
HEADER    TEST PROTEIN
HELIX    1   1 ALA A    1  ALA A    4  1                                   4
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
    viewer.load_structure(tmp_path)

    # Default SS (assigned automatically at load time)
    print("SS after load (auto):")
    print(f"  {viewer.get_secondary_structure_summary()}")

    # Reassign using the heuristic method
    ss = viewer.assign_secondary_structure('heuristic')
    print(f"SS after heuristic: {ss}")

    # Reassign from PDB HELIX/SHEET records
    ss = viewer.assign_secondary_structure('pdb_records')
    print(f"SS after pdb_records: {ss}")

    # Try psique (may not be installed)
    try:
        ss = viewer.assign_secondary_structure('psique')
        print(f"SS after psique: {ss}")
    except ViewerError as e:
        print(f"psique not available: {e}")

    # Auto method (same priority as load)
    ss = viewer.assign_secondary_structure('auto')
    print(f"SS after auto: {ss}")
finally:
    os.unlink(tmp_path)
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
