# Structure Manager Module

Module for loading, inspecting, editing, and saving molecular structures.
Uses MDAnalysis for PDB parsing and atom selections.

- Load PDB files from disk or by PDB ID
- Query chains, residues, secondary structure
- Select atoms by criteria (protein, backbone, ligand, chain, range, etc.)
- Auto-detect molecules (protein, water, ligands)
- Edit: rename chains/residues, renumber residues, delete atoms
- Save modified structures as PDB

## Import

```python
from gatewizard.core.structure_manager import StructureManager, Selection
from gatewizard import StructureManager  # Also available at top level
```

## Class: StructureManager

Main API class for programmatic structure viewing and editing.

### Constructor

```python
StructureManager()
```

Creates a new StructureManager instance. No structure is loaded initially.

### Example 1: Create a StructureManager and inspect defaults
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_01.py"
```

---

## Loading Methods

### load_structure

```python
viewer.load_structure(filepath: str) -> dict
```

Load a PDB file from disk.

**Parameters:**
- `filepath` (str): Path to the PDB file.

**Returns:** Dictionary with keys `n_atoms`, `n_residues`, `n_chains`, `n_bonds`, `title`.

**Raises:** `StructureError` if file not found or parse fails.

### Example 2: Load a PDB structure and print summary info
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_02.py"
```

---

### load_from_pdb_id

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

### get_structure_info

```python
viewer.get_structure_info() -> dict
```

Return summary of the loaded structure.

**Returns:** Dictionary with `n_atoms`, `n_residues`, `n_chains`, `n_bonds`, `title`.

---

### get_chains

```python
viewer.get_chains() -> dict
```

Return chain IDs and their residue counts.

**Returns:** `{"A": 150, "B": 120, ...}`

---

### get_residues

```python
viewer.get_residues(chain_id: str = None) -> list
```

List residues, optionally filtered by chain.

**Parameters:**
- `chain_id` (str, optional): Filter by chain. If `None`, returns all residues.

**Returns:** List of dicts with keys `name`, `seq_id`, `chain_id`, `n_atoms`, `ss`.

---

### get_secondary_structure_summary

```python
viewer.get_secondary_structure_summary() -> dict
```

Count residues by secondary structure type.

**Returns:** `{"H": 45, "E": 30, "C": 75, ...}`

### Example 3: Query chains, residues, secondary structure
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_03.py"
```

---

## Selection Methods

### select_atoms

```python
viewer.select_atoms(selection_string: str) -> list[int]
```

Select atom indices using MDAnalysis selection language (`"protein"`, `"backbone"`, `"resname LIG"`, `"around 5 protein"`).

**Parameters:**
- `selection_string` (str): MDAnalysis selection expression.

**Returns:** Atom indices into the current structure.

See [Example 4](#example-4-select-atoms-by-criteria-protein-backbone-ligand-etc) for the related convenience API.

### select_by_criteria

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
--8<-- "tests/structure_manager_examples/structure_manager_example_04.py"
```

### Example 5: Select by chain and residue range
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_05.py"
```

---

### auto_detect_molecules

```python
viewer.auto_detect_molecules() -> list[Selection]
```

Automatically group atoms into protein, water, and individual ligand selections.

**Returns:** List of `Selection` objects with sensible defaults (representation, color scheme).

### Example 6: Auto-detect molecules (protein, water, ligands)
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_06.py"
```

---

## Edit Methods

### rename_chain

```python
viewer.rename_chain(old_chain: str, new_chain: str) -> int
```

Rename all atoms in a chain.

**Parameters:**
- `old_chain` (str): Current chain ID.
- `new_chain` (str): New chain ID (1 character).

**Returns:** Number of atoms renamed.

---

### rename_residues

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
--8<-- "tests/structure_manager_examples/structure_manager_example_07.py"
```

---

### renumber_residues

```python
viewer.renumber_residues(chain_id: str, start: int, end: int, new_start: int) -> int
```

Renumber residues sequentially from `new_start`.

**Returns:** Number of atoms renumbered.

### Example 8: Renumber residues
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_08.py"
```

---

### rename_chain_by_indices

```python
viewer.rename_chain_by_indices(indices: list[int], new_chain: str) -> int
```

Rename the chain of only the selected atoms. Sibling atoms that share a residue id but are not in `indices` are left unchanged.

### rename_residues_by_indices

```python
viewer.rename_residues_by_indices(indices: list[int], new_name: str) -> int
```

Rename residues for the selected atoms only (does not rename a sibling residue with the same resid).

### renumber_residues_by_indices

```python
viewer.renumber_residues_by_indices(indices: list[int], new_start: int) -> int
```

Renumber residues that contain the selected atoms, starting at `new_start`.

### delete_atoms

```python
viewer.delete_atoms(indices: list) -> int
```

Remove atoms by index. Rebuilds residues, chains, and bonds.

**Returns:** Number of atoms removed.

---

### save_pdb

```python
viewer.save_pdb(filepath: str) -> str
```

Write the current structure to a PDB file.

**Returns:** Absolute path of the saved file.

**Raises:** `StructureError` if no structure is loaded.

### Example 9: Delete atoms and save modified PDB
```python
--8<-- "tests/structure_manager_examples/structure_manager_example_09.py"
```

---

## Coordinate Transformations

### rotate_atoms

```python
viewer.rotate_atoms(angle_degrees: float, axis: str,
                    indices: list = None, center: str = 'selection') -> int
```

Rotate atoms around a Cartesian axis.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `angle_degrees` | float | — | Rotation angle in degrees |
| `axis` | str | — | `'x'`, `'y'`, or `'z'` |
| `indices` | list[int] | `None` | Atom indices to rotate. `None` = all atoms |
| `center` | str | `'selection'` | `'selection'` rotates around the centroid of the affected atoms; `'origin'` rotates around (0, 0, 0) |

**Returns:** Number of atoms rotated.

---

### translate_atoms

```python
viewer.translate_atoms(displacement: list, indices: list = None) -> int
```

Translate atoms by a displacement vector.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `displacement` | list[float] | — | `[dx, dy, dz]` in angstroms |
| `indices` | list[int] | `None` | Atom indices to move. `None` = all atoms |

**Returns:** Number of atoms translated.

---

### center_atoms

```python
viewer.center_atoms(indices: list = None) -> numpy.ndarray
```

Move the entire structure so that the centroid of the selected atoms is at the
origin.  The shift is always applied to **all** atoms.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `indices` | list[int] | `None` | Atoms whose centroid defines the shift. `None` = all atoms |

**Returns:** The displacement applied (old centroid position) as a NumPy array.

---

### align_to_axis

```python
viewer.align_to_axis(primary_indices: list, target_axis: str = 'z',
                     secondary_indices: list = None,
                     secondary_axis: str = None,
                     apply_to: list = None) -> int
```

Align a selection's principal direction to a reference axis using SVD.  The
first singular vector fitted through the primary atom positions is rotated onto
the target axis.  An optional secondary alignment adds a rotation around the
primary axis so that the centroid of the secondary selection projects onto the
secondary axis.

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `primary_indices` | list[int] | — | Atoms whose principal direction defines the alignment vector |
| `target_axis` | str | `'z'` | `'x'`, `'y'`, or `'z'` |
| `secondary_indices` | list[int] | `None` | Atoms for the secondary axis alignment |
| `secondary_axis` | str | `None` | `'x'`, `'y'`, or `'z'`; must differ from `target_axis` |
| `apply_to` | list[int] | `None` | Atom indices to actually transform. `None` = all atoms |

**Returns:** Number of atoms transformed.

!!! tip "Use case: ion channel alignment"
    For channels it is common to select the filter ions (`name K`) and align
    them to the Z-axis, then add the selectivity filter residues as secondary
    alignment to the X-axis.  This orients the pore along Z with a known
    reference direction.

### Example 12: Rotate atoms around an axis

```python
--8<-- "tests/structure_manager_examples/structure_manager_example_12.py"
```

**Expected output:**
```
Loaded 8 atoms

Rotated 8 atoms 90° around Z:
  Atom 1 before: (0.0, 0.0, 0.0)
  Atom 1 after:  (4.5, -4.0, 0.0)

Rotated 4 atoms (residue 1) 45° around X

Rotated 8 atoms 180° around Y (origin):
  Atom 1: (0.0, 0.0, 0.0)
```

### Example 13: Translate and center structure

```python
--8<-- "tests/structure_manager_examples/structure_manager_example_13.py"
```

**Expected output:**
```
Translated 8 atoms by (5, -10, 0) Å:
  Atom 1 before: (10.0, 20.0, 30.0)
  Atom 1 after:  (15.0, 10.0, 30.0)

Translated 4 atoms (residue 2) by (0, 0, 5) Å

Centered structure:
  Centroid before: (14.2, 20.2, 30.0)
  Shift applied:   (14.2, 20.2, 30.0)
  Centroid after:  (0.0000, 0.0000, 0.0000)

Centered on residue 1:
  Residue 1 centroid: (0.0000, 0.0000, 0.0000)
```

### Example 14: Align structure to an axis

```python
--8<-- "tests/structure_manager_examples/structure_manager_example_14.py"
```

**Expected output:**
```
Before alignment (axis spans):
  X: 12.50 Å
  Y: 1.30 Å
  Z: 0.30 Å
  Principal axis: X (largest span)

Aligned 12 atoms to Z-axis (axis spans):
  X: 0.30 Å
  Y: 1.28 Å
  Z: 12.51 Å
  Principal axis: Z (largest span)

Aligned CA atoms to Y-axis, transformed all 12 atoms:
  X: 1.29 Å
  Y: 12.50 Å
  Z: 0.40 Å
```

### Example 15: Align with primary and secondary axes

```python
--8<-- "tests/structure_manager_examples/structure_manager_example_15.py"
```

**Expected output:**
```
Before alignment:
  Chain A spans: X=14.0, Y=0.4, Z=0.3
  Chain B centroid: (6.0, 4.1, 0.4)

Aligned 10 atoms (primary → Z, secondary → X):
  Chain A spans: X=0.30, Y=0.31, Z=14.00
  Chain B centroid: (10.07, 0.83, -0.73)
  Chain A now mostly along Z (Z span >> X, Y).
  Chain B centroid now has largest offset along X.

Aligned only chain A (8 atoms), chain B unchanged:
  Chain B moved: False
```

---

## Secondary Structure Assignment

### apply_mempro_orientation

```python
viewer.apply_mempro_orientation(oriented_pdb: str) -> int
```

Apply a MemPrO rigid-body orientation to the loaded structure. The oriented PDB supplies the target protein pose; every atom in the current structure (protein, lipids, ligands, water) receives the same transform.

**Parameters:**
- `oriented_pdb` (str): Path to a MemPrO `oriented_rank_*.pdb` file.

**Returns:** Number of atoms transformed.

See also [`compute_orientation_transform`](mempro.md#compute_orientation_transform) and [MemPrO Example 12](mempro.md#example-12-parse-existing-results-and-load-into-structuremanager).

### assign_secondary_structure

```python
viewer.assign_secondary_structure(method: str = 'auto') -> dict
```

Reassign secondary structure using a specific method.

**Parameters:**
- `method` (str): Assignment method. One of:
    - `'auto'` – PDB HELIX/SHEET records → psique → heuristic (default, same as initial load).
    - `'psique'` – Use the psique tool (raises `StructureError` if psique is not available).
    - `'heuristic'` – CA-angle heuristic (always available).
    - `'pdb_records'` – Only read HELIX/SHEET from the PDB file (raises `StructureError` if none found).

**Returns:** Updated secondary structure summary `{"H": n, "E": n, ...}`.

**Raises:** `StructureError` if the requested method is not available or fails.

### Example 11: Reassign secondary structure (psique, heuristic, pdb_records)

Uses backbone atoms from residues 1–20 of PDB 2MVJ (an alpha-helical region)
so that `psique` can compute secondary structure from geometry. No HELIX/SHEET
header records are included — this forces the `auto` method to fall through to
`psique` instead of just reading PDB records.

```python
--8<-- "tests/structure_manager_examples/structure_manager_example_11.py"
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
--8<-- "tests/structure_manager_examples/structure_manager_example_10.py"
```

---

### parse_pdb

```python
from gatewizard.core.structure_manager import parse_pdb

struct = parse_pdb(filepath: str) -> ProteinStructure
```

Parse a PDB file into a `ProteinStructure` without constructing a `StructureManager`.

### assign_secondary_structure_map

```python
from gatewizard.core.structure_manager import assign_secondary_structure_map

ss = assign_secondary_structure_map(filepath: str, method: str = "auto") -> dict
```

Return a per-residue secondary-structure map. `method` is `auto`, `psique`, `pdb_records`, or `heuristic`. Used by callers that do not need a full `StructureManager`.

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


