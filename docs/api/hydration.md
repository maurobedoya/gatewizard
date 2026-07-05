# Packmol Hydration Module

Fill internal protein cavities and pores with TIP3P waters using standalone **PACKMOL** (AmberTools). Designed for heavy-atom structures from Visualize using **heavy-atom-safe inflated exclusion** radii.

- Check PACKMOL availability
- Detect hydrogen status on the protein
- Estimate free volume inside a 3D box
- Build and run PACKMOL input for cavity hydration
- Optional custom PACKMOL input (advanced)

## Import

```python
from gatewizard.tools.packmol_hydration import (
    check_packmol_available,
    detect_hydrogen_status,
    estimate_cavity_volume,
    build_hydrate_inp_text,
    prepare_hydration_job,
    hydrate_cavity,
    preview_hydrate_inp,
    run_custom_packmol,
    PackmolHydrationError,
)
```

---

## Function: check_packmol_available()

```python
check_packmol_available() -> dict
```

**Returns:** `{ "available": bool, "version": str|None, "resolved_path": str|None }`

### Example 1: Check PACKMOL availability

```python
from gatewizard.tools.packmol_hydration import check_packmol_available

info = check_packmol_available()
print(f"PACKMOL available: {info['available']}")
if info["resolved_path"]:
    print(f"Path: {info['resolved_path']}")
if info["version"]:
    print(f"Version: {info['version']}")
```

---

## Function: detect_hydrogen_status()

```python
detect_hydrogen_status(
    pdb_file: str | None = None,
    atoms: Sequence[Atom] | None = None,
) -> str
```

**Returns:** `"full"`, `"partial"`, or `"none"` (protein heavy-atom vs H ratio).

### Example 2: Detect hydrogen status on 6RV3_AB

```python
from pathlib import Path

from gatewizard.tools.packmol_hydration import detect_hydrogen_status

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
status = detect_hydrogen_status(pdb_file=PDB_FILE)
print(f"Hydrogen status for 6RV3_AB: {status}")
```

---

## Function: estimate_cavity_volume()

```python
estimate_cavity_volume(
    pdb_file: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    solute_radius: float | None = None,
    exclusion_mode: str | None = None,
    grid_spacing: float | None = None,
    atom_indices: Sequence[int] | None = None,
) -> VolumeEstimate
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | `str` | — | Input PDB path |
| `box_min` | `(x, y, z)` | — | Box minimum corner (Å) |
| `box_max` | `(x, y, z)` | — | Box maximum corner (Å) |
| `solute_radius` | `float` | auto | PACKMOL solute `radius`; default 2.5 Å (heavy-atom) or 1.5 Å (explicit) |
| `exclusion_mode` | `str` | auto | `"heavy_atom_safe"` or `"explicit"` |
| `grid_spacing` | `float` | `0.5` | Volume grid spacing (Å) |
| `atom_indices` | `list[int]` | `None` | Optional subset of atoms for exclusion |

**Returns:** `VolumeEstimate` dataclass (also `.as_dict()`).

### Example 3: Estimate cavity volume with a fixed box

```python
from pathlib import Path

from gatewizard.tools.packmol_hydration import estimate_cavity_volume

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

result = estimate_cavity_volume(
    pdb_file=PDB_FILE,
    box_min=BOX_MIN,
    box_max=BOX_MAX,
)
print(f"Box volume: {result.box_volume_A3:.1f} Å³")
print(f"Free volume: {result.free_volume_A3:.1f} Å³")
print(f"Suggested waters: {result.suggested_waters}")
print(f"Exclusion mode: {result.exclusion_mode}")
print(f"Hydrogen status: {result.hydrogen_status}")
```

---

## Function: build_hydrate_inp_text()

```python
build_hydrate_inp_text(
    protein_path: str,
    tip3p_path: str,
    output_pdb: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: int,
    solute_radius: float = 2.5,
    tolerance: float = 2.0,
    nloop: int = 20,
) -> str
```

**Returns:** PACKMOL input file contents as a string.

### Example 4: Build PACKMOL input text (preview)

```python
from pathlib import Path

from gatewizard.tools.packmol_hydration import build_hydrate_inp_text, estimate_cavity_volume

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

vol = estimate_cavity_volume(PDB_FILE, BOX_MIN, BOX_MAX)
n_waters = max(1, min(vol.suggested_waters, 50))

inp = build_hydrate_inp_text(
    protein_path=PDB_FILE,
    tip3p_path="TIP3P.pdb",
    output_pdb="6RV3_AB_hydrated.pdb",
    box_min=BOX_MIN,
    box_max=BOX_MAX,
    n_waters=n_waters,
    solute_radius=2.5,
)
print(inp)
```

---

## Function: prepare_hydration_job()

```python
prepare_hydration_job(
    pdb_file: str,
    job_dir: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: int,
    ...
) -> dict
```

Creates `job_dir`, copies the protein and `TIP3P.pdb`, writes `packmol.inp`.

### Example 5: Prepare hydration job files

```python
import tempfile
from pathlib import Path

from gatewizard.tools.packmol_hydration import estimate_cavity_volume, prepare_hydration_job

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

vol = estimate_cavity_volume(PDB_FILE, BOX_MIN, BOX_MAX)
n_waters = max(1, min(vol.suggested_waters, 50))

with tempfile.TemporaryDirectory() as tmp:
    job = prepare_hydration_job(
        pdb_file=PDB_FILE,
        job_dir=tmp,
        box_min=BOX_MIN,
        box_max=BOX_MAX,
        n_waters=n_waters,
    )
    print(f"Job dir: {job['job_dir']}")
    print(f"Input file: {job['packmol_inp_path']}")
    print(f"Output PDB name: {job['output_pdb_name']}")
    for name in ("packmol.inp", "TIP3P.pdb", Path(PDB_FILE).name):
        path = Path(tmp) / name
        print(f"  {name}: {'OK' if path.is_file() else 'MISSING'}")
```

---

## Function: hydrate_cavity()

```python
hydrate_cavity(
    pdb_file: str,
    working_dir: str,
    output_folder_name: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: int | None = None,
    ...
) -> HydrationJobResult
```

Prepares the job under `working_dir/output_folder_name`, runs PACKMOL, returns paths and log.

### Example 6: Full cavity hydration workflow

```python
import tempfile
from pathlib import Path

from gatewizard.tools.packmol_hydration import (
    check_packmol_available,
    estimate_cavity_volume,
    hydrate_cavity,
)

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

if check_packmol_available()["available"]:
    vol = estimate_cavity_volume(PDB_FILE, BOX_MIN, BOX_MAX)
    n_waters = max(1, min(vol.suggested_waters, 20))
    with tempfile.TemporaryDirectory() as tmp:
        result = hydrate_cavity(
            pdb_file=PDB_FILE,
            working_dir=tmp,
            output_folder_name="hydration_6RV3_AB",
            box_min=BOX_MIN,
            box_max=BOX_MAX,
            n_waters=n_waters,
        )
        print(f"Success: {result.success}")
        print(f"Output: {result.output_pdb}")
        print(f"Message: {result.message}")
else:
    print("PACKMOL not installed; skipping hydration run")
```

---

## Complete workflow

Typical GateWizard usage:

1. Load a heavy-atom PDB in **Visualize**
2. Define a hydration box around an internal cavity
3. Run **Example 3** to estimate volume (or use the GUI)
4. Run **Example 6** (or the GUI) to fill with TIP3P waters into `hydration_{pdb}/`
5. Load the hydrated PDB and continue with **Preparation → Builder**

When no hydrogens are present, the module uses **heavy-atom-safe** inflated radii and a default solute `radius` of **2.5 Å** so waters do not clash with polar H added later in Builder/tleap.

---

## Function: run_custom_packmol()

Run user-supplied PACKMOL input text (see GUI custom tab). Returns `{ job_dir, packmol_log, output_pdb, success, message }`.
