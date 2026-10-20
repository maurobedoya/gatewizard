# MemPrO Module

Module for orienting membrane proteins using [MemPrO](https://github.com/pstansfeld/MemPrO).
MemPrO positions proteins correctly in the membrane prior to system building with packmol-memgen.

- Run MemPrO orientation from Python
- Parse ranked orientation results
- Access oriented PDB files by rank
- Build command lines without execution
- Check tool availability

## Import

```python
from gatewizard.core.mempro import MemPrO, OrientationResult, MemProError
```

## Class: MemPrO

Main class for orienting membrane proteins. Wraps the `mempro` CLI tool.

### Constructor

```python
MemPrO()
```

**Parameters:** None

**Returns:** `MemPrO` instance

### Example 1: Create a MemPrO instance and check availability
```python
--8<-- "tests/mempro_examples/mempro_example_01.py"
```

---

## Static Methods

### is_available

Check whether the `mempro` executable is on PATH.

```python
MemPrO.is_available() -> bool
```

**Returns:** `True` if `mempro` is found, `False` otherwise.

### Example 2: Check if MemPrO is installed
```python
--8<-- "tests/mempro_examples/mempro_example_02.py"
```

---

## Core Methods

### run

Run MemPrO orientation on a PDB file.

```python
run(
    pdb_file: str,
    output_dir: Optional[str] = None,
    n_cpus: Optional[int] = None,
    n_iters: int = 150,
    grid_size: int = 36,
    dual_membrane: bool = False,
    peripheral: bool = False,
    use_weights: bool = False,
    flip: bool = False,
    membrane_thickness: Optional[float] = None,
    seed: Optional[int] = None,
    device: Optional[str] = None,
    single_precision: bool = False,
    rank: Optional[str] = None,
    predict_pg: bool = False,
    pg_guess: Optional[float] = None,
    charge: Optional[float] = None,
    charge_outer: Optional[float] = None,
    outer_membrane_thickness: Optional[float] = None,
    optimize_thickness: bool = False,
    transmembrane_residues: Optional[str] = None,
    write_bfactors: bool = False,
    curvature: bool = False,
    itp_file: Optional[str] = None,
    insane_path: Optional[str] = None,
    build_system: Optional[int] = None,
    build_arguments: Optional[str] = None,
    additional_residues: Optional[str] = None,
    additional_residues_itp: Optional[str] = None,
    residue_cg_file: Optional[str] = None,
    extra_args: Optional[List[str]] = None,
) -> List[OrientationResult]
```

**Parameters:**

A value that matches the MemPrO default is left off the command.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | `str` | — | Path to the input PDB file |
| `output_dir` | `str` | `None` | Output directory name (default: `Orient`) |
| `n_cpus` | `int` | `None` | CPU cores (`-nc`). Omit to use MemPrO's default |
| `n_iters` | `int` | `150` | Minimisation iterations (`-ni`) |
| `grid_size` | `int` | `36` | Starting configurations (`-ng`) |
| `dual_membrane` | `bool` | `False` | Dual membrane (`-dm`) |
| `peripheral` | `bool` | `False` | Peripheral protein (`-pr`). Not with dual membrane |
| `use_weights` | `bool` | `False` | Weight by B-factors (`-w`). Not with transmembrane residues |
| `flip` | `bool` | `False` | Flip in Z after orientation (`-flip`) |
| `membrane_thickness` | `float` | `None` | Inner membrane thickness in Å (`-mt`, default 28) |
| `seed` | `int` | `None` | Random seed (`-seed`). Omit for an unseeded run. The same seed, protein, and CPU count repeat |
| `device` | `str` | `None` | `auto`, `cpu`, or `gpu` (`-dev`). `auto` is omitted |
| `single_precision` | `bool` | `False` | 32-bit floats (`-sp`) |
| `rank` | `str` | `None` | `auto`, `h` (hits), or `p` (potential). `auto` is omitted |
| `predict_pg` | `bool` | `False` | Peptidoglycan layer (`-pg`). Needs dual membrane |
| `pg_guess` | `float` | `None` | Peptidoglycan position guess in Å (`-pg_guess`) |
| `charge` | `float` | `None` | Inner membrane charge (`-ch`, default 0) |
| `charge_outer` | `float` | `None` | Outer membrane charge (`-ch_o`). Needs dual membrane |
| `outer_membrane_thickness` | `float` | `None` | Outer thickness in Å (`-mt_o`, default 24). Needs dual membrane |
| `optimize_thickness` | `bool` | `False` | Optimise thickness (`-mt_opt`). Not with curvature |
| `transmembrane_residues` | `str` | `None` | Ranges such as `10-40,50-60` (`-tm`). Not with B-factor weights |
| `write_bfactors` | `bool` | `False` | Write the potential into B-factors (`-wb`). Not with curvature |
| `curvature` | `bool` | `False` | Curvature minimisation (`-c`). Not with dual membrane |
| `itp_file` | `str` | `None` | Martini topology (`-itp`) |
| `insane_path` | `str` | `None` | Path to insane (`-insane`) |
| `build_system` | `int` | `None` | Build a coarse-grained system for ranks below this (`-bd`). `0` is omitted |
| `build_arguments` | `str` | `None` | Extra insane arguments (`-bd_args`) |
| `additional_residues` | `str` | `None` | Comma-separated residue names (`-res`) |
| `additional_residues_itp` | `str` | `None` | Topology for those residues (`-res_itp`) |
| `residue_cg_file` | `str` | `None` | Folder of coarse-grained residue files (`-res_cg`) |
| `extra_args` | `list` | `None` | Additional CLI arguments passed verbatim |

**Returns:** List of `OrientationResult` sorted by rank.

**Raises:**

- `MemProError` — If `mempro` is not installed or execution fails.
- `FileNotFoundError` — If the input PDB file does not exist.

**Output Structure:**

```
Orient/
├── orientation.txt        # Summary of all ranked orientations
├── Rank_1/
│   └── oriented_rank_1.pdb
├── Rank_2/
│   └── oriented_rank_2.pdb
└── ...
```

### Example 3: Run MemPrO with default settings
```python
--8<-- "tests/mempro_examples/mempro_example_03.py"
```

### Example 4: Run with custom parameters
```python
--8<-- "tests/mempro_examples/mempro_example_04.py"
```

### Example 5: Run with dual membrane mode
```python
--8<-- "tests/mempro_examples/mempro_example_05.py"
```

### Example 6: Run with peripheral mode
```python
--8<-- "tests/mempro_examples/mempro_example_06.py"
```

---

### parse_results

Parse MemPrO results from an existing Orient directory (static method).

```python
MemPrO.parse_results(orient_dir: str) -> List[OrientationResult]
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `orient_dir` | `str` | Yes | Path to the Orient output directory |

**Returns:** List of `OrientationResult` sorted by rank.

**Raises:** `MemProError` — If `orientation.txt` is not found or cannot be parsed.

### Example 7: Parse results from an existing Orient directory
```python
--8<-- "tests/mempro_examples/mempro_example_07.py"
```

---

### get_oriented_pdb

Get the path to an oriented PDB file for a specific rank (static method).

```python
MemPrO.get_oriented_pdb(orient_dir: str, rank: int = 1) -> str
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `orient_dir` | `str` | — | Path to the Orient output directory |
| `rank` | `int` | `1` | Desired rank number |

**Returns:** Absolute path to the oriented PDB file.

**Raises:** `FileNotFoundError` — If the PDB for the rank does not exist.

### Example 8: Get the best oriented PDB
```python
--8<-- "tests/mempro_examples/mempro_example_08.py"
```

---

### build_command

Build the MemPrO command line without executing it.

```python
build_command(
    pdb_file: str,
    output_dir: Optional[str] = None,
    ...
) -> List[str]
```

Accepts the same parameters as `run()`.

**Returns:** List of command-line tokens.

### Example 9: Build and inspect command
```python
--8<-- "tests/mempro_examples/mempro_example_09.py"
```

---

## Class: OrientationResult

Data class representing a single ranked orientation from a MemPrO run.

### Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `rank` | `int` | Rank number (1 = best) |
| `relative_potential` | `float` | Relative potential score |
| `hits_pct` | `float` | Percentage of hits |
| `rerank_potential` | `float` | Re-ranked potential |
| `rerank_depth` | `float` | Re-rank minima depth |
| `rerank_value` | `float` | Re-rank value |
| `pdb_path` | `str` | Path to the oriented PDB file |

### Example 10: Inspect OrientationResult attributes
```python
--8<-- "tests/mempro_examples/mempro_example_10.py"
```

---

### compute_orientation_transform

```python
from gatewizard.core.mempro import compute_orientation_transform

R, t = compute_orientation_transform(source_pdb, oriented_pdb)
```

Compute the rigid transform that maps protein atoms in `source_pdb` onto a MemPrO `oriented_rank_*.pdb`. Matching uses standard amino-acid atoms keyed by `(chain, resid, atom_name)`.

**Returns:** Rotation matrix `R` and translation `t` such that `oriented ≈ R @ source + t`.

**Raises:** `MemProError` if fewer than three atom pairs match.

### apply_orientation_transform

```python
from gatewizard.core.mempro import apply_orientation_transform

out = apply_orientation_transform(source_pdb, oriented_pdb, output_pdb)
```

Apply the MemPrO transform to every atom in `source_pdb` (protein, ligands, water) and write `output_pdb`. Dummy MemPrO atoms are not copied.

See [Example 12](#example-12-parse-existing-results-and-load-into-structuremanager) and [`StructureManager.apply_mempro_orientation`](structure_manager.md#apply_mempro_orientation).

## Class: MemProError

Custom exception for MemPrO-related errors.

```python
from gatewizard.core.mempro import MemPrO, MemProError

mp = MemPrO()
try:
    results = mp.run("nonexistent.pdb")
except FileNotFoundError:
    print("PDB file not found")
except MemProError as e:
    print(f"MemPrO error: {e}")
```

---

## Complete Workflow

### Example 11: Full orientation workflow
```python
--8<-- "tests/mempro_examples/mempro_example_11.py"
```

### Example 12: Parse existing results and load into StructureManager
```python
--8<-- "tests/mempro_examples/mempro_example_12.py"
```
