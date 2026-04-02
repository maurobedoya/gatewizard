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
from gatewizard.core.mempro import MemPrO

mp = MemPrO()
print(f"MemPrO available: {MemPrO.is_available()}")
```

---

## Static Methods

### Method: is_available()

Check whether the `mempro` executable is on PATH.

```python
MemPrO.is_available() -> bool
```

**Returns:** `True` if `mempro` is found, `False` otherwise.

### Example 2: Check if MemPrO is installed
```python
from gatewizard.core.mempro import MemPrO

if MemPrO.is_available():
    print("MemPrO is installed and ready to use")
else:
    print("Install MemPrO: pip install git+https://github.com/pstansfeld/MemPrO.git")
```

---

## Core Methods

### Method: run()

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
    extra_args: Optional[List[str]] = None,
) -> List[OrientationResult]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | `str` | — | Path to the input PDB file |
| `output_dir` | `str` | `None` | Output directory name (default: `Orient`) |
| `n_cpus` | `int` | `None` | Number of CPU cores (default: all available) |
| `n_iters` | `int` | `150` | Number of minimisation iterations |
| `grid_size` | `int` | `36` | Number of starting configurations |
| `dual_membrane` | `bool` | `False` | Enable dual membrane orientation |
| `peripheral` | `bool` | `False` | Enable peripheral protein orientation |
| `use_weights` | `bool` | `False` | Use B-factors to weight orientation |
| `flip` | `bool` | `False` | Flip protein in Z-axis after orientation |
| `membrane_thickness` | `float` | `None` | Initial membrane thickness in Å (default: 28) |
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
from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run("protein.pdb")
    for r in results:
        print(f"Rank {r.rank}: potential={r.relative_potential:.2f}, hits={r.hits_pct:.1f}%")
        print(f"  PDB: {r.pdb_path}")
```

### Example 4: Run with custom parameters
```python
from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run(
        "protein.pdb",
        output_dir="my_orient",
        n_cpus=4,
        n_iters=200,
        grid_size=72,
    )
    print(f"Found {len(results)} orientations")
```

### Example 5: Run with dual membrane mode
```python
from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run("protein.pdb", dual_membrane=True)
    print(f"Dual membrane: {len(results)} orientations")
```

### Example 6: Run with peripheral mode
```python
from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run("protein.pdb", peripheral=True)
    print(f"Peripheral: {len(results)} orientations")
```

---

### Method: parse_results()

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
from gatewizard.core.mempro import MemPrO

results = MemPrO.parse_results("Orient")
for r in results:
    print(f"Rank {r.rank}: potential={r.relative_potential:.2f}, "
          f"hits={r.hits_pct:.1f}%, pdb={r.pdb_path}")
```

---

### Method: get_oriented_pdb()

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
from gatewizard.core.mempro import MemPrO

try:
    pdb = MemPrO.get_oriented_pdb("Orient", rank=1)
    print(f"Best orientation: {pdb}")
except FileNotFoundError as e:
    print(f"Not found: {e}")
```

---

### Method: build_command()

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
from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    cmd = mp.build_command("protein.pdb", n_cpus=4, dual_membrane=True)
    print("Command:", " ".join(cmd))
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
from gatewizard.core.mempro import MemPrO

results = MemPrO.parse_results("Orient")
if results:
    best = results[0]
    print(f"Rank: {best.rank}")
    print(f"Relative potential: {best.relative_potential}")
    print(f"Hits: {best.hits_pct}%")
    print(f"Re-rank potential: {best.rerank_potential}")
    print(f"Re-rank depth: {best.rerank_depth}")
    print(f"Re-rank value: {best.rerank_value}")
    print(f"PDB path: {best.pdb_path}")
```

---

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
from gatewizard.core.mempro import MemPrO, MemProError

mp = MemPrO()

# Check availability
if not MemPrO.is_available():
    print("Install: pip install git+https://github.com/pstansfeld/MemPrO.git")
else:
    # Run orientation
    try:
        results = mp.run("protein.pdb", n_cpus=4)
        print(f"Found {len(results)} orientations\n")

        # Display results table
        print(f"{'Rank':>4}  {'Potential':>10}  {'Hits%':>6}  {'Re-rank':>10}")
        print("-" * 36)
        for r in results:
            print(f"{r.rank:>4d}  {r.relative_potential:>10.3f}  "
                  f"{r.hits_pct:>5.1f}%  {r.rerank_potential:>10.3f}")

        # Load the best orientation
        if results and results[0].pdb_path:
            print(f"\nBest PDB: {results[0].pdb_path}")
    except MemProError as e:
        print(f"Error: {e}")
```

### Example 12: Parse existing results and load into viewer
```python
from gatewizard.core.mempro import MemPrO
from gatewizard.core.viewer import MolecularViewer

# Parse pre-computed results
results = MemPrO.parse_results("Orient")

# Load the best orientation into the viewer
if results and results[0].pdb_path:
    viewer = MolecularViewer()
    info = viewer.load_structure(results[0].pdb_path)
    print(f"Loaded rank 1: {info['n_atoms']} atoms, {info['n_chains']} chains")
```
