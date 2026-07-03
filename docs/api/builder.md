# Builder Module

Module for preparing membrane protein systems with automated lipid bilayer construction, solvation, and parametrization. This module includes:

- Membrane system building with packmol-memgen
- Custom lipid composition (asymmetric membranes supported)
- Force field selection and validation
- Salt addition and ionic strength control
- Two-stage Propka workflow support
- Automated parametrization with pdb4amber and tleap

## Import

```python
from gatewizard.core.builder import Builder
from gatewizard.tools.force_fields import ForceFieldManager
```

## Class: Builder

Main class for building membrane protein systems with complete control over lipid composition, force fields, and system parameters.

### Constructor

```python
Builder()
```

**Parameters:** None

**Returns:** `Builder` instance

**Default Configuration:**

| Parameter | Default | Description |
|-----------|---------|-------------|
| `water_model` | `"tip3p"` | Water model for solvation |
| `protein_ff` | `"ff19SB"` | Protein force field |
| `lipid_ff` | `"lipid21"` | Lipid force field |
| `preoriented` | `True` | Protein already oriented for membrane insertion |
| `parametrize` | `True` | Run parametrization with tleap |
| `salt_concentration` | `0.15` | Salt concentration in M (molar) |
| `cation` | `"K+"` | Cation type for salt |
| `anion` | `"Cl-"` | Anion type for salt |
| `dist` | `12` | Minimum solute-to-box-boundary distance in Å |
| `dist_wat` | `26` | Water layer thickness in Å |
| `notprotonate` | `False` | Skip protonation (preserve residue names) |

### Example 1: Basic Configuration
```python
from gatewizard.core.builder import Builder

builder = Builder()
print(f"Default water model: {builder.config['water_model']}")
print(f"Default protein force field: {builder.config['protein_ff']}")
print(f"Default lipid parameters: {builder.config['lipid_ff']}")
```

---

## Configuration Methods

### Method: set_configuration()

Update configuration parameters for system building.

```python
set_configuration(**kwargs)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `water_model` | `str` | `"tip3p"` | Water model (tip3p, tip4p, spce, opc) |
| `protein_ff` | `str` | `"ff19SB"` | Protein force field (ff14SB, ff19SB) |
| `lipid_ff` | `str` | `"lipid21"` | Lipid force field (lipid17, lipid21) |
| `preoriented` | `bool` | `True` | Whether protein is pre-oriented |
| `parametrize` | `bool` | `True` | Run parametrization after packing |
| `salt_concentration` | `float` | `0.15` | Salt concentration in M |
| `cation` | `str` | `"K+"` | Cation type (Na+, K+, etc.) |
| `anion` | `str` | `"Cl-"` | Anion type (Cl-, Br-, etc.) |
| `dist` | `float` | `12` | Minimum solute-to-box-boundary distance in Å |
| `dist_wat` | `float` | `26` | Water layer thickness in Å |
| `notprotonate` | `bool` | `False` | Skip protonation during parametrization |
| `add_salt` | `bool` | `True` | Whether to add salt to system |

**Returns:** None

### Example 2: Custom Configuration
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Configure for specific system
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    salt_concentration=0.15,
    cation="Na+",
    anion="Cl-",
    dist=12,  # Minimum distance to box boundaries
    dist_wat=20.0,  # Larger water layer
    preoriented=True
)
```

### Example 3: Available Water Models
```python
from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get available water models
water_models = ff_manager.get_water_models()

print("\nAvailable Water Models:")
for water in water_models:
    print(f"  - {water}")

print(f"\nTotal: {len(water_models)} water models available")
```

Similarly can be done for lipids:

### Example 4: Available Lipid Models
```python
from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get all available lipids
lipids = ff_manager.get_available_lipids()

print(f"Total available lipids: {len(lipids)}\n")

for lipid in lipids:
     print(f"  - {lipid}")

print(f"\nTotal: {len(lipids)} lipid models available")
```


### Example 5: Available Protein Force Fields
```python
from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get available protein force fields
protein_ffs = ff_manager.get_protein_force_fields()

print("\nAvailable Protein Force Fields:")
for protein_ff in protein_ffs:
    print(f"  - {protein_ff}")

print(f"\nTotal: {len(protein_ffs)} protein force fields available")
```

### Example 6: Available Lipid Force Fields
```python
from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get available lipid force fields
lipid_ffs = ff_manager.get_lipid_force_fields()

print("\nAvailable Lipid Force Fields:")
for lipid_ff in lipid_ffs:
    print(f"  - {lipid_ff}")

print(f"\nTotal: {len(lipid_ffs)} lipid force fields available")
```

---

## Validation Methods

Before preparing systems, you can validate your inputs and check force field compatibility.

### Method: validate_system_inputs()

Validate all inputs before system preparation.

```python
validate_system_inputs(
    pdb_file: str,
    upper_lipids: List[str],
    lower_lipids: List[str],
    lipid_ratios: str = "",
    **kwargs
) -> Tuple[bool, str]
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pdb_file` | `str` | Yes | Path to input PDB file |
| `upper_lipids` | `List[str]` | Yes | List of lipids for upper leaflet |
| `lower_lipids` | `List[str]` | Yes | List of lipids for lower leaflet |
| `lipid_ratios` | `str` | No | Lipid ratios string |
| `**kwargs` | `Any` | No | Additional parameters to validate |

**Returns:** `Tuple[bool, str]`

- `bool`: Validation status
- `str`: Error message (empty if valid)

**Validation Checks:**

1. PDB file exists and is readable
2. All lipid types are valid/available
3. Lipid ratios match number of lipids
4. Ratios are positive numbers
5. Force field combinations are compatible
6. Ion types are valid

### Example 7: Input Validation
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Validate inputs before preparation
valid, msg = builder.validate_system_inputs(
    pdb_file="protein_protonated_prepared.pdb",
    upper_lipids=["POPC", "POPE"],
    lower_lipids=["POPC", "POPE"],
    lipid_ratios="7:3//7:3",
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21"
)

if valid:
    if "WARNING" in msg:
        print(f"[!] Inputs valid with warning: {msg}")
        print("You may proceed at your own risk")
    else:
        print("[OK] All inputs are valid, proceed with preparation")
    # Now call prepare_system()
else:
    print(f"[ERROR] Validation failed: {msg}")
    # Fix issues before proceeding
```

### Example 8: Force Field Validation
```python
from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()
valid, message, is_warning = ff_manager.validate_combination("tip3p", "ff14SB", "lipid21")

if valid:
    if is_warning:
        print(f"[!] Warning: {message}")
    else:
        print("[OK] Force field combination is compatible")
else:
    print(f"[ERROR] Incompatible: {message}")
```

---

## Core Preparation Methods

### Method: prepare_system()

Prepare a complete membrane protein system with lipid bilayer, water, and ions.

```python
prepare_system(
    pdb_file: str,
    working_dir: str,
    upper_lipids: List[str],
    lower_lipids: List[str],
    lipid_ratios: str = "",
    **kwargs
) -> Tuple[bool, str, Optional[Path]]
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `pdb_file` | `str` | Yes | Path to input PDB file (oriented protein) |
| `working_dir` | `str` | Yes | Working directory for output files |
| `upper_lipids` | `List[str]` | Yes | List of lipid types for upper leaflet |
| `lower_lipids` | `List[str]` | Yes | List of lipid types for lower leaflet |
| `lipid_ratios` | `str` | No | Lipid molar ratios (format: "ratio1:ratio2//ratio3:ratio4") |
| `output_folder_name` | `str` | No | Custom output folder name (kwarg) |
| `wait` | `bool` | No | Block until the job completes or errors (default `False`) |
| `wait_timeout` | `float` | No | Maximum seconds to wait when `wait=True` (`None` = unlimited) |
| `wait_poll_interval` | `float` | No | Seconds between status checks (default `5.0`) |
| `wait_verbose` | `bool` | No | Print elapsed-time progress while waiting (default `True`) |
| `**kwargs` | `Any` | No | Additional configuration parameters (override defaults) |

**Returns:** `Tuple[bool, str, Optional[Path]]`

- `bool`: Success status
- `str`: Status message
- `Optional[Path]`: Job directory path (None if failed)

**Lipid Ratios Format:**

The `lipid_ratios` string specifies molar ratios for each leaflet:

- Format: `"upper_ratio1:upper_ratio2//lower_ratio1:lower_ratio2"`
- Ratios must match the order of lipids in the lipid lists
- Ratios are normalized automatically (don't need to sum to 1.0)
- Use `//` to separate upper and lower leaflet ratios

**Available Lipids:**

| Category | Lipid Types |
|----------|------------|
| **Phosphatidylcholine (PC)** | DHPC, DLPC, DMPC, DPPC, DSPC, DOPC, POPC, PLPC, SOPC |
| **Phosphatidylethanolamine (PE)** | DHPE, DLPE, DMPE, DPPE, DSPE, DOPE, POPE, PLPE, SOPE |
| **Phosphatidylserine (PS)** | DHPS, DLPS, DMPS, DPPS, DSPS, DOPS, POPS, PLPS, SOPS |
| **Phosphatidylglycerol (PG)** | DHPG, DLPG, DMPG, DPPG, DSPG, DOPG, POPG, PLPG, SOPG |
| **Phosphatidic Acid (PA)** | DHPA, DLPA, DMPA, DPPA, DSPA, DOPA, POPA, PLPA, SOPA |
| **Sphingomyelins** | PSM, ASM, LSM, MSM, HSM, OSM |
| **Cholesterol** | CHL1 |
| **Other** | CARDIOLIPIN, PIP, PIP2 |

**Output Files:**

In `{output_folder_name}/` directory:

- `bilayer_*.pdb` - Packed system (protein in membrane)
- `bilayer_*_lipid.pdb` - Lipid bilayer with CRYST1 info (for equilibration)
- `system.prmtop` - AMBER topology file (if parametrize=True)
- `system.inpcrd` - AMBER coordinate file (if parametrize=True)
- `system_solv.pdb` - Solvated system in PDB format
- `system_4amber.pdb` - PDB prepared by pdb4amber
- `tleap.in` - tleap input script
- `logs/preparation.log` - Main execution log
- `logs/parametrization.log` - Parametrization log
- `status.json` - Job status tracking
- `run_preparation.sh` - Execution script

**Raises:**

- `FileNotFoundError` - If input PDB file doesn't exist
- `BuilderError` - If validation fails or system preparation fails

**Note:** By default, system preparation runs in the background and the method returns immediately.
Pass `wait=True` to block until the job finishes — useful for scripting sequential preparations.
For asynchronous monitoring, use the log files or the `JobMonitor` class (see examples below).

### Example 9: Simple Symmetric Membrane
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Configure system
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    salt_concentration=0.15,
    cation="K+",
    anion="Cl-"
)

# Prepare system with 100% POPC (symmetric)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",  # 100% POPC both leaflets
    output_folder_name="popc_membrane"
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor progress: {job_dir / 'preparation.log'}")
    print(f"  When complete, files will be at:")
    print(f"    - Topology: {job_dir / 'system.prmtop'}")
    print(f"    - Coordinates: {job_dir / 'system.inpcrd'}")
else:
    print(f"✗ Preparation failed: {message}")
```

### Example 10: Monitoring Job Progress
```python
from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path

# Create monitor for your working directory
monitor = JobMonitor(working_directory=Path("./systems"))

# Scan for jobs
monitor.scan_for_jobs(force=True)

# Get active jobs
active_jobs = monitor.get_active_jobs()

if active_jobs:
    print(f"✓ Found {len(active_jobs)} active job(s)")
    
    for job_id, job_info in active_jobs.items():
        print(f"Job: {job_info.job_dir.name}")
        print(f"  Status: {job_info.status.value}")
        print(f"  Progress: {job_info.progress:.1f}%")
        print(f"  Current step: {job_info.current_step}")
        print(f"  Elapsed: {job_info.elapsed_time:.1f}s")
        
        # Show completed steps
        if job_info.steps_completed:
            print(f"  Completed steps:")
            for step in job_info.steps_completed:
                print(f"    ✓ {step}")
else:
    print("No active jobs found")

# Check for completed jobs
completed_jobs = monitor.get_completed_jobs()
if completed_jobs:
    print(f"✓ Found {len(completed_jobs)} completed job(s)")
    for job_id, job_info in completed_jobs.items():
        print(f"  - {job_info.job_dir.name}: {job_info.status.value}")
```

### Example 11: Asymmetric Membrane with Multiple Lipids
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Configure for asymmetric membrane
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    # Note: When using anionic lipids (POPS, POPG, etc.), the system may require
    # a higher salt concentration for neutralization. If you get an error like:
    # "The concentration of ions required to neutralize the system is higher than
    # the concentration specified", increase salt_concentration or use add_salt=False.
    # POPS is anionic (-1 charge), so 50% POPS in lower leaflet adds significant
    # negative charge requiring more cations for neutralization.
    salt_concentration=0.5,  # Increased from default 0.15 M due to anionic lipids
    dist_wat=20.0  # Thicker water layer
)

# Upper leaflet: 70% POPC + 30% cholesterol
# Lower leaflet: 50% POPE + 50% POPS (anionic)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC", "CHL1"],
    lower_lipids=["POPE", "POPS"],
    lipid_ratios="7:3//5:5",  # Ratios normalized automatically
    output_folder_name="asymmetric_membrane"
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Upper leaflet: 70% POPC, 30% CHL1")
    print(f"  Lower leaflet: 50% POPE, 50% POPS")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")
```

### Example 12: Complex Composition (Plasma Membrane Mimic)
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Plasma membrane-like composition
# Upper: PC-rich with cholesterol
# Lower: PE/PS-rich (cytoplasmic side)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC", "POPE", "CHL1"],
    lower_lipids=["POPE", "POPS", "CHL1"],
    lipid_ratios="5:2:3//4:4:2",  # Upper: 50% POPC, 20% POPE, 30% CHL1
                                   # Lower: 40% POPE, 40% POPS, 20% CHL1
    output_folder_name="plasma_membrane",
    salt_concentration=0.5,
    cation="Na+",  # Use sodium instead of potassium
    dist_wat=25.0  # Extra water for large protein
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")
```

### Example 13: Packing Only (No Parametrization)
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Only pack the system, don't parametrize
# Useful for visual inspection before parametrization
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="packed_only",
    parametrize=False  # Skip parametrization
)

if success:
    print(f"✓ Packing started in background (no parametrization)")
    print(f"  {message}")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
    print(f"  When complete, inspect: {job_dir / 'bilayer_*.pdb'}")
else:
    print(f"✗ Preparation failed: {message}")
```

### Example 14: Custom Salt Concentration
```python
from gatewizard.core.builder import Builder

builder = Builder()

# High salt concentration for ionic strength studies
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC", "POPS"],
    lower_lipids=["POPC", "POPS"],
    lipid_ratios="8:2//8:2",  # 80% POPC, 20% POPS
    output_folder_name="high_salt",
    salt_concentration=2.0,  # 2000 mM (high salt)
    cation="Na+",
    anion="Cl-"
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  High salt: 2.0 M NaCl")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")
```

### Example 15: No Salt (Charge Neutralization Only)
```python
from gatewizard.core.builder import Builder

builder = Builder()

# Neutralize system charges only (no extra salt)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="no_salt",
    salt_concentration=0.0,  # Only neutralize, no extra salt
    add_salt=True  # Still add ions for neutralization
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Neutralization only, no extra salt")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")
```

---

### Method: wait_for_completion()

Block until a preparation job completes or fails by polling `status.json`.

```python
wait_for_completion(
    job_dir,
    poll_interval: float = 5.0,
    timeout: Optional[float] = None,
    verbose: bool = True,
) -> Tuple[bool, str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `job_dir` | `Path \| str` | — | Path to the job directory |
| `poll_interval` | `float` | `5.0` | Seconds between status checks |
| `timeout` | `float \| None` | `None` | Maximum seconds to wait (`None` = no limit) |
| `verbose` | `bool` | `True` | Print elapsed-time progress to stdout |

**Returns:** `Tuple[bool, str]` — (success, message)

### Example 25: Blocking (synchronous) Preparation with wait=True
```python
"""
Builder Example 25: Blocking (synchronous) preparation with wait=True

Demonstrates how to run prepare_system() in blocking mode so that the
Python script does not continue until the job finishes (or errors).
This is especially useful for scripting multiple sequential preparations.

Two approaches are shown:
  A) Inline ``wait=True`` inside prepare_system()
  B) Launch in background + call wait_for_completion() later

NOTE: This example requires AmberTools and packmol-memgen.
      It uses the real test PDB file tests/2MVJ_2ligs.pdb.
"""

from gatewizard.core.builder import Builder
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    parametrize_all_ligands,
)

pdb_file = "tests/2MVJ_2ligs.pdb"
working_dir = "./systems"

# ── Step 1: Detect and parametrize ligands (same as example 23) ──────
print("Step 1: Detecting ligands...")
ligands = detect_ligands(pdb_file)
for lig in ligands:
    print(f"  Found: {lig.name} ({lig.num_atoms} atoms, {lig.formula})")

print("\nStep 2: Parametrizing ligands...")
ligand_results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=f"{working_dir}/ligand_params",
    charges={"AAA": 0, "BBB": 0},
    charge_method="bcc",
)
print(f"  Parametrized: {list(ligand_results.keys())}")

# ── Step 3: Configure builder ────────────────────────────────────────
builder = Builder()
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    preoriented=True,
    parametrize=True,
    salt_concentration=0.15,
    dist=12,
    dist_wat=26,
    notprotonate=True,
    ligand_params=ligand_results,
)

# ── Approach A: inline wait ──────────────────────────────────────────
# With wait=True the call blocks until the job finishes.
# prepare_system returns (success, message, job_dir) where
#   success reflects the *final* outcome, not just "launched OK".
#
# NOTE: The actual preparation is long-running.  Set a short timeout
#       so the example returns quickly during tests; remove or increase
#       the timeout for real production runs.
print("\n── Approach A: prepare_system with wait=True ──")
success_a, message_a, job_dir_a = builder.prepare_system(
    pdb_file=pdb_file,
    working_dir=working_dir,
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="membrane_2MVJ_wait",
    wait=True,                # ← block until done or error
    wait_timeout=10,          # short timeout for demo (use 3600+ for real runs)
    wait_poll_interval=2,     # check every 2 s
    wait_verbose=True,        # print elapsed time
)
print(f"Result : {success_a}")
print(f"Message: {message_a}")

# ── Approach B: launch then wait later ───────────────────────────────
print("\n── Approach B: launch + wait_for_completion() ──")
success_b, message_b, job_dir_b = builder.prepare_system(
    pdb_file=pdb_file,
    working_dir=working_dir,
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="membrane_2MVJ_bg",
)

if success_b and job_dir_b is not None:
    # Do other work here …
    print("Doing other work while system builds …")

    # Then block until that specific job is done
    completed, wait_msg = builder.wait_for_completion(
        job_dir_b,
        poll_interval=2,
        timeout=10,          # short timeout for demo
        verbose=True,
    )
    print(f"Completed: {completed}")
    print(f"Message  : {wait_msg}")
else:
    print(f"Not launched: {message_b}")

print("\nExample 25 finished.")
```

---

## Job Monitoring

Since system preparation runs in the background, GateWizard provides the `JobMonitor` class to track progress, check status, and manage running jobs.

### Class: JobMonitor

Monitor and track system preparation jobs.

### Example 16: JobMonitor class
```python
from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path

monitor = JobMonitor(working_directory=Path("./systems"))
```

**Constructor Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `working_directory` | `Path` or `str` | `Path.cwd()` | Directory to monitor for jobs |

**Returns:** `JobMonitor` instance

---

### Method: scan_for_jobs()

Scan the working directory for preparation jobs.

```python
monitor.scan_for_jobs(force=False)
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `force` | `bool` | `False` | Force scan even if within scan interval (3 seconds) |

**Returns:** None (updates internal job list)

**Description:** Scans the working directory for subdirectories containing `status.json` files, which indicate active or completed preparation jobs.

---

### Method: get_active_jobs()

Get all currently running jobs.

```python
active_jobs = monitor.get_active_jobs()
```

**Returns:** `Dict[str, JobInfo]`

- Dictionary mapping job IDs to `JobInfo` objects
- Only includes jobs with status `RUNNING`

---

### Method: get_completed_jobs()

Get all completed or errored jobs.

```python
completed_jobs = monitor.get_completed_jobs()
```

**Returns:** `Dict[str, JobInfo]`

- Dictionary mapping job IDs to `JobInfo` objects
- Includes jobs with status `COMPLETED` or `ERROR`

---

### Method: get_job()

Get specific job by ID.

```python
job = monitor.get_job(job_id)
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `job_id` | `str` | Yes | Job identifier (usually the full job directory path) |

**Returns:** `Optional[JobInfo]`

- `JobInfo` object if job exists
- `None` if job not found

---

### Method: refresh_job()

Force refresh of a specific job's status.

```python
success = monitor.refresh_job(job_id)
```

**Parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `job_id` | `str` | Yes | Job identifier to refresh |

**Returns:** `bool`

- `True` if job was refreshed successfully
- `False` if job not found or refresh failed

---

### JobInfo Object

Each job is represented by a `JobInfo` object with the following attributes:

| Attribute | Type | Description |
|-----------|------|-------------|
| `job_id` | `str` | Unique job identifier |
| `job_dir` | `Path` | Job directory path |
| `status` | `JobStatus` | Current status (RUNNING, COMPLETED, ERROR, UNKNOWN) |
| `progress` | `float` | Progress percentage (0.0 - 100.0) |
| `current_step` | `str` | Name of current processing step |
| `steps_completed` | `List[str]` | List of completed step names |
| `elapsed_time` | `float` | Elapsed time in seconds |
| `start_time` | `datetime` | Job start timestamp |
| `end_time` | `Optional[datetime]` | Job end timestamp (None if running) |

**JobStatus Values:**

- `RUNNING` - Job is currently executing
- `COMPLETED` - Job finished successfully
- `ERROR` - Job encountered an error
- `UNKNOWN` - Status cannot be determined

---

### Example 17: Real-time Progress Tracking

```python
from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path
import time

# Start a system preparation (runs in background)
from gatewizard.core.builder import Builder

builder = Builder()
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1"
)

if success:
    print(f"✓ Job started: {job_dir}")

    # Monitor progress in real-time
    monitor = JobMonitor(working_directory=Path("./systems"))

    while True:
        monitor.scan_for_jobs(force=True)
        active_jobs = monitor.get_active_jobs()

        if not active_jobs:
            # Job completed
            completed = monitor.get_completed_jobs()
            if completed:
                for job_id, job_info in completed.items():
                    if str(job_dir) in job_id:
                        print(f"\n✓ Job completed: {job_info.status.value}")
                        print(f"  Total time: {job_info.elapsed_time:.1f}s")
                        break
            break

        # Show progress
        for job_id, job_info in active_jobs.items():
            if str(job_dir) in job_id:
                print(f"\rProgress: {job_info.progress:.1f}% - {job_info.current_step}", end="", flush=True)

        time.sleep(2)  # Check every 2 seconds
```

---

### Example 18: Monitoring Batch Job Management

```python
from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path

# Monitor all jobs in a directory
monitor = JobMonitor(working_directory=Path("./systems"))
monitor.scan_for_jobs(force=True)

# Check active jobs
active_jobs = monitor.get_active_jobs()
print(f"\n{'='*60}")
print(f"ACTIVE JOBS: {len(active_jobs)}")
print(f"{'='*60}")

for job_id, job_info in active_jobs.items():
    print(f"\n📁 {job_info.job_dir.name}")
    print(f"   Status: {job_info.status.value}")
    print(f"   Progress: {job_info.progress:.1f}%")
    print(f"   Current: {job_info.current_step}")
    print(f"   Runtime: {job_info.elapsed_time:.0f}s")

    if job_info.steps_completed:
        print(f"   Completed steps:")
        for step in job_info.steps_completed:
            print(f"     ✓ {step}")

# Check completed jobs
completed_jobs = monitor.get_completed_jobs()
print(f"\n{'='*60}")
print(f"COMPLETED JOBS: {len(completed_jobs)}")
print(f"{'='*60}")

for job_id, job_info in completed_jobs.items():
    status_icon = "✓" if job_info.status.value == "completed" else "✗"
    print(f"{status_icon} {job_info.job_dir.name}: {job_info.status.value} ({job_info.elapsed_time:.0f}s)")
```

---

### Monitoring Tips

**Real-time monitoring:**

- Use `scan_for_jobs(force=True)` to bypass the 3-second scan interval
- Check status every 1-2 seconds for responsive progress updates
- Use `flush=True` in print statements for immediate display

**Long-running jobs:**

- Check status periodically (every 30-60 seconds)
- Log progress to file for later review
- Use `elapsed_time` to estimate remaining time

**Multiple jobs:**

- Monitor entire working directory for batch processing
- Filter jobs by status for targeted monitoring
- Use job_dir.name for user-friendly display

**Debugging:**

- Check `job_dir / 'logs/preparation.log'` for detailed output
- Check `job_dir / 'logs/parametrization.log'` if parametrization fails
- Examine `status.json` for detailed job state

**Manual log inspection:**
```bash
# Follow preparation log in real-time
tail -f systems/popc_membrane/logs/preparation.log

# Check for errors
grep -i error systems/popc_membrane/logs/*.log

# View job status
cat systems/popc_membrane/status.json
```

---

## Tips and Best Practices

### Lipid Selection

**For simple systems:**

- Use 100% POPC for both leaflets (symmetric membrane)
- Well-characterized, stable, and widely used
- Good starting point for method development

**For realistic membranes:**

- Include cholesterol (10-30%) for membrane stability
- Use PE lipids (POPE) for inner leaflet (more physiological)
- Add PS lipids (POPS) for negative charge (cytoplasmic side)

**For specialized studies:**

- Use asymmetric compositions (different upper/lower)
- Match experimental lipid compositions when available
- Consider protein's native membrane environment

### Force Field Selection

**Default recommendation (most stable):**

- Water: tip3p
- Protein: ff14SB
- Lipid: lipid21

**For latest parameters:**

- Water: opc
- Protein: ff19SB
- Lipid: lipid21

### Salt Concentration

**Physiological (default):**

- 0.15 M (150 mM) NaCl or KCl
- Mimics intracellular or extracellular conditions

**High salt studies:**

- 0.5-1.0 M for ionic strength effects
- May require longer equilibration

**Charge neutralization only:**

- Set `salt_concentration=0.0`
- System will still be neutralized with minimum ions

### Water Layer Thickness

**Default (26 Å):**

- Sufficient for most membrane proteins
- ~3-4 water layers above/below membrane

**Boundary distance (12 Å default):**

- Set `dist` to control the minimum distance from the solute extents to box boundaries
- packmol-memgen may choose a larger real distance because it uses the worst-case scenario

**Large proteins or complexes:**

- Increase to 20-25 Å
- Ensures adequate solvation

**Minimize system size:**

- Decrease to 15 Å (minimum recommended)
- Faster simulations but ensure full solvation

### Protein Orientation

**Pre-oriented proteins:**

- Use [OPM database](https://opm.phar.umich.edu/)
- Protein already positioned in membrane plane
- Set `preoriented=True` (default)

**Non-oriented proteins:**

- Let packmol-memgen use MEMEMBED for orientation
- Set `preoriented=False`
- Adds extra time to preparation

### Troubleshooting

**"PDB file not found":**

- Check absolute/relative paths
- Verify file permissions
- Ensure PDB file is valid format

**"Packing failed":**

- Protein may be too large for default settings
- Try increasing `dist_wat` parameter
- Check protein orientation (OPM database)
- Verify protein is pre-oriented if using `preoriented=True`

**"Parametrization failed":**

- Check logs/parametrization.log for details
- May have non-standard residues
- Verify protonation states (use Propka workflow)
- Check for missing atoms or unusual residue names

**"Invalid lipid ratios":**

- Ensure ratios match number of lipids
- Format: "ratio1:ratio2//ratio3:ratio4"
- Ratios don't need to sum to 1.0 (auto-normalized)
- Example: "7:3//5:5" for 2 lipids per leaflet

**"Force field compatibility error":**

- Use ForceFieldManager to check combinations
- Not all force fields are compatible
- Use recommended combinations for stability

### Performance Tips

**Faster preparation:**

- Use symmetric membranes (same upper/lower)
- Fewer lipid types (1-2 per leaflet)
- Pre-oriented proteins
- Smaller water layers (15-20 Å)

**More realistic systems:**

- Asymmetric membranes (different leaflets)
- Multiple lipid types (3-4 per leaflet)
- Include cholesterol
- Larger water layers (20-25 Å)
- Use Propka workflow for correct protonation

---

## Ligand Parametrization

Module for detecting, extracting, and parametrizing non-standard (ligand) residues found in PDB files.
Uses the AMBER/GAFF workflow (antechamber → parmchk2 → tleap) to generate force field parameters
that integrate with the Builder's packmol-memgen and final tleap steps.
Supports both GAFF and GAFF2 atom type sets, with recommended pairings
per the AMBER manual: **gaff/bcc** and **gaff2/abcg2**.

### Import

```python
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    extract_ligand_pdb,
    parametrize_ligand,
    parametrize_all_ligands,
    get_ligand_2d_image,
    get_ligand_2d_image_from_pdb_lines,
    build_ligand_param_args,
    build_tleap_ligand_lines,
    LigandInfo,
    LigandParametrizationError,
    STANDARD_RESIDUES,
    CHARGE_METHODS,
    ATOM_TYPES,
    DEFAULT_CHARGE_METHOD,
    DEFAULT_ATOM_TYPE,
    RECOMMENDED_COMBOS,
    NON_RECOMMENDED_COMBOS,
)
```

### Class: LigandInfo

Information container for a detected ligand residue.

| Attribute | Type | Description |
|-----------|------|-------------|
| `name` | `str` | 3-letter residue name (e.g., `"AAA"`) |
| `chain` | `str` | Chain identifier |
| `res_id` | `int` | Residue sequence number |
| `num_atoms` | `int` | Number of atoms in the ligand |
| `elements` | `Dict[str, int]` | Element counts (e.g., `{'C': 9, 'H': 10, 'O': 1}`) |
| `pdb_lines` | `List[str]` | Raw HETATM lines from the PDB file |
| `formula` | `str` (property) | Molecular formula string (e.g., `"C9H10O"`) |

**Methods:**

| Method | Returns | Description |
|--------|---------|-------------|
| `to_dict()` | `Dict[str, Any]` | Convert to serializable dictionary |

### Constants

| Constant | Type | Description |
|----------|------|-------------|
| `STANDARD_RESIDUES` | `set` | Residue names excluded from detection (amino acids, water, ions, lipids, capping groups) |
| `CHARGE_METHODS` | `dict` | Antechamber charge methods: `'bcc'` (AM1-BCC), `'abcg2'` (ABCG2), `'gas'` (Gasteiger), `'mul'` (Mulliken), `'cm2'` (CM2), `'rc'` (Read-in), `'resp'` (RESP — requires Gaussian), `'esp'` (ESP — requires Gaussian) |
| `ATOM_TYPES` | `dict` | Antechamber atom type sets: `'gaff2'` → `'GAFF2'`, `'gaff'` → `'GAFF'` |
| `DEFAULT_CHARGE_METHOD` | `str` | `'bcc'` |
| `DEFAULT_ATOM_TYPE` | `str` | `'gaff2'` |
| `RECOMMENDED_COMBOS` | `set` | Recommended (atom_type, charge_method) pairings per AMBER manual: `{('gaff', 'bcc'), ('gaff2', 'abcg2')}` |
| `NON_RECOMMENDED_COMBOS` | `set` | Non-recommended pairings (warning shown in GUI): `{('gaff2', 'bcc'), ('gaff', 'abcg2')}` |

!!! warning "Recommended Atom Type / Charge Method Pairings"
    Per the AMBER manual, the efficient charge models `bcc` (AM1-BCC) and `abcg2` (ABCG2)
    should be paired with specific atom type sets:

    | Atom Type | Charge Method | Status |
    |-----------|---------------|--------|
    | `gaff` | `bcc` | **Recommended** ✓ |
    | `gaff2` | `abcg2` | **Recommended** ✓ |
    | `gaff2` | `bcc` | Not recommended ✗ |
    | `gaff` | `abcg2` | Not recommended ✗ |

    The GUI shows an amber warning label when a non-recommended combination is selected,
    and displays a confirmation dialog before proceeding with parametrization.

!!! info "External QM Software Requirements"
    The `resp` and `esp` charge methods require **Gaussian** (external quantum-mechanics
    software, not included in AmberTools). All other methods (`bcc`, `abcg2`, `gas`,
    `mul`, `cm2`) use **sqm**, which is bundled with AmberTools.

### Function: detect_ligands()

Detect non-standard (ligand) residues in a PDB file by scanning HETATM records.

```python
detect_ligands(pdb_file: str) -> List[LigandInfo]
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_file` | `str` | Path to the PDB file to analyze |

**Returns:** `List[LigandInfo]` — one entry per unique ligand residue name

**Raises:** `LigandParametrizationError` if the file cannot be read

### Example 19: Detect Ligands

```python
"""
Builder Example 19: Detect ligands in a PDB file

Demonstrates how to detect non-standard residues (ligands)
in a PDB file using the ligand parametrization tools.
"""

from gatewizard.tools.ligand_parametrization import detect_ligands

# Detect ligands in a PDB file with two ligands (AAA and BBB)
pdb_file = "tests/2MVJ_2ligs.pdb"
ligands = detect_ligands(pdb_file)

print(f"Detected {len(ligands)} ligand(s) in {pdb_file}:")
print(f"{'='*60}")

for lig in ligands:
    print(f"\nLigand: {lig.name}")
    print(f"  Chain: {lig.chain}")
    print(f"  Residue ID: {lig.res_id}")
    print(f"  Number of atoms: {lig.num_atoms}")
    print(f"  Molecular formula: {lig.formula}")
    print(f"  Elements: {lig.elements}")
    print(f"  As dict: {lig.to_dict()}")
```

---

### Function: extract_ligand_pdb()

Extract a single ligand from a PDB file into its own file.

```python
extract_ligand_pdb(pdb_file: str, ligand_name: str, output_dir: str) -> str
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_file` | `str` | Path to the source PDB file |
| `ligand_name` | `str` | 3-letter residue name of the ligand |
| `output_dir` | `str` | Directory to write the extracted PDB file |

**Returns:** `str` — path to the extracted PDB file (`output_dir/LIGAND.pdb`)

**Raises:** `LigandParametrizationError` if the ligand is not found or extraction fails

### Example 20: Extract a Ligand

```python
"""
Builder Example 20: Extract a ligand from a PDB file

Demonstrates how to extract a specific ligand into its own PDB file
for individual parametrization.
"""

from pathlib import Path
from gatewizard.tools.ligand_parametrization import detect_ligands, extract_ligand_pdb

pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_extraction"

# First detect ligands
ligands = detect_ligands(pdb_file)
print(f"Detected ligands: {[l.name for l in ligands]}")

# Extract each ligand to its own subdirectory
for lig in ligands:
    lig_dir = str(Path(output_dir) / lig.name)
    extracted_pdb = extract_ligand_pdb(pdb_file, lig.name, lig_dir)

    # Verify extraction
    with open(extracted_pdb) as f:
        lines = [l for l in f if l.startswith("HETATM")]

    print(f"\nExtracted {lig.name}:")
    print(f"  Output: {extracted_pdb}")
    print(f"  Atoms: {len(lines)}")
    print(f"  First line: {lines[0].strip()[:60]}...")

print(f"\nAll extracted ligands saved in: {output_dir}")
```

---

### Function: parametrize_ligand()

Parametrize a single ligand using the AMBER/GAFF workflow (antechamber → parmchk2 → tleap).

```python
parametrize_ligand(
    ligand_pdb: str,
    ligand_name: str,
    output_dir: str,
    charge: int = 0,
    charge_method: str = 'bcc',
    multiplicity: int = 1,
    atom_type: str = 'gaff2',
) -> Dict[str, str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ligand_pdb` | `str` | — | Path to the ligand PDB file |
| `ligand_name` | `str` | — | 3-letter residue name (must match PDB) |
| `output_dir` | `str` | — | Directory for output files |
| `charge` | `int` | `0` | Net charge of the ligand |
| `charge_method` | `str` | `'bcc'` | Charge method (see `CHARGE_METHODS`) |
| `multiplicity` | `int` | `1` | Spin multiplicity |
| `atom_type` | `str` | `'gaff2'` | Atom type set — `'gaff2'` or `'gaff'` (see `ATOM_TYPES`) |

**Returns:** Dictionary with paths to generated files:

| Key | Description |
|-----|-------------|
| `'mol2'` | Typed MOL2 file (GAFF/GAFF2 atom types + charges) |
| `'frcmod'` | Force field modification file (missing parameters) |
| `'lib'` | Residue library file for tleap |
| `'prmtop'` | AMBER topology file |
| `'inpcrd'` | AMBER coordinate file |

**Raises:** `LigandParametrizationError` if any step fails

!!! warning "Critical: Variable Name Must Match Residue Name"
    The tleap variable name **must** match the residue name in the PDB.
    For example, `AAA = loadmol2 AAA.mol2` — using a generic name like
    `mol = loadmol2 AAA.mol2` will cause `saveoff` to store the unit under the
    wrong name and tleap will fail to recognize the residue when loading the full system PDB.

---

### Function: parametrize_all_ligands()

Detect and parametrize all ligands in a PDB file in one call.

```python
parametrize_all_ligands(
    pdb_file: str,
    output_dir: str,
    charges: Optional[Dict[str, int]] = None,
    charge_method: str = 'bcc',
    atom_type: str = 'gaff2',
) -> Dict[str, Dict[str, str]]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_file` | `str` | — | PDB file containing ligands |
| `output_dir` | `str` | — | Base directory (each ligand gets a subdirectory) |
| `charges` | `Dict[str, int]` or `None` | `None` | Map of ligand name → net charge (default: 0 for all) |
| `charge_method` | `str` | `'bcc'` | Charge method for antechamber |
| `atom_type` | `str` | `'gaff2'` | Atom type set — `'gaff2'` or `'gaff'` (see `ATOM_TYPES`) |

**Returns:** `Dict[str, Dict[str, str]]` — maps ligand names to their file paths (same structure as `parametrize_ligand()`)

### Example 22: Parametrize All Ligands

```python
"""
Builder Example 22: Parametrize all ligands in a PDB file

Demonstrates the full ligand parametrization workflow:
1. Detect ligands from PDB
2. Extract each ligand
3. Run antechamber + parmchk2 + tleap for each
4. Collect .frcmod and .lib files

NOTE: This example requires AmberTools (antechamber, parmchk2, tleap)
to be installed and accessible in the PATH.
"""

from pathlib import Path
from gatewizard.tools.ligand_parametrization import (
    parametrize_all_ligands,
    build_ligand_param_args,
    build_tleap_ligand_lines,
)

pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_params"

# Set charges for each ligand (default is 0 if not specified)
charges = {
    'AAA': 0,
    'BBB': 0,
}

print(f"Output directory: {output_dir}")

# Parametrize all ligands
results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=output_dir,
    charges=charges,
    charge_method='bcc',  # AM1-BCC charges
    atom_type='gaff2',    # GAFF2 atom types (recommended with abcg2; bcc also works)
)

print(f"\nParametrized {len(results)} ligand(s):")
for name, files in results.items():
    print(f"\n  {name}:")
    for file_type, file_path in files.items():
        exists = Path(file_path).exists()
        print(f"    {file_type}: {file_path} ({'exists' if exists else 'MISSING'})")

# Now these can be passed to Builder.prepare_system()
print("\n\npackmol-memgen arguments:")
pmm_args = build_ligand_param_args(results)
print(f"  {' '.join(pmm_args)}")

print("\ntleap ligand lines:")
tleap_lines = build_tleap_ligand_lines(results)
print(tleap_lines)

print(f"\nAll parameter files saved in: {output_dir}")
```

---

### Function: get_ligand_2d_image()

Generate a publication-quality 2D molecular structure image using RDKit.

```python
get_ligand_2d_image(
    ligand_pdb_or_mol2: str,
    output_image: str,
    width: int = 400,
    height: int = 300,
    *,
    remove_nonpolar_h: bool = True,
    remove_all_h: bool = False,
    dpi: int = 150,
    bond_line_width: float = 2.5,
    atom_label_font_size: int = 0,
    background_color: tuple = (0.11, 0.11, 0.11, 1.0),
    padding: float = 0.15,
    kekulize: bool = True,
    wedge_bonds: bool = True,
    atom_palette: dict | None = None,
    highlight_atoms: list[int] | None = None,
    highlight_color: tuple = (1.0, 0.8, 0.0, 0.3),
    transparent_background: bool = False,
) -> Optional[str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `ligand_pdb_or_mol2` | `str` | — | Path to PDB or MOL2 file |
| `output_image` | `str` | — | Path for the output PNG |
| `width` | `int` | `400` | Image width in pixels (before DPI scaling) |
| `height` | `int` | `300` | Image height in pixels (before DPI scaling) |
| `remove_nonpolar_h` | `bool` | `True` | Remove non-polar hydrogens (C-H) for cleaner figures while keeping polar H (O-H, N-H, etc.) |
| `remove_all_h` | `bool` | `False` | Remove **all** explicit hydrogens. Overrides `remove_nonpolar_h` |
| `dpi` | `int` | `150` | DPI multiplier. Use `300` for print-quality images |
| `bond_line_width` | `float` | `2.5` | Thickness of bond lines |
| `atom_label_font_size` | `int` | `0` | Font size for atom labels (0 = auto) |
| `background_color` | `tuple` | `(0.11, 0.11, 0.11, 1.0)` | RGBA background colour (0–1 range) |
| `padding` | `float` | `0.15` | Fractional padding around the molecule (0–1) |
| `kekulize` | `bool` | `True` | Draw aromatic bonds as alternating single/double |
| `wedge_bonds` | `bool` | `True` | Draw stereo wedge/dash bonds |
| `atom_palette` | `dict` | `None` | Custom `{atomic_number: (r, g, b)}` colour map. Uses built-in dark palette when `None` |
| `highlight_atoms` | `list[int]` | `None` | 0-based atom indices to highlight |
| `highlight_color` | `tuple` | `(1.0, 0.8, 0.0, 0.3)` | RGBA colour for highlights |
| `transparent_background` | `bool` | `False` | Produce a transparent PNG |

**Returns:** Path to the generated PNG image, or `None` if generation failed

**Built-in palettes:**

| Palette | Import | Best for |
|---------|--------|----------|
| `_DEFAULT_DARK_PALETTE` | (used automatically) | Dark backgrounds |
| `LIGHT_PALETTE` | `from gatewizard.tools.ligand_parametrization import LIGHT_PALETTE` | White / light backgrounds |

### Function: get_ligand_2d_image_from_pdb_lines()

Generate a 2D image directly from PDB HETATM lines (no file required).
All keyword arguments are forwarded to `get_ligand_2d_image`.

```python
get_ligand_2d_image_from_pdb_lines(
    pdb_lines: List[str],
    output_image: str,
    width: int = 400,
    height: int = 300,
    **kwargs,
) -> Optional[str]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_lines` | `List[str]` | — | HETATM lines for the ligand |
| `output_image` | `str` | — | Path for the output PNG |
| `width` | `int` | `400` | Image width |
| `height` | `int` | `300` | Image height |
| `**kwargs` | | | All keyword options from `get_ligand_2d_image` |

**Returns:** Path to the generated PNG image, or `None` if failed

### Example 24: Generate 2D Structure Images

```python
"""
Builder Example 24: Generate 2D structure images of ligands

Demonstrates generating 2D molecular structure images from PDB files
using RDKit, including options for hydrogen removal, DPI control,
custom colour palettes, and transparent backgrounds.
"""

from pathlib import Path
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    get_ligand_2d_image_from_pdb_lines,
    get_ligand_2d_image,
    extract_ligand_pdb,
    LIGHT_PALETTE,
)

pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_images"
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Detect ligands
ligands = detect_ligands(pdb_file)
print(f"Detected {len(ligands)} ligands\n")

for lig in ligands:
    # --- Style 1: Default (dark background, non-polar H removed) --------
    img = str(Path(output_dir) / f"{lig.name}_default.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=400, height=300,
        remove_nonpolar_h=True,           # cleaner look (default)
    )
    if result:
        print(f"{lig.name} default : {Path(result).stat().st_size:>6} bytes")

    # --- Style 2: All hydrogens removed ----------------------------------
    img = str(Path(output_dir) / f"{lig.name}_no_h.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=400, height=300,
        remove_all_h=True,                # skeleton only
    )
    if result:
        print(f"{lig.name} no-H    : {Path(result).stat().st_size:>6} bytes")

    # --- Style 3: High-DPI for publication (white background) ------------
    img = str(Path(output_dir) / f"{lig.name}_hires.png")
    lig_dir = str(Path(output_dir) / lig.name)
    extracted_pdb = extract_ligand_pdb(pdb_file, lig.name, lig_dir)
    result = get_ligand_2d_image(
        extracted_pdb, img,
        width=800, height=600,
        dpi=300,                           # high DPI
        remove_all_h=True,
        background_color=(1, 1, 1, 1),     # white
        atom_palette=LIGHT_PALETTE,        # colours for light background
        bond_line_width=1.5,
    )
    if result:
        print(f"{lig.name} hi-res  : {Path(result).stat().st_size:>6} bytes")

    # --- Style 4: Transparent background ---------------------------------
    img = str(Path(output_dir) / f"{lig.name}_transparent.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=400, height=300,
        remove_nonpolar_h=True,
        transparent_background=True,
    )
    if result:
        print(f"{lig.name} transp. : {Path(result).stat().st_size:>6} bytes")

    # --- Style 5: All hydrogens visible, thicker bonds -------------------
    img = str(Path(output_dir) / f"{lig.name}_all_h.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=500, height=400,
        remove_nonpolar_h=False,
        remove_all_h=False,
        bond_line_width=3.5,
        padding=0.2,
    )
    if result:
        print(f"{lig.name} all-H   : {Path(result).stat().st_size:>6} bytes")

    print()

print(f"All images saved in: {output_dir}")
```

---

### Function: build_ligand_param_args()

Build `--ligand_param` command-line arguments for packmol-memgen.

```python
build_ligand_param_args(
    ligand_files: Dict[str, Dict[str, str]]
) -> List[str]
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `ligand_files` | `Dict[str, Dict[str, str]]` | Ligand name → file paths (as returned by `parametrize_ligand`) |

**Returns:** List of CLI arguments — each ligand produces `['--ligand_param', 'path.frcmod:path.lib']`

!!! info "One Flag Per Ligand"
    packmol-memgen requires a **separate** `--ligand_param` flag per ligand.
    Combining multiple ligands into a single flag will not work.

### Function: build_tleap_ligand_lines()

Build tleap input lines to load GAFF/GAFF2 and ligand parameters. Insert these **before** `loadPDB`.

```python
build_tleap_ligand_lines(
    ligand_files: Dict[str, Dict[str, str]],
    atom_type: str = 'gaff2',
) -> str
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `ligand_files` | `Dict[str, Dict[str, str]]` | Ligand name → file paths |
| `atom_type` | `str` | Atom type set — `'gaff2'` or `'gaff'` (default: `'gaff2'`). Determines the `leaprc` source line |

**Returns:** Multi-line string with tleap commands (`source leaprc.gaff2` or `source leaprc.gaff`, `loadamberparams`, `loadoff`)

### Example 21: Build packmol-memgen and tleap Arguments

```python
"""
Builder Example 21: Build packmol-memgen and tleap commands with ligand parameters

Demonstrates how to construct the packmol-memgen --ligand_param arguments
and tleap input lines for systems with ligands.
"""

from gatewizard.tools.ligand_parametrization import (
    build_ligand_param_args,
    build_tleap_ligand_lines,
)

# Simulated parametrization results (paths to .frcmod and .lib files)
ligand_files = {
    'AAA': {
        'frcmod': 'ligand_params/AAA/AAA.frcmod',
        'lib': 'ligand_params/AAA/AAA.lib',
        'mol2': 'ligand_params/AAA/AAA.mol2',
    },
    'BBB': {
        'frcmod': 'ligand_params/BBB/BBB.frcmod',
        'lib': 'ligand_params/BBB/BBB.lib',
        'mol2': 'ligand_params/BBB/BBB.mol2',
    },
}

# Build packmol-memgen arguments
# Each ligand gets its own --ligand_param flag (CANNOT combine in one flag)
pmm_args = build_ligand_param_args(ligand_files)
print("packmol-memgen arguments:")
for i in range(0, len(pmm_args), 2):
    print(f"  {pmm_args[i]} {pmm_args[i+1]}")

print()

# Build tleap input lines
# These go BEFORE loadPDB in the tleap input
tleap_lines = build_tleap_ligand_lines(ligand_files)
print("tleap input lines:")
print(tleap_lines)

# With GAFF atom types instead of GAFF2
tleap_lines_gaff = build_tleap_ligand_lines(ligand_files, atom_type='gaff')
print("\ntleap input lines (GAFF):")
print(tleap_lines_gaff)
```

---

### Builder Integration: `ligand_params` Config Key

When ligands are parametrized (either through the GUI or programmatically), the Builder uses the `ligand_params` configuration key to integrate them into the build process.

**Config Structure:**

```python
config['ligand_params'] = {
    'AAA': {
        'frcmod': '/path/to/AAA/AAA.frcmod',
        'lib': '/path/to/AAA/AAA.lib',
    },
    'BBB': {
        'frcmod': '/path/to/BBB/BBB.frcmod',
        'lib': '/path/to/BBB/BBB.lib',
    },
}
```

**What happens during build:**

1. **packmol-memgen** receives `--ligand_param frcmod:lib` for each ligand, plus `--gaff2` (or `--gaff` depending on the atom type used)
2. **tleap parametrization** loads `source leaprc.gaff2` (or `source leaprc.gaff`), then `loadamberparams` and `loadoff` for each ligand before `loadPDB`
3. The atom type is automatically extracted from the parametrized ligand results, so matching `leaprc` is always used

### Example 23: Full Membrane System with Ligand Parametrization

```python
"""
Builder Example 23: Full membrane system with ligand parametrization

Demonstrates setting up a complete membrane system that includes
non-standard ligands. The ligand .frcmod/.lib files are passed
to both packmol-memgen and the final tleap parametrization.

NOTE: This example requires AmberTools and packmol-memgen.
"""

from gatewizard.core.builder import Builder
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    parametrize_all_ligands,
)

# Create builder
builder = Builder()

pdb_file = "tests/2MVJ_2ligs.pdb"
working_dir = "./systems"

# Step 1: Detect and parametrize ligands
print("Step 1: Detecting ligands...")
ligands = detect_ligands(pdb_file)
for lig in ligands:
    print(f"  Found: {lig.name} ({lig.num_atoms} atoms, {lig.formula})")

print("\nStep 2: Parametrizing ligands...")
ligand_results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=f"{working_dir}/ligand_params",
    charges={'AAA': 0, 'BBB': 0},
    charge_method='bcc',
    atom_type='gaff2',    # GAFF2 atom types
)

print(f"  Parametrized: {list(ligand_results.keys())}")

# Step 3: Configure builder with ligand parameters
builder.set_configuration(
    water_model='tip3p',
    protein_ff='ff14SB',
    lipid_ff='lipid21',
    preoriented=True,
    parametrize=True,
    salt_concentration=0.15,
    dist=12,
    dist_wat=26,
    notprotonate=True,
    ligand_params=ligand_results,  # Pass parametrized ligand files
)

print("\nStep 3: Builder configured with ligand parameters")
print(f"  Ligands in config: {list(builder.config['ligand_params'].keys())}")

# Step 4: Prepare system
success, message, job_dir = builder.prepare_system(
    pdb_file=pdb_file,
    working_dir=working_dir,
    upper_lipids=['POPC'],
    lower_lipids=['POPC'],
    lipid_ratios='1.0//1.0',
)
print(f"\nResult: {message}")

# The builder will:
# 1. Add --ligand_param AAA/AAA.frcmod:AAA/AAA.lib to packmol-memgen
# 2. Add --ligand_param BBB/BBB.frcmod:BBB/BBB.lib to packmol-memgen
# 3. Add --gaff2 to packmol-memgen
# 4. Load ligand .frcmod and .lib in the tleap parametrization step
print("\nWorkflow complete. The builder will pass ligand params to:")
print("  - packmol-memgen (--ligand_param flags)")
print("  - tleap (loadamberparams/loadoff commands)")
```

### Example 26: Atom Type Selection and Recommended Pairings

```python
"""
Builder Example 26: Atom type selection and recommended pairings

Demonstrates how to choose between GAFF and GAFF2 atom types,
check recommended pairings per the AMBER manual, and use
the ABCG2 charge method with GAFF2.

NOTE: This example requires AmberTools (antechamber, parmchk2, tleap)
to be installed and accessible in the PATH.
"""

from gatewizard.tools.ligand_parametrization import (
    parametrize_all_ligands,
    build_tleap_ligand_lines,
    ATOM_TYPES,
    CHARGE_METHODS,
    DEFAULT_ATOM_TYPE,
    DEFAULT_CHARGE_METHOD,
    RECOMMENDED_COMBOS,
    NON_RECOMMENDED_COMBOS,
)

# ── Available options ────────────────────────────────────────────────
print("Available atom types:")
for key, label in ATOM_TYPES.items():
    default = " (default)" if key == DEFAULT_ATOM_TYPE else ""
    print(f"  {key}: {label}{default}")

print("\nAvailable charge methods:")
for key, label in CHARGE_METHODS.items():
    default = " (default)" if key == DEFAULT_CHARGE_METHOD else ""
    print(f"  {key}: {label}{default}")

# ── Recommended pairings ─────────────────────────────────────────────
print("\nRecommended pairings (AMBER manual):")
for at, cm in sorted(RECOMMENDED_COMBOS):
    print(f"  {at} + {cm}  ✓")

print("\nNon-recommended pairings (will show warning):")
for at, cm in sorted(NON_RECOMMENDED_COMBOS):
    print(f"  {at} + {cm}  ✗")

# ── Check a pairing before parametrizing ─────────────────────────────
atom_type = 'gaff2'
charge_method = 'abcg2'

if (atom_type, charge_method) in RECOMMENDED_COMBOS:
    print(f"\n{atom_type}/{charge_method} is a recommended pairing.")
elif (atom_type, charge_method) in NON_RECOMMENDED_COMBOS:
    print(f"\nWARNING: {atom_type}/{charge_method} is NOT recommended.")
else:
    print(f"\n{atom_type}/{charge_method} has no specific recommendation.")

# ── Parametrize with gaff2/abcg2 (recommended) ──────────────────────
pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_params_gaff2_abcg2"

results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=output_dir,
    charges={'AAA': 0, 'BBB': 0},
    charge_method='abcg2',   # ABCG2 charges
    atom_type='gaff2',       # GAFF2 atom types (recommended with abcg2)
)

print(f"\nParametrized {len(results)} ligand(s) with gaff2/abcg2:")
for name, files in results.items():
    print(f"  {name}: {files.get('frcmod', 'N/A')}")

# ── tleap lines reflect the chosen atom type ─────────────────────────
tleap_gaff2 = build_tleap_ligand_lines(results, atom_type='gaff2')
print(f"\ntleap lines (GAFF2):\n{tleap_gaff2}")

tleap_gaff = build_tleap_ligand_lines(results, atom_type='gaff')
print(f"\ntleap lines (GAFF):\n{tleap_gaff}")
```

### Ligand Parametrization Troubleshooting

**"Antechamber failed":**

- Check that AMBER/AmberTools is installed and in `$PATH`
- Verify ligand PDB has proper HETATM records
- Try a different charge method (e.g., `gas` is faster and may work when `bcc` fails)
- Check `logs/antechamber.log` for detailed error messages

**"No ligands detected":**

- Ligands must use HETATM records, not ATOM
- Residue names must NOT be in `STANDARD_RESIDUES` set
- Water (HOH, WAT), ions (NA, CL), and lipids (POPC, etc.) are excluded by design

**"tleap: unknown residue name":**

- The tleap variable name must match the 3-letter residue name exactly
- Ensure `.lib` file was generated with the correct residue name in `saveoff`
- Check that `source leaprc.gaff2` (or `source leaprc.gaff`) is loaded before `loadamberparams`

**"Non-recommended combination" warning (gaff2/bcc or gaff/abcg2):**

- The AMBER manual recommends `gaff/bcc` and `gaff2/abcg2` pairings
- Using `gaff2/bcc` or `gaff/abcg2` may produce suboptimal parameters
- The GUI shows an amber warning label and a confirmation dialog
- To silence the warning, switch to a recommended pairing
- Check `RECOMMENDED_COMBOS` and `NON_RECOMMENDED_COMBOS` constants for the full list

**"resp/esp: Gaussian not found":**

- The `resp` and `esp` charge methods require Gaussian (external QM software)
- Gaussian is **not** included in AmberTools
- Use `bcc` or `abcg2` instead (both use the built-in `sqm` engine)

**"RDKit not available":**

- Install with `conda install -c conda-forge rdkit`
- RDKit is a required dependency for ligand 2D visualization

---

## See Also

- [Preparation Module](preparation.md) - Protonation state analysis
- [Equilibration Module](equilibration.md) - MD equilibration protocols
- [User Guide - Builder Tab](../user-guide.md#builder-tab-system-building) - GUI workflow
- [Examples on GitHub](https://github.com/maurobedoya/gatewizard/tree/main/examples) - Complete workflow examples
