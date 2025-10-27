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
from gatewizard.core.system_builder import SystemBuilder
from gatewizard.tools.force_fields import ForceFieldManager
```

## Class: SystemBuilder

Main class for building membrane protein systems with complete control over lipid composition, force fields, and system parameters.

### Constructor

```python
SystemBuilder()
```

**Parameters:** None

**Returns:** `SystemBuilder` instance

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
| `dist_wat` | `17.5` | Water layer thickness in Å |
| `notprotonate` | `False` | Skip protonation (preserve residue names) |

### Example 1: Basic Configuration
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()
print(f"Default water model: {builder.config['water_model']}")
print(f"Default protein force field: {builder.config['protein_ff']}")
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
| `dist_wat` | `float` | `17.5` | Water layer thickness in Å |
| `notprotonate` | `bool` | `False` | Skip protonation during parametrization |
| `add_salt` | `bool` | `True` | Whether to add salt to system |

**Returns:** None

### Example 2: Custom Configuration
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

# Configure for specific system
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    salt_concentration=0.15,
    cation="Na+",
    anion="Cl-",
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

### Example 4: Available Protein Force Fields
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

### Example 5: Available Lipid Force Fields
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

### Example 6: Input Validation
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

# Validate inputs before preparation
valid, error_msg = builder.validate_system_inputs(
    pdb_file="protein_protonated_prepared.pdb",
    upper_lipids=["POPC", "POPE"],
    lower_lipids=["POPC", "POPE"],
    lipid_ratios="7:3//7:3",
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21"
)

if valid:
    print("✓ All inputs are valid, proceeding with preparation")
    # Now call prepare_system()
else:
    print(f"✗ Validation failed: {error_msg}")
    # Fix issues before proceeding
```

### Example 7: Force Field Validation
```python
from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()
valid, message = ff_manager.validate_combination("tip3p", "ff14SB", "lipid21")

if valid:
    print("✓ Force field combination is compatible")
else:
    print(f"✗ Incompatible: {message}")
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
- `SystemBuilderError` - If validation fails or system preparation fails

**Note:** System preparation runs in the background. The method returns immediately after launching the job. Monitor progress using the log files or the JobMonitor class (see monitoring example below).

### Example 8: Simple Symmetric Membrane
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

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

### Example 9: Monitoring Job Progress
```python
from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path
import time

# Create monitor for your working directory
monitor = JobMonitor(working_directory=Path("./systems"))

# Scan for jobs
monitor.scan_for_jobs(force=True)

# Get active jobs
active_jobs = monitor.get_active_jobs()

if active_jobs:
    print(f"✓ Found {len(active_jobs)} active job(s)")
    
    for job_id, job_info in active_jobs.items():
        print(f"\nJob: {job_info.job_dir.name}")
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
    print(f"\n✓ Found {len(completed_jobs)} completed job(s)")
    for job_id, job_info in completed_jobs.items():
        print(f"  - {job_info.job_dir.name}: {job_info.status.value}")
```

### Example 10: Asymmetric Membrane with Multiple Lipids
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

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

### Example 11: Complex Composition (Plasma Membrane Mimic)
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

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

### Example 12: Packing Only (No Parametrization)
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

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

### Example 13: Custom Salt Concentration
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

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

### Example 14: No Salt (Charge Neutralization Only)
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()

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

## Job Monitoring

Since system preparation runs in the background, GateWizard provides the `JobMonitor` class to track progress, check status, and manage running jobs.

### Class: JobMonitor

Monitor and track system preparation jobs.

### Example 15: JobMonitor class
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

### Example 16: Real-time Progress Tracking

```python
from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path
import time

# Start a system preparation (runs in background)
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()
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

### Example 17: Monitoring Batch Job Management

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

**Default (17.5 Å):**

- Sufficient for most membrane proteins
- ~3-4 water layers above/below membrane

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
- Smaller water layers (15-17.5 Å)

**More realistic systems:**

- Asymmetric membranes (different leaflets)
- Multiple lipid types (3-4 per leaflet)
- Include cholesterol
- Larger water layers (20-25 Å)
- Use Propka workflow for correct protonation

---

## See Also

- [Preparation Module](preparation.md) - Protonation state analysis
- [Equilibration Module](equilibration.md) - MD equilibration protocols
- [User Guide - Builder Tab](../user-guide.md#builder-tab-system-building) - GUI workflow
- [Examples on GitHub](https://github.com/maurobedoya/gatewizard/tree/main/examples) - Complete workflow examples
