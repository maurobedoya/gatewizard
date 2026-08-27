# Equilibration Module

Module for setting up NAMD, OpenMM, GROMACS, and Amber equilibration protocols for membrane protein systems. Generates configuration files, restraint files, and run scripts for multi-stage equilibration simulations using AMBER force fields.

The three engine managers share a similar user-facing API, with engine-specific file formats and run scripts.

Systems without a protein (for example a packed bilayer) omit `protein_backbone` / `protein_sidechain` restraints and protein COM/rotation. Detection uses the system PDB or GRO (`pdb_has_protein`). Lipid restraints are unchanged.

## Import

```python
from gatewizard.tools.equilibration import (
    NAMDEquilibrationManager,
    OpenMMEquilibrationManager,
    GROMACSEquilibrationManager,
    AmberEquilibrationManager,
)
```

## Remote cluster helpers (Slurm-first)

Programmatic helpers for probing modules, rendering batch scripts, and recording remote `execution` metadata live in `gatewizard.utils.cluster` (used by the GUI). Typical flow: generate inputs locally (produces `run_equilibration.sh` for local runs and `run_equilibration_cluster.sh` for Slurm), then wrap the cluster runner with `render_batch_script(...)` (`run_command` defaults to `bash run_equilibration_cluster.sh`; workdir strategy e.g. `scratch_job_id`).

```python
from gatewizard.utils.cluster import BatchScriptRequest, render_batch_script

script = render_batch_script(
    BatchScriptRequest(
        job_name="eq_popc",
        cpus=8,
        gpus=1,
        modules=["md/namd/3.0b6+cuda"],
        workdir_strategy="scratch_job_id",
        scratch_root="$SCRATCH_DIR",
    )
)
```

See the user guide section **Running on HPC Clusters** for the GUI workflow (Settings → Clusters, Equilibration Run target).

## Class: NAMDEquilibrationManager

Manager for NAMD equilibration simulations with support for multi-stage protocols, flexible restraints, and ensemble control.

### Constructor

```python
NAMDEquilibrationManager(working_dir: Path, namd_executable: str = "namd3")
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `working_dir` | `Path` | Required | Working directory containing system files |
| `namd_executable` | `str` | `"namd3"` | NAMD executable name or path |

**Returns:** `NAMDEquilibrationManager` instance

---

## Quick Start

### Example 1: Automatic File Detection (Simplest)
```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to folder with system files
system_folder = Path("popc_membrane")

# Define equilibration stages
stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 310.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    }
]

# Setup with automatic file detection (no system_files needed!)
# scheme_type is auto-detected from the 'ensemble' field in stages
manager = NAMDEquilibrationManager(system_folder)
result = manager.setup_namd_equilibration(
    stage_params_list=stages,
    output_name="equilibration_example_01"
)

# Note that output folder is created inside the system_folder

print(f"Setup complete: {result['namd_dir']}")
# Run with: cd {result['namd_dir']} && ./run_equilibration.sh
```

### Example 2: Explicit File Paths (Alternative)
```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to folder with system files
system_folder = Path("popc_membrane")

# Explicitly define system files (if auto-detection doesn't work)
system_files = {
    'prmtop': str(system_folder / 'system.prmtop'),
    'inpcrd': str(system_folder / 'system.inpcrd'),
    'pdb': str(system_folder / 'system.pdb'),
    'bilayer_pdb': str(system_folder / 'bilayer_protein_protonated_prepared_lipid.pdb')
}

# Define equilibration stages
stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 310.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    }
]

# Setup with explicit file paths
# scheme_type is auto-detected from the 'ensemble' field in stages
manager = NAMDEquilibrationManager(system_folder)
result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    output_name="equilibration_example_02"
)

print(f"Setup complete: {result['namd_dir']}")
# Run with: cd {result['namd_dir']} && ./run_equilibration.sh
```

---

## Helper Methods

### Method: find_system_files()

Automatically detect system files in the working directory. Useful for verifying which files will be used before running setup.

```python
find_system_files() -> Optional[Dict[str, str]]
```

**Parameters:** None

**Returns:** `Optional[Dict[str, str]]`

- Dictionary with detected file paths, or `None` if required files not found
- Keys: `'prmtop'`, `'inpcrd'`, `'pdb'`, `'bilayer_pdb'`

**Detection Priority:**

1. **Topology:** `*.prmtop` files
2. **Coordinates:** `*.inpcrd` → `*.crd` → `*.rst`
3. **System PDB:** `system.pdb` → any non-bilayer `*.pdb`
4. **Bilayer PDB:** `bilayer*_lipid.pdb` (with CRYST1) → `bilayer_*.pdb` → `*_bilayer.pdb`

### Example 3:
```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

manager = NAMDEquilibrationManager(Path("popc_membrane"))

# Check which files will be used
system_files = manager.find_system_files()
if system_files:
    print("Detected files:")
    for key, path in system_files.items():
        print(f"  {key}: {Path(path).name}")
    
    # Now run setup with auto-detection
    #result = manager.setup_namd_equilibration(
    #    stage_params_list=stages
    #)
else:
    print("Required files not found - please check working directory")

# Output:
# Detected files:
#   prmtop: system.prmtop
#   inpcrd: system.inpcrd
#   pdb: system.pdb
#   bilayer_pdb: bilayer_protein_protonated_prepared_lipid.pdb

```

---

## Core Method

### Method: setup_namd_equilibration()

Complete equilibration setup with automatic file generation. This method:

1. Auto-detects system files (if not provided)
2. Creates output directory structure
3. Copies system files (prmtop, inpcrd, pdb) to NAMD directory
4. Generates NAMD configuration files for each stage
5. Creates restraint files based on constraints
6. Generates executable run script
7. Creates protocol summary JSON

```python
setup_namd_equilibration(
    system_files: Optional[Dict[str, str]] = None,
    stage_params_list: List[Dict[str, Any]] = None,
    output_name: str = "equilibration",
    scheme_type: str = "NPT",
    namd_executable: str = "namd3"
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `system_files` | `Optional[Dict[str, str]]` | `None` | System file paths (see below). If `None`, auto-detects files in working_dir |
| `stage_params_list` | `List[Dict[str, Any]]` | Required | List of stage parameter dictionaries (see below). Each stage must include `ensemble` field |
| `output_name` | `str` | `"equilibration"` | Base name for output folder |
| `scheme_type` | `Optional[str]` | `None` | (Optional) Equilibration scheme (NVT, NPT, NPAT, NPgT). If `None`, auto-detected from first stage's `ensemble` field |
| `namd_executable` | `str` | `"namd3"` | NAMD executable name/path |

**Returns:** `Dict[str, Any]` with keys:

- `namd_dir` (`Path`): NAMD output directory path
- `config_files` (`List[Path]`): Generated .conf files
- `run_script` (`Path`): Executable run script path
- `summary_file` (`Path`): Protocol summary JSON path

### System Files Dictionary

**Auto-Detection (Recommended):**

If `system_files=None` (default), the method will automatically search the working directory for:

- `*.prmtop` - AMBER topology file
- `*.inpcrd`, `*.crd`, or `*.rst` - AMBER coordinate file
- `system.pdb` - System PDB file
- `bilayer*_lipid.pdb` - Bilayer PDB with CRYST1 record (highest priority)

**Manual Specification (Alternative):**

If auto-detection fails or you need specific files, provide a dictionary with these keys:

| Key | Type | Description |
|-----|------|-------------|
| `prmtop` | `str` | Path to AMBER topology file (system.prmtop) |
| `inpcrd` | `str` | Path to AMBER coordinate file (system.inpcrd) |
| `pdb` | `str` | Path to system PDB file (system.pdb) |
| `bilayer_pdb` | `str` | **REQUIRED** Path to bilayer PDB with CRYST1 record |

**Important: CRYST1 Box Dimensions**

The `bilayer_pdb` file **must contain** a CRYST1 record that defines the periodic box dimensions. This file is typically:

- Named `bilayer_*_lipid.pdb` (from packmol-memgen --parametrize)
- Generated during system preparation
- Contains the header: `CRYST1   70.335   70.833   85.067  90.00  90.00  90.00 P 1`

> **Warning:**
> Without proper CRYST1 information, NAMD will use incorrect box dimensions, leading to simulation failures.

### Stage Parameters Dictionary

Each stage in `stage_params_list` must be a dictionary with the following structure:

**Required Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `name` | `str` | Stage name (e.g., "Equilibration 1"). **Note:** Names are optional and used primarily for logging and user feedback. The actual configuration file names are generated sequentially (step1, step2, etc.) regardless of the name you provide |
| `time_ns` | `float` | Simulation time in nanoseconds |
| `steps` | `int` | Number of MD steps |
| `ensemble` | `str` | Ensemble type (NVT, NPT, NPAT, NPgT) |
| `temperature` | `float` | Temperature in Kelvin |
| `timestep` | `float` | Timestep in femtoseconds (1.0 or 2.0) |
| `constraints` | `Dict[str, float]` | Restraint forces (see below) |

**Optional Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `minimize_steps` | `int` | `0` | Minimization steps before MD (first stage only) |
| `pressure` | `float` | `1.0` | Pressure in bar (NPT, NPAT, NPgT only) |
| `surface_tension` | `float` | `0.0` | Surface tension in dyn/cm (NPAT only) |
| `dcd_freq` | `int` | `5000` | DCD output frequency (steps) |
| `use_gpu` | `bool` | `True` | Enable GPU acceleration |
| `cpu_cores` | `int` | `1` | Number of CPU cores |
| `gpu_id` | `int` | `0` | GPU device ID |
| `num_gpus` | `int` | `1` | Number of GPUs to use |
| `custom_template` | `str` | `None` | Explicitly specify CHARMM-GUI template file (e.g., `'step6.3_equilibration.inp'`) |

**Constraints Dictionary:**

Restraint forces in kcal/mol/Å². Applied to atomic groups:

| Key | Description | Typical Values |
|-----|-------------|----------------|
| `protein_backbone` | Protein backbone atoms (CA, C, N, O) | 10.0 → 5.0 → 1.0 → 0.0 |
| `protein_sidechain` | Protein sidechain heavy atoms | 5.0 → 2.5 → 0.5 → 0.0 |
| `lipid_head` | Lipid head group atoms (P, O) | 5.0 → 2.5 → 1.0 → 0.0 |
| `lipid_tail` | Lipid tail carbon atoms | 5.0 → 2.5 → 0.5 → 0.0 |
| `water` | Water molecules (usually unrestrained) | 0.0 |
| `ions` | Ion atoms (K+, Cl-, etc.) | 10.0 → 0.0 |
| `other` | Other molecules/ligands | 0.0 |

---

## Ensemble Types

| Ensemble | Full Name | Control | Use Case |
|----------|-----------|---------|----------|
| **NVT** | Canonical | Temperature | Initial heating, fixed box |
| **NPT** | Isothermal-isobaric | Temp + Pressure (isotropic) | General equilibration |
| **NPAT** | Constant surface tension | Temp + Pressure (anisotropic) + Surface tension | Membrane systems (recommended) |
| **NPgT** | Constant surface tension | Temp + Pressure + Surface tension | Alternative membrane ensemble |

---

## Advanced Features

### Per-Stage Ensemble Selection

Each stage can use a **different ensemble** by specifying the `ensemble` field. This allows flexible equilibration protocols:

```python
stages = [
    {'name': 'Equilibration 1', 'ensemble': 'NVT', ...},   # Heat with NVT
    {'name': 'Equilibration 2', 'ensemble': 'NPT', ...},   # Equilibrate with NPT
    {'name': 'Equilibration 3', 'ensemble': 'NPAT', ...},  # Relax membrane with NPAT
]
```

**How it works:**

- The `scheme_type` parameter is auto-detected from the first stage's `ensemble` field
- Each stage uses its own `ensemble` to select the appropriate CHARMM-GUI template
- When a stage's ensemble differs from the protocol default, a **warning is logged**
- Templates are automatically loaded from the correct ensemble folder (01_NVT, 02_NPT, 03_NPAT, 04_NPgT)

**Example: Mixed Ensemble Protocol**

```python
stages = [
    {
        'name': 'Initial Heating',
        'ensemble': 'NVT',      # Uses NVT template
        'time_ns': 0.25,
        'temperature': 310.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {'protein_backbone': 10.0, ...}
    },
    {
        'name': 'Pressure Equilibration',
        'ensemble': 'NPT',      # Uses NPT template (warning logged)
        'time_ns': 0.5,
        'temperature': 310.15,
        'pressure': 1.0,
        'timestep': 1.0,
        'constraints': {'protein_backbone': 5.0, ...}
    },
    {
        'name': 'Membrane Relaxation',
        'ensemble': 'NPAT',     # Uses NPAT template (warning logged)
        'time_ns': 1.0,
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {'protein_backbone': 1.0, ...}
    }
]

# Setup automatically detects scheme_type='NVT' from first stage
# Stages 2-3 use different ensembles with appropriate templates
result = manager.setup_namd_equilibration(stage_params_list=stages)
```

**Output:**
```
INFO - Auto-detected scheme_type from stages: NVT
INFO - Generated: step1_equilibration.conf
WARNING - Stage 2 (Pressure Equilibration) uses ensemble 'NPT' but protocol default is 'NVT'. Using 'NPT' template.
INFO - Generated: step2_equilibration.conf
WARNING - Stage 3 (Membrane Relaxation) uses ensemble 'NPAT' but protocol default is 'NVT'. Using 'NPAT' template.
INFO - Generated: step3_equilibration.conf
```

### Custom Template Selection

For advanced control, explicitly specify which CHARMM-GUI template to use with the `custom_template` parameter:

```python
stages = [
    {
        'name': 'Strong Restraints',
        'ensemble': 'NPT',
        'custom_template': 'step6.1_equilibration.inp',  # Use template from stage 1
        'time_ns': 0.25,
        'constraints': {'protein_backbone': 10.0, ...}
    },
    {
        'name': 'Medium Restraints',
        'ensemble': 'NPT',
        'custom_template': 'step6.3_equilibration.inp',  # Use template from stage 3
        'time_ns': 0.5,
        'constraints': {'protein_backbone': 5.0, ...}
    },
    {
        'name': 'Light Restraints',
        'ensemble': 'NPAT',
        'custom_template': 'step6.5_equilibration.inp',  # Use template from stage 5
        'time_ns': 1.0,
        'constraints': {'protein_backbone': 1.0, ...}
    }
]
```

**Use cases for custom templates:**

- Skip intermediate equilibration stages (use step6.5 directly)
- Repeat a specific stage with different parameters
- Mix templates from different equilibration phases
- Test different template configurations

**Available templates per ensemble:**

- `step6.1_equilibration.inp` through `step6.6_equilibration.inp` (6 equilibration stages)
- `step7_production.inp` (production stage)

### MDAnalysis Atom Counts and Selection Editing (GUI)

When using the Equilibration GUI frame, the constraint section for each stage now displays:

- **Atom count labels** next to each constraint entry, showing how many atoms match the selection
- **Gear button** (⚙) on each constraint row to edit the MDAnalysis selection string
- **Add Selection button** (+) at the bottom of each stage to add custom named selections
- **Auto-detect ligands** when an input folder is selected — non-standard residues are automatically added as `ligand_<RESNAME>` entries

**Gear Button (Selection Editor):**

Clicking the gear icon opens a modal dialog where you can:

- View the current MDAnalysis selection string
- Edit the selection expression
- Click **Test** to count matching atoms in the loaded PDB
- Click **Apply** to save the modified selection

**Add Selection Button:**

Clicking the **+** button opens a dialog to:

- Enter a custom name (e.g., `drug_A`)
- Provide an MDAnalysis selection string (e.g., `resname LIG and around 5 protein`)
- Set the default restraint force (kcal/mol/Å²)
- Test the selection before adding

New selections are added to **all stages** in the protocol.

---

### Custom Stage Names

Stage names can be anything - they don't need to follow the "Equilibration N" convention:

```python
stages = [
    {'name': 'Initial Heating', ...},           # Maps to step1
    {'name': 'Pressure Equilibration', ...},    # Maps to step2
    {'name': 'Membrane Relaxation', ...},       # Maps to step3
    {'name': 'Production Preparation', ...},    # Maps to step4
]
```

**How it works:**

- Custom names are automatically mapped to sequential step numbers (step1, step2, step3, ...)
- Config files use sequential naming: `step1_equilibration.conf`, `step2_equilibration.conf`, etc.
- Restart files are properly chained: `step2` loads from `step1`, `step3` loads from `step2`, etc.
- `inputname` variables are correctly set for all stages

---

## Working Examples

### Example 4: Three-Stage Protocol

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to system folder
work_dir = Path(__file__).parent / "popc_membrane"

stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 303.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 2',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 303.15,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 5.0,
            'protein_sidechain': 2.5,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 3',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 2.5,
            'protein_sidechain': 1.0,
            'lipid_head': 1.0,
            'lipid_tail': 1.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Auto-detect files and setup
# scheme_type auto-detected from stages (NVT from first stage)
# Stage 3 uses NPT ensemble - warning will be logged
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    stage_params_list=stages,
    output_name="equilibration_example_04",
    namd_executable="namd3"
)

print(f"\n✓ Setup complete!")
print(f"  Config files: {len(result['config_files'])}")
print(f"  Run script: {result['run_script'].name}")
```

### Example 5: Custom Four-Stage Protocol

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to system folder
work_dir = Path(__file__).parent / "popc_membrane"
system_files = {
    'prmtop': str(work_dir / 'system.prmtop'),
    'inpcrd': str(work_dir / 'system.inpcrd'),
    'pdb': str(work_dir / 'system.pdb'),
    'bilayer_pdb': str(work_dir / 'bilayer_protein_protonated_prepared_lipid.pdb')

custom_protocol = [
    {
        'name': 'Initial Equilibration',
        'time_ns': 0.25,
        'steps': 250000,
        'ensemble': 'NVT',
        'temperature': 310.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 5.0,
            'lipid_tail': 5.0,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    },
    {
        'name': 'Pressure Equilibration',
        'time_ns': 0.5,
        'steps': 500000,
        'ensemble': 'NPT',
        'temperature': 310.15,
        'pressure': 1.0,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 5.0,
            'protein_sidechain': 2.5,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Membrane Relaxation',
        'time_ns': 1.0,
        'steps': 500000,
        'ensemble': 'NPAT',
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 2.0,
            'protein_sidechain': 1.0,
            'lipid_head': 1.0,
            'lipid_tail': 0.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Production Preparation',
        'time_ns': 2.0,
        'steps': 1000000,
        'ensemble': 'NPAT',
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.5,
            'protein_sidechain': 0.0,
            'lipid_head': 0.0,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Auto-detect and setup
# scheme_type auto-detected from stages (NVT from first stage)
# Stages 2-4 use different ensembles - warnings will be logged
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    stage_params_list=custom_protocol,
    output_name="equilibration_example_05",
    namd_executable="namd3"
)

print(f"\n✓ Setup complete!")
print(f"  Total stages: {len(custom_protocol)}")
print(f"  Total time: {sum(s['time_ns'] for s in custom_protocol):.1f} ns")
print(f"\nTo run:")
print(f"  cd {result['namd_dir']}")
print(f"  ./run_equilibration.sh")
```

### Example 6: Complete CHARMM-GUI 7-Stage Protocol in NPT ensemble

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# System folder
work_dir = Path(__file__).parent / "popc_membrane"

stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'minimize_steps': 10000,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 2',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 5.0,
            'protein_sidechain': 2.5,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 3',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 2.5,
            'protein_sidechain': 1.0,
            'lipid_head': 1.0,
            'lipid_tail': 1.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 4',
        'time_ns': 0.5,
        'steps': 250000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 1.0,
            'protein_sidechain': 0.5,
            'lipid_head': 0.5,
            'lipid_tail': 0.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 5',
        'time_ns': 0.5,
        'steps': 250000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.5,
            'protein_sidechain': 0.1,
            'lipid_head': 0.1,
            'lipid_tail': 0.1,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 6',
        'time_ns': 0.5,
        'steps': 250000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.1,
            'protein_sidechain': 0.0,
            'lipid_head': 0.0,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Production',
        'time_ns': 10.0,
        'steps': 5000000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.0,
            'protein_sidechain': 0.0,
            'lipid_head': 0.0,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Auto-detect and setup
# scheme_type auto-detected from first stage's ensemble
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    stage_params_list=stages,
    output_name="equilibration_example_06",
    namd_executable="namd3"
)

print(f"\n✓ Complete! Generated {len(result['config_files'])} configuration files")
print(f"  Total equilibration: {sum(s['time_ns'] for s in stages[:-1]):.3f} ns")
print(f"  Production: {stages[-1]['time_ns']:.1f} ns")
```

### Example 7: Custom Template Selection

This example demonstrates explicit template control for advanced workflows:

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# System folder
work_dir = Path(__file__).parent / "popc_membrane"

system_files = {
    'prmtop': str(work_dir / 'system.prmtop'),
    'inpcrd': str(work_dir / 'system.inpcrd'),
    'pdb': str(work_dir / 'system.pdb'),
    'bilayer_pdb': str(work_dir / 'bilayer_protein_protonated_prepared_lipid.pdb')
}

# Use explicit templates to skip intermediate stages
stages = [
    {
        'name': 'Strong Restraints Phase',
        'ensemble': 'NPT',
        'custom_template': 'step6.1_equilibration.inp',  # Use stage 1 template
        'time_ns': 0.25,
        'timestep': 1.0,
        'temperature': 310.15,
        'pressure': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 5.0,
            'lipid_tail': 5.0,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    },
    {
        'name': 'Medium Restraints Phase',
        'ensemble': 'NPT',
        'custom_template': 'step6.3_equilibration.inp',  # Skip to stage 3 template
        'time_ns': 0.5,
        'timestep': 1.0,
        'temperature': 310.15,
        'pressure': 1.0,
        'constraints': {
            'protein_backbone': 5.0,
            'protein_sidechain': 2.5,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Light Restraints Phase',
        'ensemble': 'NPAT',
        'custom_template': 'step6.5_equilibration.inp',  # Use stage 5 template
        'time_ns': 1.0,
        'timestep': 2.0,
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'constraints': {
            'protein_backbone': 1.0,
            'protein_sidechain': 0.5,
            'lipid_head': 0.5,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Setup with custom template selection
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    output_name="equilibration_example_07",
    namd_executable="namd3"
)

print(f"\n✓ Setup complete with custom templates!")
print(f"  Stage 1: Using template step6.1")
print(f"  Stage 2: Using template step6.3 (skipped step6.2)")
print(f"  Stage 3: Using template step6.5 (skipped step6.4)")
```

**When to use custom templates:**

- **Skip intermediate stages**: Jump from step6.1 → step6.3 → step6.5
- **Repeat a stage**: Use step6.2 multiple times with different parameters
- **Mix templates**: Combine templates from different equilibration phases
- **Test protocols**: Experiment with different template combinations

---

### Example 8: MDAnalysis Selections for Restraints

This example demonstrates the MDAnalysis-based selection system for precise atom counting and restraint generation, including auto-detection of non-standard residues (ligands, ions):

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to the system folder
work_dir = Path("popc_membrane")
system_pdb = work_dir / "bilayer_protein_protonated_prepared_lipid.pdb"

manager = NAMDEquilibrationManager(work_dir)

# 1. Inspect default selections and atom counts
for name, sel in NAMDEquilibrationManager.DEFAULT_SELECTIONS.items():
    count = NAMDEquilibrationManager.count_selection_atoms(str(system_pdb), sel)
    print(f"  {name:25s}  →  {count:>7d} atoms")

# 2. Auto-detect ligands / non-standard residues
all_sels = NAMDEquilibrationManager.get_default_selections(str(system_pdb))
for name, sel in all_sels.items():
    if name.startswith("ligand_"):
        count = NAMDEquilibrationManager.count_selection_atoms(str(system_pdb), sel)
        print(f"  {name:25s}  →  {count:>7d} atoms  |  {sel}")

# 3. Count all selections at once
counts = NAMDEquilibrationManager.count_all_selections(str(system_pdb))

# 4. Generate restraints PDB via MDAnalysis selections
output_file = work_dir / "namd" / "restraints" / "step1_restraints.pdb"
selections_with_forces = {
    "protein_backbone":  ("protein and backbone", 10.0),
    "protein_sidechain": ("protein and not backbone", 5.0),
    "lipid_head":        (NAMDEquilibrationManager.DEFAULT_SELECTIONS["lipid_head"], 2.5),
    "lipid_tail":        (NAMDEquilibrationManager.DEFAULT_SELECTIONS["lipid_tail"], 2.5),
    "water":             (NAMDEquilibrationManager.DEFAULT_SELECTIONS["water"], 0.0),
    "ions":              (NAMDEquilibrationManager.DEFAULT_SELECTIONS["ions"], 10.0),
}
# Add any auto-detected ligand with force 1.0
for name, sel in all_sels.items():
    if name.startswith("ligand_"):
        selections_with_forces[name] = (sel, 1.0)

manager.generate_restraints_file_mda(
    system_pdb, selections_with_forces, output_file,
    stage_name="Equilibration 1",
)

# 5. Or use the high-level API with selections parameter
constraints = {"protein_backbone": 10.0, "protein_sidechain": 5.0, "lipid_head": 2.5,
               "lipid_tail": 2.5, "water": 0.0, "ions": 10.0}
selections = {name: sel for name, (sel, _) in selections_with_forces.items()}

manager.generate_restraints_file(
    system_pdb, constraints, output_file,
    stage_name="Eq1", selections=selections,
)
```

**Output:**
```
  protein_backbone            →      204 atoms
  protein_sidechain           →      488 atoms
  lipid_head                  →     2904 atoms
  lipid_tail                  →    13310 atoms
  water                       →    14685 atoms
  ions                        →        0 atoms
  other                       →       21 atoms
  ligand_Cl-                  →       10 atoms  |  resname Cl-
  ligand_K+                   →       11 atoms  |  resname K+
```

---

### Example 9: Custom Restraints — Three Levels of Customisation

This example demonstrates three progressive levels of customisation for NAMD restraints using `EquilibrationStage.replace()`:

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager, EquilibrationStage

work_dir = Path(__file__).parent / "popc_membrane"
system_files = {
    "prmtop": str(work_dir / "system.prmtop"),
    "inpcrd": str(work_dir / "system.inpcrd"),
    "pdb": str(work_dir / "system.pdb"),
    "bilayer_pdb": str(work_dir / "bilayer_protein_protonated_prepared_lipid.pdb"),
}

WORK_DIR = work_dir
manager = NAMDEquilibrationManager(working_dir=WORK_DIR)

# ---------------------------------------------------------------------------
# Level 1 — Override a single force constant key (no MDAnalysis needed)
# ---------------------------------------------------------------------------
# Turn off sidechain restraints entirely; keep default backbone + lipid forces.

print("=== Level 1: Override protein_sidechain to 0 ===")
stages_l1 = [
    s.replace(constraints={**s.constraints, "protein_sidechain": 0.0})
    for s in NAMDEquilibrationManager.get_default_stage_params()
]

result_l1 = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages_l1,
    output_name="level1_no_sc",
)
print(f"Output: {result_l1['namd_dir']}")

# ---------------------------------------------------------------------------
# Level 2 — Override selections for standard categories (MDAnalysis)
# ---------------------------------------------------------------------------
# Useful when your PSF uses non-standard segment names or residue types.

print("\n=== Level 2: Custom selections for standard categories ===")
stages_l2 = NAMDEquilibrationManager.get_default_stage_params()

result_l2 = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages_l2,
    output_name="level2_custom_sel",
    selections={
        "protein_backbone": "backbone",
        "protein_sidechain": "protein and not backbone",
        "lipid_head": "resname POPC and name P O11 O12 O13 O14",
        "lipid_tail": "resname POPC and not (name P O11 O12 O13 O14 N)",
    },
)
print(f"Output: {result_l2['namd_dir']}")

# ---------------------------------------------------------------------------
# Level 3 — Full MDAnalysis control with a custom atom category
# ---------------------------------------------------------------------------
# Restrain ions in the first 3 stages at 10 kcal/mol/Å², then release.
# Replace "ions" with "ligand_ABC" and "resname ABC" for a real ligand system.

print("\n=== Level 3: Custom ion restraints (demonstrates ligand-style) ===")
stages_l3 = NAMDEquilibrationManager.get_default_stage_params()

# Apply 10 kcal/mol/Å² to ions in stages 1-3; zero thereafter
stages_l3_dicts = []
for i, s in enumerate(stages_l3):
    ion_force = 10.0 if i < 3 else 0.0
    new_constraints = {**s.constraints, "custom_ions": ion_force}
    stages_l3_dicts.append(s.replace(constraints=new_constraints).to_dict())

result_l3 = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages_l3_dicts,
    output_name="level3_custom_ions",
    selections={
        "custom_ions": "resname SOD CLA POT",  # MDAnalysis selection
    },
)
print(f"Output:           {result_l3['namd_dir']}")
```

**When to use each level:**

- **Level 1** (`replace(constraints=...)`): Change only force constants. No MDAnalysis needed — fastest approach.
- **Level 2** (`selections=` parameter): Override which atoms are selected for standard categories (backbone, sidechain, lipid head/tail) using your own MDAnalysis selection strings.
- **Level 3** (custom keys): Add entirely new atom categories (ligands, ions, cofactors) not covered by the default keys. Set the selection string and per-stage force schedule.

---

## Best Practices

### Restraint Progression

**Recommended restraint schedule** (kcal/mol/Å²):

| Component | Stage 1 | Stage 2 | Stage 3 | Stage 4 | Production |
|-----------|---------|---------|---------|---------|------------|
| Protein backbone | 10.0 | 5.0 | 2.5 | 1.0 | 0.0 |
| Protein sidechain | 5.0 | 2.5 | 1.0 | 0.5 | 0.0 |
| Lipid heads | 5.0 | 2.5 | 1.0 | 0.5 | 0.0 |
| Lipid tails | 5.0 | 2.5 | 1.0 | 0.0 | 0.0 |
| Ions | 10.0 | 0.0 | 0.0 | 0.0 | 0.0 |
| Water | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

**Guidelines:**

- Start with **strong restraints** (10.0) on protein backbone and ions
- **Gradually reduce** restraints over 3-6 stages
- Always keep **water unrestrained** (0.0)
- Release **ions early** (after stage 1)
- Release **lipid tails** before heads
- Final stage should have **minimal or zero** restraints

### Timestep Progression

| Stage | Timestep | Restraints | Notes |
|-------|----------|------------|-------|
| 1-2 | 1.0 fs | Strong (10.0-5.0) | Initial equilibration, constrained |
| 3-4 | 1.0-2.0 fs | Medium (2.5-1.0) | Transition to larger timestep |
| 5+ | 2.0 fs | Light (< 1.0) | Standard production timestep |

**Rules:**

- Use **1.0 fs** with strong restraints (> 5.0 kcal/mol/Å²)
- Switch to **2.0 fs** when restraints < 2.5 kcal/mol/Å²
- NAMD supports **2.0 fs** with SHAKE/SETTLE for bonds to hydrogen
- Never use > 2.0 fs for equilibration

### Stage Lengths

**Minimum recommended times:**

| Ensemble | Stage Purpose | Min Time | Typical Time |
|----------|---------------|----------|--------------|
| NVT | Initial heating | 0.1 ns | 0.125-0.25 ns |
| NPT | Pressure equilibration | 0.1 ns | 0.125-0.5 ns |
| NPAT | Membrane equilibration | 0.5 ns | 0.5-2.0 ns |
| Production | Data collection | 10 ns | 50-500 ns |

**Total equilibration time:**

- **Minimum:** 1.0 ns (3-4 stages)
- **Standard:** 2.0 ns (6 stages, CHARMM-GUI protocol)
- **Thorough:** 3-5 ns (custom research protocols)

### Ensemble Selection

**For membrane protein systems:**

**Start with NVT** (0.1-0.25 ns)

   - Heat system from 0 K to target temperature
   - Strong restraints on protein and ions
   - Fixed box dimensions

**Switch to NPT** (0.1-0.5 ns)

   - Equilibrate pressure (isotropic)
   - Medium restraints
   - Box can change size uniformly

**Finish with NPAT** (0.5-2.0 ns)

   - Equilibrate membrane (anisotropic pressure)
   - Light restraints
   - Lateral box dimensions independent from Z-axis
   - **Recommended for membrane systems**

### Minimization

- Use **minimize_steps** only in **first stage**
- Typical value: **5000-10000 steps**
- Removes bad contacts/clashes
- Always minimize before MD

### GPU Acceleration

```python
# Enable GPU acceleration (default)
stage = {
    'use_gpu': True,
    'gpu_id': 0,      # First GPU
    'num_gpus': 1,    # Single GPU
    'cpu_cores': 1    # Minimal CPU usage with GPU. Recommended > 4 depending on the system´s size.
}
```

---

## Output Structure

```
{output_name}/
└── namd/
    ├── system.prmtop              # Copied AMBER topology
    ├── system.inpcrd              # Copied AMBER coordinates  
    ├── system.pdb                 # Copied system PDB (for structure/restraints)
    ├── bilayer_*_lipid.pdb        # Copied bilayer PDB (for CRYST1 box info only)
    ├── step1.conf                 # Stage 1 config
    ├── step2.conf                 # Stage 2 config
    ├── ...
    ├── run_equilibration.sh       # Executable run script
    ├── protocol_summary.json      # Protocol metadata
    ├── restraints/
        ├── step1_equilibration_restraints.pdb   # Stage 1 restraints
        ├── step2_equilibration_restraints.pdb   # Stage 2 restraints
    └── logs/
        ├── step1.log              # Stage 1 output
        ├── step2.log              # Stage 2 output
        └── ...
```

**Generated Files:**

| File | Purpose |
|------|---------|
| `system.pdb` | System structure for NAMD simulations and restraints |
| `bilayer_*_lipid.pdb` | Original bilayer file (kept for reference, CRYST1 info read from this) |
| `stepN.conf` | NAMD configuration for stage N |
| `stepN_restraints.txt` | Harmonic restraints for stage N |
| `run_equilibration.sh` | Bash script to run all stages sequentially |
| `protocol_summary.json` | Protocol metadata (stages, parameters, files) |
| `logs/stepN.log` | NAMD output log for stage N |

**Important Notes:**

- `system.pdb` is the actual system file used for simulations and restraint generation
- `bilayer_*_lipid.pdb` is kept with its original name for reference and CRYST1 box dimensions
- Both files are copied to the output directory but serve different purposes

**Run Script Features:**

- Sequential execution of all stages
- Automatic dependency management (restart files)
- Progress tracking
- Error handling
- Timing information

---

## Internal Methods (Advanced)

### Class Attribute: DEFAULT_SELECTIONS

Default MDAnalysis selection strings for the seven standard restraint categories.
These selections are used when MDAnalysis-based restraint generation is enabled.

```python
NAMDEquilibrationManager.DEFAULT_SELECTIONS = {
    'protein_backbone':  'protein and backbone',
    'protein_sidechain': 'protein and not backbone',
    'lipid_head':        '(resname POPC POPE POPS DPPC ...) and (name P O11 O12 ...)',
    'lipid_tail':        '(resname POPC POPE POPS DPPC ...) and not (name P O11 ...)',
    'water':             'resname TIP3 HOH WAT SOL TIP4 SPC T3P T4P',
    'ions':              'resname NA CL K CA MG ZN FE CU SOD CLA POT CAL MAG ZIN IRN COP',
    'other':             'not (protein or lipids or water or ions)',
}
```

Each value is a valid [MDAnalysis selection string](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).
Users can override any selection through the GUI gear button or the API `selections` parameter.

---

### Static Method: count_selection_atoms()

Count atoms matching an MDAnalysis selection expression.

```python
NAMDEquilibrationManager.count_selection_atoms(pdb_path: str, selection: str) -> int
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_path` | `str` | Path to a PDB file |
| `selection` | `str` | MDAnalysis selection string |

**Returns:** `int` — number of matching atoms (0 if selection is invalid or MDAnalysis unavailable)

```python
count = NAMDEquilibrationManager.count_selection_atoms(
    "system.pdb", "protein and backbone"
)
print(f"Backbone atoms: {count}")
```

---

### Static Method: get_default_selections()

Build the default selection dict, auto-detecting extra ligands / non-standard residues.

```python
NAMDEquilibrationManager.get_default_selections(pdb_path: str) -> Dict[str, str]
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pdb_path` | `str` | Path to a PDB file |

**Returns:** `Dict[str, str]` — `{category_name: mda_selection_string, ...}` including any
auto-detected `ligand_<RESNAME>` entries.

The seven standard categories are always present. Any residue that falls into the *other*
category is additionally split into individual `ligand_<RESNAME>` entries so users can
assign per-ligand restraint forces.

```python
sels = NAMDEquilibrationManager.get_default_selections("system.pdb")
for name, sel in sels.items():
    if name.startswith("ligand_"):
        print(f"  Detected: {name} → {sel}")
# Output:
#   Detected: ligand_Cl- → resname Cl-
#   Detected: ligand_K+  → resname K+
```

---

### Static Method: count_all_selections()

Count atoms for every selection in a dictionary in one call.

```python
NAMDEquilibrationManager.count_all_selections(
    pdb_path: str,
    selections: Optional[Dict[str, str]] = None
) -> Dict[str, int]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_path` | `str` | Required | Path to a PDB file |
| `selections` | `Optional[Dict[str, str]]` | `None` | `{name: mda_selection_string}`. If `None`, uses `get_default_selections()` (with auto-detected ligands) |

**Returns:** `Dict[str, int]` — `{name: atom_count, ...}`

```python
counts = NAMDEquilibrationManager.count_all_selections("system.pdb")
for name, n in counts.items():
    print(f"  {name:25s}  {n:>7d} atoms")
```

---

### Method: generate_restraints_file_mda()

Generate a restraints PDB using MDAnalysis selections instead of the built-in heuristic.
Each entry maps a category name to a `(mda_selection_string, force)` tuple.
For every ATOM/HETATM line the **first matching** selection determines the B-factor.
Atoms matching no selection get B-factor 0.0.

```python
generate_restraints_file_mda(
    system_pdb: Path,
    selections_with_forces: Dict[str, Tuple[str, float]],
    output_file: Path,
    stage_name: str = ""
) -> None
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `system_pdb` | `Path` | Path to the system PDB file |
| `selections_with_forces` | `Dict[str, Tuple[str, float]]` | `{name: (selection_string, force), ...}` |
| `output_file` | `Path` | Destination path for the restraints PDB |
| `stage_name` | `str` | Label for log messages (optional) |

```python
selections_with_forces = {
    "protein_backbone":  ("protein and backbone", 10.0),
    "protein_sidechain": ("protein and not backbone", 5.0),
    "lipid_head":        (DEFAULT_SELS["lipid_head"], 2.5),
    "lipid_tail":        (DEFAULT_SELS["lipid_tail"], 2.5),
    "water":             (DEFAULT_SELS["water"], 0.0),
    "ions":              (DEFAULT_SELS["ions"], 10.0),
    "ligand_Cl-":        ("resname Cl-", 1.0),
}

manager.generate_restraints_file_mda(
    system_pdb, selections_with_forces, output_file,
    stage_name="Equilibration 1",
)
```

---

### Method: generate_restraints_file()

Generates restraint PDB files with B-factors encoding restraint forces for each atom type.
Now supports optional MDAnalysis selections for precise atom classification.

```python
generate_restraints_file(
    system_pdb: Path,
    constraints: Dict[str, float],
    output_file: Path,
    stage_name: str = "",
    selections: Optional[Dict[str, str]] = None
) -> None
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `system_pdb` | `Path` | Required | Path to the **system.pdb** file (full system: protein + lipids + water + ions) |
| `constraints` | `Dict[str, float]` | Required | Dictionary of restraint forces (kcal/mol/Å²) for each atom type |
| `output_file` | `Path` | Required | Path for output restraints PDB file |
| `stage_name` | `str` | `""` | Stage name for logging (optional) |
| `selections` | `Optional[Dict[str, str]]` | `None` | `{constraint_name: mda_selection_string}`. When provided, MDAnalysis is used instead of the built-in heuristic |

**Behavior:**

When `selections` is `None` (legacy mode):

A. Classifies each atom by type using built-in residue/atom-name heuristics:

   - `protein_backbone`: CA, C, N, O atoms in protein residues
   - `protein_sidechain`: Heavy sidechain atoms in protein residues
   - `lipid_head`: P, O atoms in lipid head groups
   - `lipid_tail`: C atoms in lipid tails
   - `water`: H2O molecules (TIP3, HOH, WAT, SOL)
   - `ions`: Na+, Cl-, K+, Ca2+, Mg2+, etc.
   - `other`: Ligands and other molecules

B. Assigns B-factor values based on constraints dictionary

C. Writes restraints PDB file with modified B-factors

When `selections` is provided (MDAnalysis mode):

A. Pairs each constraint name with its MDAnalysis selection string

B. Delegates to `generate_restraints_file_mda()` for precise MDAnalysis-based classification

C. Useful for custom selections, auto-detected ligands, or non-standard residues

```python
# Legacy mode (heuristic-based)
manager.generate_restraints_file(system_pdb, constraints, output_file)

# MDAnalysis mode (selection-based)
selections = {
    "protein_backbone": "protein and backbone",
    "protein_sidechain": "protein and not backbone",
    "ligand_Cl-": "resname Cl-",
}
manager.generate_restraints_file(
    system_pdb, constraints, output_file,
    selections=selections,
)
```

**Example Output Log:**
```
INFO - Generated restraints file: step1_equilibration_restraints.pdb
INFO - Source PDB: system.pdb
INFO - Stage: Equilibration 1
INFO - Total atoms processed: 45678
INFO -   protein_backbone: 1234 atoms, force = 10.0 kcal/mol/Å²
INFO -   protein_sidechain: 2345 atoms, force = 5.0 kcal/mol/Å²
INFO -   lipid_head: 512 atoms, force = 5.0 kcal/mol/Å²
INFO -   lipid_tail: 3456 atoms, force = 5.0 kcal/mol/Å²
INFO -   ions: 89 atoms, force = 10.0 kcal/mol/Å²
```

---

### Method: generate_charmm_gui_config_file()

Generates NAMD configuration files using CHARMM-GUI templates with GateWizard customizations.

```python
generate_charmm_gui_config_file(
    stage_name: str,
    stage_params: Dict[str, Any],
    stage_index: int,
    system_files: Dict[str, str],
    scheme_type: str,
    previous_stage_name: Optional[str] = None,
    all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None,
    force_scheme_type: bool = False
) -> str
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `stage_name` | `str` | Name of the equilibration stage |
| `stage_params` | `Dict[str, Any]` | Stage parameters dictionary |
| `stage_index` | `int` | Stage index (0-based) |
| `system_files` | `Dict[str, str]` | System file paths |
| `scheme_type` | `str` | CHARMM-GUI scheme type (NVT, NPT, NPAT, NPgT) |
| `previous_stage_name` | `Optional[str]` | Previous stage name for restart files |
| `all_stage_settings` | `Optional[Dict]` | All stages settings for context |
| `force_scheme_type` | `bool` | **[NEW]** If True, always use `scheme_type` for all stages (GUI mode) |

**Returns:** NAMD configuration file content as string

**Behavior:**

**API Mode (force_scheme_type=False, default):**

- Each stage can specify its own `ensemble` field
- Template folder is selected based on the stage's ensemble value
- Allows mixing NVT/NPT/NPAT templates across stages
- Warning logged when stage ensemble differs from protocol scheme_type

**GUI Mode (force_scheme_type=True):**

- ALL stages use the global `scheme_type` for template selection
- Individual stage `ensemble` values are ignored for template selection
- Ensures CHARMM-GUI protocol consistency when selected from GUI dropdown
- Example: If "NPT" scheme selected, all stages use NPT templates (02_NPT folder)

**Template Selection Logic:**

1. If `custom_template` specified → use that template file
2. If `force_scheme_type=True` → use `scheme_type` for folder selection
3. If `force_scheme_type=False` → use stage's `ensemble` for folder selection
4. Templates loaded from: `{scheme_folder}/{template_file}` under `equilibration/namd/`

**Scheme to Folder Mapping:**

- `NVT` → `01_NVT/`
- `NPT` → `02_NPT/`
- `NPAT` → `03_NPAT/`
- `NPgT` → `04_NPgT/`

**Example (API Mode):**
```python
# Stage 1 uses NVT template, Stage 2 uses NPT template
config1 = manager.generate_charmm_gui_config_file(
    stage_name="Heating",
    stage_params={'ensemble': 'NVT', ...},
    stage_index=0,
    system_files=files,
    scheme_type='NVT',
    force_scheme_type=False  # API mode - use stage's ensemble
)

config2 = manager.generate_charmm_gui_config_file(
    stage_name="Pressure Eq",
    stage_params={'ensemble': 'NPT', ...},  # Different ensemble
    stage_index=1,
    system_files=files,
    scheme_type='NVT',  # Protocol default is NVT
    force_scheme_type=False  # Stage uses NPT template (warning logged)
)
```

**Example (GUI Mode):**
```python
# User selected "NPT" scheme in GUI dropdown
# ALL stages use NPT templates regardless of individual ensemble values
config1 = manager.generate_charmm_gui_config_file(
    stage_name="Stage 1",
    stage_params={'ensemble': 'NVT', ...},  # Ignored for template selection
    stage_index=0,
    system_files=files,
    scheme_type='NPT',  # From GUI dropdown
    force_scheme_type=True  # GUI mode - enforce scheme_type
)
# Uses template from: 02_NPT/step6.1_equilibration.inp

config2 = manager.generate_charmm_gui_config_file(
    stage_name="Stage 2",
    stage_params={'ensemble': 'NVT', ...},  # Ignored for template selection
    stage_index=1,
    system_files=files,
    scheme_type='NPT',  # From GUI dropdown
    force_scheme_type=True  # GUI mode - enforce scheme_type
)
# Uses template from: 02_NPT/step6.2_equilibration.inp
```

**Rationale for force_scheme_type:**

The `force_scheme_type` parameter was added to support two different use cases:

1. **GUI Users:** Select a CHARMM-GUI scheme (NPT, NPAT, etc.) and expect all stages to use templates from that scheme folder for consistency with CHARMM-GUI protocols.

2. **API Users:** Have flexibility to mix different ensemble templates across stages for custom equilibration protocols (e.g., NVT → NPT → NPAT progression).

---

## Common Issues

### Missing CRYST1 Record

**Error:** Box dimensions not found in PDB file

**Solution:** Ensure `bilayer_pdb` contains CRYST1 record:
```bash
head bilayer_protein_protonated_prepared_lipid.pdb
# Should show: CRYST1   70.335   70.833   85.067  90.00  90.00  90.00 P 1
```

**Correct file:**

- From packmol-memgen `--parametrize` step
- Named `bilayer_*_lipid.pdb`
- Contains proper PDB header with CRYST1
- Not an intermediate Packmol file

### Restraint File Errors

**Error:** Cannot generate restraint file OR restraints only applied to protein

**Solution:** Ensure **system.pdb** is used (not protein.pdb):

- The restraints generation requires the **full system.pdb** file
- This file must contain: protein + lipids + water + ions
- Do NOT use protein.pdb (protein-only) for restraint generation
- Check that system.pdb exists in the working directory or input folder
- Verify atom names match AMBER topology
- Verify segnames/chains are present (PROA, MEMB, SOLV)

**Common mistake:** Using protein.pdb instead of system.pdb results in:

- No restraints applied to lipids (lipid_head, lipid_tail)
- No restraints applied to ions
- Membrane collapse during equilibration
- System instability

**Fix:** The GUI now automatically uses system.pdb from the output directory where files are copied. If you see restraints being generated with protein.pdb, check:

1. Is system.pdb present in the input folder?
2. Is system.pdb copied to the output directory?
3. Check the log for "Using PDB file for restraints generation: ..."

### Timestep Too Large

**Error:** NAMD crashes with "Atoms moving too fast"

**Solution:** Reduce timestep or increase restraints:

- Use 1.0 fs with restraints > 5.0 kcal/mol/Å²
- Switch to 2.0 fs only when restraints < 2.5 kcal/mol/Å²
- Add minimization step if starting from poor geometry

### Ensemble Switching

**Error:** System unstable when switching NPT → NPAT

**Solution:** Ensure pressure is equilibrated first:

- Run NPT for at least 0.25 ns before NPAT
- Check box dimensions are stable in NPT
- Use medium restraints during transition
- Don't switch ensemble and reduce restraints simultaneously

---

---

## Class: OpenMMEquilibrationManager

Manager for OpenMM equilibration simulations using CHARMM-GUI templates and the AMBER force field.
Uses the same user-facing API as `NAMDEquilibrationManager`.

### Constructor

```python
OpenMMEquilibrationManager(working_dir: Path)
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `working_dir` | `Path` | Required | Working directory containing system files |

---

### Example 1: Single NVT Stage

```python
from pathlib import Path
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

# Point to folder with system files
work_dir = Path(__file__).parent / "popc_membrane"

# Define a single NVT equilibration stage
stages = [
    {
        "name": "Equilibration 1",
        "time_ns": 0.125,
        "ensemble": "NVT",
        "temperature": 310.15,
        "timestep": 1.0,
        "minimize_steps": 5000,
        "constraints": {
            "protein_backbone": 10.0,
            "protein_sidechain": 5.0,
            "lipid_head": 2.5,
            "lipid_tail": 0.0,
        },
    }
]

# Setup with automatic file detection (no system_files needed!)
# scheme_type is auto-detected from the 'ensemble' field in stages
manager = OpenMMEquilibrationManager(work_dir)
result = manager.setup_openmm_equilibration(
    stage_params_list=stages,
    output_name="openmm_example_01",
)

print(f"Setup complete: {result['openmm_dir']}")
# Run with: cd {result['openmm_dir']} && bash run_equilibration.sh
```

---

### Example 2: Full CHARMM-GUI Protocol (7 stages)

```python
from pathlib import Path
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

# Point to folder with system files
work_dir = Path(__file__).parent / "popc_membrane"

# Full 6-stage NPT membrane equilibration protocol (CHARMM-GUI style)
# Gradual relaxation of restraints following the standard protocol
stages = [
    {
        "name": "Equilibration 1 - NVT with strong restraints",
        "time_ns": 0.125,
        "ensemble": "NVT",
        "temperature": 303.15,
        "timestep": 1.0,
        "minimize_steps": 5000,
        "constraints": {
            "protein_backbone": 10.0,
            "protein_sidechain": 5.0,
            "lipid_head": 2.5,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 2 - NVT relaxing restraints",
        "time_ns": 0.125,
        "ensemble": "NVT",
        "temperature": 303.15,
        "timestep": 1.0,
        "constraints": {
            "protein_backbone": 5.0,
            "protein_sidechain": 2.5,
            "lipid_head": 1.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 3 - NPT with pressure coupling",
        "time_ns": 0.125,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 1.0,
        "constraints": {
            "protein_backbone": 2.5,
            "protein_sidechain": 1.0,
            "lipid_head": 0.5,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 4 - NPT further relaxing",
        "time_ns": 0.25,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 1.0,
            "protein_sidechain": 0.5,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 5 - NPT backbone only",
        "time_ns": 0.25,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 0.5,
            "protein_sidechain": 0.0,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 6 - NPT light backbone restraints",
        "time_ns": 0.5,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 0.1,
            "protein_sidechain": 0.0,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Production - NPT unrestrained",
        "time_ns": 50.0,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 0.0,
            "protein_sidechain": 0.0,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
]

# Setup with automatic file detection
# scheme_type is auto-detected from 'ensemble' field of first stage (NVT -> 01_NVT)
# Note: Mixed ensembles (NVT stages 1-2, NPT stages 3-7) are handled automatically.
#       The scheme_type controls which pressure coupling templates are used for
#       stages 3+ — pass scheme_type="NPT" explicitly if needed.
manager = OpenMMEquilibrationManager(work_dir)
result = manager.setup_openmm_equilibration(
    stage_params_list=stages,
    output_name="openmm_example_02",
    scheme_type="NPT",
)

print(f"Setup complete: {result['openmm_dir']}")
print(f"Config files: {len(result['config_files'])}")
print(f"Run script: {result['run_script'].name}")
# Run with: cd {result['openmm_dir']} && bash run_equilibration.sh
```

---

### Example 3: Custom Ligand Restraints

```python
from pathlib import Path
from gatewizard.tools.equilibration import (
    OpenMMEquilibrationManager,
    EquilibrationStage,
)

WORK_DIR = Path("openmm_ligand_restraints")
WORK_DIR.mkdir(exist_ok=True)

manager = OpenMMEquilibrationManager(working_dir=WORK_DIR)

system_files = {
    "prmtop": "system.prmtop",
    "inpcrd": "system.inpcrd",
    "pdb": "system.pdb",
}

# --- Standard protein + lipid restraints (auto-detected) ---
print("=== Example 1: Standard protein/lipid restraints ===")
stages = OpenMMEquilibrationManager.get_default_stage_params()

result = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    output_name="standard_restraints",
)
print(f"OpenMM dir:      {result['openmm_dir']}")
print(f"Restraint files: {result['restraint_files']}")
# → restraint_files["prot_pos"]  = Path(".../restraints/prot_pos.txt")
# → restraint_files["lipid_pos"] = Path(".../restraints/lipid_pos.txt")  (if lipid forces > 0)
# → restraint_files["custom_pos"] = None

# --- Add ligand ABC restraints in stages 1-3 ---
print("\n=== Example 2: Ligand ABC restraints in stages 1-3 ===")
raw_stages = OpenMMEquilibrationManager.get_default_stage_params()
stage_objs = [EquilibrationStage(**s) for s in raw_stages]

# Apply 5 kcal/mol/Å² to ligand ABC in the first 3 stages; zero thereafter
stage_dicts = []
for i, s in enumerate(stage_objs):
    ligand_force = 5.0 if i < 3 else 0.0
    new_constraints = {**s.constraints, "ligand_ABC": ligand_force}
    stage_dicts.append(s.replace(constraints=new_constraints).to_dict())

result2 = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=stage_dicts,
    output_name="ligand_ABC_restraints",
    selections={
        "ligand_ABC": "resname ABC",  # MDAnalysis selection string
    },
)
print(f"OpenMM dir:      {result2['openmm_dir']}")
print(f"Restraint files: {result2['restraint_files']}")
# → restraint_files["custom_pos"] = Path(".../restraints/custom_pos.txt")
#   custom_pos.txt force = 5.0 kcal/mol/Å² × 418.4 = 2092.0 kJ/mol/nm²

# --- Custom backbone taper + ligand restraints ---
print("\n=== Example 3: Custom backbone taper + ligand ABC ===")
raw_stages = OpenMMEquilibrationManager.get_default_stage_params()

# Apply a linear backbone taper and add the ligand
bb_schedule = [10.0, 5.0, 2.5, 1.0, 0.5, 0.0]
sc_schedule = [5.0, 2.5, 1.0, 0.5, 0.0, 0.0]
lig_schedule = [5.0, 5.0, 5.0, 0.0, 0.0, 0.0]

stage_dicts3 = []
for i, s in enumerate(raw_stages):
    s["constraints"]["protein_backbone"] = bb_schedule[i]
    s["constraints"]["protein_sidechain"] = sc_schedule[i]
    s["constraints"]["ligand_ABC"] = lig_schedule[i]
    stage_dicts3.append(s)

result3 = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=stage_dicts3,
    output_name="taper_plus_ligand",
    selections={"ligand_ABC": "resname ABC"},
)
print(f"OpenMM dir:      {result3['openmm_dir']}")
print(f"Config files:    {[p.name for p in result3['config_files']]}")
print(f"Restraint files: {result3['restraint_files']}")
```

---

### Ensemble Types (OpenMM)

| Ensemble | Folder | Pressure Coupling | Use Case |
|----------|--------|------------------|----------|
| **NVT** | `01_NVT` | None (`pcouple = no`) | Initial heating / fixed box |
| **NPT** | `02_NPT` | `MonteCarloMembraneBarostat` | General membrane equilibration |
| **NPAT** | `03_NPAT` | `MonteCarloAnisotropicBarostat` (Z-free) | Anisotropic pressure control |
| **NPgT** | `04_NPgT` | `MonteCarloMembraneBarostat` + surface tension | Membrane with explicit surface tension |

---

### Important Notes

#### Force constant schedule

The CHARMM-GUI template files contain a built-in force constant schedule
(`fc_bb`, `fc_sc`, `fc_lpos`) that decreases progressively across the six equilibration steps
(4000 → 2000 → 1000 → 500 → 200 → 50 kJ/mol/nm²). These values are **hardcoded in the templates**
and are not controlled by the user's `constraints` dict.

The `constraints` dict controls **which atom types are included in the restraint index files**:

- A non-zero value → that atom type is written to `prot_pos.txt` or `lipid_pos.txt`
- A zero value → that atom type is excluded from the restraint files

#### Units

| Parameter | NAMD | OpenMM template | User API |
|-----------|------|-----------------|----------|
| Force constants | kcal/mol/Å² (B-factor) | kJ/mol/nm² (hardcoded in template) | kcal/mol/Å² (controls inclusion only) |
| Timestep | fs | ps (converted internally) | fs |
| Temperature | K | K | K |

#### Lipid dihedral restraints

`fc_ldih = 0` is hardcoded in all templates. OpenMM dihedral restraints require equilibrium
dihedral angles from the structure, which GateWizard does not compute. Only positional
restraints (`prot_pos.txt`, `lipid_pos.txt`) are used.

---

### Core Method: setup_openmm_equilibration()

```python
setup_openmm_equilibration(
    system_files: Optional[Dict[str, str]] = None,
    stage_params_list: Optional[List[Dict[str, Any]]] = None,
    output_name: str = "equilibration",
    scheme_type: Optional[str] = None,
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `system_files` | `Optional[Dict[str, str]]` | `None` | System file paths. If `None`, auto-detects in `working_dir` |
| `stage_params_list` | `List[Dict[str, Any]]` | Required | List of stage parameter dicts |
| `output_name` | `str` | `"equilibration"` | Name for output folder |
| `scheme_type` | `Optional[str]` | `None` | Equilibration scheme (NVT, NPT, NPAT, NPgT). Auto-detected from first stage if `None` |

**Returns:** `Dict[str, Any]` with keys:

- `openmm_dir` (`Path`): Output directory path
- `config_files` (`List[Path]`): Generated `.inp` files
- `run_script` (`Path`): Bash run script path

**What it does:**

1. Auto-detects system files if not provided
2. Creates `{working_dir}/{output_name}/` output directory
3. Copies system files (`prmtop`, `inpcrd`, `pdb`) into the output directory
4. Copies the 7 CHARMM-GUI Python helper scripts into the output directory
5. Generates `.inp` config files for each stage from the appropriate template
6. Generates restraint index files (`prot_pos.txt`, `lipid_pos.txt`, `dihe.txt`) in `restraints/`
7. Generates `run_equilibration.sh` bash script

---

### Stage Parameters (OpenMM)

Each entry in `stage_params_list` accepts:

| Parameter | Type | Required | Default | Description |
|-----------|------|----------|---------|-------------|
| `name` | `str` | Yes | — | Stage label (for logging) |
| `time_ns` | `float` | Yes | — | Simulation time in nanoseconds |
| `ensemble` | `str` | Yes | — | NVT, NPT, NPAT, or NPgT |
| `temperature` | `float` | Yes | — | Temperature in Kelvin |
| `timestep` | `float` | No | `2.0` | Timestep in femtoseconds |
| `minimize_steps` | `int` | No | `0` | Minimization steps (first stage only) |
| `dcd_freq` | `int` | No | `5000` | DCD output frequency (steps) |
| `constraints` | `Dict[str, float]` | No | all zero | Restraint inclusion dict (kcal/mol/Å²) |

**Constraints keys:** `protein_backbone`, `protein_sidechain`, `lipid_head`, `lipid_tail`

---

### Output Structure (OpenMM)

```
{output_name}/
├── system.prmtop              # Copied AMBER topology
├── system.inpcrd              # Copied AMBER coordinates
├── system.pdb                 # Copied system PDB
├── openmm_run.py              # CHARMM-GUI runner script
├── omm_readinputs.py          # Config parser
├── omm_readparams.py          # File I/O
├── omm_restraints.py          # Restraint forces
├── omm_barostat.py            # Pressure coupling
├── omm_vfswitch.py            # vdW force switching
├── omm_rewrap.py              # Coordinate rewrapping
├── step6.1_equilibration.inp  # Stage 1 config
├── step6.2_equilibration.inp  # Stage 2 config
├── ...
├── step7_production.inp       # Production config
├── run_equilibration.sh       # Bash run script
└── restraints/
    ├── prot_pos.txt           # Protein atom indices (BB / SC labels)
    ├── lipid_pos.txt          # Lipid head-group atom indices
    └── dihe.txt               # Dihedral restraints (always empty)
```

---

### How to Run (OpenMM)

```bash
cd {output_name}/
bash run_equilibration.sh

# Override Python interpreter:
PYTHON=python3 bash run_equilibration.sh
```

The run script chains all stages sequentially. Each stage uses `-irst <prev>.rst` to restart
from the previous stage's checkpoint, except the first stage which reads from the AMBER
coordinate file.

---

---

## EquilibrationStage Dataclass

`EquilibrationStage` is a typed, immutable-friendly dataclass that wraps the
stage-parameter dictionary used by both engines.  It is the recommended way to
build and manipulate individual stages.

### Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `name` | `str` | `"Stage"` | Human-readable label |
| `ensemble` | `str` | `"NPT"` | Thermodynamic ensemble (`NVT`, `NPT`, `NPAT`, `NPgT`) |
| `time_ns` | `float` | `0.125` | Simulation time per stage (ns) |
| `temperature` | `float` | `310.15` | Temperature (K) |
| `timestep` | `float` | `2.0` | Integration timestep (fs) |
| `pressure` | `float` | `1.0` | Target pressure (bar) |
| `minimize_steps` | `int` | `0` | Energy minimization steps (first stage only) |
| `dcd_freq` | `int` | `5000` | Trajectory save frequency (steps) |
| `constraints` | `Dict[str, float]` | `{}` | Force constants (kcal/mol/Å²) per category |

### Methods

```python
# Replace specific fields → new object (original unchanged)
new_stage = stage.replace(time_ns=0.5, temperature=298.15)

# Convert to plain dict for the manager API
d = stage.to_dict()
```

### Example

```python
from gatewizard.tools.equilibration import EquilibrationStage

s = EquilibrationStage(
    name="NVT Heating",
    ensemble="NVT",
    time_ns=0.25,
    temperature=310.15,
    timestep=2.0,
    minimize_steps=5000,
    constraints={"protein_backbone": 10.0, "protein_sidechain": 5.0},
)
# Create a NPT variant at different temperature
s_npt = s.replace(ensemble="NPT", name="NPT Stage 1", minimize_steps=0)
```

---

## Default Equilibration Protocol

Both engines expose `get_default_stage_params()` to generate a standard
CHARMM-GUI-compatible six-stage schedule with gradually decreasing positional
restraints.

### NAMD — zero-config defaults

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

manager = NAMDEquilibrationManager(working_dir=Path("/tmp/work"))

result = manager.setup_namd_equilibration(
    system_files={
        "psf":  "system.psf",
        "pdb":  "system.pdb",
        "prm":  "par_all36m_prot.prm",
    },
    # stage_params_list omitted → uses 6-stage default
)
print(result["namd_dir"])   # Path to output folder
```

### OpenMM — zero-config defaults

```python
from pathlib import Path
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

manager = OpenMMEquilibrationManager(working_dir=Path("/tmp/work"))

result = manager.setup_openmm_equilibration(
    system_files={
        "prmtop": "system.prmtop",
        "inpcrd": "system.inpcrd",
        "pdb":    "system.pdb",
    },
    # stage_params_list omitted → uses 6-stage CHARMM-GUI default
)
print(result["openmm_dir"])
```

### Retrieving the default stages

```python
# NAMD
stages = NAMDEquilibrationManager.get_default_stage_params()

# OpenMM  (NPT, NPAT, or NPgT scheme)
stages = OpenMMEquilibrationManager.get_default_stage_params(scheme_type="NPT")
```

---

## Modifying Stage Parameters

### Direct attribute assignment

```python
stages = OpenMMEquilibrationManager.get_default_stage_params()

# Run all stages at 298 K instead of 310.15 K
for s in stages:
    s["temperature"] = 298.15

# Double the simulation time for stage 6
stages[5]["time_ns"] = 0.5
```

### Using `EquilibrationStage.replace()` for immutable updates

```python
from gatewizard.tools.equilibration import EquilibrationStage

raw = OpenMMEquilibrationManager.get_default_stage_params()
stages = [EquilibrationStage(**s) for s in raw]

# Build a taper schedule: halve backbone restraint each stage
bb_forces = [10.0, 5.0, 2.5, 1.0, 0.5, 0.0]
stages = [
    s.replace(constraints={**s.constraints, "protein_backbone": bb_forces[i]})
    for i, s in enumerate(stages)
]

result = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=[s.to_dict() for s in stages],
)
```

---

## Custom Restraints — NAMD

Three levels of customisation are available, from simplest to most flexible.

### Level 1 — Override a single force constant key

Use this when you only need to change one category (e.g. turn off sidechain
restraints entirely):

```python
stages = NAMDEquilibrationManager.get_default_stage_params()
for s in stages:
    s["constraints"]["protein_sidechain"] = 0.0  # no SC restraints

result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
)
```

### Level 2 — Auto-detected MDAnalysis selections

GateWizard can auto-detect standard selections from your PDB.  Pass a custom
`selections` dict to override specific categories:

```python
result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    selections={
        "protein_backbone": "backbone and segid PROA PROB",
        "protein_sidechain": "protein and not backbone and segid PROA PROB",
        "lipid_head":        "resname POPC and name P O11 O12 O13 O14",
    },
)
```

Supported keys: `protein_backbone`, `protein_sidechain`, `lipid_head`,
`lipid_tail`, `water`, `ions`, `other`.

### Level 3 — Full MDAnalysis control with ligands

Add a custom `ligand_<RESNAME>` or `other` key to the `constraints` dict, and
provide a matching MDAnalysis selection string:

```python
stages = NAMDEquilibrationManager.get_default_stage_params()
for s in stages[:3]:                            # restrain ligand only in stages 1-3
    s["constraints"]["ligand_ABC"] = 5.0        # kcal/mol/Å²

result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    selections={
        "ligand_ABC": "resname ABC",            # MDAnalysis selection string
    },
)
# Generates: restraints/ligand_ABC_restraints.txt (PDB-format NAMD restraint file)
```

> **Requires MDAnalysis.**  Install with `conda install -c conda-forge mdanalysis`.

---

## Custom Restraints — OpenMM

OpenMM restraints use atom-index files consumed by `omm_restraints.py` at
runtime.

### Standard protein / lipid restraints

These are enabled automatically when the `constraints` dict contains non-zero
forces for the standard keys:

```python
stages = OpenMMEquilibrationManager.get_default_stage_params()
# stages already contain protein_backbone, protein_sidechain, lipid_head, etc.

result = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=stages,
)
# Generates:
#   restraints/prot_pos.txt   (BB / SC labels)
#   restraints/lipid_pos.txt  (lipid atom indices)
```

GateWizard writes `rest = yes` automatically for any stage where at least one
force constant is > 0, and injects the correct `fc_bb`, `fc_sc`, and `fc_lpos`
values (in kJ/mol/nm²) into the CHARMM-GUI template.

### Ligand / cofactor restraints (`custom_pos.txt`)

Add a `ligand_<RESNAME>` (or any non-standard) key to the `constraints` dict
and supply a matching MDAnalysis selection:

```python
stages = OpenMMEquilibrationManager.get_default_stage_params()

# Restrain ligand ABC at 5 kcal/mol/Å² in the first three stages
for s in stages[:3]:
    s["constraints"]["ligand_ABC"] = 5.0

result = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    selections={
        "ligand_ABC": "resname ABC",   # MDAnalysis selection string
    },
)
# Generates:
#   restraints/custom_pos.txt  (per-atom force constants in kJ/mol/nm²)
```

The `custom_pos.txt` file is read at runtime by the GateWizard extension block
inside `omm_restraints.py` via a `CustomExternalForce` with a periodic harmonic
potential (`k*periodicdistance(x,y,z,x0,y0,z0)²`).

### How force constants are determined

| Key | Restraint file | Force constant source |
|-----|---------------|----------------------|
| `protein_backbone` | `prot_pos.txt` (BB) | Per-stage `fc_bb` in the `.inp` config |
| `protein_sidechain` | `prot_pos.txt` (SC) | Per-stage `fc_sc` in the `.inp` config |
| `lipid_head` / `lipid_tail` | `lipid_pos.txt` | Per-stage `fc_lpos` in the `.inp` config |
| Any other key | `custom_pos.txt` | Fixed at **maximum across all stages**; baked into the file |

> **Note — `custom_pos.txt` limitation**: force constants in `custom_pos.txt`
> are fixed at the maximum value found across all stages.  Per-stage scaling of
> custom forces is not supported in the current implementation.  If you need
> different forces per stage, generate the files manually and place them in the
> `restraints/` directory before running.

### Using pre-computed atom indices (PSF + PDB)

MDAnalysis can load topology files for more accurate selections.  If you have
already computed the 0-based atom indices externally, pass them as a literal
index selection:

```python
selections = {
    "ligand_ABC": "index 1234 1235 1236 1237",  # 0-based, space-separated
}
```

Or generate the `restraints/` files manually and simply omit `selections` from
the call.

> **Requires MDAnalysis** for custom categories.
> Install with `conda install -c conda-forge mdanalysis`.

### Return value

`setup_openmm_equilibration()` now returns a `"restraint_files"` key:

```python
result = manager.setup_openmm_equilibration(...)
print(result["restraint_files"])
# {"prot_pos": PosixPath("...restraints/prot_pos.txt"),
#  "lipid_pos": None,
#  "custom_pos": PosixPath("...restraints/custom_pos.txt")}
```

---

## Class: GROMACSEquilibrationManager

Manager for GROMACS equilibration simulations using CHARMM-GUI-style MDP templates.
Mirrors the NAMD and OpenMM manager APIs.

- Accepts GROMACS-native (``.gro`` + ``topol.top``) **or** AMBER (``.prmtop`` + ``.inpcrd``) files.
  AMBER files are automatically converted using **ParmEd**.
- Position restraints are handled via MDP ``define`` macros (`POSRES_FC_BB`,
  `POSRES_FC_SC`, `POSRES_FC_LIPID`) — identical to the CHARMM-GUI scheme.
- Force constants are specified in **kcal/mol/Å²** and converted to kJ/mol/nm² internally.
- COM restraints via GROMACS Colvars (``distanceX/Y/Z`` + optional ``orientation``).

### Import

```python
from gatewizard.tools.equilibration import GROMACSEquilibrationManager
```

### Constructor

```python
GROMACSEquilibrationManager(working_dir: Path, gmx_executable: str = "gmx")
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `working_dir` | `Path` | Required | Directory containing system files |
| `gmx_executable` | `str` | `"gmx"` | GROMACS executable name or full path |

---

### `get_default_stage_params` (static)

```python
GROMACSEquilibrationManager.get_default_stage_params(
    scheme_type: str = "NPT",
    temperature: float = 310.15,
    include_production: bool = False,
) -> List[EquilibrationStage]
```

Returns 7 (or 8) CHARMM-GUI-style equilibration stages with gradually decreasing
positional restraints.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `scheme_type` | `str` | `"NPT"` | Ensemble: NVT, NPT, NPAT, or NPgT |
| `temperature` | `float` | `310.15` | Temperature in Kelvin |
| `include_production` | `bool` | `False` | Append a 50 ns unrestrained production stage |

**Returns:** `List[EquilibrationStage]`

**Example:**
```python
stages = GROMACSEquilibrationManager.get_default_stage_params("NPT", temperature=300.0)
stages[-1].time_ns = 5.0   # extend last equilibration stage
```

---

### `setup_gromacs_equilibration`

```python
manager.setup_gromacs_equilibration(
    system_files: Optional[Dict[str, str]] = None,
    stage_params_list: Optional[List] = None,
    output_name: str = "equilibration",
    scheme_type: Optional[str] = None,
    selections: Optional[Dict[str, str]] = None,
    gmx_executable: str = "gmx",
    add_com_restraint: bool = False,
    com_restraint_k: float = 10.0,
    add_rotation_restraint: bool = False,
    rotation_restraint_k: float = 2000.0,
) -> Dict[str, Any]
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `system_files` | `dict` | Auto-detect | Keys: `gro`+`top` (GROMACS) **or** `prmtop`+`inpcrd`+`pdb`+`bilayer_pdb` (AMBER) |
| `stage_params_list` | `list` | 7-stage default | List of `EquilibrationStage` objects or dicts |
| `output_name` | `str` | `"equilibration"` | Subdirectory name under `working_dir` |
| `scheme_type` | `str` | From first stage | Ensemble: NVT, NPT, NPAT, NPgT |
| `gmx_executable` | `str` | `"gmx"` | GROMACS executable |
| `add_com_restraint` | `bool` | `False` | Generate Colvars COM restraint file |
| `com_restraint_k` | `float` | `10.0` | COM force constant in kcal/mol/Å² |
| `add_rotation_restraint` | `bool` | `False` | Add orientation restraint in Colvars |
| `rotation_restraint_k` | `float` | `2000.0` | Rotation force constant in kcal/mol/Å² |

**Returns:**

```python
{
    "gromacs_dir":   Path,          # output directory
    "mdp_files":     List[Path],    # generated MDP files
    "run_script":    Path,          # run_equilibration.sh
    "system_files":  dict,          # detected input files
    "posres_files":  dict,          # backbone/sidechain/lipid .itp files
    "gro":           Path,
    "top":           Path,
    "ndx":           Path,
    "com_colvars":   Path | None,   # com_restraint.dat (if requested)
}
```

**Example — minimal (AMBER files, auto-detected):**
```python
from pathlib import Path
from gatewizard.tools.equilibration import GROMACSEquilibrationManager

manager = GROMACSEquilibrationManager(Path("popc_membrane"))
stages = GROMACSEquilibrationManager.get_default_stage_params("NPT",
                                                               include_production=True)
stages[-1].time_ns = 100.0     # 100 ns production
result = manager.setup_gromacs_equilibration(stage_params_list=stages)
# cd result["gromacs_dir"] && bash run_equilibration.sh
```

**Example — with COM restraint:**
```python
result = manager.setup_gromacs_equilibration(
    stage_params_list=stages,
    add_com_restraint=True,
    com_restraint_k=5.0,           # kcal/mol/Å²
    add_rotation_restraint=True,
    rotation_restraint_k=500.0,
)
# Activate Colvars in the desired MDP files:
#   colvars-active         = yes
#   colvars-configfile     = com_restraint.dat
```

---

### `generate_mdp_file`

```python
manager.generate_mdp_file(
    stage_name: str,
    stage_params: Dict[str, Any],
    stage_index: int,      # 0=minimization, 1-6=equilibration, 7=production
    scheme_type: str,
) -> str
```

Returns the contents of an MDP file for a single stage, generated by substituting
runtime parameters (temperature, timestep, nsteps, force constants) into the
appropriate template from `equilibration/gromacs/{ensemble}/`.

---

### `convert_from_amber`

```python
manager.convert_from_amber(
    prmtop: Path,
    inpcrd: Path,
    output_dir: Path,
    bilayer_pdb: Optional[Path] = None,
) -> Dict[str, Path]   # {"gro": ..., "top": ...}
```

Converts AMBER topology + coordinates to GROMACS format using **ParmEd**.
When the prmtop has no box information (common in membrane systems), the box
is read from the CRYST1 record of `bilayer_pdb`.

**Requires:** `conda install -c conda-forge parmed`

---

### `generate_com_colvars_config`

```python
manager.generate_com_colvars_config(
    pdb_path: Path,
    output_file: Path,
    com_restraint_k: float = 10.0,
    add_rotation_restraint: bool = False,
    rotation_restraint_k: float = 2000.0,
    selection: str = "name CA",
) -> Optional[Path]
```

Generates a Colvars configuration file that restrains the geometric centre of
the selected atoms to its initial position.  Uses ``distanceX/Y/Z`` CVs for
translation and an ``orientation`` CV for rotation.

When COM restraints are enabled through ``setup_namd_equilibration()`` or
``setup_gromacs_equilibration()``, GateWizard writes the corresponding colvars
activation block directly into the generated input files:

- NAMD: ``colvars on`` + ``colvarsConfig com_restraint.col``
- GROMACS: ``colvars-active = yes`` + ``colvars-configfile = com_restraint.dat``

**Requires:** `conda install -c conda-forge mdanalysis`

---

## Class: AmberEquilibrationManager

Manager for Amber (`pmemd` / `sander`) equilibration using GateWizard mdin
templates under `equilibration/amber/{01_NVT,02_NPT,03_NPAT,04_NPgT}/`.

Amber `01_NVT` is **true NVT** (constant volume) for all stages through
production: no `barostat` / `ntp` / surface-tension controls. Early heating in
the other ensembles remains NVT-like; from `step3` onward those packs enable
their ensemble-specific pressure settings (`02_NPT` semi-isotropic, `03_NPAT`
anisotropic Z-scaling, `04_NPgT` surface tension). Positional restraints use Amber
`ntr=1` with a GROUP block generated from MDAnalysis selections (no dihedral /
`nmropt` / `DISANG`). Stages include a separate minimization
(`step0_minimization.mdin`) plus six equilibration stages and optional production,
matching the GROMACS stage layout. Timestep ladder: **1.0 fs through equilibration 4**,
then **2.0 fs**.

Generated inputs are stamped with the GateWizard API version, local generation
time (with timezone), and the shared equilibration **templates version** so runs
remain traceable across releases.

```python
from pathlib import Path
from gatewizard.tools.equilibration import AmberEquilibrationManager

manager = AmberEquilibrationManager(Path("popc_membrane"), amber_executable="pmemd.cuda")
stages = AmberEquilibrationManager.get_default_stage_params("NPT", include_production=True)
result = manager.setup_amber_equilibration(
    stage_params_list=stages,
    amber_executable="pmemd.cuda",
)
# cd result["amber_dir"] && bash run_equilibration.sh
# Resume unfinished stages: RESUME=1 bash run_equilibration.sh
```

Executable preference when several are found: `pmemd.cuda` → `pmemd` →
`pmemd.MPI` / `pmemd.cuda.MPI` → `sander`. Minimization prefers a CPU binary when
CUDA was selected (CHARMM-GUI caution). GPU device selection uses
`CUDA_VISIBLE_DEVICES`.

### Output Structure (Amber)

```
equilibration/
├── system.prmtop / system.inpcrd / system.pdb
├── step0_minimization.mdin … step6_equilibration.mdin
├── step7_production.mdin
└── run_equilibration.sh
```

Stage outputs: `*.mdout`, `*.rst7`, `*.nc` (dynamics), `*.mdinfo`.

---

## COM Restraints (all engines)

GateWizard supports **centre-of-mass / centre-of-geometry restraints** for
NAMD, OpenMM, and GROMACS.  Unlike per-atom positional restraints, COM restraints act on
the *centroid* of the selected group — they prevent rigid-body translation (and
optionally rotation) of the protein without introducing bias on individual atoms.
Amber equilibration currently uses GROUP positional restraints only (no COM colvars).

### Motivation

In membrane simulations, proteins can undergo slow lateral drift or rotation
within the membrane plane over long equilibration runs.  COM restraints counteract
this while still allowing internal conformational flexibility.

### API summary

| Engine | Parameter | Effect |
|--------|-----------|--------|
| **NAMD** | `add_com_restraint=True` in `setup_namd_equilibration()` | Writes `com_restraint.col` and injects `colvars on` / `colvarsConfig ...` into the generated `.conf` files |
| **OpenMM** | `add_com_restraint=True` in `setup_openmm_equilibration()` | Writes `com_restraint_params.json`; `omm_restraints.py` reads it automatically |
| **GROMACS** | `add_com_restraint=True` in `setup_gromacs_equilibration()` | Writes `com_restraint.dat` and injects `colvars-active = yes` / `colvars-configfile = ...` into the generated `.mdp` files |

All engines use kcal/mol/Å² for `com_restraint_k` and `rotation_restraint_k`.

### NAMD example

```python
from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

manager = NAMDEquilibrationManager(Path("popc_membrane"))
result = manager.setup_namd_equilibration(
    add_com_restraint=True,
    com_restraint_k=5.0,
    add_rotation_restraint=True,
    rotation_restraint_k=200.0,
)
# The generated .conf files already include the colvars activation block.
```

### OpenMM example

```python
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

manager = OpenMMEquilibrationManager(Path("popc_membrane"))
result = manager.setup_openmm_equilibration(
    add_com_restraint=True,
    com_restraint_k=5.0,
)
# com_restraint_params.json is auto-loaded by omm_restraints.py
```

### GROMACS example

```python
from gatewizard.tools.equilibration import GROMACSEquilibrationManager

manager = GROMACSEquilibrationManager(Path("popc_membrane"))
stages = GROMACSEquilibrationManager.get_default_stage_params("NPT")
result = manager.setup_gromacs_equilibration(
    stage_params_list=stages,
    add_com_restraint=True,
    com_restraint_k=5.0,
)
# The generated .mdp files already include the colvars activation block.
```

---

## See Also

- [User Guide](../user-guide.md) - Complete usage guide
- [Examples](https://github.com/maurobedoya/gatewizard/tree/main/tests/analysis_examples) - Working code examples

- [Troubleshooting](../troubleshooting.md) - Common issues