# Equilibration Module

Module for setting up NAMD equilibration protocols for membrane protein systems. Generates configuration files, restraint files, and run scripts for multi-stage equilibration simulations using AMBER force fields.

## Import

```python
from gatewizard.tools.equilibration import NAMDEquilibrationManager
```

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

### Method: generate_restraints_file()

Generates restraint PDB files with B-factors encoding restraint forces for each atom type.

```python
generate_restraints_file(
    system_pdb: Path,
    constraints: Dict[str, float],
    output_file: Path,
    stage_name: str = ""
) -> None
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `system_pdb` | `Path` | Path to the **system.pdb** file (full system: protein + lipids + water + ions) |
| `constraints` | `Dict[str, float]` | Dictionary of restraint forces (kcal/mol/Å²) for each atom type |
| `output_file` | `Path` | Path for output restraints PDB file |
| `stage_name` | `str` | Stage name for logging (optional) |

**Behavior:**

Reads the **system.pdb** file (NOT protein.pdb - includes all atoms)

A. Classifies each atom by type:

   - `protein_backbone`: CA, C, N, O atoms in protein residues
   - `protein_sidechain`: Heavy sidechain atoms in protein residues
   - `lipid_head`: P, O atoms in lipid head groups
   - `lipid_tail`: C atoms in lipid tails
   - `water`: H2O molecules (TIP3, HOH, WAT, SOL)
   - `ions`: Na+, Cl-, K+, Ca2+, Mg2+, etc.
   - `other`: Ligands and other molecules

B. Assigns B-factor values based on constraints dictionary

C. Writes restraints PDB file with modified B-factors

D. Logs detailed statistics: total atoms, counts per type, force values


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
4. Templates loaded from: `charmm_gui_templates/{scheme_folder}/{template_file}`

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

## See Also

- [User Guide](../user-guide.md) - Complete usage guide
- [Examples](../../tests/analysis_examples/) - Working code examples
- [Troubleshooting](../troubleshooting.md) - Common issues