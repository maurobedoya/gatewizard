# Equilibration Examples# Equilibration Examples# Equilibration Examples# Equilibration Examples# Equilibration Examples# Equilibration Examples - Simplified API



This directory contains working examples demonstrating the equilibration API with various features and protocols.



## Available ExamplesSimplified NAMD equilibration setup examples using the one-call API.



### Basic Examples



1. **equilibration_example_01.py** - Auto-detection from folder (NVT single stage)## Examples OverviewThis directory contains examples demonstrating the simplified NAMD equilibration API.

   - Demonstrates automatic file detection without specifying system_files

   - Uses single NVT equilibration stage

   - Shows scheme_type auto-detection from ensemble field

### Example 01: From System Folder ⭐ RECOMMENDED

2. **equilibration_example_02.py** - Explicit file paths (NPT single stage)

   - Demonstrates explicit system_files specification**File:** `equilibration_example_01_from_folder.py`

   - Uses single NPT equilibration stage

   - Alternative approach when auto-detection doesn't work## Examples OverviewThis directory contains examples demonstrating NAMD equilibration setup using the Gatewizard API.



3. **equilibration_example_03.py** - Helper method demoBest starting point - demonstrates:

   - Shows using find_system_files() to check detected files

   - Useful for debugging file detection issues- Reading files from System Builder output folder (`popc_membrane/`)



### Advanced Examples- Using bilayer PDB file for CRYST1 box dimensions



4. **equilibration_example_04.py** - Mixed ensemble protocol (NVT → NVT → NPT)- 3-stage equilibration (NVT → NPT → NPAT)All examples demonstrate the one-call `setup_namd_equilibration()` method that replicates the GUI workflow.

   - Demonstrates per-stage ensemble selection

   - Stage 3 switches to NPT ensemble- Complete workflow from folder to equilibration setup

   - Shows warning system for ensemble changes

   - **Features**: Custom stage names automatically mapped to sequential steps



5. **equilibration_example_05.py** - Multi-ensemble progressive protocol### Example 02: Basic Single Stage

   - Custom stage names: "Initial Equilibration", "Pressure Equilibration", etc.

   - Progressive ensemble changes: NVT → NPT → NPAT → NPAT**File:** `equilibration_example_02_basic.py`### Example 01: Setup from System Folder ⭐ RECOMMENDED## ExamplesThis directory contains examples demonstrating NAMD equilibration setup using the Gatewizard API.This directory contains examples of using the simplified NAMD equilibration API.

   - Shows proper restart file chaining with custom names

   - **Features**: Per-stage ensemble detection, custom naming, restart chaining



6. **equilibration_example_06.py** - Complete 7-stage CHARMM-GUI protocolMinimal example:**File:** `equilibration_example_01_from_folder.py`

   - Full 6-stage equilibration + production

   - Standard CHARMM-GUI naming convention- Single equilibration stage

   - Complete restraint progression schedule

- NPT ensemble

### Special Features

- Basic restraints with minimization

7. **equilibration_example_custom_template.py** - Explicit template control

   - Demonstrates using `custom_template` parameter**Best starting point** - demonstrates:

   - Skip intermediate stages (step6.1 → step6.3 → step6.5)

   - Mix templates from different equilibration phases### Example 03: Three-Stage Protocol

   - **Features**: Advanced template selection for custom workflows

**File:** `equilibration_example_03_three_stage.py`- Automatic file discovery from system folder1. **equilibration_example_01_basic.py** - Single stage equilibration (minimal example)

8. **test_auto_detect.py** - Testing helper for auto-detection

   - Minimal example for testing file detection

   - Useful for debugging

Standard CHARMM-GUI workflow:- Reading CRYST1 box dimensions from PDB

## Key Features Demonstrated

- Progressive restraint reduction

### 1. Auto-Detection

- **Examples 01, 03, 04, 05, 06**: Automatic detection of system files- NVT → NPT ensemble progression- Using existing system preparation output (popc_membrane/)2. **equilibration_example_02_three_stage.py** - Three-stage equilibration protocol

- No need to specify prmtop, inpcrd, pdb, bilayer_pdb paths

- Just point to the folder containing system files- 3 equilibration stages



### 2. Scheme Type Auto-Detection- 3-stage equilibration (NVT → NPT → NPAT)

- **All examples**: scheme_type automatically extracted from first stage's ensemble

- No need to specify scheme_type parameter explicitly### Example 04: Custom Four-Stage

- Backwards compatible - can still specify if needed

**File:** `equilibration_example_04_custom.py`3. **equilibration_example_03_custom.py** - Custom four-stage protocol with gradual relaxation## Examples## New Simplified API

### 3. Per-Stage Ensemble Selection

- **Examples 04, 05**: Each stage can use different ensemble (NVT, NPT, NPAT, NPgT)

- System warns when stage ensemble differs from protocol default

- Automatically selects correct CHARMM-GUI template for each ensembleAdvanced custom protocol:This example uses the complete `popc_membrane/` folder with all system files.



### 4. Custom Stage Names- Custom ensemble progression (NVT → NPT → NPAT)

- **Example 05**: Use descriptive names like "Initial Equilibration"

- Automatically mapped to sequential step numbers (step1, step2, step3)- Timestep increase (1.0 → 2.0 fs)4. **equilibration_example_04_complete.py** - Complete CHARMM-GUI 6-stage + production protocol

- Restart files properly chained regardless of naming

- Gradual restraint reduction

### 5. Custom Template Selection

- **Example custom_template**: Explicitly specify which CHARMM-GUI template to use- 4 stages, 3.75 ns total### Example 02: Basic Single Stage

- Use `custom_template: 'step6.3_equilibration.inp'` in stage parameters

- Skip stages, repeat stages, or mix templates



## Running Examples### Example 05: Complete CHARMM-GUI Protocol**File:** `equilibration_example_02_basic.py`



### Basic Usage**File:** `equilibration_example_05_complete.py`



```bash

cd tests/equilibration_examples

python equilibration_example_01.pyFull production workflow:

```

- 6 equilibration stagesMinimal example showing:## Requirements

### Check Generated Files

- 1 production stage (10 ns)

Each example creates output in a folder structure:

```- NPAT ensemble for membrane systems- Single equilibration stage

equilibration_<name>/

└── namd/- Standard CHARMM-GUI restraint reduction

    ├── step1_equilibration.conf

    ├── step2_equilibration.conf- 11.5 ns total simulation time- NPT ensemble1. **equilibration_example_01_basic.py** - Single stage equilibration (minimal example)The `setup_namd_equilibration()` method provides a complete one-call solution that replicates the GUI workflow:

    ├── step3_equilibration.conf

    ├── run_equilibration.sh

    ├── protocol_summary.json

    ├── system.prmtop---- Basic restraints with minimization

    ├── system.inpcrd

    ├── system.pdb

    ├── bilayer_protein_protonated_prepared_lipid.pdb

    └── restraints/## CRITICAL: Bilayer PDB Requirement- **PDB file required for CRYST1 box info**- System files must be present in this directory:

        ├── step1_equilibration_restraints.pdb

        ├── step2_equilibration_restraints.pdb

        └── step3_equilibration_restraints.pdb

```**All examples require a bilayer PDB file** with CRYST1 box dimensions.



### Run NAMD Simulation



```bash### Why Bilayer PDB is Required### Example 03: Three-Stage Protocol  - `system.prmtop` (AMBER topology)2. **equilibration_example_02_three_stage.py** - Three-stage equilibration protocol

cd equilibration_<name>/namd

./run_equilibration.sh

```

1. **CRYST1 record** provides accurate periodic box dimensions**File:** `equilibration_example_03_three_stage.py`

## Example Comparison

2. **NAMD requires** correct box size for MD simulations

| Example | Stages | Ensembles | Features |

|---------|--------|-----------|----------|3. **System Builder generates** `bilayer_*_lipid.pdb` with CRYST1  - `system.inpcrd` (AMBER coordinates)

| 01 | 1 | NVT | Auto-detection, single stage |

| 02 | 1 | NPT | Explicit files, single stage |4. **system.pdb does NOT contain CRYST1** - only the bilayer PDB has it

| 03 | 1 | - | Helper method demo |

| 04 | 3 | NVT,NVT,NPT | Mixed ensembles, warnings |Standard protocol demonstrating:

| 05 | 4 | NVT,NPT,NPAT,NPAT | Custom names, progressive protocol |

| 06 | 7 | NPAT (all) | Complete CHARMM-GUI protocol |### CRYST1 Record Example

| custom_template | 3 | NPT,NPT,NPAT | Explicit template control |

- Progressive restraint reduction  - `system.pdb` (System PDB file)3. **equilibration_example_03_custom.py** - Custom four-stage protocol with gradual relaxation```python

## Common Parameters

```

All examples use these common stage parameters:

CRYST1   72.440   72.700   85.450  90.00  90.00  90.00 P 1           1- NVT → NPT ensemble progression

```python

{         |        |        |        |      |      |

    'name': 'Stage Name',              # Can be anything

    'ensemble': 'NPT',                 # NVT, NPT, NPAT, or NPgT         a (Å)   b (Å)    c (Å)    α      β      γ- Typical CHARMM-GUI workflow

    'time_ns': 0.125,                  # Simulation time in ns

    'timestep': 1.0,                   # Timestep in fs (1.0 or 2.0)```

    'temperature': 310.15,             # Temperature in K

    'pressure': 1.0,                   # Pressure in bar (NPT/NPAT/NPgT only)- **PDB file required for CRYST1 box info**

    'surface_tension': 0.0,            # Surface tension (NPAT only)

    'minimize_steps': 10000,           # Minimization steps (first stage only)### Required system_files Dictionary

    'constraints': {                   # Restraint forces in kcal/mol/Å²

        'protein_backbone': 10.0,## Running Examples4. **equilibration_example_04_complete.py** - Complete CHARMM-GUI 6-stage + production protocolfrom pathlib import Path

        'protein_sidechain': 5.0,

        'lipid_head': 2.5,```python

        'lipid_tail': 2.5,

        'water': 0.0,system_files = {### Example 04: Custom Four-Stage Protocol

        'ions': 10.0,

        'other': 0.0    'prmtop': 'path/to/system.prmtop',

    },

    # Optional: Explicit template selection    'inpcrd': 'path/to/system.inpcrd',**File:** `equilibration_example_04_custom.py`

    'custom_template': 'step6.3_equilibration.inp'

}    'pdb': 'path/to/system.pdb',

```

    'bilayer_pdb': 'path/to/bilayer_protein_protonated_prepared_lipid.pdb'  # REQUIRED!

## Tips

}

1. **Start with Example 01** - Simplest working example

2. **Use Example 04** for mixed ensemble protocols```Advanced protocol showing:```bashfrom gatewizard.tools.equilibration import NAMDEquilibrationManager

3. **Use Example 05** for custom stage names with progressive protocols

4. **Use Example custom_template** for advanced template control

5. **Always check logs** - Warnings indicate ensemble changes

6. **Verify files** - Check that both system.pdb and bilayer PDB exist**Without the bilayer PDB, equilibration will fail or use incorrect box dimensions.**- Custom ensemble progression (NVT → NPT → NPAT)



## Troubleshooting



### Auto-detection fails---- Timestep increase (1.0 → 2.0 fs)python equilibration_example_01_basic.py

- Use Example 03 to check which files are detected

- Use Example 02 with explicit file paths



### Custom names not working## System Files Location- Gradual restraint reduction

- Check Example 05 - custom names are mapped to step1, step2, etc.

- Verify restart file chaining: `grep "set inputname" step2_equilibration.conf`



### Mixed ensembles not workingAll examples use files from the `popc_membrane/` folder:- Longer simulation timespython equilibration_example_02_three_stage.py## Requirements

- Check Example 04 - warnings should be logged for ensemble changes

- Each stage uses its own ensemble's template folder



### Template selection issues```- **PDB file required for CRYST1 box info**

- Use `custom_template` parameter (Example custom_template)

- Available templates: step6.1 through step6.6, step7_productionpopc_membrane/



## Additional Resources├── system.prmtop                          # AMBER topologypython equilibration_example_03_custom.py



- See `docs/api/equilibration.md` for complete API documentation├── system.inpcrd                          # AMBER coordinates  

- See `ENSEMBLE_PER_STAGE_FIX.md` for technical details on per-stage ensemble detection

- See `SCHEME_TYPE_AUTO_DETECTION.md` for scheme_type auto-detection implementation├── system.pdb                             # System PDB### Example 05: Complete CHARMM-GUI Protocol


└── bilayer_protein_protonated_prepared_lipid.pdb  # CRYST1 box info (REQUIRED!)

```**File:** `equilibration_example_05_complete.py`python equilibration_example_04_complete.py# Setup



---



## Running ExamplesFull production workflow:```



```bash- 6 equilibration stages

# Run folder-based example (recommended)

python equilibration_example_01_from_folder.py- 1 production stage (10 ns)- System files must be present in this directory:manager = NAMDEquilibrationManager(Path("/path/to/work"))



# Run basic single stage- NPAT ensemble for membranes

python equilibration_example_02_basic.py

- Standard CHARMM-GUI restraint reductionEach example creates an output directory with NAMD configuration files ready to run.

# Run three-stage protocol

python equilibration_example_03_three_stage.py- **PDB file required for CRYST1 box info**



# Run custom four-stage  - `system.prmtop` (AMBER topology)

python equilibration_example_04_custom.py

## Important: PDB File Requirement

# Run complete CHARMM-GUI protocol

python equilibration_example_05_complete.py## API Usage

```

**All examples require a PDB file** for reading CRYST1 box dimensions. The system_files dictionary must include:

---

  - `system.inpcrd` (AMBER coordinates)# Define system files (source paths)

## Output Structure

```python

Each example creates:

system_files = {All examples use the simplified `setup_namd_equilibration()` method:

```

output_name/    'prmtop': 'path/to/system.prmtop',

└── namd/

    ├── system.prmtop    'inpcrd': 'path/to/system.inpcrd',  - `system.pdb` (System PDB file)system_files = {

    ├── system.inpcrd

    ├── system.pdb    'pdb': 'path/to/system.pdb'  # REQUIRED for CRYST1!

    ├── step1_equilibration.conf

    ├── step2_equilibration.conf}```python

    ├── ...

    ├── run_equilibration.sh```

    ├── protocol_summary.json

    └── restraints/from pathlib import Path    'prmtop': '/path/to/system.prmtop',

        ├── step1_equilibration_restraints.pdb

        └── ...The PDB file provides essential box dimension information that cannot be reliably extracted from topology/coordinate files alone.

```

from gatewizard.tools.equilibration import NAMDEquilibrationManager

To run equilibration:

## Requirements

```bash

cd output_name/namd/## Running Examples    'inpcrd': '/path/to/system.inpcrd',

./run_equilibration.sh

```System files must be present:



---- `system.prmtop` - AMBER topology filemanager = NAMDEquilibrationManager(Path("/work/dir"))



## API Quick Reference- `system.inpcrd` - AMBER coordinate file



```python- `system.pdb` - PDB structure file with CRYST1 record    'pdb': '/path/to/system.pdb'

from pathlib import Path

from gatewizard.tools.equilibration import NAMDEquilibrationManager



# Create managerFor Example 01, these are provided in the `popc_membrane/` folder.result = manager.setup_namd_equilibration(

manager = NAMDEquilibrationManager(working_dir=Path("work"))



# Setup equilibration

result = manager.setup_namd_equilibration(## Running Examples    system_files={```bash}

    system_files={

        'prmtop': 'system.prmtop',

        'inpcrd': 'system.inpcrd',

        'pdb': 'system.pdb',```bash        'prmtop': 'path/to/system.prmtop',

        'bilayer_pdb': 'bilayer_protein_protonated_prepared_lipid.pdb'  # REQUIRED

    },# Run folder-based example (recommended)

    stage_params_list=[stage1, stage2, ...],

    output_name="equilibration",python equilibration_example_01_from_folder.py        'inpcrd': 'path/to/system.inpcrd',python equilibration_example_01_basic.py

    scheme_type="NPAT"

)

```

# Run basic single stage        'pdb': 'path/to/system.pdb'

---

python equilibration_example_02_basic.py

## See Also

    },python equilibration_example_02_three_stage.py# Define equilibration stages

- [Equilibration API Documentation](../../docs/api/equilibration.md)

- [System Builder](../../docs/api/system-builder.md)# Run three-stage protocol


python equilibration_example_03_three_stage.py    stage_params_list=[{...}, {...}],



# Run custom four-stage protocol    output_name="equilibration",python equilibration_example_03_custom.pystages = [

python equilibration_example_04_custom.py

    scheme_type="NPT",

# Run complete CHARMM-GUI protocol

python equilibration_example_05_complete.py    namd_executable="namd3"python equilibration_example_04_complete.py    {

```

)

## Output Structure

``````        'name': 'Equilibration 1',

Each example creates an output directory:



```

output_name/See the individual examples for complete working code and the [documentation](../../docs/api/equilibration.md) for detailed API reference.        'time_ns': 0.125,

└── namd/

    ├── system.prmtop              # Copied topology

    ├── system.inpcrd              # Copied coordinatesEach example creates an output directory with NAMD configuration files ready to run.        'steps': 125000,

    ├── system.pdb                 # Copied PDB

    ├── step1_equilibration.conf   # Stage 1 config        'ensemble': 'NVT',

    ├── step2_equilibration.conf   # Stage 2 config

    ├── ...## API Usage        'temperature': 303.15,

    ├── run_equilibration.sh       # Executable run script

    ├── protocol_summary.json      # Protocol metadata        'timestep': 1.0,

    └── restraints/

        ├── step1_equilibration_restraints.pdbAll examples use the simplified `setup_namd_equilibration()` method:        'minimize_steps': 10000,

        ├── step2_equilibration_restraints.pdb

        └── ...        'constraints': {

```

```python            'protein_backbone': 10.0,

## System Folder Structure (Example 01)

from pathlib import Path            'protein_sidechain': 5.0,

The `popc_membrane/` folder contains a complete system preparation:

from gatewizard.tools.equilibration import NAMDEquilibrationManager            'lipid_head': 2.5,

```

popc_membrane/            'lipid_tail': 2.5,

├── system.prmtop              # AMBER topology

├── system.inpcrd              # AMBER coordinatesmanager = NAMDEquilibrationManager(Path("/work/dir"))            'water': 0.0,

├── system.pdb                 # System PDB with CRYST1

├── bilayer_protein_protonated_prepared_lipid.pdb  # Full bilayer PDB            'ions': 10.0,

└── ... (other preparation files)

```result = manager.setup_namd_equilibration(            'other': 0.0



## API Quick Reference    system_files={        }



```python        'prmtop': 'path/to/system.prmtop',    },

from pathlib import Path

from gatewizard.tools.equilibration import NAMDEquilibrationManager        'inpcrd': 'path/to/system.inpcrd',    # Add more stages...



# Create manager        'pdb': 'path/to/system.pdb']

manager = NAMDEquilibrationManager(working_dir=Path("path/to/work"))

    },

# Setup equilibration

result = manager.setup_namd_equilibration(    stage_params_list=[{...}, {...}],# One method call does everything!

    system_files={

        'prmtop': 'system.prmtop',    output_name="equilibration",result = manager.setup_namd_equilibration(

        'inpcrd': 'system.inpcrd',

        'pdb': 'system.pdb'  # REQUIRED    scheme_type="NPT",    system_files=system_files,

    },

    stage_params_list=[stage1, stage2, ...],    namd_executable="namd3"    stage_params_list=stages,

    output_name="equilibration",

    scheme_type="NPAT")    output_name="equilibration",

)

``````    scheme_type="NPT",



## See Also    namd_executable="namd3"



- [Equilibration API Documentation](../../docs/api/equilibration.md)See the individual examples for complete working code.)

- [System Builder Module](../../docs/api/system-builder.md)


# Ready to run!
# cd result['namd_dir']
# ./run_equilibration.sh
```

## What It Does

The `setup_namd_equilibration()` method automatically:

1. **Creates directory structure**: `equilibration/namd/` and `restraints/`
2. **Copies system files**: Copies prmtop, inpcrd, pdb to NAMD directory
3. **Generates config files**: Creates NAMD .conf files for all stages
4. **Generates restraints**: Creates restraint .pdb files for each stage
5. **Creates run script**: Generates executable `run_equilibration.sh`
6. **Creates summary**: Generates JSON summary of protocol

All files use **relative paths** for portability.

## Examples

### Quick Start
- `equilibration_example_simplified.py` - Minimal 3-stage setup
- `equilibration_example_02_new.py` - Single stage basic example

### Complete Protocols
- `equilibration_example_charmm_gui_complete.py` - Standard 6-stage + production

### Advanced
- `equilibration_example_custom.py` - Custom stage configurations

## Scheme Types

- **NVT**: Constant volume and temperature
- **NPT**: Constant pressure and temperature
- **NPAT**: Constant surface tension (membrane systems)
- **NPgT**: Constant surface tension with pressure coupling

## Stage Parameters

Each stage dictionary can include:

```python
{
    'name': 'Equilibration 1',           # Stage name
    'time_ns': 0.125,                    # Simulation time (ns)
    'steps': 125000,                     # Number of steps
    'ensemble': 'NVT',                   # NVT, NPT, NPAT, NPgT
    'temperature': 303.15,               # Temperature (K)
    'pressure': 1.0,                     # Pressure (atm) - for NPT/NPAT/NPgT
    'surface_tension': 0.0,              # Surface tension (dyn/cm) - for NPAT/NPgT
    'timestep': 1.0,                     # Timestep (fs)
    'minimize_steps': 10000,             # Minimization steps (first stage only)
    'constraints': {                     # Restraint forces (kcal/mol/Å²)
        'protein_backbone': 10.0,
        'protein_sidechain': 5.0,
        'lipid_head': 2.5,
        'lipid_tail': 2.5,
        'water': 0.0,
        'ions': 10.0,
        'other': 0.0
    }
}
```

## Return Value

The method returns a dictionary with paths:

```python
{
    'output_dir': Path('equilibration/'),
    'namd_dir': Path('equilibration/namd/'),
    'config_files': [Path('step1_equilibration.conf'), ...],
    'restraints_dir': Path('equilibration/namd/restraints/'),
    'run_script': Path('run_equilibration.sh'),
    'summary_file': Path('protocol_summary.json')
}
```

## Running Equilibration

After setup:

```bash
cd equilibration/namd/
./run_equilibration.sh
```

The script runs all stages sequentially with proper restart handling.