# GateWizard Quick Reference

Quick reference for common GateWizard API operations.

---

## Propka Analysis

### Basic Analysis
```python
from gatewizard.core.propka import PropkaAnalyzer

analyzer = PropkaAnalyzer(propka_version="3")
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)
```

### Access Results
```python
for res in residues:
    print(f"{res['residue']}{res['res_id']} Chain {res['chain']}: pKa={res['pka']:.2f}")
```

### Filter by Type
```python
acidic = [r for r in residues if r['residue'] in ['ASP', 'GLU']]
basic = [r for r in residues if r['residue'] in ['LYS', 'ARG', 'HIS']]
```

### Write PDB with Protonation States
```python
# Apply protonation states at target pH
# NOTE: This only changes residue names (ASP→ASH), does NOT add hydrogens
stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="protein_ph7.4.pdb",
    ph=7.4
)
print(f"Modified {stats['residue_changes']} residues")

# For MD simulations, add hydrogens with pdb4amber:
# pdb4amber -i protein_ph7.4.pdb -o protein_prepared.pdb

# With custom states
custom_states = {
    "HIS67": "HID",      # Delta-protonated histidine
    "ASP42_A": "ASH"     # Chain-specific override
}
stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="protein_custom.pdb",
    ph=7.4,
    custom_states=custom_states
)
```

### Generate Multiple pH States
```python
for ph in [5.0, 6.0, 7.0, 8.0]:
    stats = analyzer.apply_protonation_states(
        input_pdb="protein.pdb",
        output_pdb=f"protein_ph{ph:.1f}.pdb",
        ph=ph
    )
```

---

## System Preparation

### Basic Membrane System
```python
from gatewizard.core.system_builder import SystemBuilder

builder = SystemBuilder()
builder.set_configuration(
    salt_concentration=0.15,
    dist_wat=17.5
)

success, msg, job_dir = builder.prepare_system(
    pdb_file="protein.pdb",
    working_dir="system_prep",
    upper_lipids=["POPC", "POPE"],
    lower_lipids=["POPC", "POPE"],
    lipid_ratios="3:1//3:1"
)
```

### Check Output
```python
if success:
    topology = job_dir / "system.prmtop"
    coords = job_dir / "system.inpcrd"
    pdb = job_dir / "system_solv.pdb"
```

---

## Equilibration Setup

### Basic Protocol
```python
from gatewizard.tools.equilibration import NAMDEquilibrationManager
from pathlib import Path

manager = NAMDEquilibrationManager(
    working_dir=Path("equilibration"),
    namd_executable="namd3"
)

protocol = {
    "Equilibration 1": {
        "name": "Equilibration 1",
        "time_ns": 0.125,
        "timestep": 1.0,
        "temperature": 310.15,
        "pressure": 1.0,
        "constraints": {
            "protein_backbone": 10.0,
            "protein_sidechain": 5.0,
            "lipid_head": 2.0,
            "water": 0.0
        },
        "minimize_steps": 10000,
        "use_gpu": True,
        "cpu_cores": 4,
        "dcd_freq": 5000
    }
}
```

### Generate Configs
```python
system_files = {
    'prmtop': 'system.prmtop',
    'inpcrd': 'system.inpcrd',
    'pdb': 'system_solv.pdb'
}

for stage_key, stage_data in protocol.items():
    stage_index = list(protocol.keys()).index(stage_key)
    
    config = manager.generate_charmm_gui_config_file(
        stage_name=stage_key,
        stage_params=stage_data,
        stage_index=stage_index,
        system_files=system_files
    )
    
    config_name = manager._get_config_name(stage_key)
    with open(f"{config_name}.conf", 'w') as f:
        f.write(config)
```

### Generate Run Script
```python
script = manager.generate_run_script(protocol, "namd3")
with open("run_equilibration.sh", 'w') as f:
    f.write(script)
Path("run_equilibration.sh").chmod(0o755)
```

---

## Analysis

### Monitor Progress
```python
from gatewizard.utils.namd_analysis import get_equilibration_progress
from pathlib import Path

progress = get_equilibration_progress(Path("equilibration"))
print(f"Progress: {progress['percent_complete']:.1f}%")
```

### Parse Log File
```python
from gatewizard.utils.namd_analysis import parse_namd_log

log_data = parse_namd_log("step1_equilibration.log")
print(f"Avg Temperature: {sum(log_data['temperature'])/len(log_data['temperature']):.2f} K")
```

---

## Error Handling

### Standard Pattern
```python
from gatewizard.core.propka import PropkaError

try:
    analyzer = PropkaAnalyzer()
    pka_file = analyzer.run_analysis("protein.pdb")
except FileNotFoundError as e:
    print(f"File not found: {e}")
except PropkaError as e:
    print(f"Propka error: {e}")
except Exception as e:
    print(f"Unexpected error: {e}")
```

---

## Common Dictionary Keys

### Propka Results
```python
residue = {
    "residue": "ASP",      # Residue name
    "res_id": 42,          # Residue number
    "chain": "A",          # Chain ID
    "pka": 3.85,           # Predicted pKa
    "model_pka": 3.80      # Model pKa
}
```

### Stage Parameters
```python
stage = {
    "time_ns": 0.125,              # Simulation time
    "timestep": 1.0,               # Timestep in fs
    "temperature": 310.15,         # Temperature in K
    "pressure": 1.0,               # Pressure in atm
    "constraints": {               # Restraints (kcal/mol/Ų)
        "protein_backbone": 10.0,
        "protein_sidechain": 5.0,
        "lipid_head": 2.0,
        "water": 0.0
    },
    "use_gpu": True,
    "cpu_cores": 4,
    "dcd_freq": 5000
}
```

---

## File Locations

### After Propka
```
protein.pka                        # Full output
protein_summary_of_prediction.txt  # Summary section
```

### After System Preparation
```
job_directory/
  ├── system.prmtop          # Topology
  ├── system.inpcrd          # Coordinates
  ├── system_solv.pdb        # Solvated structure
  └── step3_input.pdb        # Pre-solvation
```

### After Equilibration Setup
```
equilibration/
  ├── step1.conf                  # Config for stage 1
  ├── step2.conf                  # Config for stage 2
  ├── step7_production.conf       # Production config
  ├── run_equilibration.sh        # Run script
  └── protocol_summary.json       # Protocol definition
```

---

## Tips

### Performance
- Use `use_gpu=True` for equilibration stages
- Increase `cpu_cores` for larger systems
- Use `timestep=1.0` for heavily restrained stages
- Use `timestep=2.0` for production

### File Management
- Use `pathlib.Path` for file operations
- Create output directories before writing
- Use absolute paths when possible

### Debugging
- Check log files: `tail -f step*.log`
- Verify file existence before processing
- Use try-except for robust error handling
