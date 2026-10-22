# GateWizard Quick Reference

This page is the Python library. The desktop app is the [User Guide](user-guide.md).

The usual order is prepare the structure, build the membrane, write the equilibration, then analyze. Paste the blocks into the same Python session, in order. Each block reads a file the block before it wrote for PDB 6RV2. The linked module page has every argument.

## Open a session

Create the `gatewizard` environment once, from the [Installation](installation.md) page. Every new terminal then needs that environment turned on. The prompt shows `(gatewizard)` when it is on.

```bash
mkdir -p ~/gatewizard-work
cd ~/gatewizard-work
conda activate gatewizard
python -c "import gatewizard; print(gatewizard.__version__)"
```

Work in that folder. The first script saves `6RV2.pdb` in the current directory. Start Python and paste the next block:

```bash
python
```

Or save the block as `start_6rv2.py` and run `python start_6rv2.py`. Leave Python with `exit()`.

---

## Start from a PDB

Download [6RV2](https://www.rcsb.org/structure/6RV2) (TASK-1), keep chains A and B, and fill their short gaps. The crystal deposits two dimers, A+B and C+D. A and B are the copy with the short gaps, so C and D are left out. Later blocks can use `6RV2_AB.pdb`. The [study case](study-cases/prepare-systems.md#6rv2-task-1) runs the full membrane build. This script does not protonate or pack. Chain selection and the download are in the [Structure Manager](api/structure_manager.md). Listing and filling gaps are in [Preparation](api/preparation.md). Filling needs `pdbfixer`, which `pip install "gatewizard[full]"` provides.

```python
from gatewizard.core.structure_manager import StructureManager
from gatewizard.core.short_gaps import classify_short_gaps, fill_short_gaps

viewer = StructureManager()
info = viewer.load_from_pdb_id("6RV2", output_dir=".")
print(info["n_atoms"], "atoms,", info["n_chains"], "chains:", info["chains"])
print(info["title"])

report = classify_short_gaps("6RV2.pdb", chains=["A", "B"])
for gap in report["filled"]:
    print(f"fill {gap['count']} residues on chain {gap['chain']} after {gap['after_resid']}")
for gap in report["skipped"]:
    print(f"leave open ({gap['reason']}): chain {gap['chain']}")

filled = fill_short_gaps("6RV2.pdb", "6RV2_AB.pdb", chains=["A", "B"])
print("wrote", filled["output_file"])
print("built", len(filled["built"]), "gaps")
```

The two gaps that can be built are 3 residues each, at residue 149 on chain A and chain B. Longer stretches that exist only in SEQRES are reported and left open.

`6RV2_AB.pdb` still has the crystal phosphatidylcholine (`PC1`), decyl maltoside (`DMU`), and cholesterol hemisuccinate (`Y01`). Keep the protein and the potassium ions named `K` before PropKa. That is the same choice as the [study case](study-cases/prepare-systems.md#6rv2-task-1). Do this while the residue names are still the deposited ones. After PropKa, names such as ASH are no longer selected as protein.

```python
from gatewizard.core.structure_manager import StructureManager

solute = StructureManager()
solute.load_structure("6RV2_AB.pdb")
drop = solute.select_atoms("not protein and not resname K")
solute.delete_atoms(drop)
solute.save_pdb("6RV2_solute.pdb")
```

---

## Propka Analysis

Predict pKa values for `6RV2_solute.pdb` and write protonation names. Every argument is in [Preparation](api/preparation.md).

### Basic Analysis
```python
from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager(propka_version="3")
pka_file = analyzer.run_analysis("6RV2_solute.pdb")
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
# Apply protonation states at target pH.
# This only changes residue names (ASP→ASH). It does not add hydrogens.
stats = analyzer.apply_protonation_states(
    input_pdb="6RV2_solute.pdb",
    output_pdb="6RV2_ph7.4.pdb",
    ph=7.4
)
print(f"Modified {stats['residue_changes']} residues")

# Optional: override one residue from the list printed above.
# The key is RESNAME + number, or RESNAME + number + _ + chain.
custom_states = {}
for res in residues:
    if res["residue"] == "ASP" and res["chain"] == "A":
        custom_states[f"ASP{res['res_id']}_A"] = "ASH"
        break
if custom_states:
    stats = analyzer.apply_protonation_states(
        input_pdb="6RV2_solute.pdb",
        output_pdb="6RV2_custom.pdb",
        ph=7.4,
        custom_states=custom_states
    )
```

The membrane step uses `6RV2_ph7.4.pdb`. `6RV2_custom.pdb` is only the override.

### Generate Multiple pH States
```python
for ph in [5.0, 6.0, 7.0, 8.0]:
    stats = analyzer.apply_protonation_states(
        input_pdb="6RV2_solute.pdb",
        output_pdb=f"6RV2_ph{ph:.1f}.pdb",
        ph=ph
    )
```

---

## System Preparation

Pack a membrane around `6RV2_ph7.4.pdb` and write the Amber topology. Details are in the [Builder](api/builder.md). Packing this protein takes about 45 minutes. `wait=True` blocks until that job finishes.

### Basic Membrane System
```python
from gatewizard.core.builder import Builder

builder = Builder()
builder.set_configuration(
    salt_concentration=0.15,
    dist=12,
    dist_wat=26
)

success, msg, job_dir = builder.prepare_system(
    pdb_file="6RV2_ph7.4.pdb",
    working_dir="6RV2_system",
    upper_lipids=["POPC", "POPE"],
    lower_lipids=["POPC", "POPE"],
    lipid_ratios="3:1//3:1",
    wait=True,
)
print(success, msg)
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

Write an OpenMM equilibration for the packed 6RV2 system. OpenMM is part of `pip install "gatewizard[full]"`, which is the install on the [Installation](installation.md) page, so there is no separate program to put on `PATH`. The files are the ones under `6RV2_system/`. Every stage field is on the [Equilibration](api/equilibration.md) page.

The default membrane protocol is the real one: Equilibration 1–5 are 0.125, 0.5, 0.25, 0.5, and 1 ns, Equilibration 6 is 47.625 ns, and production is 50 ns. The block below is only for this tutorial. It sets Equilibration 6 and production to 1 ns each, like the longer of the short stages, so the run is a few nanoseconds. Delete those two assignments to keep 47.625 ns and 50 ns.

```python
from pathlib import Path
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

job_dir = Path("6RV2_system")
system_files = {
    "prmtop": str(job_dir / "system.prmtop"),
    "inpcrd": str(job_dir / "system.inpcrd"),
    "pdb": str(job_dir / "system_solv.pdb"),
}
stages = OpenMMEquilibrationManager.get_default_stage_params(
    "NPT", include_production=True
)
for stage in stages:
    if stage.name == "Equilibration 6":
        stage.time_ns = 1.0  # tutorial only; the real stage is 47.625 ns
    if stage.name == "Production":
        stage.time_ns = 1.0  # tutorial only; the real stage is 50 ns
manager = OpenMMEquilibrationManager(Path("."))
result = manager.setup_openmm_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    output_name="6RV2_equilibration",
)
print(result["run_script"])
```

Run it from that folder:

```bash
cd 6RV2_equilibration
bash run_equilibration.sh
```

---

## Analysis

Read the OpenMM log after `bash run_equilibration.sh` has been run in `6RV2_equilibration`. The [Analysis](api/analysis.md) page covers trajectories and the other engines.

```python
from pathlib import Path
from gatewizard.utils.openmm_analysis import get_equilibration_progress

progress = get_equilibration_progress(Path("6RV2_equilibration"))
stage = progress["equilibration_1"]
print(stage.status)
if stage.timing:
    print(stage.timing.steps_completed, "of", stage.timing.total_steps, "steps")
```

---

## Error Handling

Wrap a call so a missing file or a preparation error prints a message instead of stopping the session.

### Standard Pattern
```python
from gatewizard.core.preparation import PreparationError, PreparationManager

try:
    analyzer = PreparationManager()
    pka_file = analyzer.run_analysis("6RV2_solute.pdb")
except FileNotFoundError as e:
    print(f"File not found: {e}")
except PreparationError as e:
    print(f"Preparation error: {e}")
except Exception as e:
    print(f"Unexpected error: {e}")
```

---

## Common Dictionary Keys

These are the fields of each residue printed from the 6RV2 PropKa list, and of one equilibration stage.

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
    "constraints": {               # Restraints (kcal/mol/Å²)
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

These are the files the calls above leave on disk.

### After Propka
```
6RV2_solute.pdb                    # Protein plus potassium, detergents removed
6RV2_solute.pka                    # Full PropKa output
6RV2_solute_summary_of_prediction.txt
6RV2_ph7.4.pdb                     # Names at pH 7.4; the membrane step uses this
```

### After System Preparation
```
6RV2_system/
  ├── system.prmtop          # Topology
  ├── system.inpcrd          # Coordinates
  ├── system_solv.pdb        # Solvated structure
  └── step3_input.pdb        # Pre-solvation
```

### After Equilibration Setup
```
6RV2_equilibration/
  ├── step1_equilibration.inp     # Config for stage 1
  ├── step1_equilibration.log     # Written when that stage is run
  ├── openmm_run.py               # OpenMM runner
  ├── run_equilibration.sh        # Run script
  ├── system.prmtop
  ├── system.inpcrd
  └── system.pdb
```

---

## Tips

Small choices that keep a first run from stalling or writing into the wrong folder.

### Performance
- Use `use_gpu=True` for equilibration stages
- Increase `cpu_cores` for larger systems
- Use `timestep=1.0` for heavily restrained stages
- Use `timestep=2.0` for production

### File Management
- Stay in `~/gatewizard-work` so the `6RV2_...` paths above keep resolving
- Use `pathlib.Path` for file operations
- Create output directories before writing

### Debugging
- Check log files: `tail -f step*.log`
- Verify file existence before processing
- Use try-except for robust error handling
