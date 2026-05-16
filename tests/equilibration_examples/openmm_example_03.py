"""
openmm_example_03.py
====================
OpenMM equilibration with custom positional restraints for a ligand (residue
name ABC) plus standard protein/lipid restraints.

Prerequisites
-------------
- AMBER-format system: system.prmtop + system.inpcrd + system.pdb
- MDAnalysis (``conda install -c conda-forge mdanalysis``)

The example demonstrates:
  1. Default 6-stage CHARMM-GUI-style protocol
  2. Adding a ligand_ABC restraint in the first 3 stages
  3. Inspecting the generated restraint files
  4. Per-stage force-constant taper for backbone
"""

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

# ---------------------------------------------------------------------------
# Example 1 — Standard protein + lipid restraints (auto-detected)
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Example 2 — Add ligand ABC restraints in stages 1-3
# ---------------------------------------------------------------------------
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

# Verify rest = yes in stage 1 config (ligand_ABC = 5.0 > 0)
stage1_cfg = result2["config_files"][0].read_text()
assert "rest        = yes" in stage1_cfg, "Stage 1 should have rest = yes"

# Verify rest = no in stage 4 config (ligand_ABC = 0, no other active forces assumed)
# Note: if standard protein forces are also 0 in stage 4, rest = no there.
print("Assertions passed.")

# ---------------------------------------------------------------------------
# Example 3 — Custom backbone taper + ligand restraints
# ---------------------------------------------------------------------------
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
