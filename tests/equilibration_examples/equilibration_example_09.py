"""
equilibration_example_09.py
===========================
NAMD custom restraints — three levels of customisation.

Prerequisites
-------------
- A prepared membrane-protein system in the ``popc_membrane/`` test folder
- MDAnalysis (``conda install -c conda-forge mdanalysis``) for levels 2 and 3

"""

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
