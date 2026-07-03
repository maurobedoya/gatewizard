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
    charges={"AAA": 0, "BBB": 0},
    charge_method="bcc",
    atom_type="gaff2",  # GAFF2 atom types
)

print(f"  Parametrized: {list(ligand_results.keys())}")

# Step 3: Configure builder with ligand parameters
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
    ligand_params=ligand_results,  # Pass parametrized ligand files
)

print("\nStep 3: Builder configured with ligand parameters")
print(f"  Ligands in config: {list(builder.config['ligand_params'].keys())}")

# Step 4: Prepare system
success, message, job_dir = builder.prepare_system(
    pdb_file=pdb_file,
    working_dir=working_dir,
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1.0//1.0",
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
