from gatewizard.core.builder import Builder

builder = Builder()

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
    print("✓ All inputs are valid, proceed with preparation")
    # Now call prepare_system()
else:
    print(f"✗ Validation failed: {error_msg}")
    # Fix issues before proceeding
