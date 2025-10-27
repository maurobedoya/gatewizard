from gatewizard.core.builder import Builder

builder = Builder()

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