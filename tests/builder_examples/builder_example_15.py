from gatewizard.core.builder import Builder

builder = Builder()

# Neutralize system charges only (no extra salt)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    output_folder_name="no_salt",
    salt_concentration=0.0,  # Only neutralize, no extra salt
    add_salt=True  # Still add ions for neutralization
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Neutralization only, no extra salt")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")