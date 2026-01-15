from gatewizard.core.builder import Builder

builder = Builder()

# Configure system
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    salt_concentration=0.15,
    cation="K+",
    anion="Cl-"
)

# Prepare system with 100% POPC (symmetric)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",  # 100% POPC both leaflets
    output_folder_name="popc_membrane"
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor progress: {job_dir / 'preparation.log'}")
    print(f"  When complete, files will be at:")
    print(f"    - Topology: {job_dir / 'system.prmtop'}")
    print(f"    - Coordinates: {job_dir / 'system.inpcrd'}")
else:
    print(f"✗ Preparation failed: {message}")