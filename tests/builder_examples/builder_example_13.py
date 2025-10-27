from gatewizard.core.builder import Builder

builder = Builder()

# High salt concentration for ionic strength studies
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC", "POPS"],
    lower_lipids=["POPC", "POPS"],
    lipid_ratios="8:2//8:2",  # 80% POPC, 20% POPS
    output_folder_name="high_salt",
    salt_concentration=2.0,  # 2000 mM (high salt)
    cation="Na+",
    anion="Cl-"
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  High salt: 2.0 M NaCl")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")