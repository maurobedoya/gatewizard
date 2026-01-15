from gatewizard.core.builder import Builder

builder = Builder()

# Plasma membrane-like composition
# Upper: PC-rich with cholesterol
# Lower: PE/PS-rich (cytoplasmic side)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC", "POPE", "CHL1"],
    lower_lipids=["POPE", "POPS", "CHL1"],
    lipid_ratios="5:2:3//4:4:2",  # Upper: 50% POPC, 20% POPE, 30% CHL1
                                   # Lower: 40% POPE, 40% POPS, 20% CHL1
    output_folder_name="plasma_membrane",
    salt_concentration=0.5,
    cation="Na+",  # Use sodium instead of potassium
    dist_wat=25.0  # Extra water for large protein
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")