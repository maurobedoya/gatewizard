from gatewizard.core.builder import Builder

builder = Builder()

# Configure for asymmetric membrane
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    # Note: When using anionic lipids (POPS, POPG, etc.), the system may require
    # a higher salt concentration for neutralization. If you get an error like:
    # "The concentration of ions required to neutralize the system is higher than
    # the concentration specified", increase salt_concentration or use add_salt=False.
    # POPS is anionic (-1 charge), so 50% POPS in lower leaflet adds significant
    # negative charge requiring more cations for neutralization.
    salt_concentration=0.5,  # Increased from default 0.15 M due to anionic lipids
    dist_wat=20.0  # Thicker water layer
)

# Upper leaflet: 70% POPC + 30% cholesterol
# Lower leaflet: 50% POPE + 50% POPS (anionic)
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC", "CHL1"],
    lower_lipids=["POPE", "POPS"],
    lipid_ratios="7:3//5:5",  # Ratios normalized automatically
    output_folder_name="asymmetric_membrane"
)

if success:
    print(f"✓ System preparation started in background")
    print(f"  {message}")
    print(f"  Upper leaflet: 70% POPC, 30% CHL1")
    print(f"  Lower leaflet: 50% POPE, 50% POPS")
    print(f"  Job directory: {job_dir}")
    print(f"  Monitor: {job_dir / 'logs/preparation.log'}")
else:
    print(f"✗ Preparation failed: {message}")
