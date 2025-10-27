from gatewizard.core.builder import Builder

builder = Builder()

# Configure for specific system
builder.set_configuration(
    water_model="tip3p",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    salt_concentration=0.15,
    cation="Na+",
    anion="Cl-",
    dist_wat=20.0,  # Larger water layer
    preoriented=True
)
