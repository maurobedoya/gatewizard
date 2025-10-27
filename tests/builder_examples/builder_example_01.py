from gatewizard.core.builder import Builder

builder = Builder()
print(f"Default water model: {builder.config['water_model']}")
print(f"Default protein force field: {builder.config['protein_ff']}")
print(f"Default lipid parameters: {builder.config['lipid_ff']}")
