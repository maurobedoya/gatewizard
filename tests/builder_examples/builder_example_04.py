from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get available protein force fields
protein_ffs = ff_manager.get_protein_force_fields()

print("\nAvailable Protein Force Fields:")
for protein_ff in protein_ffs:
    print(f"  - {protein_ff}")

print(f"\nTotal: {len(protein_ffs)} protein force fields available")
