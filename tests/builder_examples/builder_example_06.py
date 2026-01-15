from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get available lipid force fields
lipid_ffs = ff_manager.get_lipid_force_fields()

print("\nAvailable Lipid Force Fields:")
for lipid_ff in lipid_ffs:
    print(f"  - {lipid_ff}")

print(f"\nTotal: {len(lipid_ffs)} lipid force fields available")
