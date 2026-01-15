from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Get all available lipids
lipids = ff_manager.get_available_lipids()

print(f"Total available lipids: {len(lipids)}\n")

for lipid in lipids:
     print(f"  - {lipid}")

print(f"\nTotal: {len(lipids)} lipid models available")
