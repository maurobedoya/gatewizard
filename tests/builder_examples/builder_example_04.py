from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()

# Residue names from packmol-memgen --available_lipids_all
lipids = ff_manager.get_available_lipids()
catalog = {row["name"]: row for row in ff_manager.get_lipid_catalog()}

print(f"Total available lipids: {len(lipids)}\n")
for name in ("POPC", "POPE", "CHL1", "1LAC"):
    row = catalog[name]
    print(f"  {name}  charge {row['charge']}  {row['full_name']}")

print(f"\nTotal: {len(lipids)} lipid models available")
