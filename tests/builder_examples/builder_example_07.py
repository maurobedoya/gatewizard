from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()
valid, message = ff_manager.validate_combination("tip3p", "ff14SB", "lipid21")

if valid:
    print("✓ Force field combination is compatible")
else:
    print(f"✗ Incompatible: {message}")
