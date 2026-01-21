from gatewizard.tools.force_fields import ForceFieldManager

ff_manager = ForceFieldManager()
valid, message, is_warning = ff_manager.validate_combination("tip3p", "ff14SB", "lipid21")

if valid:
    if is_warning:
        print(f"[!] Warning: {message}")
    else:
        print("[OK] Force field combination is compatible")
else:
    print(f"[ERROR] Incompatible: {message}")
