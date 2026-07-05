from gatewizard.tools.packmol_hydration import check_packmol_available

info = check_packmol_available()
print(f"PACKMOL available: {info['available']}")
if info["resolved_path"]:
    print(f"Path: {info['resolved_path']}")
if info["version"]:
    print(f"Version: {info['version']}")
