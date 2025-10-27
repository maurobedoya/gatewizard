from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

manager = NAMDEquilibrationManager(Path("popc_membrane"))

# Check which files will be used
system_files = manager.find_system_files()
if system_files:
    print("Detected files:")
    for key, path in system_files.items():
        print(f"  {key}: {Path(path).name}")
    
    # Now run setup with auto-detection
    #result = manager.setup_namd_equilibration(
    #    stage_params_list=stages
    #)
else:
    print("Required files not found - please check working directory")