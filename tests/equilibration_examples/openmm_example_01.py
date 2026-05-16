from pathlib import Path
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

# Point to folder with system files
work_dir = Path(__file__).parent / "popc_membrane"

# Define a single NVT equilibration stage
stages = [
    {
        "name": "Equilibration 1",
        "time_ns": 0.125,
        "ensemble": "NVT",
        "temperature": 310.15,
        "timestep": 1.0,
        "minimize_steps": 5000,
        "constraints": {
            "protein_backbone": 10.0,
            "protein_sidechain": 5.0,
            "lipid_head": 2.5,
            "lipid_tail": 0.0,
        },
    }
]

# Setup with automatic file detection (no system_files needed!)
# scheme_type is auto-detected from the 'ensemble' field in stages
manager = OpenMMEquilibrationManager(work_dir)
result = manager.setup_openmm_equilibration(
    stage_params_list=stages,
    output_name="openmm_example_01",
)

print(f"Setup complete: {result['openmm_dir']}")
# Run with: cd {result['openmm_dir']} && bash run_equilibration.sh
