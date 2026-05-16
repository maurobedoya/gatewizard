from pathlib import Path
from gatewizard.tools.equilibration import OpenMMEquilibrationManager

# Point to folder with system files
work_dir = Path(__file__).parent / "popc_membrane"

# Full 6-stage NPT membrane equilibration protocol (CHARMM-GUI style)
# Gradual relaxation of restraints following the standard protocol
stages = [
    {
        "name": "Equilibration 1 - NVT with strong restraints",
        "time_ns": 0.125,
        "ensemble": "NVT",
        "temperature": 303.15,
        "timestep": 1.0,
        "minimize_steps": 5000,
        "constraints": {
            "protein_backbone": 10.0,
            "protein_sidechain": 5.0,
            "lipid_head": 2.5,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 2 - NVT relaxing restraints",
        "time_ns": 0.125,
        "ensemble": "NVT",
        "temperature": 303.15,
        "timestep": 1.0,
        "constraints": {
            "protein_backbone": 5.0,
            "protein_sidechain": 2.5,
            "lipid_head": 1.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 3 - NPT with pressure coupling",
        "time_ns": 0.125,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 1.0,
        "constraints": {
            "protein_backbone": 2.5,
            "protein_sidechain": 1.0,
            "lipid_head": 0.5,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 4 - NPT further relaxing",
        "time_ns": 0.25,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 1.0,
            "protein_sidechain": 0.5,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 5 - NPT backbone only",
        "time_ns": 0.25,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 0.5,
            "protein_sidechain": 0.0,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Equilibration 6 - NPT light backbone restraints",
        "time_ns": 0.5,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 0.1,
            "protein_sidechain": 0.0,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
    {
        "name": "Production - NPT unrestrained",
        "time_ns": 50.0,
        "ensemble": "NPT",
        "temperature": 303.15,
        "timestep": 2.0,
        "constraints": {
            "protein_backbone": 0.0,
            "protein_sidechain": 0.0,
            "lipid_head": 0.0,
            "lipid_tail": 0.0,
        },
    },
]

# Setup with automatic file detection
# scheme_type is auto-detected from 'ensemble' field of first stage (NVT -> 01_NVT)
# Note: Mixed ensembles (NVT stages 1-2, NPT stages 3-7) are handled automatically.
#       The scheme_type controls which pressure coupling templates are used for
#       stages 3+ — pass scheme_type="NPT" explicitly if needed.
manager = OpenMMEquilibrationManager(work_dir)
result = manager.setup_openmm_equilibration(
    stage_params_list=stages,
    output_name="openmm_example_02",
    scheme_type="NPT",
)

print(f"Setup complete: {result['openmm_dir']}")
print(f"Config files: {len(result['config_files'])}")
print(f"Run script: {result['run_script'].name}")
# Run with: cd {result['openmm_dir']} && bash run_equilibration.sh
