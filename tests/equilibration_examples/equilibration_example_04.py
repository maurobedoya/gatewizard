from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to system folder
work_dir = Path(__file__).parent / "popc_membrane"

stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 303.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 2',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 303.15,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 5.0,
            'protein_sidechain': 2.5,
            'lipid_head': 2.5,
            'lipid_tail': 2.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 3',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'timestep': 1.0,
        'constraints': {
            'protein_backbone': 2.5,
            'protein_sidechain': 1.0,
            'lipid_head': 1.0,
            'lipid_tail': 1.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Auto-detect files and setup
# scheme_type auto-detected from stages (NVT from first stage)
# Stage 3 uses NPT ensemble - warning will be logged
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    stage_params_list=stages,
    output_name="equilibration_example_04",
    namd_executable="namd3"
)

print(f"\n✓ Setup complete!")
print(f"  Config files: {len(result['config_files'])}")
print(f"  Run script: {result['run_script'].name}")
