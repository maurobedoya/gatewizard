from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# System folder
work_dir = Path(__file__).parent / "popc_membrane"

stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'minimize_steps': 10000,
        'timestep': 1.0,
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
        'ensemble': 'NPT',
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
        'surface_tension': 0.0,
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
    },
    {
        'name': 'Equilibration 4',
        'time_ns': 0.5,
        'steps': 250000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 1.0,
            'protein_sidechain': 0.5,
            'lipid_head': 0.5,
            'lipid_tail': 0.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 5',
        'time_ns': 0.5,
        'steps': 250000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.5,
            'protein_sidechain': 0.1,
            'lipid_head': 0.1,
            'lipid_tail': 0.1,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Equilibration 6',
        'time_ns': 0.5,
        'steps': 250000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.1,
            'protein_sidechain': 0.0,
            'lipid_head': 0.0,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Production',
        'time_ns': 10.0,
        'steps': 5000000,
        'ensemble': 'NPT',
        'temperature': 303.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.0,
            'protein_sidechain': 0.0,
            'lipid_head': 0.0,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Auto-detect and setup
# scheme_type auto-detected from first stage's ensemble
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    stage_params_list=stages,
    output_name="equilibration_example_06",
    namd_executable="namd3"
)

print(f"\n✓ Complete! Generated {len(result['config_files'])} configuration files")
print(f"  Total equilibration: {sum(s['time_ns'] for s in stages[:-1]):.3f} ns")
print(f"  Production: {stages[-1]['time_ns']:.1f} ns")