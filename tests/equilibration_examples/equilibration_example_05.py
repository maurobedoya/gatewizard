from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to system folder
work_dir = Path(__file__).parent / "popc_membrane"
system_files = {
    'prmtop': str(work_dir / 'system.prmtop'),
    'inpcrd': str(work_dir / 'system.inpcrd'),
    'pdb': str(work_dir / 'system.pdb'),
    'bilayer_pdb': str(work_dir / 'bilayer_protein_protonated_prepared_lipid.pdb')
}

custom_protocol = [
    {
        'name': 'Initial Equilibration',
        'time_ns': 0.25,
        'steps': 250000,
        'ensemble': 'NVT',
        'temperature': 310.15,
        'timestep': 1.0,
        'minimize_steps': 10000,
        'constraints': {
            'protein_backbone': 10.0,
            'protein_sidechain': 5.0,
            'lipid_head': 5.0,
            'lipid_tail': 5.0,
            'water': 0.0,
            'ions': 10.0,
            'other': 0.0
        }
    },
    {
        'name': 'Pressure Equilibration',
        'time_ns': 0.5,
        'steps': 500000,
        'ensemble': 'NPT',
        'temperature': 310.15,
        'pressure': 1.0,
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
        'name': 'Membrane Relaxation',
        'time_ns': 1.0,
        'steps': 500000,
        'ensemble': 'NPAT',
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 2.0,
            'protein_sidechain': 1.0,
            'lipid_head': 1.0,
            'lipid_tail': 0.5,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    },
    {
        'name': 'Production Preparation',
        'time_ns': 2.0,
        'steps': 1000000,
        'ensemble': 'NPAT',
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'timestep': 2.0,
        'constraints': {
            'protein_backbone': 0.5,
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
# scheme_type auto-detected from stages (NVT from first stage)
# Stages 2-4 use different ensembles - warnings will be logged
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    stage_params_list=custom_protocol,
    output_name="equilibration_example_05",
    namd_executable="namd3"
)

print(f"\n✓ Setup complete!")
print(f"  Total stages: {len(custom_protocol)}")
print(f"  Total time: {sum(s['time_ns'] for s in custom_protocol):.1f} ns")
print(f"\nTo run:")
print(f"  cd {result['namd_dir']}")
print(f"  ./run_equilibration.sh")