from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# Point to folder with system files
work_dir = Path(__file__).parent / "popc_membrane"

# Explicitly define system files (if auto-detection doesn't work)
system_files = {
    'prmtop': str(work_dir / 'system.prmtop'),
    'inpcrd': str(work_dir / 'system.inpcrd'),
    'pdb': str(work_dir / 'system.pdb'),
    'bilayer_pdb': str(work_dir / 'bilayer_protein_protonated_prepared_lipid.pdb')
}

# Define equilibration stages
stages = [
    {
        'name': 'Equilibration 1',
        'time_ns': 0.125,
        'steps': 125000,
        'ensemble': 'NVT',
        'temperature': 310.15,
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
    }
]

# Setup with explicit file paths
# scheme_type is auto-detected from the 'ensemble' field in stages
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    output_name="equilibration_example_02"
)

print(f"Setup complete: {result['namd_dir']}")
# Run with: cd {result['namd_dir']} && ./run_equilibration.sh
