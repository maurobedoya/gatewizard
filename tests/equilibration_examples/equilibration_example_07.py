from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# System folder
work_dir = Path(__file__).parent / "popc_membrane"

system_files = {
    'prmtop': str(work_dir / 'system.prmtop'),
    'inpcrd': str(work_dir / 'system.inpcrd'),
    'pdb': str(work_dir / 'system.pdb'),
    'bilayer_pdb': str(work_dir / 'bilayer_protein_protonated_prepared_lipid.pdb')
}

# Use explicit templates to skip intermediate stages
stages = [
    {
        'name': 'Strong Restraints Phase',
        'ensemble': 'NPT',
        'custom_template': 'step6.1_equilibration.inp',  # Use stage 1 template
        'time_ns': 0.25,
        'timestep': 1.0,
        'temperature': 310.15,
        'pressure': 1.0,
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
        'name': 'Medium Restraints Phase',
        'ensemble': 'NPT',
        'custom_template': 'step6.3_equilibration.inp',  # Skip to stage 3 template
        'time_ns': 0.5,
        'timestep': 1.0,
        'temperature': 310.15,
        'pressure': 1.0,
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
        'name': 'Light Restraints Phase',
        'ensemble': 'NPAT',
        'custom_template': 'step6.5_equilibration.inp',  # Use stage 5 template
        'time_ns': 1.0,
        'timestep': 2.0,
        'temperature': 310.15,
        'pressure': 1.0,
        'surface_tension': 0.0,
        'constraints': {
            'protein_backbone': 1.0,
            'protein_sidechain': 0.5,
            'lipid_head': 0.5,
            'lipid_tail': 0.0,
            'water': 0.0,
            'ions': 0.0,
            'other': 0.0
        }
    }
]

# Setup with custom template selection
manager = NAMDEquilibrationManager(work_dir)
result = manager.setup_namd_equilibration(
    system_files=system_files,
    stage_params_list=stages,
    output_name="equilibration_example_07",
    namd_executable="namd3"
)

print(f"\n✓ Setup complete with custom templates!")
print(f"  Stage 1: Using template step6.1")
print(f"  Stage 2: Using template step6.3 (skipped step6.2)")
print(f"  Stage 3: Using template step6.5 (skipped step6.4)")