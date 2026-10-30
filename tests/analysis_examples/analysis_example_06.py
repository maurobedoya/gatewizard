from pathlib import Path
from gatewizard.utils.trajectory_analysis import TrajectoryAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Path to topology in equilibration_folder
topology_file = data_dir / "system.pdb"

# Full default-protocol trajectories (Eq1–6 + production)
trajectory_files = [
    data_dir / "step1_equilibration.dcd",
    data_dir / "step2_equilibration.dcd",
    data_dir / "step3_equilibration.dcd",
    data_dir / "step4_equilibration.dcd",
    data_dir / "step5_equilibration.dcd",
    data_dir / "step6_equilibration.dcd",
    data_dir / "step7_production.dcd",
]

# Default protocol stage times (ns); production set to 200 ns for this dataset.
file_times = {
    "step1_equilibration.dcd": 0.125,
    "step2_equilibration.dcd": 0.5,
    "step3_equilibration.dcd": 0.25,
    "step4_equilibration.dcd": 0.5,
    "step5_equilibration.dcd": 1.0,
    "step6_equilibration.dcd": 47.625,
    "step7_production.dcd": 200.0,
}
file_strides = {name: 1 for name in file_times}

analyzer = TrajectoryAnalyzer(
    topology_file,
    trajectory_files,
    file_times=file_times,
    file_strides=file_strides,
)

# Calculate RMSF for alpha carbon atoms
data = analyzer.calculate_rmsf("protein and name CA")
print(f"Residues analyzed: {len(data['resids'])}")
print(f"RMSF range: {data['rmsf'].min():.2f} - {data['rmsf'].max():.2f} Å")
print(f"Mean RMSF: {data['rmsf'].mean():.2f} Å")
