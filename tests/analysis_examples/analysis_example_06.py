from pathlib import Path
from gatewizard.utils.namd_analysis import TrajectoryAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Path to topology in equilibration_folder
topology_file = data_dir / "system.pdb"

# Multiple trajectory files from equilibration_folder
trajectory_files = [
    data_dir / "step1_equilibration.dcd",
    data_dir / "step2_equilibration.dcd",
    data_dir / "step3_equilibration.dcd"
]

# Initialize analyzer with custom time for each file (in nanoseconds)
analyzer = TrajectoryAnalyzer(
    topology_file,
    trajectory_files,
    file_times={
        "step1_equilibration.dcd": 0.1,  # 100 ps
        "step2_equilibration.dcd": 0.1,  # 100 ps
        "step3_equilibration.dcd": 0.1   # 100 ps
    }
)

# Calculate RMSF for alpha carbon atoms
data = analyzer.calculate_rmsf("protein and name CA")
print(f"Residues analyzed: {len(data['resids'])}")
print(f"RMSF range: {data['rmsf'].min():.2f} - {data['rmsf'].max():.2f} Å")
print(f"Mean RMSF: {data['rmsf'].mean():.2f} Å")
