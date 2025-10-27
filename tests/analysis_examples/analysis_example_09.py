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
    data_dir / "step3_equilibration.dcd",
    data_dir / "step4_equilibration.dcd",
    data_dir / "step5_equilibration.dcd",  
    data_dir / "step6_equilibration.dcd",
    data_dir / "step7_production.dcd",
]

# Initialize analyzer with custom time for each file (in nanoseconds)
analyzer = TrajectoryAnalyzer(
    topology_file,
    trajectory_files,
    file_times={
        "step1_equilibration.dcd": 0.1,  # 100 ps
        "step2_equilibration.dcd": 0.1,  # 100 ps
        "step3_equilibration.dcd": 0.1,   # 100 ps
        "step4_equilibration.dcd": 0.1,  # 100 ps
        "step5_equilibration.dcd": 0.1,  # 100 ps
        "step6_equilibration.dcd": 0.1,  # 100 ps
        "step7_production.dcd": 0.1      # 100 ps
    }
)
# Calculate and plot radius of gyration of a selection
analyzer.plot_radius_of_gyration(
    selection="protein",
    bg_color="white",
    fig_bg_color="white",
    text_color="black",
    line_width=3,
    save="trajectory_analysis_example_09_rdgyr.png",
    dpi=300,
)
