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
        "step3_equilibration.dcd": 0.1,  # 100 ps
        "step4_equilibration.dcd": 0.1,  # 100 ps
        "step5_equilibration.dcd": 0.1,  # 100 ps
        "step6_equilibration.dcd": 0.1,  # 100 ps
        "step7_production.dcd": 0.1      # 100 ps
    }
)
# Calculate and plot distances between selections
analyzer.plot_distances(
    selections={
        "gate_distance": ("resid 1 and name C", "resid 28 and name C"),
        "domain_distance": ("resid 1-2 and name C", "resid 9-10 and name C")
    },
    bg_color="white",
    fig_bg_color="white",
    text_color="black",
    line_width=3,
    dpi=300,
    show_mean_lines=False,
    figsize=(5, 4),
    save="trajectory_analysis_example_08_distances.png",
)

print(f"Distances plot saved: trajectory_analysis_example_08_distances.png")
