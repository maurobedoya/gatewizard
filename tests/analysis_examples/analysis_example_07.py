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

# Plot RMSF across all trajectories
analyzer.plot_rmsf(
    selection="protein and name C",
    xaxis_type= "residue_type_number",
    residue_name_format="triple",
    label_frequency="all",
    highlight_threshold=0.6,
    highlight_color="red",
    bg_color="white",
    fig_bg_color="white",
    text_color="black",
    show_grid=False,
    line_color="#1f77b4",
    line_width=2,
    #title=" ",
    save="trajectory_analysis_example_07.png",
    dpi=300,
    # other settings...
)