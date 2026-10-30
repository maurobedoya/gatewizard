from pathlib import Path
from gatewizard.utils.trajectory_analysis import TrajectoryAnalyzer

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

# Default protocol stage times (ns); production set to 200 ns for this dataset.
# DCDs were exported with stride already applied — analyse every stored frame.
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

# Plot RMSF across all trajectories
analyzer.plot_rmsf(
    selection="protein and name CA",
    xaxis_type="residue_type_number",
    residue_name_format="triple",
    label_frequency="every_50",
    highlight_threshold=0.6,
    highlight_color="red",
    bg_color="white",
    fig_bg_color="white",
    text_color="black",
    show_grid=False,
    line_color="#1f77b4",
    line_width=2,
    # title=" ",
    save="trajectory_analysis_example_07.png",
    dpi=300,
    # other settings...
)
