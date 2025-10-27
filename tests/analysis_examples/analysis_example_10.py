from pathlib import Path
from gatewizard.utils.namd_analysis import TrajectoryAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# ============================================================================
# Dark Theme RMSD Analysis - Full Customization
# ============================================================================

topology_file = data_dir / "system.pdb"
# Multiple trajectory files from equilibration_folder
trajectory_files = [
    data_dir / "step1_equilibration.dcd",
    data_dir / "step2_equilibration.dcd",
    data_dir / "step3_equilibration.dcd",
    data_dir / "step4_equilibration.dcd",
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
    }
)

# RMSD plot with full dark theme customization
analyzer.plot_rmsd(
    selection="protein and backbone",
    reference_frame=0,
    align=True,
    distance_units="Å",
    time_units="ps",
    line_color="#00d9ff",        # Bright cyan
    line_width=2.0,
    line_style="-",
    bg_color="#1a1a2e",          # Dark blue-black
    fig_bg_color="#16213e",      # Darker border
    text_color="#eee",           # Light gray
    show_grid=True,
    xlim=None,
    ylim=None,
    title="Dark Theme RMSD - Full Customization",
    xlabel="Time (ps)",
    ylabel="RMSD (Å)",
    highlight_threshold=None,
    highlight_color="orange",
    highlight_alpha=0.2,
    show_convergence=True,
    convergence_color="#ff006e",  # Magenta
    convergence_style="--",
    convergence_width=2.0,
    hlines=None,
    hline_colors=None,
    hline_styles=None,
    hline_widths=None,
    vlines=None,
    vline_colors=None,
    vline_styles=None,
    vline_widths=None,
    save="dark_theme_rmsd_example_10.png",
    show=False,
    figsize=(12, 7),
    dpi=300
)

print(f"Dark theme RMSD plot saved: dark_theme_rmsd_example_10.png")

# ============================================================================
# Dark Theme RMSD - With Threshold and Reference Lines
# ============================================================================

analyzer.plot_rmsd(
    selection="protein and backbone",
    align=True,
    distance_units="Å",
    time_units="ps",
    line_color="#7fff00",        # Chartreuse
    line_width=1.8,
    line_style="-",
    bg_color="#0d1117",          # GitHub dark
    fig_bg_color="#010409",
    text_color="#c9d1d9",
    show_grid=True,
    title="Dark Theme RMSD - With Threshold and Reference Lines",
    xlabel="Time (ps)",
    ylabel="RMSD (Å)",
    highlight_threshold=0.25,     # Highlight regions > 0.25 Å
    highlight_color="#ff6b6b",   # Red highlight
    highlight_alpha=0.5,
    show_convergence=True,
    convergence_color="#ffd700",  # Gold
    convergence_style="-.",
    convergence_width=1.5,
    hlines=[0.1, 0.2, 0.25],      # Horizontal reference lines
    hline_colors=["#4ecdc4", "#95e1d3", "#ff6b6b"],  # Teal to red gradient
    hline_styles=[":", ":", ":"],
    hline_widths=[1.0, 1.0, 1.0],
    save="dark_theme_rmsd_with_lines_example_10.png",
    figsize=(12, 7),
    dpi=300
)

print(f"Dark theme RMSD with lines saved: dark_theme_rmsd_with_lines_example_10.png")

# ============================================================================
# Dark Theme RMSD - Minimal Style
# ============================================================================

analyzer.plot_rmsd(
    selection="protein and backbone",
    align=True,
    distance_units="Å",
    time_units="ps",
    line_color="#ffffff",        # Pure white line
    line_width=1.5,
    line_style="-",
    bg_color="#000000",          # Pure black
    fig_bg_color="#000000",
    text_color="#ffffff",
    show_grid=False,             # No grid for minimal look
    title="",                    # No title
    show_convergence=False,      # No convergence line
    save="dark_theme_rmsd_minimal_example_10.png",
    figsize=(10, 6),
    dpi=300
)

print(f"Dark theme RMSD minimal saved: dark_theme_rmsd_minimal_example_10.png")

print(f"\nAll dark theme RMSD figures created with:")
print(f"  - Custom dark backgrounds")
print(f"  - Full customization options")
print(f"  - High resolution (300 DPI)")
