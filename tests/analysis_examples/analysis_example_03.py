from pathlib import Path
from gatewizard.utils.namd_analysis import EnergyAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Multiple log files from equilibration_folder
log_files = [
    data_dir / "step1_equilibration.log",
    data_dir / "step2_equilibration.log",
    data_dir / "step3_equilibration.log",
    data_dir / "step4_equilibration.log",
    data_dir / "step5_equilibration.log",
    data_dir / "step6_equilibration.log",
    data_dir / "step7_production.log"
]

# Initialize analyzer with custom time for each file (in nanoseconds)
analyzer = EnergyAnalyzer(
    [str(f) for f in log_files],
    file_times={
        "step1_equilibration.log": 0.1,  # 100 ps
        "step2_equilibration.log": 0.1,  # 100 ps
        "step3_equilibration.log": 0.1,   # 100 ps
        "step4_equilibration.log": 0.1,  # 100 ps
        "step5_equilibration.log": 0.1,  # 100 ps
        "step6_equilibration.log": 0.1,  # 100 ps
        "step7_production.log": 0.1,      # 100 ps
    }
)

# Plot specific energy properties
analyzer.plot_properties(
    properties=["bond energy", "angle energy", "dihedral energy"],
    energy_units="kcal/mol",
    time_units="ps",
    bg_color = "#ffffff",
    fig_bg_color = "#FFFFFF",
    save="energy_analysis_example_03.png",
    dpi=300,
)

print(f"Energy properties plot complete!")
print(f"Plot saved: energy_analysis_example_03.png")
