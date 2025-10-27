from pathlib import Path
from gatewizard.utils.namd_analysis import EnergyAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Multiple log files from equilibration_folder
log_files = [
    data_dir / "step1_equilibration.log",
    data_dir / "step2_equilibration.log",
    data_dir / "step3_equilibration.log"
]

# Initialize analyzer with custom time for each file (in nanoseconds)
analyzer = EnergyAnalyzer(
    [str(f) for f in log_files],
    file_times={
        "step1_equilibration.log": 0.1,  # 100 ps
        "step2_equilibration.log": 0.1,  # 100 ps
        "step3_equilibration.log": 0.1   # 100 ps
    }
)

# Generate 4-panel energy plot
analyzer.plot_energy(target_temperature=300,    # 300 K
                     target_pressure=1.01325,   # 1.01325 atm (1 bar)
                     time_units="ps",
                     save="energy_analysis_example_01.png",
                     dpi=300)

print(f"Energy analysis complete!")
print(f"Plot saved: energy_analysis_example_01.png")
