from pathlib import Path
from gatewizard.utils.namd_analysis import EnergyAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Full default-protocol logs (Eq1–6 + production)
log_files = [
    data_dir / "step1_equilibration.log",
    data_dir / "step2_equilibration.log",
    data_dir / "step3_equilibration.log",
    data_dir / "step4_equilibration.log",
    data_dir / "step5_equilibration.log",
    data_dir / "step6_equilibration.log",
    data_dir / "step7_production.log",
]

# Default protocol stage times (ns); production set to 200 ns for this dataset.
file_times = {
    "step1_equilibration.log": 0.125,
    "step2_equilibration.log": 0.5,
    "step3_equilibration.log": 0.25,
    "step4_equilibration.log": 0.5,
    "step5_equilibration.log": 1.0,
    "step6_equilibration.log": 47.625,
    "step7_production.log": 200.0,
}

analyzer = EnergyAnalyzer(
    [str(f) for f in log_files],
    file_times=file_times,
)

# Generate 4-panel energy plot
analyzer.plot_energy(
    target_temperature=300,  # 300 K
    target_pressure=1.01325,  # 1.01325 atm (1 bar)
    time_units="ns",
    save="energy_analysis_example_01.png",
    dpi=300,
)

print(f"Energy analysis complete!")
print(f"Plot saved: energy_analysis_example_01.png")
print(f"Total simulation time: {sum(file_times.values()):.3f} ns")
