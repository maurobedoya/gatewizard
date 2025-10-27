from pathlib import Path
from gatewizard.utils.namd_analysis import EnergyAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent

# Path to log file in equilibration_folder
log_file = script_dir / "equilibration_folder" / "step1_equilibration.log"

# Initialize energy analyzer
analyzer = EnergyAnalyzer(log_file)
props = analyzer.get_available_properties()
print(f"Can plot: {', '.join(props)}")
# Output: Can plot: Total Energy, Potential Energy, Kinetic Energy, Electrostatic Energy, Van der Waals Energy, Bond Energy, Angle Energy, Dihedral Energy, Improper Energy, Temperature, Pressure, Volume
