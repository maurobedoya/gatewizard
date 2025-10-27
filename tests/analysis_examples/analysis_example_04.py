from pathlib import Path
from gatewizard.utils.namd_analysis import EnergyAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent

# Path to log file in equilibration_folder
log_file = script_dir / "equilibration_folder" / "step2_equilibration.log"

# Initialize energy analyzer
analyzer = EnergyAnalyzer(log_file)

stats = analyzer.get_statistics()

# Temperature
temp = stats['temp']
print(f"Temperature: {temp['mean']:.1f} ± {temp['std']:.1f} K")
print(f"  Range: {temp['min']:.1f} - {temp['max']:.1f} K")

# Energy
energy = stats['bond']
print(f"Bond Energy: {energy['mean']:.0f} kcal/mol")
print(f"  Final: {energy['final']:.0f} kcal/mol")
print(f"  Convergence: {abs(energy['final'] - energy['mean']):.0f} kcal/mol")

# Output: 
# Temperature: 303.6 ± 2.1 K
#   Range: 299.4 - 306.3 K
# Bond Energy: 2080 kcal/mol
#   Final: 2097 kcal/mol
#   Convergence: 17 kcal/mol