"""
analysis_example_11.py
======================
Basic OpenMM log analysis — single log file.

Prerequisites
-------------
- ``step7_production.log`` output from a StateDataReporter run
- matplotlib / numpy (``conda install matplotlib numpy``)

"""

from pathlib import Path
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer

# Get the directory where this script is located
script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Single OpenMM log file (StateDataReporter output)
log_file = data_dir / "step7_production.log"

# Initialize the analyzer
analyzer = OpenMMLogAnalyzer(log_file)

# ------------------------------------------------------------------
# 1. Get statistics for all properties
# ------------------------------------------------------------------
stats = analyzer.get_statistics()

for key, s in stats.items():
    print(
        f"  {key:20s}  mean={s['mean']:12.3f}  std={s['std']:10.3f}"
        f"  min={s['min']:12.3f}  max={s['max']:12.3f}"
    )
# → All values in kJ/mol (energies) or native units (temperature in K,
#   volume in nm³, density in g/mL)

# ------------------------------------------------------------------
# 2. 2×2 summary plot: total energy, potential + kinetic, temperature,
#    volume / density
# ------------------------------------------------------------------
analyzer.plot_energy(
    save="energy_summary_example_11.png",
    show=False,
    target_temperature=303.15,
)

print("Energy summary saved: energy_summary_example_11.png")

# ------------------------------------------------------------------
# 3. Same plot converted to kcal/mol
# ------------------------------------------------------------------
analyzer.plot_energy(
    energy_units="kcal/mol",
    save="energy_summary_kcal_example_11.png",
    show=False,
)

print("Energy summary (kcal/mol) saved: energy_summary_kcal_example_11.png")
