"""
analysis_example_13.py
======================
OpenMM plot_properties — select specific columns, custom colours.

Prerequisites
-------------
- Log file ``step7_production.log`` from a StateDataReporter run
- matplotlib / numpy

"""

from pathlib import Path
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer

script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

log_file = data_dir / "step7_production.log"
analyzer = OpenMMLogAnalyzer(log_file)

# ------------------------------------------------------------------
# 1. Plot all available properties (auto-detected)
# ------------------------------------------------------------------
analyzer.plot_properties(
    save="all_properties_example_13.png",
    show=False,
)
print("All-properties plot saved: all_properties_example_13.png")

# ------------------------------------------------------------------
# 2. Plot only energies, converted to kcal/mol
# ------------------------------------------------------------------
analyzer.plot_properties(
    properties=["potential", "kinetic", "total"],
    energy_units="kcal/mol",
    save="energies_kcal_example_13.png",
    show=False,
    colors=["#61afef", "#98c379", "#e06c75"],
)
print("Energy plots (kcal/mol) saved: energies_kcal_example_13.png")

# ------------------------------------------------------------------
# 3. Temperature and density as separate figures
# ------------------------------------------------------------------
analyzer.plot_properties(
    properties=["temp", "density"],
    separate_plots=True,
    save="temp_density_example_13",
    show=False,
)
print(
    "Separate figures saved: temp_density_example_13_temp.png, temp_density_example_13_density.png"
)
