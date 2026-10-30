"""
analysis_example_15.py
======================
Membrane thickness analysis using lipyphilic (BilayerTrajectoryAnalyzer).

Prerequisites
-------------
- ``equilibration_folder/system.pdb`` and equilibration DCD trajectories
- lipyphilic (core gatewizard dependency; installed via ``pip install -e .``)
"""

from pathlib import Path

from gatewizard.utils.lipid_bilayer_analysis import BilayerTrajectoryAnalyzer, run_bilayer_analysis

script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

topology_file = data_dir / "system.pdb"
trajectory_files = [
    data_dir / "step1_equilibration.dcd",
    data_dir / "step2_equilibration.dcd",
    data_dir / "step3_equilibration.dcd",
    data_dir / "step4_equilibration.dcd",
    data_dir / "step5_equilibration.dcd",
    data_dir / "step6_equilibration.dcd",
    data_dir / "step7_production.dcd",
]

LIPID_SEL = "resname PC and name P31"

# Default protocol stage times (ns); production set to 200 ns for this dataset.
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

# ------------------------------------------------------------------
# 1. Class API — calculate membrane thickness
# ------------------------------------------------------------------
analyzer = BilayerTrajectoryAnalyzer(
    topology_file,
    trajectory_files,
    file_times=file_times,
    file_strides=file_strides,
)
data = analyzer.calculate_membrane_thickness(lipid_sel=LIPID_SEL)

mean_thickness = float(data["thickness"].mean())
print(f"Mean membrane thickness: {mean_thickness:.1f} Å")
print(f"Frames analysed: {len(data['thickness'])}")

# ------------------------------------------------------------------
# 2. Plot membrane thickness time series
# ------------------------------------------------------------------
analyzer.plot_membrane_thickness(
    lipid_sel=LIPID_SEL,
    time_units="ns",
    save="membrane_thickness_example_15.png",
    show=False,
)
print("Plot saved: membrane_thickness_example_15.png")

# ------------------------------------------------------------------
# 3. JSON API — run_bilayer_analysis()
# ------------------------------------------------------------------
result = run_bilayer_analysis(
    topology_file,
    trajectory_files,
    analysis_type="membrane_thickness",
    lipid_sel=LIPID_SEL,
    file_times=file_times,
    file_strides=file_strides,
)
print(f"JSON API mean thickness: {result['stats']['mean']:.1f} Å")
