"""
analysis_example_14.py
======================
Area per lipid analysis using lipyphilic (BilayerTrajectoryAnalyzer).

Prerequisites
-------------
- ``equilibration_folder/system.pdb`` and equilibration DCD trajectories
- lipyphilic (core gatewizard dependency; installed via ``pip install -e .``)
"""

from pathlib import Path

from gatewizard.utils.namd_analysis import BilayerTrajectoryAnalyzer, run_bilayer_analysis

script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

topology_file = data_dir / "system.pdb"
trajectory_files = [
    data_dir / "step1_equilibration.dcd",
    data_dir / "step2_equilibration.dcd",
    data_dir / "step3_equilibration.dcd",
]

# POPC headgroup phosphorus atoms in the AMBER lipid parametrization
LIPID_SEL = "resname PC and name P31"

file_times = {
    "step1_equilibration.dcd": 0.1,  # 0.1 ns
    "step2_equilibration.dcd": 0.1,  # 0.1 ns
    "step3_equilibration.dcd": 0.1,  # 0.1 ns
}

# ------------------------------------------------------------------
# 1. Class API — calculate area per lipid
# ------------------------------------------------------------------
analyzer = BilayerTrajectoryAnalyzer(
    topology_file,
    trajectory_files,
    file_times=file_times,
)
data = analyzer.calculate_area_per_lipid(lipid_sel=LIPID_SEL)

mean_area = float(data["mean_area_per_lipid"].mean())
print(f"Mean area per lipid: {mean_area:.1f} Å²")
print(f"Lipids analysed: {len(data['resids'])}")
print(f"Frames analysed: {data['areas'].shape[1]}")

# ------------------------------------------------------------------
# 2. Plot mean area per lipid time series
# ------------------------------------------------------------------
analyzer.plot_area_per_lipid(
    lipid_sel=LIPID_SEL,
    series="mean",
    time_units="ns",
    save="area_per_lipid_example_14.png",
    show=False,
)
print("Plot saved: area_per_lipid_example_14.png")

# ------------------------------------------------------------------
# 3. JSON API — run_bilayer_analysis()
# ------------------------------------------------------------------
result = run_bilayer_analysis(
    topology_file,
    trajectory_files,
    analysis_type="area_per_lipid",
    lipid_sel=LIPID_SEL,
    file_times=file_times,
)
print(f"JSON API mean area: {result['stats']['mean']:.1f} Å²")
