"""
analysis_example_12.py
======================
Multi-stage OpenMM analysis — concatenate several log files.

Prerequisites
-------------
- Log files ``step1_equilibration.log`` … ``step7_production.log``
  (StateDataReporter output, one per equilibration stage)
- matplotlib / numpy

"""

from pathlib import Path
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer

script_dir = Path(__file__).parent
data_dir = script_dir / "equilibration_folder"

# Pass a list of log files — time axis is concatenated automatically
log_files = [
    data_dir / "step1_equilibration.log",
    data_dir / "step2_equilibration.log",
    data_dir / "step3_equilibration.log",
    data_dir / "step4_equilibration.log",
    data_dir / "step5_equilibration.log",
    data_dir / "step6_equilibration.log",
    data_dir / "step7_production.log",
]

# Provide real stage durations (ns) to override the "Time (ps)" column
# — useful when logs lack the time column or have restarted counters.
file_times = {
    "step1_equilibration.log": 0.125,
    "step2_equilibration.log": 0.125,
    "step3_equilibration.log": 0.125,
    "step4_equilibration.log": 0.25,
    "step5_equilibration.log": 0.25,
    "step6_equilibration.log": 0.5,
    "step7_production.log": 50.0,
}

analyzer = OpenMMLogAnalyzer(log_files, file_times=file_times)

# ------------------------------------------------------------------
# 1. Print key statistics
# ------------------------------------------------------------------
stats = analyzer.get_statistics()

print("=== Multi-stage OpenMM analysis ===")
for key in ("potential", "kinetic", "total", "temp", "volume", "density"):
    if key not in stats:
        continue
    s = stats[key]
    print(
        f"  {key:12s}  mean={s['mean']:12.3f}"
        f"  initial={s['initial']:12.3f}  final={s['final']:12.3f}"
    )

# ------------------------------------------------------------------
# 2. Energy summary over the full trajectory
# ------------------------------------------------------------------
analyzer.plot_energy(
    save="energy_multistage_example_12.png",
    show=False,
    title="Full equilibration + production (51.375 ns)",
    target_temperature=303.15,
)

print("Saved: energy_multistage_example_12.png")
