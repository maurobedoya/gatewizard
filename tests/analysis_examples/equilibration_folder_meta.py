"""Shared equilibration_folder trajectories, protocol times, and stride map.

Default protocol stage durations (ns) for Eq1–6; production set to 200 ns for
this dataset (base protocol ships 100 ns). DCDs were exported with stride
already applied, so analysis uses every stored frame (stride=1).
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List

DATA_DIR = Path(__file__).resolve().parent / "equilibration_folder"

TOPOLOGY = DATA_DIR / "system.pdb"

TRAJECTORY_NAMES = (
    "step1_equilibration.dcd",
    "step2_equilibration.dcd",
    "step3_equilibration.dcd",
    "step4_equilibration.dcd",
    "step5_equilibration.dcd",
    "step6_equilibration.dcd",
    "step7_production.dcd",
)

FILE_TIMES: Dict[str, float] = {
    "step1_equilibration.dcd": 0.125,
    "step2_equilibration.dcd": 0.5,
    "step3_equilibration.dcd": 0.25,
    "step4_equilibration.dcd": 0.5,
    "step5_equilibration.dcd": 1.0,
    "step6_equilibration.dcd": 47.625,
    "step7_production.dcd": 200.0,
}

FILE_STRIDES: Dict[str, int] = {name: 1 for name in TRAJECTORY_NAMES}

TOTAL_TIME_NS = sum(FILE_TIMES.values())


def trajectory_paths() -> List[Path]:
    return [DATA_DIR / name for name in TRAJECTORY_NAMES]


def data_available() -> bool:
    if not TOPOLOGY.is_file():
        return False
    return all(p.is_file() for p in trajectory_paths())
