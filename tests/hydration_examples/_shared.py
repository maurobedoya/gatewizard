"""Shared box bounds for hydration examples (6RV3_AB sub-region, Å)."""

from pathlib import Path

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")

# 20 Å cube around structure centroid (see tests/hydration_examples/README.md)
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)
