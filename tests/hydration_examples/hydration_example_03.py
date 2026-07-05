from pathlib import Path

from gatewizard.tools.packmol_hydration import estimate_cavity_volume

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

result = estimate_cavity_volume(
    pdb_file=PDB_FILE,
    box_min=BOX_MIN,
    box_max=BOX_MAX,
)
print(f"Box volume: {result.box_volume_A3:.1f} Å³")
print(f"Free volume: {result.free_volume_A3:.1f} Å³")
print(f"Suggested waters: {result.suggested_waters}")
print(f"Exclusion mode: {result.exclusion_mode}")
print(f"Hydrogen status: {result.hydrogen_status}")
