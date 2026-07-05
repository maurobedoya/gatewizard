import tempfile
from pathlib import Path

from gatewizard.tools.packmol_hydration import (
    check_packmol_available,
    estimate_cavity_volume,
    hydrate_cavity,
)

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

if check_packmol_available()["available"]:
    vol = estimate_cavity_volume(PDB_FILE, BOX_MIN, BOX_MAX)
    n_waters = max(1, min(vol.suggested_waters, 20))
    with tempfile.TemporaryDirectory() as tmp:
        result = hydrate_cavity(
            pdb_file=PDB_FILE,
            working_dir=tmp,
            output_folder_name="hydration_6RV3_AB",
            box_min=BOX_MIN,
            box_max=BOX_MAX,
            n_waters=n_waters,
        )
        print(f"Success: {result.success}")
        print(f"Output: {result.output_pdb}")
        print(f"Message: {result.message}")
else:
    print("PACKMOL not installed; skipping hydration run")
