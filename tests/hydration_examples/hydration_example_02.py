from pathlib import Path

from gatewizard.tools.packmol_hydration import detect_hydrogen_status

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
status = detect_hydrogen_status(pdb_file=PDB_FILE)
print(f"Hydrogen status for 6RV3_AB: {status}")
