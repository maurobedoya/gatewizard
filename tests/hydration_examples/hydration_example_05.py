import tempfile
from pathlib import Path

from gatewizard.tools.packmol_hydration import estimate_cavity_volume, prepare_hydration_job

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

vol = estimate_cavity_volume(PDB_FILE, BOX_MIN, BOX_MAX)
n_waters = max(1, min(vol.suggested_waters, 50))

with tempfile.TemporaryDirectory() as tmp:
    job = prepare_hydration_job(
        pdb_file=PDB_FILE,
        job_dir=tmp,
        box_min=BOX_MIN,
        box_max=BOX_MAX,
        n_waters=n_waters,
    )
    print(f"Job dir: {job['job_dir']}")
    print(f"Input file: {job['packmol_inp_path']}")
    print(f"Output PDB name: {job['output_pdb_name']}")
    for name in ("packmol.inp", "TIP3P.pdb", Path(PDB_FILE).name):
        path = Path(tmp) / name
        print(f"  {name}: {'OK' if path.is_file() else 'MISSING'}")
