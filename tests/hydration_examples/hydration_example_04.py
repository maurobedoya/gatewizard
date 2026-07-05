from pathlib import Path

from gatewizard.tools.packmol_hydration import build_hydrate_inp_text, estimate_cavity_volume

PDB_FILE = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

vol = estimate_cavity_volume(PDB_FILE, BOX_MIN, BOX_MAX)
n_waters = max(1, min(vol.suggested_waters, 50))

inp = build_hydrate_inp_text(
    protein_path=PDB_FILE,
    tip3p_path="TIP3P.pdb",
    output_pdb="6RV3_AB_hydrated.pdb",
    box_min=BOX_MIN,
    box_max=BOX_MAX,
    n_waters=n_waters,
    solute_radius=2.5,
)
print(inp)
