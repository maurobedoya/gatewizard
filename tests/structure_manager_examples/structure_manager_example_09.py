import os
import tempfile
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
HETATM    5  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

out_path = tmp_path + "_out.pdb"

try:
    viewer.load_structure(tmp_path)

    # Delete water atoms
    water_idx = viewer.select_by_criteria("Water")
    print(f"Deleting {len(water_idx)} water atoms")
    removed = viewer.delete_atoms(water_idx)
    print(f"Removed {removed} atoms")

    info = viewer.get_structure_info()
    print(f"After deletion: {info['n_atoms']} atoms")

    # Save modified structure
    saved = viewer.save_pdb(out_path)
    print(f"Saved to: {os.path.basename(saved)}")
finally:
    os.unlink(tmp_path)
    if os.path.exists(out_path):
        os.unlink(out_path)
