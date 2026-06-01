import os
import tempfile
import numpy as np
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

# Structure offset from origin so centering is visible
pdb_content = """\
ATOM      1  N   ALA A   1      10.000  20.000  30.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      11.500  20.000  30.000  1.00  0.00           C
ATOM      3  C   ALA A   1      13.000  20.000  30.000  1.00  0.00           C
ATOM      4  O   ALA A   1      13.500  21.000  30.000  1.00  0.00           O
ATOM      5  N   ALA A   2      14.500  20.000  30.000  1.00  0.00           N
ATOM      6  CA  ALA A   2      16.000  20.000  30.000  1.00  0.00           C
ATOM      7  C   ALA A   2      17.500  20.000  30.000  1.00  0.00           C
ATOM      8  O   ALA A   2      18.000  21.000  30.000  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    atoms = viewer.structure.atoms

    # --- Translate all atoms by (5, -10, 0) Å ---
    before = np.array([a.coord.copy() for a in atoms])
    n = viewer.translate_atoms([5.0, -10.0, 0.0])
    after = np.array([a.coord for a in atoms])
    print(f"Translated {n} atoms by (5, -10, 0) Å:")
    print(
        f"  Atom 1 before: ({before[0][0]:.1f}, {before[0][1]:.1f}, {before[0][2]:.1f})"
    )
    print(f"  Atom 1 after:  ({after[0][0]:.1f}, {after[0][1]:.1f}, {after[0][2]:.1f})")

    # --- Translate only residue 2 ---
    viewer.load_structure(tmp_path)
    res2 = [i for i, a in enumerate(viewer.structure.atoms) if a.res_id == 2]
    n = viewer.translate_atoms([0.0, 0.0, 5.0], indices=res2)
    print(f"\nTranslated {n} atoms (residue 2) by (0, 0, 5) Å")

    # --- Center structure at origin ---
    viewer.load_structure(tmp_path)
    before_center = np.array([a.coord for a in viewer.structure.atoms]).mean(axis=0)
    shift = viewer.center_atoms()
    after_center = np.array([a.coord for a in viewer.structure.atoms]).mean(axis=0)
    print(f"\nCentered structure:")
    print(
        f"  Centroid before: ({before_center[0]:.1f}, {before_center[1]:.1f}, {before_center[2]:.1f})"
    )
    print(f"  Shift applied:   ({shift[0]:.1f}, {shift[1]:.1f}, {shift[2]:.1f})")
    print(
        f"  Centroid after:  ({after_center[0]:.4f}, {after_center[1]:.4f}, {after_center[2]:.4f})"
    )

    # --- Center using a subset as reference, shift applied to all ---
    viewer.load_structure(tmp_path)
    res1 = [i for i, a in enumerate(viewer.structure.atoms) if a.res_id == 1]
    shift = viewer.center_atoms(indices=res1)
    res1_center = np.array([viewer.structure.atoms[i].coord for i in res1]).mean(axis=0)
    print(f"\nCentered on residue 1:")
    print(
        f"  Residue 1 centroid: ({res1_center[0]:.4f}, {res1_center[1]:.4f}, {res1_center[2]:.4f})"
    )
finally:
    os.unlink(tmp_path)
