import os
import tempfile
import numpy as np
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

# A small chain along the X-axis so rotations are easy to verify
pdb_content = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.500   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   1.000   0.000  1.00  0.00           O
ATOM      5  N   ALA A   2       4.500   0.000   0.000  1.00  0.00           N
ATOM      6  CA  ALA A   2       6.000   0.000   0.000  1.00  0.00           C
ATOM      7  C   ALA A   2       7.500   0.000   0.000  1.00  0.00           C
ATOM      8  O   ALA A   2       8.000   1.000   0.000  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    atoms = viewer.structure.atoms
    print(f"Loaded {len(atoms)} atoms")

    # --- Rotate all atoms 90° around Z (X → Y) ---
    before = np.array([a.coord.copy() for a in atoms])
    n = viewer.rotate_atoms(90, "z")
    after = np.array([a.coord for a in atoms])
    print(f"\nRotated {n} atoms 90° around Z:")
    print(
        f"  Atom 1 before: ({before[0][0]:.1f}, {before[0][1]:.1f}, {before[0][2]:.1f})"
    )
    print(f"  Atom 1 after:  ({after[0][0]:.1f}, {after[0][1]:.1f}, {after[0][2]:.1f})")

    # --- Rotate only first residue 45° around X ---
    viewer.load_structure(tmp_path)  # reload
    res1 = [i for i, a in enumerate(viewer.structure.atoms) if a.res_id == 1]
    n = viewer.rotate_atoms(45, "x", indices=res1)
    print(f"\nRotated {n} atoms (residue 1) 45° around X")

    # --- Rotate around origin instead of selection centroid ---
    viewer.load_structure(tmp_path)
    n = viewer.rotate_atoms(180, "y", center="origin")
    after180 = np.array([a.coord for a in viewer.structure.atoms])
    print(f"\nRotated {n} atoms 180° around Y (origin):")
    print(
        f"  Atom 1: ({after180[0][0]:.1f}, {after180[0][1]:.1f}, {after180[0][2]:.1f})"
    )
finally:
    os.unlink(tmp_path)
