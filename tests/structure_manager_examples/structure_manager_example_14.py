import os
import tempfile
import numpy as np
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

# Structure extended along the X-axis (principal axis ≈ X)
# We will align it so the principal axis points along Z
pdb_content = """\
ATOM      1  N   ALA A   1       0.000   0.200   0.100  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.500   0.100  -0.050  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000  -0.100   0.200  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   0.900   0.100  1.00  0.00           O
ATOM      5  N   ALA A   2       4.500   0.050  -0.100  1.00  0.00           N
ATOM      6  CA  ALA A   2       6.000  -0.200   0.050  1.00  0.00           C
ATOM      7  C   ALA A   2       7.500   0.100   0.150  1.00  0.00           C
ATOM      8  O   ALA A   2       8.000   1.100  -0.050  1.00  0.00           O
ATOM      9  N   ALA A   3       9.000  -0.050   0.000  1.00  0.00           N
ATOM     10  CA  ALA A   3      10.500   0.150   0.100  1.00  0.00           C
ATOM     11  C   ALA A   3      12.000  -0.100  -0.050  1.00  0.00           C
ATOM     12  O   ALA A   3      12.500   0.800   0.200  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    atoms = viewer.structure.atoms

    # Before alignment: measure span along each axis
    coords_before = np.array([a.coord for a in atoms])
    spans_before = coords_before.max(axis=0) - coords_before.min(axis=0)
    print("Before alignment (axis spans):")
    print(f"  X: {spans_before[0]:.2f} Å")
    print(f"  Y: {spans_before[1]:.2f} Å")
    print(f"  Z: {spans_before[2]:.2f} Å")
    print(f"  Principal axis: X (largest span)")

    # --- Align all atoms to Z-axis ---
    all_idx = list(range(len(atoms)))
    n = viewer.align_to_axis(all_idx, target_axis="z")
    coords_after = np.array([a.coord for a in atoms])
    spans_after = coords_after.max(axis=0) - coords_after.min(axis=0)
    print(f"\nAligned {n} atoms to Z-axis (axis spans):")
    print(f"  X: {spans_after[0]:.2f} Å")
    print(f"  Y: {spans_after[1]:.2f} Å")
    print(f"  Z: {spans_after[2]:.2f} Å")
    print(f"  Principal axis: Z (largest span)")

    # --- Align using only backbone CA atoms, transform all ---
    viewer.load_structure(tmp_path)
    ca_idx = [i for i, a in enumerate(viewer.structure.atoms) if a.name == "CA"]
    n = viewer.align_to_axis(ca_idx, target_axis="y")
    coords_ca = np.array([a.coord for a in viewer.structure.atoms])
    spans_ca = coords_ca.max(axis=0) - coords_ca.min(axis=0)
    print(f"\nAligned CA atoms to Y-axis, transformed all {n} atoms:")
    print(f"  X: {spans_ca[0]:.2f} Å")
    print(f"  Y: {spans_ca[1]:.2f} Å")
    print(f"  Z: {spans_ca[2]:.2f} Å")
finally:
    os.unlink(tmp_path)
