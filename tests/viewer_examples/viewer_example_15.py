import os
import tempfile
import numpy as np
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

# A structure with two chains:
#   Chain A runs along the X-axis (the "channel pore")
#   Chain B has an atom offset in the Y direction (reference for secondary axis)
# This mimics aligning a channel pore to Z with a pore-lining residue on X.
pdb_content = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   0.200   0.100  1.00  0.00           C
ATOM      3  C   ALA A   1       4.000  -0.100   0.050  1.00  0.00           C
ATOM      4  O   ALA A   1       6.000   0.100  -0.100  1.00  0.00           O
ATOM      5  N   ALA A   2       8.000  -0.050   0.200  1.00  0.00           N
ATOM      6  CA  ALA A   2      10.000   0.150  -0.050  1.00  0.00           C
ATOM      7  C   ALA A   2      12.000   0.000   0.100  1.00  0.00           C
ATOM      8  O   ALA A   2      14.000  -0.200   0.000  1.00  0.00           O
ATOM      9  N   GLY B   1       5.000   4.000   0.500  1.00  0.00           N
ATOM     10  CA  GLY B   1       7.000   4.200   0.300  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    atoms = viewer.structure.atoms

    # Primary: chain A backbone → align to Z-axis
    # Secondary: chain B atoms → align to X-axis
    chainA = [i for i, a in enumerate(atoms) if a.chain_id == 'A']
    chainB = [i for i, a in enumerate(atoms) if a.chain_id == 'B']

    print("Before alignment:")
    coordsA = np.array([atoms[i].coord for i in chainA])
    coordsB = np.array([atoms[i].coord for i in chainB])
    spansA = coordsA.max(axis=0) - coordsA.min(axis=0)
    print(f"  Chain A spans: X={spansA[0]:.1f}, Y={spansA[1]:.1f}, Z={spansA[2]:.1f}")
    print(f"  Chain B centroid: ({coordsB.mean(0)[0]:.1f}, "
          f"{coordsB.mean(0)[1]:.1f}, {coordsB.mean(0)[2]:.1f})")

    # Align with primary + secondary axes
    n = viewer.align_to_axis(
        primary_indices=chainA,
        target_axis='z',
        secondary_indices=chainB,
        secondary_axis='x',
    )

    print(f"\nAligned {n} atoms (primary → Z, secondary → X):")
    coordsA2 = np.array([atoms[i].coord for i in chainA])
    coordsB2 = np.array([atoms[i].coord for i in chainB])
    spansA2 = coordsA2.max(axis=0) - coordsA2.min(axis=0)
    print(f"  Chain A spans: X={spansA2[0]:.2f}, Y={spansA2[1]:.2f}, Z={spansA2[2]:.2f}")
    print(f"  Chain B centroid: ({coordsB2.mean(0)[0]:.2f}, "
          f"{coordsB2.mean(0)[1]:.2f}, {coordsB2.mean(0)[2]:.2f})")
    print(f"  Chain A now mostly along Z (Z span >> X, Y).")
    print(f"  Chain B centroid now has largest offset along X.")

    # --- Align only chain A, keep chain B fixed ---
    viewer.load_structure(tmp_path)
    atoms = viewer.structure.atoms
    chainA = [i for i, a in enumerate(atoms) if a.chain_id == 'A']
    chainB_fixed = [i for i, a in enumerate(atoms) if a.chain_id == 'B']
    coordsB_before = np.array([atoms[i].coord.copy() for i in chainB_fixed])

    n = viewer.align_to_axis(
        primary_indices=chainA,
        target_axis='z',
        apply_to=chainA,  # only move chain A
    )
    coordsB_after = np.array([atoms[i].coord for i in chainB_fixed])
    print(f"\nAligned only chain A ({n} atoms), chain B unchanged:")
    print(f"  Chain B moved: {not np.allclose(coordsB_before, coordsB_after)}")
finally:
    os.unlink(tmp_path)
