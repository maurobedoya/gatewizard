import os
import tempfile
from gatewizard.core.viewer import MolecularViewer, ViewerError

viewer = MolecularViewer()

# A small structure with HELIX/SHEET records
pdb_content = """\
HEADER    TEST PROTEIN
HELIX    1   1 ALA A    1  ALA A    4  1                                   4
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00  0.00           O
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Default SS (assigned automatically at load time)
    print("SS after load (auto):")
    print(f"  {viewer.get_secondary_structure_summary()}")

    # Reassign using the heuristic method
    ss = viewer.assign_secondary_structure('heuristic')
    print(f"SS after heuristic: {ss}")

    # Reassign from PDB HELIX/SHEET records
    ss = viewer.assign_secondary_structure('pdb_records')
    print(f"SS after pdb_records: {ss}")

    # Try psique (may not be installed)
    try:
        ss = viewer.assign_secondary_structure('psique')
        print(f"SS after psique: {ss}")
    except ViewerError as e:
        print(f"psique not available: {e}")

    # Auto method (same priority as load)
    ss = viewer.assign_secondary_structure('auto')
    print(f"SS after auto: {ss}")
finally:
    os.unlink(tmp_path)
