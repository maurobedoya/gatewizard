import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Rename chain A -> X
    count = viewer.rename_chain("A", "X")
    print(f"Renamed {count} atoms from chain A to X")
    chains = viewer.get_chains()
    print(f"Chains after rename: {chains}")

    # Rename residue
    count = viewer.rename_residues("X", 1, 1, "MET")
    print(f"Renamed {count} atoms to MET")
    residues = viewer.get_residues("X")
    for r in residues:
        print(f"  {r['name']} {r['seq_id']}")
finally:
    os.unlink(tmp_path)
