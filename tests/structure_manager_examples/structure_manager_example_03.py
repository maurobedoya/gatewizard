import os
import tempfile
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00  0.00           O
ATOM      9  N   ALA B   1      11.000   2.000   3.000  1.00  0.00           N
ATOM     10  CA  ALA B   1      12.000   2.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    chains = viewer.get_chains()
    print(f"Chains: {chains}")

    residues = viewer.get_residues(chain_id="A")
    print(f"Chain A residues: {len(residues)}")
    for r in residues:
        print(f"  {r['name']} {r['seq_id']} ({r['n_atoms']} atoms, SS: {r['ss']})")

    ss = viewer.get_secondary_structure_summary()
    print(f"SS summary: {ss}")
finally:
    os.unlink(tmp_path)
