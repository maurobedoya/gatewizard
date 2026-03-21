import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A  10       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A  10       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A  11       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A  11       5.000   1.000   3.000  1.00  0.00           C
ATOM      5  N   ALA A  12       7.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  ALA A  12       8.000   1.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    # Renumber residues 10-12 to start at 1
    count = viewer.renumber_residues('A', 10, 12, new_start=1)
    print(f"Renumbered {count} atoms")
    residues = viewer.get_residues('A')
    for r in residues:
        print(f"  {r['name']} {r['seq_id']}")
finally:
    os.unlink(tmp_path)
