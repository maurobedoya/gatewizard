import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
HETATM    4  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
HETATM    5  C1  LIG A 200      30.000  30.000  30.000  1.00  0.00           C
HETATM    6  C2  LIG A 200      31.000  30.000  30.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix='.pdb', mode='w', delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)
    selections = viewer.auto_detect_molecules()
    for sel in selections:
        print(f"Selection '{sel.name}': {len(sel.atom_indices)} atoms, "
              f"rep={sel.representation}, cs={sel.color_scheme}")
finally:
    os.unlink(tmp_path)
