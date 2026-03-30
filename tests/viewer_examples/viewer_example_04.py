import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       2.000   3.000   4.000  1.00  0.00           C
HETATM    6  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
HETATM    7  C1  LIG A 200      30.000  30.000  30.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    protein = viewer.select_by_criteria("Protein")
    print(f"Protein atoms: {len(protein)}")

    backbone = viewer.select_by_criteria("Backbone")
    print(f"Backbone atoms: {len(backbone)}")

    sidechain = viewer.select_by_criteria("Sidechain")
    print(f"Sidechain atoms: {len(sidechain)}")

    water = viewer.select_by_criteria("Water")
    print(f"Water atoms: {len(water)}")

    ligand = viewer.select_by_criteria("Ligand")
    print(f"Ligand atoms: {len(ligand)}")
finally:
    os.unlink(tmp_path)
