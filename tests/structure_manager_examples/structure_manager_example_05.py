import os
import tempfile
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

pdb_content = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      4  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      5  N   ALA B   1      11.000   2.000   3.000  1.00  0.00           N
ATOM      6  CA  ALA B   1      12.000   2.000   3.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    viewer.load_structure(tmp_path)

    chain_a = viewer.select_by_criteria("Chain...", "A")
    print(f"Chain A atoms: {len(chain_a)}")

    rng = viewer.select_by_criteria("Residue range...", "A:1-2")
    print(f"A:1-2 atoms: {len(rng)}")

    all_atoms = viewer.select_by_criteria("All")
    print(f"All atoms: {len(all_atoms)}")
finally:
    os.unlink(tmp_path)
