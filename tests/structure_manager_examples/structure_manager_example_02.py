import os
import tempfile
from gatewizard.core.structure_manager import StructureManager

viewer = StructureManager()

# Create a minimal PDB for testing
pdb_content = """\
HEADER    TEST PROTEIN
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

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name

try:
    info = viewer.load_structure(tmp_path)
    print(f"Loaded: {info['n_atoms']} atoms, {info['n_residues']} residues")
    print(f"Chains: {info['n_chains']}, Bonds: {info['n_bonds']}")
    print(f"Title: {info.get('title', 'N/A')}")
finally:
    os.unlink(tmp_path)
