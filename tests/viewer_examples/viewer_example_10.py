import os
import tempfile
from gatewizard.core.viewer import MolecularViewer

viewer = MolecularViewer()

pdb_content = """\
ATOM      1  N   ALA A  50       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A  50       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A  50       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A  50       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A  51       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A  51       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A  51       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A  51       6.500   2.000   3.000  1.00  0.00           O
HETATM    9  O   HOH A 300      20.000  20.000  20.000  1.00  0.00           O
HETATM   10  C1  LIG B   1      30.000  30.000  30.000  1.00  0.00           C
HETATM   11  C2  LIG B   1      31.000  30.000  30.000  1.00  0.00           C
END
"""

with tempfile.NamedTemporaryFile(suffix=".pdb", mode="w", delete=False) as f:
    f.write(pdb_content)
    tmp_path = f.name
out_path = tmp_path + "_edited.pdb"

try:
    # 1. Load
    info = viewer.load_structure(tmp_path)
    print(f"Loaded: {info['n_atoms']} atoms, {info['n_chains']} chains")

    # 2. Inspect
    print(f"Chains: {viewer.get_chains()}")
    sels = viewer.auto_detect_molecules()
    for s in sels:
        print(f"  Detected: {s.name} ({len(s.atom_indices)} atoms)")

    # 3. Edit: rename chain A -> X, renumber residues
    viewer.rename_chain("A", "X")
    viewer.renumber_residues("X", 50, 51, new_start=1)

    # 4. Delete water
    water = viewer.select_by_criteria("Water")
    viewer.delete_atoms(water)

    # 5. Save
    viewer.save_pdb(out_path)
    print(f"Saved edited structure: {os.path.basename(out_path)}")

    # 6. Verify
    viewer2 = MolecularViewer()
    info2 = viewer2.load_structure(out_path)
    print(f"Verified: {info2['n_atoms']} atoms, chains: {viewer2.get_chains()}")
    residues = viewer2.get_residues("X")
    for r in residues:
        print(f"  {r['name']} {r['seq_id']}")
finally:
    os.unlink(tmp_path)
    if os.path.exists(out_path):
        os.unlink(out_path)
