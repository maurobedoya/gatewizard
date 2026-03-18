"""
Builder Example 19: Detect ligands in a PDB file

Demonstrates how to detect non-standard residues (ligands)
in a PDB file using the ligand parametrization tools.
"""

from gatewizard.tools.ligand_parametrization import detect_ligands

# Detect ligands in a PDB file with two ligands (AAA and BBB)
pdb_file = "tests/2MVJ_2ligs.pdb"
ligands = detect_ligands(pdb_file)

print(f"Detected {len(ligands)} ligand(s) in {pdb_file}:")
print(f"{'='*60}")

for lig in ligands:
    print(f"\nLigand: {lig.name}")
    print(f"  Chain: {lig.chain}")
    print(f"  Residue ID: {lig.res_id}")
    print(f"  Number of atoms: {lig.num_atoms}")
    print(f"  Molecular formula: {lig.formula}")
    print(f"  Elements: {lig.elements}")
    print(f"  As dict: {lig.to_dict()}")
