"""
Builder Example 20: Extract a ligand from a PDB file

Demonstrates how to extract a specific ligand into its own PDB file
for individual parametrization.
"""

from pathlib import Path
from gatewizard.tools.ligand_parametrization import detect_ligands, extract_ligand_pdb

pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_extraction"

# First detect ligands
ligands = detect_ligands(pdb_file)
print(f"Detected ligands: {[l.name for l in ligands]}")

# Extract each ligand to its own subdirectory
for lig in ligands:
    lig_dir = str(Path(output_dir) / lig.name)
    extracted_pdb = extract_ligand_pdb(pdb_file, lig.name, lig_dir)

    # Verify extraction
    with open(extracted_pdb) as f:
        lines = [l for l in f if l.startswith("HETATM")]

    print(f"\nExtracted {lig.name}:")
    print(f"  Output: {extracted_pdb}")
    print(f"  Atoms: {len(lines)}")
    print(f"  First line: {lines[0].strip()[:60]}...")

print(f"\nAll extracted ligands saved in: {output_dir}")
