from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()

# First run analysis and extract summary to generate the summary file
analyzer.run_analysis("protein.pdb")
analyzer.extract_summary("protein.pka")

# Now parse the summary file
residues = analyzer.parse_summary("protein_summary_of_prediction.txt")
print(residues)

# Separate protein residues from ligand atoms
protein_residues = [r for r in residues if r['res_id'] > 0]
ligand_atoms = [r for r in residues if r['res_id'] == 0]

print(f"Found {len(protein_residues)} ionizable protein residues")
print(f"Found {len(ligand_atoms)} ionizable ligand atoms")

# Analyze ligands
for lig in ligand_atoms:
    print(f"Ligand {lig['residue']} atom {lig['atom']} ({lig['atom_type']}): "
          f"pKa = {lig['pka']:.2f}")

# Group by ligand type
from collections import defaultdict
ligands_by_type = defaultdict(list)
for lig in ligand_atoms:
    ligands_by_type[lig['residue']].append(lig)

for ligand_name, atoms in ligands_by_type.items():
    print(f"\n{ligand_name}: {len(atoms)} ionizable atoms")
    for atom in atoms:
        print(f"  {atom['atom']:6s} pKa={atom['pka']:5.2f} ({atom['atom_type']})")
