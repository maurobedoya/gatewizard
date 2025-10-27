from gatewizard.utils.protein_capping import ProteinCapper

capper = ProteinCapper()

# Add caps and get mapping
capped_file, residue_mapping = capper.remove_hydrogens_and_cap(
    input_file="protein.pdb",
    output_file="protein_capped.pdb",
    target_dir="output"
)

print(f"✓ Capped protein: {capped_file}")
print(f"✓ Residue mapping: {len(residue_mapping)} residues tracked")

# Show first few mappings
for (orig_resname, chain, orig_resid), (new_resname, new_chain, new_resid) in list(residue_mapping.items())[:5]:
    print(f"  {orig_resname}{orig_resid} (Chain {chain}) → {new_resname}{new_resid}")
# Output:
#  LEU1 (Chain A) → LEU2
#  PRO2 (Chain A) → PRO3
#  ALA3 (Chain A) → ALA4
#  GLY4 (Chain A) → GLY5
#  ILE5 (Chain A) → ILE6
