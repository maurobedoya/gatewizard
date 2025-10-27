from gatewizard.utils.protein_capping import cap_protein

capped_file, mapping = cap_protein(
    input_file="protein.pdb",
    output_file="protein_capped_convenient.pdb"
)
