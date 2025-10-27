from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
# Auto-detect and apply
num_bonds = analyzer.apply_disulfide_bonds(
    input_pdb="protein.pdb",
    output_pdb="protein_ss_auto.pdb"
)
print(f"✓ Applied {num_bonds} disulfide bonds")

# Use pre-detected bonds with custom threshold
bonds = analyzer.detect_disulfide_bonds("protein.pdb", distance_threshold=2.3)
num_bonds = analyzer.apply_disulfide_bonds(
    input_pdb="protein.pdb",
    output_pdb="protein_ss_custom.pdb",
    disulfide_bonds=bonds,
    auto_detect=False
)
