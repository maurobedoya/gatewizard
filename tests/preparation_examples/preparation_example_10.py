from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)
# Step 1: Apply protonation states
stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="protein_ph7.pdb",
    ph=7.4,
    residues=residues,
)

# Step 2: Apply disulfide bonds to protonated structure
num_bonds = analyzer.apply_disulfide_bonds(
    input_pdb="protein_ph7.pdb",
    output_pdb="protein_ph7_ss.pdb"
)
