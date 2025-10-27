from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)

# Basic usage - automatic protonation at pH 7.4
stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="protein_ph7.4.pdb",
    ph=7.4,
    residues=residues
)
print(f"Modified {stats['residue_changes']} residues ({stats['record_changes']} atoms)")

# With custom states - override specific residues
custom = {
    "ASP12": "ASH",    # Want to protonate original ASP12
    "GLU13": "GLH",    # Want to protonate original GLU13
    "HIS15": "HID"     # Want delta-protonated original HIS15
}

stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="protein_custom.pdb",
    ph=7.4,
    custom_states=custom,
    residues=residues
)
