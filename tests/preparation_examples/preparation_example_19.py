from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()

# Run analysis once
pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
summary_file = analyzer.extract_summary(pka_file, output_dir="output")
residues = analyzer.parse_summary(summary_file)
bonds = analyzer.detect_disulfide_bonds("protein.pdb")

# Generate structures for different pH values
for ph in [5.0, 6.0, 7.0, 7.4, 8.0, 9.0]:
    # Apply protonation
    stats = analyzer.apply_protonation_states(
        input_pdb="protein.pdb",
        output_pdb=f"output/protein_ph{ph:.1f}.pdb",
        ph=ph,
        residues=residues
    )
    
    # Apply disulfide bonds
    analyzer.apply_disulfide_bonds(
        input_pdb=f"output/protein_ph{ph:.1f}.pdb",
        output_pdb=f"output/protein_ph{ph:.1f}_ss.pdb",
        disulfide_bonds=bonds,
        auto_detect=False
    )
    
    print(f"pH {ph:.1f}: {stats['residue_changes']} residues modified, "
          f"{len(bonds)} S-S bonds")
