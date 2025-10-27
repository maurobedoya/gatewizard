from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)

stats = analyzer.get_residue_statistics()
for res_type, count in stats.items():
    print(f"{res_type}: {count}")
