from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()

# First run analysis to generate the .pka file
analyzer.run_analysis("protein.pdb")

# Then extract summary
summary_file = analyzer.extract_summary("protein.pka")
# Returns: "protein_summary_of_prediction.txt"