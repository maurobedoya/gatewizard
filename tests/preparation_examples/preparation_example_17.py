from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()

# Step 1: Run Propka analysis
pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
summary_file = analyzer.extract_summary(pka_file, output_dir="output")
residues = analyzer.parse_summary(summary_file)
print(f"✓ Found {len(residues)} ionizable residues")

# Step 2: Detect disulfide bonds
bonds = analyzer.detect_disulfide_bonds("protein.pdb", distance_threshold=2.5)
print(f"✓ Detected {len(bonds)} disulfide bonds")

# Step 3: Apply protonation states
stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="output/protein_ph7.pdb",
    ph=7.4,
    residues=residues
)
print(f"✓ Protonation: {stats['residue_changes']} residues, "
      f"{stats['record_changes']} PDB records changed")

# Step 4: Apply disulfide bonds
num_bonds = analyzer.apply_disulfide_bonds(
    input_pdb="output/protein_ph7.pdb",
    output_pdb="output/protein_ph7_ss.pdb",
    disulfide_bonds=bonds,
    auto_detect=False
)
print(f"✓ Applied {num_bonds} disulfide bonds (CYS → CYX)")
