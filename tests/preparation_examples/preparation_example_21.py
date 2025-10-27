from gatewizard.core.preparation import PreparationManager

# Step 1: Run standard Propka workflow
analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
summary_file = analyzer.extract_summary(pka_file, output_dir="output")
residues = analyzer.parse_summary(summary_file)
print(f"✓ Analyzed {len(residues)} ionizable residues")

# Step 2: Specify custom states to override automatic pH-based assignment
custom_states = {
    "ASP12_A": "ASH",  # Force protonated aspartate in chain A
    "HIS15_A": "HID",  # Force delta-protonated histidine (instead of HIE or HIP)
    "GLU22": "GLH"     # Force protonated in all chains
}

# Step 3: Apply protonation with custom overrides
stats = analyzer.apply_protonation_states(
    input_pdb="protein.pdb",
    output_pdb="output/protein_custom.pdb",
    ph=7.4,
    custom_states=custom_states,
    residues=residues
)
print(f"✓ Modified {stats['residue_changes']} residues ({stats['record_changes']} atoms)")
print(f"  Custom states applied: {len(custom_states)}")

# Step 4: Verify the changes
for spec, state in custom_states.items():
    print(f"  {spec} → {state}")
