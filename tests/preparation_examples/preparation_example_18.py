from gatewizard.core.preparation import PreparationManager
from gatewizard.utils.protein_capping import ProteinCapper

# Step 1: Add ACE/NME caps
capper = ProteinCapper()
capped_file, residue_mapping = capper.remove_hydrogens_and_cap(
    input_file="protein.pdb",
    output_file="protein_capped.pdb",
    target_dir="output"
)
print(f"✓ Capped protein: {len(residue_mapping)} residues tracked")

# Step 2: Run Propka on capped structure
analyzer = PreparationManager()
pka_file = analyzer.run_analysis(capped_file, output_dir="output")
summary_file = analyzer.extract_summary(pka_file, output_dir="output")
residues = analyzer.parse_summary(summary_file)

# Step 3: Detect disulfide bonds
bonds = analyzer.detect_disulfide_bonds(capped_file)

# Step 4: Apply protonation
stats = analyzer.apply_protonation_states(
    input_pdb=capped_file,
    output_pdb="output/protein_capped_ph7.pdb",
    ph=7.4,
    residues=residues
)

# Step 5: Apply disulfide bonds
num_bonds = analyzer.apply_disulfide_bonds(
    input_pdb="output/protein_capped_ph7.pdb",
    output_pdb="output/protein_capped_ph7_ss.pdb",
    disulfide_bonds=bonds,
    auto_detect=False
)

# Step 6: Run pdb4amber with automatic ACE/NME HETATM fix
result = analyzer.run_pdb4amber_with_cap_fix(
    input_pdb="output/protein_capped_ph7_ss.pdb",
    output_pdb="output/protein_prepared.pdb",
    fix_caps=True  # Automatically fix ACE/NME HETATM records
)

print(f"✓ pdb4amber completed: {result['output_file']}")
if result['hetatm_fixed'] > 0:
    print(f"✓ Fixed {result['hetatm_fixed']} HETATM records for ACE/NME caps")
print(f"✓ Final structure ready: protein_prepared.pdb")
