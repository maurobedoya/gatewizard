import re
from gatewizard.core.preparation import PreparationManager
from gatewizard.utils.protein_capping import ProteinCapper

# Step 1: Add ACE/NME caps to get residue mapping
capper = ProteinCapper()
capped_file, residue_mapping = capper.remove_hydrogens_and_cap(
    input_file="protein.pdb",
    output_file="protein_capped.pdb",
    target_dir="output"
)
print(f"✓ Capped protein with {len(residue_mapping)} residues tracked")

# Step 2: Run Propka analysis on capped structure
analyzer = PreparationManager()
pka_file = analyzer.run_analysis(capped_file, output_dir="output")
summary_file = analyzer.extract_summary(pka_file, output_dir="output")
residues = analyzer.parse_summary(summary_file)
print(f"✓ Analyzed {len(residues)} ionizable residues")

# Step 3: Define custom states using ORIGINAL (uncapped) numbering
original_custom_states = {
    "ASP12": "ASH",  # Force protonate ASP at original position 12
    "GLU13": "GLH"   # Force protonate GLU at original position 13
}

# Step 4: Translate to NEW numbering after capping (all shifted by +1)
custom_states_new = {}
for orig_spec, state in original_custom_states.items():
    match = re.match(r'([A-Z]+)(\d+)(?:_([A-Z]))?', orig_spec)
    if match:
        resname, resid, chain = match.groups()
        resid = int(resid)
        chain = chain or "A"
        
        # Look up the new numbering from mapping
        original_key = (resname, chain, resid)
        if original_key in residue_mapping:
            new_resname, new_chain, new_resid = residue_mapping[original_key]
            new_spec = f"{new_resname}{new_resid}_{new_chain}"
            custom_states_new[new_spec] = state
            print(f"  Mapped {orig_spec} → {new_spec} ({state})")
        else:
            print(f"  Warning: {orig_spec} not found in mapping")

# Step 5: Apply protonation states using the TRANSLATED custom states
stats = analyzer.apply_protonation_states(
    input_pdb=capped_file,
    output_pdb="output/protein_capped_custom.pdb",
    ph=7.4,
    custom_states=custom_states_new,  # Use translated states
    residues=residues
)
print(f"✓ Applied protonation: {stats['residue_changes']} residues modified")

# Step 6: Run pdb4amber with automatic ACE/NME HETATM fix
result = analyzer.run_pdb4amber_with_cap_fix(
    input_pdb="output/protein_capped_custom.pdb",
    output_pdb="output/protein_prepared.pdb",
    fix_caps=True  # Automatically fix ACE/NME HETATM records
)

print(f"✓ pdb4amber completed: {result['output_file']}")
if result['hetatm_fixed'] > 0:
    print(f"✓ Fixed {result['hetatm_fixed']} HETATM records for ACE/NME caps")
print(f"✓ Final structure ready: protein_prepared.pdb")
