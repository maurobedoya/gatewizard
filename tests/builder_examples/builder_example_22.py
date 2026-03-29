"""
Builder Example 22: Parametrize all ligands in a PDB file

Demonstrates the full ligand parametrization workflow:
1. Detect ligands from PDB
2. Extract each ligand
3. Run antechamber + parmchk2 + tleap for each
4. Collect .frcmod and .lib files

NOTE: This example requires AmberTools (antechamber, parmchk2, tleap)
to be installed and accessible in the PATH.
"""

from pathlib import Path
from gatewizard.tools.ligand_parametrization import (
    parametrize_all_ligands,
    build_ligand_param_args,
    build_tleap_ligand_lines,
)

pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_params"

# Set charges for each ligand (default is 0 if not specified)
charges = {
    'AAA': 0,
    'BBB': 0,
}

print(f"Output directory: {output_dir}")

# Parametrize all ligands
results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=output_dir,
    charges=charges,
    charge_method='bcc',  # AM1-BCC charges
    atom_type='gaff2',    # GAFF2 atom types (recommended with abcg2; bcc also works)
)

print(f"\nParametrized {len(results)} ligand(s):")
for name, files in results.items():
    print(f"\n  {name}:")
    for file_type, file_path in files.items():
        exists = Path(file_path).exists()
        print(f"    {file_type}: {file_path} ({'exists' if exists else 'MISSING'})")

# Now these can be passed to Builder.prepare_system()
print("\n\npackmol-memgen arguments:")
pmm_args = build_ligand_param_args(results)
print(f"  {' '.join(pmm_args)}")

print("\ntleap ligand lines:")
tleap_lines = build_tleap_ligand_lines(results)
print(tleap_lines)

print(f"\nAll parameter files saved in: {output_dir}")
