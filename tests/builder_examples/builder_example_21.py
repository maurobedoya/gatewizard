"""
Builder Example 21: Build packmol-memgen and tleap commands with ligand parameters

Demonstrates how to construct the packmol-memgen --ligand_param arguments
and tleap input lines for systems with ligands.
"""

from gatewizard.tools.ligand_parametrization import (
    build_ligand_param_args,
    build_tleap_ligand_lines,
)

# Simulated parametrization results (paths to .frcmod and .lib files)
ligand_files = {
    'AAA': {
        'frcmod': 'ligand_params/AAA/AAA.frcmod',
        'lib': 'ligand_params/AAA/AAA.lib',
        'mol2': 'ligand_params/AAA/AAA.mol2',
    },
    'BBB': {
        'frcmod': 'ligand_params/BBB/BBB.frcmod',
        'lib': 'ligand_params/BBB/BBB.lib',
        'mol2': 'ligand_params/BBB/BBB.mol2',
    },
}

# Build packmol-memgen arguments
# Each ligand gets its own --ligand_param flag (CANNOT combine in one flag)
pmm_args = build_ligand_param_args(ligand_files)
print("packmol-memgen arguments:")
for i in range(0, len(pmm_args), 2):
    print(f"  {pmm_args[i]} {pmm_args[i+1]}")

print()

# Build tleap input lines
# These go BEFORE loadPDB in the tleap input
tleap_lines = build_tleap_ligand_lines(ligand_files)
print("tleap input lines:")
print(tleap_lines)
