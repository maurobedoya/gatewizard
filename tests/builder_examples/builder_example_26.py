"""
Builder Example 26: Atom type selection and recommended pairings

Demonstrates how to choose between GAFF and GAFF2 atom types,
check recommended pairings per the AMBER manual, and use
the ABCG2 charge method with GAFF2.

NOTE: This example requires AmberTools (antechamber, parmchk2, tleap)
to be installed and accessible in the PATH.
"""

from gatewizard.tools.ligand_parametrization import (
    parametrize_all_ligands,
    build_tleap_ligand_lines,
    ATOM_TYPES,
    CHARGE_METHODS,
    DEFAULT_ATOM_TYPE,
    DEFAULT_CHARGE_METHOD,
    RECOMMENDED_COMBOS,
    NON_RECOMMENDED_COMBOS,
)

# ── Available options ────────────────────────────────────────────────
print("Available atom types:")
for key, label in ATOM_TYPES.items():
    default = " (default)" if key == DEFAULT_ATOM_TYPE else ""
    print(f"  {key}: {label}{default}")

print("\nAvailable charge methods:")
for key, label in CHARGE_METHODS.items():
    default = " (default)" if key == DEFAULT_CHARGE_METHOD else ""
    print(f"  {key}: {label}{default}")

# ── Recommended pairings ─────────────────────────────────────────────
print("\nRecommended pairings (AMBER manual):")
for at, cm in sorted(RECOMMENDED_COMBOS):
    print(f"  {at} + {cm}  ✓")

print("\nNon-recommended pairings (will show warning):")
for at, cm in sorted(NON_RECOMMENDED_COMBOS):
    print(f"  {at} + {cm}  ✗")

# ── Check a pairing before parametrizing ─────────────────────────────
atom_type = 'gaff2'
charge_method = 'abcg2'

if (atom_type, charge_method) in RECOMMENDED_COMBOS:
    print(f"\n{atom_type}/{charge_method} is a recommended pairing.")
elif (atom_type, charge_method) in NON_RECOMMENDED_COMBOS:
    print(f"\nWARNING: {atom_type}/{charge_method} is NOT recommended.")
else:
    print(f"\n{atom_type}/{charge_method} has no specific recommendation.")

# ── Parametrize with gaff2/abcg2 (recommended) ──────────────────────
pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_params_gaff2_abcg2"

results = parametrize_all_ligands(
    pdb_file=pdb_file,
    output_dir=output_dir,
    charges={'AAA': 0, 'BBB': 0},
    charge_method='abcg2',   # ABCG2 charges
    atom_type='gaff2',       # GAFF2 atom types (recommended with abcg2)
)

print(f"\nParametrized {len(results)} ligand(s) with gaff2/abcg2:")
for name, files in results.items():
    print(f"  {name}: {files.get('frcmod', 'N/A')}")

# ── tleap lines reflect the chosen atom type ─────────────────────────
tleap_gaff2 = build_tleap_ligand_lines(results, atom_type='gaff2')
print(f"\ntleap lines (GAFF2):\n{tleap_gaff2}")

tleap_gaff = build_tleap_ligand_lines(results, atom_type='gaff')
print(f"\ntleap lines (GAFF):\n{tleap_gaff}")
