from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
# Detect with default 2.5 Å threshold
bonds = analyzer.detect_disulfide_bonds("protein.pdb")
print(f"Found {len(bonds)} disulfide bonds:")

for bond in bonds:
    (res1_name, res1_id), (res2_name, res2_id) = bond
    print(f"  {res1_name}{res1_id} ↔ {res2_name}{res2_id}")
# Output:
# Found 1 disulfide bonds:
#   CYS18 ↔ CYS34

# Use stricter threshold
bonds_strict = analyzer.detect_disulfide_bonds("protein.pdb", distance_threshold=2.0)
print(f"Strict (2.0 Å): {len(bonds_strict)} bonds")
