from gatewizard.core.preparation import PreparationManager

# Step 1: Run Propka analysis
analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb", output_dir="output")
summary_file = analyzer.extract_summary(pka_file, output_dir="output")
residues = analyzer.parse_summary(summary_file)

# Step 2: Define expected model pKa values
expected_pka = {
    'ASP': 3.9, 'GLU': 4.3, 'HIS': 6.0,
    'LYS': 10.5, 'ARG': 12.5, 'CYS': 8.3, 'TYR': 10.1
}

# Step 3: Find residues with significant pKa shifts (>1.0 units)
print("Residues with unusual pKa shifts (>1.0 units):")
print("=" * 70)

shifted_residues = []
for res in residues:
    # Only analyze protein residues (skip ligands)
    if res['res_id'] == 0:
        continue
    
    res_name = res['residue']
    if res_name in expected_pka:
        pka_diff = abs(res['pka'] - expected_pka[res_name])
        if pka_diff > 1.0:
            shifted_residues.append({
                'id': f"{res_name}{res['res_id']}_{res['chain']}",
                'type': res_name,
                'pka': res['pka'],
                'expected': expected_pka[res_name],
                'shift': res['pka'] - expected_pka[res_name]
            })
            print(f"{res_name}{res['res_id']:3d} (Chain {res['chain']}): "
                  f"pKa={res['pka']:5.2f} (expected ~{expected_pka[res_name]:.1f}, "
                  f"shift={res['pka'] - expected_pka[res_name]:+.1f})")

print(f"\nSummary: Found {len(shifted_residues)} residues with significant pKa shifts")

# Step 4: Categorize by shift direction
upshifted = [r for r in shifted_residues if r['shift'] > 1.0]
downshifted = [r for r in shifted_residues if r['shift'] < -1.0]

print(f"  Upshifted (more basic): {len(upshifted)}")
print(f"  Downshifted (more acidic): {len(downshifted)}")

# Step 5: Find extreme cases (>2.0 units shift)
extreme_shifts = [r for r in shifted_residues if abs(r['shift']) > 2.0]
if extreme_shifts:
    print(f"\n⚠ Extreme shifts (>2.0 units): {len(extreme_shifts)}")
    for r in extreme_shifts:
        print(f"  {r['id']}: {r['shift']:+.2f} units")
