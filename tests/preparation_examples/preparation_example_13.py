from gatewizard.core.preparation import PreparationManager
import matplotlib.pyplot as plt

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)

# Get titration curves
curves = analyzer.get_ph_titration_curve(ph_range=(0, 14), ph_step=0.5)

# State encoding for plotting (1 = protonated, 0 = deprotonated)
state_map = {'ASH': 1, 'ASP': 0, 'GLH': 1, 'GLU': 0, 
             'HIP': 1, 'HIE': 0, 'HID': 0, 'LYS': 1, 'LYN': 0,
             'CYS': 0, 'CYM': 0, 'TYR': 0, 'TYM': 0, 'ARG': 1}

# Option 1: Select specific residues by name
selected_residues = ['ASP12', 'GLU13', 'LYS16']
filtered_curves = {k: v for k, v in curves.items() if k in selected_residues}

# Option 2: Select by residue type (e.g., only histidines)
# filtered_curves = {k: v for k, v in curves.items() if k.startswith('HIS')}

# Option 3: Select residues with pKa near physiological pH (6-8)
# pka_dict = {f"{r['residue']}{r['res_id']}_{r['chain']}": r['pka'] 
#             for r in residues if r['res_id'] > 0}
# filtered_curves = {k: v for k, v in curves.items() 
#                    if k in pka_dict and 6 <= pka_dict[k] <= 8}

# Plot selected curves
fig, ax = plt.subplots(figsize=(10, 6))

for residue_id, curve in filtered_curves.items():
    ph_values = [ph for ph, _ in curve]
    protonation = [state_map.get(state, 0.5) for _, state in curve]
    ax.plot(ph_values, protonation, marker='o', label=residue_id, 
            linewidth=2, markersize=6, alpha=0.8)

ax.axvline(x=7.4, color='gray', linestyle='--', linewidth=1.5, 
           alpha=0.5, label='pH 7.4 (physiological)')
ax.set_xlabel('pH', fontsize=12)
ax.set_ylabel('Protonation State (1=protonated, 0=deprotonated)', fontsize=12)
ax.set_title('Protein Titration Curves', fontsize=14)
ax.set_ylim(-0.1, 1.1)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)
plt.tight_layout()
plt.savefig('propka_titration_curves.png', dpi=600, bbox_inches='tight')
#plt.show()

print(f"Plotted {len(filtered_curves)} out of {len(curves)} total residues")
