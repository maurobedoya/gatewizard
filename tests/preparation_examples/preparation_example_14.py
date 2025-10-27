from gatewizard.core.preparation import PreparationManager
import matplotlib.pyplot as plt
import numpy as np

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)

# Filter protein residues only
protein_residues = [r for r in residues if r['res_id'] > 0]

# Set target pH for protonation state analysis
target_ph = 5.0  # Change this to analyze different pH

# Group by residue type and calculate protonation states
residue_data = {}
for res in protein_residues:
    res_type = res['residue']
    if res_type not in residue_data:
        residue_data[res_type] = {'pka': [], 'protonated': [], 'state': []}
    
    residue_data[res_type]['pka'].append(res['pka'])
    
    # Determine protonation state at target pH
    state = analyzer.get_default_protonation_state(res, ph=target_ph)
    residue_data[res_type]['state'].append(state)
    
    # Check if protonated (pH < pKa => protonated)
    is_protonated = target_ph < res['pka']
    residue_data[res_type]['protonated'].append(is_protonated)

# Prepare data for plotting
res_names = []
pka_values_list = []
protonation_list = []
colors_list = []

color_map = {'ASP': '#e74c3c', 'GLU': "#f92f94", 'HIS': "#72e183", 
             'LYS': "#60cae7", 'CYS': '#9b59b6', 'TYR': "#5b76f8", 
             'ARG': '#f39c12', 'N+': '#95a5a6', 'C-': '#34495e'}

for res_type in sorted(residue_data.keys()):
    res_names.append(res_type)
    pka_values_list.append(residue_data[res_type]['pka'])
    protonation_list.append(residue_data[res_type]['protonated'])
    colors_list.append(color_map.get(res_type, '#7f8c8d'))

# Create figure with two subplots
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), 
                                gridspec_kw={'height_ratios': [2, 1]})

# --- Top plot: pKa distribution with protonation states ---
bp = ax1.boxplot(pka_values_list, labels=res_names, patch_artist=True,
                 showmeans=True, meanline=False,
                 medianprops={'color': 'black', 'linewidth': 2},
                 meanprops={'marker': 'D', 'markerfacecolor': 'red', 
                           'markeredgecolor': 'darkred', 'markersize': 6})

# Color each box
for patch, color in zip(bp['boxes'], colors_list):
    patch.set_facecolor(color)
    patch.set_alpha(0.6)

# Add individual points colored by protonation state
for i, (res_type, pka_vals, prot_states) in enumerate(zip(res_names, pka_values_list, protonation_list)):
    x = np.random.normal(i+1, 0.04, size=len(pka_vals))
    # Protonated = filled circles, Deprotonated = open circles
    for xi, pka, is_prot in zip(x, pka_vals, prot_states):
        if is_prot:
            ax1.scatter(xi, pka, s=50, color=colors_list[i], alpha=0.7, 
                       edgecolors='black', linewidth=1, marker='o', zorder=3)
        else:
            ax1.scatter(xi, pka, s=50, facecolors='none', edgecolors=colors_list[i], 
                       alpha=0.7, linewidth=2, marker='o', zorder=3)

# Add pH reference line
ax1.axhline(y=target_ph, color='red', linestyle='--', linewidth=2.5, 
           alpha=0.8, label=f'pH {target_ph} (target)', zorder=2)

ax1.set_ylabel('pKa', fontsize=13, fontweight='bold')
ax1.set_title(f'pKa Distribution and Protonation States at pH {target_ph}', 
              fontsize=15, fontweight='bold')
ax1.grid(True, alpha=0.3, axis='y', linestyle=':')
ax1.legend(loc='upper right', fontsize=11)
ax1.set_ylim(0, 20)

# Add legend for protonation states
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', 
           markersize=10, markeredgecolor='black', linewidth=1, 
           label='Protonated (pH < pKa)'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='none', 
           markersize=10, markeredgecolor='gray', linewidth=2, 
           label='Deprotonated (pH ≥ pKa)')
]
ax1.legend(handles=legend_elements, loc='upper left', fontsize=10)

# --- Bottom plot: Protonation fraction bar chart ---
protonation_fractions = []
for prot_states in protonation_list:
    if len(prot_states) > 0:
        fraction = sum(prot_states) / len(prot_states)
    else:
        fraction = 0
    protonation_fractions.append(fraction)

bars = ax2.bar(range(1, len(res_names)+1), protonation_fractions, 
               color=colors_list, alpha=0.7, edgecolor='black', linewidth=1.5)

# Add percentage labels on bars
for i, (bar, frac) in enumerate(zip(bars, protonation_fractions)):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.02,
             f'{frac*100:.0f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')

ax2.set_xlabel('Residue Type', fontsize=13, fontweight='bold')
ax2.set_ylabel('Protonated Fraction', fontsize=13, fontweight='bold')
ax2.set_title(f'Protonation State Distribution at pH {target_ph}', fontsize=14, fontweight='bold')
ax2.set_xticks(range(1, len(res_names)+1))
ax2.set_xticklabels(res_names)
ax2.set_ylim(0, 1.1)
ax2.axhline(y=0.5, color='gray', linestyle=':', linewidth=1, alpha=0.5)
ax2.grid(True, alpha=0.3, axis='y', linestyle=':')

plt.tight_layout()
plt.savefig('propka_pka_distribution.png', dpi=600, bbox_inches='tight')
#plt.show()

# Print statistics
print(f"\npKa Statistics and Protonation States at pH {target_ph}:")
print("=" * 80)
for res_type, data in sorted(residue_data.items()):
    pka_vals = data['pka']
    n_prot = sum(data['protonated'])
    n_total = len(pka_vals)
    prot_frac = n_prot / n_total if n_total > 0 else 0
    
    print(f"\n{res_type:4s}: n={n_total:2d}  "
          f"pKa range=[{min(pka_vals):5.2f}, {max(pka_vals):5.2f}]  "
          f"mean={np.mean(pka_vals):5.2f}")
    print(f"      Protonated: {n_prot}/{n_total} ({prot_frac*100:.0f}%)  "
          f"Deprotonated: {n_total-n_prot}/{n_total} ({(1-prot_frac)*100:.0f}%)")
    print(f"      States: {', '.join(set(data['state']))}")
