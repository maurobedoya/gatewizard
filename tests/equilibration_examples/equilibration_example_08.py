#!/usr/bin/env python3
"""
Equilibration Example 08 — MDAnalysis Selections for Restraints

Demonstrates:
  1. Listing default MDAnalysis selections and their atom counts.
  2. Auto-detecting extra ligand residues in the PDB as named selections.
  3. Using custom MDAnalysis selection strings when generating restraint
     PDB files so that the B-factor column is filled via MDAnalysis
     rather than the built-in heuristic.

Requirements:
  - MDAnalysis >= 2.0
  - A system PDB file (the example uses the bundled popc_membrane data).
"""

from pathlib import Path
from gatewizard.tools.equilibration import NAMDEquilibrationManager

# ── 1. Point to the test system ──────────────────────────────────────────
work_dir = Path(__file__).parent / "popc_membrane"
system_pdb = work_dir / "bilayer_protein_protonated_prepared_lipid.pdb"

manager = NAMDEquilibrationManager(work_dir)

# ── 2. Inspect the default selections ────────────────────────────────────
print("=== Default MDAnalysis selections ===")
for name, sel in NAMDEquilibrationManager.DEFAULT_SELECTIONS.items():
    count = NAMDEquilibrationManager.count_selection_atoms(str(system_pdb), sel)
    print(f"  {name:25s}  →  {count:>7d} atoms")

# ── 3. Auto-detect ligands / non-standard residues ──────────────────────
print("\n=== Auto-detected selections (includes ligands) ===")
all_sels = NAMDEquilibrationManager.get_default_selections(str(system_pdb))
for name, sel in all_sels.items():
    if name.startswith("ligand_"):
        count = NAMDEquilibrationManager.count_selection_atoms(str(system_pdb), sel)
        print(f"  {name:25s}  →  {count:>7d} atoms  |  {sel}")

# ── 4. Count all selections at once ─────────────────────────────────────
print("\n=== Full atom count summary ===")
counts = NAMDEquilibrationManager.count_all_selections(str(system_pdb))
for name, count in counts.items():
    print(f"  {name:25s}  →  {count:>7d} atoms")

# ── 5. Generate a restraints PDB using MDAnalysis selections ─────────────
output_dir = work_dir / "equilibration_example_08" / "namd" / "restraints"
output_dir.mkdir(parents=True, exist_ok=True)

# Build {name: (mda_selection, force)} dict
selections_with_forces = {
    "protein_backbone": ("protein and backbone", 10.0),
    "protein_sidechain": ("protein and not backbone", 5.0),
    "lipid_head": (NAMDEquilibrationManager.DEFAULT_SELECTIONS["lipid_head"], 2.5),
    "lipid_tail": (NAMDEquilibrationManager.DEFAULT_SELECTIONS["lipid_tail"], 2.5),
    "water": (NAMDEquilibrationManager.DEFAULT_SELECTIONS["water"], 0.0),
    "ions": (NAMDEquilibrationManager.DEFAULT_SELECTIONS["ions"], 10.0),
}

# Add any auto-detected ligand with force 1.0
for name, sel in all_sels.items():
    if name.startswith("ligand_"):
        selections_with_forces[name] = (sel, 1.0)

output_file = output_dir / "step1_equilibration_restraints.pdb"
manager.generate_restraints_file_mda(
    system_pdb, selections_with_forces, output_file, stage_name="Equilibration 1"
)
print(f"\n✓ Restraints PDB written: {output_file}")

# ── 6. Use the high-level API (generate_restraints_file with selections) ─
constraints = {
    "protein_backbone": 10.0,
    "protein_sidechain": 5.0,
    "lipid_head": 2.5,
    "lipid_tail": 2.5,
    "water": 0.0,
    "ions": 10.0,
}
selections = {name: sel for name, (sel, _) in selections_with_forces.items()}

output_file2 = output_dir / "step1_via_high_level.pdb"
manager.generate_restraints_file(
    system_pdb, constraints, output_file2, stage_name="Eq1 (high-level)",
    selections=selections,
)
print(f"✓ High-level restraints PDB written: {output_file2}")

print("\nDone!")
