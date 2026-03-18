"""
Builder Example 24: Generate 2D structure images of ligands

Demonstrates generating 2D molecular structure images from PDB files
using RDKit, including options for hydrogen removal, DPI control,
custom colour palettes, and transparent backgrounds.
"""

from pathlib import Path
from gatewizard.tools.ligand_parametrization import (
    detect_ligands,
    get_ligand_2d_image_from_pdb_lines,
    get_ligand_2d_image,
    extract_ligand_pdb,
    LIGHT_PALETTE,
)

pdb_file = "tests/2MVJ_2ligs.pdb"
output_dir = "./systems/ligand_images"
Path(output_dir).mkdir(parents=True, exist_ok=True)

# Detect ligands
ligands = detect_ligands(pdb_file)
print(f"Detected {len(ligands)} ligands\n")

for lig in ligands:
    # --- Style 1: Default (dark background, non-polar H removed) --------
    img = str(Path(output_dir) / f"{lig.name}_default.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=400, height=300,
        remove_nonpolar_h=True,           # cleaner look (default)
    )
    if result:
        print(f"{lig.name} default : {Path(result).stat().st_size:>6} bytes")

    # --- Style 2: All hydrogens removed ----------------------------------
    img = str(Path(output_dir) / f"{lig.name}_no_h.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=400, height=300,
        remove_all_h=True,                # skeleton only
    )
    if result:
        print(f"{lig.name} no-H    : {Path(result).stat().st_size:>6} bytes")

    # --- Style 3: High-DPI for publication (white background) ------------
    img = str(Path(output_dir) / f"{lig.name}_hires.png")
    lig_dir = str(Path(output_dir) / lig.name)
    extracted_pdb = extract_ligand_pdb(pdb_file, lig.name, lig_dir)
    result = get_ligand_2d_image(
        extracted_pdb, img,
        width=800, height=600,
        dpi=300,                           # high DPI
        remove_all_h=True,
        background_color=(1, 1, 1, 1),     # white
        atom_palette=LIGHT_PALETTE,        # colours for light background
        bond_line_width=1.5,
    )
    if result:
        print(f"{lig.name} hi-res  : {Path(result).stat().st_size:>6} bytes")

    # --- Style 4: Transparent background ---------------------------------
    img = str(Path(output_dir) / f"{lig.name}_transparent.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=400, height=300,
        remove_nonpolar_h=True,
        transparent_background=True,
    )
    if result:
        print(f"{lig.name} transp. : {Path(result).stat().st_size:>6} bytes")

    # --- Style 5: All hydrogens visible, thicker bonds -------------------
    img = str(Path(output_dir) / f"{lig.name}_all_h.png")
    result = get_ligand_2d_image_from_pdb_lines(
        lig.pdb_lines, img, width=500, height=400,
        remove_nonpolar_h=False,
        remove_all_h=False,
        bond_line_width=3.5,
        padding=0.2,
    )
    if result:
        print(f"{lig.name} all-H   : {Path(result).stat().st_size:>6} bytes")

    print()

print(f"All images saved in: {output_dir}")
