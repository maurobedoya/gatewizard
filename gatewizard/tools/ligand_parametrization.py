# gatewizard/tools/ligand_parametrization.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Ligand parametrization tools for membrane systems.

This module provides functionality to:
- Detect non-standard (ligand) residues in PDB files
- Parametrize ligands using antechamber + parmchk2 + tleap (GAFF2)
- Generate 2D structure images of ligands using RDKit
- Produce .frcmod and .lib files for use with packmol-memgen and tleap

The workflow follows AMBER conventions:
1. Extract ligand from PDB
2. Run antechamber for atom typing and charge assignment
3. Run parmchk2 for missing parameter generation
4. Run tleap to generate .lib (library) file
5. Pass .frcmod/.lib to packmol-memgen and final tleap parametrization
"""

import os
import json
import shutil
import subprocess
import tempfile
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional, Tuple

from gatewizard.utils.logger import get_logger
from gatewizard.utils.helpers import get_clean_env

logger = get_logger(__name__)

# Standard amino acid residue names (not ligands)
STANDARD_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Alternative protonation states
    'HIE', 'HID', 'HIP', 'ASH', 'GLH', 'CYX', 'CYM', 'LYN', 'TYM',
    # Capping groups
    'ACE', 'NME', 'NHE',
    # Common solvent/ions
    'WAT', 'HOH', 'TIP', 'TIP3', 'TP3', 'SOL',
    'NA', 'NA+', 'CL', 'CL-', 'K', 'K+', 'MG', 'CA', 'ZN', 'FE',
    # Lipids (common AMBER lipid residue names)
    'POPC', 'POPE', 'POPS', 'DPPC', 'DPPE', 'DMPC', 'DOPC', 'CHL1',
    'PA', 'PC', 'PE', 'OL', 'MY', 'ST', 'AR',
}

# Charge methods supported by antechamber
CHARGE_METHODS = {
    'bcc': 'AM1-BCC (recommended)',
    'resp': 'RESP Fitting',
    'cm2': 'CM2 Charges',
    'mul': 'Mulliken Charges',
    'rc': 'RC Charges',
    'esp': 'ESP Fitting',
    'gas': 'Gasteiger Charges',
}

DEFAULT_CHARGE_METHOD = 'bcc'


class LigandParametrizationError(Exception):
    """Custom exception for ligand parametrization errors."""
    pass


class LigandInfo:
    """Information about a detected ligand in a PDB file."""

    def __init__(self, name: str, chain: str, res_id: int, num_atoms: int,
                 elements: Dict[str, int], pdb_lines: List[str]):
        self.name = name
        self.chain = chain
        self.res_id = res_id
        self.num_atoms = num_atoms
        self.elements = elements  # e.g., {'C': 10, 'H': 12, 'O': 2, 'N': 1}
        self.pdb_lines = pdb_lines

    @property
    def formula(self) -> str:
        """Return molecular formula string."""
        # Standard ordering: C, H, then alphabetical
        parts = []
        for elem in ['C', 'H']:
            if elem in self.elements:
                count = self.elements[elem]
                parts.append(f"{elem}{count}" if count > 1 else elem)
        for elem in sorted(self.elements.keys()):
            if elem not in ['C', 'H']:
                count = self.elements[elem]
                parts.append(f"{elem}{count}" if count > 1 else elem)
        return ''.join(parts)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return {
            'name': self.name,
            'chain': self.chain,
            'res_id': self.res_id,
            'num_atoms': self.num_atoms,
            'elements': self.elements,
            'formula': self.formula,
        }


def detect_ligands(pdb_file: str) -> List[LigandInfo]:
    """
    Detect non-standard (ligand) residues in a PDB file.

    Scans HETATM records for residues that are not standard amino acids,
    water, ions, or common lipid components.

    Args:
        pdb_file: Path to the PDB file to analyze

    Returns:
        List of LigandInfo objects for each unique ligand found

    Raises:
        LigandParametrizationError: If the PDB file cannot be read

    Example:
        >>> ligands = detect_ligands("system.pdb")
        >>> for lig in ligands:
        ...     print(f"{lig.name}: {lig.num_atoms} atoms, formula={lig.formula}")
    """
    pdb_path = Path(pdb_file)
    if not pdb_path.exists():
        raise LigandParametrizationError(f"PDB file not found: {pdb_file}")

    # Group HETATM lines by (residue_name, chain, res_id)
    ligand_groups: Dict[Tuple[str, str, int], List[str]] = {}

    try:
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    chain = line[21:22].strip()
                    try:
                        res_id = int(line[22:26].strip())
                    except ValueError:
                        continue

                    if res_name.upper() not in STANDARD_RESIDUES:
                        key = (res_name, chain, res_id)
                        if key not in ligand_groups:
                            ligand_groups[key] = []
                        ligand_groups[key].append(line)
    except Exception as e:
        raise LigandParametrizationError(f"Error reading PDB file: {e}")

    # Build LigandInfo objects
    ligands = []
    seen_names = set()

    for (res_name, chain, res_id), pdb_lines in ligand_groups.items():
        # Count elements
        elements: Dict[str, int] = {}
        for line in pdb_lines:
            # Element is in columns 77-78 of PDB format
            element = line[76:78].strip()
            if not element:
                # Fallback: guess from atom name
                atom_name = line[12:16].strip()
                element = ''.join(c for c in atom_name if c.isalpha())[:2]
                if len(element) > 1:
                    element = element[0].upper() + element[1].lower()
                else:
                    element = element.upper()
            elements[element] = elements.get(element, 0) + 1

        # Only add unique ligand names (avoid duplicates from multiple chains)
        lig_key = res_name
        if lig_key not in seen_names:
            seen_names.add(lig_key)
            ligands.append(LigandInfo(
                name=res_name,
                chain=chain,
                res_id=res_id,
                num_atoms=len(pdb_lines),
                elements=elements,
                pdb_lines=pdb_lines
            ))

    logger.info(f"Detected {len(ligands)} ligand(s) in {pdb_file}: "
                f"{[l.name for l in ligands]}")
    return ligands


def extract_ligand_pdb(pdb_file: str, ligand_name: str, output_dir: str) -> str:
    """
    Extract a ligand from a PDB file into its own PDB file.

    Args:
        pdb_file: Path to the source PDB file
        ligand_name: 3-letter residue name of the ligand
        output_dir: Directory to write the extracted PDB

    Returns:
        Path to the extracted ligand PDB file

    Raises:
        LigandParametrizationError: If extraction fails
    """
    pdb_path = Path(pdb_file).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    output_pdb = out_dir / f"{ligand_name}.pdb"

    try:
        lines = []
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith('HETATM'):
                    res_name = line[17:20].strip()
                    if res_name == ligand_name:
                        lines.append(line)

        if not lines:
            raise LigandParametrizationError(
                f"Ligand '{ligand_name}' not found in {pdb_file}")

        with open(output_pdb, 'w') as f:
            for line in lines:
                f.write(line)
            f.write("END\n")

        logger.info(f"Extracted ligand {ligand_name} ({len(lines)} atoms) to {output_pdb}")
        return str(output_pdb)

    except LigandParametrizationError:
        raise
    except Exception as e:
        raise LigandParametrizationError(f"Error extracting ligand: {e}")


def parametrize_ligand(
    ligand_pdb: str,
    ligand_name: str,
    output_dir: str,
    charge: int = 0,
    charge_method: str = DEFAULT_CHARGE_METHOD,
    multiplicity: int = 1,
) -> Dict[str, str]:
    """
    Parametrize a single ligand using antechamber + parmchk2 + tleap.

    This follows the AMBER/GAFF2 workflow:
    1. antechamber: atom typing (GAFF2) and charge assignment
    2. parmchk2: missing parameter generation
    3. tleap: generate .lib file with residue library

    CRITICAL: The tleap variable name MUST match the residue name
    (e.g., ``AAA = loadmol2 AAA.mol2``), otherwise ``saveoff`` will store
    the unit under the wrong name and tleap will fail to recognize the
    residue when loading the full system PDB.

    Args:
        ligand_pdb: Path to the ligand PDB file
        ligand_name: 3-letter residue name (must match PDB residue name)
        output_dir: Directory for output files
        charge: Net charge of the ligand (default: 0)
        charge_method: Charge calculation method (default: 'bcc')
        multiplicity: Spin multiplicity (default: 1)

    Returns:
        Dictionary with paths to generated files:
        - 'mol2': path to typed MOL2 file
        - 'frcmod': path to force field modification file
        - 'lib': path to residue library file
        - 'prmtop': path to topology file
        - 'inpcrd': path to coordinate file

    Raises:
        LigandParametrizationError: If any step fails

    Example:
        >>> files = parametrize_ligand("AAA.pdb", "AAA", "./AAA/", charge=0)
        >>> print(files['frcmod'])  # AAA/AAA.frcmod
        >>> print(files['lib'])     # AAA/AAA.lib
    """
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Resolve input path so it works regardless of subprocess cwd
    ligand_pdb = str(Path(ligand_pdb).resolve())

    # Create logs directory
    log_dir = out_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    # Define output file paths
    mol2_file = out_dir / f"{ligand_name}.mol2"
    frcmod_file = out_dir / f"{ligand_name}.frcmod"
    lib_file = out_dir / f"{ligand_name}.lib"
    prmtop_file = out_dir / f"{ligand_name}.prmtop"
    inpcrd_file = out_dir / f"{ligand_name}.inpcrd"

    # Status tracking
    status_file = out_dir / "status.json"
    status = {
        "ligand_name": ligand_name,
        "charge": charge,
        "charge_method": charge_method,
        "status": "running",
        "current_step": "antechamber",
        "steps_completed": [],
        "start_time": datetime.now().isoformat(),
        "error": None,
    }
    _write_status(status_file, status)

    try:
        # Step 1: Antechamber - atom typing and charge assignment
        logger.info(f"Running antechamber for {ligand_name} "
                     f"(charge={charge}, method={charge_method})")

        antechamber_cmd = [
            'antechamber',
            '-i', str(ligand_pdb),
            '-fi', 'pdb',
            '-o', str(mol2_file),
            '-fo', 'mol2',
            '-c', charge_method,
            '-nc', str(charge),
            '-m', str(multiplicity),
            '-rn', ligand_name,
            '-s', '2',
            '-at', 'gaff2',
        ]

        result = subprocess.run(
            antechamber_cmd,
            capture_output=True,
            text=True,
            cwd=str(out_dir),
            timeout=600,  # 10 min timeout
            env=get_clean_env(),
        )

        # Save log
        with open(log_dir / "antechamber.log", 'w') as f:
            f.write(f"COMMAND: {' '.join(antechamber_cmd)}\n\n")
            f.write(f"STDOUT:\n{result.stdout}\n\n")
            f.write(f"STDERR:\n{result.stderr}\n")

        if result.returncode != 0:
            raise LigandParametrizationError(
                f"Antechamber failed for {ligand_name}: {result.stderr}")

        if not mol2_file.exists():
            raise LigandParametrizationError(
                f"Antechamber did not produce {mol2_file.name}")

        status['steps_completed'].append('antechamber')
        status['current_step'] = 'parmchk2'
        _write_status(status_file, status)
        logger.info(f"Antechamber completed for {ligand_name}")

        # Step 2: Parmchk2 - missing parameter generation
        logger.info(f"Running parmchk2 for {ligand_name}")

        parmchk_cmd = [
            'parmchk2',
            '-i', str(mol2_file),
            '-f', 'mol2',
            '-o', str(frcmod_file),
            '-s', 'gaff2',
        ]

        result = subprocess.run(
            parmchk_cmd,
            capture_output=True,
            text=True,
            cwd=str(out_dir),
            timeout=120,
            env=get_clean_env(),
        )

        with open(log_dir / "parmchk2.log", 'w') as f:
            f.write(f"COMMAND: {' '.join(parmchk_cmd)}\n\n")
            f.write(f"STDOUT:\n{result.stdout}\n\n")
            f.write(f"STDERR:\n{result.stderr}\n")

        if result.returncode != 0:
            raise LigandParametrizationError(
                f"Parmchk2 failed for {ligand_name}: {result.stderr}")

        if not frcmod_file.exists():
            raise LigandParametrizationError(
                f"Parmchk2 did not produce {frcmod_file.name}")

        status['steps_completed'].append('parmchk2')
        status['current_step'] = 'tleap'
        _write_status(status_file, status)
        logger.info(f"Parmchk2 completed for {ligand_name}")

        # Step 3: tleap - generate .lib file
        # CRITICAL: variable name must match residue name
        logger.info(f"Running tleap for {ligand_name}")

        tleap_input = out_dir / "tleap.in"
        tleap_content = f"""source leaprc.gaff2
loadamberparams {ligand_name}.frcmod
{ligand_name} = loadmol2 {ligand_name}.mol2
check {ligand_name}
saveoff {ligand_name} {ligand_name}.lib
saveamberparm {ligand_name} {ligand_name}.prmtop {ligand_name}.inpcrd
quit
"""
        with open(tleap_input, 'w') as f:
            f.write(tleap_content)

        result = subprocess.run(
            ['tleap', '-f', 'tleap.in'],
            capture_output=True,
            text=True,
            cwd=str(out_dir),
            timeout=120,
            env=get_clean_env(),
        )

        with open(log_dir / "tleap.log", 'w') as f:
            f.write(f"COMMAND: tleap -f tleap.in\n\n")
            f.write(f"INPUT:\n{tleap_content}\n\n")
            f.write(f"STDOUT:\n{result.stdout}\n\n")
            f.write(f"STDERR:\n{result.stderr}\n")

        if result.returncode != 0:
            raise LigandParametrizationError(
                f"tleap failed for {ligand_name}: {result.stderr}")

        if not lib_file.exists():
            raise LigandParametrizationError(
                f"tleap did not produce {lib_file.name}")

        status['steps_completed'].append('tleap')
        status['current_step'] = 'completed'
        status['status'] = 'completed'
        status['end_time'] = datetime.now().isoformat()
        _write_status(status_file, status)
        logger.info(f"Ligand {ligand_name} parametrization completed successfully")

        return {
            'mol2': str(mol2_file),
            'frcmod': str(frcmod_file),
            'lib': str(lib_file),
            'prmtop': str(prmtop_file),
            'inpcrd': str(inpcrd_file),
        }

    except LigandParametrizationError:
        status['status'] = 'error'
        status['error'] = str(status.get('error', ''))
        status['end_time'] = datetime.now().isoformat()
        _write_status(status_file, status)
        raise
    except subprocess.TimeoutExpired as e:
        msg = f"Timeout during parametrization of {ligand_name} at step {status['current_step']}"
        logger.error(msg)
        status['status'] = 'error'
        status['error'] = msg
        status['end_time'] = datetime.now().isoformat()
        _write_status(status_file, status)
        raise LigandParametrizationError(msg) from e
    except Exception as e:
        msg = f"Unexpected error during parametrization of {ligand_name}: {e}"
        logger.error(msg, exc_info=True)
        status['status'] = 'error'
        status['error'] = msg
        status['end_time'] = datetime.now().isoformat()
        _write_status(status_file, status)
        raise LigandParametrizationError(msg) from e


def parametrize_all_ligands(
    pdb_file: str,
    output_dir: str,
    charges: Optional[Dict[str, int]] = None,
    charge_method: str = DEFAULT_CHARGE_METHOD,
) -> Dict[str, Dict[str, str]]:
    """
    Detect and parametrize all ligands in a PDB file.

    Args:
        pdb_file: Path to the PDB file containing ligands
        output_dir: Base output directory (each ligand gets a subdirectory)
        charges: Dictionary mapping ligand names to net charges
                 (default: 0 for all)
        charge_method: Charge method for antechamber (default: 'bcc')

    Returns:
        Dictionary mapping ligand names to their output file paths

    Raises:
        LigandParametrizationError: If detection or parametrization fails

    Example:
        >>> results = parametrize_all_ligands(
        ...     "system.pdb", "./ligands/",
        ...     charges={"AAA": 0, "BBB": -1},
        ...     charge_method="bcc"
        ... )
        >>> for name, files in results.items():
        ...     print(f"{name}: {files['frcmod']}")
    """
    if charges is None:
        charges = {}

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Detect ligands
    ligands = detect_ligands(pdb_file)

    if not ligands:
        logger.info("No ligands detected in PDB file")
        return {}

    results = {}

    for ligand in ligands:
        lig_dir = out_dir / ligand.name
        charge = charges.get(ligand.name, 0)

        # Extract ligand PDB
        lig_pdb = extract_ligand_pdb(pdb_file, ligand.name, str(lig_dir))

        # Parametrize
        files = parametrize_ligand(
            ligand_pdb=lig_pdb,
            ligand_name=ligand.name,
            output_dir=str(lig_dir),
            charge=charge,
            charge_method=charge_method,
        )

        results[ligand.name] = files

    return results


def get_ligand_2d_image(
    ligand_pdb_or_mol2: str,
    output_image: str,
    width: int = 400,
    height: int = 300,
    *,
    remove_nonpolar_h: bool = True,
    remove_all_h: bool = False,
    dpi: int = 150,
    bond_line_width: float = 2.5,
    atom_label_font_size: int = 0,
    background_color: Tuple = (0.11, 0.11, 0.11, 1.0),
    padding: float = 0.15,
    kekulize: bool = True,
    wedge_bonds: bool = True,
    atom_palette: Optional[Dict[int, Tuple]] = None,
    highlight_atoms: Optional[List[int]] = None,
    highlight_color: Tuple = (1.0, 0.8, 0.0, 0.3),
    transparent_background: bool = False,
) -> Optional[str]:
    """
    Generate a 2D structure image of a ligand using RDKit.

    Produces a publication-quality 2D depiction from a PDB or MOL2 file
    with extensive customisation options.

    Args:
        ligand_pdb_or_mol2: Path to a PDB or MOL2 file containing the ligand.
        output_image: Path to the output PNG image.
        width: Image width in pixels (before DPI scaling).
        height: Image height in pixels (before DPI scaling).
        remove_nonpolar_h: Remove non-polar hydrogens (C-H) for a cleaner
            look while keeping polar ones (O-H, N-H, S-H …).  Default ``True``.
        remove_all_h: Remove **all** explicit hydrogens.  Overrides
            *remove_nonpolar_h* when ``True``.  Default ``False``.
        dpi: Dots-per-inch multiplier applied to *width* / *height*.
            A value of 300 doubles the pixel count relative to the
            default 150 and produces sharper images for print.
        bond_line_width: Thickness of bond lines in the drawing.
        atom_label_font_size: Font size for atom labels.  ``0`` lets
            RDKit choose automatically.
        background_color: RGBA tuple ``(r, g, b, a)`` in 0-1 range.
            Ignored when *transparent_background* is ``True``.
        padding: Fractional padding around the molecule (0-1).
        kekulize: Draw aromatic bonds as alternating single/double
            instead of dashed circles.
        wedge_bonds: Draw stereo wedge/dash bonds.
        atom_palette: Custom mapping of atomic number → RGB tuple.
            Falls back to a built-in dark-background palette when ``None``.
        highlight_atoms: List of 0-based atom indices to highlight.
        highlight_color: RGBA colour for the atom highlight circles.
        transparent_background: Use a transparent PNG background.

    Returns:
        Path to the generated image, or ``None`` if generation failed.

    Example:
        >>> # Minimal usage
        >>> get_ligand_2d_image("AAA.pdb", "AAA_2d.png")

        >>> # Publication-quality, no hydrogens, high DPI
        >>> get_ligand_2d_image(
        ...     "AAA.pdb", "AAA_hires.png",
        ...     width=800, height=600,
        ...     remove_all_h=True, dpi=300,
        ...     background_color=(1, 1, 1, 1),  # white
        ...     bond_line_width=1.5,
        ... )

        >>> # Keep only polar H, transparent background
        >>> get_ligand_2d_image(
        ...     "AAA.mol2", "AAA_trans.png",
        ...     remove_nonpolar_h=True,
        ...     transparent_background=True,
        ... )
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, AllChem, rdDepictor, rdmolops

        file_path = Path(ligand_pdb_or_mol2)

        # Load molecule based on file type
        mol = None
        if file_path.suffix.lower() == '.mol2':
            mol = Chem.MolFromMol2File(str(file_path), removeHs=False)
        elif file_path.suffix.lower() == '.pdb':
            mol = Chem.MolFromPDBFile(str(file_path), removeHs=False)

        if mol is None:
            logger.warning(f"Could not load molecule from {ligand_pdb_or_mol2}")
            return None

        # --- Hydrogen handling ---
        if remove_all_h:
            mol = rdmolops.RemoveAllHs(mol)
        elif remove_nonpolar_h:
            mol = _remove_nonpolar_hydrogens(mol)

        # --- Kekulize ---
        if kekulize:
            try:
                Chem.Kekulize(mol, clearAromaticFlags=False)
            except Exception:
                pass  # some molecules cannot be kekulized

        # --- 2D coordinates ---
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol, clearConfs=True)

        # --- DPI scaling ---
        scale = dpi / 150.0
        render_w = int(width * scale)
        render_h = int(height * scale)

        # --- Drawer ---
        drawer = Draw.rdMolDraw2D.MolDraw2DCairo(render_w, render_h)
        opts = drawer.drawOptions()
        opts.clearBackground = not transparent_background

        if not transparent_background:
            opts.backgroundColour = background_color

        opts.bondLineWidth = bond_line_width * scale
        opts.padding = padding
        opts.addStereoAnnotation = wedge_bonds

        if atom_label_font_size > 0:
            opts.maxFontSize = int(atom_label_font_size * scale)
            opts.minFontSize = int(atom_label_font_size * scale * 0.6)

        # --- Atom palette ---
        palette = atom_palette or _DEFAULT_DARK_PALETTE
        try:
            opts.setAtomPalette(palette)
        except Exception:
            pass  # older RDKit versions

        # --- Highlight ---
        highlight_atom_list = highlight_atoms or []
        highlight_atom_colors = {}
        if highlight_atom_list:
            for idx in highlight_atom_list:
                highlight_atom_colors[idx] = highlight_color

        drawer.DrawMolecule(
            mol,
            highlightAtoms=highlight_atom_list,
            highlightAtomColors=highlight_atom_colors if highlight_atom_list else {},
        )
        drawer.FinishDrawing()

        # --- Write to file ---
        png_data = drawer.GetDrawingText()

        # If transparent background requested, use PIL to strip bg
        if transparent_background:
            try:
                from PIL import Image as PILImage
                import io
                img = PILImage.open(io.BytesIO(png_data)).convert("RGBA")
                # Replace background colour pixels with transparent
                data = img.getdata()
                bg = tuple(int(c * 255) for c in background_color[:3]) + (255,)
                new_data = []
                for item in data:
                    # Treat near-black as the drawn background
                    if item[0] < 35 and item[1] < 35 and item[2] < 35:
                        new_data.append((0, 0, 0, 0))
                    else:
                        new_data.append(item)
                img.putdata(new_data)
                img.save(str(output_image))
            except ImportError:
                # Fallback without transparency
                with open(str(output_image), 'wb') as f:
                    f.write(png_data)
        else:
            with open(str(output_image), 'wb') as f:
                f.write(png_data)

        logger.info(f"Generated 2D image: {output_image} "
                     f"({render_w}x{render_h}px, dpi={dpi})")
        return str(output_image)

    except Exception as e:
        logger.warning(f"Error generating 2D image: {e}")
        return None


# Built-in colour palette optimised for dark backgrounds
_DEFAULT_DARK_PALETTE: Dict[int, Tuple] = {
    6:  (0.9, 0.9, 0.9),    # C:  light gray
    7:  (0.3, 0.5, 1.0),    # N:  blue
    8:  (1.0, 0.3, 0.3),    # O:  red
    1:  (0.7, 0.7, 0.7),    # H:  gray
    16: (1.0, 1.0, 0.0),    # S:  yellow
    15: (1.0, 0.5, 0.0),    # P:  orange
    9:  (0.0, 1.0, 0.0),    # F:  green
    17: (0.0, 0.8, 0.0),    # Cl: green
    35: (0.6, 0.2, 0.2),    # Br: brown
    53: (0.5, 0.0, 0.5),    # I:  purple
}

# Light-background palette (convenience for users)
LIGHT_PALETTE: Dict[int, Tuple] = {
    6:  (0.2, 0.2, 0.2),    # C:  dark gray
    7:  (0.0, 0.0, 0.8),    # N:  blue
    8:  (0.8, 0.0, 0.0),    # O:  red
    1:  (0.5, 0.5, 0.5),    # H:  gray
    16: (0.8, 0.6, 0.0),    # S:  yellow-brown
    15: (0.8, 0.3, 0.0),    # P:  orange
    9:  (0.0, 0.6, 0.0),    # F:  green
    17: (0.0, 0.5, 0.0),    # Cl: green
    35: (0.5, 0.1, 0.1),    # Br: brown
    53: (0.4, 0.0, 0.4),    # I:  purple
}


def _remove_nonpolar_hydrogens(mol):
    """Remove hydrogens bonded only to carbon, keeping polar H (O-H, N-H, etc.)."""
    from rdkit import Chem

    remove_ids = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 1:
            continue
        # Check neighbour: if bonded to carbon, mark for removal
        neighbours = atom.GetNeighbors()
        if neighbours and all(n.GetAtomicNum() == 6 for n in neighbours):
            remove_ids.append(atom.GetIdx())

    if not remove_ids:
        return mol

    em = Chem.RWMol(mol)
    # Remove in reverse order to keep indices valid
    for idx in sorted(remove_ids, reverse=True):
        em.RemoveAtom(idx)
    return em.GetMol()


def get_ligand_2d_image_from_pdb_lines(
    pdb_lines: List[str],
    output_image: str,
    width: int = 400,
    height: int = 300,
    **kwargs,
) -> Optional[str]:
    """
    Generate a 2D structure image from PDB HETATM lines.

    Creates a temporary PDB file from the lines and generates a 2D image.
    All keyword arguments are forwarded to :func:`get_ligand_2d_image`.

    Args:
        pdb_lines: List of PDB HETATM lines for the ligand.
        output_image: Path to the output PNG image.
        width: Image width in pixels.
        height: Image height in pixels.
        **kwargs: Additional options forwarded to ``get_ligand_2d_image``
            (e.g. *remove_nonpolar_h*, *dpi*, *background_color*, etc.).

    Returns:
        Path to the generated image, or ``None`` if failed.
    """
    try:
        with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb',
                                          delete=False) as tmp:
            for line in pdb_lines:
                tmp.write(line)
            tmp.write("END\n")
            tmp_path = tmp.name

        result = get_ligand_2d_image(
            tmp_path, output_image, width, height, **kwargs
        )
        return result
    except Exception as e:
        logger.warning(f"Error generating 2D image from PDB lines: {e}")
        return None
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass


def build_ligand_param_args(
    ligand_files: Dict[str, Dict[str, str]]
) -> List[str]:
    """
    Build ``--ligand_param`` arguments for packmol-memgen.

    Each ligand requires a separate ``--ligand_param frcmod:lib`` flag.

    Args:
        ligand_files: Dict mapping ligand name to file paths
                      (as returned by parametrize_ligand)

    Returns:
        List of command-line arguments for packmol-memgen

    Example:
        >>> args = build_ligand_param_args({
        ...     'AAA': {'frcmod': 'AAA/AAA.frcmod', 'lib': 'AAA/AAA.lib'},
        ...     'BBB': {'frcmod': 'BBB/BBB.frcmod', 'lib': 'BBB/BBB.lib'},
        ... })
        >>> print(args)
        ['--ligand_param', 'AAA/AAA.frcmod:AAA/AAA.lib',
         '--ligand_param', 'BBB/BBB.frcmod:BBB/BBB.lib']
    """
    args = []
    for name, files in ligand_files.items():
        frcmod = files.get('frcmod', '')
        lib = files.get('lib', '')
        if frcmod and lib:
            args.extend(['--ligand_param', f"{frcmod}:{lib}"])
    return args


def build_tleap_ligand_lines(
    ligand_files: Dict[str, Dict[str, str]]
) -> str:
    """
    Build tleap input lines to load ligand parameters.

    These lines should be inserted BEFORE ``loadPDB`` in the tleap input.

    Args:
        ligand_files: Dict mapping ligand name to file paths

    Returns:
        Multi-line string with tleap commands

    Example:
        >>> lines = build_tleap_ligand_lines({
        ...     'AAA': {'frcmod': 'AAA/AAA.frcmod', 'lib': 'AAA/AAA.lib'},
        ... })
        >>> print(lines)
        # Load GAFF2 and ligand parameters
        source leaprc.gaff2
        loadamberparams AAA/AAA.frcmod
        loadoff AAA/AAA.lib
    """
    if not ligand_files:
        return ""

    lines = ["# Load GAFF2 and ligand parameters", "source leaprc.gaff2"]

    for name, files in ligand_files.items():
        frcmod = files.get('frcmod', '')
        lib = files.get('lib', '')
        if frcmod:
            lines.append(f"loadamberparams {frcmod}")
        if lib:
            lines.append(f"loadoff {lib}")

    return "\n".join(lines)


def _write_status(status_file: Path, status: Dict[str, Any]) -> None:
    """Write status atomically."""
    tmp = str(status_file) + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(status, f, indent=2)
    os.replace(tmp, str(status_file))
