# gatewizard/core/structure_manager.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Structure manager and editor using MDAnalysis.

Provides a programmatic API for loading, inspecting, selecting, editing,
and saving molecular structures.  Secondary structure is assigned via
PDB HELIX/SHEET records, the psique program, or a CA-angle heuristic.
"""

import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import psique

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------


class StructureError(Exception):
    """Error raised by StructureManager operations."""


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Van der Waals radii (Å). Bondi 1964 / Mantina 2009.
VDW_RADII: Dict[str, float] = {
    "H": 1.20,
    "HE": 1.40,
    "LI": 1.82,
    "BE": 1.53,
    "B": 1.92,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.35,
    "NE": 1.54,
    "NA": 2.27,
    "MG": 1.73,
    "AL": 1.84,
    "SI": 2.10,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "AR": 1.88,
    "K": 2.75,
    "CA": 2.31,
    "SC": 2.11,
    "TI": 2.10,
    "V": 2.05,
    "CR": 2.00,
    "MN": 1.80,
    "FE": 1.94,
    "CO": 1.95,
    "NI": 1.63,
    "CU": 1.40,
    "ZN": 1.39,
    "GA": 1.87,
    "GE": 2.11,
    "AS": 1.85,
    "SE": 1.90,
    "BR": 1.85,
    "KR": 2.02,
    "RB": 3.03,
    "SR": 2.49,
    "Y": 2.40,
    "ZR": 2.30,
    "NB": 2.15,
    "MO": 2.10,
    "TC": 2.10,
    "RU": 2.05,
    "RH": 2.00,
    "PD": 1.63,
    "AG": 1.72,
    "CD": 1.58,
    "IN": 1.93,
    "SN": 2.17,
    "SB": 2.06,
    "TE": 2.06,
    "I": 1.98,
    "XE": 2.16,
    "CS": 3.43,
    "BA": 2.68,
    "LA": 2.50,
    "HF": 2.25,
    "TA": 2.20,
    "W": 2.15,
    "RE": 2.10,
    "OS": 2.00,
    "IR": 2.00,
    "PT": 1.75,
    "AU": 1.66,
    "HG": 1.55,
    "TL": 1.96,
    "PB": 2.02,
    "BI": 2.07,
    "U": 1.86,
    "DEFAULT": 1.70,
}

# Covalent radii (Å). Cordero 2008 / Pyykkö 2009.
COVALENT_RADII: Dict[str, float] = {
    "H": 0.31,
    "HE": 0.46,
    "LI": 1.28,
    "BE": 0.96,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "NE": 0.58,
    "NA": 1.66,
    "MG": 1.41,
    "AL": 1.21,
    "SI": 1.11,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "AR": 1.06,
    "K": 2.03,
    "CA": 1.76,
    "SC": 1.70,
    "TI": 1.60,
    "V": 1.53,
    "CR": 1.39,
    "MN": 1.39,
    "FE": 1.32,
    "CO": 1.26,
    "NI": 1.24,
    "CU": 1.32,
    "ZN": 1.22,
    "GA": 1.22,
    "GE": 1.20,
    "AS": 1.19,
    "SE": 1.20,
    "BR": 1.20,
    "KR": 1.16,
    "RB": 2.20,
    "SR": 1.95,
    "Y": 1.90,
    "ZR": 1.75,
    "NB": 1.64,
    "MO": 1.54,
    "TC": 1.47,
    "RU": 1.46,
    "RH": 1.42,
    "PD": 1.39,
    "AG": 1.45,
    "CD": 1.44,
    "IN": 1.42,
    "SN": 1.39,
    "SB": 1.39,
    "TE": 1.38,
    "I": 1.39,
    "XE": 1.40,
    "CS": 2.44,
    "BA": 2.15,
    "LA": 2.07,
    "HF": 1.75,
    "TA": 1.70,
    "W": 1.62,
    "RE": 1.51,
    "OS": 1.44,
    "IR": 1.41,
    "PT": 1.36,
    "AU": 1.36,
    "HG": 1.32,
    "TL": 1.45,
    "PB": 1.46,
    "BI": 1.48,
    "U": 1.96,
    "DEFAULT": 1.50,
}

# Element colours – Jmol / CPK scheme (RGB 0-255).
ELEMENT_COLORS: Dict[str, Tuple[int, int, int]] = {
    "H": (255, 255, 255),
    "HE": (217, 255, 255),
    "LI": (204, 128, 255),
    "BE": (194, 255, 0),
    "B": (255, 181, 181),
    "C": (144, 144, 144),
    "N": (48, 80, 248),
    "O": (255, 13, 13),
    "F": (144, 224, 80),
    "NE": (179, 227, 245),
    "NA": (171, 92, 242),
    "MG": (138, 255, 0),
    "AL": (191, 166, 166),
    "SI": (240, 200, 160),
    "P": (255, 128, 0),
    "S": (255, 255, 48),
    "CL": (31, 240, 31),
    "AR": (128, 209, 227),
    "K": (143, 64, 212),
    "CA": (61, 255, 0),
    "SC": (230, 230, 230),
    "TI": (191, 194, 199),
    "V": (166, 166, 171),
    "CR": (138, 153, 199),
    "MN": (156, 122, 199),
    "FE": (224, 102, 51),
    "CO": (240, 144, 160),
    "NI": (80, 208, 80),
    "CU": (200, 128, 51),
    "ZN": (125, 128, 176),
    "GA": (194, 143, 143),
    "GE": (102, 143, 143),
    "AS": (189, 128, 227),
    "SE": (255, 161, 0),
    "BR": (166, 41, 41),
    "KR": (92, 184, 209),
    "RB": (112, 46, 176),
    "SR": (0, 255, 0),
    "Y": (148, 255, 255),
    "ZR": (148, 224, 224),
    "NB": (115, 194, 201),
    "MO": (84, 181, 181),
    "TC": (59, 158, 158),
    "RU": (36, 143, 143),
    "RH": (10, 125, 140),
    "PD": (0, 105, 133),
    "AG": (192, 192, 192),
    "CD": (255, 217, 143),
    "IN": (166, 117, 115),
    "SN": (102, 128, 128),
    "SB": (158, 99, 181),
    "TE": (212, 122, 0),
    "I": (148, 0, 148),
    "XE": (66, 158, 176),
    "CS": (87, 23, 143),
    "BA": (0, 201, 0),
    "LA": (112, 212, 255),
    "HF": (77, 194, 255),
    "TA": (77, 166, 255),
    "W": (33, 148, 214),
    "RE": (38, 125, 171),
    "OS": (38, 102, 150),
    "IR": (23, 84, 135),
    "PT": (208, 208, 224),
    "AU": (255, 209, 35),
    "HG": (184, 184, 208),
    "TL": (166, 84, 77),
    "PB": (87, 89, 97),
    "BI": (158, 79, 181),
    "U": (0, 143, 56),
    "DEFAULT": (255, 20, 147),
}

SS_COLORS: Dict[str, Tuple[int, int, int]] = {
    "H": (180, 141, 218),  # alpha helix – lavender
    "G": (123, 63, 181),  # 3₁₀-helix – medium violet
    "I": (61, 26, 110),  # pi helix – deep indigo
    "PP": (249, 199, 79),  # polyproline – golden yellow
    "E": (33, 150, 166),  # beta sheet – deep teal
    "C": (232, 232, 232),  # coil – light gray
    "T": (181, 213, 200),  # turn – soft sage green
    "DEFAULT": (200, 200, 200),
}

HELIX_SS_TYPES = {"H", "G", "I", "PP"}

_HELIX_CLASS_TO_SS = {
    1: "H",
    2: "H",
    3: "I",
    4: "H",
    5: "G",
    6: "H",
    7: "H",
    8: "H",
    9: "PP",
    10: "PP",
    11: "G",
    13: "I",
}

SS_LABELS: Dict[str, str] = {
    "H": "Alpha helix",
    "G": "3-10 helix",
    "I": "Pi helix",
    "PP": "Polyproline",
    "E": "Sheet",
    "C": "Coil",
    "T": "Turn",
}

CHAIN_PALETTE: List[Tuple[int, int, int]] = [
    (230, 25, 75),
    (60, 180, 75),
    (255, 225, 25),
    (0, 130, 200),
    (245, 130, 48),
    (145, 30, 180),
    (70, 240, 240),
    (240, 50, 230),
    (210, 245, 60),
    (250, 190, 212),
    (0, 128, 128),
    (220, 190, 255),
    (170, 110, 40),
    (128, 0, 0),
    (0, 0, 128),
    (128, 128, 0),
]

# Residue nature classification (7-category).
RESIDUE_NATURE: Dict[str, str] = {
    # Acidic (negatively charged)
    "ASP": "acidic",
    "GLU": "acidic",
    # Basic (positively charged)
    "ARG": "basic",
    "LYS": "basic",
    "HIS": "basic",
    # Polar uncharged
    "SER": "polar",
    "THR": "polar",
    "ASN": "polar",
    "GLN": "polar",
    # Hydrophobic aliphatic
    "ALA": "aliphatic",
    "VAL": "aliphatic",
    "LEU": "aliphatic",
    "ILE": "aliphatic",
    "MET": "aliphatic",
    # Hydrophobic aromatic
    "PHE": "aromatic",
    "TRP": "aromatic",
    "TYR": "aromatic",
    # Special
    "CYS": "special",
    "GLY": "special",
    "PRO": "special",
}

RESIDUE_NATURE_COLORS: Dict[str, Tuple[int, int, int]] = {
    "acidic": (220, 60, 60),  # Red
    "basic": (70, 100, 220),  # Blue
    "polar": (60, 180, 75),  # Green
    "aliphatic": (230, 200, 50),  # Yellow
    "aromatic": (240, 150, 50),  # Orange
    "special": (170, 80, 200),  # Purple
    "other": (180, 180, 180),  # Gray
}

RESIDUE_NATURE_LABELS: Dict[str, str] = {
    "acidic": "Acidic",
    "basic": "Basic",
    "polar": "Polar",
    "aliphatic": "Aliphatic",
    "aromatic": "Aromatic",
    "special": "Special",
    "other": "Other",
}

BACKBONE_NAMES = {"CA", "C", "N", "O", "OXT"}
AA_NAMES = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "MSE",
    "SEC",
    "PYL",
}


# ---------------------------------------------------------------------------
# Internal data classes
# ---------------------------------------------------------------------------


class Atom:
    """Lightweight atom container."""

    __slots__ = (
        "serial",
        "name",
        "element",
        "coord",
        "res_name",
        "res_id",
        "chain_id",
        "bfactor",
        "occupancy",
    )

    def __init__(
        self,
        serial,
        name,
        element,
        coord,
        res_name,
        res_id,
        chain_id,
        bfactor=0.0,
        occupancy=1.0,
    ):
        self.serial = serial
        self.name = name
        self.element = element.upper().strip()
        self.coord = np.asarray(coord, dtype=np.float64)
        self.res_name = res_name
        self.res_id = res_id
        self.chain_id = chain_id
        self.bfactor = bfactor
        self.occupancy = occupancy


class Residue:
    """Lightweight residue container."""

    __slots__ = ("name", "seq_id", "chain_id", "atoms", "ss", "ca_coord", "o_coord")

    def __init__(self, name, seq_id, chain_id):
        self.name = name
        self.seq_id = seq_id
        self.chain_id = chain_id
        self.atoms: List[Atom] = []
        self.ss = "C"
        self.ca_coord = None
        self.o_coord = None

    def add_atom(self, atom: Atom):
        self.atoms.append(atom)
        if atom.name == "CA":
            self.ca_coord = atom.coord
        elif atom.name == "O":
            self.o_coord = atom.coord


class ProteinStructure:
    """In-memory protein structure with atoms, residues, chains, and bonds."""

    def __init__(self):
        self.atoms: List[Atom] = []
        self.residues: List[Residue] = []
        self.chains: Dict[str, List[Residue]] = {}
        self.bonds: List[Tuple[int, int]] = []
        self.title: str = ""
        self.source_file: Optional[str] = None

    # -- bond detection ------------------------------------------------

    def refresh_residue_coords(self):
        """Sync residue ca_coord / o_coord with current atom coordinates."""
        for res in self.residues:
            res.ca_coord = None
            res.o_coord = None
            for a in res.atoms:
                if a.name == "CA":
                    res.ca_coord = a.coord
                elif a.name == "O":
                    res.o_coord = a.coord

    def build_bonds(self, cutoff_factor: float = 1.3):
        """Detect covalent bonds from distances and covalent radii."""
        self.refresh_residue_coords()
        coords = np.array([a.coord for a in self.atoms])
        n = len(coords)
        if n == 0:
            return
        grid: Dict[Tuple[int, int, int], List[int]] = defaultdict(list)
        cell = 3.5
        for i, c in enumerate(coords):
            key = (int(c[0] // cell), int(c[1] // cell), int(c[2] // cell))
            grid[key].append(i)
        seen: set = set()
        bonds = []
        for key, indices in grid.items():
            nbrs: List[int] = []
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        nk = (key[0] + dx, key[1] + dy, key[2] + dz)
                        if nk in grid:
                            nbrs.extend(grid[nk])
            for i in indices:
                ai = self.atoms[i]
                ri = COVALENT_RADII.get(ai.element, COVALENT_RADII["DEFAULT"])
                for j in nbrs:
                    if j <= i:
                        continue
                    pair = (i, j)
                    if pair in seen:
                        continue
                    aj = self.atoms[j]
                    rj = COVALENT_RADII.get(aj.element, COVALENT_RADII["DEFAULT"])
                    d = np.linalg.norm(coords[i] - coords[j])
                    if 0.4 < d < (ri + rj) * cutoff_factor:
                        bonds.append(pair)
                        seen.add(pair)
        self.bonds = bonds

    # -- CA-angle heuristic (fallback) ---------------------------------

    def assign_secondary_structure_heuristic(self):
        """Assign SS from CA-CA angles (simple fallback)."""
        for chain_id, residues in self.chains.items():
            ca = [r.ca_coord for r in residues]
            n = len(ca)
            for i, r in enumerate(residues):
                if i < 2 or i >= n - 2:
                    r.ss = "C"
                    continue
                cds = [ca[j] for j in range(i - 2, i + 3)]
                if any(c is None for c in cds):
                    r.ss = "C"
                    continue
                v1 = cds[1] - cds[0]
                v2 = cds[2] - cds[1]
                v3 = cds[3] - cds[2]
                a1 = _vec_angle(v1, v2)
                a2 = _vec_angle(v2, v3)
                if 80 < a1 < 120 and 80 < a2 < 120:
                    r.ss = "H"
                elif a1 > 150 or a2 > 150:
                    r.ss = "E"
                else:
                    r.ss = "C"

    # -- PDB writer ----------------------------------------------------

    def write_pdb(self, filepath: Union[str, Path]):
        """Write the structure to a PDB file."""
        filepath = Path(filepath)
        filepath.parent.mkdir(parents=True, exist_ok=True)
        with open(filepath, "w") as f:
            if self.title:
                f.write(f"TITLE     {self.title}\n")
            serial = 1
            for atom in self.atoms:
                rec = "HETATM" if atom.res_name not in AA_NAMES else "ATOM  "
                # Clamp each field to its PDB column width to prevent column
                # overflow that would corrupt the fixed-width coordinate columns
                # and cause parsers (especially MDAnalysis) to fail.
                name = str(atom.name)[:4]
                res_name = str(atom.res_name)[:3]
                chain_id = str(atom.chain_id)[:1]
                # res_id must fit in 4 chars (-999 … 9999); wrap large values
                res_id = atom.res_id % 10000
                # serial must fit in 5 chars (1–99999); wrap if needed
                ser = serial % 100000 or 100000
                f.write(
                    f"{rec}{ser:5d} {name:<4s} {res_name:>3s} "
                    f"{chain_id:1s}{res_id:4d}    "
                    f"{atom.coord[0]:8.3f}{atom.coord[1]:8.3f}{atom.coord[2]:8.3f}"
                    f"{atom.occupancy:6.2f}{atom.bfactor:6.2f}"
                    f"          {atom.element:>2s}\n"
                )
                serial += 1
            f.write("END\n")

    # -- helpers -------------------------------------------------------

    def _rebuild_residues_and_chains(self):
        """Rebuild residue/chain dicts from the atoms list."""
        self.residues.clear()
        self.chains.clear()
        cur_key = None
        cur_res = None
        for atom in self.atoms:
            rk = (atom.chain_id, atom.res_id, atom.res_name)
            if rk != cur_key:
                cur_key = rk
                cur_res = Residue(atom.res_name, atom.res_id, atom.chain_id)
                self.residues.append(cur_res)
                self.chains.setdefault(atom.chain_id, []).append(cur_res)
            cur_res.add_atom(atom)


def _vec_angle(v1, v2):
    c = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
    return np.degrees(np.arccos(np.clip(c, -1, 1)))


# ---------------------------------------------------------------------------
# Coordinate transformation helpers
# ---------------------------------------------------------------------------


def _axis_rotation_matrix(axis: str, angle_rad: float) -> np.ndarray:
    """Return a 3×3 rotation matrix for rotation around *axis* by *angle_rad*."""
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    if axis == "x":
        return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
    if axis == "y":
        return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    # z
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def _rotation_matrix_from_vectors(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Return 3×3 rotation matrix that aligns unit vector *v1* to *v2*."""
    a = v1 / (np.linalg.norm(v1) + 1e-12)
    b = v2 / (np.linalg.norm(v2) + 1e-12)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if c > 1.0 - 1e-8:  # already aligned
        return np.eye(3)
    if c < -1.0 + 1e-8:  # opposite direction — 180° rotation
        perp = np.array([1, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1, 0])
        perp = perp - np.dot(perp, a) * a
        perp /= np.linalg.norm(perp)
        return 2.0 * np.outer(perp, perp) - np.eye(3)
    s = np.linalg.norm(v)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / (s * s))


# ---------------------------------------------------------------------------
# Secondary-structure assignment helpers (PDB records → psique → heuristic)
# ---------------------------------------------------------------------------


def _read_ss_from_pdb_records(filepath: str) -> Optional[Dict]:
    """Read HELIX/SHEET records from PDB file."""
    ss_map: Dict[Tuple[str, int], str] = {}
    found = False
    with open(filepath, "r") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "HELIX":
                try:
                    chain = line[19]
                    start = int(line[21:25])
                    end = int(line[33:37])
                    helix_class = int(line[38:40]) if len(line) >= 40 else 1
                    ss_code = _HELIX_CLASS_TO_SS.get(helix_class, "H")
                    if chain.strip():
                        for seq in range(start, end + 1):
                            ss_map[(chain, seq)] = ss_code
                        found = True
                except (ValueError, IndexError):
                    continue
            elif rec == "SHEET":
                try:
                    chain = line[21]
                    start = int(line[22:26])
                    end = int(line[33:37])
                    if chain.strip():
                        for seq in range(start, end + 1):
                            ss_map[(chain, seq)] = "E"
                        found = True
                except (ValueError, IndexError):
                    continue
    return ss_map if found else None


def _assign_ss_psique(filepath: str) -> Optional[Dict]:
    """Run PSIQUE on PDB.

    Returns
    -------
    dict or None
        SS mapping, or ``None`` if PSIQUE is unavailable or produced no SS records.
    """
    try:
        ss_map: Dict[Tuple[str, int], str] = {}
        for ss in psique.assign(filepath):
            for i in range(ss.start.number, ss.end.number + 1):
                ss_map[(ss.start.chain, i)] = ss.kind.value
        return ss_map or None
    except Exception as exc:
        logger.warning("PSIQUE secondary structure assignment failed: %s", exc)
        return None


def _apply_ss_map(struct: ProteinStructure, ss_map: Dict[Tuple[str, int], str]):
    """Write an SS mapping onto *struct* residues (default ``'C'``)."""
    for r in struct.residues:
        r.ss = ss_map.get((r.chain_id, r.seq_id), "C")


def _assign_secondary_structure(
    struct: ProteinStructure, filepath: Optional[str] = None
):
    """Assign SS using best available method:
    1) PSIQUE external tool
    2) PDB HELIX/SHEET records
    3) CA-angle heuristic (fallback)
    """
    if filepath:
        ss_map = _assign_ss_psique(filepath)
        if ss_map:
            _apply_ss_map(struct, ss_map)
            return
        ss_map = _read_ss_from_pdb_records(filepath)
        if ss_map:
            _apply_ss_map(struct, ss_map)
            return
    struct.assign_secondary_structure_heuristic()


def assign_secondary_structure_map(
    filepath: str,
    method: str = "auto",
) -> Dict[Tuple[str, int], str]:
    """Assign secondary structure and return ``{(chain_id, resid): ss_code}``.

    Parameters
    ----------
    filepath : str
        Path to a coordinate file readable by MDAnalysis (typically PDB).
    method : str
        Assignment method. ``'auto'`` tries PSIQUE, then PDB HELIX/SHEET records,
        then the CA-angle heuristic. Other values match
        :meth:`StructureManager.assign_secondary_structure`.

    Returns
    -------
    dict
        Mapping from ``(chain_id, resid)`` to SS code (``'H'``, ``'E'``, ``'C'``, …).
    """
    filepath = str(filepath)
    method = method.lower()
    struct = _load_structure_from_mdanalysis(filepath, build_bonds=False)
    if method == "auto":
        _assign_secondary_structure(struct, filepath=filepath)
        return {(r.chain_id, r.seq_id): r.ss for r in struct.residues}

    if method == "psique":
        ss_map = _assign_ss_psique(filepath)
        if not ss_map:
            raise StructureError(
                "PSIQUE produced no secondary structure assignments "
                "for this structure (too few residues?)."
            )
        _apply_ss_map(struct, ss_map)
    elif method == "heuristic":
        struct.assign_secondary_structure_heuristic()
    elif method == "pdb_records":
        ss_map = _read_ss_from_pdb_records(filepath)
        if not ss_map:
            raise StructureError("No HELIX/SHEET records found in PDB file")
        _apply_ss_map(struct, ss_map)
    else:
        raise StructureError(
            f"Unknown method '{method}'. "
            "Use 'auto', 'psique', 'heuristic', or 'pdb_records'."
        )
    return {(r.chain_id, r.seq_id): r.ss for r in struct.residues}


# ---------------------------------------------------------------------------
# MDAnalysis-based PDB parser
# ---------------------------------------------------------------------------


def _load_structure_from_mdanalysis(
    filepath: str, *, build_bonds: bool = True
) -> ProteinStructure:
    """Parse a PDB (or MDAnalysis-readable) file into a ProteinStructure."""
    import MDAnalysis as mda
    import warnings

    struct = ProteinStructure()
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        u = mda.Universe(filepath)

    # Extract title from raw PDB header if available
    _header_text = ""
    with open(filepath, "r") as fh:
        for line in fh:
            if line.startswith("TITLE"):
                struct.title += line[10:].strip() + " "
            elif line.startswith("HEADER"):
                _header_text = line[10:50].strip()
            elif line.startswith(("ATOM", "HETATM")):
                break
    struct.title = struct.title.strip()
    if not struct.title and _header_text:
        struct.title = _header_text

    cur_key = None
    cur_res: Optional[Residue] = None
    for ag_atom in u.atoms:
        # MDAnalysis uses segids for chain identification.
        # CHARMM/MemPrO files populate segid (e.g. "PROT", "MEMB");
        # standard PDB files (e.g. packmol-memgen) leave segid empty and
        # store the chain letter in chainID (col 22).  Fall back accordingly.
        chain_id = (ag_atom.segid or "").strip()
        if len(chain_id) > 1:
            # Long segid (CHARMM style) – prefer the single-char chainID
            try:
                chain_id = ag_atom.chainID
            except AttributeError:
                chain_id = chain_id[0]
        if not chain_id:
            # segid was empty – try the PDB chainID column
            try:
                chain_id = (ag_atom.chainID or "").strip()
            except AttributeError:
                pass
        chain_id = chain_id or "A"
        element = (
            ag_atom.element
            if hasattr(ag_atom, "element") and ag_atom.element
            else ag_atom.name[0]
        )
        atom = Atom(
            serial=int(ag_atom.id),
            name=ag_atom.name,
            element=element.upper().strip(),
            coord=ag_atom.position.copy(),
            res_name=ag_atom.resname,
            res_id=int(ag_atom.resid),
            chain_id=chain_id,
            bfactor=(
                float(ag_atom.tempfactor) if hasattr(ag_atom, "tempfactor") else 0.0
            ),
            occupancy=(
                float(ag_atom.occupancy) if hasattr(ag_atom, "occupancy") else 1.0
            ),
        )
        struct.atoms.append(atom)
        rk = (chain_id, atom.res_id, atom.res_name)
        if rk != cur_key:
            cur_key = rk
            cur_res = Residue(atom.res_name, atom.res_id, chain_id)
            struct.residues.append(cur_res)
            struct.chains.setdefault(chain_id, []).append(cur_res)
        cur_res.add_atom(atom)

    struct.source_file = filepath
    if build_bonds:
        struct.build_bonds()
    return struct


def _parse_with_mdanalysis(filepath: str) -> ProteinStructure:
    """Parse a PDB file using MDAnalysis."""
    struct = _load_structure_from_mdanalysis(filepath, build_bonds=True)
    _assign_secondary_structure(struct, filepath=filepath)
    return struct


# ---------------------------------------------------------------------------
# Fallback manual PDB parser (no external deps)
# ---------------------------------------------------------------------------


def _parse_pdb_manual(filepath: str) -> ProteinStructure:
    """Parse PDB with no library dependency (fallback)."""
    struct = ProteinStructure()
    cur_key = None
    cur_res: Optional[Residue] = None
    _header_text = ""
    with open(filepath, "r") as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec == "TITLE":
                struct.title += line[10:].strip() + " "
            elif rec == "HEADER":
                _header_text = line[10:50].strip()
            if rec not in ("ATOM", "HETATM"):
                continue
            try:
                serial = int(line[6:11])
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[21] if line[21] != " " else "A"
                res_id = int(line[22:26])
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                occ = float(line[54:60]) if len(line) >= 60 else 1.0
                bf = float(line[60:66]) if len(line) >= 66 else 0.0
                elem = line[76:78].strip() if len(line) >= 78 else atom_name[0]
            except (ValueError, IndexError):
                continue
            atom = Atom(
                serial, atom_name, elem, (x, y, z), res_name, res_id, chain_id, bf, occ
            )
            struct.atoms.append(atom)
            rk = (chain_id, res_id, res_name)
            if rk != cur_key:
                cur_key = rk
                cur_res = Residue(res_name, res_id, chain_id)
                struct.residues.append(cur_res)
                struct.chains.setdefault(chain_id, []).append(cur_res)
            cur_res.add_atom(atom)
    struct.title = struct.title.strip()
    if not struct.title and _header_text:
        struct.title = _header_text
    struct.source_file = filepath
    struct.build_bonds()
    _assign_secondary_structure(struct, filepath=filepath)
    return struct


def parse_pdb(filepath: str) -> ProteinStructure:
    """Parse a PDB file.  Prefers MDAnalysis; falls back to manual parser."""
    try:
        return _parse_with_mdanalysis(filepath)
    except Exception:
        logger.debug("MDAnalysis not available or failed, using manual parser")
        return _parse_pdb_manual(filepath)


# ---------------------------------------------------------------------------
# Selection helper
# ---------------------------------------------------------------------------


class Selection:
    """Describes a visual selection: atom indices + rendering parameters."""

    def __init__(
        self,
        name: str,
        atom_indices: List[int],
        *,
        representation: str = "ball_stick",
        color_scheme: str = "element",
        uniform_color: Optional[Tuple[int, int, int]] = None,
        visible: bool = True,
        carbon_color: Optional[Tuple[int, int, int]] = None,
        quality: int = 3,
        opacity: float = 0.5,
        surface_resolution: int = 64,
        surface_radius: float = 0.12,
        atom_scale: float = 1.0,
        bond_radius: float = 0.15,
        ball_scale: float = 0.3,
        stick_radius: float = 0.2,
        backbone_radius: float = 0.3,
        helix_width: float = 3.25,
        sheet_width: float = 2.5,
        coil_width: float = 0.5,
        ambient: float = 0.2,
        diffuse: float = 0.8,
        specular: float = 0.05,
        specular_power: float = 1.0,
        ss_colors: Optional[Dict[str, Tuple[int, int, int]]] = None,
        criteria: str = "",
        criteria_extra: str = "",
    ):
        self.name = name
        self.atom_indices = atom_indices
        self.representation = representation
        self.color_scheme = color_scheme
        self.uniform_color = uniform_color
        self.visible = visible
        self.carbon_color = carbon_color
        self.quality = quality
        self.opacity = opacity
        self.surface_resolution = surface_resolution
        self.surface_radius = surface_radius
        self.atom_scale = atom_scale
        self.bond_radius = bond_radius
        self.ball_scale = ball_scale
        self.stick_radius = stick_radius
        self.backbone_radius = backbone_radius
        self.helix_width = helix_width
        self.sheet_width = sheet_width
        self.coil_width = coil_width
        self.ambient = ambient
        self.diffuse = diffuse
        self.specular = specular
        self.specular_power = specular_power
        self.ss_colors = ss_colors if ss_colors else dict(SS_COLORS)
        self.criteria = criteria
        self.criteria_extra = criteria_extra
        self.actors: list = []


# ---------------------------------------------------------------------------
# Public API class
# ---------------------------------------------------------------------------


class StructureManager:
    """
    Programmatic API for loading, inspecting, editing and saving
    molecular structures.

    Uses MDAnalysis for PDB parsing and atom selections.  The returned
    data are lightweight ``ProteinStructure`` / ``Atom`` objects that
    can be used programmatically or exported to PDB files.

    Examples
    --------
    >>> from gatewizard.core.structure_manager import StructureManager
    >>> sm = StructureManager()
    >>> sm.load_structure("protein.pdb")
    >>> info = sm.get_structure_info()
    >>> print(info['n_atoms'], info['n_residues'])
    """

    def __init__(self):
        self.structure: Optional[ProteinStructure] = None
        self.selections: List[Selection] = []
        self._filepath: Optional[str] = None

    # -- loading -------------------------------------------------------

    def load_structure(self, filepath: Union[str, Path]) -> Dict[str, Any]:
        """Load a PDB file and return summary info.

        Parameters
        ----------
        filepath : str or Path
            Path to a PDB file.

        Returns
        -------
        dict
            Keys: ``n_atoms``, ``n_residues``, ``n_chains``, ``n_bonds``,
            ``chains``, ``title``.
        """
        filepath = str(Path(filepath).resolve())
        if not os.path.isfile(filepath):
            raise StructureError(f"File not found: {filepath}")
        self.structure = parse_pdb(filepath)
        self._filepath = filepath
        self.selections.clear()
        self.selections.append(
            Selection(
                "All",
                list(range(len(self.structure.atoms))),
                representation="vdw",
                color_scheme="element",
            )
        )
        return self.get_structure_info()

    def load_from_pdb_id(
        self, pdb_id: str, output_dir: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """Download a PDB from RCSB and load it.

        Parameters
        ----------
        pdb_id : str
            Four-letter PDB identifier (e.g. ``"1CRN"``).
        output_dir : str or Path, optional
            Directory to save the downloaded file.  Defaults to a temp dir.

        Returns
        -------
        dict
            Same as :meth:`load_structure`.
        """
        import requests as _requests

        pdb_id = pdb_id.strip().upper()
        if not re.match(r"^[0-9A-Z]{4}$", pdb_id):
            raise StructureError(f"Invalid PDB ID: {pdb_id}")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        resp = _requests.get(url, timeout=30)
        resp.raise_for_status()
        if output_dir:
            out = Path(output_dir)
            out.mkdir(parents=True, exist_ok=True)
            pdb_path = out / f"{pdb_id}.pdb"
        else:
            pdb_path = Path(tempfile.mkdtemp()) / f"{pdb_id}.pdb"
        pdb_path.write_text(resp.text)
        return self.load_structure(pdb_path)

    # -- inspection ----------------------------------------------------

    def get_structure_info(self) -> Dict[str, Any]:
        """Return summary info about the loaded structure."""
        self._require_structure()
        s = self.structure
        return {
            "n_atoms": len(s.atoms),
            "n_residues": len(s.residues),
            "n_chains": len(s.chains),
            "n_bonds": len(s.bonds),
            "chains": sorted(s.chains.keys()),
            "title": s.title,
            "source_file": s.source_file,
        }

    def get_chains(self) -> Dict[str, int]:
        """Return chain IDs with their residue counts."""
        self._require_structure()
        return {ch: len(res) for ch, res in self.structure.chains.items()}

    def get_residues(self, chain_id: Optional[str] = None) -> List[Dict[str, Any]]:
        """Return residue info, optionally filtered by chain.

        Parameters
        ----------
        chain_id : str, optional
            Filter to a specific chain.

        Returns
        -------
        list of dict
            Each dict has ``name``, ``seq_id``, ``chain_id``, ``n_atoms``, ``ss``.
        """
        self._require_structure()
        result = []
        for r in self.structure.residues:
            if chain_id and r.chain_id != chain_id:
                continue
            result.append(
                {
                    "name": r.name,
                    "seq_id": r.seq_id,
                    "chain_id": r.chain_id,
                    "n_atoms": len(r.atoms),
                    "ss": r.ss,
                }
            )
        return result

    def get_secondary_structure_summary(self) -> Dict[str, int]:
        """Count residues per secondary structure type."""
        self._require_structure()
        counts: Dict[str, int] = defaultdict(int)
        for r in self.structure.residues:
            counts[r.ss] += 1
        return dict(counts)

    def assign_secondary_structure(self, method: str = "auto") -> Dict[str, int]:
        """Reassign secondary structure using a specific method.

        Parameters
        ----------
        method : str
            Assignment method. One of:

            - ``'auto'`` – PSIQUE → PDB HELIX/SHEET records → heuristic
              (default, same as initial load).
            - ``'psique'`` – Use the PSIQUE program.
            - ``'heuristic'`` – CA-angle heuristic (always available).
            - ``'pdb_records'`` – Only read HELIX/SHEET from the PDB file
              (raises ``StructureError`` if none found).

        Returns
        -------
        dict
            Updated secondary structure summary ``{"H": n, "E": n, ...}``.

        Raises
        ------
        StructureError
            If the requested method is not available or fails.
        """
        self._require_structure()
        method = method.lower()
        if method == "auto":
            _assign_secondary_structure(self.structure, filepath=self._filepath)
        elif method == "psique":
            if not self._filepath:
                raise StructureError("No PDB file path – cannot run PSIQUE")
            ss_map = _assign_ss_psique(self._filepath)
            if ss_map is None:
                raise StructureError(
                    "PSIQUE produced no secondary structure assignments "
                    "for this structure (too few residues?)."
                )
            _apply_ss_map(self.structure, ss_map)
        elif method == "heuristic":
            self.structure.assign_secondary_structure_heuristic()
        elif method == "pdb_records":
            if not self._filepath:
                raise StructureError("No PDB file path available")
            ss_map = _read_ss_from_pdb_records(self._filepath)
            if not ss_map:
                raise StructureError("No HELIX/SHEET records found in PDB file")
            _apply_ss_map(self.structure, ss_map)
        else:
            raise StructureError(
                f"Unknown method '{method}'. "
                "Use 'auto', 'psique', 'heuristic', or 'pdb_records'."
            )
        return self.get_secondary_structure_summary()

    # -- atom selection (MDAnalysis syntax) ----------------------------

    def select_atoms(self, selection_string: str) -> List[int]:
        """Select atom indices using MDAnalysis selection language.

        Parameters
        ----------
        selection_string : str
            MDAnalysis selection expression, e.g. ``"protein"``,
            ``"backbone"``, ``"resname LIG"``, ``"around 5 protein"``.

        Returns
        -------
        list of int
            Atom indices into the current structure.
        """
        self._require_structure()
        import MDAnalysis as mda
        import warnings

        # Build a Universe from current atoms so selections reflect edits
        filepath = self._get_current_pdb_path()
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            u = mda.Universe(filepath)
        ag = u.select_atoms(selection_string)
        # Map MDAnalysis indices to our internal indices
        return [int(i) for i in ag.indices]

    def select_by_criteria(self, criteria: str, extra: str = "") -> List[int]:
        """Convenience method for common selections without MDAnalysis syntax.

        Parameters
        ----------
        criteria : str
            One of ``'All'``, ``'Protein'``, ``'Backbone'``, ``'Sidechain'``,
            ``'Water'``, ``'Ligand'``, ``'Chain'``, ``'Residue range'``.
        extra : str
            Additional info for ``'Chain'`` (chain ID) or
            ``'Residue range'`` (e.g. ``'A:10-50'``).

        Returns
        -------
        list of int
        """
        self._require_structure()
        atoms = self.structure.atoms
        if criteria == "All":
            return list(range(len(atoms)))
        if criteria == "Protein":
            return [i for i, a in enumerate(atoms) if a.res_name in AA_NAMES]
        if criteria == "Backbone":
            return [
                i
                for i, a in enumerate(atoms)
                if a.res_name in AA_NAMES and a.name in BACKBONE_NAMES
            ]
        if criteria == "Sidechain":
            return [
                i
                for i, a in enumerate(atoms)
                if a.res_name in AA_NAMES and a.name not in BACKBONE_NAMES
            ]
        if criteria == "Water":
            return [
                i for i, a in enumerate(atoms) if a.res_name in ("HOH", "WAT", "TIP")
            ]
        if criteria == "Ligand":
            return [
                i
                for i, a in enumerate(atoms)
                if a.res_name not in AA_NAMES
                and a.res_name not in ("HOH", "WAT", "TIP")
            ]
        if criteria in ("Chain", "Chain..."):
            ch = extra.strip().upper()
            return [i for i, a in enumerate(atoms) if a.chain_id == ch]
        if criteria in ("Residue range", "Residue range..."):
            return self._parse_range_text(extra)
        return []

    def auto_detect_molecules(self) -> List[Selection]:
        """Auto-create selections by molecule type (protein, water, ligands).

        Returns
        -------
        list of Selection
        """
        self._require_structure()
        groups: Dict[str, List[int]] = defaultdict(list)
        for i, a in enumerate(self.structure.atoms):
            if a.res_name in AA_NAMES:
                groups["Protein"].append(i)
            elif a.res_name in ("HOH", "WAT", "TIP"):
                groups["Water"].append(i)
            else:
                groups[a.res_name].append(i)
        self.selections.clear()
        color_idx = 0
        for name, indices in groups.items():
            if name == "Protein":
                sel = Selection(
                    name, indices, representation="tube_ss", color_scheme="ss"
                )
            elif name == "Water":
                sel = Selection(
                    name,
                    indices,
                    representation="vdw",
                    color_scheme="element",
                    visible=False,
                )
            else:
                c = CHAIN_PALETTE[color_idx % len(CHAIN_PALETTE)]
                color_idx += 1
                sel = Selection(
                    name,
                    indices,
                    representation="vdw",
                    color_scheme="element",
                    carbon_color=c,
                )
            self.selections.append(sel)
        return list(self.selections)

    # -- editing -------------------------------------------------------

    def rename_chain(self, old_chain: str, new_chain: str) -> int:
        """Rename all atoms/residues in *old_chain* to *new_chain*.

        Parameters
        ----------
        old_chain : str
            Current chain ID (1 character).
        new_chain : str
            New chain ID (1 character).

        Returns
        -------
        int
            Number of atoms renamed.
        """
        self._require_structure()
        if len(new_chain) != 1:
            raise StructureError("Chain ID must be 1 character")
        count = 0
        for atom in self.structure.atoms:
            if atom.chain_id == old_chain:
                atom.chain_id = new_chain
                count += 1
        for res in self.structure.residues:
            if res.chain_id == old_chain:
                res.chain_id = new_chain
        self.structure.chains.clear()
        for res in self.structure.residues:
            self.structure.chains.setdefault(res.chain_id, []).append(res)
        logger.info(f"Renamed chain {old_chain} -> {new_chain} ({count} atoms)")
        return count

    def rename_residues(
        self, chain_id: str, start: int, end: int, new_name: str
    ) -> int:
        """Rename residues in a chain within a range.

        Parameters
        ----------
        chain_id : str
            Chain identifier.
        start : int
            First residue number (inclusive).
        end : int
            Last residue number (inclusive).
        new_name : str
            New residue name (max 3 characters).

        Returns
        -------
        int
            Number of atoms affected.
        """
        self._require_structure()
        new_name = new_name.strip().upper()
        chain_id = chain_id.strip().upper()
        count = 0
        for atom in self.structure.atoms:
            if atom.chain_id == chain_id and start <= atom.res_id <= end:
                atom.res_name = new_name
                count += 1
        for res in self.structure.residues:
            if res.chain_id == chain_id and start <= res.seq_id <= end:
                res.name = new_name
        logger.info(
            f"Renamed residues {chain_id}:{start}-{end} -> {new_name} ({count} atoms)"
        )
        return count

    def renumber_residues(
        self, chain_id: str, start: int, end: int, new_start: int = 1
    ) -> int:
        """Renumber residues in a chain from *new_start*.

        Parameters
        ----------
        chain_id : str
        start : int
            Current first residue number.
        end : int
            Current last residue number.
        new_start : int
            New starting number (default 1).

        Returns
        -------
        int
            Number of atoms renumbered.
        """
        self._require_structure()
        chain_id = chain_id.strip().upper()
        old_ids = sorted(
            set(
                a.res_id
                for a in self.structure.atoms
                if a.chain_id == chain_id and start <= a.res_id <= end
            )
        )
        remap = {old: new_start + i for i, old in enumerate(old_ids)}
        count = 0
        for atom in self.structure.atoms:
            if atom.chain_id == chain_id and atom.res_id in remap:
                atom.res_id = remap[atom.res_id]
                count += 1
        for res in self.structure.residues:
            if res.chain_id == chain_id and res.seq_id in remap:
                res.seq_id = remap[res.seq_id]
        logger.info(
            f"Renumbered {len(remap)} residues in chain {chain_id} ({count} atoms)"
        )
        return count

    def delete_atoms(self, indices: List[int]) -> int:
        """Delete atoms by index and rebuild residues/chains/bonds.

        Parameters
        ----------
        indices : list of int
            Atom indices to remove.

        Returns
        -------
        int
            Number of atoms deleted.
        """
        self._require_structure()
        to_remove = set(indices)
        n_before = len(self.structure.atoms)
        self.structure.atoms = [
            a for i, a in enumerate(self.structure.atoms) if i not in to_remove
        ]
        self.structure._rebuild_residues_and_chains()
        self.structure.build_bonds()
        _assign_secondary_structure(self.structure, filepath=self._filepath)
        # Reset selections to All
        self.selections.clear()
        self.selections.append(
            Selection(
                "All",
                list(range(len(self.structure.atoms))),
                representation="vdw",
                color_scheme="element",
            )
        )
        deleted = n_before - len(self.structure.atoms)
        logger.info(f"Deleted {deleted} atoms")
        return deleted

    def rename_chain_by_indices(self, indices: List[int], new_chain: str) -> int:
        """Rename the chain ID of atoms specified by index list.

        Parameters
        ----------
        indices : list of int
            Atom indices (0-based position in ``structure.atoms``).
        new_chain : str
            New chain ID (1 character).

        Returns
        -------
        int
            Number of atoms renamed.
        """
        self._require_structure()
        new_chain = new_chain.strip().upper()
        if len(new_chain) != 1:
            raise StructureError("Chain ID must be 1 character")
        idx_set = set(indices)
        count = 0
        for i, atom in enumerate(self.structure.atoms):
            if i in idx_set:
                atom.chain_id = new_chain
                count += 1
        # Update residues whose chain_id no longer matches their atoms: use the
        # chain_id of the first atom in the residue as the authoritative value.
        atom_by_pos = {i: a for i, a in enumerate(self.structure.atoms)}
        atom_positions: Dict[Tuple[str, int], List[int]] = {}
        for i, a in enumerate(self.structure.atoms):
            key = (a.chain_id, a.res_id)
            atom_positions.setdefault(key, []).append(i)
        for res in self.structure.residues:
            # Find atoms that belong to this residue by original chain/res_id match
            for i, a in enumerate(self.structure.atoms):
                if a.res_id == res.seq_id and i in idx_set:
                    res.chain_id = new_chain
                    break
        self.structure.chains.clear()
        for res in self.structure.residues:
            self.structure.chains.setdefault(res.chain_id, []).append(res)
        logger.info(f"Renamed chain for {count} atoms (by indices) -> {new_chain}")
        return count

    def rename_residues_by_indices(self, indices: List[int], new_name: str) -> int:
        """Rename the residue name of atoms specified by index list.

        Parameters
        ----------
        indices : list of int
            Atom indices (0-based position in ``structure.atoms``).
        new_name : str
            New residue name (max 4 characters).

        Returns
        -------
        int
            Number of atoms renamed.
        """
        self._require_structure()
        new_name = new_name.strip().upper()
        idx_set = set(indices)
        # Collect (chain_id, res_id) pairs that are in the selection
        sel_keys: set = set()
        count = 0
        for i, atom in enumerate(self.structure.atoms):
            if i in idx_set:
                atom.res_name = new_name
                sel_keys.add((atom.chain_id, atom.res_id))
                count += 1
        for res in self.structure.residues:
            if (res.chain_id, res.seq_id) in sel_keys:
                res.name = new_name
        logger.info(f"Renamed residues for {count} atoms (by indices) -> {new_name}")
        return count

    def renumber_residues_by_indices(
        self, indices: List[int], new_start: int = 1
    ) -> int:
        """Renumber residues that contain atoms in *indices*, starting from *new_start*.

        The unique (chain_id, res_id) pairs found among the selected atoms are
        sorted by res_id and assigned new sequential IDs beginning at
        *new_start*.

        Parameters
        ----------
        indices : list of int
            Atom indices (0-based position in ``structure.atoms``).
        new_start : int
            Starting residue number for the renumbered residues.

        Returns
        -------
        int
            Number of atoms whose res_id was updated.
        """
        self._require_structure()
        idx_set = set(indices)
        pairs = sorted(
            {
                (a.chain_id, a.res_id)
                for i, a in enumerate(self.structure.atoms)
                if i in idx_set
            },
            key=lambda p: (p[0], p[1]),
        )
        remap: Dict[Tuple[str, int], int] = {
            pair: new_start + j for j, pair in enumerate(pairs)
        }
        count = 0
        for atom in self.structure.atoms:
            key = (atom.chain_id, atom.res_id)
            if key in remap:
                atom.res_id = remap[key]
                count += 1
        for res in self.structure.residues:
            key = (res.chain_id, res.seq_id)
            if key in remap:
                res.seq_id = remap[key]
        logger.info(
            f"Renumbered {len(pairs)} residues (by indices), new_start={new_start} ({count} atoms)"
        )
        return count

    def apply_mempro_orientation(self, oriented_pdb: str) -> int:
        """Apply a MemPrO rigid-body orientation to the loaded structure.

        The oriented PDB supplies the target protein pose (plus MemPrO dummy
        atoms).  All atoms in the currently loaded structure — protein, lipids,
        ligands, water, etc. — receive the same transform.

        Parameters
        ----------
        oriented_pdb : str
            Path to a MemPro ``oriented_rank_*.pdb`` file.

        Returns
        -------
        int
            Number of atoms transformed.
        """
        from gatewizard.core.mempro import compute_orientation_transform

        self._require_structure()
        if not self._filepath:
            raise StructureError(
                "No source PDB path available — cannot match atoms for MemPro orientation"
            )
        oriented_path = str(Path(oriented_pdb).resolve())
        if not Path(oriented_path).is_file():
            raise StructureError(f"Oriented PDB not found: {oriented_path}")

        R, t = compute_orientation_transform(self._filepath, oriented_path)
        atoms = self.structure.atoms
        for atom in atoms:
            atom.coord = R @ atom.coord + t
        self.structure.build_bonds()
        self._reassign_ss()
        logger.info(
            f"Applied MemPro orientation from {oriented_path} to {len(atoms)} atoms"
        )
        return len(atoms)

    # -- coordinate transformations ------------------------------------

    def _reassign_ss(self):
        """Re-assign secondary structure after coordinate changes."""
        import tempfile, os

        fd, tmp = tempfile.mkstemp(suffix=".pdb")
        os.close(fd)
        try:
            self.structure.write_pdb(tmp)
            _assign_secondary_structure(self.structure, filepath=tmp)
        finally:
            try:
                os.unlink(tmp)
            except OSError:
                pass

    def rotate_atoms(
        self,
        angle_degrees: float,
        axis: str,
        indices: Optional[List[int]] = None,
        center: str = "selection",
    ) -> int:
        """Rotate atoms around an axis.

        Parameters
        ----------
        angle_degrees : float
            Rotation angle in degrees.
        axis : str
            ``'x'``, ``'y'``, or ``'z'``.
        indices : list of int, optional
            Atom indices to rotate.  ``None`` means all atoms.
        center : str
            ``'selection'`` rotates around the centroid of the affected
            atoms; ``'origin'`` rotates around (0, 0, 0).

        Returns
        -------
        int
            Number of atoms rotated.
        """
        self._require_structure()
        atoms = self.structure.atoms
        if indices is None:
            indices = list(range(len(atoms)))
        if not indices:
            return 0
        coords = np.array([atoms[i].coord for i in indices])
        pivot = coords.mean(axis=0) if center == "selection" else np.zeros(3)
        R = _axis_rotation_matrix(axis, np.radians(angle_degrees))
        for i in indices:
            atoms[i].coord = R @ (atoms[i].coord - pivot) + pivot
        self.structure.build_bonds()
        self._reassign_ss()
        logger.info(f"Rotated {len(indices)} atoms by {angle_degrees}° around {axis}")
        return len(indices)

    def translate_atoms(
        self, displacement: List[float], indices: Optional[List[int]] = None
    ) -> int:
        """Translate atoms by a displacement vector.

        Parameters
        ----------
        displacement : list of float
            ``[dx, dy, dz]`` in angstroms.
        indices : list of int, optional
            Atom indices to translate.  ``None`` means all atoms.

        Returns
        -------
        int
            Number of atoms translated.
        """
        self._require_structure()
        atoms = self.structure.atoms
        if indices is None:
            indices = list(range(len(atoms)))
        if not indices:
            return 0
        d = np.asarray(displacement, dtype=float)
        for i in indices:
            atoms[i].coord = atoms[i].coord + d
        self.structure.build_bonds()
        self._reassign_ss()
        logger.info(f"Translated {len(indices)} atoms by {displacement}")
        return len(indices)

    def center_atoms(self, indices: Optional[List[int]] = None) -> np.ndarray:
        """Move atoms so that their centroid is at the origin.

        Parameters
        ----------
        indices : list of int, optional
            Atom indices whose centroid defines the shift.  ``None`` means
            all atoms.  The shift is always applied to **all** atoms so
            that the structure stays intact.

        Returns
        -------
        numpy.ndarray
            The displacement applied (old centroid position).
        """
        self._require_structure()
        atoms = self.structure.atoms
        ref = indices if indices else list(range(len(atoms)))
        centroid = np.array([atoms[i].coord for i in ref]).mean(axis=0)
        for a in atoms:
            a.coord = a.coord - centroid
        self.structure.build_bonds()
        self._reassign_ss()
        logger.info(f"Centered structure (shift {centroid})")
        return centroid

    def align_to_axis(
        self,
        primary_indices: List[int],
        target_axis: str = "z",
        secondary_indices: Optional[List[int]] = None,
        secondary_axis: Optional[str] = None,
        apply_to: Optional[List[int]] = None,
    ) -> int:
        """Align a selection's principal direction to a reference axis.

        The principal direction is the first singular vector (SVD) fitted
        through the selected atom positions.  If *secondary_indices* and
        *secondary_axis* are given, a secondary rotation around the primary
        axis is applied so that the centroid of the secondary selection
        projects onto the secondary axis in the plane perpendicular to the
        primary.

        Parameters
        ----------
        primary_indices : list of int
            Atoms whose principal direction defines the alignment vector.
        target_axis : str
            ``'x'``, ``'y'``, or ``'z'``.
        secondary_indices : list of int, optional
            Atoms for secondary axis alignment.
        secondary_axis : str, optional
            ``'x'``, ``'y'``, or ``'z'``; must differ from *target_axis*.
        apply_to : list of int, optional
            Atom indices to actually transform.  ``None`` means all atoms.

        Returns
        -------
        int
            Number of atoms transformed.
        """
        self._require_structure()
        atoms = self.structure.atoms
        if apply_to is None:
            apply_to = list(range(len(atoms)))
        if not primary_indices or not apply_to:
            return 0

        axis_idx = {"x": 0, "y": 1, "z": 2}

        # --- compute centroid of atoms being transformed as pivot ---
        pivot = np.array([atoms[i].coord for i in apply_to]).mean(axis=0)

        # --- primary alignment: SVD principal direction → target axis ---
        pri_coords = np.array([atoms[i].coord for i in primary_indices])
        pri_centered = pri_coords - pri_coords.mean(axis=0)
        _, _, Vt = np.linalg.svd(pri_centered, full_matrices=False)
        principal = Vt[0]
        # ensure positive projection on target so direction is consistent
        tidx = axis_idx[target_axis]
        if principal[tidx] < 0:
            principal = -principal
        target_vec = np.zeros(3)
        target_vec[tidx] = 1.0
        R1 = _rotation_matrix_from_vectors(principal, target_vec)

        for i in apply_to:
            atoms[i].coord = R1 @ (atoms[i].coord - pivot) + pivot

        # --- secondary alignment (rotation around the primary axis) ---
        if secondary_indices and secondary_axis:
            sidx = axis_idx[secondary_axis]
            sec_coords = np.array([atoms[i].coord for i in secondary_indices])
            sec_center = sec_coords.mean(axis=0)
            direction = sec_center - pivot
            # project onto plane perpendicular to target axis
            direction[tidx] = 0.0
            norm = np.linalg.norm(direction)
            if norm > 1e-8:
                direction /= norm
                sec_target = np.zeros(3)
                sec_target[sidx] = 1.0
                # angle between projection and secondary target
                cos_a = np.clip(np.dot(direction, sec_target), -1, 1)
                cross = np.cross(direction, sec_target)
                sign = 1.0 if np.dot(cross, target_vec) >= 0 else -1.0
                angle = sign * np.arccos(cos_a)
                R2 = _axis_rotation_matrix(target_axis, angle)
                for i in apply_to:
                    atoms[i].coord = R2 @ (atoms[i].coord - pivot) + pivot

        self.structure.build_bonds()
        self._reassign_ss()
        logger.info(
            f"Aligned {len(primary_indices)} atoms to {target_axis}-axis "
            f"({len(apply_to)} atoms transformed)"
        )
        return len(apply_to)

    # -- saving --------------------------------------------------------

    def save_pdb(self, filepath: Union[str, Path]) -> str:
        """Save the current structure to a PDB file.

        Parameters
        ----------
        filepath : str or Path

        Returns
        -------
        str
            Absolute path to the saved file.
        """
        self._require_structure()
        filepath = Path(filepath).resolve()
        self.structure.write_pdb(filepath)
        logger.info(f"Saved PDB to {filepath}")
        return str(filepath)

    # -- private helpers -----------------------------------------------

    def _require_structure(self):
        if self.structure is None:
            raise StructureError("No structure loaded. Call load_structure() first.")

    def _get_current_pdb_path(self) -> str:
        """Return a file path to the current structure (writing a temp if needed)."""
        if self._filepath and os.path.isfile(self._filepath):
            return self._filepath
        fd, tmp = tempfile.mkstemp(suffix=".pdb")
        os.close(fd)
        self.structure.write_pdb(tmp)
        return tmp

    def _parse_range_text(self, text: str) -> List[int]:
        indices: List[int] = []
        for part in text.split(","):
            part = part.strip()
            if ":" not in part:
                continue
            ch, rng = part.split(":", 1)
            ch = ch.strip().upper()
            try:
                if "-" in rng:
                    lo, hi = int(rng.split("-")[0]), int(rng.split("-")[1])
                else:
                    lo = hi = int(rng.strip())
            except ValueError:
                continue
            for i, a in enumerate(self.structure.atoms):
                if a.chain_id == ch and lo <= a.res_id <= hi:
                    indices.append(i)
        return indices
