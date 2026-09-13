# gatewizard/utils/peptide_residues.py
# SPDX-License-Identifier: MIT
"""Shared L/D / peptide-cap residue recognition for prep, visualize, builder, eq, analysis.

Standard MDA ``protein`` / AA_NAMES miss D-amino acids and formyl/ethanolamine
caps (e.g. gramicidin 1jno: FVA, DLE, DVA, ETA). This module is the single
source of truth for treating those residues as part of a peptide polymer
rather than free GAFF ligands.
"""

from __future__ import annotations

from typing import Dict, FrozenSet, Iterable, Optional

# ---------------------------------------------------------------------------
# Standard L-amino acids (+ common protonation / MSE variants)
# ---------------------------------------------------------------------------

STANDARD_L_AMINO_ACIDS: FrozenSet[str] = frozenset(
    {
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
        # Protonation / Amber variants
        "HID",
        "HIE",
        "HIP",
        "HSD",
        "HSE",
        "HSP",
        "ASH",
        "GLH",
        "LYN",
        "CYX",
        "CYM",
        "TYM",
        "ARN",
        "SEP",
        "T2P",
        "COO",
    }
)

# Common PDB CCD / Amber-style D-amino acid residue names
D_AMINO_ACIDS: FrozenSet[str] = frozenset(
    {
        "DAL",  # D-Ala
        "DAR",  # D-Arg
        "DAS",  # D-Asp
        "DCY",  # D-Cys
        "DGL",  # D-Glu / D-Gln variously
        "DGN",  # D-Gln
        "DHI",  # D-His
        "DIL",  # D-Ile
        "DLE",  # D-Leu (gramicidin)
        "DLY",  # D-Lys
        "DPN",  # D-Phe
        "DPR",  # D-Pro
        "DSN",  # D-Ser
        "DTH",  # D-Thr
        "DTR",  # D-Trp
        "DTY",  # D-Tyr
        "DVA",  # D-Val (gramicidin)
        "MED",  # N-methyl-D-something (keep as polymer when linked)
    }
)

# Terminal / peptide caps (standard Amber + native peptide termini)
PEPTIDE_CAPS: FrozenSet[str] = frozenset(
    {
        "ACE",
        "NME",
        "NMA",
        "NHE",
        "FVA",  # N-formyl-L-valine (gramicidin)
        "FOR",  # formyl
        "ETA",  # ethanolamine C-terminus (gramicidin)
        "NH2",
    }
)

# Caps that mean "already capped — do not add ACE/NME"
NATIVE_PEPTIDE_CAPS: FrozenSet[str] = frozenset(
    {"ACE", "NME", "NMA", "NHE", "FVA", "FOR", "ETA"}
)

# Markers that classify the whole polymer as Peptide (not Protein).
# ACE/NME alone stay Protein; formyl / ethanolamine / D-aa flip to Peptide.
PEPTIDE_CLASS_MARKERS: FrozenSet[str] = frozenset(D_AMINO_ACIDS | {"FVA", "FOR", "ETA", "NH2"})

# Full polymer set used for Visualize / Prep / restraints / strip-H
PEPTIDE_POLYMER_RESIDUES: FrozenSet[str] = frozenset(
    STANDARD_L_AMINO_ACIDS | D_AMINO_ACIDS | PEPTIDE_CAPS
)

# Amber / tleap-oriented renames for common CCD peptide names.
# Identity entries document names that are already Amber-friendly.
# FVA/ETA often need a supporting leaprc / frcmod; we keep CCD names but
# classify them as polymer (not GAFF ligands). Callers may warn.
AMBER_PEPTIDE_RESNAME_MAP: Dict[str, str] = {
    "DLE": "DLE",
    "DVA": "DVA",
    "DAL": "DAL",
    "DAR": "DAR",
    "DAS": "DAS",
    "DPN": "DPN",
    "DPR": "DPR",
    "DSN": "DSN",
    "DTH": "DTH",
    "DTR": "DTR",
    "DTY": "DTY",
    "DLY": "DLY",
    "DIL": "DIL",
    "DHI": "DHI",
    "DGL": "DGL",
    "DGN": "DGN",
    "DCY": "DCY",
    "FVA": "FVA",
    "FOR": "FOR",
    "ETA": "ETA",
    "ACE": "ACE",
    "NME": "NME",
    "NMA": "NME",
    "NHE": "NHE",
}

# Residue names that commonly lack stock ff19SB parameters (warn in Builder).
AMBER_UNSUPPORTED_PEPTIDE_WARN: FrozenSet[str] = frozenset({"FVA", "FOR", "ETA"})

# CCD D-amino acid → stock L Amber residue names for tleap loadPDB.
# Chirality comes from coordinates; do not use for PropKa rename.
D_TO_L_AMBER_MAP: Dict[str, str] = {
    "DAL": "ALA",
    "DAR": "ARG",
    "DAS": "ASP",
    "DCY": "CYS",
    "DGL": "GLU",
    "DGN": "GLN",
    "DHI": "HIS",
    "DIL": "ILE",
    "DLE": "LEU",
    "DLY": "LYS",
    "DPN": "PHE",
    "DPR": "PRO",
    "DSN": "SER",
    "DTH": "THR",
    "DTR": "TRP",
    "DTY": "TYR",
    "DVA": "VAL",
}

# Polymer-cap roles for GAFF residue libraries (head/tail connect atoms).
PEPTIDE_CAP_CONNECT: Dict[str, Dict[str, Optional[str]]] = {
    "FVA": {"role": "n_term", "head_atom": None, "tail_atom": "C"},
    "FOR": {"role": "n_term", "head_atom": None, "tail_atom": "C"},
    "ETA": {"role": "c_term", "head_atom": "N", "tail_atom": None},
}


def normalize_resname(res_name: Optional[str]) -> str:
    return (res_name or "").strip().upper()


def is_peptide_polymer_residue(res_name: Optional[str]) -> bool:
    """True if *res_name* should be treated as part of a peptide/protein polymer."""
    return normalize_resname(res_name) in PEPTIDE_POLYMER_RESIDUES


def is_native_peptide_cap(res_name: Optional[str]) -> bool:
    return normalize_resname(res_name) in NATIVE_PEPTIDE_CAPS


def is_d_amino_acid(res_name: Optional[str]) -> bool:
    return normalize_resname(res_name) in D_AMINO_ACIDS


def is_peptide_class_marker(res_name: Optional[str]) -> bool:
    """True if *res_name* marks a polymer as Peptide (vs standard Protein)."""
    return normalize_resname(res_name) in PEPTIDE_CLASS_MARKERS


def classify_polymer_kind(res_names: Iterable[Optional[str]]) -> str:
    """Return ``'peptide'`` or ``'protein'`` for a set of polymer residue names.

    Any D-aa / formyl / ethanolamine marker → peptide; otherwise protein.
    """
    for name in res_names:
        if is_peptide_class_marker(name):
            return "peptide"
    return "protein"


def amber_resname_for(res_name: Optional[str]) -> str:
    """Map a PDB CCD peptide name to an Amber-oriented residue name."""
    key = normalize_resname(res_name)
    return AMBER_PEPTIDE_RESNAME_MAP.get(key, key)


def peptide_polymer_resnames_sorted() -> list[str]:
    return sorted(PEPTIDE_POLYMER_RESIDUES)


def mda_peptide_or_protein_selection() -> str:
    """MDAnalysis selection covering standard protein + peptide polymer residues.

    Prefer this over bare ``protein`` when D-aa / formyl / ETA may be present.
    """
    names = " ".join(peptide_polymer_resnames_sorted())
    return f"(protein or resname {names})"


def mda_peptide_polymer_only_selection() -> str:
    """Selection by residue name only (no MDA protein keyword)."""
    names = " ".join(peptide_polymer_resnames_sorted())
    return f"resname {names}"


def default_peptide_exclude_selection() -> str:
    """Default APL / EVAPL exclude set for peptide-in-membrane systems."""
    return mda_peptide_or_protein_selection()


def polymer_resnames_in_pdb_lines(lines: Iterable[str]) -> FrozenSet[str]:
    found: set[str] = set()
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
            name = normalize_resname(line[17:20])
            if name in PEPTIDE_POLYMER_RESIDUES:
                found.add(name)
    return frozenset(found)


def detect_peptide_caps_in_pdb(pdb_file: str) -> list[str]:
    """Return sorted unique native/standard peptide cap names in a PDB."""
    found: set[str] = set()
    try:
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                    name = normalize_resname(line[17:20])
                    if name in NATIVE_PEPTIDE_CAPS:
                        found.add(name)
    except OSError:
        return []
    return sorted(found)


def amber_unsupported_peptide_names_in_pdb(pdb_file: str) -> list[str]:
    """Return sorted peptide polymer names that likely need extra Amber libs."""
    found: set[str] = set()
    try:
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                    name = normalize_resname(line[17:20])
                    if name in AMBER_UNSUPPORTED_PEPTIDE_WARN:
                        found.add(name)
    except OSError:
        return []
    return sorted(found)


def d_amino_acid_names_in_pdb(pdb_file: str) -> list[str]:
    """Return sorted unique D-amino acid residue names present in a PDB."""
    found: set[str] = set()
    try:
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                    name = normalize_resname(line[17:20])
                    if name in D_AMINO_ACIDS and name in D_TO_L_AMBER_MAP:
                        found.add(name)
    except OSError:
        return []
    return sorted(found)


def pdb_has_peptide_hetatm_markers(pdb_file: str) -> bool:
    """True if PDB has D-aa or GAFF-needed caps as polymer markers (keepligs hint)."""
    markers = D_AMINO_ACIDS | AMBER_UNSUPPORTED_PEPTIDE_WARN
    try:
        with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                    if normalize_resname(line[17:20]) in markers:
                        return True
    except OSError:
        return False
    return False


def detect_peptide_caps_needing_gaff(pdb_file: str) -> list[dict]:
    """Return structured GAFF-needed polymer caps: ``[{name, role}, ...]``."""
    present = set(amber_unsupported_peptide_names_in_pdb(pdb_file))
    out: list[dict] = []
    for name in sorted(present):
        meta = PEPTIDE_CAP_CONNECT.get(name, {})
        out.append(
            {
                "name": name,
                "role": meta.get("role") or "unknown",
                "head_atom": meta.get("head_atom"),
                "tail_atom": meta.get("tail_atom"),
            }
        )
    return out


def remap_d_amino_acids_for_tleap(input_pdb: str, output_pdb: str) -> dict:
    """Remap CCD D-amino acid names to L Amber partners for tleap loadPDB.

    Keeps FVA/FOR/ETA unchanged. Returns the same stats dict as
    :func:`remap_pdb_peptide_resnames`.
    """
    return remap_pdb_peptide_resnames(
        input_pdb, output_pdb, mapping=D_TO_L_AMBER_MAP
    )


def remap_pdb_peptide_resnames(
    input_pdb: str,
    output_pdb: str,
    *,
    mapping: Optional[Dict[str, str]] = None,
) -> dict:
    """Rewrite ATOM/HETATM residue names using *mapping* (default Amber map).

    Returns counts of renamed records / unique residue keys.
    """
    map_table = mapping if mapping is not None else AMBER_PEPTIDE_RESNAME_MAP
    renamed_records = 0
    renamed_keys: set[str] = set()
    out_lines: list[str] = []
    with open(input_pdb, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                old = normalize_resname(line[17:20])
                new = map_table.get(old, old)
                if new != old:
                    line = line[:17] + f"{new:<3}" + line[20:]
                    renamed_records += 1
                    chain = line[21:22].strip() if len(line) > 21 else ""
                    try:
                        resid = int(line[22:26].strip())
                    except ValueError:
                        resid = 0
                    renamed_keys.add(f"{old}->{new}:{chain}{resid}")
            out_lines.append(line)
    with open(output_pdb, "w", encoding="utf-8") as handle:
        handle.writelines(out_lines)
    return {
        "record_changes": renamed_records,
        "residue_changes": len(renamed_keys),
        "renamed": sorted(renamed_keys),
    }
