# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""Report missing atoms, gaps, alternate conformations, and chain differences.

The likelihood of an alternate conformation is the PDB occupancy, the column
just before the B-factor. The altLoc letter is column 17. The most likely
letter is the highest occupancy. A tie keeps the earlier letter.
"""

from __future__ import annotations

import os
import shutil
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np

from gatewizard.utils.helpers import resolve_pdb_chain_id
from gatewizard.utils.peptide_residues import PEPTIDE_POLYMER_RESIDUES
from gatewizard.utils.trajectory_tools import _AMBER_LIPID_RESNAMES

_HEAVY: Dict[str, frozenset] = {
    "ALA": frozenset({"N", "CA", "C", "O", "CB"}),
    "ARG": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"}),
    "ASN": frozenset({"N", "CA", "C", "O", "CB", "CG", "OD1", "ND2"}),
    "ASP": frozenset({"N", "CA", "C", "O", "CB", "CG", "OD1", "OD2"}),
    "CYS": frozenset({"N", "CA", "C", "O", "CB", "SG"}),
    "GLN": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "NE2"}),
    "GLU": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD", "OE1", "OE2"}),
    "GLY": frozenset({"N", "CA", "C", "O"}),
    "HIS": frozenset({"N", "CA", "C", "O", "CB", "CG", "ND1", "CD2", "CE1", "NE2"}),
    "ILE": frozenset({"N", "CA", "C", "O", "CB", "CG1", "CG2", "CD1"}),
    "LEU": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"}),
    "LYS": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"}),
    "MET": frozenset({"N", "CA", "C", "O", "CB", "CG", "SD", "CE"}),
    "PHE": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"}),
    "PRO": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD"}),
    "SER": frozenset({"N", "CA", "C", "O", "CB", "OG"}),
    "THR": frozenset({"N", "CA", "C", "O", "CB", "OG1", "CG2"}),
    "TRP": frozenset(
        {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"}
    ),
    "TYR": frozenset({"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"}),
    "VAL": frozenset({"N", "CA", "C", "O", "CB", "CG1", "CG2"}),
    "MSE": frozenset({"N", "CA", "C", "O", "CB", "CG", "SE", "CE"}),
    "SEC": frozenset({"N", "CA", "C", "O", "CB", "SE"}),
}

_HEAVY_PARENT: Dict[str, str] = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "HSD": "HIS",
    "HSE": "HIS",
    "HSP": "HIS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
    "CYX": "CYS",
    "CYM": "CYS",
    "TYM": "TYR",
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

_WATERS = frozenset({"HOH", "WAT", "DOD", "TIP", "TIP3", "TIP4", "SOL", "H2O"})
_IONS = frozenset(
    {
        "NA",
        "CL",
        "K",
        "MG",
        "ZN",
        "CA",
        "FE",
        "MN",
        "CU",
        "NI",
        "CO",
        "ION",
        "NA+",
        "CL-",
        "BR",
        "IOD",
        "CS",
        "RB",
        "LI",
        "CD",
        "HG",
    }
)
_NUCLEOTIDES = frozenset(
    {
        "A",
        "C",
        "G",
        "U",
        "I",
        "DA",
        "DC",
        "DG",
        "DT",
        "DI",
        "ADE",
        "CYT",
        "GUA",
        "THY",
        "URA",
        "URI",
        "RA",
        "RC",
        "RG",
        "RU",
    }
)
_KNOWN = frozenset(PEPTIDE_POLYMER_RESIDUES) | _WATERS | _IONS | _NUCLEOTIDES | frozenset(
    _AMBER_LIPID_RESNAMES
)


def _is_ion_resname(name: str) -> bool:
    """True for K, K+, Cl-, NA+, and the same ion with a charge written on the name."""
    token = str(name or "").upper().strip()
    if token in _IONS:
        return True
    stem = token.rstrip("+-")
    while stem and stem[-1].isdigit():
        stem = stem[:-1]
    return bool(stem) and stem in _IONS


def _is_known_resname(name: str) -> bool:
    token = str(name or "").upper().strip()
    return token in _KNOWN or _is_ion_resname(token)
_BACKBONE_BREAK_A = 2.0
_CHAIN_IDENTITY = 0.8
# Two copies occupy one site when some heavy atoms sit on top of each other
# and at least half of the smaller copy is within a short contact of the other.
_OVERLAP_CONTACT_A = 1.5
_OVERLAP_NEAR_A = 2.5
_OVERLAP_FRACTION = 0.5
_DUPLICATE_POSE_A = 0.05
_MAX_CROSS_COPIES = 40


def _heavy_names(resname: str) -> Optional[frozenset]:
    name = resname.upper()
    parent = _HEAVY_PARENT.get(name, name)
    return _HEAVY.get(parent)


def _is_hydrogen(name: str, element: str) -> bool:
    if element == "H":
        return True
    if element:
        return False
    return name.startswith("H") and name not in {"HG"}


def _selection(chain: str, resid: int, atom: str = "") -> str:
    parts = []
    if chain:
        parts.append(f"chainID {chain}")
    parts.append(f"resid {int(resid)}")
    if atom:
        parts.append(f"name {atom}")
    return " and ".join(parts)


def _altloc_selection(resname: str, chain: str, resid: int) -> str:
    """Highlight one residue, including a ligand that reuses a protein number."""
    parts = []
    if resname:
        parts.append(f"resname {resname}")
    if chain:
        parts.append(f"chainID {chain}")
    parts.append(f"resid {int(resid)}")
    return " and ".join(parts)


def _or_selection(chain: str, resids: Sequence[int]) -> str:
    body = " or ".join(f"resid {int(resid)}" for resid in resids)
    if chain:
        return f"chainID {chain} and ({body})"
    return body


def _count_models(path: str) -> int:
    count = 0
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("MODEL"):
                    count += 1
                elif line.startswith(("data_", "loop_", "_atom_site.")):
                    return 0
    except OSError:
        return 0
    return count


def _parse_seqres(path: str) -> Dict[str, List[str]]:
    chains: Dict[str, List[str]] = {}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("SEQRES") or len(line) < 20:
                    continue
                chain = line[11:12].strip() or "_"
                names = [part.upper() for part in line[19:].split() if part.strip()]
                chains.setdefault(chain, []).extend(names)
    except OSError:
        return {}
    return chains


def _atom_field(atom: Any, name: str, default: str = "") -> str:
    try:
        value = getattr(atom, name)
    except Exception:
        return default
    if value is None:
        return default
    return str(value).strip()


def _load_atoms(path: str, topology: str | None = None) -> Tuple[Any, List[Dict[str, Any]], int]:
    import MDAnalysis as mda

    top = str(topology or "").strip()
    if top and os.path.isfile(top):
        universe = mda.Universe(top, path)
    else:
        universe = mda.Universe(path)
    n_frames = int(getattr(universe.trajectory, "n_frames", 1) or 1)
    if n_frames > 1:
        universe.trajectory[0]
    rows: List[Dict[str, Any]] = []
    for index, atom in enumerate(universe.atoms):
        # Same letter the viewer stores. A prmtop segment named SYSTEM with an
        # empty chain column becomes S, which is what the 3D view selects.
        chain = resolve_pdb_chain_id(_atom_field(atom, "segid"), _atom_field(atom, "chainID"))
        name = _atom_field(atom, "name").upper()
        element = _atom_field(atom, "element").upper()
        if element in {"", "DUMMY"}:
            element = ""
        try:
            occupancy = float(atom.occupancy)
        except Exception:
            occupancy = 1.0
        try:
            bfactor = float(atom.tempfactor)
        except Exception:
            bfactor = 0.0
        try:
            resid = int(atom.resid)
        except Exception:
            resid = 0
        icode = _atom_field(atom, "icode")
        altloc = _atom_field(atom, "altLoc").upper()
        resname = _atom_field(atom, "resname").upper()
        try:
            coord = np.asarray(atom.position, dtype=float)
        except Exception:
            coord = np.zeros(3)
        rows.append(
            {
                "index": index,
                "chain": chain,
                "resid": resid,
                "icode": icode,
                "resname": resname,
                "name": name,
                "element": element,
                "altloc": altloc,
                "occupancy": occupancy,
                "bfactor": bfactor,
                "coord": coord,
            }
        )
    return universe, rows, n_frames


def _most_likely(options: Sequence[Mapping[str, Any]]) -> str:
    ordered = sorted(
        options,
        key=lambda item: (-float(item["occupancy"]), str(item["altloc"])),
    )
    return str(ordered[0]["altloc"]) if ordered else ""


def _residues_by_chain(rows: Sequence[Mapping[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    grouped: Dict[Tuple[str, int, str, str], Dict[str, Any]] = {}
    order: List[Tuple[str, int, str, str]] = []
    for row in rows:
        if row["resname"] not in PEPTIDE_POLYMER_RESIDUES:
            continue
        key = (str(row["chain"]), int(row["resid"]), str(row["icode"]), str(row["resname"]))
        slot = grouped.get(key)
        if slot is None:
            slot = {
                "chain": key[0],
                "resid": key[1],
                "icode": key[2],
                "resname": key[3],
                "names": set(),
                "c": None,
                "n": None,
            }
            grouped[key] = slot
            order.append(key)
        if not _is_hydrogen(str(row["name"]), str(row["element"])):
            slot["names"].add(str(row["name"]))
        if row["name"] == "C" and (not row["altloc"] or slot["c"] is None):
            slot["c"] = row["coord"]
        if row["name"] == "N" and (not row["altloc"] or slot["n"] is None):
            slot["n"] = row["coord"]
    chains: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for key in order:
        chains[key[0]].append(grouped[key])
    for chain in chains:
        chains[chain].sort(key=lambda item: (int(item["resid"]), str(item["icode"])))
    return chains


def _missing_atoms(chains: Mapping[str, Sequence[Mapping[str, Any]]]) -> List[Dict[str, Any]]:
    missing: List[Dict[str, Any]] = []
    for residues in chains.values():
        for residue in residues:
            expected = _heavy_names(str(residue["resname"]))
            if not expected:
                continue
            absent = sorted(expected - set(residue["names"]))
            for atom in absent:
                missing.append(
                    {
                        "chain": residue["chain"],
                        "resid": residue["resid"],
                        "icode": residue["icode"],
                        "resname": residue["resname"],
                        "atom": atom,
                        "selection": _selection(str(residue["chain"]), int(residue["resid"]), atom),
                    }
                )
    return missing


def _numbering_gaps(chains: Mapping[str, Sequence[Mapping[str, Any]]]) -> List[Dict[str, Any]]:
    gaps: List[Dict[str, Any]] = []
    for residues in chains.values():
        for prev, nxt in zip(residues, residues[1:]):
            if prev["icode"] or nxt["icode"]:
                continue
            span = int(nxt["resid"]) - int(prev["resid"])
            if span <= 1:
                continue
            missing_ids = list(range(int(prev["resid"]) + 1, int(nxt["resid"])))
            gaps.append(
                {
                    "chain": prev["chain"],
                    "after_resid": prev["resid"],
                    "before_resid": nxt["resid"],
                    "count": len(missing_ids),
                    "missing": missing_ids[:30],
                    "selection": _or_selection(str(prev["chain"]), [int(prev["resid"]), int(nxt["resid"])]),
                }
            )
    return gaps


def _seqres_gaps(path: str, chains: Mapping[str, Sequence[Mapping[str, Any]]]) -> List[Dict[str, Any]]:
    seqres = _parse_seqres(path)
    found: List[Dict[str, Any]] = []
    for chain, expected in seqres.items():
        built = [str(residue["resname"]) for residue in chains.get(chain, [])]
        if not expected or expected == built:
            continue
        matcher = SequenceMatcher(a=expected, b=built, autojunk=False)
        for tag, i1, i2, _j1, _j2 in matcher.get_opcodes():
            if tag not in {"delete", "replace"} or i2 <= i1:
                continue
            names = expected[i1:i2]
            found.append(
                {
                    "chain": "" if chain == "_" else chain,
                    "source": "seqres",
                    "count": len(names),
                    "resnames": names[:30],
                    "selection": _selection(chain, int(chains[chain][0]["resid"])) if chains.get(chain) else "",
                }
            )
    return found


def _backbone_breaks(chains: Mapping[str, Sequence[Mapping[str, Any]]]) -> List[Dict[str, Any]]:
    breaks: List[Dict[str, Any]] = []
    for residues in chains.values():
        for prev, nxt in zip(residues, residues[1:]):
            if prev["icode"] or nxt["icode"]:
                continue
            if int(nxt["resid"]) - int(prev["resid"]) != 1:
                continue
            if prev["c"] is None or nxt["n"] is None:
                continue
            distance = float(np.linalg.norm(np.asarray(nxt["n"]) - np.asarray(prev["c"])))
            if distance <= _BACKBONE_BREAK_A:
                continue
            breaks.append(
                {
                    "chain": prev["chain"],
                    "resid_a": prev["resid"],
                    "resname_a": prev["resname"],
                    "resid_b": nxt["resid"],
                    "resname_b": nxt["resname"],
                    "distance": round(distance, 2),
                    "selection": _or_selection(str(prev["chain"]), [int(prev["resid"]), int(nxt["resid"])]),
                }
            )
    return breaks


def _altloc_groups(rows: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """One entry per residue that has two or more altLoc letters."""
    grouped: Dict[Tuple[str, int, str, str], List[Mapping[str, Any]]] = defaultdict(list)
    for row in rows:
        if not row["altloc"]:
            continue
        key = (str(row["chain"]), int(row["resid"]), str(row["icode"]), str(row["resname"]))
        grouped[key].append(row)
    groups: List[Dict[str, Any]] = []
    for key, members in grouped.items():
        by_letter: Dict[str, List[Mapping[str, Any]]] = defaultdict(list)
        for member in members:
            by_letter[str(member["altloc"])].append(member)
        if len(by_letter) < 2:
            continue
        options = []
        for letter, atoms in sorted(by_letter.items()):
            counts: Dict[float, int] = defaultdict(int)
            for atom in atoms:
                counts[round(float(atom["occupancy"]), 2)] += 1
            occupancy = max(counts.items(), key=lambda item: (item[1], item[0]))[0]
            options.append(
                {
                    "altloc": letter,
                    "occupancy": occupancy,
                    "atoms": len(atoms),
                }
            )
        default = _most_likely(options)
        chain, resid, icode, resname = key
        groups.append(
            {
                "id": f"{chain}|{resid}|{icode}|{resname}",
                "chain": chain,
                "resid": resid,
                "icode": icode,
                "resname": resname,
                "default": default,
                "options": [{**option, "likely": option["altloc"] == default} for option in options],
                "selection": _altloc_selection(resname, chain, resid),
            }
        )
    groups.sort(key=lambda item: (str(item["chain"]), int(item["resid"]), str(item["icode"]), str(item["resname"])))
    return groups


def _row_coord(row: Mapping[str, Any]) -> Optional[np.ndarray]:
    coord = row.get("coord")
    if coord is None:
        return None
    arr = np.asarray(coord, dtype=float)
    if arr.shape != (3,) or not np.isfinite(arr).all():
        return None
    return arr


def _heavy_rows(rows: Sequence[Mapping[str, Any]]) -> List[Mapping[str, Any]]:
    heavy: List[Mapping[str, Any]] = []
    for row in rows:
        if _row_coord(row) is None:
            continue
        if _is_hydrogen(str(row["name"]), str(row.get("element") or "")):
            continue
        heavy.append(row)
    return heavy


def _poses_overlap(left: Sequence[Mapping[str, Any]], right: Sequence[Mapping[str, Any]]) -> bool:
    """True when two atom sets occupy one site and are not the same coordinates twice.

    A real alternate conformation puts part of one pose on top of the other.
    An exact second copy of the same coordinates is a duplicate record.
    """
    heavy_left = _heavy_rows(left)
    heavy_right = _heavy_rows(right)
    if not heavy_left or not heavy_right:
        return False
    xyz_left = np.stack([_row_coord(row) for row in heavy_left])
    xyz_right = np.stack([_row_coord(row) for row in heavy_right])
    distances = np.linalg.norm(xyz_left[:, None, :] - xyz_right[None, :, :], axis=2)
    if float(distances.min()) > _OVERLAP_CONTACT_A:
        return False
    if len(heavy_left) <= len(heavy_right):
        near = distances.min(axis=1)
    else:
        near = distances.min(axis=0)
    if float(np.mean(near <= _OVERLAP_NEAR_A)) < _OVERLAP_FRACTION:
        return False
    by_name = {str(row["name"]): _row_coord(row) for row in heavy_left}
    paired = [
        float(np.linalg.norm(by_name[str(row["name"])] - _row_coord(row)))
        for row in heavy_right
        if str(row["name"]) in by_name
    ]
    if paired and float(np.median(paired)) < _DUPLICATE_POSE_A:
        return False
    return True


def _occurrence_blocks(members: Sequence[Mapping[str, Any]]) -> List[List[Mapping[str, Any]]]:
    """Split one residue into poses when an atom name is stored more than once."""
    counts: Dict[str, int] = defaultdict(int)
    for row in members:
        counts[str(row["name"])] += 1
    if not counts or max(counts.values()) < 2:
        return []
    seen: Dict[str, int] = defaultdict(int)
    blocks: Dict[int, List[Mapping[str, Any]]] = defaultdict(list)
    for row in members:
        name = str(row["name"])
        seen[name] += 1
        blocks[seen[name]].append(row)
    return [blocks[index] for index in sorted(blocks)]


def _cluster_overlapping(blocks: Sequence[Sequence[Mapping[str, Any]]]) -> List[List[int]]:
    parent = list(range(len(blocks)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for left in range(len(blocks)):
        for right in range(left + 1, len(blocks)):
            if _poses_overlap(blocks[left], blocks[right]):
                parent[find(left)] = find(right)
    clusters: Dict[int, List[int]] = defaultdict(list)
    for index in range(len(blocks)):
        clusters[find(index)].append(index)
    return [members for members in clusters.values() if len(members) >= 2]


def _option_from_rows(letter: str, rows: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    counts: Dict[float, int] = defaultdict(int)
    for row in rows:
        counts[round(float(row.get("occupancy") or 0.0), 2)] += 1
    occupancy = max(counts.items(), key=lambda item: (item[1], item[0]))[0] if counts else 1.0
    return {"altloc": letter, "occupancy": occupancy, "atoms": len(rows)}


def _cross_residue_eligible(resname: str) -> bool:
    name = resname.upper()
    if name in PEPTIDE_POLYMER_RESIDUES or name in _WATERS or _is_ion_resname(name) or name in _NUCLEOTIDES:
        return False
    if name in _AMBER_LIPID_RESNAMES:
        return False
    return True


def _copies_selection(options: Sequence[Mapping[str, Any]]) -> str:
    parts = [
        f"({_altloc_selection(str(option['resname']), str(option['chain']), int(option['resid']))})"
        for option in options
    ]
    return " or ".join(parts)


def _inferred_overlap_groups(
    rows: Sequence[Mapping[str, Any]],
    labeled: Iterable[Tuple[str, int, str, str]],
) -> List[Dict[str, Any]]:
    """Alternate conformations that have no altLoc letter.

    One residue with a repeated atom name is included when those copies
    overlap. Separate residues of the same ligand name are included when
    their heavy atoms occupy one site. Protein, water, ions, nucleotides,
    and Amber lipids are not paired across residue numbers.
    """
    labeled_keys = set(labeled)
    grouped: Dict[Tuple[str, int, str, str], List[Mapping[str, Any]]] = defaultdict(list)
    order: List[Tuple[str, int, str, str]] = []
    for row in rows:
        key = (str(row["chain"]), int(row["resid"]), str(row["icode"]), str(row["resname"]))
        if key not in grouped:
            order.append(key)
        grouped[key].append(row)
    groups: List[Dict[str, Any]] = []
    for key in order:
        if key in labeled_keys:
            continue
        members = grouped[key]
        letters = {str(row["altloc"]) for row in members if row["altloc"]}
        if len(letters) >= 2:
            continue
        blocks = _occurrence_blocks(members)
        for cluster in _cluster_overlapping(blocks):
            chosen = [blocks[index] for index in cluster]
            options = []
            for offset, block in enumerate(chosen):
                letter = chr(ord("A") + offset)
                option = _option_from_rows(letter, block)
                option.update(
                    {
                        "chain": key[0],
                        "resid": key[1],
                        "icode": key[2],
                        "resname": key[3],
                    }
                )
                options.append(option)
            default = _most_likely(options)
            chain, resid, icode, resname = key
            groups.append(
                {
                    "id": f"{chain}|{resid}|{icode}|{resname}",
                    "chain": chain,
                    "resid": resid,
                    "icode": icode,
                    "resname": resname,
                    "default": default,
                    "inferred": True,
                    "mode": "occurrence",
                    "options": [{**option, "likely": option["altloc"] == default} for option in options],
                    "selection": _altloc_selection(resname, chain, resid),
                }
            )
    by_name: Dict[str, List[Tuple[str, int, str, str]]] = defaultdict(list)
    for key in order:
        if _cross_residue_eligible(key[3]):
            by_name[key[3]].append(key)
    for resname, keys in by_name.items():
        if len(keys) < 2 or len(keys) > _MAX_CROSS_COPIES:
            continue
        blocks = [grouped[key] for key in keys]
        for cluster in _cluster_overlapping(blocks):
            chosen_keys = [keys[index] for index in cluster]
            options = []
            for offset, key in enumerate(chosen_keys):
                letter = chr(ord("A") + offset)
                option = _option_from_rows(letter, grouped[key])
                option.update({"chain": key[0], "resid": key[1], "icode": key[2], "resname": key[3]})
                options.append(option)
            default = _most_likely(options)
            first = next(option for option in options if option["altloc"] == default)
            ident = "|".join(f"{option['chain']}{option['resid']}{option['icode']}" for option in options)
            groups.append(
                {
                    "id": f"overlap|{resname}|{ident}",
                    "chain": first["chain"],
                    "resid": first["resid"],
                    "icode": first["icode"],
                    "resname": resname,
                    "default": default,
                    "inferred": True,
                    "mode": "copies",
                    "label": " and ".join(
                        f"{option['resname']} {option['chain']}{option['resid']}{option['icode']} ({option['altloc']})"
                        for option in options
                    ),
                    "options": [{**option, "likely": option["altloc"] == default} for option in options],
                    "selection": _copies_selection(options),
                }
            )
    return groups


def _alternate_groups(rows: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """AltLoc letters, plus overlapping atoms that were stored without a letter."""
    labeled = _altloc_groups(rows)
    inferred = _inferred_overlap_groups(rows, (_residue_key(group) for group in labeled))
    groups = labeled + inferred
    groups.sort(
        key=lambda item: (
            str(item["chain"]),
            int(item["resid"]),
            str(item["icode"]),
            str(item["resname"]),
            str(item["id"]),
        )
    )
    return groups


def _chain_pairs(chains: Mapping[str, Sequence[Mapping[str, Any]]]) -> List[Dict[str, Any]]:
    names = [chain for chain, residues in chains.items() if residues]
    pairs: List[Dict[str, Any]] = []
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            seq_a = list(chains[left])
            seq_b = list(chains[right])
            names_a = [str(item["resname"]) for item in seq_a]
            names_b = [str(item["resname"]) for item in seq_b]
            longest = max(len(names_a), len(names_b), 1)
            length_close = abs(len(names_a) - len(names_b)) <= max(2, round(0.1 * longest))
            matcher = SequenceMatcher(a=names_a, b=names_b, autojunk=False)
            matches = sum(i2 - i1 for tag, i1, i2, _j1, _j2 in matcher.get_opcodes() if tag == "equal")
            identity = matches / longest
            if identity < _CHAIN_IDENTITY and not length_close:
                continue
            substitutions = []
            length_notes = []
            for tag, i1, i2, j1, j2 in matcher.get_opcodes():
                if tag == "equal":
                    continue
                if tag == "replace" and (i2 - i1) == (j2 - j1):
                    for offset in range(i2 - i1):
                        residue_a = seq_a[i1 + offset]
                        residue_b = seq_b[j1 + offset]
                        substitutions.append(
                            {
                                "resname_a": residue_a["resname"],
                                "resid_a": residue_a["resid"],
                                "resname_b": residue_b["resname"],
                                "resid_b": residue_b["resid"],
                                "selection": _selection(left, int(residue_a["resid"])),
                            }
                        )
                    continue
                length_notes.append(
                    {
                        "chain_a": left,
                        "from_resid_a": seq_a[i1]["resid"] if i1 < len(seq_a) else None,
                        "to_resid_a": seq_a[i2 - 1]["resid"] if i1 < i2 else None,
                        "chain_b": right,
                        "from_resid_b": seq_b[j1]["resid"] if j1 < len(seq_b) else None,
                        "to_resid_b": seq_b[j2 - 1]["resid"] if j1 < j2 else None,
                        "count_a": i2 - i1,
                        "count_b": j2 - j1,
                    }
                )
            if len(names_a) == len(names_b):
                length_text = "Same length."
            elif len(names_a) > len(names_b):
                length_text = f"Chain {left} is {len(names_a) - len(names_b)} residues longer."
            else:
                length_text = f"Chain {right} is {len(names_b) - len(names_a)} residues longer."
            pairs.append(
                {
                    "chain_a": left,
                    "chain_b": right,
                    "length_a": len(names_a),
                    "length_b": len(names_b),
                    "identity": round(identity, 3),
                    "same": not substitutions and not length_notes and len(names_a) == len(names_b),
                    "length_note": length_text,
                    "substitutions": substitutions,
                    "length_differences": length_notes,
                    "selection": _selection(left, int(seq_a[0]["resid"])) if seq_a else "",
                }
            )
    return pairs


def _other_issues(
    path: str,
    rows: Sequence[Mapping[str, Any]],
    chains: Mapping[str, Sequence[Mapping[str, Any]]],
    n_frames: int,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    issues: List[Dict[str, Any]] = []
    model_count = max(_count_models(path), n_frames if n_frames > 1 else 0)
    if model_count > 1:
        issues.append(
            {
                "kind": "models",
                "message": f"{model_count} models are in the file. This check uses the first one.",
                "selection": "all",
            }
        )
    seen: Dict[Tuple[Any, ...], int] = {}
    hydrogen_count = 0
    heavy_count = 0
    unknown: Dict[str, Dict[str, Any]] = {}
    for row in rows:
        protein = row["resname"] in PEPTIDE_POLYMER_RESIDUES
        if protein and _is_hydrogen(str(row["name"]), str(row["element"])):
            hydrogen_count += 1
        elif protein:
            heavy_count += 1
        key = (row["chain"], row["resid"], row["icode"], row["name"], row["altloc"])
        seen[key] = seen.get(key, 0) + 1
        if row["altloc"] == "" and float(row["occupancy"]) <= 0.0:
            issues.append(
                {
                    "kind": "occupancy",
                    "message": (
                        f"{row['resname']} {row['chain']}{row['resid']} {row['name']} has occupancy 0."
                    ),
                    "selection": _selection(str(row["chain"]), int(row["resid"]), str(row["name"])),
                }
            )
        if row["resname"] and not _is_known_resname(str(row["resname"])):
            slot = unknown.setdefault(
                str(row["resname"]),
                {"count": 0, "chains": set()},
            )
            slot["count"] += 1
            if row["chain"]:
                slot["chains"].add(str(row["chain"]))
    for key, count in seen.items():
        if count < 2:
            continue
        chain, resid, _icode, atom, altloc = key
        letter = f" altLoc {altloc}" if altloc else ""
        issues.append(
            {
                "kind": "duplicate",
                "message": f"{chain}{resid} {atom}{letter} is recorded {count} times.",
                "selection": _selection(str(chain), int(resid), str(atom)),
            }
        )
    for resname, slot in sorted(unknown.items()):
        chains_text = ", ".join(sorted(slot["chains"])) or "no chain id"
        issues.append(
            {
                "kind": "unknown",
                "message": f"{resname} is not a standard residue ({slot['count']} atoms, {chains_text}).",
                "selection": f"resname {resname}",
            }
        )
    if heavy_count and hydrogen_count == 0:
        issues.append(
            {
                "kind": "hydrogens",
                "message": "No hydrogens are present on the protein. That is normal for many PDB files.",
                "selection": "protein",
            }
        )
    elif heavy_count and hydrogen_count:
        residues_with_h = {
            (row["chain"], row["resid"], row["icode"])
            for row in rows
            if row["resname"] in PEPTIDE_POLYMER_RESIDUES
            and _is_hydrogen(str(row["name"]), str(row["element"]))
        }
        protein_residues = sum(len(residues) for residues in chains.values())
        if residues_with_h and len(residues_with_h) < protein_residues:
            issues.append(
                {
                    "kind": "hydrogens",
                    "message": (
                        f"Hydrogens are present on {len(residues_with_h)} of {protein_residues} residues."
                    ),
                    "selection": "byres element H",
                }
            )
    hydrogens = {
        "protein_heavy": heavy_count,
        "protein_hydrogen": hydrogen_count,
    }
    return issues, hydrogens


def diagnose_structure(path: str, topology: str | None = None) -> Dict[str, Any]:
    """Report missing atoms, gaps, alternate conformations, and chain differences.

    PDB and mmCIF are read with MDAnalysis. GRO files still report missing
    atoms and gaps. They have no altLoc. Only the first model is inspected.
    A residue with no altLoc letter is still reported when two copies of its
    heavy atoms occupy the same site.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Structure not found: {path}")
    _universe, rows, n_frames = _load_atoms(path, topology)
    chains = _residues_by_chain(rows)
    missing_atoms = _missing_atoms(chains)
    numbering = _numbering_gaps(chains)
    seqres = _seqres_gaps(path, chains)
    missing_residues = numbering + seqres
    breaks = _backbone_breaks(chains)
    altlocs = _alternate_groups(_altloc_rows(path, rows))
    pairs = _chain_pairs(chains)
    issues, hydrogens = _other_issues(path, rows, chains, n_frames)
    residue_count = sum(item.get("count", 1) for item in missing_residues)
    return {
        "path": os.path.abspath(path),
        "summary": {
            "missing_atoms": len(missing_atoms),
            "missing_residues": residue_count,
            "gaps": len(breaks),
            "altlocs": len(altlocs),
            "chain_pairs": len(pairs),
            "chain_differences": sum(1 for pair in pairs if not pair["same"]),
            "other": len(issues),
        },
        "hydrogens": hydrogens,
        "missing_atoms": missing_atoms,
        "missing_residues": missing_residues,
        "gaps": breaks,
        "altlocs": altlocs,
        "chain_pairs": pairs,
        "issues": issues,
    }


def _altloc_rows(path: str, fallback: Sequence[Mapping[str, Any]]) -> List[Mapping[str, Any]]:
    """AltLoc rows from PDB ``ATOM`` and ``HETATM`` lines.

    A PDB line is the source for ligands and other hetero groups. MDAnalysis
    sometimes leaves their altLoc empty. GRO and mmCIF have no PDB coordinate
    lines, so those files keep the MDAnalysis rows. Only the first model is read.
    """
    rows: List[Dict[str, Any]] = []
    saw_coordinate = False
    try:
        handle = open(path, "r", encoding="utf-8", errors="replace")
    except OSError:
        return list(fallback)
    with handle:
        for line in handle:
            if line.startswith("MODEL") and saw_coordinate:
                break
            if line.startswith("ENDMDL"):
                break
            parsed = _pdb_coordinate(line)
            if parsed is None:
                continue
            saw_coordinate = True
            rows.append(parsed)
    if not saw_coordinate:
        return list(fallback)
    return rows


def _pdb_coordinate(line: str) -> Optional[Dict[str, Any]]:
    """One ``ATOM`` or ``HETATM`` line, or None for every other record."""
    if len(line) < 27 or not line.startswith(("ATOM", "HETATM")):
        return None
    try:
        resid = int(line[22:26])
    except ValueError:
        return None
    try:
        occupancy = float(line[54:60]) if len(line) >= 60 else 1.0
    except ValueError:
        occupancy = 1.0
    coord = None
    if len(line) >= 54:
        try:
            coord = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
        except ValueError:
            coord = None
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    return {
        "chain": line[21].strip(),
        "resid": resid,
        "icode": line[26].strip(),
        "resname": line[17:20].strip().upper(),
        "name": line[12:16].strip().upper(),
        "element": element,
        "altloc": line[16].strip().upper(),
        "occupancy": occupancy,
        "coord": coord,
    }


def _pdb_altloc_fields(line: str) -> Optional[Tuple[str, int, str, str, str, str]]:
    """Chain, residue number, insertion code, residue name, atom name, and altLoc.

    Returned only when the line is ``ATOM``, ``HETATM``, or ``ANISOU`` and the
    altLoc column is set. Other records stay in the file untouched.
    """
    if len(line) < 27 or not line.startswith(("ATOM", "HETATM", "ANISOU")):
        return None
    letter = line[16].strip().upper()
    if not letter:
        return None
    try:
        resid = int(line[22:26])
    except ValueError:
        return None
    return (
        line[21].strip(),
        resid,
        line[26].strip(),
        line[17:20].strip().upper(),
        line[12:16].strip().upper(),
        letter,
    )


def _without_altloc(line: str) -> str:
    """Blank the altLoc column and set occupancy to 1.00."""
    if len(line) < 60 or not line.startswith(("ATOM", "HETATM")):
        if line.startswith("ANISOU") and len(line) > 16:
            return line[:16] + " " + line[17:]
        return line
    return line[:16] + " " + line[17:54] + "  1.00" + line[60:]


def _residue_key(group: Mapping[str, Any]) -> Tuple[str, int, str, str]:
    return (
        str(group["chain"]),
        int(group["resid"]),
        str(group["icode"]),
        str(group["resname"]).upper(),
    )


def _choice_targets(
    choices: Optional[Iterable[Mapping[str, Any]]],
    groups: Sequence[Mapping[str, Any]],
) -> Tuple[Dict[Tuple[str, int, str, str, str], str], Dict[Tuple[str, int, str, str], str]]:
    """Atom choices win over a residue-wide letter.

    A choice with ``resname`` applies to that residue. A choice without one
    applies only when a single residue has that chain, number, and insertion
    code, so a ligand that reuses a protein number is left alone.
    """
    places: Dict[Tuple[str, int, str], List[Tuple[str, int, str, str]]] = defaultdict(list)
    for group in groups:
        key = _residue_key(group)
        places[key[:3]].append(key)
    by_atom: Dict[Tuple[str, int, str, str, str], str] = {}
    by_residue: Dict[Tuple[str, int, str, str], str] = {}
    for choice in choices or []:
        altloc = str(choice.get("altloc") or "").strip().upper()
        if not altloc:
            continue
        try:
            resid = int(choice.get("resid"))
        except (TypeError, ValueError):
            continue
        place = (
            str(choice.get("chain") or "").strip(),
            resid,
            str(choice.get("icode") or "").strip(),
        )
        resname = str(choice.get("resname") or "").strip().upper()
        if resname:
            keys = [(*place, resname)]
        else:
            keys = places.get(place) or []
            if len(keys) != 1:
                continue
        atom = str(choice.get("atom") or "").strip().upper()
        for key in keys:
            if atom:
                by_atom[(*key, atom)] = altloc
            else:
                by_residue[key] = altloc
    return by_atom, by_residue


def _pdb_atom_identity(line: str) -> Optional[Dict[str, Any]]:
    """Chain, residue, atom name, and altLoc for ``ATOM``, ``HETATM``, and ``ANISOU``."""
    if len(line) < 27 or not line.startswith(("ATOM", "HETATM", "ANISOU")):
        return None
    try:
        resid = int(line[22:26])
    except ValueError:
        return None
    return {
        "record": line[:6].strip(),
        "chain": line[21].strip(),
        "resid": resid,
        "icode": line[26].strip(),
        "resname": line[17:20].strip().upper(),
        "name": line[12:16].strip().upper(),
        "altloc": line[16].strip().upper(),
    }


def _conformer_plan(
    groups: Sequence[Mapping[str, Any]],
    choices: Optional[Iterable[Mapping[str, Any]]],
) -> Dict[str, Any]:
    """Letters to keep for altLoc rows, repeated atom names, and overlapping copies."""
    by_atom, by_residue = _choice_targets(choices, groups)
    altloc: Dict[Tuple[str, int, str, str], str] = {}
    occurrence: Dict[Tuple[str, int, str, str], str] = {}
    drop: set[Tuple[str, int, str, str]] = set()
    for group in groups:
        key = _residue_key(group)
        letter = str(by_residue.get(key) or group["default"]).upper()
        mode = str(group.get("mode") or "")
        if mode == "copies":
            for option in group["options"]:
                copy = (
                    str(option["chain"]),
                    int(option["resid"]),
                    str(option["icode"]),
                    str(option["resname"]).upper(),
                )
                if str(option["altloc"]).upper() != letter:
                    drop.add(copy)
            continue
        if mode == "occurrence":
            occurrence[key] = letter
            continue
        altloc[key] = letter
    return {"altloc": altloc, "by_atom": by_atom, "occurrence": occurrence, "drop": drop}


def _filter_structure_lines(
    lines: Iterable[str],
    plan: Mapping[str, Any],
) -> Tuple[List[str], int, int]:
    """Drop conformations that were not chosen. SEQRES, LINK, and CONECT stay."""
    raw = [line if line.endswith("\n") else line + "\n" for line in lines]
    occurrence: Mapping[Tuple[str, int, str, str], str] = plan["occurrence"]
    totals: Dict[Tuple[Tuple[str, int, str, str], str], int] = defaultdict(int)
    identities: List[Optional[Dict[str, Any]]] = []
    for line in raw:
        ident = _pdb_atom_identity(line)
        identities.append(ident)
        if ident is None or ident["record"] not in {"ATOM", "HETATM"}:
            continue
        residue = (ident["chain"], ident["resid"], ident["icode"], ident["resname"])
        if residue in occurrence:
            totals[(residue, ident["name"])] += 1
    seen: Dict[Tuple[Tuple[str, int, str, str], str], int] = defaultdict(int)
    last_kept: Dict[Tuple[Tuple[str, int, str, str], str], bool] = {}
    kept_lines: List[str] = []
    kept_atoms = 0
    removed = 0
    altloc: Mapping[Tuple[str, int, str, str], str] = plan["altloc"]
    by_atom: Mapping[Tuple[str, int, str, str, str], str] = plan["by_atom"]
    drop: set[Tuple[str, int, str, str]] = plan["drop"]
    for line, ident in zip(raw, identities):
        if ident is None:
            kept_lines.append(line)
            continue
        residue = (ident["chain"], ident["resid"], ident["icode"], ident["resname"])
        coordinate = ident["record"] in {"ATOM", "HETATM"}
        if residue in drop:
            if coordinate:
                removed += 1
            continue
        if residue in occurrence:
            atom_key = (residue, ident["name"])
            if ident["record"] == "ANISOU":
                if last_kept.get(atom_key, True):
                    kept_lines.append(line)
                continue
            seen[atom_key] += 1
            shared = totals[atom_key] == 1
            letter = chr(ord("A") + seen[atom_key] - 1)
            keep = shared or letter == occurrence[residue]
            last_kept[atom_key] = keep
            if not keep:
                removed += 1
                continue
            kept_atoms += 1
            kept_lines.append(_without_altloc(line))
            continue
        if not ident["altloc"]:
            if coordinate:
                kept_atoms += 1
            kept_lines.append(line)
            continue
        wanted = by_atom.get((*residue, ident["name"])) or altloc.get(residue)
        if wanted and ident["altloc"] != str(wanted).upper():
            if coordinate:
                removed += 1
            continue
        if coordinate:
            kept_atoms += 1
        kept_lines.append(_without_altloc(line) if coordinate else line[:16] + " " + line[17:])
    return kept_lines, kept_atoms, removed


def single_conformer_pdb_lines(lines: Sequence[str]) -> Tuple[List[str], int]:
    """Keep one pose of a ligand that has altLoc letters or overlapping copies.

    The kept pose is the higher occupancy. A tie keeps the earlier letter.
    A ligand with one pose is returned unchanged. The count is how many poses
    were present.
    """
    rows = []
    for line in lines:
        parsed = _pdb_coordinate(line)
        if parsed is not None:
            rows.append(parsed)
    groups = _alternate_groups(rows)
    poses = max((len(group["options"]) for group in groups), default=1)
    if not groups:
        return list(lines), 1
    kept, _atoms, _removed = _filter_structure_lines(lines, _conformer_plan(groups, None))
    return kept, poses


def _write_filtered_pdb(path: str, destination: str, lines: Sequence[str]) -> None:
    Path(destination).parent.mkdir(parents=True, exist_ok=True)
    with open(destination, "w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(lines)


def apply_altloc_choices(
    path: str,
    choices: Optional[Sequence[Mapping[str, Any]]] = None,
    output_path: Optional[str] = None,
) -> Dict[str, Any]:
    """Write a PDB that keeps one altLoc for each disordered atom.

    A choice is ``{chain, resid, icode, resname, altloc}`` for the whole
    residue, or the same fields plus ``atom`` for one atom. ``resname``
    distinguishes a ligand from a protein residue that uses the same number.
    Without ``resname``, the letter is used only when one residue has that
    chain, number, and insertion code. A residue without a choice keeps the
    highest occupancy. ``ATOM`` and ``HETATM`` stay in their original records.
    SEQRES, LINK, and CONECT stay in the file. The original file is not changed.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Structure not found: {path}")
    _universe, rows, _n_frames = _load_atoms(path)
    groups = _alternate_groups(_altloc_rows(path, rows))
    plan = _conformer_plan(groups, choices)
    destination = os.path.abspath(
        output_path or str(Path(path).with_name(f"{Path(path).stem}_checked.pdb"))
    )
    if not groups and not plan["by_atom"]:
        if os.path.abspath(path) != destination:
            Path(destination).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, destination)
        return {
            "output_file": destination,
            "atoms_kept": len(rows),
            "atoms_removed": 0,
        }
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        source_lines = handle.readlines()
    lines, kept_atoms, removed = _filter_structure_lines(source_lines, plan)
    if kept_atoms <= 0:
        raise ValueError("No atoms left after applying alternate conformations.")
    _write_filtered_pdb(path, destination, lines)
    return {
        "output_file": destination,
        "atoms_kept": kept_atoms,
        "atoms_removed": removed,
    }


def keep_most_likely_altlocs(path: str, output_path: str) -> Dict[str, Any]:
    """Keep the higher-occupancy altLoc of each residue.

    Protein, ligands, waters, and other hetero groups are separate residues.
    A file with no alternate conformations is copied unchanged, including
    SEQRES and ANISOU. A tie keeps the earlier letter, the same rule as
    :func:`apply_altloc_choices` with no explicit choices. Other records
    stay in the file.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Structure not found: {path}")
    _universe, rows, _n_frames = _load_atoms(path)
    groups = _alternate_groups(_altloc_rows(path, rows))
    destination = os.path.abspath(output_path)
    if not groups:
        if os.path.abspath(path) != destination:
            Path(destination).parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, destination)
        return {
            "output_file": destination,
            "atoms_kept": len(rows),
            "atoms_removed": 0,
            "residues": 0,
        }
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        source_lines = handle.readlines()
    lines, kept_atoms, removed = _filter_structure_lines(source_lines, _conformer_plan(groups, None))
    _write_filtered_pdb(path, destination, lines)
    return {
        "output_file": destination,
        "atoms_kept": kept_atoms,
        "atoms_removed": removed,
        "residues": len(groups),
    }
