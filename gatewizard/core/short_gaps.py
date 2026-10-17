"""Fill protein numbering gaps of 1 to 6 residues with PDBFixer.

Longer gaps are reported and left open. The sequence comes from SEQRES.
A gap with no SEQRES sequence is skipped. A sequence that is not a
standard amino acid is skipped.
"""

from __future__ import annotations

import os
import shutil
import tempfile
from difflib import SequenceMatcher
from typing import Any, Dict, List, Optional, Sequence, Tuple

from gatewizard.core.structure_checker import (
    _load_atoms,
    _numbering_gaps,
    _parse_seqres,
    _residues_by_chain,
)

SHORT_GAP_MAX = 6
# PDBFixer can build these. A gap that names anything else, such as HSK, stays open.
_STANDARD_AMINO_ACIDS = frozenset(
    {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }
)
GapKey = Tuple[str, int, int]


def classify_short_gaps(
    path: str, chains: Optional[Sequence[str]] = None
) -> Dict[str, Any]:
    """Mark each numbering gap as fill (1–6 residues with a SEQRES sequence) or skip."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Structure not found: {path}")
    _universe, rows, _n_frames = _load_atoms(path)
    by_chain = _residues_by_chain(rows)
    allowed = {str(chain).strip() for chain in chains or [] if str(chain).strip()}
    sequences = _sequences_for_gaps(path, by_chain)
    gaps: List[Dict[str, Any]] = []
    for gap in _numbering_gaps(by_chain):
        chain = str(gap["chain"])
        if allowed and chain not in allowed:
            continue
        count = int(gap["count"])
        after_resid = int(gap["after_resid"])
        before_resid = int(gap["before_resid"])
        sequence = sequences.get((chain, after_resid, before_resid)) or []
        if count > SHORT_GAP_MAX:
            action, reason = "skip", "longer than 6 residues"
        elif len(sequence) != count:
            action, reason = "skip", "no SEQRES sequence for this gap"
        elif any(name not in _STANDARD_AMINO_ACIDS for name in sequence):
            action, reason = "skip", "sequence is not a standard amino acid"
        else:
            action, reason = "fill", ""
        gaps.append(
            {
                "chain": chain,
                "after_resid": after_resid,
                "before_resid": before_resid,
                "count": count,
                "sequence": sequence,
                "action": action,
                "reason": reason,
            }
        )
    return {
        "gaps": gaps,
        "filled": [gap for gap in gaps if gap["action"] == "fill"],
        "skipped": [gap for gap in gaps if gap["action"] == "skip"],
    }


def fill_short_gaps(
    input_pdb: str,
    output_pdb: str,
    chains: Optional[Sequence[str]] = None,
    sequences: Optional[Sequence[Mapping[str, Any]]] = None,
) -> Dict[str, Any]:
    """Build every 1–6 residue gap that has a SEQRES sequence.

    The built residue atoms are inserted into a copy of the original PDB.
    Every other line stays as it was. Atom serials after each insertion
    increase by the number of atoms added there, and CONECT serials follow.
    Chains not listed in ``chains`` are dropped from the written model when
    ``chains`` is given. Longer gaps stay missing.
    """
    report = classify_short_gaps(input_pdb, chains)
    _apply_sequence_overrides(report, sequences)
    to_fill = report["filled"]
    if not to_fill:
        if os.path.abspath(output_pdb) != os.path.abspath(input_pdb):
            shutil.copyfile(input_pdb, output_pdb)
        report["output_file"] = output_pdb
        report["built"] = []
        return report
    try:
        from openmm.app import PDBFile
        from pdbfixer import PDBFixer
    except ImportError as exc:
        raise RuntimeError(
            "pdbfixer is not installed. Install it in the GateWizard environment: pip install pdbfixer"
        ) from exc

    fixer = PDBFixer(filename=input_pdb)
    allowed = {str(chain).strip() for chain in chains or [] if str(chain).strip()}
    if allowed:
        drop = [index for index, chain in enumerate(fixer.topology.chains()) if chain.id not in allowed]
        if drop:
            fixer.removeChains(drop)
    fixer.findMissingResidues()
    wanted = {
        (gap["chain"], gap["after_resid"], gap["before_resid"]): list(gap["sequence"])
        for gap in to_fill
    }
    kept: Dict[Tuple[int, int], List[str]] = {}
    built: List[Dict[str, Any]] = []
    chain_list = list(fixer.topology.chains())
    seen: set[Tuple[str, int, int]] = set()
    for (chain_index, residue_index), names in list(fixer.missingResidues.items()):
        if residue_index <= 0 or residue_index >= len(list(chain_list[chain_index].residues())):
            continue
        residues = list(chain_list[chain_index].residues())
        after = residues[residue_index - 1]
        before = residues[residue_index]
        try:
            after_id = int(str(after.id).strip())
            before_id = int(str(before.id).strip())
        except ValueError:
            continue
        key = (chain_list[chain_index].id, after_id, before_id)
        sequence = wanted.get(key)
        if sequence is None:
            continue
        kept[(chain_index, residue_index)] = list(sequence)
        seen.add(key)
        built.append(
            {
                "chain": key[0],
                "after_resid": after_id,
                "before_resid": before_id,
                "count": len(sequence),
                "sequence": list(sequence),
            }
        )
    for gap in to_fill:
        key = (gap["chain"], gap["after_resid"], gap["before_resid"])
        if key in seen:
            continue
        placed = _place_missing_gap(chain_list, gap)
        if placed is None:
            continue
        chain_index, residue_index = placed
        sequence = list(gap["sequence"])
        kept[(chain_index, residue_index)] = sequence
        built.append(
            {
                "chain": gap["chain"],
                "after_resid": gap["after_resid"],
                "before_resid": gap["before_resid"],
                "count": len(sequence),
                "sequence": sequence,
            }
        )
    fixer.missingResidues = kept
    if not kept:
        shutil.copyfile(input_pdb, output_pdb)
        report["output_file"] = output_pdb
        report["built"] = []
        return report
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    try:
        fixer.addMissingHydrogens(7.0)
    except Exception as exc:
        report["hydrogen_warning"] = str(exc)
    fd, built_path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    try:
        with open(built_path, "w", encoding="utf-8") as handle:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        _splice_new_residues(input_pdb, built_path, output_pdb, built, chains)
    finally:
        if os.path.isfile(built_path):
            os.remove(built_path)
    report["output_file"] = os.path.abspath(output_pdb)
    report["built"] = built
    return report


def _atom_identity(line: str) -> tuple[str, int, str] | None:
    if not line.startswith(("ATOM", "HETATM")) or len(line) < 26:
        return None
    chain = line[21:22]
    try:
        resid = int(line[22:26])
    except ValueError:
        return None
    return chain, resid, line[12:16]


def _is_hydrogen(line: str) -> bool:
    element = line[76:78].strip() if len(line) >= 78 else ""
    if element:
        return element.upper() == "H"
    name = line[12:16].strip()
    return name.startswith(("H", "h")) or (bool(name) and name[0].isdigit() and name[1:2] in {"H", "h"})


def _serial_of(line: str) -> int | None:
    if len(line) < 11:
        return None
    text = line[6:11].strip()
    if text.isdigit():
        return int(text)
    return None


def _with_serial(line: str, serial: int) -> str:
    body = line[11:] if len(line) > 11 else ""
    if not body.endswith("\n"):
        body = body.rstrip("\n") + "\n"
    return f"{line[:6]}{serial:5d}{body}"


def _keep_chain_record(line: str, allowed: set[str]) -> bool:
    if not allowed:
        return True
    if line.startswith(("ATOM", "HETATM")):
        return len(line) > 21 and line[21:22] in allowed
    if line.startswith("SEQRES"):
        return len(line) > 11 and line[11:12] in allowed
    if line.startswith("SSBOND"):
        chains = {line[15:16], line[29:30]} if len(line) > 30 else set()
        return bool(chains & allowed) and chains <= allowed
    if line.startswith("LINK"):
        chains = set()
        if len(line) > 21:
            chains.add(line[21:22])
        if len(line) > 51:
            chains.add(line[51:52])
        return bool(chains & allowed) and chains <= allowed
    return True


def _conect_serials(line: str) -> list[int]:
    body = line[6:].strip()
    numbers: list[int] = []
    for start in range(0, len(body), 5):
        chunk = body[start : start + 5].strip()
        if not chunk:
            continue
        if not chunk.isdigit():
            break
        numbers.append(int(chunk))
    return numbers


def _format_conect(serials: list[int]) -> str:
    return "CONECT" + "".join(f"{serial:5d}" for serial in serials) + "\n"


def _splice_new_residues(
    original_pdb: str,
    built_pdb: str,
    output_pdb: str,
    built: Sequence[Mapping[str, Any]],
    chains: Optional[Sequence[str]],
) -> None:
    """Insert heavy atoms of the new residues into a copy of ``original_pdb``."""
    allowed = {str(chain).strip() for chain in chains or [] if str(chain).strip()}
    wanted: Dict[tuple[str, int], None] = {}
    for gap in built:
        chain = str(gap["chain"])
        after = int(gap["after_resid"])
        before = int(gap["before_resid"])
        for resid in range(after + 1, before):
            wanted[(chain, resid)] = None
    new_lines: Dict[tuple[str, int], List[str]] = {}
    with open(built_pdb, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            identity = _atom_identity(line)
            if identity is None or (identity[0], identity[1]) not in wanted:
                continue
            if _is_hydrogen(line):
                continue
            new_lines.setdefault((identity[0], identity[1]), []).append(line if line.endswith("\n") else line + "\n")

    inserts: Dict[int, List[str]] = {}
    with open(original_pdb, encoding="utf-8", errors="replace") as handle:
        source_lines = handle.readlines()
    last_index: Dict[tuple[str, int], int] = {}
    for index, line in enumerate(source_lines):
        if not _keep_chain_record(line, allowed):
            continue
        identity = _atom_identity(line)
        if identity is not None:
            last_index[(identity[0], identity[1])] = index
    for gap in built:
        chain = str(gap["chain"])
        anchor = last_index.get((chain, int(gap["after_resid"])))
        if anchor is None:
            continue
        block: List[str] = []
        for resid in range(int(gap["after_resid"]) + 1, int(gap["before_resid"])):
            block.extend(new_lines.get((chain, resid), []))
        if block:
            inserts.setdefault(anchor, []).extend(block)

    merged: List[str] = []
    inserted_flags: List[bool] = []
    for index, line in enumerate(source_lines):
        if not _keep_chain_record(line, allowed):
            continue
        merged.append(line if line.endswith("\n") else line + "\n")
        inserted_flags.append(False)
        for extra in inserts.get(index, []):
            merged.append(extra)
            inserted_flags.append(True)

    old_to_new: Dict[int, int] = {}
    added_before = 0
    previous_serial = 0
    for index, line in enumerate(merged):
        if not line.startswith(("ATOM", "HETATM", "TER", "ANISOU")):
            continue
        if inserted_flags[index]:
            added_before += 1
            previous_serial += 1
            merged[index] = _with_serial(line, previous_serial)
            continue
        old = _serial_of(line)
        if old is None:
            continue
        new_serial = old + added_before
        previous_serial = new_serial
        old_to_new[old] = new_serial
        if new_serial != old:
            merged[index] = _with_serial(line, new_serial)

    rewritten: List[str] = []
    for line in merged:
        if not line.startswith("CONECT"):
            rewritten.append(line)
            continue
        serials = _conect_serials(line)
        if not serials or serials[0] not in old_to_new:
            continue
        mapped = [old_to_new[serial] for serial in serials if serial in old_to_new]
        if len(mapped) < 2:
            continue
        rewritten.append(_format_conect(mapped))

    with open(output_pdb, "w", encoding="utf-8") as handle:
        handle.writelines(rewritten)


def _apply_sequence_overrides(
    report: Dict[str, Any], sequences: Optional[Sequence[Mapping[str, Any]]]
) -> None:
    """Use a caller sequence for a 1–6 residue gap instead of SEQRES."""
    if not sequences:
        return
    by_gap: Dict[GapKey, List[str]] = {}
    for item in sequences:
        chain = str(item.get("chain") or "").strip()
        try:
            after_resid = int(item.get("after_resid"))
            before_resid = int(item.get("before_resid"))
        except (TypeError, ValueError):
            continue
        names = [str(name).strip().upper() for name in item.get("sequence") or [] if str(name).strip()]
        if chain and names:
            by_gap[(chain, after_resid, before_resid)] = names
    if not by_gap:
        return
    for gap in report["gaps"]:
        key = (str(gap["chain"]), int(gap["after_resid"]), int(gap["before_resid"]))
        names = by_gap.get(key)
        if not names or int(gap["count"]) > SHORT_GAP_MAX:
            continue
        if len(names) != int(gap["count"]):
            gap["action"] = "skip"
            gap["reason"] = "sequence length does not match the gap"
            gap["sequence"] = names
            continue
        gap["sequence"] = names
        gap["action"] = "fill"
        gap["reason"] = ""
    report["filled"] = [gap for gap in report["gaps"] if gap["action"] == "fill"]
    report["skipped"] = [gap for gap in report["gaps"] if gap["action"] == "skip"]


def _place_missing_gap(chain_list: Sequence[Any], gap: Mapping[str, Any]) -> Optional[Tuple[int, int]]:
    """Residue index PDBFixer uses for a hole it did not find in SEQRES."""
    for chain_index, chain in enumerate(chain_list):
        if chain.id != gap["chain"]:
            continue
        residues = list(chain.residues())
        for residue_index, residue in enumerate(residues):
            if residue_index <= 0:
                continue
            try:
                before_id = int(str(residue.id).strip())
                after_id = int(str(residues[residue_index - 1].id).strip())
            except ValueError:
                continue
            if after_id == int(gap["after_resid"]) and before_id == int(gap["before_resid"]):
                return chain_index, residue_index
    return None


def _sequences_for_gaps(path: str, by_chain: Dict[str, List[Dict[str, Any]]]) -> Dict[GapKey, List[str]]:
    seqres = _parse_seqres(path)
    found: Dict[GapKey, List[str]] = {}
    for chain, residues in by_chain.items():
        expected = seqres.get(chain) or []
        if not expected or not residues:
            continue
        built = [str(residue["resname"]) for residue in residues]
        matcher = SequenceMatcher(a=expected, b=built, autojunk=False)
        for tag, i1, i2, j1, j2 in matcher.get_opcodes():
            if tag != "delete" or i2 <= i1 or j1 != j2 or j1 <= 0 or j1 >= len(residues):
                continue
            names = [str(name) for name in expected[i1:i2]]
            after = residues[j1 - 1]
            before = residues[j1]
            if int(before["resid"]) - int(after["resid"]) - 1 != len(names):
                continue
            found[(str(chain), int(after["resid"]), int(before["resid"]))] = names
    return found
