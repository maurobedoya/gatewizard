"""Restore original PDB residue numbers after pdb4amber / capping."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

from gatewizard.utils.logger import get_logger
from gatewizard.utils.peptide_residues import is_peptide_polymer_residue
from gatewizard.utils.residue_mapping import (
    ResidueChainKey,
    align_mapping_keys_to_pdb,
    compose_amber_to_original,
    parse_residue_mapping_by_chain,
    parse_residue_mapping_rows,
)

logger = get_logger(__name__)

CAP_RESNAMES = frozenset({"ACE", "NME", "NMA"})

PdbPath = Union[str, Path]
ResidueSnapshot = Tuple[str, int, str, str]  # chain, resid, icode, resname
# (amber_chain, amber_resid) → (original_chain, original_resid); None = new cap
RestoreTarget = Optional[Tuple[str, int]]


def find_pdb4amber_renum_path(output_pdb: PdbPath) -> Optional[Path]:
    """Return ``{stem}_renum.txt`` next to a pdb4amber output, if it exists."""
    p = Path(output_pdb)
    candidates = [
        p.with_name(f"{p.stem}_renum.txt"),
        p.with_name(f"{p.stem}renum.txt"),
    ]
    for cand in candidates:
        if cand.is_file():
            return cand
    return None


def discover_cap_mapping_path(input_pdb: PdbPath) -> Optional[Path]:
    """Sibling ``*_gatewizard_residue_mapping.txt`` written by ProteinCapper."""
    p = Path(input_pdb)
    names = [
        f"{p.stem}_gatewizard_residue_mapping.txt",
        f"{p.stem}_gatewizard_to_original.txt",
    ]
    if p.stem.endswith("_capped"):
        base = p.stem[: -len("_capped")]
        names.append(f"{base}_capped_gatewizard_residue_mapping.txt")
        names.append(f"{base}_gatewizard_residue_mapping.txt")
    for name in names:
        cand = p.parent / name
        if cand.is_file():
            return cand
    return None


def discover_original_pdb(input_pdb: PdbPath) -> Optional[Path]:
    """Uncapped sibling when *input_pdb* is ``*_capped.pdb``."""
    p = Path(input_pdb)
    if not p.stem.endswith("_capped"):
        return None
    sibling = p.with_name(f"{p.stem[: -len('_capped')]}.pdb")
    if sibling.is_file() and sibling.resolve() != p.resolve():
        return sibling
    return None


def snapshot_pdb_residues(pdb_path: PdbPath) -> List[ResidueSnapshot]:
    """Unique residues in file order: ``(chain, resid, icode, resname)``."""
    seen: List[ResidueSnapshot] = []
    last: Optional[ResidueSnapshot] = None
    path = Path(pdb_path)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            rec = _parse_pdb_atom(line)
            if rec is None:
                continue
            key = (rec["chain"], rec["resid"], rec["icode"], rec["resname"])
            if key != last:
                seen.append(key)
                last = key
    return seen


def restore_original_residue_numbers(
    pdb_path: PdbPath,
    remum_path: Optional[PdbPath] = None,
    output_pdb: Optional[PdbPath] = None,
    cap_mapping_path: Optional[PdbPath] = None,
    original_pdb: Optional[PdbPath] = None,
    assign_new_caps: bool = True,
    pre_amber_snapshot: Optional[Sequence[ResidueSnapshot]] = None,
) -> Dict[str, Any]:
    """
    Rewrite ATOM/HETATM resSeq (and chain) back to pre-pipeline numbers.

    When *assign_new_caps* is True, ACE/NME do not reuse protein ids:
    ACE is N-terminus − 1 (0 if the chain starts at 1), NME is C-terminus + 1,
    skipping ids already used on that chain.

    Returns a dict with ``residue_numbers_preserved``, ``output_file``,
    ``remum_path``, ``cap_mapping_path``, and ``cap_assignments``.
    """
    src = Path(pdb_path)
    if not src.is_file():
        raise FileNotFoundError(f"PDB file not found: {src}")

    remum_file = Path(remum_path) if remum_path else None
    cap_file = Path(cap_mapping_path) if cap_mapping_path else None
    orig_file = Path(original_pdb) if original_pdb else None

    remum_map = (
        parse_residue_mapping_by_chain(remum_file) if remum_file and remum_file.is_file() else {}
    )
    cap_map = (
        parse_residue_mapping_by_chain(cap_file) if cap_file and cap_file.is_file() else {}
    )
    pdb_keys = [
        (chain, resid) for chain, resid, _icode, _name in snapshot_pdb_residues(src)
    ]
    remum_map = align_mapping_keys_to_pdb(remum_map, pdb_keys)

    amber_to_original = _build_restore_map(
        remum_map=remum_map,
        cap_map=cap_map,
        original_pdb=orig_file,
        remum_path=remum_file,
        pre_amber_snapshot=pre_amber_snapshot,
        output_pdb=src,
    )

    out_path = Path(output_pdb) if output_pdb else src
    rewritten, cap_assignments = _rewrite_pdb_resids(
        src,
        out_path,
        amber_to_original,
        assign_new_caps=assign_new_caps,
        original_pdb=orig_file,
    )

    logger.info(
        "Restored original residue numbers in %s (%s atom records; caps=%s)",
        out_path,
        rewritten,
        cap_assignments,
    )
    return {
        "residue_numbers_preserved": True,
        "output_file": str(out_path),
        "remum_path": str(remum_file) if remum_file and remum_file.is_file() else None,
        "cap_mapping_path": str(cap_file) if cap_file and cap_file.is_file() else None,
        "cap_assignments": {
            f"{chain}:{resname}": resid for (chain, resname), resid in cap_assignments.items()
        },
        "atoms_rewritten": rewritten,
    }


def _build_restore_map(
    remum_map: Mapping[ResidueChainKey, Optional[int]],
    cap_map: Mapping[ResidueChainKey, Optional[int]],
    original_pdb: Optional[Path],
    remum_path: Optional[Path],
    pre_amber_snapshot: Optional[Sequence[ResidueSnapshot]],
    output_pdb: Path,
) -> Dict[ResidueChainKey, RestoreTarget]:
    original_res: List[ResidueSnapshot] = []
    if original_pdb and original_pdb.is_file():
        original_res = [
            snap for snap in snapshot_pdb_residues(original_pdb) if snap[3] not in CAP_RESNAMES
        ]
    elif pre_amber_snapshot:
        original_res = [
            snap for snap in pre_amber_snapshot if snap[3] not in CAP_RESNAMES
        ]
    out_res = [
        snap for snap in snapshot_pdb_residues(output_pdb) if snap[3] not in CAP_RESNAMES
    ]
    snapshot_map = _pair_snapshots_to_targets(out_res, original_res) if original_res else {}
    if snapshot_map and _restore_coverage(snapshot_map, out_res) >= 0.9:
        return snapshot_map

    if remum_map and cap_map:
        return _resid_map_to_targets(compose_amber_to_original(remum_map, cap_map))

    if remum_map and original_pdb and original_pdb.is_file():
        paired = _pair_remum_to_original_pdb(remum_path, original_pdb)
        if paired and _restore_coverage(_resid_map_to_targets(paired), out_res) >= 0.9:
            return _resid_map_to_targets(paired)

    if remum_map:
        return _resid_map_to_targets(remum_map)

    if snapshot_map:
        return snapshot_map

    logger.warning(
        "No remum / cap mapping / original PDB to restore residue numbers from %s",
        output_pdb,
    )
    return {}


def _resid_map_to_targets(
    resid_map: Mapping[ResidueChainKey, Optional[int]],
) -> Dict[ResidueChainKey, RestoreTarget]:
    mapped: Dict[ResidueChainKey, RestoreTarget] = {}
    for (chain, final_id), orig_id in resid_map.items():
        if orig_id is None:
            mapped[(chain, final_id)] = None
        else:
            mapped[(chain, final_id)] = (chain, orig_id)
    return mapped


def _restore_coverage(
    mapped: Mapping[ResidueChainKey, RestoreTarget],
    output_res: Sequence[ResidueSnapshot],
) -> float:
    if not output_res:
        return 0.0
    hits = sum(1 for chain, resid, _icode, _name in output_res if (chain, resid) in mapped)
    return hits / len(output_res)


def _split_polymer_hetero(
    snaps: Sequence[ResidueSnapshot],
) -> Tuple[List[ResidueSnapshot], List[ResidueSnapshot]]:
    polymer: List[ResidueSnapshot] = []
    hetero: List[ResidueSnapshot] = []
    for snap in snaps:
        if _is_polymer_non_cap(snap[3]):
            polymer.append(snap)
        else:
            hetero.append(snap)
    return polymer, hetero


def _zip_snapshot_targets(
    source: Sequence[ResidueSnapshot],
    dest: Sequence[ResidueSnapshot],
) -> Dict[ResidueChainKey, RestoreTarget]:
    mapped: Dict[ResidueChainKey, RestoreTarget] = {}
    for src_snap, dst_snap in zip(source, dest):
        mapped[(src_snap[0], src_snap[1])] = (dst_snap[0], dst_snap[1])
    return mapped


def _pair_snapshots_to_targets(
    output_res: Sequence[ResidueSnapshot],
    original_res: Sequence[ResidueSnapshot],
) -> Dict[ResidueChainKey, RestoreTarget]:
    """
    Map Amber output residues back to the original PDB.

    pdb4amber / tleap often merge chains onto a blank chain and sequential
    ids. When per-chain counts no longer match, pair polymer then hetero
    in file order so A/B 1–261 is restored from a 1–522 blank chain.
    """
    if not output_res or not original_res:
        return {}

    out_by_chain: Dict[str, List[ResidueSnapshot]] = {}
    orig_by_chain: Dict[str, List[ResidueSnapshot]] = {}
    for snap in output_res:
        out_by_chain.setdefault(snap[0], []).append(snap)
    for snap in original_res:
        orig_by_chain.setdefault(snap[0], []).append(snap)

    per_chain_ok = bool(out_by_chain) and all(
        chain in orig_by_chain and len(orig_by_chain[chain]) == len(snaps)
        for chain, snaps in out_by_chain.items()
    )
    if per_chain_ok:
        mapped: Dict[ResidueChainKey, RestoreTarget] = {}
        for chain, snaps in out_by_chain.items():
            mapped.update(_zip_snapshot_targets(snaps, orig_by_chain[chain]))
        return mapped

    out_polymer, out_hetero = _split_polymer_hetero(output_res)
    orig_polymer, orig_hetero = _split_polymer_hetero(original_res)
    mapped = {}
    if out_polymer and len(out_polymer) == len(orig_polymer):
        mapped.update(_zip_snapshot_targets(out_polymer, orig_polymer))
    if out_hetero and len(out_hetero) == len(orig_hetero):
        mapped.update(_zip_snapshot_targets(out_hetero, orig_hetero))
    return mapped


def _pair_remum_to_original_pdb(
    remum_path: Optional[Path],
    original_pdb: Path,
) -> Dict[ResidueChainKey, Optional[int]]:
    """Order-match remum protein/hetero rows to the pre-cap PDB (keeps gaps)."""
    if remum_path is None or not remum_path.is_file():
        return {}
    text = remum_path.read_text(encoding="utf-8", errors="replace")
    rows = parse_residue_mapping_rows(text, path=remum_path)
    orig_by_chain: Dict[str, List[ResidueSnapshot]] = {}
    for snap in snapshot_pdb_residues(original_pdb):
        if snap[3] in CAP_RESNAMES:
            continue
        orig_by_chain.setdefault(snap[0], []).append(snap)

    remum_by_chain: Dict[str, List] = {}
    cap_finals: Dict[ResidueChainKey, None] = {}
    for row in rows:
        key = (row.chain, row.final_resid)
        if row.final_resname in CAP_RESNAMES or row.original_resname in CAP_RESNAMES:
            cap_finals[key] = None
            continue
        remum_by_chain.setdefault(row.chain, []).append(row)

    paired: Dict[ResidueChainKey, Optional[int]] = {key: None for key in cap_finals}
    for chain, remum_rows in remum_by_chain.items():
        originals = orig_by_chain.get(chain, [])
        remum_rows = sorted(remum_rows, key=lambda r: r.final_resid)
        if originals and len(originals) != len(remum_rows):
            logger.debug(
                "Residue count mismatch chain %s: remum %s vs original %s; "
                "using remum originals",
                chain,
                len(remum_rows),
                len(originals),
            )
            for row in remum_rows:
                paired[(row.chain, row.final_resid)] = row.original_resid
            continue
        for row, snap in zip(remum_rows, originals):
            paired[(row.chain, row.final_resid)] = snap[1]
    return paired


def _pair_snapshots_to_map(
    output_res: Sequence[ResidueSnapshot],
    original_res: Sequence[ResidueSnapshot],
) -> Dict[ResidueChainKey, Optional[int]]:
    mapped: Dict[ResidueChainKey, Optional[int]] = {}
    for key, target in _pair_snapshots_to_targets(output_res, original_res).items():
        mapped[key] = None if target is None else target[1]
    return mapped


def _rewrite_pdb_resids(
    src: Path,
    dest: Path,
    amber_to_original: Mapping[ResidueChainKey, RestoreTarget],
    *,
    assign_new_caps: bool,
    original_pdb: Optional[Path],
) -> Tuple[int, Dict[Tuple[str, str], int]]:
    icode_by_orig = _icode_lookup(original_pdb)

    with src.open("r", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()

    residue_plan: Dict[ResidueChainKey, Tuple[str, int]] = {}
    cap_keys: List[Tuple[ResidueChainKey, str]] = []
    polymer_originals: Dict[str, List[int]] = {}
    used: Dict[str, set] = {}

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        rec = _parse_pdb_atom(line)
        if rec is None:
            continue
        key: ResidueChainKey = (rec["chain"], rec["resid"])
        resname = rec["resname"]
        if key in residue_plan or any(ck == key for ck, _ in cap_keys):
            continue
        if assign_new_caps and resname in CAP_RESNAMES:
            cap_keys.append((key, resname))
            continue
        target = amber_to_original.get(key)
        if target is None:
            dest_chain, new_id = rec["chain"], rec["resid"]
        else:
            dest_chain, new_id = target
        residue_plan[key] = (dest_chain, new_id)
        used.setdefault(dest_chain, set()).add(new_id)
        if _is_polymer_non_cap(resname):
            polymer_originals.setdefault(dest_chain, []).append(new_id)

    cap_assignments: Dict[Tuple[str, str], int] = {}
    if assign_new_caps:
        for key, resname in cap_keys:
            dest_chain = key[0]
            chain_used = used.setdefault(dest_chain, set())
            polymer_ids = polymer_originals.get(dest_chain, [])
            if resname == "ACE":
                start = (min(polymer_ids) - 1) if polymer_ids else 0
                new_id = _next_free_down(chain_used, start)
            else:
                start = (max(polymer_ids) + 1) if polymer_ids else 1
                new_id = _next_free_up(chain_used, start)
            residue_plan[key] = (dest_chain, new_id)
            chain_used.add(new_id)
            cap_assignments[(dest_chain, resname)] = new_id

    rewritten = 0
    out_lines: List[str] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM", "TER")):
            out_lines.append(line)
            continue
        rec = _parse_pdb_atom(line)
        if rec is None:
            out_lines.append(line)
            continue
        key = (rec["chain"], rec["resid"])
        planned = residue_plan.get(key)
        if planned is None:
            out_lines.append(line)
            continue
        dest_chain, new_id = planned
        icode = icode_by_orig.get((dest_chain, new_id), rec["icode"])
        new_line = _rewrite_pdb_resid_line(
            line, new_id, icode=icode, chain=dest_chain
        )
        if new_line != line:
            rewritten += 1
        out_lines.append(new_line)

    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(out_lines)
    return rewritten, cap_assignments


def _icode_lookup(original_pdb: Optional[Path]) -> Dict[Tuple[str, int], str]:
    if original_pdb is None or not original_pdb.is_file():
        return {}
    lookup: Dict[Tuple[str, int], str] = {}
    for chain, resid, icode, _resname in snapshot_pdb_residues(original_pdb):
        lookup[(chain, resid)] = icode
    return lookup


def _is_polymer_non_cap(resname: str) -> bool:
    if resname in CAP_RESNAMES:
        return False
    return is_peptide_polymer_residue(resname)


def _next_free_down(used: Iterable[int], start: int) -> int:
    used_set = set(used)
    n = start
    while n in used_set:
        n -= 1
    return n


def _next_free_up(used: Iterable[int], start: int) -> int:
    used_set = set(used)
    n = start
    while n in used_set:
        n += 1
    return n


def _parse_pdb_atom(line: str) -> Optional[Dict[str, Any]]:
    if len(line) < 22:
        return None
    padded = line.rstrip("\n").ljust(27)
    resname = padded[17:20].strip()
    chain = padded[21:22]
    resid_s = padded[22:26].strip()
    icode = padded[26:27]
    if not resid_s:
        return None
    try:
        resid = int(resid_s)
    except ValueError:
        return None
    return {
        "resname": resname,
        "chain": chain if chain.strip() else " ",
        "resid": resid,
        "icode": icode,
    }


def _rewrite_pdb_resid_line(
    line: str,
    new_resid: int,
    *,
    icode: Optional[str] = None,
    chain: Optional[str] = None,
) -> str:
    nl = "\n" if line.endswith("\n") else ""
    body = line[:-1] if nl else line
    if len(body) < 27:
        body = body.ljust(27)
    resid_s = f"{int(new_resid):4d}"
    if len(resid_s) > 4:
        raise ValueError(f"Residue number {new_resid} does not fit PDB resSeq columns")
    ic = icode if icode is not None else (body[26:27] if len(body) > 26 else " ")
    if not ic:
        ic = " "
    ch = chain if chain is not None else (body[21:22] if len(body) > 21 else " ")
    if not ch:
        ch = " "
    else:
        ch = ch[0]
    return body[:21] + ch + resid_s + ic + body[27:] + nl


def stamp_atom_lines_from_snapshots(
    lines: Sequence[str],
    targets: Sequence[ResidueSnapshot],
) -> List[str]:
    """
    Rewrite ATOM/HETATM/TER chain + resid in residue order from *targets*.

    Used after tleap ``savePdb``, which sequentializes ids and blanks chains.
    Residue count must match; otherwise *lines* are returned unchanged.
    """
    seen: List[Tuple[str, int, str, str]] = []
    last: Optional[Tuple[str, int, str, str]] = None
    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        rec = _parse_pdb_atom(line)
        if rec is None:
            continue
        key = (rec["chain"], rec["resid"], rec["icode"], rec["resname"])
        if key != last:
            seen.append(key)
            last = key
    if len(seen) != len(targets):
        logger.warning(
            "Cannot stamp residue ids: tleap wrote %s residues, original had %s",
            len(seen),
            len(targets),
        )
        return list(lines)

    snap_i = -1
    last_key: Optional[Tuple[str, int, str, str]] = None
    stamped: List[str] = []
    for line in lines:
        if not line.startswith(("ATOM", "HETATM", "TER")):
            stamped.append(line)
            continue
        rec = _parse_pdb_atom(line)
        if rec is None:
            stamped.append(line)
            continue
        key = (rec["chain"], rec["resid"], rec["icode"], rec["resname"])
        if line.startswith(("ATOM", "HETATM")) and key != last_key:
            snap_i += 1
            last_key = key
        if snap_i < 0 or snap_i >= len(targets):
            stamped.append(line)
            continue
        ochain, oresid, oicode, _oname = targets[snap_i]
        stamped.append(
            _rewrite_pdb_resid_line(line, oresid, icode=oicode, chain=ochain)
        )
    return stamped
