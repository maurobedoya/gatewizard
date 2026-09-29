"""
Parse GateWizard / pdb4amber residue renumbering files.

Supported files (same 5-column layout, column meaning differs)::

1. ``*_gatewizard_residue_mapping.txt`` (capping / PropKa prep)::

       # Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID
       VAL A     21    VAL   2

2. pdb4amber / tleap ``*_renum.txt`` — column order varies by tool version::

       ALA A     10    ALA    1    # original then final (AmberTools pdb4amber)
       VAL A      2    VAL   21    # final then original (some tleap / GW files)

Caps may use ``-`` for a missing original id. The returned map is always
``final_resid → original_resid``. Detect order from the header or from
which column is systematically larger; do not trust the filename alone.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Mapping, NamedTuple, Optional, Sequence, Tuple, Union

_MAPPING_LINE = re.compile(
    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$"
)
# pdb4amber drops the chain column when tleap blanked it: "MET 1 MET 1"
_MAPPING_LINE_NO_CHAIN = re.compile(
    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$"
)

# (chain, final_resid) → original resid; ``None`` means a new cap (``-``).
ResidueChainKey = Tuple[str, int]


class ResidueMappingRow(NamedTuple):
    """One 5-column mapping line, ids already oriented as final → original."""

    original_resname: str
    chain: str
    original_resid: Optional[int]
    final_resname: str
    final_resid: int


def _detect_column_order(text: str, path: Optional[Union[str, Path]] = None) -> str:
    """
    Return ``'original_final'`` or ``'final_original'``.
    """
    name = Path(path).name.lower() if path else ""
    if "gatewizard_residue_mapping" in name or "gatewizard_to_original" in name:
        return "original_final"

    header = "\n".join(
        line for line in text.splitlines()[:12] if line.strip().startswith("#")
    ).upper()
    if "ORIGINAL" in header and "FINAL" in header:
        return "original_final"

    # Content heuristic on pairs where the two ids differ.
    orig_gt = 0
    final_gt = 0
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = _MAPPING_LINE.match(line)
        if not match:
            continue
        a_s, b_s = match.group(3), match.group(5)
        if a_s in {"-", "."} or b_s in {"-", "."}:
            continue
        try:
            a, b = int(a_s), int(b_s)
        except ValueError:
            continue
        if a == b:
            continue
        if a > b:
            orig_gt += 1
        else:
            final_gt += 1
    if orig_gt > final_gt:
        return "original_final"
    if final_gt > orig_gt:
        return "final_original"
    # Identity / empty: either order yields the same map.
    return "original_final"


def parse_residue_mapping_file(
    path: Union[str, Path],
) -> Dict[int, int]:
    """
    Load ``final_resid → original_resid`` from a renumbering / mapping file.

    Entries whose original id is ``-`` (ACE/NME caps) are skipped.
    If the same final id appears more than once, the last line wins.
    """
    p = Path(path).expanduser().resolve()
    text = p.read_text(encoding="utf-8", errors="replace")
    return parse_residue_mapping_text(text, path=p)


def parse_residue_mapping_rows(
    text: str,
    *,
    path: Optional[Union[str, Path]] = None,
    column_order: Optional[str] = None,
) -> List[ResidueMappingRow]:
    """Parse mapping text into rows. Caps keep ``original_resid=None``."""
    order = column_order or _detect_column_order(text, path)
    rows: List[ResidueMappingRow] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = _MAPPING_LINE.match(line)
        if match:
            orig_name, chain, id_a_s, final_name, id_b_s = match.groups()
        else:
            match = _MAPPING_LINE_NO_CHAIN.match(line)
            if not match:
                continue
            orig_name, id_a_s, final_name, id_b_s = match.groups()
            chain = " "
        if order == "final_original":
            final_id_s, original_id_s = id_a_s, id_b_s
            # remum: first resname is the final/Amber name
            original_resname, final_resname = final_name, orig_name
        else:
            original_id_s, final_id_s = id_a_s, id_b_s
            original_resname, final_resname = orig_name, final_name
        if original_id_s in {"-", "."}:
            original_id: Optional[int] = None
        else:
            try:
                original_id = int(original_id_s)
            except ValueError:
                continue
        try:
            final_id = int(final_id_s)
        except ValueError:
            continue
        rows.append(
            ResidueMappingRow(
                original_resname=original_resname,
                chain=chain,
                original_resid=original_id,
                final_resname=final_resname,
                final_resid=final_id,
            )
        )
    return rows


def parse_residue_mapping_text(
    text: str,
    *,
    path: Optional[Union[str, Path]] = None,
    column_order: Optional[str] = None,
) -> Dict[int, int]:
    """Parse mapping text into ``final_resid → original_resid``."""
    final_to_original: Dict[int, int] = {}
    for row in parse_residue_mapping_rows(text, path=path, column_order=column_order):
        if row.original_resid is None:
            continue
        final_to_original[row.final_resid] = row.original_resid
    return final_to_original


def parse_residue_mapping_by_chain(
    path: Union[str, Path],
) -> Dict[ResidueChainKey, Optional[int]]:
    """
    Load ``(chain, final_resid) → original_resid`` from a mapping file.

    Cap lines with original ``-`` are kept as ``None``. Two chains that share
    the same final id stay distinct.
    """
    p = Path(path).expanduser().resolve()
    text = p.read_text(encoding="utf-8", errors="replace")
    return parse_residue_mapping_text_by_chain(text, path=p)


def parse_residue_mapping_text_by_chain(
    text: str,
    *,
    path: Optional[Union[str, Path]] = None,
    column_order: Optional[str] = None,
) -> Dict[ResidueChainKey, Optional[int]]:
    """Parse mapping text into ``(chain, final_resid) → original_resid``."""
    by_chain: Dict[ResidueChainKey, Optional[int]] = {}
    for row in parse_residue_mapping_rows(text, path=path, column_order=column_order):
        by_chain[(row.chain, row.final_resid)] = row.original_resid
    return by_chain


def align_mapping_keys_to_pdb(
    by_chain: Mapping[ResidueChainKey, Optional[int]],
    pdb_keys: Sequence[ResidueChainKey],
) -> Dict[ResidueChainKey, Optional[int]]:
    """
    If a remum was parsed with the columns swapped, flip it so keys match the PDB.

    Compares ``(chain, final)`` hits vs ``(chain, original)`` hits against *pdb_keys*.
    """
    if not by_chain:
        return {}
    key_set = set(pdb_keys)
    current = dict(by_chain)
    hits = sum(1 for key in current if key in key_set)
    flipped: Dict[ResidueChainKey, Optional[int]] = {}
    for (chain, final_id), orig_id in current.items():
        if orig_id is None:
            flipped[(chain, final_id)] = None
            continue
        flipped[(chain, orig_id)] = final_id
    flip_hits = sum(1 for key in flipped if key in key_set)
    if flip_hits > hits:
        return flipped
    return current


def compose_amber_to_original(
    remum_by_chain: Mapping[ResidueChainKey, Optional[int]],
    cap_by_chain: Optional[Mapping[ResidueChainKey, Optional[int]]] = None,
) -> Dict[ResidueChainKey, Optional[int]]:
    """
    Map Amber (pdb4amber) ids to pre-cap originals.

    ``remum_by_chain`` is ``(chain, amber_id) → id_in_pdb4amber_input``.
    ``cap_by_chain`` is ``(chain, capped_id) → original_id`` (``None`` for ACE/NME).
    """
    if not cap_by_chain:
        return dict(remum_by_chain)
    composed: Dict[ResidueChainKey, Optional[int]] = {}
    for (chain, amber_id), capped_id in remum_by_chain.items():
        if capped_id is None:
            composed[(chain, amber_id)] = None
            continue
        cap_key = (chain, capped_id)
        if cap_key in cap_by_chain:
            composed[(chain, amber_id)] = cap_by_chain[cap_key]
        else:
            composed[(chain, amber_id)] = capped_id
    return composed


def should_apply_residue_mapping(
    resids: Sequence[Union[int, float]],
    final_to_original: Mapping[int, int],
) -> bool:
    """
    True when ``resids`` look like GateWizard finals (renumbered from ~1),
    not Charmm/PDB originals that already match the mapping targets.

    Overlay of Charmm (starts ~21) + GateWizard (starts at 2) can share one
    mapping file: only the renumbered series is remapped.
    """
    if not final_to_original:
        return False
    nums = [int(r) for r in resids if r is not None and r == r]  # drop NaN
    if not nums:
        return False
    start = min(nums)
    if start > 2:
        return False
    mapped = final_to_original.get(start)
    return mapped is not None and mapped != start


def remap_resids_to_original(
    resids: Sequence[Union[int, float]],
    final_to_original: Mapping[int, int],
    *,
    force: bool = False,
) -> List[int]:
    """
    Map topology (final) residue ids to original PDB numbers.

    Unmapped ids are left unchanged. When ``force`` is False:
    - GateWizard finals (start ≤ 2) are remapped
    - Already-original axes are left alone
    - Over-remapped full-length series (same length as the map) are rebuilt
      from sorted finals (Y order is still residue order)
    """
    nums = [int(r) for r in resids]
    if not nums or not final_to_original:
        return list(nums)
    if force:
        return [int(final_to_original.get(r, r)) for r in nums]
    if should_apply_residue_mapping(nums, final_to_original):
        return [int(final_to_original.get(r, r)) for r in nums]

    originals = list(final_to_original.values())
    min_orig = min(originals)
    orig_set = set(originals)
    start = min(nums)
    orig_hits = sum(1 for n in nums if n in orig_set)
    looks_original = start > 2 and (orig_hits >= len(nums) * 0.75 or start >= min_orig)
    if not looks_original:
        return list(nums)

    mono = all(nums[i] >= nums[i - 1] for i in range(1, len(nums)))
    if start == min_orig and mono and len(set(nums)) == len(nums):
        return list(nums)

    finals = sorted(final_to_original.keys())
    if len(nums) == len(finals):
        return [int(final_to_original[f]) for f in finals]
    return list(nums)


def remap_residue_type_labels(
    labels: Sequence[str],
    resids: Sequence[Union[int, float]],
    final_to_original: Mapping[int, int],
    *,
    force: bool = False,
) -> List[str]:
    """
    Rewrite ``ALA2``-style labels using original residue numbers.

    Prefers parallel ``resids`` when lengths match; otherwise parses a trailing
    integer from each label.
    """
    label_list = [str(x) for x in labels]
    if not label_list:
        return []
    use_resids = len(resids) == len(label_list)
    topology_ids: List[int] = []
    if use_resids:
        topology_ids = [int(r) for r in resids]
    else:
        for lab in label_list:
            m = re.search(r"(\d+)\s*$", lab)
            topology_ids.append(int(m.group(1)) if m else -1)

    if not force and not should_apply_residue_mapping(
        [r for r in topology_ids if r >= 0], final_to_original
    ):
        return list(label_list)

    out: List[str] = []
    for lab, rid in zip(label_list, topology_ids):
        if rid < 0 or rid not in final_to_original:
            out.append(lab)
            continue
        orig = final_to_original[rid]
        m = re.match(r"^(.*?)(\d+)\s*$", lab)
        if m:
            out.append(f"{m.group(1)}{orig}")
        else:
            out.append(str(orig))
    return out
