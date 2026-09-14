"""
Parse GateWizard / pdb4amber residue renumbering files.

Supported files (same 5-column layout, column meaning differs)::

1. ``*_gatewizard_residue_mapping.txt`` (capping / PropKa prep)::

       # Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID
       VAL A     21    VAL   2

2. ``*_protonated_renum.txt`` / ``system_for_tleap_renum.txt`` (pdb4amber)::

       VAL A      2    VAL  21
       # first id = new/final, second id = old/original

Caps may use ``-`` for a missing original id. The returned map is always
``final_resid → original_resid``.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence, Union

_MAPPING_LINE = re.compile(
    r"^\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$"
)


def _detect_column_order(text: str, path: Optional[Union[str, Path]] = None) -> str:
    """
    Return ``'original_final'`` or ``'final_original'``.
    """
    name = Path(path).name.lower() if path else ""
    if "gatewizard_residue_mapping" in name or "gatewizard_to_original" in name:
        return "original_final"
    if name.endswith("_renum.txt") or name.endswith("renum.txt") or "_renum." in name:
        return "final_original"

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


def parse_residue_mapping_text(
    text: str,
    *,
    path: Optional[Union[str, Path]] = None,
    column_order: Optional[str] = None,
) -> Dict[int, int]:
    """Parse mapping text into ``final_resid → original_resid``."""
    order = column_order or _detect_column_order(text, path)
    final_to_original: Dict[int, int] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        match = _MAPPING_LINE.match(line)
        if not match:
            continue
        id_a_s, id_b_s = match.group(3), match.group(5)
        if order == "final_original":
            final_id_s, original_id_s = id_a_s, id_b_s
        else:
            original_id_s, final_id_s = id_a_s, id_b_s
        if original_id_s in {"-", "."}:
            continue
        try:
            original_id = int(original_id_s)
            final_id = int(final_id_s)
        except ValueError:
            continue
        final_to_original[final_id] = original_id
    return final_to_original


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
