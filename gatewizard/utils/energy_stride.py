"""Shared per-file stride helpers for energetic (log) analysis."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union


def lookup_file_map(file_map: Optional[Dict[str, Any]], path: Path) -> Any:
    """Look up a per-file value using basename, with case-insensitive fallback."""
    if not file_map:
        return None
    name = path.name
    if name in file_map:
        return file_map[name]
    # Also try full path / as-stored keys
    sp = str(path)
    if sp in file_map:
        return file_map[sp]
    name_lower = name.lower()
    for key, val in file_map.items():
        key_base = Path(str(key)).name.lower()
        if key_base == name_lower or str(key).lower() == name_lower:
            return val
    return None


def energy_keep_indices(
    n_points: int,
    log_files: Sequence[Union[str, Path]],
    file_ranges: Dict[str, tuple],
    file_strides: Optional[Dict[str, int]],
) -> Optional[List[int]]:
    """
    Build indices to keep after applying per-file strides.

    Returns None when no striding is needed (all strides ≤ 1 or empty map).
    ``file_ranges`` values are ``(start_idx, end_idx, ...)``.
    """
    if n_points <= 0 or not file_strides:
        return None
    if not any(max(1, int(v or 1)) > 1 for v in file_strides.values()):
        return None

    keep: List[int] = []
    for log_file in log_files:
        path = Path(log_file)
        rng = file_ranges.get(str(path))
        if rng is None:
            rng = file_ranges.get(path.name)
        if rng is None:
            # case-insensitive basename
            for key, val in file_ranges.items():
                if Path(key).name.lower() == path.name.lower():
                    rng = val
                    break
        if rng is None:
            continue
        start_idx = int(rng[0])
        end_idx = int(rng[1])
        stride = max(1, int(lookup_file_map(file_strides, path) or 1))
        keep.extend(range(start_idx, end_idx, stride))

    if not keep or len(keep) >= n_points:
        return None
    return keep


def apply_energy_stride_to_result(
    result: Dict[str, Any],
    log_files: Sequence[Union[str, Path]],
    file_ranges: Dict[str, tuple],
    file_strides: Optional[Dict[str, int]],
) -> Dict[str, Any]:
    """Subsample energetic analysis result arrays using per-file strides."""
    import numpy as np

    x = np.asarray(result.get("x") or [], dtype=float)
    keep = energy_keep_indices(len(x), log_files, file_ranges, file_strides)
    if keep is None:
        return result

    idx = np.asarray(keep, dtype=int)
    new_series = []
    for s in result.get("series") or []:
        y = np.asarray(s.get("y") or [], dtype=float)
        if len(y) != len(x):
            new_series.append(s)
            continue
        y2 = y[idx]
        new_series.append({**s, "y": y2.tolist()})
        # refresh statistics for this series key when present
        key = s.get("key")
        if key and isinstance(result.get("statistics"), dict) and key in result["statistics"]:
            valid = y2[np.isfinite(y2)]
            if len(valid):
                result["statistics"][key] = {
                    "mean": float(np.mean(valid)),
                    "std": float(np.std(valid)),
                    "min": float(np.min(valid)),
                    "max": float(np.max(valid)),
                    "initial": float(valid[0]),
                    "final": float(valid[-1]),
                }

    out = {**result, "x": x[idx].tolist(), "series": new_series}
    out["stride_applied"] = True
    out["n_points_before_stride"] = int(len(x))
    out["n_points_after_stride"] = int(len(idx))
    return out
