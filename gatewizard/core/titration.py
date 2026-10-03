"""Titration figure data shared by the library and the GUI.

A residue is protonated at the target pH when pH < pKa.
``exclude_ids`` drop residues from panels A and B only.
``curve_ids`` chooses panel C. ``None`` keeps every curve.

Cysteines in a disulfide bond are omitted from panels A and B. PropKa writes
pKa 99.99 for those residues as a cysteine-bridge placeholder, not a titration
constant.
"""

from __future__ import annotations

import io
import math
import os
from collections import defaultdict
from typing import Any, Iterable

# Amber state → 1 protonated, 0 deprotonated. Unknown states stay at 0.5.
STATE_TO_Y = {
    "ASH": 1.0,
    "ASP": 0.0,
    "GLH": 1.0,
    "GLU": 0.0,
    "HIP": 1.0,
    "HIE": 0.0,
    "HID": 0.0,
    "LYS": 1.0,
    "LYN": 0.0,
    "CYS": 0.0,
    "CYM": 0.0,
    "CYX": 0.0,
    "TYR": 0.0,
    "TYM": 0.0,
    "ARG": 1.0,
}

TYPE_ORDER = ["ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG"]


def protonated_at_ph(pka: float, target_ph: float) -> bool:
    return float(target_ph) < float(pka)


def state_protonation(state: str) -> float:
    return float(STATE_TO_Y.get(str(state or "").strip().upper(), 0.5))


def _ph_grid(ph_min: float, ph_max: float, ph_step: float) -> list[float]:
    lo = float(ph_min)
    hi = float(ph_max)
    step = float(ph_step)
    if not (hi > lo):
        raise ValueError("pH max must be greater than pH min.")
    if step <= 0:
        raise ValueError("pH step must be positive.")
    if step < 0.05:
        step = 0.05
    values: list[float] = []
    current = lo
    while current <= hi + 1e-9 and len(values) < 400:
        values.append(round(current, 6))
        current += step
    if len(values) < 2:
        raise ValueError("pH range does not contain two points.")
    return values


def _curve_for_residue(manager: Any, residue: dict, ph_values: list[float]) -> list[float]:
    ys: list[float] = []
    for ph in ph_values:
        state = manager.get_default_protonation_state(residue, ph)
        ys.append(state_protonation(state))
    return ys


def _as_tokens(values: Iterable[str] | None) -> set[str]:
    if not values:
        return set()
    return {str(item).strip().upper() for item in values if str(item).strip()}


def residue_keys(row: dict) -> set[str]:
    """Ids a caller may use to name this residue: CYS77, CYS77:A, or the figure id."""
    residue = str(row.get("residue") or "").strip().upper()
    number = row.get("res_id")
    chain = str(row.get("chain") or "").strip().upper()
    base = f"{residue}{number}"
    keys = {
        str(row.get("id") or "").strip().upper(),
        str(row.get("label") or "").strip().upper(),
        base,
    }
    if chain and chain != "_":
        keys.add(f"{base}:{chain}")
    return {key for key in keys if key}


def _matches(row: dict, tokens: set[str]) -> bool:
    if not tokens:
        return False
    return bool(residue_keys(row) & tokens)


def assemble_titration(
    manager: Any,
    residues: list[dict],
    *,
    target_ph: float,
    ph_min: float,
    ph_max: float,
    ph_step: float,
) -> dict:
    """Residues and curves for the titration figure.

    get_ph_titration_curve is called so the pH grid matches the manager. Per-residue
    curves are then sampled with get_default_protonation_state so two chains that
    share a residue number stay separate. Ligands and non-positive residue numbers
    are dropped. A chain suffix is added to the label only when that residue number
    appears more than once.
    """
    ph_values = _ph_grid(ph_min, ph_max, ph_step)
    manager.get_ph_titration_curve((ph_values[0], ph_values[-1]), ph_step)

    protein = []
    for res in residues or []:
        try:
            res_id = int(res.get("res_id") or 0)
        except (TypeError, ValueError):
            continue
        if res_id <= 0:
            continue
        try:
            pka = float(res.get("pka"))
        except (TypeError, ValueError):
            continue
        protein.append(
            {
                "residue": str(res.get("residue") or "").strip(),
                "res_id": res_id,
                "chain": str(res.get("chain") or "").strip(),
                "pka": pka,
                "raw": res,
            }
        )

    label_counts: dict[str, int] = defaultdict(int)
    for row in protein:
        label_counts[f"{row['residue']}{row['res_id']}"] += 1

    out_residues = []
    curves = []
    for row in protein:
        base = f"{row['residue']}{row['res_id']}"
        chain = row["chain"]
        label = f"{base}:{chain}" if chain and label_counts[base] > 1 else base
        rid = f"{base}:{chain or '_'}"
        state = manager.get_default_protonation_state(row["raw"], target_ph)
        protonated = protonated_at_ph(row["pka"], target_ph)
        out_residues.append(
            {
                "id": rid,
                "label": label,
                "residue": row["residue"],
                "res_id": row["res_id"],
                "chain": chain,
                "pka": row["pka"],
                "state": str(state or ""),
                "protonated": bool(protonated),
            }
        )
        curves.append(
            {
                "id": rid,
                "label": label,
                "residue": row["residue"],
                "res_id": row["res_id"],
                "chain": chain,
                "pka": row["pka"],
                "ph": ph_values,
                "y": _curve_for_residue(manager, row["raw"], ph_values),
            }
        )

    return {
        "target_ph": float(target_ph),
        "ph_min": float(ph_values[0]),
        "ph_max": float(ph_values[-1]),
        "ph_step": float(ph_step),
        "residues": out_residues,
        "curves": curves,
    }


def titration_figure(
    payload: dict,
    *,
    exclude_ids: Iterable[str] | None = None,
    curve_ids: Iterable[str] | None = None,
    type_filter: Iterable[str] | None = None,
) -> dict:
    """Panels A/B/C from an assembled titration payload.

    ``exclude_ids`` omit residues from the pKa and protonation groups (A and B).
    Curves stay available so panel C can still show them.
    ``curve_ids`` is ``None`` for every curve, or a list of ids/labels to keep.
    ``type_filter`` limits A, B, and C to those residue names. Empty means all types.
    The token ``__none__`` hides every type.
    """
    residues = list(payload.get("residues") or [])
    curves = list(payload.get("curves") or [])
    target_ph = float(payload.get("target_ph", 7.0))
    excluded = _as_tokens(exclude_ids)
    types_wanted = _as_tokens(type_filter)
    hide_all = "__NONE__" in types_wanted
    if hide_all:
        types_wanted = set()

    def type_ok(name: str) -> bool:
        if hide_all:
            return False
        if not types_wanted:
            return True
        return str(name).strip().upper() in types_wanted

    kept = [
        row
        for row in residues
        if type_ok(row.get("residue")) and not _matches(row, excluded)
    ]
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in kept:
        grouped[str(row["residue"])].append(row)

    names = [name for name in TYPE_ORDER if name in grouped]
    names.extend(sorted(name for name in grouped if name not in TYPE_ORDER))
    types = []
    for name in names:
        rows = grouped[name]
        protonated = sum(1 for row in rows if protonated_at_ph(row["pka"], target_ph))
        types.append(
            {
                "name": name,
                "pka": [float(row["pka"]) for row in rows],
                "fraction": protonated / len(rows) if rows else 0.0,
                "points": [
                    {
                        "id": row["id"],
                        "pka": float(row["pka"]),
                        "protonated": protonated_at_ph(row["pka"], target_ph),
                    }
                    for row in rows
                ],
            }
        )

    if curve_ids is None:
        selected = [row for row in curves if type_ok(row.get("residue"))]
    else:
        wanted = _as_tokens(curve_ids)
        selected = [
            row
            for row in curves
            if type_ok(row.get("residue")) and _matches(row, wanted)
        ]

    return {
        "target_ph": target_ph,
        "types": types,
        "curves": selected,
        "excluded_ids": sorted(
            row["id"] for row in residues if _matches(row, excluded)
        ),
    }


# PropKa's cysteine-bridge placeholder. It is not a physical pKa.
_BRIDGE_PKA = 99.99
_BRIDGE_NAMES = {"CYS", "CYX"}


def disulfide_bridges(pdb_file: str, distance_threshold: float = 2.5) -> list[dict]:
    """Cysteine pairs that form a disulfide in a PDB file.

    A pair is kept when two SG atoms of CYS or CYX lie within
    ``distance_threshold`` angstroms, or when an SSBOND record names them.
    Only the first MODEL is read. Each bond is
    ``residue_a``, ``chain_a``, ``res_id_a`` and the same fields for ``b``,
    plus ``distance`` (angstroms, or None for an SSBOND with no coordinates)
    and ``source`` (``distance`` or ``ssbond``).
    """
    if not pdb_file or not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"PDB file {pdb_file} does not exist.")

    sulfurs: list[dict] = []
    records: list[tuple] = []
    seen_model = False
    with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("MODEL"):
                if seen_model:
                    break
                seen_model = True
                continue
            if line.startswith("ENDMDL") and seen_model:
                break
            if line.startswith("SSBOND") and len(line) >= 35:
                try:
                    records.append(
                        (
                            line[11:14].strip().upper() or "CYS",
                            line[15:16].strip().upper(),
                            int(line[17:21]),
                            line[25:28].strip().upper() or "CYS",
                            line[29:30].strip().upper(),
                            int(line[31:35]),
                        )
                    )
                except ValueError:
                    continue
                continue
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                continue
            if line[12:16].strip().upper() != "SG":
                continue
            name = line[17:20].strip().upper()
            if name not in _BRIDGE_NAMES:
                continue
            try:
                sulfurs.append(
                    {
                        "residue": name,
                        "chain": line[21:22].strip().upper(),
                        "res_id": int(line[22:26]),
                        "coords": (
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ),
                    }
                )
            except ValueError:
                continue

    bonds: dict[tuple, dict] = {}

    def _store(left: dict, right: dict, distance: float | None, source: str) -> None:
        key = tuple(
            sorted(
                (
                    (left["residue"], left["chain"], int(left["res_id"])),
                    (right["residue"], right["chain"], int(right["res_id"])),
                )
            )
        )
        if key[0] == key[1]:
            return
        current = bonds.get(key)
        if current and current.get("distance") is not None and distance is None:
            return
        a, b = key
        bonds[key] = {
            "residue_a": a[0],
            "chain_a": a[1],
            "res_id_a": a[2],
            "residue_b": b[0],
            "chain_b": b[1],
            "res_id_b": b[2],
            "distance": None if distance is None else round(float(distance), 3),
            "source": source if not current else current["source"],
        }
        if distance is not None:
            bonds[key]["distance"] = round(float(distance), 3)
            bonds[key]["source"] = "distance"

    cutoff = float(distance_threshold)
    for i, left in enumerate(sulfurs):
        for right in sulfurs[i + 1 :]:
            x1, y1, z1 = left["coords"]
            x2, y2, z2 = right["coords"]
            distance = ((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2) ** 0.5
            if distance <= cutoff:
                _store(left, right, distance, "distance")

    for name_a, chain_a, num_a, name_b, chain_b, num_b in records:
        _store(
            {"residue": name_a, "chain": chain_a, "res_id": num_a},
            {"residue": name_b, "chain": chain_b, "res_id": num_b},
            None,
            "ssbond",
        )

    return list(bonds.values())


def _site_tokens(residue: str, res_id: int, chain: str) -> set[str]:
    tokens = set()
    chain_id = str(chain or "").strip().upper()
    for name in _BRIDGE_NAMES:
        base = f"{name}{int(res_id)}"
        if chain_id:
            tokens.add(f"{base}:{chain_id}")
        else:
            tokens.add(base)
    return tokens


def _site_label(residue: str, res_id: int, chain: str) -> str:
    base = f"{str(residue or 'CYS').strip().upper()}{int(res_id)}"
    chain_id = str(chain or "").strip().upper()
    return f"{base}:{chain_id}" if chain_id else base


def disulfide_report(residues: list[dict], bridges: list[dict]) -> dict:
    """Residues to drop from panels A and B, and the warning that explains why.

    A residue is hidden when it is one side of a detected disulfide, or when
    PropKa gave that cysteine the 99.99 cysteine-bridge placeholder.
    ``warning`` is empty when nothing is hidden.
    """
    rows = list(residues or [])
    hidden: list[dict] = []
    used: set[str] = set()

    def _take(row: dict, partner: str, distance: float | None, source: str) -> None:
        rid = str(row.get("id") or "")
        if not rid or rid in used:
            return
        used.add(rid)
        hidden.append(
            {
                "id": rid,
                "label": str(row.get("label") or rid),
                "pka": float(row.get("pka")),
                "partner": partner,
                "distance": distance,
                "source": source,
            }
        )

    for bond in bridges or []:
        tokens_a = _site_tokens(bond["residue_a"], bond["res_id_a"], bond["chain_a"])
        tokens_b = _site_tokens(bond["residue_b"], bond["res_id_b"], bond["chain_b"])
        rows_a = [row for row in rows if _matches(row, tokens_a)]
        rows_b = [row for row in rows if _matches(row, tokens_b)]
        label_a = (
            str(rows_a[0].get("label") or "")
            if rows_a
            else _site_label(bond["residue_a"], bond["res_id_a"], bond["chain_a"])
        )
        label_b = (
            str(rows_b[0].get("label") or "")
            if rows_b
            else _site_label(bond["residue_b"], bond["res_id_b"], bond["chain_b"])
        )
        distance = bond.get("distance")
        source = str(bond.get("source") or "distance")
        ids_a = {id(row) for row in rows_a}
        for row in rows_a:
            _take(row, label_b, distance, source)
        for row in rows_b:
            if id(row) in ids_a:
                continue
            _take(row, label_a, distance, source)

    for row in rows:
        if str(row.get("residue") or "").strip().upper() not in _BRIDGE_NAMES:
            continue
        try:
            pka = float(row.get("pka"))
        except (TypeError, ValueError):
            continue
        if abs(pka - _BRIDGE_PKA) > 0.011:
            continue
        _take(row, "", None, "propka")

    return {
        "bridges": list(bridges or []),
        "hidden": hidden,
        "warning": _disulfide_warning_text(hidden),
    }


def _disulfide_warning_text(hidden: list[dict]) -> str:
    if not hidden:
        return ""
    lines = ["Omitted from panels A and B:"]
    seen_pairs: set[tuple[str, str]] = set()
    for item in hidden:
        partner = str(item.get("partner") or "").strip()
        distance = item.get("distance")
        label = str(item.get("label") or item.get("id"))
        if partner:
            key = tuple(sorted((label, partner)))
            if key in seen_pairs:
                continue
            seen_pairs.add(key)
            if distance is not None:
                lines.append(f"- {label} and {partner} (SG–SG {float(distance):.2f} Å)")
            else:
                lines.append(f"- {label} and {partner} (SSBOND record)")
        else:
            lines.append(f"- {label} (cysteine bridge)")
    lines.extend(
        [
            "",
            "These cysteines form disulfide bonds, so they were left out of panels A and B. "
            "PropKa does not calculate their pKa: if the sulfur is in a bridge it writes "
            "99.99 and stops. That marker is not used in these panels.",
        ]
    )
    return "\n".join(lines)


def run_titration(
    path: str,
    *,
    target_ph: float = 7.0,
    ph_min: float = 0.0,
    ph_max: float = 14.0,
    ph_step: float = 0.5,
    exclude_ids: Iterable[str] | None = None,
    curve_ids: Iterable[str] | None = None,
    type_filter: Iterable[str] | None = None,
    manager: Any | None = None,
) -> dict:
    """Run PropKa and return residues, every curve, and the filtered figure.

    Cysteines in a disulfide (SG–SG distance or an SSBOND record), and any
    cysteine PropKa marked with the 99.99 bridge placeholder, are added to
    ``exclude_ids`` so panels A and B omit them. ``disulfide_warning`` explains
    why. The full residue list still includes them.
    """
    if manager is None:
        from gatewizard.core.preparation import PreparationManager

        manager = PreparationManager(propka_version="3")
    manager.run_analysis(path)
    summary = manager.extract_summary(manager.last_analysis_file)
    residues = manager.parse_summary(summary)
    payload = assemble_titration(
        manager,
        residues,
        target_ph=target_ph,
        ph_min=ph_min,
        ph_max=ph_max,
        ph_step=ph_step,
    )
    payload["pdb_path"] = path
    report = disulfide_report(payload["residues"], disulfide_bridges(path))
    payload["disulfides"] = report["bridges"]
    payload["disulfide_hidden"] = report["hidden"]
    payload["disulfide_warning"] = report["warning"]
    locked = [row["id"] for row in report["hidden"]]
    payload["figure"] = titration_figure(
        payload,
        exclude_ids=[*(exclude_ids or []), *locked],
        curve_ids=curve_ids,
        type_filter=type_filter,
    )
    return payload


_MPL_LINESTYLE = {
    "solid": "-",
    "dashed": "--",
    "dotted": ":",
    "dashdot": "-.",
}


def _as_float(value: Any, default: float) -> float:
    try:
        n = float(value)
    except (TypeError, ValueError):
        return default
    return n


def render_titration_png(spec: dict[str, Any]) -> bytes:
    """Draw panels A (boxes), B (bars), and C (curves) into one PNG."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np

    style = spec.get("style") if isinstance(spec.get("style"), dict) else {}
    show = spec.get("show") if isinstance(spec.get("show"), dict) else {}
    show_a = show.get("a", True) is not False
    show_b = show.get("b", True) is not False
    show_c = show.get("c", True) is not False
    if not (show_a or show_b or show_c):
        raise ValueError("At least one titration panel must be visible.")

    titles = spec.get("titles") if isinstance(spec.get("titles"), dict) else {}
    labels = spec.get("labels") if isinstance(spec.get("labels"), dict) else {}
    letters = spec.get("letters") if isinstance(spec.get("letters"), dict) else {}
    types = spec.get("types") if isinstance(spec.get("types"), list) else []
    curves = spec.get("curves") if isinstance(spec.get("curves"), list) else []
    ref_lines = spec.get("reference_lines") if isinstance(spec.get("reference_lines"), list) else []
    target_ph = _as_float(spec.get("target_ph"), 7.0)
    dpi = int(_as_float(style.get("dpi"), 200))
    text = str(style.get("text_color") or "#222222")
    plot_bg = str(style.get("plot_bg") or "white")
    fig_bg = str(style.get("fig_bg") or plot_bg or "white")
    font = str(style.get("font_family") or "sans-serif").split(",")[0].strip() or "sans-serif"

    layout = spec.get("layout") if isinstance(spec.get("layout"), dict) else {}
    markers = spec.get("markers") if isinstance(spec.get("markers"), dict) else {}
    aspect_a = max(0.3, _as_float(layout.get("aspect_a"), 2.0))
    aspect_b = max(0.3, _as_float(layout.get("aspect_b"), 3.2))
    gap_px = max(0.0, _as_float(layout.get("gap"), 4.0))
    pka_marker = max(1.2, _as_float(markers.get("pka"), 5.0))
    curve_marker = max(1.2, _as_float(markers.get("curve"), 4.0))
    spine_w = max(0.4, _as_float(style.get("spine_width"), 1.0))
    line_w = max(0.4, _as_float(style.get("line_width"), 1.6))
    _line_styles = {"dashed": "--", "dotted": ":", "dashdot": "-."}
    mpl_ls = _line_styles.get(str(style.get("line_style") or "solid"), "-")
    fonts = spec.get("fonts") if isinstance(spec.get("fonts"), dict) else {}

    side_by_side = show_c and (show_a or show_b)
    pixels = layout.get("pixels") if isinstance(layout.get("pixels"), dict) else {}

    def _panel_px(key: str) -> tuple[float, float] | None:
        row = pixels.get(key) if isinstance(pixels.get(key), dict) else None
        if not row:
            return None
        w = _as_float(row.get("w"), 0.0)
        h = _as_float(row.get("h"), 0.0)
        if w < 8 or h < 8:
            return None
        return w, h

    pa, pb, pc = _panel_px("a"), _panel_px("b"), _panel_px("c")
    col_gap = max(0.0, _as_float(layout.get("column_gap"), 12.0))
    use_pixels = side_by_side and (pa or pb) and pc
    if use_pixels:
        left_w = (pa or pb or (1.0, 1.0))[0]
        right_w = pc[0] if pc else left_w
        ha = pa[1] if pa else 0.001
        hb = pb[1] if pb else 0.001
        fig_w = max(4.0, (left_w + col_gap + right_w) / 100.0)
        fig_h = max(3.0, ((pa[1] if pa else 0.0) + gap_px + (pb[1] if pb else 0.0)) / 100.0)
    else:
        height_a = max(0.25, _as_float(layout.get("height_a"), 1.0 / aspect_a))
        height_b = max(0.25, _as_float(layout.get("height_b"), 1.0 / aspect_b))
        left_share = min(80.0, max(20.0, _as_float(layout.get("left_share"), 52.0)))
        left_w, right_w = left_share, max(15.0, 100.0 - left_share)
        ha, hb = height_a, height_b
        fig_w = 11.2 if side_by_side else 7.2
        fig_h = 7.2 if side_by_side else (5.6 if show_a and show_b else 4.4)
    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    if fig_bg and fig_bg.lower() != "none":
        fig.patch.set_facecolor(fig_bg)

    extra_l = _as_float(style.get("extra_left"), 0.0) / max(fig_w * dpi, 1)
    extra_r = _as_float(style.get("extra_right"), 0.0) / max(fig_w * dpi, 1)
    extra_t = _as_float(style.get("extra_top"), 0.0) / max(fig_h * dpi, 1)
    extra_b = _as_float(style.get("extra_bottom"), 0.0) / max(fig_h * dpi, 1)
    avg_w = max((left_w + right_w) / 2.0, 1.0)
    wspace = col_gap / avg_w if use_pixels else 0.22
    tight = 0.012 if use_pixels else 0.08
    # Figure inches are screen CSS pixels / 100, so N px is N * 0.72 pt.
    font_scale = 72.0 / 100.0 if use_pixels else 1.0
    if use_pixels:
        pka_marker *= font_scale
        curve_marker *= font_scale

    def _font_pt(raw: Any, default_px: float) -> float:
        try:
            n = float(raw)
        except (TypeError, ValueError):
            n = float(default_px)
        if not (n > 0):
            n = float(default_px)
        return n * font_scale

    def _face(key: str) -> tuple[float, float, float, bool]:
        row = fonts.get(key) if isinstance(fonts.get(key), dict) else {}
        bold = row.get("axis_font_bold") is True
        return (
            _font_pt(row.get("axis_font_size"), 11),
            _font_pt(row.get("title_font_size"), 13),
            _font_pt(row.get("legend_font_size"), 11),
            bold,
        )

    def _tick_values(lo: float, hi: float, step_raw: Any, count_raw: Any) -> list[float]:
        step = _as_float(step_raw, 0.0)
        if step > 0 and hi > lo:
            start = math.ceil((lo - 1e-12) / step) * step
            values: list[float] = []
            v = start
            while v <= hi + 1e-12 and len(values) < 30:
                values.append(v)
                v += step
            if len(values) >= 2:
                return values
        count = int(_as_float(count_raw, 5) or 5)
        count = max(2, min(20, count))
        if hi == lo:
            return [lo]
        return [lo + (hi - lo) * i / (count - 1) for i in range(count)]

    if side_by_side:
        gs = fig.add_gridspec(
            3,
            2,
            height_ratios=[ha if show_a else 0.001, max(gap_px, 0.001), hb if show_b else 0.001],
            width_ratios=[max(0.2, left_w), max(0.2, right_w)],
            wspace=wspace,
            hspace=0.0,
            left=min(0.22, tight + max(0.0, extra_l)),
            right=max(0.78, (0.995 if use_pixels else 0.98) - max(0.0, extra_r)),
            top=max(0.78, (0.995 if use_pixels else 0.94) - max(0.0, extra_t)),
            bottom=min(0.22, tight + max(0.0, extra_b)),
        )
        slot_a = gs[0, 0] if show_b else gs[:, 0]
        slot_b = gs[2, 0] if show_a else gs[:, 0]
        slot_c = gs[:, 1]
    elif show_a and show_b:
        gs = fig.add_gridspec(
            3,
            1,
            height_ratios=[ha, max(gap_px, 0.001), hb],
            hspace=0.0,
            left=0.1,
            right=0.98,
            top=0.96,
            bottom=0.06,
        )
        slot_a, slot_b, slot_c = gs[0, 0], gs[2, 0], None
    else:
        gs = fig.add_gridspec(1, 1, left=0.1, right=0.98, top=0.94, bottom=0.08)
        only = gs[0, 0]
        slot_a = only if show_a else None
        slot_b = only if show_b else None
        slot_c = only if show_c else None

    def _letter(ax, key: str) -> None:
        letter = str(letters.get(key) or "").strip()
        if not letter:
            return
        _, title_pt, _, _ = _face(key)
        ax.text(
            -0.08,
            1.06,
            letter,
            transform=ax.transAxes,
            fontsize=title_pt + 3 * font_scale,
            fontweight="bold",
            color=text,
            fontfamily=font,
            clip_on=False,
        )

    def _style_ax(ax) -> None:
        ax.set_facecolor(plot_bg if plot_bg.lower() != "none" else "white")
        ax.tick_params(colors=text)
        ax.xaxis.label.set_color(text)
        ax.yaxis.label.set_color(text)
        ax.title.set_color(text)
        spine_on = {
            "left": style.get("spine_left", True) is not False,
            "bottom": style.get("spine_bottom", True) is not False,
            "top": style.get("spine_top", False) is True,
            "right": style.get("spine_right", False) is True,
        }
        for name, spine in ax.spines.items():
            spine.set_color(text)
            spine.set_linewidth(spine_w)
            spine.set_linestyle(mpl_ls)
            spine.set_visible(spine_on.get(name, False))
        if style.get("show_grid", True) is not False:
            ax.grid(True, axis="y", linestyle=":", alpha=0.35, color=text)

    def _refs(ax, axis: str) -> None:
        for line in ref_lines:
            if not isinstance(line, dict):
                continue
            if str(line.get("axis") or "y") != axis:
                continue
            value = _as_float(line.get("value"), float("nan"))
            if value != value:
                continue
            color = str(line.get("color") or "#e74c3c")
            ls = _MPL_LINESTYLE.get(str(line.get("style") or "dashed"), "--")
            lw = _as_float(line.get("width"), 1.5)
            label = str(line.get("label") or "")
            if axis == "y":
                ax.axhline(value, color=color, linestyle=ls, linewidth=lw, label=label or None, zorder=2)
            else:
                ax.axvline(value, color=color, linestyle=ls, linewidth=lw, label=label or None, zorder=2)

    if show_a and slot_a is not None:
        ax1 = fig.add_subplot(slot_a)
        names = [str(t.get("name") or "") for t in types]
        pka_lists = [list(t.get("pka") or []) for t in types]
        colors = [str(t.get("color") or "#7f8c8d") for t in types]
        if names and any(pka_lists):
            box_kwargs = dict(
                patch_artist=True,
                showmeans=True,
                medianprops={"color": "black", "linewidth": 1.4},
                meanprops={
                    "marker": "D",
                    "markerfacecolor": "#c0392b",
                    "markeredgecolor": "#7b241c",
                    "markersize": 5,
                },
            )
            try:
                bp = ax1.boxplot(pka_lists, tick_labels=names, **box_kwargs)
            except TypeError:
                bp = ax1.boxplot(pka_lists, labels=names, **box_kwargs)
            for patch, color in zip(bp["boxes"], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.55)
            rng = np.random.default_rng(1)
            for i, group in enumerate(types):
                points = group.get("points") or []
                if not points:
                    continue
                xs = rng.normal(i + 1, 0.04, size=len(points))
                for x, point in zip(xs, points):
                    filled = bool(point.get("protonated"))
                    y = _as_float(point.get("pka"), float("nan"))
                    if y != y:
                        continue
                    if filled:
                        ax1.scatter(
                            x, y, s=pka_marker ** 2, color=colors[i], edgecolors="black", linewidths=0.6, zorder=3
                        )
                    else:
                        ax1.scatter(
                            x,
                            y,
                            s=pka_marker ** 2,
                            facecolors="none",
                            edgecolors=colors[i],
                            linewidths=1.1,
                            zorder=3,
                        )
        axis_pt, title_pt, legend_pt, bold = _face("a")
        weight = "bold" if bold else "normal"
        ax1.set_ylabel(str(labels.get("a_y") or "pKa"), fontfamily=font, fontsize=axis_pt, fontweight=weight)
        ax1.set_title(
            str(titles.get("a") or f"pKa distribution at pH {target_ph:g}"),
            fontfamily=font,
            fontsize=title_pt,
            fontweight=weight,
        )
        y0 = spec.get("a_ymin")
        y1 = spec.get("a_ymax")
        if y0 is not None and y1 is not None and str(y0) != "" and str(y1) != "":
            ax1.set_ylim(_as_float(y0, 0), _as_float(y1, 20))
        _refs(ax1, "y")
        handles = [
            Line2D([0], [0], marker="o", color="none", markerfacecolor="gray", markeredgecolor="black", label="Protonated (pH < pKa)"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="none", markeredgecolor="gray", label="Deprotonated (pH ≥ pKa)"),
        ]
        ax1.legend(handles=handles, loc="upper left", fontsize=legend_pt, frameon=False)
        _style_ax(ax1)
        ax1.tick_params(labelsize=axis_pt)
        for lbl in list(ax1.get_xticklabels()) + list(ax1.get_yticklabels()):
            lbl.set_fontweight(weight)
        _letter(ax1, "a")

    if show_b and slot_b is not None:
        ax2 = fig.add_subplot(slot_b)
        names = [str(t.get("name") or "") for t in types]
        fracs = [_as_float(t.get("fraction"), 0.0) for t in types]
        colors = [str(t.get("color") or "#7f8c8d") for t in types]
        if names:
            bars = ax2.bar(
                range(len(names)),
                fracs,
                color=colors,
                alpha=0.75,
                edgecolor="black",
                linewidth=0.8,
            )
            ax2.set_xticks(range(len(names)))
            ax2.set_xticklabels(names)
            for bar, frac in zip(bars, fracs):
                ax2.text(
                    bar.get_x() + bar.get_width() / 2,
                    bar.get_height() + 0.02,
                    f"{frac * 100:.0f}%",
                    ha="center",
                    va="bottom",
                    fontsize=_face("b")[2],
                    fontweight="bold",
                    color=text,
                    fontfamily=font,
                )
        ax2.set_ylim(0, 1.15)
        by0, by1 = spec.get("b_ymin"), spec.get("b_ymax")
        if by0 is not None and by1 is not None and str(by0) != "" and str(by1) != "":
            ax2.set_ylim(_as_float(by0, 0), _as_float(by1, 1.15))
        axis_pt, title_pt, _, bold = _face("b")
        weight = "bold" if bold else "normal"
        ax2.set_ylabel(
            str(labels.get("b_y") or "Protonated fraction"),
            fontfamily=font,
            fontsize=axis_pt,
            fontweight=weight,
        )
        ax2.set_xlabel(
            str(labels.get("b_x") or "Residue type"),
            fontfamily=font,
            fontsize=axis_pt,
            fontweight=weight,
        )
        ax2.set_title(
            str(titles.get("b") or f"Protonation state distribution at pH {target_ph:g}"),
            fontfamily=font,
            fontsize=title_pt,
            fontweight=weight,
        )
        _refs(ax2, "y")
        _style_ax(ax2)
        ax2.tick_params(labelsize=axis_pt)
        for lbl in list(ax2.get_xticklabels()) + list(ax2.get_yticklabels()):
            lbl.set_fontweight(weight)
        _letter(ax2, "b")

    if show_c and slot_c is not None:
        ax3 = fig.add_subplot(slot_c)
        for curve in curves:
            ph = [ _as_float(v, float("nan")) for v in (curve.get("ph") or []) ]
            yy = [ _as_float(v, float("nan")) for v in (curve.get("y") or []) ]
            n = min(len(ph), len(yy))
            if n < 2:
                continue
            ax3.plot(
                ph[:n],
                yy[:n],
                color=str(curve.get("color") or "#333333"),
                linewidth=line_w * font_scale,
                linestyle=mpl_ls,
                marker="o",
                markersize=curve_marker,
                label=str(curve.get("label") or curve.get("id") or ""),
            )
        axis_pt, title_pt, legend_pt, bold = _face("c")
        weight = "bold" if bold else "normal"
        ax3.set_xlabel(str(labels.get("c_x") or "pH"), fontfamily=font, fontsize=axis_pt, fontweight=weight)
        ax3.set_ylabel(
            str(labels.get("c_y") or "Protonation state (1=protonated, 0=deprotonated)"),
            fontfamily=font,
            fontsize=axis_pt,
            fontweight=weight,
        )
        ax3.set_title(
            str(titles.get("c") or "Protein titration curves"),
            fontfamily=font,
            fontsize=title_pt,
            fontweight=weight,
        )
        ax3.set_ylim(-0.05, 1.08)
        y0 = spec.get("c_ymin")
        y1 = spec.get("c_ymax")
        if y0 is not None and y1 is not None and str(y0) != "" and str(y1) != "":
            lo = _as_float(y0, 0.0)
            hi = _as_float(y1, 1.0)
            if hi < lo:
                lo, hi = hi, lo
            pad = _as_float(spec.get("c_ypad"), 0.05)
            if pad < 0:
                pad = 0.0
            ax3.set_ylim(lo - pad, hi + pad)
            ax3.set_yticks(_tick_values(lo, hi, spec.get("c_ytick_step"), spec.get("c_ytick_count")))
            decimals = spec.get("c_ytick_decimals")
            if decimals is not None and str(decimals).strip() != "":
                places = max(0, min(6, int(_as_float(decimals, 2))))
                ax3.yaxis.set_major_formatter(plt.FormatStrFormatter(f"%.{places}f"))
        x0, x1 = spec.get("c_xmin"), spec.get("c_xmax")
        if x0 is not None and x1 is not None and str(x0) != "" and str(x1) != "":
            ax3.set_xlim(_as_float(x0, 0), _as_float(x1, 14))
        _refs(ax3, "x")
        if curves:
            ax3.legend(
                fontsize=legend_pt,
                frameon=style.get("legend_frame_c", True) is not False,
                loc="best",
            )
        _style_ax(ax3)
        ax3.tick_params(labelsize=axis_pt)
        for lbl in list(ax3.get_xticklabels()) + list(ax3.get_yticklabels()):
            lbl.set_fontweight(weight)
        _letter(ax3, "c")

    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=dpi, facecolor=fig.get_facecolor())
    plt.close(fig)
    return buf.getvalue()
