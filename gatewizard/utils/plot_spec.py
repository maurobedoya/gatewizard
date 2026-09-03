"""
Shared PlotSpec contract for energetic analysis plots (API + GUI export).
"""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

PLOT_SPEC_VERSION = 1

DEFAULT_GLOBAL_STYLE: Dict[str, Any] = {
    "time_units": "ns",
    "energy_units": "kcal/mol",
    "pressure_units": "atm",
    "temperature_units": "K",
    "volume_units": "Å³",
    "plot_bg": "#2b2b2b",
    "fig_bg": "#212121",
    "text_color": "Auto",
    "grid_color": None,
    "show_grid": True,
    "figsize": [10.0, 6.0],
    "dpi": 300,
    "font_family": "sans-serif",
    "xlabel": None,
    "ylabel": None,
    "title": None,
    "xlim": None,
    "ylim": None,
    "show_ticks": True,
    "tick_length": 4.0,
    "tick_width": 1.0,
    "spine_width": 1.0,
    "show_spine_left": True,
    "show_spine_bottom": True,
    "show_spine_top": False,
    "show_spine_right": False,
    "extra_left": 0.0,
    "extra_right": 0.0,
    "extra_top": 0.0,
    "extra_bottom": 0.0,
}

DEFAULT_LINE_COLORS = [
    "#61afef",
    "#98c379",
    "#e06c75",
    "#e5c07b",
    "#c678dd",
    "#56b6c2",
    "#d19a66",
    "#abb2bf",
    "blue",
    "red",
    "green",
    "orange",
]


def _as_pair(value: Any) -> Optional[Tuple[float, float]]:
    if value is None:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 2:
        try:
            return float(value[0]), float(value[1])
        except (TypeError, ValueError):
            return None
    return None


def _as_float(value: Any, default: Optional[float] = None) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def normalize_reference_lines(raw: Any) -> List[Dict[str, Any]]:
    """Normalize ``[{axis, value, color, width, style, label}, ...]`` (also hlines/vlines)."""
    out: List[Dict[str, Any]] = []
    if not isinstance(raw, list):
        return out
    for item in raw:
        axis = "y"
        style = "dashed"
        color = "#888888"
        width = 1.2
        label = ""
        value: Optional[float] = None
        if isinstance(item, dict):
            value = _as_float(item.get("value"))
            ax = str(item.get("axis") or "y").lower()
            axis = ax if ax in ("x", "y") else "y"
            st = str(item.get("style") or "dashed").lower()
            style = st if st in ("solid", "dashed", "dotted", "dashdot", "--", ":", "-.", "-") else "dashed"
            color = str(item.get("color") or color)
            width = _as_float(item.get("width"), 1.2) or 1.2
            label = str(item.get("label") or "")
        else:
            value = _as_float(item)
        if value is None:
            continue
        out.append(
            {
                "axis": axis,
                "value": value,
                "color": color,
                "width": width if width > 0 else 1.2,
                "style": style,
                "label": label,
            }
        )
    return out


def grid_spec_slices(
    n_panels: int,
    cols: int,
    last_row_align: str = "start",
) -> Tuple[List[Tuple[int, int, int]], int, int]:
    """Slices into a (rows × cols*2) micro-grid so a short last row can align.

    ``last_row_align`` is ``start``, ``center``, or ``end`` (right).
    Returns ``(slices, rows, micro_cols)`` where each slice is ``(row, c0, c1)``.
    """
    n_cols = max(1, int(cols))
    n = max(0, int(n_panels))
    rows = max(1, (n + n_cols - 1) // n_cols if n else 1)
    micro_cols = n_cols * 2
    slices: List[Tuple[int, int, int]] = []
    full = n // n_cols
    rem = n % n_cols
    for i in range(full * n_cols):
        r = i // n_cols
        c = i % n_cols
        slices.append((r, c * 2, c * 2 + 2))
    if rem:
        r = full
        if last_row_align == "center":
            start = n_cols - rem
        elif last_row_align == "end":
            start = 2 * (n_cols - rem)
        else:
            start = 0
        for k in range(rem):
            c0 = start + k * 2
            slices.append((r, c0, c0 + 2))
    return slices, rows, micro_cols


def _clamp_legend_fontsize(raw: Any, default: float = 8.0) -> float:
    """Map GUI px (often 10–40) to matplotlib points so figure legends stay readable."""
    n = _as_float(raw, default)
    if n is None:
        n = default
    if n > 14:
        n = max(6.0, min(11.0, round(n * 0.28 + 5.2)))
    return max(6.0, min(11.0, float(n)))


def _normalize_legend(raw: Any) -> Dict[str, Any]:
    src = raw if isinstance(raw, dict) else {}
    mode = str(src.get("mode") or "").lower()
    if mode not in ("each", "one", "outside", "none", ""):
        mode = ""
    loc = str(src.get("loc") or "bottom").lower()
    if loc not in ("top", "bottom", "left", "right"):
        loc = "bottom"
    entries = str(src.get("entries") or "sets").lower()
    if entries not in ("sets", "roles", "both"):
        entries = "sets"
    fontsize = _clamp_legend_fontsize(src.get("fontsize"), 8.0)
    ncol = int(src.get("ncol") or 1)
    ncol = max(1, min(ncol, 8))
    cell = int(src.get("cell") or 0)
    return {
        "mode": mode,
        "cell": max(0, cell),
        "loc": loc,
        "entries": entries,
        "fontsize": fontsize,
        "ncol": ncol,
        "title": str(src.get("title") or ""),
    }


def normalize_plot_spec(spec: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Return a validated PlotSpec dict with defaults filled in."""
    src = deepcopy(spec or {})
    global_style = {**DEFAULT_GLOBAL_STYLE, **(src.get("global") or {})}
    panels_in = src.get("panels") or []
    panels: List[Dict[str, Any]] = []
    for i, panel in enumerate(panels_in):
        if not isinstance(panel, dict):
            continue
        p = {
            "key": str(panel.get("key") or panel.get("name") or f"series_{i}"),
            "title": panel.get("title"),
            "ylabel": panel.get("ylabel"),
            "line_color": panel.get("line_color")
            or panel.get("color")
            or DEFAULT_LINE_COLORS[i % len(DEFAULT_LINE_COLORS)],
            "xlim": _as_pair(panel.get("xlim")),
            "ylim": _as_pair(panel.get("ylim")),
            "show_grid": panel.get("show_grid"),
        }
        if panel.get("name"):
            p["name"] = panel["name"]
        # Multi-set compare panels list every series drawn on the subplot
        # (one property across sets, or one set across properties).
        series_keys = panel.get("series_keys")
        if isinstance(series_keys, list) and series_keys:
            p["series_keys"] = [str(k) for k in series_keys if k is not None and str(k)]
        if "show_xlabel" in panel:
            p["show_xlabel"] = bool(panel.get("show_xlabel"))
        if "show_ylabel" in panel:
            p["show_ylabel"] = bool(panel.get("show_ylabel"))
        if "show_ticks" in panel:
            p["show_ticks"] = bool(panel.get("show_ticks"))
        if "show_xticklabels" in panel:
            p["show_xticklabels"] = bool(panel.get("show_xticklabels"))
        if "show_yticklabels" in panel:
            p["show_yticklabels"] = bool(panel.get("show_yticklabels"))
        for spine_key in (
            "show_spine_left",
            "show_spine_bottom",
            "show_spine_top",
            "show_spine_right",
        ):
            if spine_key in panel:
                p[spine_key] = bool(panel.get(spine_key))
        for num_key in ("tick_length", "tick_width", "spine_width"):
            if panel.get(num_key) is not None:
                p[num_key] = panel.get(num_key)
        if panel.get("linewidth") is not None:
            p["linewidth"] = panel.get("linewidth")
        if panel.get("linestyle"):
            p["linestyle"] = str(panel.get("linestyle"))
        if "show_legend" in panel:
            p["show_legend"] = bool(panel.get("show_legend"))
        if panel.get("legend_loc"):
            p["legend_loc"] = str(panel.get("legend_loc"))
        if panel.get("legend_fontsize") is not None:
            p["legend_fontsize"] = _clamp_legend_fontsize(panel.get("legend_fontsize"))
        refs = normalize_reference_lines(panel.get("reference_lines"))
        if refs:
            p["reference_lines"] = refs
        panels.append(p)

    layout = str(src.get("layout") or "overlay").lower()
    if layout not in ("overlay", "grid"):
        layout = "overlay"

    try:
        cols = int(src.get("cols") or 2)
    except (TypeError, ValueError):
        cols = 2
    cols = max(1, min(cols, 8))

    try:
        rows_raw = src.get("rows")
        rows = int(rows_raw) if rows_raw not in (None, "") else None
    except (TypeError, ValueError):
        rows = None
    if rows is not None:
        rows = max(1, min(rows, 16))

    last_row_align = str(src.get("last_row_align") or "start").lower()
    if last_row_align not in ("start", "center", "end"):
        last_row_align = "start"

    wspace = _as_float(src.get("wspace"))
    hspace = _as_float(src.get("hspace"))
    cell_aspect = _as_float(src.get("cell_aspect"))

    ref_lines = normalize_reference_lines(src.get("reference_lines"))
    g_refs = normalize_reference_lines(global_style.get("reference_lines"))
    hlines = normalize_reference_lines(
        [{"axis": "y", "value": v} if not isinstance(v, dict) else {**v, "axis": "y"} for v in (global_style.get("hlines") or [])]
        if isinstance(global_style.get("hlines"), list)
        else []
    )
    vlines = normalize_reference_lines(
        [{"axis": "x", "value": v} if not isinstance(v, dict) else {**v, "axis": "x"} for v in (global_style.get("vlines") or [])]
        if isinstance(global_style.get("vlines"), list)
        else []
    )
    merged_refs = ref_lines + g_refs + hlines + vlines

    out: Dict[str, Any] = {
        "version": int(src.get("version") or PLOT_SPEC_VERSION),
        "layout": layout,
        "cols": cols,
        "sync_x": bool(src.get("sync_x", True)),
        "last_row_align": last_row_align,
        "legend": _normalize_legend(src.get("legend")),
        "global": global_style,
        "panels": panels,
    }
    if rows is not None:
        out["rows"] = rows
    if wspace is not None:
        out["wspace"] = wspace
    if hspace is not None:
        out["hspace"] = hspace
    if cell_aspect is not None and cell_aspect > 0:
        out["cell_aspect"] = cell_aspect
    if merged_refs:
        out["reference_lines"] = merged_refs
    return out


def build_plot_spec_from_series(
    series: Sequence[Dict[str, Any]],
    *,
    layout: str = "overlay",
    cols: int = 2,
    sync_x: bool = True,
    global_style: Optional[Dict[str, Any]] = None,
    panel_overrides: Optional[Dict[str, Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    """Build PlotSpec from energetic JSON ``{x, series: [{name, key, unit, y}]}``."""
    overrides = panel_overrides or {}
    panels: List[Dict[str, Any]] = []
    for i, item in enumerate(series):
        if not isinstance(item, dict):
            continue
        key = str(item.get("key") or item.get("name") or f"series_{i}")
        name = str(item.get("name") or key)
        unit = item.get("unit") or ""
        ylabel = f"{name} ({unit})" if unit else name
        base = {
            "key": key,
            "name": name,
            "title": name,
            "ylabel": ylabel,
            "line_color": DEFAULT_LINE_COLORS[i % len(DEFAULT_LINE_COLORS)],
            "xlim": None,
            "ylim": None,
            "show_grid": None,
        }
        if key in overrides:
            base.update({k: v for k, v in overrides[key].items() if v is not None})
        elif name in overrides:
            base.update({k: v for k, v in overrides[name].items() if v is not None})
        panels.append(base)

    if layout == "grid" and len(panels) == 1:
        layout = "overlay"

    g = {**DEFAULT_GLOBAL_STYLE, **(global_style or {})}
    if not g.get("xlabel"):
        g["xlabel"] = f"Time ({g.get('time_units', 'ns')})"

    return normalize_plot_spec(
        {
            "layout": layout,
            "cols": cols,
            "sync_x": sync_x,
            "global": g,
            "panels": panels,
        }
    )


def plot_spec_from_plot_properties_kwargs(
    properties: Sequence[str],
    *,
    separate_plots: bool = False,
    line_colors: Optional[Sequence[str]] = None,
    **style: Any,
) -> Dict[str, Any]:
    """Map legacy ``plot_properties`` kwargs to PlotSpec."""
    panels: List[Dict[str, Any]] = []
    colors = list(line_colors or DEFAULT_LINE_COLORS)
    for i, prop in enumerate(properties):
        panels.append(
            {
                "key": str(prop),
                "name": str(prop),
                "title": style.get("title") or str(prop),
                "ylabel": style.get("ylabel"),
                "line_color": colors[i % len(colors)],
                "xlim": _as_pair(style.get("xlim")),
                "ylim": _as_pair(style.get("ylim")),
            }
        )
    layout = "grid" if separate_plots else "overlay"
    if separate_plots and len(panels) == 1:
        layout = "overlay"
    global_style = {
        "time_units": style.get("time_units", "ns"),
        "energy_units": style.get("energy_units", "kcal/mol"),
        "pressure_units": style.get("pressure_units", "atm"),
        "temperature_units": style.get("temperature_units", "K"),
        "volume_units": style.get("volume_units", "Å³"),
        "plot_bg": style.get("bg_color", DEFAULT_GLOBAL_STYLE["plot_bg"]),
        "fig_bg": style.get("fig_bg_color", DEFAULT_GLOBAL_STYLE["fig_bg"]),
        "text_color": style.get("text_color", "Auto"),
        "grid_color": style.get("grid_color"),
        "show_grid": style.get("show_grid", True),
        "figsize": list(style.get("figsize") or DEFAULT_GLOBAL_STYLE["figsize"]),
        "dpi": style.get("dpi", 300),
        "font_family": style.get("font_family", "sans-serif"),
        "xlabel": style.get("xlabel"),
        "title": style.get("title"),
        "xlim": _as_pair(style.get("xlim")),
    }
    return normalize_plot_spec(
        {
            "layout": layout,
            "cols": 2,
            "sync_x": True,
            "global": global_style,
            "panels": panels,
        }
    )


def panel_effective_limits(
    spec: Dict[str, Any], panel: Dict[str, Any]
) -> Tuple[Optional[Tuple[float, float]], Optional[Tuple[float, float]]]:
    """Resolve xlim/ylim for one panel (panel override → global)."""
    spec = normalize_plot_spec(spec)
    g = spec["global"]
    xlim = _as_pair(panel.get("xlim")) or _as_pair(g.get("xlim"))
    ylim = _as_pair(panel.get("ylim")) or _as_pair(g.get("ylim"))
    return xlim, ylim


def union_axis_limits(
    pairs: Sequence[Optional[Tuple[float, float]]],
) -> Optional[Tuple[float, float]]:
    """Combine several (min, max) windows into one spanning window."""
    vals = [p for p in pairs if p is not None]
    if not vals:
        return None
    return (min(p[0] for p in vals), max(p[1] for p in vals))


def panel_show_grid(spec: Dict[str, Any], panel: Dict[str, Any]) -> bool:
    spec = normalize_plot_spec(spec)
    if panel.get("show_grid") is not None:
        return bool(panel.get("show_grid"))
    return bool(spec["global"].get("show_grid", True))
