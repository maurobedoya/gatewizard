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
    "title": None,
    "xlim": None,
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
        panels.append(p)

    layout = str(src.get("layout") or "overlay").lower()
    if layout not in ("overlay", "grid"):
        layout = "overlay"

    cols = int(src.get("cols") or 2)
    cols = max(1, min(cols, 4))

    return {
        "version": int(src.get("version") or PLOT_SPEC_VERSION),
        "layout": layout,
        "cols": cols,
        "sync_x": bool(src.get("sync_x", True)),
        "global": global_style,
        "panels": panels,
    }


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
    ylim = _as_pair(panel.get("ylim"))
    return xlim, ylim


def panel_show_grid(spec: Dict[str, Any], panel: Dict[str, Any]) -> bool:
    spec = normalize_plot_spec(spec)
    if panel.get("show_grid") is not None:
        return bool(panel.get("show_grid"))
    return bool(spec["global"].get("show_grid", True))
