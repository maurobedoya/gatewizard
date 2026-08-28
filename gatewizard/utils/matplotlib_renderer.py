"""
Matplotlib renderer for energetic PlotSpec (shared by API and GUI export).
"""

from __future__ import annotations

import io
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from gatewizard.utils.plot_spec import (
    DEFAULT_LINE_COLORS,
    grid_spec_slices,
    normalize_plot_spec,
    panel_effective_limits,
    panel_show_grid,
    union_axis_limits,
)

logger = logging.getLogger(__name__)


def _configure_matplotlib_for_headless() -> None:
    """Use Agg so FastAPI/CLI export never touches Tk (no main loop in worker threads)."""
    import matplotlib

    if matplotlib.get_backend().lower() != "agg":
        matplotlib.use("Agg", force=True)


def _pyplot():
    _configure_matplotlib_for_headless()
    import matplotlib.pyplot as plt

    return plt


def _auto_text_color(bg_color: str, text_color: str) -> str:
    if text_color and text_color != "Auto":
        return text_color
    if bg_color == "none":
        return "black"
    try:
        hex_color = str(bg_color).lstrip("#")
        if len(hex_color) != 6:
            return "white"
        r, g, b = (
            int(hex_color[0:2], 16),
            int(hex_color[2:4], 16),
            int(hex_color[4:6], 16),
        )
        luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
        return "black" if luminance > 0.5 else "white"
    except Exception:
        return "white"


_LEGEND_LOCS = {
    "upper left": "upper left",
    "upper right": "upper right",
    "lower left": "lower left",
    "lower right": "lower right",
    "best": "best",
    "top-left": "upper left",
    "top-right": "upper right",
    "bottom-left": "lower left",
    "bottom-right": "lower right",
}


def _legend_fontsize(raw, default: float = 8.0) -> float:
    try:
        n = float(raw)
    except (TypeError, ValueError):
        n = default
    if n > 14:
        n = max(6.0, min(11.0, round(n * 0.28 + 5.2)))
    return max(6.0, min(11.0, n))


def _legend_loc(raw) -> str:
    key = str(raw or "best").replace("_", " ").strip().lower()
    return _LEGEND_LOCS.get(key, "best")


def _style_axes(
    ax,
    *,
    text_color: str,
    grid_color: str,
    show_grid: bool,
    bg_color: str,
    spec: Optional[Dict[str, Any]] = None,
    panel: Optional[Dict[str, Any]] = None,
) -> None:
    g = (spec or {}).get("global") or {}
    p = panel or {}

    def _flag(key: str, default: bool) -> bool:
        if key in p and p.get(key) is not None:
            return bool(p.get(key))
        if key in g and g.get(key) is not None:
            return bool(g.get(key))
        return default

    def _num(key: str, default: float, lo: float, hi: float) -> float:
        raw = p.get(key)
        if raw is None:
            raw = g.get(key)
        try:
            n = float(raw)
        except (TypeError, ValueError):
            n = default
        return max(lo, min(hi, n))

    spine_w = _num("spine_width", 1.0, 0.2, 8.0)
    for side, default in (
        ("left", True),
        ("bottom", True),
        ("top", False),
        ("right", False),
    ):
        vis = _flag(f"show_spine_{side}", default)
        sp = ax.spines[side]
        sp.set_visible(vis)
        sp.set_color(text_color)
        sp.set_linewidth(spine_w)

    tick_len = _num("tick_length", 4.0, 0.0, 16.0)
    tick_w = _num("tick_width", 1.0, 0.2, 8.0)
    show_ticks = _flag("show_ticks", True)
    ax.tick_params(
        colors=text_color,
        length=0 if not show_ticks else tick_len,
        width=tick_w,
    )
    for lbl in ax.get_xticklabels():
        lbl.set_horizontalalignment("center")
    if bg_color != "none":
        ax.set_facecolor(bg_color)
    if show_grid:
        ax.grid(True, alpha=0.3, color=grid_color, linewidth=0.5)


def _apply_extra_figure_margins(fig, g: Dict[str, Any]) -> None:
    """Shift the subplot box by GUI extra_* values (CSS px ≈ 1/96 in)."""
    w, h = fig.get_size_inches()

    def _frac(key: str, dim: float) -> float:
        try:
            n = float(g.get(key) or 0)
        except (TypeError, ValueError):
            n = 0.0
        if dim <= 0:
            return 0.0
        return max(-0.25, min(0.4, (n / 96.0) / dim))

    dl = _frac("extra_left", w)
    dr = _frac("extra_right", w)
    dt = _frac("extra_top", h)
    db = _frac("extra_bottom", h)
    if dl == dr == dt == db == 0:
        return
    sp = fig.subplotpars
    left = min(max(0.02, sp.left + dl), 0.85)
    right = max(min(0.98, sp.right - dr), left + 0.1)
    bottom = min(max(0.02, sp.bottom + db), 0.85)
    top = max(min(0.98, sp.top - dt), bottom + 0.1)
    fig.subplots_adjust(left=left, right=right, top=top, bottom=bottom)


def _series_lookup(data: Dict[str, Any]) -> Dict[str, Dict[str, Any]]:
    out: Dict[str, Dict[str, Any]] = {}
    x = data.get("x") or []
    for item in data.get("series") or []:
        if not isinstance(item, dict):
            continue
        key = str(item.get("key") or item.get("name") or "")
        name = str(item.get("name") or key)
        out[key] = item
        out[name] = item
        out[key.lower()] = item
        out[name.lower()] = item
        item.setdefault("x", x)
    return out


def _panel_series(panel: Dict[str, Any], lookup: Dict[str, Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    key = str(panel.get("key") or "")
    for candidate in (key, panel.get("name"), key.lower(), str(panel.get("name") or "").lower()):
        if candidate and candidate in lookup:
            return lookup[candidate]
    return None


def _panel_series_list(
    panel: Dict[str, Any], lookup: Dict[str, Dict[str, Any]]
) -> List[Dict[str, Any]]:
    """All series for a grid panel (multi-set compare on one property or one set)."""
    keys = panel.get("series_keys")
    if isinstance(keys, list) and keys:
        out: List[Dict[str, Any]] = []
        for raw_key in keys:
            key = str(raw_key or "")
            for candidate in (key, key.lower()):
                if candidate and candidate in lookup:
                    out.append(lookup[candidate])
                    break
        if out:
            return out
    single = _panel_series(panel, lookup)
    return [single] if single else []


_LINESTYLES = {
    "solid": "-",
    "-": "-",
    "dashed": "--",
    "--": "--",
    "dotted": ":",
    ":": ":",
    "dashdot": "-.",
    "-.": "-.",
}


def _series_line_color(
    series: Dict[str, Any],
    panel: Dict[str, Any],
    index: int,
) -> str:
    color = series.get("color") or series.get("line_color")
    if color:
        return str(color)
    panel_color = panel.get("line_color")
    if panel_color:
        return str(panel_color)
    return DEFAULT_LINE_COLORS[index % len(DEFAULT_LINE_COLORS)]


def _series_linewidth(series: Dict[str, Any], panel: Dict[str, Any]) -> float:
    raw = series.get("linewidth")
    if raw is None:
        raw = series.get("line_width")
    if raw is None:
        raw = panel.get("linewidth") or panel.get("line_width")
    try:
        n = float(raw)
    except (TypeError, ValueError):
        n = 1.5
    return max(0.4, min(12.0, n))


_DASH_UNITS = {
    "dashed": (4.0, 2.5),
    "--": (4.0, 2.5),
    "dotted": (1.0, 2.0),
    ":": (1.0, 2.0),
    "dashdot": (4.0, 1.8, 1.0, 1.8),
    "-.": (4.0, 1.8, 1.0, 1.8),
}


def _style_name(raw: Any) -> str:
    return str(raw or "solid").lower()


def scaled_linestyle(style: str, linewidth: float):
    """Named dash pattern whose on/off lengths grow with ``linewidth`` (points)."""
    key = _style_name(style)
    if key in ("solid", "-", ""):
        return "-"
    pattern = _DASH_UNITS.get(key)
    if not pattern:
        return _LINESTYLES.get(key, "-")
    try:
        lw = float(linewidth)
    except (TypeError, ValueError):
        lw = 1.5
    lw = max(0.4, lw)
    return (0, tuple(round(part * lw, 3) for part in pattern))


def _series_style_name(series: Dict[str, Any], panel: Dict[str, Any]) -> str:
    return _style_name(
        series.get("linestyle")
        or series.get("line_style")
        or panel.get("linestyle")
        or panel.get("line_style")
        or "solid"
    )


def _line_cap_kwargs(linestyle) -> Dict[str, str]:
    if linestyle == "-" or linestyle is None:
        return {}
    return {"dash_capstyle": "butt", "solid_capstyle": "butt"}


def _series_plot_style(series: Dict[str, Any], panel: Dict[str, Any]) -> Dict[str, Any]:
    lw = _series_linewidth(series, panel)
    ls = scaled_linestyle(_series_style_name(series, panel), lw)
    out: Dict[str, Any] = {"linewidth": lw, "linestyle": ls}
    out.update(_line_cap_kwargs(ls))
    return out


def _draw_reference_lines(ax, lines: Optional[Sequence[Dict[str, Any]]], text_color: str) -> None:
    if not lines:
        return
    for line in lines:
        if not isinstance(line, dict):
            continue
        try:
            value = float(line.get("value"))
        except (TypeError, ValueError):
            continue
        color = str(line.get("color") or "#888888")
        try:
            width = float(line.get("width") or 1.2)
        except (TypeError, ValueError):
            width = 1.2
        ls = scaled_linestyle(str(line.get("style") or "dashed"), width)
        caps = _line_cap_kwargs(ls)
        axis = str(line.get("axis") or "y").lower()
        label = str(line.get("label") or "")
        if axis == "x":
            ax.axvline(value, color=color, linewidth=width, linestyle=ls, zorder=2, **caps)
            if label:
                ax.text(
                    value,
                    0.98,
                    label,
                    transform=ax.get_xaxis_transform(),
                    color=text_color,
                    fontsize=8,
                    ha="left",
                    va="top",
                )
        else:
            ax.axhline(value, color=color, linewidth=width, linestyle=ls, zorder=2, **caps)
            if label:
                ax.text(
                    0.02,
                    value,
                    label,
                    transform=ax.get_yaxis_transform(),
                    color=text_color,
                    fontsize=8,
                    ha="left",
                    va="bottom",
                )


def _panel_show_legend(spec: Dict[str, Any], panel: Dict[str, Any], index: int, n_labels: int) -> bool:
    legend = spec.get("legend") or {}
    mode = str(legend.get("mode") or "")
    if mode in ("none", "outside"):
        return False
    if "show_legend" in panel:
        show = bool(panel.get("show_legend"))
        if mode == "one":
            return show and index == int(legend.get("cell") or 0)
        return show
    if mode == "one":
        return index == int(legend.get("cell") or 0) and n_labels > 0
    return n_labels > 1


def _figure_legend_entries(
    axes,
    panels: Sequence[Dict[str, Any]],
    lookup: Dict[str, Dict[str, Any]],
    entries: str,
):
    handles = []
    labels = []
    seen = set()

    def add(handle, label, key):
        if key in seen or handle is None:
            return
        seen.add(key)
        handles.append(handle)
        labels.append(label)

    def walk(kind: str) -> None:
        for ax, panel in zip(axes, panels):
            hlist, llist = ax.get_legend_handles_labels()
            series_list = _panel_series_list(panel, lookup)
            for i, (hi, li) in enumerate(zip(hlist, llist)):
                series = series_list[i] if i < len(series_list) else {}
                if kind == "sets":
                    key = str(series.get("set_id") or series.get("set_name") or li)
                    lab = str(series.get("set_name") or li)
                    add(hi, lab, f"set:{key}")
                else:
                    role = str(series.get("series_role") or "mean")
                    lab = {"mean": "Mean", "upper": "Upper", "lower": "Lower"}.get(role, li)
                    add(hi, lab, f"role:{role}")

    if entries == "roles":
        walk("roles")
    elif entries == "both":
        walk("sets")
        walk("roles")
    else:
        walk("sets")
    return handles, labels


def render_energetic(
    data: Dict[str, Any],
    spec: Dict[str, Any],
):
    """Render energetic data with PlotSpec. Returns matplotlib Figure."""
    try:
        plt = _pyplot()
        import numpy as np
    except ImportError as exc:
        raise ImportError("matplotlib and numpy are required for plotting") from exc

    spec = normalize_plot_spec(spec)
    g = spec["global"]
    lookup = _series_lookup(data)
    panels = spec["panels"]
    if not panels:
        raise ValueError("PlotSpec has no panels")

    text_color = _auto_text_color(g.get("plot_bg", "#2b2b2b"), g.get("text_color", "Auto"))
    grid_color = g.get("grid_color") or text_color
    figsize = tuple(g.get("figsize") or (10, 6))
    bg_color = g.get("plot_bg", "#2b2b2b")
    fig_bg = g.get("fig_bg", "#212121")

    layout = spec["layout"]
    if layout == "overlay" or len(panels) == 1:
        fig, ax = plt.subplots(figsize=figsize)
        if fig_bg != "none":
            fig.patch.set_facecolor(fig_bg)
        lines = []
        labels = []
        for i, panel in enumerate(panels):
            series = _panel_series(panel, lookup)
            if not series:
                logger.warning("Panel %s: no matching series", panel.get("key"))
                continue
            xs = np.asarray(series.get("x") or data.get("x") or [], dtype=float)
            ys = np.asarray(series.get("y") or [], dtype=float)
            n = min(len(xs), len(ys))
            if n == 0:
                continue
            color = _series_line_color(series, panel, i)
            name = str(series.get("name") or panel.get("name") or panel.get("key"))
            unit = series.get("unit") or ""
            label = f"{name} ({unit})" if unit else name
            use_marker = len(panels) > 1 and n <= 80
            line = ax.plot(
                xs[:n],
                ys[:n],
                color=color,
                marker="o" if use_marker else None,
                markersize=2,
                label=label,
                **_series_plot_style(series, panel),
            )
            lines.extend(line)
            labels.append(label)
        _style_axes(
            ax,
            text_color=text_color,
            grid_color=grid_color,
            show_grid=bool(g.get("show_grid", True)),
            bg_color=bg_color,
            spec=spec,
        )
        ax.set_xlabel(g.get("xlabel") or f"Time ({g.get('time_units', 'ns')})", color=text_color)
        ylabels = {
            str(p.get("ylabel") or "").strip()
            for p in panels
            if str(p.get("ylabel") or "").strip()
        }
        if len(ylabels) == 1:
            ylabel = next(iter(ylabels))
        elif len(panels) == 1:
            ylabel = panels[0].get("ylabel")
        else:
            ylabel = g.get("ylabel") or "Value"
        ax.set_ylabel(ylabel or "Value", color=text_color)
        ax.set_title(g.get("title") or "", color=text_color, fontweight="bold")
        if len(labels) > 1:
            legend = ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
            if legend:
                plt.setp(legend.get_texts(), color=text_color)
        # Overlay shares one axis — span all panel windows so secondary series
        # (e.g. APL upper/lower leaflets) are not clipped to the first series.
        xlim = union_axis_limits(panel_effective_limits(spec, p)[0] for p in panels)
        ylim = union_axis_limits(panel_effective_limits(spec, p)[1] for p in panels)
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        _draw_reference_lines(ax, spec.get("reference_lines"), text_color)
        plt.tight_layout()
        _apply_extra_figure_margins(fig, g)
        return fig

    # grid layout — GridSpec so a short last row can be centered (no leftover hidden axes)
    cols = int(spec.get("cols") or 2)
    n = len(panels)
    last_row_align = str(spec.get("last_row_align") or "start")
    slices, rows, micro = grid_spec_slices(n, cols, last_row_align)
    legend_cfg = spec.get("legend") or {}
    legend_mode = str(legend_cfg.get("mode") or "")
    legend_loc = str(legend_cfg.get("loc") or "bottom")
    outside = legend_mode == "outside"
    extra_row = 1 if outside and legend_loc in ("top", "bottom") else 0
    extra_col = 1 if outside and legend_loc in ("left", "right") else 0

    fig_w, fig_h = figsize[0], figsize[1] * max(rows, 1) * 0.55
    if spec.get("cell_aspect") and figsize[1]:
        fig_w, fig_h = figsize[0], figsize[1]

    wspace = spec.get("wspace")
    hspace = spec.get("hspace")
    if wspace is None:
        wspace = 0.25
    if hspace is None:
        hspace = 0.35

    gs_rows = rows + extra_row
    gs_cols = micro + extra_col
    height_ratios = None
    width_ratios = None
    panel_row_offset = 0
    panel_col_offset = 0
    legend_row = None
    legend_col = None
    if extra_row:
        strip = 0.22
        if legend_loc == "bottom":
            height_ratios = [1] * rows + [strip]
            legend_row = rows
        else:
            height_ratios = [strip] + [1] * rows
            panel_row_offset = 1
            legend_row = 0
    if extra_col:
        strip = 0.22
        if legend_loc == "right":
            width_ratios = [1] * micro + [strip]
            legend_col = micro
        else:
            width_ratios = [strip] + [1] * micro
            panel_col_offset = 1
            legend_col = 0

    fig = plt.figure(figsize=(fig_w, fig_h))
    if fig_bg != "none":
        fig.patch.set_facecolor(fig_bg)
    gs = fig.add_gridspec(
        gs_rows,
        gs_cols,
        wspace=float(wspace),
        hspace=float(hspace),
        height_ratios=height_ratios,
        width_ratios=width_ratios,
    )

    shared_xlim = None
    if spec.get("sync_x"):
        shared_xlim, _ = panel_effective_limits(spec, panels[0])
        if shared_xlim is None:
            shared_xlim = tuple(g.get("xlim") or ()) if g.get("xlim") else None

    axes = []
    for i, panel in enumerate(panels):
        r, c0, c1 = slices[i]
        ax = fig.add_subplot(
            gs[
                r + panel_row_offset,
                (c0 + panel_col_offset) : (c1 + panel_col_offset),
            ]
        )
        axes.append(ax)
        series_list = _panel_series_list(panel, lookup)
        legend_labels: List[str] = []
        for j, series in enumerate(series_list):
            xs = np.asarray(series.get("x") or data.get("x") or [], dtype=float)
            ys = np.asarray(series.get("y") or [], dtype=float)
            npts = min(len(xs), len(ys))
            if not npts:
                continue
            color = _series_line_color(series, panel, j)
            name = str(series.get("name") or panel.get("name") or panel.get("key"))
            unit = series.get("unit") or ""
            label = f"{name} ({unit})" if unit else name
            ax.plot(
                xs[:npts],
                ys[:npts],
                color=color,
                label=label,
                **_series_plot_style(series, panel),
            )
            legend_labels.append(label)
        _style_axes(
            ax,
            text_color=text_color,
            grid_color=grid_color,
            show_grid=panel_show_grid(spec, panel),
            bg_color=bg_color,
            spec=spec,
            panel=panel,
        )
        name = str(panel.get("title") or panel.get("name") or panel.get("key"))
        unit = (series_list[0] if series_list else {}).get("unit") or ""
        show_xlabel = panel.get("show_xlabel", True)
        show_ylabel = panel.get("show_ylabel", True)
        show_ticks = panel.get("show_ticks", True)
        if show_xlabel:
            ax.set_xlabel(g.get("xlabel") or f"Time ({g.get('time_units', 'ns')})", color=text_color)
        else:
            ax.set_xlabel("")
        if show_ylabel:
            ax.set_ylabel(panel.get("ylabel") or (f"{name} ({unit})" if unit else name), color=text_color)
        else:
            ax.set_ylabel("")
        ax.set_title(name, color=text_color, fontweight="bold")
        if not show_ticks:
            ax.tick_params(length=0)
        if panel.get("show_xticklabels") is False:
            ax.tick_params(labelbottom=False)
        if panel.get("show_yticklabels") is False:
            ax.tick_params(labelleft=False)
        if _panel_show_legend(spec, panel, i, len(legend_labels)):
            legend = ax.legend(
                fontsize=_legend_fontsize(
                    panel.get("legend_fontsize") or legend_cfg.get("fontsize") or 8
                ),
                loc=_legend_loc(panel.get("legend_loc") or "best"),
            )
            if legend:
                plt.setp(legend.get_texts(), color=text_color)
        xlim, ylim = panel_effective_limits(spec, panel)
        if spec.get("sync_x") and shared_xlim:
            ax.set_xlim(shared_xlim)
        elif xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        refs = list(spec.get("reference_lines") or []) + list(panel.get("reference_lines") or [])
        _draw_reference_lines(ax, refs, text_color)

    if outside and axes:
        handles, labels = _figure_legend_entries(
            axes, panels, lookup, str(legend_cfg.get("entries") or "sets")
        )
        if handles:
            if extra_row and legend_row is not None:
                lax = fig.add_subplot(gs[legend_row, :])
            else:
                lax = fig.add_subplot(gs[:, legend_col])
            lax.axis("off")
            leg = lax.legend(
                handles,
                labels,
                loc="center",
                ncol=int(legend_cfg.get("ncol") or 1),
                fontsize=_legend_fontsize(legend_cfg.get("fontsize") or 8),
                title=legend_cfg.get("title") or None,
                frameon=False,
            )
            if leg:
                plt.setp(leg.get_texts(), color=text_color)
                if leg.get_title():
                    plt.setp(leg.get_title(), color=text_color)

    if g.get("title"):
        fig.suptitle(g["title"], color=text_color, fontweight="bold")
    # GridSpec already applies wspace/hspace; tight_layout fights a centered last row.
    _apply_extra_figure_margins(fig, g)
    return fig


def render_energetic_to_bytes(
    data: Dict[str, Any],
    spec: Dict[str, Any],
    *,
    fmt: str = "png",
    dpi: Optional[int] = None,
) -> bytes:
    """Render PlotSpec to PNG/SVG bytes."""
    plt = _pyplot()

    spec = normalize_plot_spec(spec)
    fig = render_energetic(data, spec)
    buf = io.BytesIO()
    try:
        fig.savefig(
            buf,
            format=fmt,
            dpi=dpi or int(spec["global"].get("dpi") or 300),
            bbox_inches="tight",
            facecolor=fig.get_facecolor(),
        )
        return buf.getvalue()
    finally:
        plt.close(fig)


def render_energetic_to_path(
    data: Dict[str, Any],
    spec: Dict[str, Any],
    path: Union[str, Path],
    *,
    fmt: Optional[str] = None,
) -> Path:
    """Save rendered figure to disk. Grid + separate_plots saves one file per panel when path is prefix."""
    plt = _pyplot()

    spec = normalize_plot_spec(spec)
    path = Path(path)
    layout = spec["layout"]
    panels = spec["panels"]
    dpi = int(spec["global"].get("dpi") or 300)

    if layout == "grid" and len(panels) > 1 and not path.suffix:
        prefix = str(path)
        if not prefix.endswith("_") and not prefix.endswith("/"):
            prefix = prefix + "_"
        for i, panel in enumerate(panels):
            single = normalize_plot_spec({**spec, "layout": "overlay", "panels": [panel]})
            fig = render_energetic(data, single)
            safe = str(panel.get("key") or panel.get("name") or i).lower().replace(" ", "_")
            out = Path(f"{prefix}{safe}.png")
            fig.savefig(out, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            logger.info("Plot saved: %s", out)
        return path

    fig = render_energetic(data, spec)
    out_path = path if path.suffix else path.with_suffix(".png")
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    logger.info("Plot saved: %s", out_path)
    return out_path
