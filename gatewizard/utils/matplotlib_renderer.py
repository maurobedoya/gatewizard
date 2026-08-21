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


def _style_axes(ax, *, text_color: str, grid_color: str, show_grid: bool, bg_color: str) -> None:
    ax.tick_params(colors=text_color)
    for spine in ax.spines.values():
        spine.set_color(text_color)
    if bg_color != "none":
        ax.set_facecolor(bg_color)
    if show_grid:
        ax.grid(True, alpha=0.3, color=grid_color, linewidth=0.5)


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
                linewidth=1.5,
                marker="o" if use_marker else None,
                markersize=2,
                label=label,
            )
            lines.extend(line)
            labels.append(label)
        _style_axes(
            ax,
            text_color=text_color,
            grid_color=grid_color,
            show_grid=bool(g.get("show_grid", True)),
            bg_color=bg_color,
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
        plt.tight_layout()
        return fig

    # grid layout — one property per panel
    cols = int(spec.get("cols") or 2)
    n = len(panels)
    rows = (n + cols - 1) // cols
    fig, axes_arr = plt.subplots(rows, cols, figsize=(figsize[0], figsize[1] * rows * 0.55), squeeze=False)
    if fig_bg != "none":
        fig.patch.set_facecolor(fig_bg)

    shared_xlim = None
    if spec.get("sync_x"):
        shared_xlim, _ = panel_effective_limits(spec, panels[0])
        if shared_xlim is None:
            shared_xlim = tuple(g.get("xlim") or ()) if g.get("xlim") else None

    for i, panel in enumerate(panels):
        r, c = divmod(i, cols)
        ax = axes_arr[r][c]
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
            ax.plot(xs[:npts], ys[:npts], color=color, linewidth=1.5, label=label)
            legend_labels.append(label)
        _style_axes(
            ax,
            text_color=text_color,
            grid_color=grid_color,
            show_grid=panel_show_grid(spec, panel),
            bg_color=bg_color,
        )
        name = str(panel.get("title") or panel.get("name") or panel.get("key"))
        unit = (series_list[0] if series_list else {}).get("unit") or ""
        ax.set_xlabel(g.get("xlabel") or f"Time ({g.get('time_units', 'ns')})", color=text_color)
        ax.set_ylabel(panel.get("ylabel") or (f"{name} ({unit})" if unit else name), color=text_color)
        ax.set_title(name, color=text_color, fontweight="bold")
        if len(legend_labels) > 1:
            legend = ax.legend(fontsize=8, loc="best")
            if legend:
                plt.setp(legend.get_texts(), color=text_color)
        xlim, ylim = panel_effective_limits(spec, panel)
        if spec.get("sync_x") and shared_xlim:
            ax.set_xlim(shared_xlim)
        elif xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

    for j in range(len(panels), rows * cols):
        r, c = divmod(j, cols)
        axes_arr[r][c].set_visible(False)

    if g.get("title"):
        fig.suptitle(g["title"], color=text_color, fontweight="bold")
    plt.tight_layout()
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
