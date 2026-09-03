"""Tests for PlotSpec and matplotlib energetic renderer."""

from __future__ import annotations

import pytest

from gatewizard.utils.plot_spec import (
    build_plot_spec_from_series,
    grid_spec_slices,
    normalize_plot_spec,
    normalize_reference_lines,
    panel_effective_limits,
    plot_spec_from_plot_properties_kwargs,
)
from gatewizard.utils.matplotlib_renderer import (
    render_energetic,
    render_energetic_to_bytes,
    scaled_linestyle,
)


@pytest.fixture
def sample_data():
    return {
        "x": [0.0, 1.0, 2.0, 3.0],
        "series": [
            {"key": "temp", "name": "Temperature", "unit": "K", "y": [300.0, 301.0, 302.0, 303.0]},
            {"key": "press", "name": "Pressure", "unit": "atm", "y": [1.0, 1.01, 0.99, 1.0]},
        ],
    }


def test_normalize_plot_spec_defaults():
    spec = normalize_plot_spec({"layout": "grid", "panels": [{"key": "TEMP"}]})
    assert spec["layout"] == "grid"
    assert spec["cols"] == 2
    assert spec["panels"][0]["line_color"]


def test_build_plot_spec_from_series_overlay(sample_data):
    spec = build_plot_spec_from_series(sample_data["series"], layout="overlay")
    assert spec["layout"] == "overlay"
    assert len(spec["panels"]) == 2


def test_plot_spec_from_plot_properties_kwargs():
    spec = plot_spec_from_plot_properties_kwargs(
        ["Temperature", "Pressure"],
        separate_plots=True,
        time_units="ns",
        bg_color="#111111",
    )
    assert spec["layout"] == "grid"


@pytest.mark.parametrize("layout", ["overlay", "grid"])
def test_render_energetic_smoke(sample_data, layout):
    pytest.importorskip("matplotlib")
    spec = build_plot_spec_from_series(sample_data["series"], layout=layout)
    fig = render_energetic(sample_data, spec)
    assert fig is not None
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_render_energetic_to_bytes(sample_data):
    pytest.importorskip("matplotlib")
    spec = build_plot_spec_from_series(sample_data["series"], layout="overlay")
    png = render_energetic_to_bytes(sample_data, spec, fmt="png", dpi=100)
    assert png[:8] == b"\x89PNG\r\n\x1a\n"


def test_render_overlay_uses_union_ylim_across_series():
    """APL-style overlay: mean + leaflets must not clip to the first series ylim."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    data = {
        "x": [0.0, 1.0, 2.0, 3.0],
        "series": [
            {"key": "mean", "name": "Mean", "y": [64.0, 65.0, 64.5, 65.2], "x": [0.0, 1.0, 2.0, 3.0]},
            {
                "key": "upper",
                "name": "Upper leaflet",
                "y": [70.0, 78.0, 72.0, 80.0],
                "x": [0.0, 1.0, 2.0, 3.0],
            },
            {
                "key": "lower",
                "name": "Lower leaflet",
                "y": [50.0, 52.0, 48.0, 55.0],
                "x": [0.0, 1.0, 2.0, 3.0],
            },
        ],
    }
    spec = {
        "layout": "overlay",
        "global": {"title": "Area per Lipid", "ylabel": "Area (Å²)"},
        "panels": [
            {"key": "mean", "name": "Mean", "ylim": [63.0, 66.0]},
            {"key": "upper", "name": "Upper leaflet", "ylim": [69.0, 81.0]},
            {"key": "lower", "name": "Lower leaflet", "ylim": [47.0, 56.0]},
        ],
    }
    fig = render_energetic(data, spec)
    try:
        ax = fig.axes[0]
        assert len(ax.lines) == 3
        y0, y1 = ax.get_ylim()
        assert y0 == pytest.approx(47.0)
        assert y1 == pytest.approx(81.0)
    finally:
        plt.close(fig)


def test_panel_effective_limits_falls_back_to_global_ylim():
    spec = normalize_plot_spec(
        {
            "layout": "overlay",
            "global": {"ylim": [10.0, 90.0]},
            "panels": [{"key": "a"}],
        }
    )
    xlim, ylim = panel_effective_limits(spec, spec["panels"][0])
    assert ylim == (10.0, 90.0)


def test_normalize_preserves_series_keys():
    spec = normalize_plot_spec(
        {
            "layout": "grid",
            "panels": [
                {
                    "key": "set-a",
                    "title": "Set A",
                    "series_keys": ["set-a:TEMP", "set-a:TOTAL"],
                }
            ],
        }
    )
    assert spec["panels"][0]["series_keys"] == ["set-a:TEMP", "set-a:TOTAL"]


def test_render_grid_by_set_with_series_keys_and_sync_x():
    """GUI energetic Pub PNG (compare: one panel per set) uses series_keys + sync_x."""
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    data = {
        "x": [0.0, 1.0, 2.0, 3.0],
        "series": [
            {
                "key": "nvt:TOTAL",
                "name": "TOTAL",
                "unit": "kJ/mol",
                "x": [0.0, 1.0, 2.0, 3.0],
                "y": [-200000.0, -199000.0, -198500.0, -198000.0],
            },
            {
                "key": "nvt:TEMP",
                "name": "TEMP",
                "unit": "K",
                "x": [0.0, 1.0, 2.0, 3.0],
                "y": [300.0, 301.0, 299.0, 300.5],
            },
            {
                "key": "npt:TOTAL",
                "name": "TOTAL",
                "unit": "kJ/mol",
                "x": [0.0, 1.0, 2.0, 3.0],
                "y": [-180000.0, -179000.0, -178500.0, -178000.0],
            },
            {
                "key": "npt:TEMP",
                "name": "TEMP",
                "unit": "K",
                "x": [0.0, 1.0, 2.0, 3.0],
                "y": [300.0, 302.0, 301.0, 300.0],
            },
        ],
    }
    spec = {
        "layout": "grid",
        "cols": 2,
        "sync_x": True,
        "global": {
            "title": "OpenMM Energetic Analysis",
            "xlabel": "Time (ns)",
            "plot_bg": "#2b2b2b",
            "fig_bg": "#212121",
            "text_color": "#cccccc",
        },
        "panels": [
            {
                "key": "nvt",
                "title": "OpenMM - NVT",
                "ylabel": "Value",
                "series_keys": ["nvt:TOTAL", "nvt:TEMP"],
                "xlim": [0.0, 3.0],
                "ylim": [-210000.0, 500.0],
            },
            {
                "key": "npt",
                "title": "OpenMM - NPT",
                "ylabel": "Value",
                "series_keys": ["npt:TOTAL", "npt:TEMP"],
                "xlim": [0.0, 3.0],
                "ylim": [-190000.0, 500.0],
            },
        ],
    }
    fig = render_energetic(data, spec)
    try:
        plotted = [ax for ax in fig.axes if ax.get_visible() and ax.lines]
        assert len(plotted) == 2
        for ax in plotted:
            assert len(ax.lines) == 2
            x0, x1 = ax.get_xlim()
            # sync_x must use time (xlim), not energy (ylim)
            assert x0 == pytest.approx(0.0)
            assert x1 == pytest.approx(3.0)
            assert x1 - x0 < 100  # not the ±1e5 energy window
    finally:
        plt.close(fig)


def test_normalize_plot_spec_raises_cols_clamp_and_legend():
    spec = normalize_plot_spec(
        {
            "layout": "grid",
            "cols": 8,
            "rows": 3,
            "last_row_align": "center",
            "legend": {"mode": "outside", "loc": "bottom", "entries": "sets", "ncol": 2},
            "reference_lines": [{"axis": "y", "value": 60.0, "style": "dashed"}],
            "panels": [
                {
                    "key": "p0",
                    "series_keys": ["a:mean", "b:mean"],
                    "show_xlabel": False,
                    "show_xticklabels": False,
                    "show_yticklabels": False,
                    "linewidth": 3.5,
                    "linestyle": "dashed",
                    "show_legend": True,
                }
            ],
        }
    )
    assert spec["cols"] == 8
    assert spec["rows"] == 3
    assert spec["last_row_align"] == "center"
    assert spec["legend"]["mode"] == "outside"
    assert spec["legend"]["ncol"] == 2
    assert spec["panels"][0]["series_keys"] == ["a:mean", "b:mean"]
    assert spec["panels"][0]["show_xlabel"] is False
    assert spec["panels"][0]["show_xticklabels"] is False
    assert spec["panels"][0]["show_yticklabels"] is False
    assert spec["panels"][0]["linewidth"] == 3.5
    assert spec["panels"][0]["linestyle"] == "dashed"
    assert spec["reference_lines"][0]["value"] == 60.0


def test_normalize_legend_clamps_gui_pixel_fontsize():
    spec = normalize_plot_spec(
        {
            "layout": "grid",
            "legend": {"mode": "outside", "fontsize": 40},
            "panels": [{"key": "p0", "legend_fontsize": 30, "legend_loc": "upper left"}],
        }
    )
    assert spec["legend"]["fontsize"] == 11
    assert spec["panels"][0]["legend_fontsize"] == 11
    assert spec["panels"][0]["legend_loc"] == "upper left"


def test_grid_spec_slices_centers_8_in_3_cols():
    slices, rows, micro = grid_spec_slices(8, 3, "center")
    assert micro == 6
    assert rows == 3
    assert len(slices) == 8
    assert slices[6] == (2, 1, 3)
    assert slices[7] == (2, 3, 5)


def test_grid_spec_slices_end_aligns_7_in_4_cols():
    slices, rows, micro = grid_spec_slices(7, 4, "end")
    assert micro == 8
    assert rows == 2
    assert len(slices) == 7
    assert slices[4] == (1, 2, 4)
    assert slices[6] == (1, 6, 8)


def test_normalize_reference_lines_accepts_hlines():
    lines = normalize_reference_lines([{"axis": "x", "value": 10, "style": "dotted"}])
    assert lines[0]["axis"] == "x"
    assert lines[0]["style"] == "dotted"


def test_render_grid_grouped_series_keys_and_reference_line():
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    data = {
        "x": [0.0, 50.0, 100.0],
        "series": [
            {
                "key": "s1:mean",
                "name": "Set 1",
                "set_id": "s1",
                "set_name": "Set 1",
                "series_role": "mean",
                "x": [0.0, 50.0, 100.0],
                "y": [1.0, 1.1, 1.2],
                "color": "#f59e0b",
            },
            {
                "key": "s3:mean",
                "name": "Set 3",
                "set_id": "s3",
                "set_name": "Set 3",
                "series_role": "mean",
                "x": [0.0, 50.0, 100.0],
                "y": [2.0, 2.1, 2.2],
                "color": "#22c55e",
            },
            {
                "key": "s2:mean",
                "name": "Set 2",
                "set_id": "s2",
                "set_name": "Set 2",
                "series_role": "mean",
                "x": [0.0, 50.0, 100.0],
                "y": [0.5, 0.6, 0.7],
                "color": "#38bdf8",
            },
        ],
    }
    spec = {
        "layout": "grid",
        "cols": 2,
        "sync_x": False,
        "legend": {"mode": "each"},
        "reference_lines": [{"axis": "y", "value": 1.5, "style": "dashed", "color": "#888888"}],
        "global": {"xlabel": "Time (ns)", "ylabel": "RMSD (Å)", "title": "RMSD"},
        "panels": [
            {
                "key": "cell-0",
                "title": "1, 3, 2",
                "series_keys": ["s1:mean", "s3:mean", "s2:mean"],
                "show_legend": True,
            },
            {
                "key": "cell-1",
                "title": "Set 2",
                "series_keys": ["s2:mean"],
                "show_legend": False,
            },
        ],
    }
    fig = render_energetic(data, spec)
    try:
        plotted = [ax for ax in fig.axes if ax.get_visible() and ax.lines]
        assert len(plotted) == 2
        # cell 1: three series + one hline
        assert len(plotted[0].lines) == 4
        labels = [ln.get_label() for ln in plotted[0].lines if not ln.get_label().startswith("_")]
        assert labels[:3] == ["Set 1", "Set 3", "Set 2"]
        assert plotted[0].get_legend() is not None
        assert plotted[1].get_legend() is None
    finally:
        plt.close(fig)


def test_render_8_panels_3_cols_center_last_row():
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    series = []
    panels = []
    for i in range(8):
        key = f"s{i}"
        series.append(
            {
                "key": key,
                "name": key,
                "x": [0.0, 1.0],
                "y": [float(i), float(i) + 0.5],
            }
        )
        panels.append({"key": key, "title": key, "series_keys": [key]})
    data = {"x": [0.0, 1.0], "series": series}
    spec = {
        "layout": "grid",
        "cols": 3,
        "last_row_align": "center",
        "legend": {"mode": "none"},
        "global": {"xlabel": "t", "title": "grid"},
        "panels": panels,
    }
    fig = render_energetic(data, spec)
    try:
        plotted = [ax for ax in fig.axes if ax.get_visible() and ax.lines]
        assert len(plotted) == 8
        # leftover 9th slot must not exist as a hidden empty axes
        assert len(plotted) == 8
        xs = [ax.get_position().x0 for ax in plotted]
        # last two panels sit inward of the first-column x0
        first_col = min(xs[:3])
        assert plotted[6].get_position().x0 > first_col + 0.02
        assert plotted[7].get_position().x0 > plotted[6].get_position().x0
    finally:
        plt.close(fig)


def test_render_outside_legend_and_one_cell_mode():
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    data = {
        "x": [0.0, 1.0],
        "series": [
            {
                "key": "a",
                "name": "A",
                "set_id": "a",
                "set_name": "A",
                "series_role": "mean",
                "x": [0.0, 1.0],
                "y": [1.0, 2.0],
            },
            {
                "key": "b",
                "name": "B",
                "set_id": "b",
                "set_name": "B",
                "series_role": "mean",
                "x": [0.0, 1.0],
                "y": [2.0, 3.0],
            },
        ],
    }
    spec = {
        "layout": "grid",
        "cols": 2,
        "legend": {"mode": "outside", "loc": "bottom", "entries": "sets", "title": "Sets"},
        "global": {"xlabel": "t"},
        "panels": [
            {"key": "p0", "title": "P0", "series_keys": ["a"]},
            {"key": "p1", "title": "P1", "series_keys": ["b"]},
        ],
    }
    fig = render_energetic(data, spec)
    try:
        plot_axes = [ax for ax in fig.axes if ax.lines]
        assert len(plot_axes) == 2
        assert plot_axes[0].get_legend() is None
        assert plot_axes[1].get_legend() is None
        fig_legs = [c for c in fig.axes if c.get_legend() is not None]
        assert len(fig_legs) >= 1
    finally:
        plt.close(fig)

    spec["legend"] = {"mode": "one", "cell": 1}
    fig = render_energetic(data, spec)
    try:
        plot_axes = [ax for ax in fig.axes if ax.lines]
        assert plot_axes[0].get_legend() is None
        assert plot_axes[1].get_legend() is not None
    finally:
        plt.close(fig)


def test_scaled_linestyle_grows_with_linewidth():
    assert scaled_linestyle("solid", 3.0) == "-"
    thin = scaled_linestyle("dashed", 1.0)
    thick = scaled_linestyle("dashed", 3.0)
    assert thin[0] == 0
    assert thick[1][0] == pytest.approx(thin[1][0] * 3)
    assert thick[1][1] == pytest.approx(thin[1][1] * 3)
    dotted = scaled_linestyle("dotted", 2.0)
    assert dotted[1][1] > dotted[1][0]


def test_normalize_plot_spec_keeps_axis_chrome():
    spec = normalize_plot_spec(
        {
            "layout": "overlay",
            "global": {
                "tick_length": 12,
                "tick_width": 2.5,
                "spine_width": 0.6,
                "show_spine_top": True,
                "show_spine_right": True,
                "show_spine_left": False,
            },
            "panels": [
                {
                    "key": "temp",
                    "show_spine_bottom": False,
                    "tick_width": 3.0,
                }
            ],
        }
    )
    g = spec["global"]
    assert g["show_spine_bottom"] is True
    assert g["show_spine_top"] is True
    assert g["show_spine_right"] is True
    assert g["show_spine_left"] is False
    assert g["tick_length"] == 12
    assert g["tick_width"] == 2.5
    assert spec["panels"][0]["show_spine_bottom"] is False
    assert spec["panels"][0]["tick_width"] == 3.0


def test_render_default_hides_top_right_spines(sample_data):
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    spec = normalize_plot_spec(
        {"layout": "overlay", "panels": [{"key": "temp", "name": "Temperature"}]}
    )
    fig = render_energetic(sample_data, spec)
    try:
        ax = fig.axes[0]
        assert ax.spines["left"].get_visible() is True
        assert ax.spines["bottom"].get_visible() is True
        assert ax.spines["top"].get_visible() is False
        assert ax.spines["right"].get_visible() is False
        fig.canvas.draw()
        labels = [lbl for lbl in ax.get_xticklabels() if lbl.get_text()]
        assert labels
        assert all(lbl.get_ha() == "center" for lbl in labels)
    finally:
        plt.close(fig)


def test_render_honors_spine_and_tick_chrome(sample_data):
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    spec = normalize_plot_spec(
        {
            "layout": "overlay",
            "global": {
                "show_spine_left": False,
                "show_spine_bottom": True,
                "show_spine_top": True,
                "show_spine_right": True,
                "spine_width": 2.5,
                "tick_width": 2.0,
                "tick_length": 8.0,
                "show_ticks": True,
            },
            "panels": [{"key": "temp", "name": "Temperature"}],
        }
    )
    fig = render_energetic(sample_data, spec)
    try:
        ax = fig.axes[0]
        assert ax.spines["left"].get_visible() is False
        assert ax.spines["top"].get_visible() is True
        assert ax.spines["right"].get_visible() is True
        assert ax.spines["bottom"].get_linewidth() == pytest.approx(2.5)
        tick = ax.xaxis.majorTicks[0]
        assert tick.tick1line.get_markeredgewidth() == pytest.approx(2.0)
        assert tick.tick1line.get_markersize() == pytest.approx(8.0)
    finally:
        plt.close(fig)


def test_normalize_plot_spec_extra_margins_default_to_zero():
    spec = normalize_plot_spec({"layout": "overlay", "panels": [{"key": "temp"}]})
    g = spec["global"]
    assert g["extra_left"] == 0
    assert g["extra_right"] == 0
    assert g["extra_top"] == 0
    assert g["extra_bottom"] == 0
    spec2 = normalize_plot_spec(
        {
            "layout": "overlay",
            "global": {"extra_left": -15, "extra_top": 8},
            "panels": [{"key": "temp"}],
        }
    )
    assert spec2["global"]["extra_left"] == -15
    assert spec2["global"]["extra_top"] == 8


def test_render_extra_left_margin_shifts_overlay(sample_data):
    pytest.importorskip("matplotlib")
    import matplotlib.pyplot as plt

    base = normalize_plot_spec(
        {"layout": "overlay", "panels": [{"key": "temp", "name": "Temperature"}]}
    )
    padded = normalize_plot_spec(
        {
            "layout": "overlay",
            "global": {"extra_left": 40},
            "panels": [{"key": "temp", "name": "Temperature"}],
        }
    )
    fig0 = render_energetic(sample_data, base)
    fig1 = render_energetic(sample_data, padded)
    try:
        assert fig1.axes[0].get_position().x0 > fig0.axes[0].get_position().x0 + 0.01
    finally:
        plt.close(fig0)
        plt.close(fig1)
