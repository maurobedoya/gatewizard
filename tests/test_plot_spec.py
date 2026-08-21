"""Tests for PlotSpec and matplotlib energetic renderer."""

from __future__ import annotations

import pytest

from gatewizard.utils.plot_spec import (
    build_plot_spec_from_series,
    normalize_plot_spec,
    panel_effective_limits,
    plot_spec_from_plot_properties_kwargs,
)
from gatewizard.utils.matplotlib_renderer import render_energetic, render_energetic_to_bytes


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
