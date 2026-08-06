"""Tests for PlotSpec and matplotlib energetic renderer."""

from __future__ import annotations

import pytest

from gatewizard.utils.plot_spec import (
    build_plot_spec_from_series,
    normalize_plot_spec,
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
