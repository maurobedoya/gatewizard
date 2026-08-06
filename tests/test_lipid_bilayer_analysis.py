#!/usr/bin/env python3
"""
Tests for lipid bilayer analysis (area per lipid and membrane thickness).

Uses the equilibration_folder membrane-protein test system from analysis examples.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

pytest.importorskip("lipyphilic")

from gatewizard.utils.lipid_bilayer_analysis import (  # noqa: E402
    BilayerTrajectoryAnalyzer,
    run_bilayer_analysis,
)
from gatewizard.utils.namd_analysis import (  # noqa: E402
    BilayerTrajectoryAnalyzer as ExportedBilayerAnalyzer,
    run_bilayer_analysis as exported_run_bilayer_analysis,
)

LIPID_SEL = "resname PC and name P31"


@pytest.fixture(scope="module")
def equilibration_bilayer_data():
    data_dir = Path(__file__).parent / "analysis_examples" / "equilibration_folder"
    topology = data_dir / "system.pdb"
    trajectories = [
        data_dir / "step1_equilibration.dcd",
        data_dir / "step2_equilibration.dcd",
        data_dir / "step3_equilibration.dcd",
    ]
    if not topology.exists() or not all(t.exists() for t in trajectories):
        pytest.skip("equilibration_folder test data not found")

    file_times = {
        "step1_equilibration.dcd": 0.1,
        "step2_equilibration.dcd": 0.1,
        "step3_equilibration.dcd": 0.1,
    }
    return topology, trajectories, file_times


class TestBilayerAnalysisExports:
    def test_reexported_from_namd_analysis(self):
        assert ExportedBilayerAnalyzer is BilayerTrajectoryAnalyzer
        assert exported_run_bilayer_analysis is run_bilayer_analysis


class TestAreaPerLipid:
    def test_calculate_area_per_lipid(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        data = analyzer.calculate_area_per_lipid(lipid_sel=LIPID_SEL)

        assert len(data["resids"]) == 122
        assert data["areas"].shape[0] == 122
        assert data["areas"].shape[1] >= 1
        mean_area = float(data["mean_area_per_lipid"].mean())
        assert 40.0 < mean_area < 120.0

    def test_run_bilayer_analysis_area_per_lipid(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        result = run_bilayer_analysis(
            topology,
            trajectories,
            analysis_type="area_per_lipid",
            lipid_sel=LIPID_SEL,
            file_times=file_times,
        )

        assert result["analysis_type"] == "area_per_lipid"
        assert result["x_label"] == "Time (ns)"
        assert result["y_label"] == "Area per lipid (Å²)"
        assert len(result["lipid_resids"]) == 122
        assert len(result["per_lipid_areas"]) == 122
        assert 40.0 < result["stats"]["mean"] < 120.0
        assert "mean" in result["stats"]


class TestMembraneThickness:
    def test_calculate_membrane_thickness(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        data = analyzer.calculate_membrane_thickness(lipid_sel=LIPID_SEL)

        assert len(data["thickness"]) >= 1
        mean_thickness = float(data["thickness"].mean())
        assert 30.0 < mean_thickness < 60.0

    def test_run_bilayer_analysis_membrane_thickness(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        result = run_bilayer_analysis(
            topology,
            trajectories,
            analysis_type="membrane_thickness",
            lipid_sel=LIPID_SEL,
            file_times=file_times,
        )

        assert result["analysis_type"] == "membrane_thickness"
        assert result["y_label"] == "Membrane thickness (Å)"
        assert 30.0 < result["stats"]["mean"] < 60.0
        assert result["n_bins"] == 1
        assert "std" in result["stats"]

    def test_bilayer_time_axis_spans_file_times(self, equilibration_bilayer_data):
        """Regression: x-axis must honor file_times, not fall back to 0.01 ns/frame."""
        topology, trajectories, file_times = equilibration_bilayer_data
        expected_total = sum(file_times.values())
        result = run_bilayer_analysis(
            topology,
            trajectories,
            analysis_type="membrane_thickness",
            lipid_sel=LIPID_SEL,
            file_times=file_times,
        )
        x = result["x"]
        assert len(x) == len(result["y"])
        assert x[-1] == pytest.approx(expected_total, rel=0.02, abs=0.05)

    def test_pbc_straddling_thickness_is_folded_to_bilayer_gap(self):
        """Water-gap path Lz−d must not be reported as thickness."""
        import numpy as np

        from gatewizard.utils.lipid_bilayer_analysis import (
            _correct_pbc_straddling_thickness,
        )

        box_z = 135.0
        water_gap = np.array([100.0, 100.3, 99.8])
        fixed = _correct_pbc_straddling_thickness(water_gap, box_z)
        np.testing.assert_allclose(fixed, box_z - water_gap)
        true_d = np.array([35.0, 34.7, 35.2])
        np.testing.assert_allclose(
            _correct_pbc_straddling_thickness(true_d, box_z), true_d
        )


class TestBilayerAnalysisValidation:
    def test_unsupported_analysis_type(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        with pytest.raises(ValueError, match="Unsupported bilayer analysis type"):
            run_bilayer_analysis(
                topology,
                trajectories,
                analysis_type="flip_flop",
                lipid_sel=LIPID_SEL,
                file_times=file_times,
            )
