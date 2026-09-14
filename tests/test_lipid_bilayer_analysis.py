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
pytest.importorskip("freud")

from gatewizard.utils.fatslim_apl import FatslimError  # noqa: E402
from gatewizard.utils.lipid_bilayer_analysis import (  # noqa: E402
    BilayerTrajectoryAnalyzer,
    _clip_polygon_halfplane,
    _evapl_clip_areas,
    _gridmat_assign_areas,
    _gridmat_build_xy_grid,
    _lipyphilic_areaperlipid_supports_exclude,
    _polygon_area_xy,
    _resolve_apl_method,
    _voronoi_atom_areas,
    _vtmc_assign_areas,
    run_bilayer_analysis,
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
    def test_canonical_import(self):
        assert BilayerTrajectoryAnalyzer is not None
        assert callable(run_bilayer_analysis)


class TestVoronoiAplHelpers:
    def test_voronoi_areas_tile_box(self):
        import numpy as np

        # 2×2 grid in a 20×20 box → each cell ≈ 100 Å²
        pos = np.array(
            [
                [5.0, 5.0, 0.0],
                [15.0, 5.0, 0.0],
                [5.0, 15.0, 0.0],
                [15.0, 15.0, 0.0],
            ],
            dtype=float,
        )
        areas = _voronoi_atom_areas(pos, 20.0, 20.0)
        assert areas.shape == (4,)
        np.testing.assert_allclose(np.sum(areas), 400.0, rtol=1e-5)
        np.testing.assert_allclose(areas, 100.0, rtol=1e-4)

    def test_gridmat_assign_four_lipids(self):
        import numpy as np

        grid = _gridmat_build_xy_grid(20.0, 20.0, 2, conserve_ratio=False)
        lipid_xy = np.array(
            [[5.0, 5.0], [15.0, 5.0], [5.0, 15.0], [15.0, 15.0]],
            dtype=float,
        )
        lipid_res = np.array([0, 1, 2, 3], dtype=int)
        areas = _gridmat_assign_areas(grid, lipid_xy, lipid_res, np.empty((0, 2)), 20.0, 20.0)
        assert len(areas) == 4
        np.testing.assert_allclose(sum(areas.values()), 400.0, rtol=1e-4)
        np.testing.assert_allclose(list(areas.values()), 100.0, rtol=1e-4)

    def test_vtmc_assign_four_lipids_no_protein(self):
        import numpy as np

        lipid_xy = np.array(
            [[5.0, 5.0], [15.0, 5.0], [5.0, 15.0], [15.0, 15.0]],
            dtype=float,
        )
        lipid_res = np.array([0, 1, 2, 3], dtype=int)
        rng = np.random.default_rng(0)
        areas = _vtmc_assign_areas(
            lipid_xy, lipid_res, np.empty((0, 2)), 20.0, 20.0, 40_000, 1.7, rng
        )
        assert len(areas) == 4
        np.testing.assert_allclose(sum(areas.values()), 400.0, rtol=1e-3)
        # Equal Voronoi tiles → ~100 Å² each (MC noise).
        for a in areas.values():
            assert 70.0 < a < 130.0

    def test_vtmc_protein_disk_reduces_area(self):
        import numpy as np

        lipid_xy = np.array(
            [[5.0, 5.0], [15.0, 5.0], [5.0, 15.0], [15.0, 15.0]],
            dtype=float,
        )
        lipid_res = np.array([0, 1, 2, 3], dtype=int)
        # Dense protein cluster in one quadrant.
        protein = np.array([[5.0, 5.0], [4.5, 5.0], [5.0, 4.5], [5.5, 5.0]], dtype=float)
        rng = np.random.default_rng(1)
        areas = _vtmc_assign_areas(
            lipid_xy, lipid_res, protein, 20.0, 20.0, 40_000, 2.0, rng
        )
        assert sum(areas.values()) < 400.0 - 5.0
        assert areas[0] < areas[1]

    def test_evapl_halfplane_keeps_ref_side(self):
        import numpy as np

        square = np.array(
            [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]],
            dtype=float,
        )
        clipped = _clip_polygon_halfplane(square, [2.0, 5.0], [8.0, 5.0])
        assert clipped.shape[0] >= 3
        area = _polygon_area_xy(clipped)
        assert 40.0 < area < 60.0
        # All remaining vertices stay on the lipid side of x=5.
        assert np.all(clipped[:, 0] <= 5.0 + 1e-9)

    def test_evapl_clip_four_lipids_protein_one_cell(self):
        import numpy as np

        lipid = np.array(
            [
                [5.0, 5.0, 0.0],
                [15.0, 5.0, 0.0],
                [5.0, 15.0, 0.0],
                [15.0, 15.0, 0.0],
            ],
            dtype=float,
        )
        empty = _evapl_clip_areas(lipid, np.empty((0, 3)), 20.0, 20.0)
        np.testing.assert_allclose(empty, 100.0, rtol=1e-4)

        # Protein cluster inside the first lipid's cell only.
        protein = np.array(
            [[6.0, 5.0, 0.0], [5.5, 5.2, 0.0], [5.2, 4.7, 0.0]],
            dtype=float,
        )
        clipped = _evapl_clip_areas(lipid, protein, 20.0, 20.0)
        assert clipped[0] < empty[0] - 1.0
        np.testing.assert_allclose(clipped[1:], empty[1:], rtol=1e-4)

    def test_evapl_multiatom_excludes_both_sides(self):
        """Two atoms on opposite sides of a lipid must both cut the cell."""
        import numpy as np

        lipid = np.array(
            [
                [5.0, 5.0, 0.0],
                [15.0, 5.0, 0.0],
                [5.0, 15.0, 0.0],
                [15.0, 15.0, 0.0],
            ],
            dtype=float,
        )
        empty = _evapl_clip_areas(lipid, np.empty((0, 3)), 20.0, 20.0)
        # COM of these two is at the lipid → old single-COM clip would barely cut.
        protein = np.array(
            [[3.0, 5.0, 0.0], [7.0, 5.0, 0.0]],
            dtype=float,
        )
        com_only = _evapl_clip_areas(
            lipid, np.array([[5.0, 5.0, 0.0]], dtype=float), 20.0, 20.0
        )
        multi = _evapl_clip_areas(lipid, protein, 20.0, 20.0)
        # Atom at the lipid site is degenerate (no cut); multi-atom must shrink more.
        assert multi[0] < empty[0] - 10.0
        assert multi[0] < com_only[0] - 5.0
        np.testing.assert_allclose(multi[1:], empty[1:], rtol=1e-4)

    def test_evapl_atom_on_lipid_side_of_com_is_clipped(self):
        """An atom closer to the lipid than the COM must still remove area."""
        import numpy as np
        from gatewizard.utils.lipid_bilayer_analysis import (
            _clip_polygon_halfplane,
            _normalize_xy_positions,
            _unwrap_xy_to_ref,
            _voronoi_compute,
        )

        lipid = np.array(
            [
                [5.0, 5.0, 0.0],
                [15.0, 5.0, 0.0],
                [5.0, 15.0, 0.0],
                [15.0, 15.0, 0.0],
            ],
            dtype=float,
        )
        # Far atom + near atom: COM sits between lipid and far atom; near atom
        # is on the lipid side of the COM bisector.
        near = np.array([5.8, 5.0, 0.0])
        far = np.array([8.5, 5.0, 0.0])
        protein = np.vstack([near, far])
        areas = _evapl_clip_areas(lipid, protein, 20.0, 20.0)

        pos = _normalize_xy_positions(lipid)
        base = _voronoi_compute(pos, 20.0, 20.0)
        ref = pos[0, :2]
        verts = _unwrap_xy_to_ref(
            np.asarray(base.polytopes[0], dtype=float)[:, :2], ref, 20.0, 20.0
        )
        com = 0.5 * (near[:2] + far[:2])
        com_clip = _clip_polygon_halfplane(verts, ref, com)
        # Near atom is on the lipid side of the COM cut → still inside COM polygon.
        mid = 0.5 * (ref + com)
        normal = ref - com
        assert float(np.dot(near[:2] - mid, normal)) >= -1e-9
        # Multi-atom EVAPL area must be strictly smaller than COM-only area.
        com_area = _polygon_area_xy(com_clip)
        assert areas[0] < com_area - 1.0

    def test_evapl_soluble_head_outside_z_range_ignored(self):
        """Atoms far above the headgroup Z-range must not change APL."""
        import numpy as np

        lipid = np.array(
            [
                [5.0, 5.0, 0.0],
                [15.0, 5.0, 0.0],
                [5.0, 15.0, 0.0],
                [15.0, 15.0, 0.0],
            ],
            dtype=float,
        )
        empty = _evapl_clip_areas(lipid, np.empty((0, 3)), 20.0, 20.0)
        # Frame path gates by Z before calling _evapl_clip_areas; simulate that
        # by not passing out-of-range atoms (helper returns empty).
        from gatewizard.utils.lipid_bilayer_analysis import _exclude_atoms_in_z_range

        class _FakeAtomGroup:
            def __init__(self, positions, names=None):
                self.positions = np.asarray(positions, dtype=float)
                self._names = names or ["C"] * len(positions)

            def __len__(self):
                return len(self.positions)

            def __getitem__(self, mask):
                idx = np.where(mask)[0] if mask.dtype == bool else mask
                return _FakeAtomGroup(self.positions[idx], [self._names[i] for i in idx])

            def select_atoms(self, sel):
                if "H*" in sel:
                    return _FakeAtomGroup(self.positions, self._names)
                return self

            def wrap(self, inplace=True):
                return None

        class _FakeUniverse:
            def select_atoms(self, sel):
                # Large XY footprint far above phosphates (z≈30 vs headgroup z∈[0,0]).
                return _FakeAtomGroup(
                    [[5.0, 5.0, 30.0], [6.0, 5.0, 31.0], [4.0, 6.0, 29.0], [7.0, 4.0, 32.0]]
                )

        gated = _exclude_atoms_in_z_range(_FakeUniverse(), "protein", 0.0, 0.0, keep_z=True)
        assert gated.shape[0] == 0
        clipped = _evapl_clip_areas(lipid, gated, 20.0, 20.0)
        np.testing.assert_allclose(clipped, empty, rtol=1e-5)

        # Same XY but inside the headgroup Z-range must reduce area.
        in_range = np.array([[5.0, 5.5, 0.0], [6.0, 5.0, 0.0]], dtype=float)
        reduced = _evapl_clip_areas(lipid, in_range, 20.0, 20.0)
        assert reduced[0] < empty[0] - 1.0

    def test_unknown_apl_method_is_rejected(self):
        with pytest.raises(ValueError, match="Unsupported apl_method"):
            _resolve_apl_method("not_a_real_method", "protein")

    def test_resolve_apl_method_fatslim_is_default(self):
        assert _resolve_apl_method("auto", "protein") == "fatslim"
        assert _resolve_apl_method("auto", "") == "fatslim"
        assert _resolve_apl_method(None, "protein") == "fatslim"
        assert _resolve_apl_method("fatslim", "protein") == "fatslim"
        assert _resolve_apl_method("fatslim_cli", "") == "fatslim"
        assert _resolve_apl_method("evapl", "protein") == "evapl"
        assert _resolve_apl_method("lipyphilic", "") == "lipyphilic"

    def test_exclude_sites_reduce_lipid_share(self):
        import numpy as np

        lipid = np.array(
            [
                [5.0, 5.0, 0.0],
                [15.0, 5.0, 0.0],
                [5.0, 15.0, 0.0],
                [15.0, 15.0, 0.0],
            ],
            dtype=float,
        )
        lipid_only = _voronoi_atom_areas(lipid, 20.0, 20.0)
        with_exclude = _voronoi_atom_areas(
            np.vstack([lipid, [[10.0, 10.0, 0.0]]]),
            20.0,
            20.0,
        )
        lipid_share = with_exclude[:4]
        assert float(np.sum(lipid_share)) < float(np.sum(lipid_only)) - 1.0
        np.testing.assert_allclose(np.sum(with_exclude), 400.0, rtol=1e-5)


class TestAreaPerLipid:
    def test_default_fatslim_missing_raises(self, equilibration_bilayer_data, monkeypatch):
        monkeypatch.setattr(
            "gatewizard.utils.fatslim_apl.resolve_fatslim_executable", lambda: None
        )
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        with pytest.raises(FatslimError, match="GATEWIZARD_FATSLIM|install_fatslim_env"):
            analyzer.calculate_area_per_lipid(
                lipid_sel=LIPID_SEL, exclude_sel="", start=0, stop=2
            )

    def test_calculate_area_per_lipid(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        data = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL, exclude_sel="", apl_method="lipyphilic", start=0, stop=2
        )

        assert len(data["resids"]) == 122
        assert data["areas"].shape[0] == 122
        assert data["areas"].shape[1] >= 1
        mean_area = float(data["mean_area_per_lipid"].mean())
        assert 40.0 < mean_area < 120.0

    def test_lipyphilic_method_uses_areaperlipid(self, equilibration_bilayer_data, monkeypatch):
        """apl_method=lipyphilic must call lipyphilic.AreaPerLipid, not freud frame loop."""
        import lipyphilic.analysis.area_per_lipid as apl_mod

        calls = []
        Real = apl_mod.AreaPerLipid

        class Spy(Real):
            def __init__(self, *args, **kwargs):
                calls.append(kwargs)
                super().__init__(*args, **kwargs)

        monkeypatch.setattr(apl_mod, "AreaPerLipid", Spy)
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        data = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL, exclude_sel="", apl_method="lipyphilic", start=0, stop=2
        )
        assert calls, "AreaPerLipid was not constructed"
        assert float(data["mean_area_per_lipid"].mean()) > 40.0

    def test_lipyphilic_exclude_requires_upstream_support(self, equilibration_bilayer_data):
        if _lipyphilic_areaperlipid_supports_exclude():
            pytest.skip("installed lipyphilic already supports exclude_sel")
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        with pytest.raises(ValueError, match="exclude support|requirements-lipyphilic-git|fatslim|evapl"):
            analyzer.calculate_area_per_lipid(
                lipid_sel=LIPID_SEL,
                exclude_sel="protein",
                apl_method="lipyphilic",
                start=0,
                stop=2,
            )

    @pytest.mark.skipif(
        not _lipyphilic_areaperlipid_supports_exclude(),
        reason="needs lipyphilic from requirements-lipyphilic-git.txt (PR #164)",
    )
    def test_lipyphilic_exclude_lowers_mean_when_supported(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        no_excl = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL, exclude_sel="", apl_method="lipyphilic", start=0, stop=2
        )
        with_prot = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL,
            exclude_sel="protein",
            apl_method="lipyphilic",
            exclude_cutoff=30.0,
            start=0,
            stop=2,
        )
        mean_no = float(no_excl["mean_area_per_lipid"].mean())
        mean_prot = float(with_prot["mean_area_per_lipid"].mean())
        assert mean_prot < mean_no - 0.5

    def test_evapl_protein_exclude_lowers_mean_apl(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        lipyphilic = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL, exclude_sel="", apl_method="lipyphilic", start=0, stop=2
        )
        analyzer2 = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        evapl = analyzer2.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL,
            exclude_sel="protein",
            apl_method="evapl",
            start=0,
            stop=2,
        )
        mean_box = float(lipyphilic["mean_area_per_lipid"].mean())
        mean_evapl = float(evapl["mean_area_per_lipid"].mean())
        assert mean_evapl < mean_box - 0.5
        assert 40.0 < mean_evapl < 120.0

    def test_gridmat_apl_on_equilibration(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        no_excl = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL, exclude_sel="", apl_method="gridmat", start=0, stop=2
        )
        with_prot = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL,
            exclude_sel="protein",
            apl_method="gridmat",
            gridmat_precision=13.0,
            start=0,
            stop=2,
        )
        mean_no = float(no_excl["mean_area_per_lipid"].mean())
        mean_prot = float(with_prot["mean_area_per_lipid"].mean())
        assert 40.0 < mean_no < 120.0
        assert mean_prot < mean_no

    def test_vtmc_apl_on_equilibration(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        analyzer = BilayerTrajectoryAnalyzer(topology, trajectories, file_times=file_times)
        no_excl = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL,
            exclude_sel="",
            apl_method="vtmc",
            vtmc_n_samples=20_000,
            start=0,
            stop=2,
        )
        with_prot = analyzer.calculate_area_per_lipid(
            lipid_sel=LIPID_SEL,
            exclude_sel="protein",
            apl_method="vtmc",
            vtmc_n_samples=20_000,
            vtmc_protein_radius=1.7,
            start=0,
            stop=2,
        )
        mean_no = float(no_excl["mean_area_per_lipid"].mean())
        mean_prot = float(with_prot["mean_area_per_lipid"].mean())
        assert 40.0 < mean_no < 120.0
        assert mean_prot < mean_no

    def test_run_bilayer_analysis_area_per_lipid(self, equilibration_bilayer_data):
        topology, trajectories, file_times = equilibration_bilayer_data
        result = run_bilayer_analysis(
            topology,
            trajectories,
            analysis_type="area_per_lipid",
            lipid_sel=LIPID_SEL,
            exclude_sel="",
            apl_method="lipyphilic",
            file_times=file_times,
            start=0,
            stop=2,
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
