#!/usr/bin/env python3
"""Unit tests for the optional FATSLiM CLI APL bridge."""

import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.utils.fatslim_apl import (  # noqa: E402
    FATSLIM_ENV_VAR,
    FatslimError,
    frame_index_chunks,
    is_available,
    normalize_fatslim_jobs,
    parse_apl_raw_csv,
    require_fatslim,
    resolve_fatslim_executable,
    run_fatslim_apl,
)
from gatewizard.utils.lipid_bilayer_analysis import _resolve_apl_method  # noqa: E402


def test_resolve_apl_defaults_to_fatslim():
    assert _resolve_apl_method("auto", "protein") == "fatslim"
    assert _resolve_apl_method("", None) == "fatslim"


def test_frame_index_chunks_half_open():
    assert frame_index_chunks(100, 1) == [(0, 100)]
    assert frame_index_chunks(10, 3) == [(0, 4), (4, 7), (7, 10)]
    assert normalize_fatslim_jobs(0) == 1
    assert normalize_fatslim_jobs(8) == 8


def test_parse_apl_raw_csv(tmp_path):
    csv_path = tmp_path / "apl_raw0.csv"
    csv_path.write_text(
        "resid,leaflet,x,y,z,apl\n"
        "1,upper leaflet,1.0,2.0,3.0,0.65\n"
        "2,lower leaflet,1.1,2.1,3.1,0.70\n",
        encoding="utf-8",
    )
    mapping = parse_apl_raw_csv(csv_path)
    assert mapping[1] == pytest.approx(65.0)
    assert mapping[2] == pytest.approx(70.0)


def test_require_fatslim_missing_raises(monkeypatch):
    monkeypatch.delenv(FATSLIM_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "gatewizard.utils.fatslim_apl.resolve_fatslim_executable", lambda: None
    )
    with pytest.raises(FatslimError, match="install_fatslim_env|GATEWIZARD_FATSLIM"):
        require_fatslim()


def test_resolve_respects_gatewizard_fatslim_env(monkeypatch, tmp_path):
    fake = tmp_path / "fatslim"
    fake.write_text("#!/bin/sh\necho ok\n", encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setenv(FATSLIM_ENV_VAR, str(fake))
    assert resolve_fatslim_executable() == str(fake.resolve())


def test_run_fatslim_apl_mocked(tmp_path, monkeypatch):
    """Export + subprocess + CSV parse without a real FATSLiM binary."""
    mda = pytest.importorskip("MDAnalysis")

    # Minimal fake universe via MDA test trajectory if available; otherwise skip.
    try:
        from MDAnalysisTests.datafiles import GRO, XTC
    except Exception:
        pytest.skip("MDAnalysisTests not installed")

    u = mda.Universe(GRO, XTC)
    lipid_sel = "resname SOL"  # may not exist
    # Prefer a selection that exists on the standard GRO
    if len(u.select_atoms(lipid_sel)) == 0:
        lipid_sel = "name OW"
    if len(u.select_atoms(lipid_sel)) == 0:
        lipid_sel = "all"
    membrane = u.select_atoms(lipid_sel)
    if membrane.n_residues < 2:
        pytest.skip("fixture topology unsuitable for APL mock")

    n_frames = min(2, len(u.trajectory))
    frame_indices = list(range(n_frames))

    def fake_run(cmd, capture_output=True, text=True, check=False):
        assert "--nthreads" in cmd
        assert "--begin-frame" in cmd
        assert "--end-frame" in cmd
        begin = int(cmd[cmd.index("--begin-frame") + 1])
        end = int(cmd[cmd.index("--end-frame") + 1])
        export_prefix = Path(cmd[cmd.index("--export-apl-raw") + 1])
        work = export_prefix.parent
        resids = membrane.residues.resids
        for i in range(end - begin):
            lines = ["resid,leaflet,x,y,z,apl\n"]
            for r in resids:
                lines.append(f"{int(r)},upper leaflet,0,0,0,0.50\n")
            (work / f"{export_prefix.name}{i}.csv").write_text(
                "".join(lines), encoding="utf-8"
            )
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(
        "gatewizard.utils.fatslim_apl.require_fatslim", lambda: "/fake/fatslim"
    )
    with patch("gatewizard.utils.fatslim_apl.subprocess.run", side_effect=fake_run):
        out = run_fatslim_apl(
            u,
            lipid_sel=lipid_sel,
            exclude_sel="",
            frame_indices=frame_indices,
            work_dir=tmp_path / "export",
            fatslim_nthreads=1,
            fatslim_jobs=1,
        )

    areas = np.asarray(out["areas"])
    assert areas.shape[0] == membrane.n_residues
    assert areas.shape[1] == n_frames
    assert np.nanmean(areas) == pytest.approx(50.0)
    assert out["meta"]["nthreads"] == 1
    assert out["meta"]["jobs"] == 1


def test_run_fatslim_apl_jobs_mocked(tmp_path, monkeypatch):
    """Two frame chunks write CSVs under separate chunk_* dirs."""
    mda = pytest.importorskip("MDAnalysis")
    try:
        from MDAnalysisTests.datafiles import GRO, XTC
    except Exception:
        pytest.skip("MDAnalysisTests not installed")

    u = mda.Universe(GRO, XTC)
    lipid_sel = "name OW" if len(u.select_atoms("name OW")) else "all"
    membrane = u.select_atoms(lipid_sel)
    if membrane.n_residues < 1 or len(u.trajectory) < 2:
        pytest.skip("fixture topology unsuitable for APL mock")

    n_frames = min(4, len(u.trajectory))
    frame_indices = list(range(n_frames))
    seen_ranges = []

    def fake_run(cmd, capture_output=True, text=True, check=False):
        begin = int(cmd[cmd.index("--begin-frame") + 1])
        end = int(cmd[cmd.index("--end-frame") + 1])
        nthreads = int(cmd[cmd.index("--nthreads") + 1])
        assert nthreads == 2
        seen_ranges.append((begin, end))
        export_prefix = Path(cmd[cmd.index("--export-apl-raw") + 1])
        work = export_prefix.parent
        resids = membrane.residues.resids
        for i in range(end - begin):
            lines = ["resid,leaflet,x,y,z,apl\n"]
            for r in resids:
                lines.append(f"{int(r)},upper leaflet,0,0,0,0.40\n")
            (work / f"{export_prefix.name}{i}.csv").write_text(
                "".join(lines), encoding="utf-8"
            )
        return MagicMock(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(
        "gatewizard.utils.fatslim_apl.require_fatslim", lambda: "/fake/fatslim"
    )
    with patch("gatewizard.utils.fatslim_apl.subprocess.run", side_effect=fake_run):
        out = run_fatslim_apl(
            u,
            lipid_sel=lipid_sel,
            exclude_sel="",
            frame_indices=frame_indices,
            work_dir=tmp_path / "export2",
            fatslim_nthreads=2,
            fatslim_jobs=2,
        )

    assert sorted(seen_ranges) == frame_index_chunks(n_frames, 2)
    areas = np.asarray(out["areas"])
    assert areas.shape == (membrane.n_residues, n_frames)
    assert out["meta"]["jobs"] == 2
    assert np.nanmean(areas) == pytest.approx(40.0)


@pytest.mark.skipif(not is_available(), reason="fatslim not installed")
def test_fatslim_live_smoke_if_installed():
    """Optional live check when companion env is on PATH / GATEWIZARD_FATSLIM."""
    assert require_fatslim()
