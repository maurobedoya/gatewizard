#!/usr/bin/env python3
"""Unit tests for the optional GridMAT-MD.pl APL bridge."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.utils.gridmat_md_apl import (  # noqa: E402
    GRIDMAT_MD_ENV_VAR,
    GridmatMdError,
    normalize_gridmat_md_jobs,
    parse_gridmat_ave_apl,
    require_gridmat_md,
    resolve_gridmat_md_script,
)
from gatewizard.utils.lipid_bilayer_analysis import _resolve_apl_method  # noqa: E402


def test_resolve_separates_gridmat_and_gridmat_md():
    assert _resolve_apl_method("gridmat") == "gridmat"
    assert _resolve_apl_method("grid") == "gridmat"
    assert _resolve_apl_method("gridmat_md") == "gridmat_md"
    assert _resolve_apl_method("gridmat_pl") == "gridmat_md"
    assert _resolve_apl_method("gridmat_cli") == "gridmat_md"


def test_normalize_gridmat_md_jobs():
    assert normalize_gridmat_md_jobs(0) >= 1
    assert normalize_gridmat_md_jobs(4) == 4


def test_parse_gridmat_ave_apl(tmp_path):
    path = tmp_path / "out.top_areas.dat"
    path.write_text("some header\nAve APL = 62.5\n", encoding="utf-8")
    assert parse_gridmat_ave_apl(path) == pytest.approx(62.5)


def test_require_missing_raises(monkeypatch):
    monkeypatch.delenv(GRIDMAT_MD_ENV_VAR, raising=False)
    monkeypatch.setattr(
        "gatewizard.utils.gridmat_md_apl.resolve_gridmat_md_script", lambda: None
    )
    with pytest.raises(GridmatMdError, match="GATEWIZARD_GRIDMAT_MD|install_gridmat"):
        require_gridmat_md()


def test_resolve_respects_env(monkeypatch, tmp_path):
    script = tmp_path / "GridMAT-MD.pl"
    script.write_text("#!/usr/bin/env perl\n", encoding="utf-8")
    monkeypatch.setenv(GRIDMAT_MD_ENV_VAR, str(script))
    assert resolve_gridmat_md_script() == str(script.resolve())
