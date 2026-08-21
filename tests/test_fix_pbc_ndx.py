"""Unit tests for Fix PBC multi-group index merging."""

from pathlib import Path

import pytest

from gatewizard.utils.trajectory_tools import (
    _normalize_group_name_list,
    _prepare_gromacs_fix_groups,
    _recommend_center_groups,
    build_compound_ndx,
    list_ndx_groups,
    parse_ndx_group_atoms,
)


SAMPLE_NDX = """\
[ System ]
   1   2   3   4   5   6   7   8   9  10
[ Protein ]
   1   2   3
[ PA ]
   4   5
[ PC ]
   6   7
[ OL ]
   8   9  10
[ MEMB ]
   4   5   6   7   8   9  10
[ System_dup ]
   1   2
"""


@pytest.fixture()
def ndx_path(tmp_path: Path) -> Path:
    p = tmp_path / "index.ndx"
    p.write_text(SAMPLE_NDX, encoding="utf-8")
    return p


def test_parse_ndx_group_atoms(ndx_path: Path):
    assert parse_ndx_group_atoms(ndx_path, "PA") == [4, 5]
    assert parse_ndx_group_atoms(ndx_path, "OL") == [8, 9, 10]
    assert parse_ndx_group_atoms(ndx_path, "Missing") == []


def test_build_compound_ndx_union(ndx_path: Path, tmp_path: Path):
    dest = tmp_path / "merged.ndx"
    info = build_compound_ndx(ndx_path, ["PA", "PC", "OL"], "GW_CENTER", dest)
    assert info["compound_name"] == "GW_CENTER"
    assert info["source_groups"] == ["PA", "PC", "OL"]
    assert info["n_atoms"] == 7  # unique 4..10
    atoms = parse_ndx_group_atoms(dest, "GW_CENTER")
    assert atoms == [4, 5, 6, 7, 8, 9, 10]
    # Original groups preserved
    assert parse_ndx_group_atoms(dest, "Protein") == [1, 2, 3]
    names = {g["name"] for g in list_ndx_groups(dest)}
    assert "GW_CENTER" in names
    assert "PA" in names


def test_build_compound_ndx_dedupes_overlap(ndx_path: Path, tmp_path: Path):
    dest = tmp_path / "overlap.ndx"
    info = build_compound_ndx(ndx_path, ["MEMB", "PA"], "GW_CENTER", dest)
    assert info["n_atoms"] == 7  # MEMB already covers PA


def test_build_compound_ndx_missing_group(ndx_path: Path, tmp_path: Path):
    with pytest.raises(ValueError, match="not found"):
        build_compound_ndx(ndx_path, ["PA", "NOPE"], "GW_CENTER", tmp_path / "x.ndx")


def test_normalize_group_name_list_prefers_names():
    assert _normalize_group_name_list(["PA", "PC"], "OL") == ["PA", "PC"]
    assert _normalize_group_name_list(None, "Protein") == ["Protein"]
    assert _normalize_group_name_list(["PA", "PA", ""], None) == ["PA"]


def test_prepare_multi_center_resolves_gw_center(ndx_path: Path, tmp_path: Path):
    work = tmp_path / "work"
    work.mkdir()
    resolved = _prepare_gromacs_fix_groups(
        ndx=ndx_path,
        tpr=None,
        gmx="gmx",
        work_dir=work,
        center_groups=["PA", "PC", "OL"],
        output_group="System",
    )
    assert resolved["center_group"] == "GW_CENTER"
    assert resolved["center_sources"] == ["PA", "PC", "OL"]
    assert resolved["center_label"] == "GW_CENTER = PA+PC+OL"
    assert resolved["output_group"] == "System"
    assert resolved["ndx"] is not None
    assert Path(resolved["ndx"]).is_file()
    # trjconv stdin would use compound name
    assert parse_ndx_group_atoms(resolved["ndx"], "GW_CENTER") == [4, 5, 6, 7, 8, 9, 10]


def test_prepare_single_group_keeps_name(ndx_path: Path, tmp_path: Path):
    work = tmp_path / "work"
    work.mkdir()
    resolved = _prepare_gromacs_fix_groups(
        ndx=ndx_path,
        tpr=None,
        gmx="gmx",
        work_dir=work,
        center_group="SOLU_MEMB",  # absent → still uses explicit singular
        output_groups=["System"],
    )
    # Explicit singular is kept even if not in ndx (trjconv will fail later if invalid)
    assert resolved["center_group"] == "SOLU_MEMB"
    assert resolved["ndx"] == ndx_path


def test_recommend_split_lipids_without_solu_memb():
    groups = [
        {"name": "System"},
        {"name": "Protein"},
        {"name": "PA"},
        {"name": "PC"},
        {"name": "OL"},
        {"name": "MEMB"},
    ]
    recommended, recommended_groups = _recommend_center_groups(groups)
    assert "PA" in recommended_groups
    assert "PC" in recommended_groups
    assert "OL" in recommended_groups
    assert "Protein" in recommended_groups
    assert recommended == recommended_groups[0]


def test_recommend_solu_memb_wins():
    groups = [
        {"name": "SOLU_MEMB"},
        {"name": "PA"},
        {"name": "PC"},
        {"name": "OL"},
    ]
    recommended, recommended_groups = _recommend_center_groups(groups)
    assert recommended == "SOLU_MEMB"
    assert recommended_groups == ["SOLU_MEMB"]
