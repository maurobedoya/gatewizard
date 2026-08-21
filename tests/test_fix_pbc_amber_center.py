"""Amber Fix PBC defaults to protein + detected lipid residues."""

from pathlib import Path

from gatewizard.utils.trajectory_tools import (
    detect_pbc_engine,
    ordered_lipid_resnames,
    recommend_cpptraj_center_selection,
)


def _write_prmtop(path: Path, labels_blob: str) -> None:
    path.write_text(
        "%VERSION  VERSION_STAMP = V0001.000\n"
        "%FLAG TITLE\n"
        "%FORMAT(20a4)\n"
        "dummy\n"
        "%FLAG RESIDUE_LABEL\n"
        "%FORMAT(20a4)\n"
        f"{labels_blob}\n"
        "%FLAG POINTERS\n"
        "%FORMAT(10I8)\n"
        "       0\n",
        encoding="utf-8",
    )


def test_ordered_lipid_resnames_prefers_lipid21_parts():
    assert ordered_lipid_resnames(["WAT", "ALA", "OL", "PC", "PA", "NA"]) == [
        "PA",
        "PC",
        "OL",
    ]


def test_recommend_protein_only_without_lipids():
    sel, lipids = recommend_cpptraj_center_selection(residue_names=["ALA", "WAT", "NA"])
    assert sel == "protein"
    assert lipids == []


def test_recommend_protein_and_lipid21_popc():
    sel, lipids = recommend_cpptraj_center_selection(
        residue_names=["NMET", "ALA", "PA", "PC", "OL", "WAT"]
    )
    assert sel == "protein or resname PA PC OL"
    assert lipids == ["PA", "PC", "OL"]


def test_recommend_from_prmtop(tmp_path):
    # 4-character Amber fields: ALA, PA, PC, OL, WAT
    blob = "ALA PA  PC  OL  WAT "
    top = tmp_path / "system.prmtop"
    _write_prmtop(top, blob)
    sel, lipids = recommend_cpptraj_center_selection(top)
    assert lipids == ["PA", "PC", "OL"]
    assert sel == "protein or resname PA PC OL"


def test_detect_pbc_engine_returns_amber_membrane_default(tmp_path):
    top = tmp_path / "system.prmtop"
    _write_prmtop(top, "ALA PA  PC  OL  WAT ")
    traj = tmp_path / "prod.nc"
    traj.write_bytes(b"")
    info = detect_pbc_engine(str(top), [str(traj)], engine_hint="amber")
    assert info["engine"] == "amber"
    assert info["lipid_resnames"] == ["PA", "PC", "OL"]
    assert info["recommended_center_selection"] == "protein or resname PA PC OL"
