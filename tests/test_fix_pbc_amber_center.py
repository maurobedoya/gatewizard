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
    assert sel.startswith("protein or resname ")
    assert "POPC" in sel
    assert "PA" in sel
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


def test_recommend_from_psf_charmm_popc(tmp_path):
    psf = tmp_path / "system.psf"
    psf.write_text(
        "PSF\n\n"
        "       1 !NTITLE\n"
        " REMARKS test\n\n"
        "       4 !NATOM\n"
        "       1 PROT 1    MET  N    NH3   -0.300000       14.0070           0\n"
        "       2 MEMB 1    POPC N    NDL    0.620000       14.0070           0\n"
        "       3 MEMB 1    POPC C    CTL    0.000000       12.0110           0\n"
        "       4 WAT  1    TIP3 OH2  OT    -0.834000       15.9994           0\n"
        "\n"
        "       0 !NBOND: bonds\n",
        encoding="utf-8",
    )
    sel, lipids = recommend_cpptraj_center_selection(psf)
    assert lipids == ["POPC"]
    assert sel == "protein or resname POPC"
    traj = tmp_path / "prod.dcd"
    traj.write_bytes(b"")
    info = detect_pbc_engine(str(psf), [str(traj)], engine_hint="namd")
    assert info["engine"] == "namd"
    assert info["lipid_resnames"] == ["POPC"]
    assert info["recommended_center_selection"] == "protein or resname POPC"


def test_recommend_from_pdb_openmm(tmp_path):
    pdb = tmp_path / "system.pdb"
    pdb.write_text(
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  N   POPC B   1       1.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      3  P   POPC B   1       2.000   0.000   0.000  1.00  0.00           P\n"
        "ATOM      4  OH2 TIP3 C   1       3.000   0.000   0.000  1.00  0.00           O\n",
        encoding="utf-8",
    )
    sel, lipids = recommend_cpptraj_center_selection(pdb)
    assert lipids == ["POPC"]
    assert sel == "protein or resname POPC"
    info = detect_pbc_engine(str(pdb), [], engine_hint="openmm")
    assert info["engine"] == "openmm"
    assert info["recommended_center_selection"] == "protein or resname POPC"
