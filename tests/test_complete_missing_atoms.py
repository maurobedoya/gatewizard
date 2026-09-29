"""tleap completion of missing protein heavy atoms."""

import shutil
from pathlib import Path

import pytest

from gatewizard.core.preparation import complete_missing_heavy_atoms


def _atom(serial: int, name: str, resname: str, chain: str, resid: int, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  {name:<3s} {resname:>3s} {chain}{resid:4d}    "
        f"{x:8.3f}{x:8.3f}{x:8.3f}  1.00  0.00           {name[0]}\n"
    )


@pytest.mark.skipif(shutil.which("tleap") is None, reason="tleap not installed")
def test_complete_missing_met_sidechain_keeps_ligand(tmp_path: Path):
    src = tmp_path / "stub.pdb"
    src.write_text(
        _atom(1, "N", "MET", "A", 1, 1.0)
        + _atom(2, "CA", "MET", "A", 1, 2.0)
        + _atom(3, "C", "MET", "A", 1, 3.0)
        + _atom(4, "O", "MET", "A", 1, 4.0)
        + _atom(5, "CB", "MET", "A", 1, 5.0)
        + "HETATM    6  C1  Y01 A 900       9.000   9.000   9.000  1.00  0.00           C\n"
        + "END\n",
        encoding="utf-8",
    )
    out = tmp_path / "done.pdb"
    info = complete_missing_heavy_atoms(str(src), str(out))
    assert info["atoms_added"] > 0
    text = out.read_text(encoding="utf-8")
    names = {
        line[12:16].strip()
        for line in text.splitlines()
        if line.startswith("ATOM") and line[17:20].strip() == "MET"
    }
    assert {"CG", "SD", "CE"} <= names
    assert "Y01" in text
    assert "900" in text
    met_line = next(
        line
        for line in text.splitlines()
        if line.startswith("ATOM") and line[17:20].strip() == "MET"
    )
    assert met_line[21:22] == "A"
    assert int(met_line[22:26]) == 1


@pytest.mark.skipif(shutil.which("tleap") is None, reason="tleap not installed")
def test_complete_missing_ash_adds_carboxylic_proton(tmp_path: Path):
    """ASH must go to tleap as ASH; renaming ASP after completion leaves HD2 off."""
    src = tmp_path / "ash.pdb"
    src.write_text(
        _atom(1, "N", "ASH", "A", 42, 1.0)
        + _atom(2, "CA", "ASH", "A", 42, 2.0)
        + _atom(3, "C", "ASH", "A", 42, 3.0)
        + _atom(4, "O", "ASH", "A", 42, 4.0)
        + _atom(5, "CB", "ASH", "A", 42, 5.0)
        + _atom(6, "CG", "ASH", "A", 42, 6.0)
        + _atom(7, "OD1", "ASH", "A", 42, 7.0)
        + _atom(8, "OD2", "ASH", "A", 42, 8.0)
        + "END\n",
        encoding="utf-8",
    )
    out = tmp_path / "ash_done.pdb"
    complete_missing_heavy_atoms(str(src), str(out))
    names = {
        line[12:16].strip()
        for line in out.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM") and line[17:20].strip() == "ASH"
    }
    assert "HD2" in names


@pytest.mark.skipif(shutil.which("tleap") is None, reason="tleap not installed")
def test_complete_then_rename_asp_to_ash_misses_hd2(tmp_path: Path):
    """Regression: tleap on ASP then a name-only ASH change does not add HD2."""
    src = tmp_path / "asp.pdb"
    src.write_text(
        _atom(1, "N", "ASP", "A", 42, 1.0)
        + _atom(2, "CA", "ASP", "A", 42, 2.0)
        + _atom(3, "C", "ASP", "A", 42, 3.0)
        + _atom(4, "O", "ASP", "A", 42, 4.0)
        + _atom(5, "CB", "ASP", "A", 42, 5.0)
        + _atom(6, "CG", "ASP", "A", 42, 6.0)
        + _atom(7, "OD1", "ASP", "A", 42, 7.0)
        + _atom(8, "OD2", "ASP", "A", 42, 8.0)
        + "END\n",
        encoding="utf-8",
    )
    completed = tmp_path / "asp_done.pdb"
    complete_missing_heavy_atoms(str(src), str(completed))
    renamed = tmp_path / "ash_renamed.pdb"
    text = completed.read_text(encoding="utf-8").replace("ASP", "ASH")
    renamed.write_text(text, encoding="utf-8")
    names = {
        line[12:16].strip()
        for line in renamed.read_text(encoding="utf-8").splitlines()
        if line.startswith("ATOM") and line[17:20].strip() == "ASH"
    }
    assert "HD2" not in names
