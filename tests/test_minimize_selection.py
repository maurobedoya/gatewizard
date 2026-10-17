"""A selected residue moves. Every other atom stays where it was."""

from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pytest

from gatewizard.core.minimize_selection import minimize_selection


def _atom(serial, name, resname, resid, x, y, z, chain="A"):
    return (
        f"ATOM  {serial:5d} {name} {resname:>3} {chain}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
    )


def _write_strained_pdb(path: Path) -> None:
    lines = [
        _atom(1, " N  ", "ALA", 1, 0.0, 0.0, 0.0),
        _atom(2, " CA ", "ALA", 1, 5.0, 0.0, 0.0),
        _atom(3, " C  ", "ALA", 1, 6.2, 1.0, 0.0),
        _atom(4, " O  ", "ALA", 1, 7.2, 0.8, 0.7),
        _atom(5, " CB ", "ALA", 1, 5.2, -1.4, 0.4),
        _atom(6, " N  ", "GLY", 2, 1.2, 3.0, 0.2),
        _atom(7, " CA ", "GLY", 2, 2.4, 3.6, 0.1),
        _atom(8, " C  ", "GLY", 2, 3.2, 2.8, 1.1),
        _atom(9, " O  ", "GLY", 2, 2.8, 1.6, 1.3),
        _atom(10, "K   ", "K", 10, 0.0, 0.0, 40.0).replace("ATOM  ", "HETATM", 1),
        _atom(15, " C  ", "FOO", 4, 0.5, 0.5, 0.5),
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


def _pos(path: Path, resid: int, name: str) -> np.ndarray:
    universe = mda.Universe(str(path))
    atom = universe.select_atoms(f"resid {resid} and name {name}")[0]
    return np.asarray(atom.position, dtype=np.float64)


def test_selected_residue_moves_and_the_shell_stays(tmp_path: Path):
    pytest.importorskip("openmm")
    pdb = tmp_path / "strain.pdb"
    out = tmp_path / "min.pdb"
    _write_strained_pdb(pdb)
    before_ca = _pos(pdb, 1, "CA")
    before_shell = _pos(pdb, 2, "CA")
    before_far = _pos(pdb, 10, "K")
    report = minimize_selection(str(pdb), str(out), selection="resid 1", max_iterations=50)
    assert report["atoms_moved"] > 0
    assert report["skipped"] == []
    assert float(np.linalg.norm(_pos(out, 1, "CA") - before_ca)) > 0.1
    assert float(np.linalg.norm(_pos(out, 2, "CA") - before_shell)) < 1e-4
    assert float(np.linalg.norm(_pos(out, 10, "K") - before_far)) < 1e-4


def test_unknown_residue_is_relaxed_with_the_protein(tmp_path: Path):
    pytest.importorskip("openmm")
    pdb = tmp_path / "strain.pdb"
    out = tmp_path / "min.pdb"
    _write_strained_pdb(pdb)
    report = minimize_selection(str(pdb), str(out), selection="resid 1 or resid 4", max_iterations=20)
    assert report["skipped"] == []
    assert any(item["resname"] == "FOO" and item["reason"] == "Relaxed sterically" for item in report["relaxed"])
    assert out.is_file()


def test_potassium_ion_can_be_minimized(tmp_path: Path):
    pytest.importorskip("openmm")
    pdb = tmp_path / "strain.pdb"
    out = tmp_path / "min.pdb"
    _write_strained_pdb(pdb)
    text = pdb.read_text(encoding="utf-8")
    ion = _atom(16, "K   ", "K", 309, 3.0, -4.0, 3.0).replace("ATOM  ", "HETATM", 1)
    pdb.write_text(text.replace("END\n", ion + "END\n"), encoding="utf-8")
    before_shell = _pos(pdb, 1, "CA")
    before_far = _pos(pdb, 10, "K")
    report = minimize_selection(str(pdb), str(out), selection="resid 309", max_iterations=50)
    assert report["skipped"] == []
    assert report["movable"] == 1
    assert out.is_file()
    assert float(np.linalg.norm(_pos(out, 1, "CA") - before_shell)) < 1e-4
    assert float(np.linalg.norm(_pos(out, 10, "K") - before_far)) < 1e-4


def test_clashing_unknown_atom_moves_and_alanine_stays(tmp_path: Path):
    pdb = tmp_path / "strain.pdb"
    out = tmp_path / "min.pdb"
    _write_strained_pdb(pdb)
    before_foo = _pos(pdb, 4, "C")
    before_ala = _pos(pdb, 1, "N")
    before_ca = _pos(pdb, 1, "CA")
    report = minimize_selection(str(pdb), str(out), selection="resid 4", max_iterations=50)
    assert report["skipped"] == []
    assert report["relaxed"][0]["resname"] == "FOO"
    after_foo = _pos(out, 4, "C")
    assert float(np.linalg.norm(after_foo - before_ala)) > float(np.linalg.norm(before_foo - before_ala))
    assert float(np.linalg.norm(_pos(out, 1, "CA") - before_ca)) < 1e-4


def test_water_hydrogen_bonds_keep_their_length(tmp_path: Path):
    pdb = tmp_path / "water.pdb"
    out = tmp_path / "min.pdb"
    lines = [
        _atom(1, " N  ", "ALA", 1, 0.0, 0.0, 0.0),
        _atom(2, " CA ", "ALA", 1, 1.4, 0.2, 0.0),
        _atom(3, " C  ", "ALA", 1, 2.2, 1.2, 0.4),
        _atom(4, " O  ", "ALA", 1, 1.8, 2.3, 0.2),
        _atom(5, " CB ", "ALA", 1, 1.6, -1.0, 0.8),
        _atom(6, " O  ", "HOH", 5, 2.000, 0.000, 0.000),
        _atom(7, " H1 ", "HOH", 5, 3.179, 0.000, 0.000),
        _atom(8, " H2 ", "HOH", 5, 1.700, 1.141, 0.000),
        "CONECT    1    2\n",
        "END\n",
    ]
    pdb.write_text("".join(lines), encoding="utf-8")
    oh1 = float(np.linalg.norm(_pos(pdb, 5, "O") - _pos(pdb, 5, "H1")))
    oh2 = float(np.linalg.norm(_pos(pdb, 5, "O") - _pos(pdb, 5, "H2")))
    before_o = _pos(pdb, 5, "O")
    before_n = _pos(pdb, 1, "N")
    before_ca = _pos(pdb, 1, "CA")
    minimize_selection(str(pdb), str(out), selection="resid 5", max_iterations=50)
    assert abs(float(np.linalg.norm(_pos(out, 5, "O") - _pos(out, 5, "H1"))) - oh1) < 0.02
    assert abs(float(np.linalg.norm(_pos(out, 5, "O") - _pos(out, 5, "H2"))) - oh2) < 0.02
    assert abs(_water_angle(out) - _water_angle(pdb)) < 1.0
    assert float(np.linalg.norm(_pos(out, 5, "O") - before_n)) > float(np.linalg.norm(before_o - before_n))
    assert float(np.linalg.norm(_pos(out, 1, "CA") - before_ca)) < 1e-4


def _water_angle(path: Path) -> float:
    oxygen = _pos(path, 5, "O")
    first = _pos(path, 5, "H1") - oxygen
    second = _pos(path, 5, "H2") - oxygen
    cosine = float(np.dot(first, second) / (np.linalg.norm(first) * np.linalg.norm(second)))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def test_rmsd_cap_limits_the_unknown_atom(tmp_path: Path):
    pdb = tmp_path / "strain.pdb"
    out = tmp_path / "min.pdb"
    _write_strained_pdb(pdb)
    before = _pos(pdb, 4, "C")
    minimize_selection(str(pdb), str(out), selection="resid 4", max_iterations=50, max_rmsd=0.05)
    assert float(np.linalg.norm(_pos(out, 4, "C") - before)) <= 0.05 + 1e-6
