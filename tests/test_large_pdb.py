"""Amber column repair for packed membranes that exceed PDB field widths."""

from pathlib import Path

from gatewizard.utils.large_pdb import (
    finalize_for_tleap,
    needs_amber_pdb_rewrite,
    read_cryst1_box,
    repair_saved_pdb,
    rewrite_if_needed,
    rewrite_pdb_for_amber,
)


def _atom(serial: str, name: str, resname: str, resid: str, x: str, y: str, z: str) -> str:
    return f"ATOM  {serial} {name} {resname} {resid} {x} {y} {z}  1.00  0.00\n"


def test_overflow_residue_and_spilled_chloride_are_realigned(tmp_path: Path):
    source = tmp_path / "packed.pdb"
    source.write_text(
        "CRYST1  100.000   80.000   90.000  90.00  90.00  90.00 P 1           1\n"
        "ATOM  541040  O   WAT  100000     58.371 115.505  51.092  1.00  0.00\n"
        "HETATM 272291Cl- CL     12      10.000  20.000  30.000  1.00  0.00\n"
        "ATOM      3 1HG  OL      4      1.000   2.000-3.000  1.00  0.00\n"
        "END\n",
        encoding="utf-8",
    )
    assert needs_amber_pdb_rewrite(source)
    dest = tmp_path / "aligned.pdb"
    stats = rewrite_pdb_for_amber(source, dest)
    assert stats["atoms"] == 3
    text = dest.read_text(encoding="utf-8")
    assert "100000" not in text
    assert "Cl-" in text
    assert "1Cl-" not in text
    assert "1HG" in text
    assert not needs_amber_pdb_rewrite(dest)
    box = read_cryst1_box(dest)
    assert box is not None
    assert box[0] == 100.0
    assert box[2] == 90.0


def test_repeated_residue_numbers_stay_unique_for_tleap(tmp_path: Path):
    source = tmp_path / "chains.pdb"
    lines = ["CRYST1   70.000   70.000   80.000  90.00  90.00  90.00 P 1           1\n"]
    # Two cysteines numbered 30, with another residue between them so they
    # stay separate. The second 30 has to take a free number.
    rows = (
        (1, "A", "CYX", 30, 0.0),
        (2, "A", "ALA", 5, 5.0),
        (3, "B", "CYX", 30, 2.0),
    )
    for serial, chain, resname, resid, x in rows:
        lines.append(
            f"ATOM  {serial:5d}  SG  {resname} {chain}{resid:4d}    "
            f"{x:8.3f}   0.000   0.000  1.00  0.00\n"
        )
    lines.append("END\n")
    source.write_text("".join(lines), encoding="utf-8")
    dest = tmp_path / "system_for_tleap.pdb"
    rewrite_pdb_for_amber(source, dest)
    fragment = tmp_path / "tleap_geometry.in"
    finalize_for_tleap(dest, fragment)
    text = fragment.read_text(encoding="utf-8")
    assert "bond system.30.SG system." in text
    assert "set system box { 70.000 70.000 80.000 }" in text
    # The second SG must not keep residue 30, or tleap would bond the first one twice.
    body = dest.read_text(encoding="utf-8")
    assert body.count("  30") == 1


def test_readable_pdb_is_left_unchanged(tmp_path: Path):
    source = tmp_path / "small.pdb"
    source.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00\nEND\n",
        encoding="utf-8",
    )
    dest = tmp_path / "bilayer_for_pdb4amber.pdb"
    dest.write_text("stale\n", encoding="utf-8")
    assert rewrite_if_needed(source, dest) is False
    assert not dest.exists()


def test_saved_pdb_regains_cryst1_and_columns(tmp_path: Path):
    source = tmp_path / "system_for_tleap.pdb"
    source.write_text(
        "CRYST1  243.816  243.691  172.096  90.00  90.00  90.00 P 1           1\n"
        "ATOM      1  O   WAT A   1      58.371 115.505  51.092  1.00  0.00\n"
        "END\n",
        encoding="utf-8",
    )
    saved = tmp_path / "system.pdb"
    saved.write_text(
        "ATOM  541040  O   WAT  100000     58.371 115.505  51.092  1.00  0.00\nEND\n",
        encoding="utf-8",
    )
    info = repair_saved_pdb(saved, source)
    assert info["rewritten"] is True
    assert info["cryst1_copied"] is True
    assert not needs_amber_pdb_rewrite(saved)
    assert read_cryst1_box(saved)[:3] == (243.816, 243.691, 172.096)
