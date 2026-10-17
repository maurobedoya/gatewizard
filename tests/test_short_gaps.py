"""Short numbering gaps are filled only when they are 1 to 6 residues and SEQRES has the sequence."""

from pathlib import Path
import urllib.request

import pytest

from gatewizard.core.short_gaps import SHORT_GAP_MAX, classify_short_gaps, fill_short_gaps


def _write_gapped_pdb(path: Path) -> None:
    # Residues 1, 2, then 6 (three missing), and 1 then 9 on chain B (seven missing).
    lines = [
        "SEQRES   1 A    6  ALA GLY ALA ALA ALA GLY",
        "SEQRES   1 B    9  ALA ALA ALA ALA ALA ALA ALA ALA VAL",
        "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N",
        "ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C",
        "ATOM      3  C   ALA A   1       2.000   0.000   0.000  1.00  0.00           C",
        "ATOM      4  N   GLY A   2       3.000   0.000   0.000  1.00  0.00           N",
        "ATOM      5  CA  GLY A   2       4.000   0.000   0.000  1.00  0.00           C",
        "ATOM      6  C   GLY A   2       5.000   0.000   0.000  1.00  0.00           C",
        "ATOM      7  N   GLY A   6       6.000   0.000   0.000  1.00  0.00           N",
        "ATOM      8  CA  GLY A   6       7.000   0.000   0.000  1.00  0.00           C",
        "ATOM      9  C   GLY A   6       8.000   0.000   0.000  1.00  0.00           C",
        "ATOM     10  N   ALA B   1       0.000   5.000   0.000  1.00  0.00           N",
        "ATOM     11  CA  ALA B   1       1.000   5.000   0.000  1.00  0.00           C",
        "ATOM     12  C   ALA B   1       2.000   5.000   0.000  1.00  0.00           C",
        "ATOM     13  N   VAL B   9       3.000   5.000   0.000  1.00  0.00           N",
        "ATOM     14  CA  VAL B   9       4.000   5.000   0.000  1.00  0.00           C",
        "ATOM     15  C   VAL B   9       5.000   5.000   0.000  1.00  0.00           C",
        "END",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_seven_residue_gap_is_not_filled(tmp_path: Path):
    pdb = tmp_path / "gap.pdb"
    _write_gapped_pdb(pdb)
    report = classify_short_gaps(str(pdb))
    by_chain = {gap["chain"]: gap for gap in report["gaps"]}
    assert by_chain["A"]["count"] == 3
    assert by_chain["A"]["action"] == "fill"
    assert by_chain["A"]["sequence"] == ["ALA", "ALA", "ALA"]
    assert by_chain["B"]["count"] == 7
    assert by_chain["B"]["action"] == "skip"
    assert "longer than 6" in by_chain["B"]["reason"]
    assert SHORT_GAP_MAX == 6


def test_pdbfixer_builds_the_short_gap_and_leaves_the_long_one(tmp_path: Path):
    pytest.importorskip("pdbfixer")
    pdb = tmp_path / "gap.pdb"
    out = tmp_path / "filled.pdb"
    _write_gapped_pdb(pdb)
    report = fill_short_gaps(str(pdb), str(out), chains=["A", "B"])
    text = out.read_text(encoding="utf-8")
    assert "ALA A   3" in text or "ALA A   4" in text
    assert report["built"]
    assert all(item["count"] <= 6 for item in report["built"])
    assert any(gap["chain"] == "B" and gap["action"] == "skip" for gap in report["skipped"])


def test_custom_sequence_replaces_seqres(tmp_path: Path):
    pytest.importorskip("pdbfixer")
    pdb = tmp_path / "gap.pdb"
    out = tmp_path / "filled.pdb"
    _write_gapped_pdb(pdb)
    report = fill_short_gaps(
        str(pdb),
        str(out),
        chains=["A"],
        sequences=[
            {
                "chain": "A",
                "after_resid": 2,
                "before_resid": 6,
                "sequence": ["GLY", "GLY", "GLY"],
            }
        ],
    )
    text = out.read_text(encoding="utf-8")
    assert "GLY A   3" in text or "GLY A   4" in text
    assert report["built"][0]["sequence"] == ["GLY", "GLY", "GLY"]


def test_6rv2_chains_a_and_b_fill_only_short_gaps(tmp_path: Path):
    pytest.importorskip("pdbfixer")
    pdb = tmp_path / "6RV2.pdb"
    urllib.request.urlretrieve("https://files.rcsb.org/download/6RV2.pdb", pdb)
    out = tmp_path / "6RV2_AB.pdb"
    before = classify_short_gaps(str(pdb), chains=["A", "B"])
    report = fill_short_gaps(str(pdb), str(out), chains=["A", "B"])
    text = out.read_text(encoding="utf-8")
    assert "ATOM" in text
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
            assert line[21] in {"A", "B"}
    assert len(before["filled"]) == 2
    assert len(report["built"]) == 2
    assert all(item["count"] == 3 for item in report["built"])
    assert "MET A 149" in text and "MET B 149" in text
    for gap in before["skipped"]:
        if gap["count"] > 6:
            assert "longer than 6" in gap["reason"]
