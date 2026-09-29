"""Tests for restoring original residue numbers after pdb4amber / capping."""

from pathlib import Path

from gatewizard.core.preparation import PreparationManager
from gatewizard.utils.residue_restore import restore_original_residue_numbers


def _atom(serial: int, name: str, resname: str, chain: str, resid: int, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  {name:<3s} {resname:>3s} {chain}{resid:4d}    "
        f"{x:8.3f}{x:8.3f}{x:8.3f}  1.00  0.00           {name[0]}\n"
    )


def _het(serial: int, name: str, resname: str, chain: str, resid: int, x: float) -> str:
    return (
        f"HETATM{serial:5d}  {name:<3s} {resname:>3s} {chain}{resid:4d}    "
        f"{x:8.3f}{x:8.3f}{x:8.3f}  1.00  0.00           {name[0]}\n"
    )


def _write_pdb(path: Path, lines: list[str]) -> Path:
    path.write_text("".join(lines) + "END\n", encoding="utf-8")
    return path


def _resids(path: Path) -> list[tuple[str, str, int]]:
    out = []
    last = None
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[17:20].strip(), line[21:22], int(line[22:26]))
        if key != last:
            out.append(key)
            last = key
    return out


def test_restore_gap_two_chains_keeps_hetero(tmp_path: Path):
    sequential = _write_pdb(
        tmp_path / "seq.pdb",
        [
            _atom(1, "CA", "ALA", "A", 1, 1.0),
            _atom(2, "CA", "GLY", "A", 2, 2.0),
            _het(3, "C1", "LIG", "A", 3, 3.0),
            _atom(4, "CA", "VAL", "B", 4, 4.0),
            _atom(5, "CA", "SER", "B", 5, 5.0),
        ],
    )
    remum = tmp_path / "seq_renum.txt"
    remum.write_text(
        "ALA A      1    ALA    10\n"
        "GLY A      2    GLY    20\n"
        "LIG A      3    LIG   900\n"
        "VAL B      4    VAL    10\n"
        "SER B      5    SER    20\n",
        encoding="utf-8",
    )
    out = tmp_path / "restored.pdb"
    result = restore_original_residue_numbers(
        sequential, remum_path=remum, output_pdb=out
    )
    assert result["residue_numbers_preserved"] is True
    assert _resids(out) == [
        ("ALA", "A", 10),
        ("GLY", "A", 20),
        ("LIG", "A", 900),
        ("VAL", "B", 10),
        ("SER", "B", 20),
    ]
    hetero = [
        line
        for line in out.read_text(encoding="utf-8").splitlines()
        if line.startswith("HETATM")
    ]
    assert len(hetero) == 1


def test_restore_caps_protein_starts_at_1(tmp_path: Path):
    sequential = _write_pdb(
        tmp_path / "capped.pdb",
        [
            _atom(1, "C", "ACE", "A", 1, 0.0),
            _atom(2, "CA", "MET", "A", 2, 1.0),
            _atom(3, "CA", "ALA", "A", 3, 2.0),
            _atom(4, "N", "NME", "A", 4, 3.0),
        ],
    )
    remum = tmp_path / "capped_renum.txt"
    remum.write_text(
        "ACE A      1    ACE     1\n"
        "MET A      2    MET     2\n"
        "ALA A      3    ALA     3\n"
        "NME A      4    NME     4\n",
        encoding="utf-8",
    )
    cap_map = tmp_path / "capped_gatewizard_residue_mapping.txt"
    cap_map.write_text(
        "# Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID\n"
        "  - A      -    ACE     1\n"
        "MET A      1    MET     2\n"
        "ALA A      2    ALA     3\n"
        "  - A      -    NME     4\n",
        encoding="utf-8",
    )
    out = tmp_path / "restored.pdb"
    result = restore_original_residue_numbers(
        sequential,
        remum_path=remum,
        output_pdb=out,
        cap_mapping_path=cap_map,
    )
    assert _resids(out) == [
        ("ACE", "A", 0),
        ("MET", "A", 1),
        ("ALA", "A", 2),
        ("NME", "A", 3),
    ]
    assert result["cap_assignments"]["A:ACE"] == 0
    assert result["cap_assignments"]["A:NME"] == 3
    keys = {(r, c, i) for r, c, i in _resids(out)}
    assert len(keys) == 4


def test_restore_caps_protein_starts_at_21(tmp_path: Path):
    sequential = _write_pdb(
        tmp_path / "capped.pdb",
        [
            _atom(1, "C", "ACE", "A", 1, 0.0),
            _atom(2, "CA", "VAL", "A", 2, 1.0),
            _atom(3, "CA", "ALA", "A", 3, 2.0),
            _atom(4, "N", "NME", "A", 4, 3.0),
        ],
    )
    remum = tmp_path / "capped_renum.txt"
    remum.write_text(
        "ACE A      1    ACE     1\n"
        "VAL A      2    VAL     2\n"
        "ALA A      3    ALA     3\n"
        "NME A      4    NME     4\n",
        encoding="utf-8",
    )
    cap_map = tmp_path / "capped_gatewizard_residue_mapping.txt"
    cap_map.write_text(
        "# Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID\n"
        "  - A      -    ACE     1\n"
        "VAL A     21    VAL     2\n"
        "ALA A     22    ALA     3\n"
        "  - A      -    NME     4\n",
        encoding="utf-8",
    )
    out = tmp_path / "restored.pdb"
    restore_original_residue_numbers(
        sequential,
        remum_path=remum,
        output_pdb=out,
        cap_mapping_path=cap_map,
    )
    assert _resids(out) == [
        ("ACE", "A", 20),
        ("VAL", "A", 21),
        ("ALA", "A", 22),
        ("NME", "A", 23),
    ]


def test_restore_cap_plus_gap_compose(tmp_path: Path):
    sequential = _write_pdb(
        tmp_path / "capped.pdb",
        [
            _atom(1, "C", "ACE", "A", 1, 0.0),
            _atom(2, "CA", "ALA", "A", 2, 1.0),
            _atom(3, "CA", "GLY", "A", 3, 2.0),
            _atom(4, "N", "NME", "A", 4, 3.0),
            _het(5, "C1", "LIG", "A", 5, 4.0),
        ],
    )
    remum = tmp_path / "capped_renum.txt"
    remum.write_text(
        "ACE A      1    ACE     1\n"
        "ALA A      2    ALA     2\n"
        "GLY A      3    GLY     3\n"
        "NME A      4    NME     4\n"
        "LIG A      5    LIG     5\n",
        encoding="utf-8",
    )
    cap_map = tmp_path / "capped_gatewizard_residue_mapping.txt"
    cap_map.write_text(
        "# Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID\n"
        "  - A      -    ACE     1\n"
        "ALA A    200    ALA     2\n"
        "GLY A    205    GLY     3\n"
        "  - A      -    NME     4\n"
        "LIG A    900    LIG     5\n",
        encoding="utf-8",
    )
    out = tmp_path / "restored.pdb"
    restore_original_residue_numbers(
        sequential,
        remum_path=remum,
        output_pdb=out,
        cap_mapping_path=cap_map,
    )
    assert _resids(out) == [
        ("ACE", "A", 199),
        ("ALA", "A", 200),
        ("GLY", "A", 205),
        ("NME", "A", 206),
        ("LIG", "A", 900),
    ]
    resnames_at_break = [r for r, _c, i in _resids(out) if i in {200, 205}]
    assert resnames_at_break == ["ALA", "GLY"]


def test_restore_merged_chains_from_original_pdb(tmp_path: Path):
    """tleap/pdb4amber blank the chain and concatenate A then B."""
    original = _write_pdb(
        tmp_path / "9g9v.pdb",
        [
            _atom(1, "CA", "MET", "A", 1, 1.0),
            _atom(2, "CA", "LYS", "A", 2, 2.0),
            _atom(3, "CA", "MET", "B", 1, 3.0),
            _atom(4, "CA", "LYS", "B", 2, 4.0),
            _het(5, "C1", "Y01", "A", 301, 5.0),
            _het(6, "K", "K", "A", 303, 6.0),
        ],
    )
    sequential = _write_pdb(
        tmp_path / "seq.pdb",
        [
            _atom(1, "CA", "MET", " ", 1, 1.0),
            _atom(2, "CA", "LYS", " ", 2, 2.0),
            _atom(3, "CA", "MET", " ", 3, 3.0),
            _atom(4, "CA", "LYS", " ", 4, 4.0),
            _het(5, "C1", "Y01", "A", 5, 5.0),
            _het(6, "K", "K", "A", 6, 6.0),
        ],
    )
    remum = tmp_path / "seq_renum.txt"
    remum.write_text(
        "MET       1    MET     1\n"
        "LYS       2    LYS     2\n"
        "MET       3    MET     3\n"
        "LYS       4    LYS     4\n"
        "Y01 A   301    Y01     5\n"
        "  K A   303      K     6\n",
        encoding="utf-8",
    )
    out = tmp_path / "restored.pdb"
    result = restore_original_residue_numbers(
        sequential,
        remum_path=remum,
        output_pdb=out,
        original_pdb=original,
    )
    assert result["residue_numbers_preserved"] is True
    assert _resids(out) == [
        ("MET", "A", 1),
        ("LYS", "A", 2),
        ("MET", "B", 1),
        ("LYS", "B", 2),
        ("Y01", "A", 301),
        ("K", "A", 303),
    ]
    assert result["atoms_rewritten"] >= 4


def test_manager_restore_delegates(tmp_path: Path):
    sequential = _write_pdb(
        tmp_path / "seq.pdb",
        [_atom(1, "CA", "ALA", "A", 1, 1.0), _atom(2, "CA", "GLY", "A", 2, 2.0)],
    )
    remum = tmp_path / "seq_renum.txt"
    remum.write_text(
        "ALA A      1    ALA    10\nGLY A      2    GLY    20\n",
        encoding="utf-8",
    )
    manager = PreparationManager()
    out = tmp_path / "out.pdb"
    info = manager.restore_original_residue_numbers(
        str(sequential), remum_path=str(remum), output_pdb=str(out)
    )
    assert info["residue_numbers_preserved"] is True
    assert _resids(out) == [("ALA", "A", 10), ("GLY", "A", 20)]
