"""Structure checker: missing atoms, gaps, altLoc, and chain differences."""

from pathlib import Path

from gatewizard.core.structure_checker import (
    apply_altloc_choices,
    diagnose_structure,
    keep_most_likely_altlocs,
)
from gatewizard.tools.ligand_parametrization import detect_ligands, extract_ligand_pdb


def _atom(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float = 0.0,
    z: float = 0.0,
    occ: float = 1.0,
    alt: str = " ",
) -> str:
    atom = name if len(name) == 4 else f" {name:<3s}"
    return (
        f"ATOM  {serial:5d} {atom}{alt:1s}{resname:>3s} {chain}"
        f"{resid:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{0.0:6.2f}           {name[0]:>2s}\n"
    )


def _sample(path: Path) -> None:
    lines = [
        _atom(1, "N", "ALA", "A", 1, 0.0),
        _atom(2, "CA", "ALA", "A", 1, 1.0),
        _atom(3, "C", "ALA", "A", 1, 2.0),
        _atom(4, "O", "ALA", "A", 1, 2.5),
        _atom(5, "CB", "ALA", "A", 1, 1.2, occ=0.70, alt="A"),
        _atom(6, "CB", "ALA", "A", 1, 1.4, occ=0.30, alt="B"),
        _atom(7, "N", "GLY", "A", 2, 3.0),
        _atom(8, "CA", "GLY", "A", 2, 4.0),
        _atom(9, "C", "GLY", "A", 2, 0.0),
        _atom(10, "O", "GLY", "A", 2, 0.5),
        _atom(11, "N", "ALA", "A", 3, 10.0),
        _atom(12, "CA", "ALA", "A", 3, 11.0),
        _atom(13, "C", "ALA", "A", 3, 12.0),
        _atom(14, "O", "ALA", "A", 3, 12.5),
        _atom(15, "N", "ALA", "A", 5, 20.0),
        _atom(16, "CA", "ALA", "A", 5, 21.0),
        _atom(17, "C", "ALA", "A", 5, 22.0),
        _atom(18, "O", "ALA", "A", 5, 22.5),
        _atom(19, "CB", "ALA", "A", 5, 21.2),
        _atom(20, "N", "ALA", "B", 1, 0.0, 5.0),
        _atom(21, "CA", "ALA", "B", 1, 1.0, 5.0),
        _atom(22, "C", "ALA", "B", 1, 2.0, 5.0),
        _atom(23, "O", "ALA", "B", 1, 2.5, 5.0),
        _atom(24, "CB", "ALA", "B", 1, 1.2, 5.0),
        _atom(25, "N", "VAL", "B", 2, 3.0, 5.0),
        _atom(26, "CA", "VAL", "B", 2, 4.0, 5.0),
        _atom(27, "C", "VAL", "B", 2, 5.0, 5.0),
        _atom(28, "O", "VAL", "B", 2, 5.5, 5.0),
        _atom(29, "CB", "VAL", "B", 2, 4.2, 5.0),
        _atom(30, "CG1", "VAL", "B", 2, 4.4, 5.0),
        _atom(31, "CG2", "VAL", "B", 2, 4.6, 5.0),
        _atom(32, "N", "ALA", "B", 3, 6.0, 5.0),
        _atom(33, "CA", "ALA", "B", 3, 7.0, 5.0),
        _atom(34, "C", "ALA", "B", 3, 8.0, 5.0),
        _atom(35, "O", "ALA", "B", 3, 8.5, 5.0),
        _atom(36, "CB", "ALA", "B", 3, 7.2, 5.0),
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


def test_diagnose_reports_missing_gap_altloc_and_chains(tmp_path: Path):
    src = tmp_path / "sample.pdb"
    _sample(src)
    report = diagnose_structure(str(src))
    missing = {(item["chain"], item["resid"], item["atom"]) for item in report["missing_atoms"]}
    assert ("A", 3, "CB") in missing
    assert ("A", 1, "CB") not in missing
    numbering = [item for item in report["missing_residues"] if item.get("count") == 1]
    assert numbering and numbering[0]["missing"] == [4]
    assert report["gaps"][0]["resid_a"] == 2
    assert report["gaps"][0]["resid_b"] == 3
    assert report["gaps"][0]["distance"] > 2
    group = report["altlocs"][0]
    assert group["resname"] == "ALA"
    assert group["resid"] == 1
    assert "atom" not in group
    assert group["default"] == "A"
    assert group["selection"] == "resname ALA and chainID A and resid 1"
    assert group["options"][0]["likely"] is True
    pair = report["chain_pairs"][0]
    assert pair["same"] is False
    assert any(item["resname_a"] == "GLY" and item["resname_b"] == "VAL" for item in pair["substitutions"])
    assert "longer" in pair["length_note"]
    assert any(item["kind"] == "hydrogens" and item["selection"] == "protein" for item in report["issues"])


def test_apply_keeps_the_most_likely_altloc_unless_chosen(tmp_path: Path):
    src = tmp_path / "sample.pdb"
    _sample(src)
    default_out = tmp_path / "default_checked.pdb"
    kept = apply_altloc_choices(str(src), output_path=str(default_out))
    assert kept["atoms_removed"] == 1
    text = default_out.read_text(encoding="utf-8")
    cb = [line for line in text.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CB" and line[21] == "A" and int(line[22:26]) == 1]
    assert len(cb) == 1
    assert cb[0][16] == " "
    chosen = tmp_path / "chosen.pdb"
    apply_altloc_choices(
        str(src),
        choices=[{"chain": "A", "resid": 1, "icode": "", "altloc": "B"}],
        output_path=str(chosen),
    )
    chosen_text = chosen.read_text(encoding="utf-8")
    chosen_cb = [
        line
        for line in chosen_text.splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CB" and line[21] == "A" and int(line[22:26]) == 1
    ]
    assert len(chosen_cb) == 1
    assert float(chosen_cb[0][54:60]) == 1.0
    again = diagnose_structure(str(chosen))
    assert again["altlocs"] == []


def test_keep_most_likely_copies_a_file_without_altlocs(tmp_path: Path):
    src = tmp_path / "plain.pdb"
    src.write_text(_atom(1, "N", "ALA", "A", 1, 0.0), encoding="utf-8")
    dest = tmp_path / "plain_out.pdb"
    kept = keep_most_likely_altlocs(str(src), str(dest))
    assert kept["residues"] == 0
    assert kept["atoms_removed"] == 0
    assert dest.read_bytes() == src.read_bytes()


def _hetatm(
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    occ: float = 1.0,
    alt: str = " ",
) -> str:
    atom = name if len(name) == 4 else f" {name:<3s}"
    return (
        f"HETATM{serial:5d} {atom}{alt:1s}{resname:>3s} {chain}"
        f"{resid:4d}    {x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{0.0:6.2f}           {name[0]:>2s}\n"
    )


def _shared_number(path: Path) -> None:
    """Protein and a ligand share chain A residue 1. A water has its own number."""
    lines = [
        "SEQRES   1 A    1  ALA\n",
        _atom(1, "N", "ALA", "A", 1, 0.0),
        _atom(2, "CA", "ALA", "A", 1, 1.0, occ=0.70, alt="A"),
        _atom(3, "CA", "ALA", "A", 1, 1.1, occ=0.30, alt="B"),
        _hetatm(4, "C1", "RET", "A", 1, 5.0, occ=0.20, alt="A"),
        _hetatm(5, "C1", "RET", "A", 1, 5.2, occ=0.80, alt="B"),
        _hetatm(6, "O", "HOH", "A", 9, 8.0, occ=0.40, alt="A"),
        _hetatm(7, "O", "HOH", "A", 9, 8.2, occ=0.60, alt="B"),
        "LINK         C1  BRET A   1                 CA  AALA A   1\n",
        "END\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")


def test_hetero_altlocs_are_their_own_residues(tmp_path: Path):
    src = tmp_path / "shared.pdb"
    _shared_number(src)
    report = diagnose_structure(str(src))
    found = {(group["resname"], group["default"], group["selection"]) for group in report["altlocs"]}
    assert ("ALA", "A", "resname ALA and chainID A and resid 1") in found
    assert ("RET", "B", "resname RET and chainID A and resid 1") in found
    assert ("HOH", "B", "resname HOH and chainID A and resid 9") in found
    dest = tmp_path / "likely.pdb"
    kept = keep_most_likely_altlocs(str(src), str(dest))
    assert kept["residues"] == 3
    assert kept["atoms_removed"] == 3
    text = dest.read_text(encoding="utf-8")
    assert "SEQRES" in text
    assert "LINK" in text
    ret = [line for line in text.splitlines() if line.startswith("HETATM") and line[17:20].strip() == "RET"]
    ala = [line for line in text.splitlines() if line.startswith("ATOM") and line[12:16].strip() == "CA"]
    hoh = [line for line in text.splitlines() if line.startswith("HETATM") and line[17:20].strip() == "HOH"]
    assert len(ret) == 1 and ret[0][16] == " " and float(ret[0][30:38]) == 5.2
    assert len(ala) == 1 and float(ala[0][30:38]) == 1.0
    assert len(hoh) == 1 and float(hoh[0][30:38]) == 8.2
    chosen = tmp_path / "ligand_a.pdb"
    apply_altloc_choices(
        str(src),
        choices=[{"chain": "A", "resid": 1, "icode": "", "resname": "RET", "altloc": "A"}],
        output_path=str(chosen),
    )
    chosen_text = chosen.read_text(encoding="utf-8")
    chosen_ret = [
        line for line in chosen_text.splitlines()
        if line.startswith("HETATM") and line[17:20].strip() == "RET"
    ]
    chosen_ala = [
        line for line in chosen_text.splitlines()
        if line.startswith("ATOM") and line[12:16].strip() == "CA"
    ]
    assert len(chosen_ret) == 1 and float(chosen_ret[0][30:38]) == 5.0
    assert len(chosen_ala) == 1 and float(chosen_ala[0][30:38]) == 1.0
    assert diagnose_structure(str(chosen))["altlocs"] == []


def test_keep_most_likely_drops_the_lower_occupancy(tmp_path: Path):
    src = tmp_path / "sample.pdb"
    _sample(src)
    dest = tmp_path / "likely.pdb"
    kept = keep_most_likely_altlocs(str(src), str(dest))
    assert kept["residues"] == 1
    assert kept["atoms_removed"] == 1
    assert diagnose_structure(str(dest))["altlocs"] == []


def _lig(path: Path, lines: list[str]) -> None:
    path.write_text("".join(lines) + "END\n", encoding="utf-8")


def test_overlapping_atoms_without_altloc_are_conformations(tmp_path: Path):
    src = tmp_path / "overlap.pdb"
    _lig(
        src,
        [
            _hetatm(1, "C1", "LIG", "A", 1, 0.0),
            _hetatm(2, "C2", "LIG", "A", 1, 1.4),
            _hetatm(3, "C1", "LIG", "A", 1, 0.4),
            _hetatm(4, "C2", "LIG", "A", 1, 1.7),
        ],
    )
    report = diagnose_structure(str(src))
    assert len(report["altlocs"]) == 1
    group = report["altlocs"][0]
    assert group["inferred"] is True
    assert group["resname"] == "LIG"
    assert group["default"] == "A"
    assert [option["altloc"] for option in group["options"]] == ["A", "B"]
    dest = tmp_path / "one.pdb"
    kept = keep_most_likely_altlocs(str(src), str(dest))
    assert kept["atoms_removed"] == 2
    coords = [
        float(line[30:38])
        for line in dest.read_text(encoding="utf-8").splitlines()
        if line.startswith("HETATM")
    ]
    assert coords == [0.0, 1.4]
    chosen = tmp_path / "pose_b.pdb"
    apply_altloc_choices(
        str(src),
        choices=[{"chain": "A", "resid": 1, "icode": "", "resname": "LIG", "altloc": "B"}],
        output_path=str(chosen),
    )
    chosen_coords = [
        float(line[30:38])
        for line in chosen.read_text(encoding="utf-8").splitlines()
        if line.startswith("HETATM")
    ]
    assert chosen_coords == [0.4, 1.7]
    ligands = detect_ligands(str(src))
    assert len(ligands) == 1
    assert ligands[0].num_atoms == 2
    assert ligands[0].conformers == 2
    extracted = extract_ligand_pdb(str(src), "LIG", str(tmp_path / "lig"))
    extracted_coords = [
        float(line[30:38])
        for line in Path(extracted).read_text(encoding="utf-8").splitlines()
        if line.startswith("HETATM")
    ]
    assert extracted_coords == [0.0, 1.4]


def test_separate_ligand_copies_that_overlap_are_one_conformation_pair(tmp_path: Path):
    src = tmp_path / "copies.pdb"
    _lig(
        src,
        [
            _hetatm(1, "C1", "LIG", "A", 1, 0.0),
            _hetatm(2, "C2", "LIG", "A", 1, 1.4),
            _hetatm(3, "C1", "LIG", "A", 2, 0.4),
            _hetatm(4, "C2", "LIG", "A", 2, 1.7),
            _hetatm(5, "C1", "LIG", "A", 9, 40.0),
            _hetatm(6, "C2", "LIG", "A", 9, 41.4),
        ],
    )
    report = diagnose_structure(str(src))
    assert len(report["altlocs"]) == 1
    group = report["altlocs"][0]
    assert group["mode"] == "copies"
    assert group["inferred"] is True
    assert "LIG A1" in group["label"]
    assert "LIG A2" in group["label"]
    dest = tmp_path / "kept.pdb"
    keep_most_likely_altlocs(str(src), str(dest))
    text = dest.read_text(encoding="utf-8")
    resids = [int(line[22:26]) for line in text.splitlines() if line.startswith("HETATM")]
    assert resids == [1, 1, 9, 9]


def test_close_potassium_ions_are_not_alternate_conformations(tmp_path: Path):
    src = tmp_path / "ions.pdb"
    path = src
    path.write_text(
        "".join(
            [
                _atom(1, "K+", "K+", "A", 1017, 0.0),
                _atom(2, "K+", "K+", "A", 1033, 0.4),
                _atom(3, "K+", "K+", "A", 1100, 40.0),
            ]
        )
        + "END\n",
        encoding="utf-8",
    )
    report = diagnose_structure(str(src))
    assert report["altlocs"] == []
    assert not any(item.get("kind") == "unknown" and "K+" in item.get("message", "") for item in report["issues"])


def test_separated_or_identical_copies_are_not_alternate_conformations(tmp_path: Path):
    far = tmp_path / "far.pdb"
    _lig(
        far,
        [
            _hetatm(1, "C1", "LIG", "A", 1, 0.0),
            _hetatm(2, "C2", "LIG", "A", 1, 1.4),
            _hetatm(3, "C1", "LIG", "A", 1, 30.0),
            _hetatm(4, "C2", "LIG", "A", 1, 31.4),
        ],
    )
    assert diagnose_structure(str(far))["altlocs"] == []
    same = tmp_path / "same.pdb"
    line = _hetatm(1, "C1", "LIG", "A", 1, 0.0)
    _lig(same, [line, line])
    report = diagnose_structure(str(same))
    assert report["altlocs"] == []
    assert any(item["kind"] == "duplicate" for item in report["issues"])
