"""Biological assembly from BIOMT / mmCIF operators."""

from __future__ import annotations

import urllib.error
import urllib.request
from pathlib import Path

import gemmi
import pytest

from gatewizard.core.assembly import (
    LARGE_ASSEMBLY_ATOM_COUNT,
    build_assembly,
    is_biological_assembly_mmcif,
    list_assemblies,
    universe_from_assembly_mmcif,
)

_DIMER = """\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.000   0.000   0.000  1.00  0.00           C
HETATM    3  O   HOH A   2       5.000   0.000   0.000  1.00  0.00           O
HETATM    4  C   LIG A   3       8.000   0.000   0.000  1.00  0.00           C
END
"""

_HETEROMER = """\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: TETRAMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A, B
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CA  GLY B   1       3.000   0.000   0.000  1.00  0.00           C
END
"""

_IDENTITY = """\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: MONOMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
END
"""

_BARE = """\
ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C
END
"""

# 1FAT's current deposit is already the tetramer (one identity operator).
# 1HHO's asymmetric unit is the αβ dimer; the assembly is the α2β2 tetramer.
# 1A34 is a 60-mer of a three-chain asymmetric unit, so the assembly file
# has 180 chains, and it is large enough to need confirmation.
_RCSB = ("1FAT", "1HHO", "1A34")


def _write(tmp_path: Path, text: str, name: str = "unit.pdb") -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _nonwater_chains(path: Path) -> int:
    structure = gemmi.read_structure(str(path))
    count = 0
    for chain in structure[0]:
        if any(str(residue.name).strip().upper() not in {"HOH", "WAT", "DOD"} for residue in chain):
            count += 1
    return count


def _download(url: str, dest: Path) -> None:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            dest.write_bytes(response.read())
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        pytest.skip(f"RCSB is not reachable: {exc}")


def test_dimer_drops_waters_and_keeps_numbered_chains(tmp_path: Path) -> None:
    path = _write(tmp_path, _DIMER)
    listed = list_assemblies(path)
    assert listed["message"] == ""
    row = listed["assemblies"][0]
    assert row["id"] == "1"
    assert row["oligomeric_count"] == 2
    assert row["atom_count"] == 6  # N, CA, LIG on each copy
    assert row["atom_count_with_waters"] == 8
    assert row["already_biological_unit"] is False
    assert "author" in row["description"].lower()

    built = build_assembly(path, "1", include_waters=False)
    assert built["built"] is True
    assert built["chains"] == ["A1", "A2"]
    assert built["n_atoms"] == 6
    assert is_biological_assembly_mmcif(built["path"])
    text = Path(built["path"]).read_text(encoding="utf-8")
    assert "HOH" not in text
    assert "LIG" in text

    universe = universe_from_assembly_mmcif(built["path"])
    assert list(universe.atoms.chainIDs) == ["A1", "A1", "A1", "A2", "A2", "A2"]
    # Copy A2 is the 180 degree operator. CA was at x=1, so it lands at x=-1.
    assert float(universe.atoms.positions[4][0]) == pytest.approx(-1.0)

    wet = build_assembly(path, "1", include_waters=True)
    assert wet["n_atoms"] == 8


_AXIS_ION = """\
REMARK 350 BIOMOLECULE: 1
REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC
REMARK 350 APPLY THE FOLLOWING TO CHAINS: A
REMARK 350   BIOMT1   1  1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   1  0.000000  1.000000  0.000000        0.00000
REMARK 350   BIOMT3   1  0.000000  0.000000  1.000000        0.00000
REMARK 350   BIOMT1   2 -1.000000  0.000000  0.000000        0.00000
REMARK 350   BIOMT2   2  0.000000 -1.000000  0.000000        0.00000
REMARK 350   BIOMT3   2  0.000000  0.000000  1.000000        0.00000
ATOM      1  CA  ALA A   1       5.000   0.000   0.000  1.00  0.00           C
HETATM    2  K     K A   2       0.000   0.000   0.000  1.00  0.00           K
END
"""


def test_symmetry_axis_ion_is_kept_once(tmp_path: Path) -> None:
    path = _write(tmp_path, _AXIS_ION)
    row = list_assemblies(path)["assemblies"][0]
    assert row["stacked_extra_atoms"] == 1
    assert row["stacked_detail"] == "K × 1"
    assert row["protein_overlap_chains"] == []

    built = build_assembly(path, "1")
    assert built["n_atoms"] == 3
    assert "Removed 1 extra atom" in built["message"]
    structure = gemmi.read_structure(built["path"])
    ions = [
        str(chain.name)
        for chain in structure[0]
        for residue in chain
        if str(residue.name).strip() == "K"
    ]
    assert ions == ["A1"]

    kept = build_assembly(path, "1", drop_overlaps=False)
    assert kept["n_atoms"] == 4
    assert kept["message"] == ""


def test_heteromer_copies_both_chains(tmp_path: Path) -> None:
    path = _write(tmp_path, _HETEROMER)
    built = build_assembly(path, "1")
    assert built["built"] is True
    assert len(built["chains"]) == 4
    assert set(built["chains"]) == {"A1", "A2", "B1", "B2"}


def test_identity_operator_does_not_duplicate(tmp_path: Path) -> None:
    path = _write(tmp_path, _IDENTITY)
    listed = list_assemblies(path)
    assert "nothing to build" in listed["message"]
    assert listed["assemblies"][0]["already_biological_unit"] is True
    built = build_assembly(path, "1")
    assert built["built"] is False
    assert built["path"] == ""
    assert "nothing to build" in built["message"]


def test_file_without_operators_reports_nothing_to_build(tmp_path: Path) -> None:
    path = _write(tmp_path, _BARE)
    listed = list_assemblies(path)
    assert listed["assemblies"] == []
    assert "nothing to build" in listed["message"]


@pytest.mark.parametrize("pdb_id", _RCSB)
def test_rcsb_assembly_chain_count(tmp_path: Path, pdb_id: str) -> None:
    pdb_path = tmp_path / f"{pdb_id}.pdb"
    official_path = tmp_path / f"{pdb_id}-assembly1.cif"
    code = pdb_id.lower()
    _download(f"https://files.rcsb.org/download/{code}.pdb", pdb_path)
    _download(f"https://files.rcsb.org/download/{code}-assembly1.cif", official_path)

    official_chains = _nonwater_chains(official_path)
    listed = list_assemblies(pdb_path)
    matches = [
        row
        for row in listed["assemblies"]
        if row["id"] != "ncs" and row["oligomeric_count"] == official_chains
    ]
    assert matches, listed
    row = matches[0]
    if pdb_id == "1A34":
        assert official_chains == 180
        assert row["atom_count"] >= LARGE_ASSEMBLY_ATOM_COUNT
        assert row["needs_confirmation"] is True
    if row["already_biological_unit"]:
        assert "nothing to build" in listed["message"]
        built = build_assembly(pdb_path, row["id"])
        assert built["built"] is False
        assert _nonwater_chains(pdb_path) == official_chains
        return
    built = build_assembly(pdb_path, row["id"], include_waters=False)
    assert built["built"] is True
    assert len(built["chains"]) == official_chains


def test_rcsb_monomer_has_nothing_to_build(tmp_path: Path) -> None:
    pdb_path = tmp_path / "1CRN.pdb"
    _download("https://files.rcsb.org/download/1crn.pdb", pdb_path)
    listed = list_assemblies(pdb_path)
    assert "nothing to build" in listed["message"]
    for row in listed["assemblies"]:
        built = build_assembly(pdb_path, row["id"])
        assert built["built"] is False
        assert "nothing to build" in built["message"]
