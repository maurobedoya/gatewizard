"""Tests for peptide polymer residue recognition (gramicidin / D-aa / formyl)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gatewizard.core.builder import Builder
from gatewizard.core.preparation import PreparationManager
from gatewizard.tools.ligand_parametrization import detect_ligands
from gatewizard.utils.peptide_residues import (
    amber_unsupported_peptide_names_in_pdb,
    classify_polymer_kind,
    d_amino_acid_names_in_pdb,
    detect_peptide_caps_in_pdb,
    detect_peptide_caps_needing_gaff,
    is_peptide_polymer_residue,
    mda_peptide_or_protein_selection,
    remap_d_amino_acids_for_tleap,
    remap_pdb_peptide_resnames,
)
from gatewizard.utils.protein_capping import detect_terminal_caps


GRAMICIDIN_LIKE = """\
HEADER    ANTIBIOTIC
ATOM      1  N   GLY A   2      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  GLY A   2      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   GLY A   2      12.000  10.000  10.000  1.00 20.00           C
ATOM      4  O   GLY A   2      12.500  11.000  10.000  1.00 20.00           O
HETATM    5  N   FVA A   1       8.000  10.000  10.000  1.00 20.00           N
HETATM    6  CA  FVA A   1       9.000  10.000  10.000  1.00 20.00           C
HETATM    7  CN  FVA A   1       7.000  10.000  10.000  1.00 20.00           C
HETATM    8  O1  FVA A   1       6.500  11.000  10.000  1.00 20.00           O
HETATM    9  N   DLE A   3      13.000  10.000  10.000  1.00 20.00           N
HETATM   10  CA  DLE A   3      14.000  10.000  10.000  1.00 20.00           C
HETATM   11  C   DLE A   3      15.000  10.000  10.000  1.00 20.00           C
HETATM   12  O   DLE A   3      15.500  11.000  10.000  1.00 20.00           O
HETATM   13  N   DVA A   4      16.000  10.000  10.000  1.00 20.00           N
HETATM   14  CA  DVA A   4      17.000  10.000  10.000  1.00 20.00           C
HETATM   15  N   ETA A   5      18.000  10.000  10.000  1.00 20.00           N
HETATM   16  CA  ETA A   5      19.000  10.000  10.000  1.00 20.00           C
HETATM   17  C   ETA A   5      20.000  10.000  10.000  1.00 20.00           C
HETATM   18  O   ETA A   5      20.500  11.000  10.000  1.00 20.00           O
END
"""

PEPTIDE_WITH_SHEET = """\
HEADER    TEST PEPTIDE
SHEET    1   A 1 GLY A   2  TRP A   4  0
ATOM      1  N   FVA A   1      0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  FVA A   1      1.500   0.000   0.000  1.00 20.00           C
ATOM      3  C   FVA A   1      2.000   1.400   0.000  1.00 20.00           C
ATOM      4  O   FVA A   1      1.300   2.400   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2      3.400   1.500   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2      4.000   2.800   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2      5.500   2.800   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2      6.100   3.800   0.000  1.00 20.00           O
ATOM      9  N   DLE A   3      6.100   1.600   0.000  1.00 20.00           N
ATOM     10  CA  DLE A   3      7.500   1.500   0.000  1.00 20.00           C
ATOM     11  C   DLE A   3      8.100   0.100   0.000  1.00 20.00           C
ATOM     12  O   DLE A   3      7.500  -1.000   0.000  1.00 20.00           O
ATOM     13  N   TRP A   4      9.400   0.100   0.000  1.00 20.00           N
ATOM     14  CA  TRP A   4     10.100  -1.200   0.000  1.00 20.00           C
ATOM     15  C   TRP A   4     11.600  -1.100   0.000  1.00 20.00           C
ATOM     16  O   TRP A   4     12.200  -0.100   0.000  1.00 20.00           O
ATOM     17  N   ETA A   5     12.200  -2.300   0.000  1.00 20.00           N
ATOM     18  CA  ETA A   5     13.600  -2.400   0.000  1.00 20.00           C
END
"""


def test_is_peptide_polymer_residue_covers_gA_names():
    for name in ("GLY", "FVA", "DLE", "DVA", "ETA", "ACE", "NME"):
        assert is_peptide_polymer_residue(name)
    assert not is_peptide_polymer_residue("LIG")
    assert not is_peptide_polymer_residue("HOH")


def test_detect_terminal_caps_includes_formyl_and_eta(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    caps = detect_terminal_caps(pdb)
    assert "FVA" in caps
    assert "ETA" in caps
    assert detect_peptide_caps_in_pdb(str(pdb)) == caps


def test_detect_ligands_skips_peptide_polymer_hets(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    ligands = detect_ligands(str(pdb))
    names = {lig.name.upper() for lig in ligands}
    assert "FVA" not in names
    assert "DLE" not in names
    assert "DVA" not in names
    assert "ETA" not in names


def test_apply_protonation_empty_list_passthrough(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    out = tmp_path / "out.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    mgr = PreparationManager()
    result = mgr.apply_protonation_states(str(pdb), str(out), 7.0, {}, [])
    assert result["residue_changes"] == 0
    assert out.is_file()
    assert "FVA" in out.read_text(encoding="utf-8")


def test_amber_unsupported_and_remap(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    out = tmp_path / "mapped.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    warn = amber_unsupported_peptide_names_in_pdb(str(pdb))
    assert "FVA" in warn
    assert "ETA" in warn
    content = GRAMICIDIN_LIKE + (
        "HETATM   19  C   NMA A   6      21.000  10.000  10.000  1.00 20.00           C\n"
    )
    pdb.write_text(content, encoding="utf-8")
    info = remap_pdb_peptide_resnames(str(pdb), str(out))
    text = out.read_text(encoding="utf-8")
    assert "NME" in text
    assert info["record_changes"] >= 1


def test_remap_d_amino_acids_for_tleap_keeps_caps(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    out = tmp_path / "dremap.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    assert d_amino_acid_names_in_pdb(str(pdb)) == ["DLE", "DVA"]
    info = remap_d_amino_acids_for_tleap(str(pdb), str(out))
    text = out.read_text(encoding="utf-8")
    assert "LEU" in text
    assert "VAL" in text
    assert "DLE" not in text
    assert "DVA" not in text
    assert "FVA" in text
    assert "ETA" in text
    assert info["residue_changes"] >= 2
    caps = detect_peptide_caps_needing_gaff(str(pdb))
    assert {c["name"] for c in caps} == {"ETA", "FVA"}
    assert {c["role"] for c in caps} == {"c_term", "n_term"}


def test_classify_polymer_kind():
    assert classify_polymer_kind(["GLY", "ALA", "TRP"]) == "protein"
    assert classify_polymer_kind(["GLY", "FVA", "DLE"]) == "peptide"
    assert classify_polymer_kind(["ACE", "ALA", "NME"]) == "protein"


def test_mda_peptide_selection_string():
    sel = mda_peptide_or_protein_selection()
    assert sel.startswith("(protein or resname ")
    assert "FVA" in sel
    assert "DLE" in sel


def test_builder_command_includes_peptide_cap_flags(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    frcmod = tmp_path / "FVA.frcmod"
    lib = tmp_path / "FVA.lib"
    frcmod.write_text("x", encoding="utf-8")
    lib.write_text("x", encoding="utf-8")
    builder = Builder()
    cmd = builder._build_command(
        pdb,
        ["POPC"],
        ["POPC"],
        "1//1",
        {
            **builder.config,
            "peptide_cap_params": {
                "FVA": {"frcmod": str(frcmod), "lib": str(lib), "charge": 0},
                "ETA": {"frcmod": str(frcmod), "lib": str(lib), "charge": -1},
            },
        },
    )
    joined = " ".join(cmd)
    assert "--ligand_param" in cmd
    assert "--gaff2" in cmd
    assert "--charge_pdb_delta" in cmd
    assert "-1" in cmd
    assert "--keepligs" in cmd
    assert "--solute" not in cmd
    assert str(frcmod) in joined


def test_builder_validate_warns_missing_peptide_caps(tmp_path: Path):
    pdb = tmp_path / "ga.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    builder = Builder()
    ok, msg = builder.validate_system_inputs(
        str(pdb),
        ["POPC"],
        ["POPC"],
        "1//1",
        parametrize=True,
        peptide_cap_params={},
    )
    assert ok
    assert "FVA" in msg or "ETA" in msg
    assert "peptide" in msg.lower() or "GAFF" in msg


def test_ensure_cap_hydrogens_adds_h_to_heavy_atom_eta(tmp_path: Path):
    pytest.importorskip("rdkit")
    from gatewizard.tools.peptide_cap_parametrization import (
        _pdb_has_hydrogen_atoms,
        ensure_cap_hydrogens,
        extract_polymer_residue_pdb,
    )

    pdb = tmp_path / "ga.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    extracted = extract_polymer_residue_pdb(str(pdb), "ETA", str(tmp_path / "caps"))
    assert not _pdb_has_hydrogen_atoms(extracted)
    prepared = ensure_cap_hydrogens(extracted, "ETA")
    assert _pdb_has_hydrogen_atoms(prepared)
    text = Path(prepared).read_text(encoding="utf-8")
    assert "ETA" in text
    assert "UNL" not in text


def test_short_peptide_ss_prefers_pdb_sheet(tmp_path: Path):
    pytest.importorskip("psique")
    from gatewizard.core.structure_manager import assign_secondary_structure_map

    pdb = tmp_path / "sheet_peptide.pdb"
    pdb.write_text(PEPTIDE_WITH_SHEET, encoding="utf-8")
    ss = assign_secondary_structure_map(str(pdb), method="auto")
    mid = [ss.get(("A", i), "C") for i in (2, 3, 4)]
    assert any(code != "C" for code in mid)


def test_auto_detect_molecules_folds_peptide_into_protein(tmp_path: Path):
    pytest.importorskip("psique")
    from gatewizard.core.structure_manager import StructureManager

    pdb = tmp_path / "ga.pdb"
    pdb.write_text(GRAMICIDIN_LIKE, encoding="utf-8")
    sm = StructureManager()
    sm.load_structure(str(pdb))
    sels = sm.auto_detect_molecules()
    labels = [s.name for s in sels]
    assert "FVA" not in labels
    assert "DLE" not in labels
    assert "DVA" not in labels
    assert "ETA" not in labels
    assert "Peptide" in labels
    assert "Protein" not in labels
    polymer = next(s for s in sels if s.name == "Peptide")
    assert len(polymer.atom_indices) == 18
