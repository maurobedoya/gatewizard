"""Tests for Maestro CT property reading."""

from __future__ import annotations

from pathlib import Path

import pytest

from gatewizard.tools.maestro import (
    HEADLINE_PROPERTY_NAMES,
    MaestroError,
    is_maestro_path,
    list_maestro_entries,
    maestro_ct_properties,
    pick_headline_property,
    read_maestro_properties,
)


def _write_mae_with_props(path: Path) -> None:
    from rdkit import Chem

    mol = Chem.AddHs(Chem.MolFromSmiles("C=O"))
    assert mol is not None
    mol.SetProp("_Name", "pose1")
    mol.SetDoubleProp("r_i_docking_score", -8.5)
    mol.SetDoubleProp("r_psp_MMGBSA_dG_Bind", -42.3)
    mol.SetIntProp("i_i_glide_lignum", 3)
    mol.SetProp("s_user_note", "ok")
    mol.SetProp("s_m_pdb_file", "X" * 600)
    writer = Chem.MaeWriter(str(path))
    writer.write(mol)
    writer.close()


@pytest.fixture()
def scored_mae(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem

    if not hasattr(Chem, "MaeWriter") or not hasattr(Chem, "MaeMolSupplier"):
        pytest.skip("RDKit Maestro I/O not available")
    mae = tmp_path / "scored.mae"
    _write_mae_with_props(mae)
    return mae


def test_is_maestro_path():
    assert is_maestro_path("a.mae")
    assert is_maestro_path("a.maegz")
    assert is_maestro_path("a.mae.gz")
    assert not is_maestro_path("a.pdb")


def test_maestro_ct_properties_from_mol(scored_mae):
    from rdkit import Chem

    mol = next(iter(Chem.MaeMolSupplier(str(scored_mae), sanitize=False, removeHs=False)))
    props = {p["name"]: p for p in maestro_ct_properties(mol)}
    assert props["r_i_docking_score"]["value"] == pytest.approx(-8.5)
    assert props["r_psp_MMGBSA_dG_Bind"]["value"] == pytest.approx(-42.3)
    assert props["i_i_glide_lignum"] == {
        "name": "i_i_glide_lignum",
        "type": "integer",
        "value": 3,
    }
    assert props["s_user_note"]["value"] == "ok"
    assert "s_m_pdb_file" not in props
    assert "_Name" not in props


def test_list_and_read_maestro_properties(scored_mae):
    info = list_maestro_entries(scored_mae)
    assert info["sourcePath"] == str(scored_mae.resolve())
    assert len(info["entries"]) == 1
    entry = info["entries"][0]
    assert entry["index"] == 0
    assert entry["title"] == "pose1"
    assert "pose1" in entry["label"]
    names = {p["name"] for p in entry["properties"]}
    assert "r_psp_MMGBSA_dG_Bind" in names
    assert "s_m_pdb_file" not in names

    rows = read_maestro_properties(scored_mae, [0])
    assert len(rows) == 1
    assert rows[0]["index"] == 0
    assert any(p["name"] == "r_i_docking_score" for p in rows[0]["properties"])


def test_read_maestro_properties_rejects_bad_index(scored_mae):
    with pytest.raises(ValueError, match="out of range"):
        read_maestro_properties(scored_mae, [99])


def test_pick_headline_property_prefers_mmgbsa():
    props = [
        {"name": "r_i_glide_gscore", "type": "real", "value": -7.0},
        {"name": "r_i_docking_score", "type": "real", "value": -8.0},
        {"name": "r_psp_MMGBSA_dG_Bind", "type": "real", "value": -42.3},
    ]
    assert pick_headline_property(props)["name"] == HEADLINE_PROPERTY_NAMES[0]
    assert pick_headline_property(props[:2])["name"] == "r_i_docking_score"
    assert pick_headline_property([]) is None


def test_list_maestro_entries_rejects_non_mae(tmp_path):
    pdb = tmp_path / "x.pdb"
    pdb.write_text("ATOM\nEND\n", encoding="utf-8")
    with pytest.raises(MaestroError, match="Not a Maestro"):
        list_maestro_entries(pdb)


def test_plain_mae_still_lists(tmp_path):
    pytest.importorskip("rdkit")
    from rdkit import Chem

    if not hasattr(Chem, "MaeWriter"):
        pytest.skip("RDKit MaeWriter not available")
    mae = tmp_path / "plain.mae"
    mol = Chem.AddHs(Chem.MolFromSmiles("C"))
    assert mol is not None
    mol.SetProp("_Name", "plain")
    writer = Chem.MaeWriter(str(mae))
    writer.write(mol)
    writer.close()
    info = list_maestro_entries(mae)
    assert len(info["entries"]) == 1
    assert isinstance(info["entries"][0]["properties"], list)
