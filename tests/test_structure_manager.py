#!/usr/bin/env python3
"""
StructureManager Test Suite

This test suite covers:
1. Core structure manager API (StructureManager class)
2. Data model classes (Atom, Residue, ProteinStructure, Selection)
3. Documentation example workflows (structure_manager_example_*.py)

The test suite automatically discovers and runs all example scripts from:
    tests/structure_manager_examples/structure_manager_example_*.py

Usage:
    # Run all tests
    pytest tests/test_structure_manager.py -v

    # Run only example tests
    pytest tests/test_structure_manager.py::TestStructureManagerExamples -v

    # Run specific example
    pytest tests/test_structure_manager.py::TestStructureManagerExamples::test_example_script[structure_manager_example_01.py] -v
"""

import pytest
import sys
import os
import tempfile
from pathlib import Path

from tests.example_runner import parametrize_example_scripts, run_example_script

sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.utils.helpers import resolve_pdb_chain_id
from gatewizard.core.structure_manager import (
    StructureManager,
    ProteinStructure,
    Atom,
    Residue,
    Selection,
    parse_pdb,
    StructureError,
    AA_NAMES,
    BACKBONE_NAMES,
    assign_secondary_structure_map,
)

MINI_PDB = """\
HEADER    TEST PROTEIN
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00  0.00           C
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00  0.00           O
HETATM    9  O   HOH A 100      20.000  20.000  20.000  1.00  0.00           O
HETATM   10  C1  LIG B   1      30.000  30.000  30.000  1.00  0.00           C
HETATM   11  C2  LIG B   1      31.000  30.000  30.000  1.00  0.00           C
END
"""

CHARMM_SEGID_PDB = """\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00  0.00           N  PROT
ATOM      2  CA  ALA A   1       2.000   2.000   3.000  1.00  0.00           C  PROT
ATOM      3  C   ALA A   1       3.000   2.000   3.000  1.00  0.00           C  PROT
ATOM      4  O   ALA A   1       3.500   3.000   3.000  1.00  0.00           O  PROT
ATOM      5  N   GLY A   2       4.000   1.000   3.000  1.00  0.00           N  PROT
ATOM      6  CA  GLY A   2       5.000   1.000   3.000  1.00  0.00           C  PROT
ATOM      7  C   GLY A   2       6.000   1.000   3.000  1.00  0.00           C  PROT
ATOM      8  O   GLY A   2       6.500   2.000   3.000  1.00  0.00           O  PROT
HETATM    9  EPW WAT A 100      20.000  20.000  20.000  1.00  0.00          EPW  SOLV
END
"""


@pytest.fixture
def mini_pdb(tmp_path):
    p = tmp_path / "test.pdb"
    p.write_text(MINI_PDB)
    return str(p)


@pytest.fixture
def viewer(mini_pdb):
    v = StructureManager()
    v.load_structure(mini_pdb)
    return v


# ============================================================================
# SECTION 1: CORE VIEWER API TESTS
# ============================================================================


class TestStructureManager:
    """Test the StructureManager class core functionality."""

    def test_create(self):
        v = StructureManager()
        assert v.structure is None

    def test_load_structure(self, mini_pdb):
        v = StructureManager()
        info = v.load_structure(mini_pdb)
        assert info["n_atoms"] == 11
        assert info["n_residues"] > 0
        assert info["n_chains"] > 0
        assert info["n_bonds"] >= 0

    def test_load_nonexistent(self):
        v = StructureManager()
        with pytest.raises(StructureError):
            v.load_structure("/nonexistent/path.pdb")

    def test_get_chains(self, viewer):
        chains = viewer.get_chains()
        assert "A" in chains
        assert "B" in chains

    def test_get_residues(self, viewer):
        residues_a = viewer.get_residues(chain_id="A")
        assert len(residues_a) >= 2
        names = [r["name"] for r in residues_a]
        assert "ALA" in names
        assert "GLY" in names

    def test_get_residues_all(self, viewer):
        all_res = viewer.get_residues()
        assert len(all_res) >= 3  # ALA, GLY, HOH or LIG at least

    def test_get_secondary_structure_summary(self, viewer):
        ss = viewer.get_secondary_structure_summary()
        assert isinstance(ss, dict)

    def test_assign_secondary_structure_map_falls_back_without_psique(
        self, mini_pdb, monkeypatch
    ):
        import gatewizard.core.structure_manager as sm

        monkeypatch.setattr(sm, "_assign_ss_psique", lambda _path: None)
        ss_map = assign_secondary_structure_map(mini_pdb, method="auto")
        assert isinstance(ss_map, dict)
        assert ss_map
        assert all(code in {"H", "E", "C", "G", "I", "T"} for code in ss_map.values())

    def test_resolve_pdb_chain_id_charmm_style(self):
        assert resolve_pdb_chain_id("PROT", "A") == "A"
        assert resolve_pdb_chain_id("PROT", "") == "P"
        assert resolve_pdb_chain_id("", "B") == "B"
        assert resolve_pdb_chain_id("", "") == "A"

    def test_assign_secondary_structure_map_charmm_chain_keys(
        self, tmp_path, monkeypatch
    ):
        import gatewizard.core.structure_manager as sm

        p = tmp_path / "charmm.pdb"
        p.write_text(CHARMM_SEGID_PDB)
        monkeypatch.setattr(sm, "_assign_ss_psique", lambda _path: None)
        ss_map = assign_secondary_structure_map(str(p), method="auto")
        assert ("A", 1) in ss_map
        assert ("PROT", 1) not in ss_map

    def test_assign_secondary_structure_map_remaps_psique_chain_keys(
        self, tmp_path, monkeypatch
    ):
        import gatewizard.core.structure_manager as sm

        p = tmp_path / "charmm.pdb"
        p.write_text(CHARMM_SEGID_PDB)
        monkeypatch.setattr(
            sm,
            "_assign_ss_psique",
            lambda _path: {("PROT", 1): "H", ("PROT", 2): "E"},
        )
        ss_map = assign_secondary_structure_map(str(p), method="auto")
        assert ss_map.get(("A", 1)) == "H"
        assert ss_map.get(("A", 2)) == "E"

    def test_assign_secondary_structure_falls_back_when_psique_all_coil(
        self, tmp_path, monkeypatch
    ):
        import gatewizard.core.structure_manager as sm

        p = tmp_path / "charmm.pdb"
        p.write_text(CHARMM_SEGID_PDB)
        monkeypatch.setattr(
            sm,
            "_assign_ss_psique",
            lambda _path: {("Z", 999): "C"},
        )
        ss_map = assign_secondary_structure_map(str(p), method="auto")
        assert isinstance(ss_map, dict)
        assert ss_map
        assert ("Z", 999) not in ss_map
        assert ("A", 1) in ss_map
        assert ("A", 2) in ss_map
        assert all(code in {"H", "E", "C"} for code in ss_map.values())

    def test_assign_ss_psique_falls_back_to_protein_only(
        self, tmp_path, monkeypatch
    ):
        import gatewizard.core.structure_manager as sm

        p = tmp_path / "charmm.pdb"
        p.write_text(CHARMM_SEGID_PDB)
        calls: list[str] = []

        def fake_psique(path):
            calls.append(path)
            if path == str(p):
                raise RuntimeError("Could not guess element")
            return {("A", 1): "H", ("A", 2): "C"}

        monkeypatch.setattr(sm, "_run_psique_assign", fake_psique)
        ss_map = sm._assign_ss_psique(str(p))
        assert ss_map == {("A", 1): "H", ("A", 2): "C"}
        assert len(calls) == 2
        assert calls[0] == str(p)
        assert calls[1] != str(p)

    def test_select_by_criteria_all(self, viewer):
        idx = viewer.select_by_criteria("All")
        assert len(idx) == 11

    def test_select_by_criteria_protein(self, viewer):
        idx = viewer.select_by_criteria("Protein")
        for i in idx:
            assert viewer.structure.atoms[i].res_name in AA_NAMES

    def test_select_by_criteria_backbone(self, viewer):
        idx = viewer.select_by_criteria("Backbone")
        for i in idx:
            a = viewer.structure.atoms[i]
            assert a.res_name in AA_NAMES
            assert a.name in BACKBONE_NAMES

    def test_select_by_criteria_water(self, viewer):
        idx = viewer.select_by_criteria("Water")
        assert len(idx) >= 1
        for i in idx:
            assert viewer.structure.atoms[i].res_name in ("HOH", "WAT", "TIP")

    def test_select_by_criteria_ligand(self, viewer):
        idx = viewer.select_by_criteria("Ligand")
        assert len(idx) >= 1

    def test_select_by_criteria_chain(self, viewer):
        idx = viewer.select_by_criteria("Chain...", "A")
        for i in idx:
            assert viewer.structure.atoms[i].chain_id == "A"

    def test_select_by_criteria_range(self, viewer):
        idx = viewer.select_by_criteria("Residue range...", "A:1-2")
        for i in idx:
            a = viewer.structure.atoms[i]
            assert a.chain_id == "A"
            assert 1 <= a.res_id <= 2

    def test_auto_detect_molecules(self, viewer):
        sels = viewer.auto_detect_molecules()
        assert len(sels) >= 1
        names = [s.name for s in sels]
        assert "Protein" in names

    def test_rename_chain(self, viewer):
        count = viewer.rename_chain("A", "X")
        assert count > 0
        chains = viewer.get_chains()
        assert "X" in chains
        assert "A" not in chains

    def test_rename_residues(self, viewer):
        count = viewer.rename_residues("A", 1, 1, "MET")
        assert count > 0
        res = viewer.get_residues("A")
        met = [r for r in res if r["name"] == "MET"]
        assert len(met) > 0

    def test_renumber_residues(self, viewer):
        count = viewer.renumber_residues("A", 1, 2, new_start=100)
        assert count > 0
        res = viewer.get_residues("A")
        seq_ids = [r["seq_id"] for r in res]
        assert 100 in seq_ids

    def test_renumber_by_indices_does_not_touch_sibling_same_resid(self, tmp_path):
        """Two ligands can share UNK/900; renumber selection must not move the other."""
        pdb = tmp_path / "dup_unk.pdb"
        pdb.write_text(
            "\n".join(
                [
                    "HETATM    1  C1  UNK X 900       0.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    2  C2  UNK X 900       1.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    3  C1  UNK X 900      10.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    4  C2  UNK X 900      11.000   0.000   0.000  1.00  0.00           C",
                    "END",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        v = StructureManager()
        v.load_structure(str(pdb))
        assert all(a.res_id == 900 for a in v.structure.atoms)
        # First ligand only
        count = v.renumber_residues_by_indices([0, 1], new_start=1)
        assert count == 2
        assert v.structure.atoms[0].res_id == 1
        assert v.structure.atoms[1].res_id == 1
        assert v.structure.atoms[2].res_id == 900
        assert v.structure.atoms[3].res_id == 900

    def test_rename_by_indices_does_not_rename_sibling_same_resid(self, tmp_path):
        pdb = tmp_path / "dup_unk.pdb"
        pdb.write_text(
            "\n".join(
                [
                    "HETATM    1  C1  UNK X 900       0.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    2  C2  UNK X 900       1.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    3  C1  UNK X 900      10.000   0.000   0.000  1.00  0.00           C",
                    "HETATM    4  C2  UNK X 900      11.000   0.000   0.000  1.00  0.00           C",
                    "END",
                    "",
                ]
            ),
            encoding="utf-8",
        )
        v = StructureManager()
        v.load_structure(str(pdb))
        count = v.rename_residues_by_indices([0, 1], "LIG")
        assert count == 2
        assert v.structure.atoms[0].res_name == "LIG"
        assert v.structure.atoms[1].res_name == "LIG"
        assert v.structure.atoms[2].res_name == "UNK"
        assert v.structure.atoms[3].res_name == "UNK"

    def test_delete_atoms(self, viewer):
        water = viewer.select_by_criteria("Water")
        n_before = len(viewer.structure.atoms)
        removed = viewer.delete_atoms(water)
        assert removed == len(water)
        assert len(viewer.structure.atoms) == n_before - removed

    def test_save_pdb(self, viewer, tmp_path):
        out = str(tmp_path / "out.pdb")
        saved = viewer.save_pdb(out)
        assert os.path.isfile(saved)
        assert os.path.getsize(saved) > 0

    def test_save_pdb_no_structure(self, tmp_path):
        v = StructureManager()
        with pytest.raises(StructureError):
            v.save_pdb(str(tmp_path / "fail.pdb"))


# ============================================================================
# SECTION 2: DATA MODEL TESTS
# ============================================================================


class TestDataModel:
    """Test data model classes."""

    def test_atom_creation(self):
        import numpy as np

        a = Atom(1, "CA", "C", (1.0, 2.0, 3.0), "ALA", 1, "A")
        assert a.name == "CA"
        assert a.element == "C"
        assert np.allclose(a.coord, (1.0, 2.0, 3.0))

    def test_residue_creation(self):
        r = Residue("ALA", 1, "A")
        assert r.name == "ALA"
        assert r.seq_id == 1
        a = Atom(1, "CA", "C", (1.0, 2.0, 3.0), "ALA", 1, "A")
        r.add_atom(a)
        assert len(r.atoms) == 1

    def test_protein_structure(self):
        s = ProteinStructure()
        assert len(s.atoms) == 0
        assert len(s.residues) == 0
        assert len(s.bonds) == 0

    def test_selection_creation(self):
        s = Selection("Test", [0, 1, 2])
        assert s.name == "Test"
        assert len(s.atom_indices) == 3
        assert s.visible is True
        assert s.representation == "ball_stick"

    def test_parse_pdb(self, mini_pdb):
        struct = parse_pdb(mini_pdb)
        assert isinstance(struct, ProteinStructure)
        assert len(struct.atoms) == 11
        assert len(struct.residues) > 0

    def test_write_pdb_roundtrip(self, mini_pdb, tmp_path):
        struct = parse_pdb(mini_pdb)
        out = str(tmp_path / "roundtrip.pdb")
        struct.write_pdb(out)
        struct2 = parse_pdb(out)
        assert len(struct2.atoms) == len(struct.atoms)

    def test_build_bonds(self, mini_pdb):
        struct = parse_pdb(mini_pdb)
        struct.build_bonds()
        assert len(struct.bonds) >= 0


# ============================================================================
# SECTION 3: EXAMPLE TESTS
# ============================================================================


class TestStructureManagerExamples:
    """Run each structure_manager_example_*.py script once."""

    @parametrize_example_scripts(
        Path(__file__).parent / "structure_manager_examples",
        "structure_manager_example_*.py",
    )
    def test_example_script(self, script):
        run_example_script(script)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
