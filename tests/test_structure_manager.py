#!/usr/bin/env python3
"""
StructureManager Test Suite

This test suite covers:
1. Core structure manager API (StructureManager class)
2. Data model classes (Atom, Residue, ProteinStructure, Selection)
3. Documentation example workflows (Examples 1-10)

The test suite automatically discovers and runs all example scripts from:
    tests/viewer_examples/viewer_example_*.py

Usage:
    # Run all tests
    pytest tests/test_viewer.py -v

    # Run only example tests
    pytest tests/test_viewer.py::TestViewerExamples -v

    # Run specific example
    pytest tests/test_viewer.py::TestViewerExamples::test_individual_examples[01] -v
"""

import pytest
import sys
import os
import tempfile
from pathlib import Path
import importlib.util

sys.path.insert(0, str(Path(__file__).parent.parent))

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


class TestViewerExamples:
    """Test viewer examples from viewer_examples directory."""

    @pytest.fixture
    def temp_dir(self, tmp_path):
        """Provide a temporary directory for test outputs."""
        orig = os.getcwd()
        os.chdir(tmp_path)
        yield tmp_path
        os.chdir(orig)

    def test_run_example_scripts(self, temp_dir):
        """Test running actual example scripts from viewer_examples directory."""
        examples_dir = Path(__file__).parent / "viewer_examples"

        if not examples_dir.exists():
            pytest.skip(f"Examples directory not found: {examples_dir}")

        example_files = sorted(examples_dir.glob("viewer_example_*.py"))

        if not example_files:
            pytest.skip("No example files found in viewer_examples directory")

        print(f"\nFound {len(example_files)} example files to test")

        failed_examples = []
        passed_examples = []

        for example_file in example_files:
            example_num = example_file.stem.split("_")[-1]

            spec = importlib.util.spec_from_file_location(
                f"viewer_example_{example_num}", example_file
            )

            if spec is None or spec.loader is None:
                failed_examples.append((example_num, "Could not load module"))
                continue

            module = importlib.util.module_from_spec(spec)

            try:
                print(f"\nTesting Example {example_num}: {example_file.name}")
                spec.loader.exec_module(module)
                print(f"  Example {example_num} executed successfully")
                passed_examples.append(example_num)
            except Exception as e:
                error_msg = f"{type(e).__name__}: {str(e)}"
                print(f"  Example {example_num} failed: {error_msg}")
                failed_examples.append((example_num, error_msg))

        print(
            f"\nSummary: {len(passed_examples)} passed, "
            f"{len(failed_examples)} failed"
        )

        if failed_examples:
            for num, error in failed_examples:
                print(f"  Failed {num}: {error}")
            pytest.fail(f"{len(failed_examples)} example(s) failed")

    @pytest.mark.parametrize("example_num", [f"{i:02d}" for i in range(1, 11)])
    def test_individual_examples(self, example_num, temp_dir):
        """Test each example individually for better pytest reporting."""
        examples_dir = Path(__file__).parent / "viewer_examples"
        example_file = examples_dir / f"viewer_example_{example_num}.py"

        if not example_file.exists():
            pytest.skip(f"Example {example_num} not found")

        spec = importlib.util.spec_from_file_location(
            f"viewer_example_{example_num}", example_file
        )

        if spec is None or spec.loader is None:
            pytest.fail(f"Could not load example {example_num}")

        module = importlib.util.module_from_spec(spec)

        try:
            spec.loader.exec_module(module)
        except Exception as e:
            pytest.fail(f"Example {example_num} failed: {type(e).__name__}: {str(e)}")


# ============================================================================
# MAIN (for running outside pytest)
# ============================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("StructureManager Test Suite - Manual Run")
    print("=" * 60)

    examples_dir = Path(__file__).parent / "viewer_examples"

    if not examples_dir.exists():
        print(f"Examples directory not found: {examples_dir}")
        sys.exit(1)

    example_files = sorted(examples_dir.glob("viewer_example_*.py"))

    if not example_files:
        print("No example files found")
        sys.exit(1)

    print(f"Found {len(example_files)} examples to run\n")

    passed = []
    failed = []

    for example_file in example_files:
        example_num = example_file.stem.split("_")[-1]
        print(f"\n{'=' * 60}")
        print(f"Example {example_num}: {example_file.name}")
        print(f"{'=' * 60}")

        spec = importlib.util.spec_from_file_location(
            f"viewer_example_{example_num}", example_file
        )

        if spec is None or spec.loader is None:
            print(f"  Could not load module")
            failed.append(example_num)
            continue

        module = importlib.util.module_from_spec(spec)
        try:
            spec.loader.exec_module(module)
            print(f"  Passed")
            passed.append(example_num)
        except Exception as e:
            print(f"  Failed: {type(e).__name__}: {e}")
            failed.append(example_num)

    print(f"\n{'=' * 60}")
    print(
        f"Results: {len(passed)} passed, {len(failed)} failed "
        f"out of {len(example_files)}"
    )
    print(f"{'=' * 60}")

    if failed:
        print(f"Failed: {', '.join(failed)}")
        sys.exit(1)
