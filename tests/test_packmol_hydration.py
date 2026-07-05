#!/usr/bin/env python3
"""
Packmol hydration test suite.

Unit tests for gatewizard.tools.packmol_hydration plus documentation example runners.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.tools.packmol_hydration import (
    build_hydrate_inp_text,
    check_packmol_available,
    detect_hydrogen_status,
    estimate_cavity_volume,
    exclusion_radius,
    hydrate_cavity,
    prepare_hydration_job,
    run_packmol,
    vdw_radius,
)
from gatewizard.core.structure_manager import Atom

PDB_6RV3 = Path(__file__).parent / "6RV3_AB.pdb"
BOX_MIN = (-13.03, -49.98, 18.67)
BOX_MAX = (6.97, -29.98, 38.67)

MINI_PDB = """\
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.389   0.000  1.00  0.00           O
END
"""


@pytest.fixture
def mini_pdb(tmp_path):
    path = tmp_path / "mini.pdb"
    path.write_text(MINI_PDB, encoding="utf-8")
    return str(path)


class TestPackmolHydrationUnit:
    def test_check_packmol_available_shape(self):
        info = check_packmol_available()
        assert "available" in info
        assert "version" in info
        assert "resolved_path" in info
        assert isinstance(info["available"], bool)

    def test_detect_hydrogen_status_none(self, mini_pdb):
        assert detect_hydrogen_status(pdb_file=mini_pdb) == "none"

    def test_detect_hydrogen_status_full(self):
        atoms = [
            Atom(1, "N", "N", (0, 0, 0), "ALA", 1, "A"),
            Atom(2, "HN", "H", (0, 0, 0), "ALA", 1, "A"),
            Atom(3, "CA", "C", (1, 0, 0), "ALA", 1, "A"),
            Atom(4, "HA", "H", (1, 0, 0), "ALA", 1, "A"),
        ]
        assert detect_hydrogen_status(atoms=atoms) == "full"

    def test_vdw_radius_heavy_atom_safe_inflation(self):
        explicit_c = vdw_radius("C", "explicit")
        safe_c = vdw_radius("C", "heavy_atom_safe")
        assert safe_c > explicit_c
        assert vdw_radius("O", "heavy_atom_safe") > vdw_radius("O", "explicit")

    def test_exclusion_radius_includes_solute(self):
        atom = Atom(1, "CA", "C", (0, 0, 0), "ALA", 1, "A")
        r = exclusion_radius(atom, "heavy_atom_safe", 2.5)
        assert r > 2.5

    def test_estimate_cavity_volume_mini(self, mini_pdb):
        result = estimate_cavity_volume(
            mini_pdb,
            box_min=(-5, -5, -5),
            box_max=(5, 5, 5),
        )
        assert result.box_volume_A3 == pytest.approx(1000.0)
        assert result.free_volume_A3 >= 0
        assert result.free_volume_A3 <= result.box_volume_A3
        assert result.hydrogen_status == "none"
        assert result.exclusion_mode == "heavy_atom_safe"
        assert result.suggested_waters >= 0

    @pytest.mark.skipif(not PDB_6RV3.is_file(), reason="6RV3_AB.pdb not found")
    def test_estimate_cavity_volume_6rv3(self):
        result = estimate_cavity_volume(
            str(PDB_6RV3),
            box_min=BOX_MIN,
            box_max=BOX_MAX,
        )
        assert result.free_volume_A3 > 0
        assert result.suggested_waters > 0

    def test_build_hydrate_inp_text(self, mini_pdb):
        text = build_hydrate_inp_text(
            protein_path=mini_pdb,
            tip3p_path="TIP3P.pdb",
            output_pdb="out.pdb",
            box_min=(0, 0, 0),
            box_max=(10, 10, 10),
            n_waters=5,
            solute_radius=2.5,
        )
        assert "tolerance 2.0" in text
        assert "radius 2.5" in text
        assert "number 5" in text
        assert "inside box" in text

    def test_prepare_hydration_job(self, mini_pdb, tmp_path):
        job = prepare_hydration_job(
            pdb_file=mini_pdb,
            job_dir=str(tmp_path / "job"),
            box_min=(0, 0, 0),
            box_max=(10, 10, 10),
            n_waters=3,
        )
        assert Path(job["packmol_inp_path"]).is_file()
        assert Path(job["tip3p_path"]).is_file()
        assert "inp_text" in job

    def test_run_packmol_mocked(self, tmp_path):
        inp = tmp_path / "packmol.inp"
        inp.write_text("tolerance 2.0\n", encoding="utf-8")
        out_pdb = tmp_path / "out.pdb"
        out_pdb.write_text("END\n", encoding="utf-8")

        def fake_run(cmd, cwd, capture_output, text, timeout, check):
            class R:
                returncode = 0
                stdout = "SUCCESS\n"
                stderr = ""

            return R()

        with patch(
            "gatewizard.tools.packmol_hydration.resolve_executable",
            return_value="/usr/bin/packmol",
        ), patch("subprocess.run", side_effect=fake_run):
            ok, log = run_packmol(str(inp), cwd=str(tmp_path))
        assert ok is True
        assert "SUCCESS" in log

    @pytest.mark.skipif(not PDB_6RV3.is_file(), reason="6RV3_AB.pdb not found")
    def test_hydrate_cavity_mocked(self, tmp_path):
        def fake_run(cmd, cwd, capture_output, text, timeout, check):
            out = Path(cwd) / "6RV3_AB_hydrated.pdb"
            out.write_text("END\n", encoding="utf-8")
            class R:
                returncode = 0
                stdout = "SUCCESS\n"
                stderr = ""
            return R()

        with patch(
            "gatewizard.tools.packmol_hydration.resolve_executable",
            return_value="/usr/bin/packmol",
        ), patch("subprocess.run", side_effect=fake_run):
            result = hydrate_cavity(
                pdb_file=str(PDB_6RV3),
                working_dir=str(tmp_path),
                output_folder_name="hydration_test",
                box_min=BOX_MIN,
                box_max=BOX_MAX,
                n_waters=5,
            )
        assert result.success is True
        assert Path(result.output_pdb).is_file()

    def test_invalid_box_raises(self, mini_pdb):
        with pytest.raises(ValueError):
            estimate_cavity_volume(mini_pdb, (5, 5, 5), (1, 1, 1))


class TestHydrationExamples:
    """Run hydration_example_*.py scripts (documentation examples)."""

    @pytest.fixture
    def temp_dir(self, tmp_path):
        orig = os.getcwd()
        os.chdir(tmp_path)
        yield tmp_path
        os.chdir(orig)

    @pytest.fixture(autouse=True)
    def cleanup_job_dirs(self):
        yield
        for pattern in ("hydration_*",):
            for target in Path(".").glob(pattern):
                if target.is_dir():
                    shutil.rmtree(target, ignore_errors=True)

    def test_run_example_scripts(self, temp_dir):
        examples_dir = Path(__file__).parent / "hydration_examples"
        if not examples_dir.exists():
            pytest.skip(f"Examples directory not found: {examples_dir}")

        example_files = sorted(examples_dir.glob("hydration_example_*.py"))
        if not example_files:
            pytest.skip("No hydration example files found")

        failed = []
        for example_file in example_files:
            example_num = example_file.stem.split("_")[-1]
            spec = importlib.util.spec_from_file_location(
                f"hydration_example_{example_num}", example_file
            )
            if spec is None or spec.loader is None:
                failed.append((example_num, "Could not load module"))
                continue
            module = importlib.util.module_from_spec(spec)
            try:
                spec.loader.exec_module(module)
            except Exception as exc:
                failed.append((example_num, f"{type(exc).__name__}: {exc}"))

        if failed:
            msg = "; ".join(f"{n}: {e}" for n, e in failed)
            pytest.fail(f"{len(failed)} example(s) failed: {msg}")

    @pytest.mark.parametrize("example_num", [f"{i:02d}" for i in range(1, 7)])
    def test_individual_examples(self, example_num, temp_dir):
        examples_dir = Path(__file__).parent / "hydration_examples"
        example_file = examples_dir / f"hydration_example_{example_num}.py"
        if not example_file.is_file():
            pytest.skip(f"Example file not found: {example_file}")

        spec = importlib.util.spec_from_file_location(
            f"hydration_example_{example_num}", example_file
        )
        assert spec is not None and spec.loader is not None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
