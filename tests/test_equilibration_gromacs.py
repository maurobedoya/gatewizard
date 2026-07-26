"""Tests for GROMACSEquilibrationManager.

These tests cover:
- Class constants and mappings
- Default stage parameters
- MDP file generation (force constant conversion, temperature, nsteps)
- Template existence (all 4 ensembles × 8 files)
- File discovery helpers
- COM colvars helpers (_build_com_colvars_config / _build_com_colvars_activation_block)
- Integration tests that require gmx executable (skipped when not available)
- AMBER → GROMACS conversion (skipped when ParmEd not available)
"""

import json
import re
import shutil
import subprocess
from pathlib import Path
from typing import Optional
from unittest.mock import MagicMock, patch

import pytest

from gatewizard.tools.equilibration import (
    EquilibrationStage,
    GROMACSEquilibrationManager,
    _build_com_colvars_activation_block,
    _build_com_colvars_config,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

POPC_DIR = Path(__file__).parent / "equilibration_examples" / "popc_membrane"
PRMTOP = POPC_DIR / "system.prmtop"
INPCRD = POPC_DIR / "system.inpcrd"
BILAYER_PDB = POPC_DIR / "bilayer_protein_protonated_prepared_lipid.pdb"

TEMPLATES_DIR = Path(__file__).parent.parent / "equilibration" / "gromacs"

GMX_AVAILABLE: bool = shutil.which("gmx") is not None
PARMED_AVAILABLE: bool = True
try:
    import parmed  # noqa: F401
except ImportError:
    PARMED_AVAILABLE = False

MDA_AVAILABLE: bool = True
try:
    import MDAnalysis  # noqa: F401
except ImportError:
    MDA_AVAILABLE = False


def _make_manager(tmp_path: Path) -> GROMACSEquilibrationManager:
    return GROMACSEquilibrationManager(tmp_path)


# ---------------------------------------------------------------------------
# Class constants
# ---------------------------------------------------------------------------


class TestClassConstants:
    def test_scheme_mapping_keys(self):
        assert set(GROMACSEquilibrationManager.SCHEME_MAPPING.keys()) == {
            "NVT",
            "NPT",
            "NPAT",
            "NPgT",
        }

    def test_scheme_mapping_values(self):
        assert GROMACSEquilibrationManager.SCHEME_MAPPING["NPT"] == "02_NPT"
        assert GROMACSEquilibrationManager.SCHEME_MAPPING["NPAT"] == "03_NPAT"

    def test_template_mapping_keys(self):
        expected = {
            "step0_minimization",
            "step1",
            "step2",
            "step3",
            "step4",
            "step5",
            "step6",
            "step7_production",
        }
        assert set(GROMACSEquilibrationManager.TEMPLATE_MAPPING.keys()) == expected

    def test_stage_index_to_key(self):
        m = GROMACSEquilibrationManager.STAGE_INDEX_TO_KEY
        assert m[0] == "step0_minimization"
        assert m[7] == "step7_production"
        for i in range(1, 7):
            assert m[i] == f"step{i}"

    def test_kcal_to_kj_constant(self):
        assert GROMACSEquilibrationManager._KCAL_TO_KJ == pytest.approx(418.4)


# ---------------------------------------------------------------------------
# Template existence
# ---------------------------------------------------------------------------


class TestTemplateExistence:
    ENSEMBLES = ["01_NVT", "02_NPT", "03_NPAT", "04_NPgT"]
    FILENAMES = [
        "step6.0_minimization.mdp",
        "step6.1_equilibration.mdp",
        "step6.2_equilibration.mdp",
        "step6.3_equilibration.mdp",
        "step6.4_equilibration.mdp",
        "step6.5_equilibration.mdp",
        "step6.6_equilibration.mdp",
        "step7_production.mdp",
    ]

    @pytest.mark.parametrize("ensemble", ENSEMBLES)
    @pytest.mark.parametrize("filename", FILENAMES)
    def test_template_exists(self, ensemble, filename):
        p = TEMPLATES_DIR / ensemble / filename
        assert p.exists(), f"Missing template: {p}"


# ---------------------------------------------------------------------------
# Default stage parameters
# ---------------------------------------------------------------------------


class TestGetDefaultStageParams:
    def test_returns_list(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params()
        assert isinstance(stages, list)

    def test_default_length(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params()
        # 1 minimization + 6 equilibration
        assert len(stages) == 7

    def test_with_production(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params(
            include_production=True
        )
        assert len(stages) == 8
        assert stages[-1].name == "Production"

    def test_first_stage_is_minimization(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params()
        assert stages[0].minimize_steps > 0

    def test_force_constants_decrease(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params()
        # backbone FC should be non-increasing across stages
        fc_bb = [s.constraints.get("protein_backbone", 0.0) for s in stages]
        for i in range(1, len(fc_bb) - 1):
            assert fc_bb[i] >= fc_bb[i + 1], (
                f"Backbone FC increased from stage {i} to {i+1}: "
                f"{fc_bb[i]} → {fc_bb[i+1]}"
            )

    def test_ensemble_propagated(self):
        for scheme in ("NVT", "NPT", "NPAT", "NPgT"):
            stages = GROMACSEquilibrationManager.get_default_stage_params(scheme)
            for s in stages:
                assert s.ensemble == scheme

    def test_invalid_scheme_raises(self):
        with pytest.raises(ValueError):
            GROMACSEquilibrationManager.get_default_stage_params("INVALID")

    def test_equilibration_4_uses_1fs_timestep(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params()
        equil_4 = next(s for s in stages if s.name == "Equilibration 4")
        assert equil_4.timestep == 1.0
        equil_5 = next(s for s in stages if s.name == "Equilibration 5")
        assert equil_5.timestep == 2.0

    def test_temperature_propagated(self):
        stages = GROMACSEquilibrationManager.get_default_stage_params(temperature=300.0)
        for s in stages:
            assert s.temperature == 300.0


# ---------------------------------------------------------------------------
# MDP file generation
# ---------------------------------------------------------------------------


class TestGenerateMdpFile:
    """Tests for MDP file content generation via GROMACSEquilibrationManager."""

    @pytest.fixture
    def manager(self, tmp_path):
        return _make_manager(tmp_path)

    def _stage_dict(self, **kwargs):
        base = {
            "name": "Test",
            "ensemble": "NPT",
            "time_ns": 0.25,
            "timestep": 2.0,
            "temperature": 300.0,
            "minimize_steps": 0,
            "constraints": {
                "protein_backbone": 5.0,
                "protein_sidechain": 2.5,
                "lipid_head": 1.0,
                "lipid_tail": 0.0,
                "water": 0.0,
                "ions": 0.0,
                "other": 0.0,
            },
        }
        base.update(kwargs)
        return base

    def test_minimization_stage_generated(self, manager):
        params = self._stage_dict(minimize_steps=5000)
        content = manager.generate_mdp_file(
            stage_name="Minimization",
            stage_params=params,
            stage_index=0,
            scheme_type="NPT",
        )
        assert (
            "minimization" in content.lower()
            or "steep" in content.lower()
            or "nsteps" in content.lower()
        )

    def test_nsteps_substitution(self, manager):
        """2 ns at 2 fs timestep → 1 000 000 steps."""
        params = self._stage_dict(time_ns=2.0, timestep=2.0, minimize_steps=0)
        content = manager.generate_mdp_file(
            stage_name="Stage 1",
            stage_params=params,
            stage_index=1,
            scheme_type="NPT",
        )
        assert "1000000" in content

    def test_temperature_substitution(self, manager):
        params = self._stage_dict(temperature=280.0)
        content = manager.generate_mdp_file(
            stage_name="Stage 1",
            stage_params=params,
            stage_index=1,
            scheme_type="NPT",
        )
        assert "280.00" in content

    def test_force_constant_conversion(self, manager):
        """5.0 kcal/mol/Å² → 2092.0 kJ/mol/nm²."""
        params = self._stage_dict(
            constraints={
                "protein_backbone": 5.0,
                "protein_sidechain": 2.5,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "water": 0.0,
                "ions": 0.0,
                "other": 0.0,
            }
        )
        content = manager.generate_mdp_file(
            stage_name="Stage 1",
            stage_params=params,
            stage_index=1,
            scheme_type="NPT",
        )
        expected_bb = f"{5.0 * 418.4:.1f}"  # "2092.0"
        assert expected_bb in content, f"Expected FC {expected_bb} in MDP define line"

    def test_zero_constraints_removes_define(self, manager):
        """When all force constants are 0, the define= line should be absent or empty."""
        params = self._stage_dict(
            constraints={
                k: 0.0
                for k in [
                    "protein_backbone",
                    "protein_sidechain",
                    "lipid_head",
                    "lipid_tail",
                    "water",
                    "ions",
                    "other",
                ]
            }
        )
        content = manager.generate_mdp_file(
            stage_name="Stage 6",
            stage_params=params,
            stage_index=6,
            scheme_type="NPT",
        )
        # Either the define line is gone or the FC values are 0.0
        if "define" in content:
            assert "POSRES_FC_BB=0.0" in content

    def test_production_has_no_define(self, manager):
        params = self._stage_dict(
            time_ns=50.0,
            constraints={
                k: 0.0
                for k in [
                    "protein_backbone",
                    "protein_sidechain",
                    "lipid_head",
                    "lipid_tail",
                    "water",
                    "ions",
                    "other",
                ]
            },
        )
        content = manager.generate_mdp_file(
            stage_name="Production",
            stage_params=params,
            stage_index=7,
            scheme_type="NPT",
        )
        # Production template should have no active POSRES define
        # (the define line should have been stripped)
        assert "POSRES_FC_BB" not in content

    def test_all_ensembles_generate(self, manager):
        for scheme in ("NVT", "NPT", "NPAT", "NPgT"):
            content = manager.generate_mdp_file(
                stage_name="Test",
                stage_params=self._stage_dict(ensemble=scheme),
                stage_index=1,
                scheme_type=scheme,
            )
            assert len(content) > 100

    def test_get_mdp_filename(self, manager):
        assert manager._get_mdp_filename(0) == "step0_minimization.mdp"
        assert manager._get_mdp_filename(1) == "step1_equilibration.mdp"
        assert manager._get_mdp_filename(6) == "step6_equilibration.mdp"
        assert manager._get_mdp_filename(7) == "step7_production.mdp"


# ---------------------------------------------------------------------------
# Run script
# ---------------------------------------------------------------------------


class TestGenerateRunScript:
    def test_run_script_created(self, tmp_path):
        manager = _make_manager(tmp_path)
        script = manager.generate_run_script(
            gromacs_dir=tmp_path,
            gro_name="system.gro",
            top_name="topol_posres.top",
            ndx_name="index.ndx",
            n_stages=6,
        )
        assert script.exists()
        text = script.read_text()
        assert "step0_minimization" in text
        assert "step7_production" in text
        assert "source /usr/local/gromacs/bin/GMXRC" not in text

    def test_run_script_optional_gmxrc(self, tmp_path):
        manager = _make_manager(tmp_path)
        gmxrc = tmp_path / "GMXRC"
        gmxrc.write_text("# test\n")
        script = manager.generate_run_script(
            gromacs_dir=tmp_path,
            gro_name="system.gro",
            top_name="topol_posres.top",
            ndx_name="index.ndx",
            n_stages=6,
            gmxrc_path=str(gmxrc),
        )
        text = script.read_text()
        assert f'source "{gmxrc}"' in text

    def test_run_script_executable(self, tmp_path):
        manager = _make_manager(tmp_path)
        script = manager.generate_run_script(
            gromacs_dir=tmp_path,
            gro_name="system.gro",
            top_name="topol.top",
            ndx_name=None,
            n_stages=6,
        )
        assert script.stat().st_mode & 0o100  # executable bit

    def test_run_script_n_stages(self, tmp_path):
        manager = _make_manager(tmp_path)
        for n in (3, 6):
            script = manager.generate_run_script(
                gromacs_dir=tmp_path,
                gro_name="system.gro",
                top_name="topol.top",
                ndx_name=None,
                n_stages=n,
            )
            text = script.read_text()
            assert f"step{n}_equilibration" in text

    def test_run_script_includes_cpu_gpu_resources(self, tmp_path):
        manager = _make_manager(tmp_path)
        script = manager.generate_run_script(
            gromacs_dir=tmp_path,
            gro_name="system.gro",
            top_name="topol.top",
            ndx_name=None,
            n_stages=2,
            cpu_cores=4,
            use_gpu=True,
            gpu_id=0,
            num_gpus=1,
        )
        text = script.read_text()
        # Dynamics: GROMACS requires -ntmpi when combining GPUs with -ntomp.
        assert "-ntmpi 1" in text
        assert "-ntomp 4" in text
        assert "-nb gpu" in text
        assert "-pme gpu" in text
        assert "-gpu_id 0" in text
        # Energy minimisation must stay on CPU (PME GPU rejects steep/cg).
        em_lines = [
            line
            for line in text.splitlines()
            if "mdrun" in line and "step0_minimization" in line
        ]
        assert em_lines
        assert "-nb cpu" in em_lines[0]
        assert "-pme cpu" in em_lines[0]
        assert "-nb gpu" not in em_lines[0]
        assert "-pme gpu" not in em_lines[0]
        assert "-ntmpi 1" in em_lines[0]
        assert "-ntomp 4" in em_lines[0]
        # GPU flags only on dynamics stages (2 eq + production).
        assert text.count("-nb gpu") >= 3

    def test_run_script_multi_gpu_id_string(self, tmp_path):
        manager = _make_manager(tmp_path)
        script = manager.generate_run_script(
            gromacs_dir=tmp_path,
            gro_name="system.gro",
            top_name="topol.top",
            ndx_name=None,
            n_stages=1,
            cpu_cores=8,
            use_gpu=True,
            gpu_id=1,
            num_gpus=2,
        )
        text = script.read_text()
        assert "-gpu_id 12" in text
        assert "-ntmpi 2" in text
        assert "-ntomp 4" in text  # 8 processors / 2 ranks

    def test_run_script_cpu_only_omits_gpu_flags(self, tmp_path):
        manager = _make_manager(tmp_path)
        script = manager.generate_run_script(
            gromacs_dir=tmp_path,
            gro_name="system.gro",
            top_name="topol.top",
            ndx_name=None,
            n_stages=1,
            cpu_cores=4,
            use_gpu=False,
        )
        text = script.read_text()
        assert "-ntmpi 1" in text
        assert "-ntomp 4" in text
        assert "-nb cpu" in text
        assert "-pme cpu" in text
        assert "-nb gpu" not in text
        assert "-gpu_id" not in text


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------


class TestFindSystemFiles:
    def test_detects_gromacs_files(self, tmp_path):
        (tmp_path / "system.gro").write_text("GROMACS GRO\n")
        (tmp_path / "topol.top").write_text("[system]\n")
        manager = _make_manager(tmp_path)
        result = manager.find_system_files()
        assert result is not None
        assert "gro" in result
        assert "top" in result

    def test_detects_amber_files(self, tmp_path):
        (tmp_path / "system.prmtop").write_text("")
        (tmp_path / "system.inpcrd").write_text("")
        manager = _make_manager(tmp_path)
        result = manager.find_system_files()
        assert result is not None
        assert "prmtop" in result
        assert "inpcrd" in result

    def test_returns_none_when_empty(self, tmp_path):
        manager = _make_manager(tmp_path)
        result = manager.find_system_files()
        assert result is None


# ---------------------------------------------------------------------------
# COM colvars helpers (module-level functions)
# ---------------------------------------------------------------------------


class TestBuildComColvarsConfig:
    @pytest.fixture
    def mock_ag(self):
        """Mock MDAnalysis AtomGroup with 3 fake atoms."""
        ag = MagicMock()
        atoms_list = []
        for i, pos in enumerate([(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), (7.0, 8.0, 9.0)]):
            a = MagicMock()
            a.index = i
            a.position = pos
            atoms_list.append(a)
        ag.atoms = atoms_list
        ag.__len__ = MagicMock(return_value=3)
        return ag

    def test_returns_string(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=1.0,
            y0=2.0,
            z0=3.0,
            com_k=10.0,
            add_rotation=False,
            rot_k=2000.0,
            ag=mock_ag,
        )
        assert isinstance(result, str)
        assert len(result) > 50

    def test_contains_center_distance_cv(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=1.0,
            y0=2.0,
            z0=3.0,
            com_k=10.0,
            add_rotation=False,
            rot_k=2000.0,
            ag=mock_ag,
        )
        assert "name center_x" in result
        assert "distanceZ {" in result
        assert "dummyAtom" in result

    def test_contains_harmonic_block(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=1.0,
            y0=2.0,
            z0=3.0,
            com_k=5.0,
            add_rotation=False,
            rot_k=2000.0,
            ag=mock_ag,
        )
        assert "harmonic" in result
        assert "colvars center_x" in result
        assert "forceConstant 5.0000" in result

    def test_atom_numbers_present(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="42 99 101",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=False,
            rot_k=0.0,
            ag=mock_ag,
        )
        assert "42 99 101" in result

    def test_rotation_cv_added_when_requested(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=True,
            rot_k=500.0,
            ag=mock_ag,
            ref_positions_file="system.pdb",
        )
        assert "name rotation" in result
        assert "orientation" in result
        assert "refPositionsFile system.pdb" in result
        assert "centers (1.0, 0.0, 0.0, 0.0)" in result

    def test_namd_rotation_uses_inline_refpositions_when_requested(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=True,
            rot_k=500.0,
            ag=mock_ag,
            engine="namd",
            rotation_ref_positions_mode="refPositions",
        )
        assert "refPositions {" in result
        assert "refPositionsFile" not in result

    def test_file_mode_requires_refpositionsfile(self, mock_ag):
        with pytest.raises(ValueError, match="requires ref_positions_file"):
            _build_com_colvars_config(
                atom_numbers="1 2 3",
                x0=0.0,
                y0=0.0,
                z0=0.0,
                com_k=1.0,
                add_rotation=True,
                rot_k=500.0,
                ag=mock_ag,
                engine="namd",
                rotation_ref_positions_mode="refPositionsFile",
                ref_positions_file=None,
            )

    def test_refpositionscol_and_value_emitted_in_file_mode(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=True,
            rot_k=500.0,
            ag=mock_ag,
            engine="namd",
            rotation_ref_positions_mode="refPositionsFile",
            ref_positions_file="system.pdb",
            ref_positions_col="B",
            ref_positions_col_value=2.0,
        )
        assert "refPositionsFile system.pdb" in result
        assert "refPositionsCol B" in result
        assert "refPositionsColValue 2" in result

    def test_refpositionscolvalue_requires_refpositionscol(self, mock_ag):
        with pytest.raises(ValueError, match="requires ref_positions_col"):
            _build_com_colvars_config(
                atom_numbers="1 2 3",
                x0=0.0,
                y0=0.0,
                z0=0.0,
                com_k=1.0,
                add_rotation=True,
                rot_k=500.0,
                ag=mock_ag,
                rotation_ref_positions_mode="refPositionsFile",
                ref_positions_file="system.pdb",
                ref_positions_col_value=1.0,
            )

    def test_gromacs_can_use_refpositionsfile_mode(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=True,
            rot_k=500.0,
            ag=mock_ag,
            engine="gromacs",
            rotation_ref_positions_mode="refPositionsFile",
            ref_positions_file="ref.pdb",
        )
        assert "refPositionsFile ref.pdb" in result

    def test_rotation_cv_absent_by_default(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=False,
            rot_k=0.0,
            ag=mock_ag,
        )
        assert "orientation" not in result

    def test_gromacs_engine_uses_hash_comments(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=False,
            rot_k=0.0,
            ag=mock_ag,
            engine="gromacs",
        )
        assert "# Colvars" in result
        assert "; Colvars" not in result

    def test_gromacs_engine_uses_direct_atomnumbers_blocks(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2 3",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=True,
            rot_k=500.0,
            ag=mock_ag,
            engine="gromacs",
            ref_positions_file="system.pdb",
        )
        assert "group1 {\n         atomNumbers {" not in result
        assert "main {\n         atomNumbers {" in result
        assert "orientation {\n      atoms {" in result
        assert "group1 {\n         atoms {" not in result

    def test_namd_engine_uses_hash(self, mock_ag):
        result = _build_com_colvars_config(
            atom_numbers="1 2",
            x0=0.0,
            y0=0.0,
            z0=0.0,
            com_k=1.0,
            add_rotation=False,
            rot_k=0.0,
            ag=mock_ag,
            engine="namd",
        )
        assert "# Colvars" in result


class TestBuildComColvarsActivationBlock:
    def test_namd_activation_block(self):
        result = _build_com_colvars_activation_block("namd", "com_restraint.col")
        assert "colvars on" in result
        assert "colvarsConfig com_restraint.col" in result

    def test_gromacs_activation_block(self):
        result = _build_com_colvars_activation_block("gromacs", "com_restraint.dat")
        assert "colvars-active         = yes" in result
        assert "colvars-configfile     = com_restraint.dat" in result


class TestComColvarsSetupPaths:
    def test_setup_gromacs_com_colvars_uses_restraints_path(
        self, tmp_path, monkeypatch
    ):
        manager = _make_manager(tmp_path)

        prmtop = tmp_path / "system.prmtop"
        inpcrd = tmp_path / "system.inpcrd"
        pdb = tmp_path / "system.pdb"
        prmtop.write_text("dummy")
        inpcrd.write_text("dummy")
        pdb.write_text(
            "ATOM      1  CA  ALA A   1      11.000  21.000  31.000  1.00 20.00           C\nEND\n"
        )

        def _fake_convert_from_amber(prmtop, inpcrd, output_dir, bilayer_pdb=None):
            gro = output_dir / "system.gro"
            top = output_dir / "topol.top"
            gro.write_text(
                "test\n1\n    1ALA     CA    1   0.000   0.000   0.000\n   1.0   1.0   1.0\n"
            )
            top.write_text("[ moleculetype ]\nProtein 3\n")
            return {"gro": gro, "top": top}

        monkeypatch.setattr(manager, "convert_from_amber", _fake_convert_from_amber)

        monkeypatch.setattr(
            manager,
            "generate_index_ndx",
            lambda gro_path, index_path: index_path.write_text("[ System ]\n1\n"),
        )
        monkeypatch.setattr(
            manager,
            "generate_mdp_file",
            lambda **kwargs: "integrator = md\n",
        )

        def _fake_run_script(
            gromacs_dir, gro_name, top_name, ndx_name, n_stages, gmx_executable
        ):
            run_script = gromacs_dir / "run_equilibration.sh"
            run_script.write_text("#!/bin/bash\n")
            return run_script

        monkeypatch.setattr(manager, "generate_run_script", _fake_run_script)

        captured = {}

        def _fake_generate_com_colvars_config(**kwargs):
            output_file = kwargs["output_file"]
            captured["output_file"] = output_file
            output_file.parent.mkdir(parents=True, exist_ok=True)
            output_file.write_text("; mock colvars\n")
            return output_file

        monkeypatch.setattr(
            manager,
            "generate_com_colvars_config",
            _fake_generate_com_colvars_config,
        )

        out_dir = tmp_path / "gmx_out"
        result = manager.setup_gromacs_equilibration(
            system_files={
                "prmtop": str(prmtop),
                "inpcrd": str(inpcrd),
                "pdb": str(pdb),
            },
            stage_params_list=[
                {
                    "name": "Equilibration 1",
                    "ensemble": "NPT",
                    "time_ns": 0.01,
                    "timestep": 1.0,
                    "temperature": 310.15,
                    "constraints": {},
                }
            ],
            output_name=str(out_dir),
            add_com_restraint=True,
        )

        expected_colvars = out_dir / "restraints" / "com_restraint.dat"
        assert captured["output_file"] == expected_colvars
        assert result["com_colvars"] == expected_colvars

        mdp_text = result["mdp_files"][0].read_text()
        assert "colvars-configfile     = restraints/com_restraint.dat" in mdp_text


# ---------------------------------------------------------------------------
# Integration: AMBER conversion (requires ParmEd + test data)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not PARMED_AVAILABLE or not PRMTOP.exists(),
    reason="ParmEd not installed or POPC test data not found",
)
class TestAmberConversion:
    def test_produces_gro_and_top(self, tmp_path):
        manager = _make_manager(tmp_path)
        result = manager.convert_from_amber(
            prmtop=PRMTOP,
            inpcrd=INPCRD,
            output_dir=tmp_path,
            bilayer_pdb=BILAYER_PDB if BILAYER_PDB.exists() else None,
        )
        assert result["gro"].exists(), "system.gro not created"
        assert result["top"].exists(), "topol.top not created"

    def test_gro_has_atoms(self, tmp_path):
        manager = _make_manager(tmp_path)
        result = manager.convert_from_amber(
            prmtop=PRMTOP,
            inpcrd=INPCRD,
            output_dir=tmp_path,
            bilayer_pdb=BILAYER_PDB if BILAYER_PDB.exists() else None,
        )
        # GRO file: line 2 is the atom count
        lines = result["gro"].read_text().splitlines()
        n_atoms = int(lines[1].strip())
        assert n_atoms > 0

    def test_topol_contains_moleculetype(self, tmp_path):
        manager = _make_manager(tmp_path)
        result = manager.convert_from_amber(
            prmtop=PRMTOP,
            inpcrd=INPCRD,
            output_dir=tmp_path,
            bilayer_pdb=BILAYER_PDB if BILAYER_PDB.exists() else None,
        )
        top_text = result["top"].read_text()
        assert "[ moleculetype ]" in top_text


# ---------------------------------------------------------------------------
# Integration: full setup (requires gmx)
# ---------------------------------------------------------------------------


@pytest.mark.skipif(
    not GMX_AVAILABLE or not PARMED_AVAILABLE or not PRMTOP.exists(),
    reason="gmx not available, ParmEd missing, or test data absent",
)
class TestSetupGromacs:
    def test_setup_creates_output_dir(self, tmp_path):
        shutil.copy2(PRMTOP, tmp_path / PRMTOP.name)
        shutil.copy2(INPCRD, tmp_path / INPCRD.name)
        if BILAYER_PDB.exists():
            shutil.copy2(BILAYER_PDB, tmp_path / BILAYER_PDB.name)

        manager = _make_manager(tmp_path)
        stages = GROMACSEquilibrationManager.get_default_stage_params("NPT")
        result = manager.setup_gromacs_equilibration(
            stage_params_list=stages,
        )
        gdir = result["gromacs_dir"]
        assert gdir.exists()

    def test_setup_creates_mdp_files(self, tmp_path):
        shutil.copy2(PRMTOP, tmp_path / PRMTOP.name)
        shutil.copy2(INPCRD, tmp_path / INPCRD.name)
        if BILAYER_PDB.exists():
            shutil.copy2(BILAYER_PDB, tmp_path / BILAYER_PDB.name)

        manager = _make_manager(tmp_path)
        stages = GROMACSEquilibrationManager.get_default_stage_params("NPT")
        result = manager.setup_gromacs_equilibration(stage_params_list=stages)
        mdp_files = result["mdp_files"]
        assert len(mdp_files) >= 7
        for f in mdp_files:
            assert f.exists()

    def test_run_script_sources_gmxrc(self, tmp_path):
        shutil.copy2(PRMTOP, tmp_path / PRMTOP.name)
        shutil.copy2(INPCRD, tmp_path / INPCRD.name)
        if BILAYER_PDB.exists():
            shutil.copy2(BILAYER_PDB, tmp_path / BILAYER_PDB.name)

        manager = _make_manager(tmp_path)
        stages = GROMACSEquilibrationManager.get_default_stage_params("NPT")
        result = manager.setup_gromacs_equilibration(stage_params_list=stages)
        script_text = result["run_script"].read_text()
        assert "source /usr/local/gromacs/bin/GMXRC" in script_text

    def test_setup_invalid_scheme_raises(self, tmp_path):
        manager = _make_manager(tmp_path)
        with pytest.raises(ValueError):
            manager.setup_gromacs_equilibration(
                system_files={
                    "gro": str(tmp_path / "x.gro"),
                    "top": str(tmp_path / "x.top"),
                },
                scheme_type="BADSCHEME",
                stage_params_list=[
                    {
                        "name": "s",
                        "ensemble": "BADSCHEME",
                        "time_ns": 0.1,
                        "timestep": 2.0,
                        "temperature": 310.15,
                        "constraints": {},
                    }
                ],
            )


# ---------------------------------------------------------------------------
# Water / ions / other / custom MDA posres
# ---------------------------------------------------------------------------


def _minimal_top_with_ions(tmp_path: Path) -> Path:
    """Write a tiny topology: Protein(2) + POPC(3) + NA(1) + CL(1) + SOL(3)."""
    top = tmp_path / "topol.top"
    top.write_text(
        """\
[ defaults ]
; nbfunc comb-rule gen-pairs fudgeLJ fudgeQQ
1 2 yes 0.5 0.8333

[ moleculetype ]
; name  nrexcl
Protein  3

[ atoms ]
; nr type resnr residue atom cgnr charge mass
1  CX  1  ALA  CA  1  0.0  12.0
2  HX  1  ALA  H   2  0.0   1.0

[ moleculetype ]
POPC  3

[ atoms ]
1  LP  1  POPC  C1  1  0.0  12.0
2  LP  1  POPC  C2  2  0.0  12.0
3  LP  1  POPC  C3  3  0.0  12.0

[ moleculetype ]
NA  1

[ atoms ]
1  Na  1  NA  NA  1  1.0  23.0

[ moleculetype ]
CL  1

[ atoms ]
1  Cl  1  CL  CL  1  -1.0  35.0

[ moleculetype ]
SOL  2

[ atoms ]
1  OW  1  SOL  OW  1  0.0  16.0
2  HW  1  SOL  HW1 2  0.0   1.0
3  HW  1  SOL  HW2 3  0.0   1.0

[ system ]
Test

[ molecules ]
Protein  1
POPC     2
NA       2
CL       2
SOL      3
"""
    )
    return top


def _minimal_structure_matching_top(tmp_path: Path) -> Path:
    """GRO with 21 atoms matching _minimal_top_with_ions (MDA-readable)."""
    gro = tmp_path / "system.gro"
    n = 21
    # residue names / atom names aligned with topology copies
    specs = []
    specs += [(1, "ALA", "CA"), (1, "ALA", "H")]
    for r in (2, 3):
        specs += [(r, "POPC", "C1"), (r, "POPC", "C2"), (r, "POPC", "C3")]
    specs += [(4, "NA", "NA"), (5, "NA", "NA")]
    specs += [(6, "CL", "CL"), (7, "CL", "CL")]
    for r in (8, 9, 10):
        specs += [(r, "SOL", "OW"), (r, "SOL", "HW1"), (r, "SOL", "HW2")]
    assert len(specs) == n
    lines = ["GateWizard test system", f"{n:5d}"]
    for i, (resid, resname, atom) in enumerate(specs, start=1):
        lines.append(
            f"{resid:5d}{resname:<5s}{atom:>5s}{i:5d}"
            f"{0.1 * i:8.3f}{0.200:8.3f}{0.300:8.3f}"
        )
    lines.append("   4.00000   4.00000   4.00000")
    gro.write_text("\n".join(lines) + "\n")
    return gro


class TestGromacsPosresHelpers:
    def test_classify_moltype(self):
        assert GROMACSEquilibrationManager._classify_moltype("SOL") == "water"
        assert GROMACSEquilibrationManager._classify_moltype("NA") == "ions"
        assert GROMACSEquilibrationManager._classify_moltype("POPC") == "lipid"
        assert GROMACSEquilibrationManager._classify_moltype("LIG") == "other"

    def test_posres_macro_name(self):
        assert GROMACSEquilibrationManager._posres_macro_name("ions") == "POSRES_FC_ION"
        assert (
            GROMACSEquilibrationManager._posres_macro_name("ligand_ABC")
            == "POSRES_FC_LIGAND_ABC"
        )

    def test_parse_topology_and_copies(self, tmp_path):
        top = _minimal_top_with_ions(tmp_path)
        defs = GROMACSEquilibrationManager._parse_moleculetype_defs(top)
        assert [n for n, _ in defs] == ["Protein", "POPC", "NA", "CL", "SOL"]
        assert [c for _, c in defs] == [2, 3, 1, 1, 3]
        _defs, copies = GROMACSEquilibrationManager._build_mol_copies(top)
        assert len(copies) == 1 + 2 + 2 + 2 + 3
        assert copies[0]["start"] == 1 and copies[0]["end"] == 2
        sol0 = next(c for c in copies if c["name"] == "SOL" and c["copy_idx"] == 0)
        assert sol0["n_atoms"] == 3

    def test_write_posre_itp(self, tmp_path):
        out = tmp_path / "posre_ions_na.itp"
        GROMACSEquilibrationManager._write_posre_itp(out, {1}, "POSRES_FC_ION")
        text = out.read_text()
        assert "[ position_restraints ]" in text
        assert "POSRES_FC_ION" in text

    def test_mdp_ions_and_water_macros(self, tmp_path):
        manager = _make_manager(tmp_path)
        params = {
            "name": "Equilibration 1",
            "ensemble": "NPT",
            "time_ns": 0.125,
            "timestep": 1.0,
            "temperature": 310.15,
            "constraints": {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "water": 1.0,
                "ions": 2.5,
                "other": 0.0,
            },
        }
        content = manager.generate_mdp_file(
            stage_name="Eq1",
            stage_params=params,
            stage_index=1,
            scheme_type="NPT",
        )
        assert f"POSRES_FC_WATER={1.0 * 418.4:.1f}" in content
        assert f"POSRES_FC_ION={2.5 * 418.4:.1f}" in content
        assert "define" in content

    def test_mdp_strips_unused_posres_macros(self, tmp_path):
        """Macros not referenced by topology ITPs must not appear in define=."""
        manager = _make_manager(tmp_path)
        params = {
            "name": "Equilibration 1",
            "ensemble": "NPT",
            "time_ns": 0.125,
            "timestep": 1.0,
            "temperature": 310.15,
            "constraints": {
                "protein_backbone": 10.0,
                "protein_sidechain": 5.0,
                "lipid_head": 2.5,
                "lipid_tail": 2.5,
                "water": 0.0,
                "ions": 10.0,
                "other": 0.0,
            },
        }
        content = manager.generate_mdp_file(
            stage_name="Eq1",
            stage_params=params,
            stage_index=1,
            scheme_type="NPT",
            used_posres_macros={
                "POSRES_FC_BB",
                "POSRES_FC_SC",
                "POSRES_FC_LIPID",
                "POSRES_FC_ION",
            },
        )
        assert "POSRES_FC_BB=" in content
        assert "POSRES_FC_ION=" in content
        assert "POSRES_FC_WATER" not in content
        assert "POSRES_FC_OTHER" not in content

    def test_mdp_custom_macro(self, tmp_path):
        manager = _make_manager(tmp_path)
        params = {
            "name": "Equilibration 1",
            "ensemble": "NPT",
            "time_ns": 0.125,
            "timestep": 1.0,
            "temperature": 310.15,
            "constraints": {
                "protein_backbone": 1.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "water": 0.0,
                "ions": 0.0,
                "other": 0.0,
                "ligand_ABC": 3.0,
            },
        }
        content = manager.generate_mdp_file(
            stage_name="Eq1",
            stage_params=params,
            stage_index=1,
            scheme_type="NPT",
        )
        assert f"POSRES_FC_LIGAND_ABC={3.0 * 418.4:.1f}" in content

    def test_add_posres_includes_ions(self, tmp_path):
        top = _minimal_top_with_ions(tmp_path)
        manager = _make_manager(tmp_path)
        ion_itp = tmp_path / "posre_ions_na.itp"
        GROMACSEquilibrationManager._write_posre_itp(ion_itp, {1}, "POSRES_FC_ION")
        bb = tmp_path / "posre_protein_backbone.itp"
        bb.write_text("[ position_restraints ]\n")
        out = manager._add_posres_to_topology(
            topol_path=top,
            posres_files={
                "backbone": bb,
                "sidechain": None,
                "lipid": None,
                "includes": {"NA": [ion_itp], "CL": [ion_itp]},
                "macros": ["POSRES_FC_ION"],
            },
            output_path=tmp_path / "topol_posres.top",
        )
        text = out.read_text()
        assert "posre_ions_na.itp" in text
        assert "posre_protein_backbone.itp" in text

    @pytest.mark.skipif(not MDA_AVAILABLE, reason="MDAnalysis not installed")
    def test_map_selection_partial_copies_warns(self, tmp_path):
        top = _minimal_top_with_ions(tmp_path)
        gro = _minimal_structure_matching_top(tmp_path)
        _defs, copies = GROMACSEquilibrationManager._build_mol_copies(top)
        manager = _make_manager(tmp_path)
        # First NA only (0-based index 8); second NA at 9 → partial
        local_by_mt, partial = manager._map_selection_to_local_indices(
            gro, "index 8", copies
        )
        assert "NA" in local_by_mt
        assert local_by_mt["NA"] == {1}
        assert partial is True

    @pytest.mark.skipif(not MDA_AVAILABLE, reason="MDAnalysis not installed")
    def test_generate_posres_mda_ions_without_gmx(self, tmp_path):
        top = _minimal_top_with_ions(tmp_path)
        gro = _minimal_structure_matching_top(tmp_path)
        ndx = tmp_path / "index.ndx"
        ndx.write_text("[ System ]\n" + " ".join(str(i) for i in range(1, 22)) + "\n")

        manager = GROMACSEquilibrationManager(
            tmp_path, gmx_executable="gmx_missing_xyz"
        )
        result = manager.generate_posres_itp(
            gro_path=gro,
            index_path=ndx,
            output_dir=tmp_path,
            top_path=top,
            pdb_path=None,
            selections=None,
            constraints_max={"ions": 5.0, "water": 1.0},
        )
        assert "POSRES_FC_ION" in result["macros"]
        assert "POSRES_FC_WATER" in result["macros"]
        assert "NA" in result["includes"]
        assert "CL" in result["includes"]
        assert "SOL" in result["includes"]
        ion_files = result["includes"]["NA"]
        assert ion_files and "POSRES_FC_ION" in ion_files[0].read_text()

    def test_pick_protein_lipid_parmed_system_names(self, tmp_path):
        """ParmEd membrane tops use system1/system2 — not Protein/POPC."""
        top = tmp_path / "topol.top"
        top.write_text(
            """\
[ moleculetype ]
system1 3
[ atoms ]
1 X 1 A CA 1 0 12
2 X 1 A CB 2 0 12
3 X 1 A CG 3 0 12
4 X 1 A CD 4 0 12
5 X 1 A CE 5 0 12

[ moleculetype ]
K+ 1
[ atoms ]
1 K 1 K+ K 1 1 39

[ moleculetype ]
system2 3
[ atoms ]
1 L 1 LIP C1 1 0 12
2 L 1 LIP C2 2 0 12
3 L 1 LIP C3 3 0 12

[ moleculetype ]
Cl- 1
[ atoms ]
1 Cl 1 Cl- Cl 1 -1 35

[ moleculetype ]
WAT 2
[ atoms ]
1 OW 1 WAT OW 1 0 16
2 HW 1 WAT HW1 2 0 1
3 HW 1 WAT HW2 3 0 1

[ system ]
x

[ molecules ]
system1  1
K+       2
system1  3
system2  10
K+       1
Cl-      2
WAT      5
"""
        )
        defs = GROMACSEquilibrationManager._parse_moleculetype_defs(top)
        mols = GROMACSEquilibrationManager._parse_molecules_section(top)
        protein, lipid = GROMACSEquilibrationManager._pick_protein_lipid_molnames(
            defs, mols
        )
        assert protein == "system1"
        assert lipid == "system2"

        _defs, copies = GROMACSEquilibrationManager._build_mol_copies(top)
        sys1 = [c for c in copies if c["name"] == "system1"]
        assert len(sys1) == 4
        assert [c["copy_idx"] for c in sys1] == [0, 1, 2, 3]
        assert sys1[0]["end"] == 5
        # first system1 (5) + 2 K+ (2) → second system1 starts at 8
        assert sys1[1]["start"] == 8
    def test_protein_posres_local_indices_multi_copy(self, tmp_path):
        """Backbone POSRES must use local indices for multi-copy system1."""
        top = tmp_path / "topol.top"
        top.write_text(
            """\
[ moleculetype ]
system1 3
[ atoms ]
1 X 1 A N  1 0 14
2 X 1 A CA 2 0 12
3 X 1 A C  3 0 12
4 X 1 A O  4 0 16
5 X 1 A CB 5 0 12

[ moleculetype ]
K+ 1
[ atoms ]
1 K 1 K+ K 1 1 39

[ moleculetype ]
system2 3
[ atoms ]
1 L 1 LIP C1 1 0 12
2 L 1 LIP C2 2 0 12

[ system ]
x
[ molecules ]
system1  2
K+       1
system2  3
"""
        )
        # Two protein copies: atoms 1-5 and 6-10. Backbone globals 2,4,7,9
        ndx = tmp_path / "index.ndx"
        ndx.write_text(
            "[ System ]\n"
            + " ".join(str(i) for i in range(1, 18))
            + "\n"
            + "[ Protein ]\n1 2 3 4 5 6 7 8 9 10\n"
            + "[ Protein-H ]\n1 2 3 4 5 6 7 8 9 10\n"
            + "[ C-alpha ]\n2 7\n"
            + "[ Backbone ]\n2 4 7 9\n"
            + "[ MainChain ]\n2 4 7 9\n"
            + "[ MainChain+Cb ]\n2 4 7 9\n"
            + "[ MainChain+H ]\n2 4 7 9\n"
            + "[ SideChain ]\n5 10\n"
        )
        gro = tmp_path / "system.gro"
        n = 2 * 5 + 1 + 3 * 2
        lines = ["test", f"{n:5d}"]
        for i in range(1, n + 1):
            lines.append(
                f"{1:5d}ALA  CA  {i:5d}   {float(i):8.3f}   {0.0:8.3f}   {0.0:8.3f}"
            )
        lines.append("  10.0 10.0 10.0")
        gro.write_text("\n".join(lines) + "\n")

        manager = GROMACSEquilibrationManager(
            tmp_path, gmx_executable="gmx_missing_xyz"
        )
        result = manager.generate_posres_itp(
            gro_path=gro,
            index_path=ndx,
            output_dir=tmp_path,
            top_path=top,
            constraints_max={"protein_backbone": 1.0, "protein_sidechain": 1.0},
        )
        assert result["backbone"] is not None
        bb = result["backbone"].read_text()
        assert "     2     1" in bb or re.search(r"^\s*2\s+1\s", bb, re.M)
        assert "     4     1" in bb or re.search(r"^\s*4\s+1\s", bb, re.M)
        # Must NOT contain globals from the second copy (7, 9)
        assert not re.search(r"^\s*7\s+1\s", bb, re.M)
        assert not re.search(r"^\s*9\s+1\s", bb, re.M)

        sc = result["sidechain"].read_text()
        assert re.search(r"^\s*5\s+1\s", sc, re.M)
        assert not re.search(r"^\s*10\s+1\s", sc, re.M)

        (tmp_path / "posre_lipid.itp").write_text("[ position_restraints ]\n")
        (tmp_path / "posre_ions_k.itp").write_text("[ position_restraints ]\n")
        out = manager._add_posres_to_topology(
            topol_path=top,
            posres_files={
                "backbone": result["backbone"],
                "sidechain": result["sidechain"],
                "lipid": tmp_path / "posre_lipid.itp",
                "includes": {"K+": [tmp_path / "posre_ions_k.itp"]},
                "macros": [],
            },
            output_path=tmp_path / "topol_posres.top",
        )
        text = out.read_text()
        # Protein includes in system1 block; lipid in system2; ions in K+
        s1 = text.split("[ moleculetype ]")[1]
        s_k = text.split("[ moleculetype ]")[2]
        s2 = text.split("[ moleculetype ]")[3]
        assert "posre_protein_backbone.itp" in s1
        assert "posre_lipid.itp" not in s_k
        assert "posre_ions_k.itp" in s_k
        assert "posre_lipid.itp" in s2
        assert "posre_protein_backbone.itp" not in s2
