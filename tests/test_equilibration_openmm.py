#!/usr/bin/env python3
"""
OpenMM Equilibration Test Suite

Tests for OpenMMEquilibrationManager: configuration generation, restraint file
generation, run script generation, and full workflow.

Usage:
    pytest tests/test_equilibration_openmm.py -v
    pytest tests/test_equilibration_openmm.py::TestOpenMMConfigGeneration -v
"""

import pytest
import sys
import os
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.tools.equilibration import OpenMMEquilibrationManager
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer

POPC_DIR = Path(__file__).parent / "equilibration_examples" / "popc_membrane"

# ============================================================================
# SECTION 1: CLASS CONSTANTS AND MAPPINGS
# ============================================================================


class TestOpenMMEquilibrationManager:
    """Test class constants and basic instance creation."""

    @pytest.fixture
    def manager(self, tmp_path):
        return OpenMMEquilibrationManager(working_dir=tmp_path)

    def test_scheme_mapping_keys(self, manager):
        assert set(manager.SCHEME_MAPPING.keys()) == {"NVT", "NPT", "NPAT", "NPgT"}

    def test_scheme_mapping_values(self, manager):
        assert manager.SCHEME_MAPPING["NVT"] == "01_NVT"
        assert manager.SCHEME_MAPPING["NPT"] == "02_NPT"
        assert manager.SCHEME_MAPPING["NPAT"] == "03_NPAT"
        assert manager.SCHEME_MAPPING["NPgT"] == "04_NPgT"

    def test_template_mapping_count(self, manager):
        assert len(manager.TEMPLATE_MAPPING) == 7

    def test_template_mapping_step1(self, manager):
        assert manager.TEMPLATE_MAPPING["step1"] == "step6.1_equilibration.inp"

    def test_template_mapping_production(self, manager):
        assert manager.TEMPLATE_MAPPING["step7_production"] == "step7_production.inp"

    def test_stage_index_to_key(self, manager):
        assert manager.STAGE_INDEX_TO_KEY[1] == "step1"
        assert manager.STAGE_INDEX_TO_KEY[6] == "step6"
        assert manager.STAGE_INDEX_TO_KEY[7] == "step7_production"

    def test_templates_dir_exists(self, manager):
        assert (
            manager.templates_dir.exists()
        ), f"Templates dir missing: {manager.templates_dir}"

    def test_scripts_dir_exists(self, manager):
        assert (
            manager.scripts_dir.exists()
        ), f"Scripts dir missing: {manager.scripts_dir}"

    def test_all_templates_present(self, manager):
        for ensemble in ("01_NVT", "02_NPT", "03_NPAT", "04_NPgT"):
            for step in (
                "step6.1_equilibration.inp",
                "step6.2_equilibration.inp",
                "step6.3_equilibration.inp",
                "step6.4_equilibration.inp",
                "step6.5_equilibration.inp",
                "step6.6_equilibration.inp",
                "step7_production.inp",
            ):
                p = manager.templates_dir / ensemble / step
                assert p.exists(), f"Template missing: {p}"

    def test_all_scripts_present(self, manager):
        for script in (
            "openmm_run.py",
            "omm_readinputs.py",
            "omm_readparams.py",
            "omm_restraints.py",
            "omm_barostat.py",
            "omm_vfswitch.py",
            "omm_rewrap.py",
        ):
            p = manager.scripts_dir / script
            assert p.exists(), f"Script missing: {p}"

    def test_config_filename_stages_1_to_6(self, manager):
        for i in range(1, 7):
            assert manager._get_config_filename(i) == f"step{i}_equilibration.inp"

    def test_config_filename_stage_7(self, manager):
        assert manager._get_config_filename(7) == "step7_production.inp"

    def test_config_filename_beyond_7(self, manager):
        assert manager._get_config_filename(8) == "step7_production.inp"


# ============================================================================
# SECTION 2: CONFIG GENERATION
# ============================================================================


class TestOpenMMConfigGeneration:
    """Test template loading and placeholder substitution."""

    @pytest.fixture
    def manager(self, tmp_path):
        return OpenMMEquilibrationManager(working_dir=tmp_path)

    @pytest.fixture
    def basic_stage_params(self):
        return {
            "name": "Test stage",
            "ensemble": "NPT",
            "time_ns": 0.125,
            "timestep": 1.0,
            "temperature": 310.15,
            "dcd_freq": 5000,
            "minimize_steps": 5000,
            "constraints": {
                "protein_backbone": 10.0,
                "protein_sidechain": 5.0,
                "lipid_head": 2.5,
            },
        }

    def test_temperature_substituted(self, manager, basic_stage_params):
        content = manager.generate_openmm_config("s1", basic_stage_params, 1, "NPT")
        assert "310.15" in content
        assert "{TEMPERATURE}" not in content

    def test_nstep_calculated_1fs(self, manager, basic_stage_params):
        # time_ns=0.125, dt=1fs=0.001ps → nstep = 0.125*1000/0.001 = 125000
        content = manager.generate_openmm_config("s1", basic_stage_params, 1, "NPT")
        assert "nstep       = 125000" in content
        assert "{NSTEP}" not in content

    def test_nstep_calculated_2fs(self, manager):
        params = {
            "time_ns": 0.5,
            "timestep": 2.0,
            "temperature": 303.15,
            "minimize_steps": 5000,
            "constraints": {},
        }
        content = manager.generate_openmm_config("s1", params, 1, "NVT")
        # 0.5*1000/0.002 = 250000
        assert "nstep       = 250000" in content

    def test_dt_substituted_1fs(self, manager, basic_stage_params):
        content = manager.generate_openmm_config("s1", basic_stage_params, 1, "NPT")
        assert "dt          = 0.001" in content
        assert "{DT}" not in content

    def test_dt_substituted_2fs(self, manager):
        params = {
            "time_ns": 0.5,
            "timestep": 2.0,
            "temperature": 303.15,
            "minimize_steps": 5000,
            "constraints": {},
        }
        content = manager.generate_openmm_config("s1", params, 1, "NVT")
        assert "dt          = 0.002" in content

    def test_rest_yes_when_constraints_active(self, manager, basic_stage_params):
        # Restraints are now enabled when constraints dict has non-zero forces
        content = manager.generate_openmm_config("s1", basic_stage_params, 1, "NPT")
        assert "rest        = yes" in content

    def test_rest_no_when_no_constraints(self, manager):
        # rest = no when all constraint forces are zero
        params = {
            "ensemble": "NPT",
            "time_ns": 0.125,
            "temperature": 310.15,
            "timestep": 2.0,
            "dcd_freq": 5000,
            "minimize_steps": 5000,
            "constraints": {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
            },
        }
        content = manager.generate_openmm_config("s1", params, 1, "NPT")
        assert "rest        = no" in content

    def test_rest_yes_for_production_when_constraints_active(
        self, manager, basic_stage_params
    ):
        content = manager.generate_openmm_config("prod", basic_stage_params, 7, "NPT")
        assert "rest        = yes" in content

    def test_rest_no_for_production_when_all_constraints_zero(self, manager):
        params = {
            "ensemble": "NPT",
            "time_ns": 0.125,
            "temperature": 310.15,
            "timestep": 2.0,
            "dcd_freq": 5000,
            "constraints": {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
            },
        }
        content = manager.generate_openmm_config("prod", params, 7, "NPT")
        assert "rest        = no" in content

    def test_rest_yes_for_ions_constraint_via_custom_restraints(
        self, manager, tmp_path
    ):
        params = {
            "ensemble": "NPT",
            "time_ns": 0.125,
            "temperature": 310.15,
            "timestep": 2.0,
            "dcd_freq": 5000,
            "constraints": {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "ions": 5.0,
            },
        }
        # Simulate a per-stage custom_pos file (the path just needs to exist for the test)
        fake_custom_file = tmp_path / "restraints" / "custom_pos_stage1.txt"
        fake_custom_file.parent.mkdir(parents=True)
        fake_custom_file.touch()
        content = manager.generate_openmm_config(
            "eq1", params, 1, "NPT", custom_pos_file=fake_custom_file
        )
        assert "rest        = yes" in content
        assert "custom_pos_file = restraints/custom_pos_stage1.txt" in content

    def test_fc_ldih_always_zero(self, manager, basic_stage_params):
        # step6.6 (stage 6) drops all lipid restraint fields per CHARMM-GUI protocol;
        # only check stages 1-5 where fc_ldih appears in the template.
        for ensemble in ("NVT", "NPT", "NPAT", "NPgT"):
            for stage_idx in range(1, 6):
                content = manager.generate_openmm_config(
                    "s", basic_stage_params, stage_idx, ensemble
                )
                assert "fc_ldih     = 0" in content
                # Should NOT have non-zero fc_ldih
                for line in content.splitlines():
                    if "fc_ldih" in line and "#" not in line.split("fc_ldih")[0]:
                        assert "= 0" in line, f"fc_ldih should be 0, got: {line}"

    def test_no_unresolved_placeholders(self, manager, basic_stage_params):
        for ensemble in ("NVT", "NPT", "NPAT", "NPgT"):
            for stage_idx in (1, 2, 3, 4, 5, 6, 7):
                content = manager.generate_openmm_config(
                    "s", basic_stage_params, stage_idx, ensemble
                )
                # Check that common placeholders are resolved
                assert "{TEMPERATURE}" not in content
                assert "{NSTEP}" not in content
                assert "{DT}" not in content
                assert "{NSTDCD}" not in content
                assert "{REST}" not in content

    def test_first_stage_has_mini_nstep(self, manager, basic_stage_params):
        content = manager.generate_openmm_config("s1", basic_stage_params, 1, "NPT")
        assert "mini_nstep  = 5000" in content
        assert "gen_vel     = yes" in content
        assert "{MINI_NSTEP}" not in content
        assert "{GEN_VEL}" not in content

    def test_second_stage_no_mini(self, manager, basic_stage_params):
        content = manager.generate_openmm_config("s2", basic_stage_params, 2, "NPT")
        assert "mini_nstep" not in content
        assert "gen_vel" not in content

    def test_nvt_no_pressure_coupling(self, manager, basic_stage_params):
        for stage_idx in range(1, 8):
            content = manager.generate_openmm_config(
                "s", basic_stage_params, stage_idx, "NVT"
            )
            assert "pcouple     = no" in content
            assert "pcouple     = yes" not in content

    def test_npt_pressure_from_step3(self, manager, basic_stage_params):
        # Steps 1 and 2: no pressure
        for stage_idx in (1, 2):
            content = manager.generate_openmm_config(
                "s", basic_stage_params, stage_idx, "NPT"
            )
            assert "pcouple     = no" in content

        # Steps 3+: membrane barostat
        for stage_idx in (3, 4, 5, 6, 7):
            content = manager.generate_openmm_config(
                "s", basic_stage_params, stage_idx, "NPT"
            )
            assert "pcouple     = yes" in content
            assert "p_type      = membrane" in content

    def test_npat_anisotropic_barostat(self, manager, basic_stage_params):
        for stage_idx in (3, 4, 5, 6, 7):
            content = manager.generate_openmm_config(
                "s", basic_stage_params, stage_idx, "NPAT"
            )
            assert "p_type      = anisotropic" in content
            assert "p_scale     = Z" in content

    def test_npgt_membrane_barostat_with_ptens(self, manager, basic_stage_params):
        for stage_idx in (3, 4, 5, 6, 7):
            content = manager.generate_openmm_config(
                "s", basic_stage_params, stage_idx, "NPgT"
            )
            assert "p_type      = membrane" in content
            assert "p_tens" in content

    def test_invalid_scheme_raises(self, manager, basic_stage_params):
        with pytest.raises((ValueError, FileNotFoundError, KeyError)):
            manager.generate_openmm_config("s1", basic_stage_params, 1, "INVALID")

    def test_dcd_freq_used(self, manager, basic_stage_params):
        params = dict(basic_stage_params, dcd_freq=10000)
        content = manager.generate_openmm_config("s1", params, 1, "NPT")
        assert "nstdcd      = 10000" in content

    def test_production_default_dcd_freq(self, manager):
        params = {
            "time_ns": 1.0,
            "timestep": 2.0,
            "temperature": 303.15,
            "constraints": {},
        }
        content = manager.generate_openmm_config("prod", params, 7, "NPT")
        assert "nstdcd      = 50000" in content


# ============================================================================
# SECTION 4: RUN SCRIPT GENERATION
# ============================================================================


class TestOpenMMRunScript:
    """Test bash run script generation."""

    @pytest.fixture
    def manager(self, tmp_path):
        return OpenMMEquilibrationManager(working_dir=tmp_path)

    def test_run_script_created(self, manager, tmp_path):
        names = ["step1_equilibration", "step2_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "system.prmtop", "system.inpcrd"
        )
        assert script.exists()
        assert script.name == "run_equilibration.sh"

    def test_run_script_executable(self, manager, tmp_path):
        names = ["step1_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "system.prmtop", "system.inpcrd"
        )
        assert os.access(script, os.X_OK)

    def test_run_script_has_shebang(self, manager, tmp_path):
        names = ["step1_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "system.prmtop", "system.inpcrd"
        )
        content = script.read_text()
        assert content.startswith("#!/bin/bash")

    def test_run_script_first_stage_no_irst(self, manager, tmp_path):
        names = ["step1_equilibration", "step2_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "system.prmtop", "system.inpcrd"
        )
        lines = script.read_text().splitlines()
        # Find the line with step1
        step1_lines = [l for l in lines if "step1_equilibration.inp" in l]
        assert step1_lines
        assert "-irst" not in step1_lines[0]

    def test_run_script_later_stages_have_irst(self, manager, tmp_path):
        names = ["step1_equilibration", "step2_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "system.prmtop", "system.inpcrd"
        )
        content = script.read_text()
        # Stage 2 should restart from stage 1
        assert "-irst step1_equilibration.rst" in content

    def test_run_script_contains_prmtop(self, manager, tmp_path):
        names = ["step1_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "custom.prmtop", "custom.inpcrd"
        )
        content = script.read_text()
        assert "custom.prmtop" in content
        assert "custom.inpcrd" in content

    def test_run_script_error_handling(self, manager, tmp_path):
        names = ["step1_equilibration", "step2_equilibration"]
        script = manager.generate_run_script(
            names, tmp_path, "system.prmtop", "system.inpcrd"
        )
        content = script.read_text()
        assert "exit 1" in content


# ============================================================================
# SECTION 5: FULL WORKFLOW
# ============================================================================


class TestOpenMMSetup:
    """Integration tests for setup_openmm_equilibration."""

    @pytest.fixture
    def manager(self, tmp_path):
        return OpenMMEquilibrationManager(working_dir=tmp_path)

    @pytest.fixture
    def mock_system_files(self, tmp_path):
        """Create minimal stub system files in tmp_path."""
        prmtop = tmp_path / "system.prmtop"
        inpcrd = tmp_path / "system.inpcrd"
        pdb = tmp_path / "system.pdb"
        prmtop.write_text("# stub prmtop")
        inpcrd.write_text("# stub inpcrd")
        # Minimal PDB with one atom (no real protein/lipid atoms)
        pdb.write_text(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
            "END\n"
        )
        return {
            "prmtop": str(prmtop),
            "inpcrd": str(inpcrd),
            "pdb": str(pdb),
        }

    @pytest.fixture
    def single_stage_nvt(self):
        return [
            {
                "name": "NVT Stage 1",
                "ensemble": "NVT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 310.15,
                "minimize_steps": 5000,
                "constraints": {
                    "protein_backbone": 10.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 0.0,
                },
            }
        ]

    def test_output_directory_created(
        self, manager, tmp_path, mock_system_files, single_stage_nvt
    ):
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=single_stage_nvt,
            output_name="test_out",
        )
        assert result["openmm_dir"].exists()
        assert result["openmm_dir"].name == "test_out"

    def test_config_file_generated(
        self, manager, tmp_path, mock_system_files, single_stage_nvt
    ):
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=single_stage_nvt,
            output_name="test_out",
        )
        assert len(result["config_files"]) == 1
        assert result["config_files"][0].exists()

    def test_run_script_generated(
        self, manager, tmp_path, mock_system_files, single_stage_nvt
    ):
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=single_stage_nvt,
            output_name="test_out",
        )
        assert result["run_script"].exists()
        assert result["run_script"].name == "run_equilibration.sh"

    def test_python_scripts_copied(
        self, manager, tmp_path, mock_system_files, single_stage_nvt
    ):
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=single_stage_nvt,
            output_name="test_out",
        )
        openmm_dir = result["openmm_dir"]
        assert (openmm_dir / "openmm_run.py").exists()
        assert (openmm_dir / "omm_readinputs.py").exists()
        assert (openmm_dir / "omm_restraints.py").exists()

    def test_system_files_copied(
        self, manager, tmp_path, mock_system_files, single_stage_nvt
    ):
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=single_stage_nvt,
            output_name="test_out",
        )
        openmm_dir = result["openmm_dir"]
        assert (openmm_dir / "system.prmtop").exists()
        assert (openmm_dir / "system.inpcrd").exists()

    def test_restraints_dir_created_when_active(
        self, manager, tmp_path, mock_system_files, single_stage_nvt
    ):
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=single_stage_nvt,
            output_name="test_out",
        )
        # restraints/ should be created because protein_backbone = 10.0 > 0
        restraints_dir = result["openmm_dir"] / "restraints"
        assert restraints_dir.exists()
        assert (restraints_dir / "prot_pos.txt").exists()

    def test_no_restraints_dir_when_all_zero(
        self, manager, tmp_path, mock_system_files
    ):
        stages = [
            {
                "ensemble": "NVT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 310.15,
                "minimize_steps": 5000,
                "constraints": {
                    "protein_backbone": 0.0,
                    "protein_sidechain": 0.0,
                    "lipid_head": 0.0,
                },
            }
        ]
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=stages,
            output_name="test_out",
        )
        restraints_dir = result["openmm_dir"] / "restraints"
        assert not restraints_dir.exists()

    def test_scheme_auto_detected_from_ensemble(
        self, manager, tmp_path, mock_system_files
    ):
        stages = [
            {
                "ensemble": "NPAT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 303.15,
                "minimize_steps": 5000,
                "constraints": {},
            }
        ]
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=stages,
            output_name="test_npat",
        )
        # Config should contain NPAT-specific anisotropic barostat from step 3 onwards
        # Stage 1 is step6.1 which has no pressure → just verify config file exists
        assert result["config_files"][0].exists()

    def test_invalid_scheme_raises(self, manager, mock_system_files):
        stages = [
            {
                "ensemble": "INVALID",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 303.15,
                "constraints": {},
            }
        ]
        with pytest.raises(ValueError, match="Unknown scheme_type"):
            manager.setup_openmm_equilibration(
                system_files=mock_system_files,
                stage_params_list=stages,
            )

    def test_empty_stage_list_uses_defaults(self, manager, mock_system_files):
        # Empty or None stage_params_list now triggers automatic default stages
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=[],
            scheme_type="NPT",
        )
        assert len(result["config_files"]) == 6

    def test_six_stage_protocol(self, manager, tmp_path, mock_system_files):
        stages = [
            {
                "ensemble": "NPT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 310.15,
                "minimize_steps": 5000,
                "constraints": {"protein_backbone": 10.0},
            },
            {
                "ensemble": "NPT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 310.15,
                "constraints": {"protein_backbone": 5.0},
            },
            {
                "ensemble": "NPT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 310.15,
                "constraints": {"protein_backbone": 2.0},
            },
            {
                "ensemble": "NPT",
                "time_ns": 0.25,
                "timestep": 2.0,
                "temperature": 310.15,
                "constraints": {"protein_backbone": 1.0},
            },
            {
                "ensemble": "NPT",
                "time_ns": 0.25,
                "timestep": 2.0,
                "temperature": 310.15,
                "constraints": {"protein_backbone": 0.5},
            },
            {
                "ensemble": "NPT",
                "time_ns": 0.25,
                "timestep": 2.0,
                "temperature": 310.15,
                "constraints": {},
            },
        ]
        result = manager.setup_openmm_equilibration(
            system_files=mock_system_files,
            stage_params_list=stages,
            output_name="six_stage",
        )
        assert len(result["config_files"]) == 6
        expected = [f"step{i}_equilibration.inp" for i in range(1, 7)]
        for cfg, exp in zip(result["config_files"], expected):
            assert cfg.name == exp

    @pytest.mark.skipif(
        not POPC_DIR.exists(), reason="popc_membrane test fixture not found"
    )
    def test_full_workflow_with_popc_system(self, tmp_path):
        """Full workflow test using the actual POPC membrane test system."""
        manager = OpenMMEquilibrationManager(working_dir=POPC_DIR)
        system_files = manager.find_system_files()
        assert (
            system_files is not None
        ), "Could not find system files in popc_membrane dir"

        stages = [
            {
                "name": "Eq1",
                "ensemble": "NPT",
                "time_ns": 0.125,
                "timestep": 1.0,
                "temperature": 310.15,
                "minimize_steps": 5000,
                "constraints": {
                    "protein_backbone": 10.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 2.5,
                },
            }
        ]

        result = manager.setup_openmm_equilibration(
            system_files=system_files,
            stage_params_list=stages,
            output_name=str(tmp_path / "openmm_test_popc"),
        )
        assert result["openmm_dir"].exists()
        config = result["config_files"][0]
        content = config.read_text()
        assert "temp        = 310.15" in content
        assert "nstep       = 125000" in content
        assert "fc_ldih     = 0" in content


# ============================================================================
# SECTION 6: FIND_SYSTEM_FILES
# ============================================================================


class TestFindSystemFiles:
    """Test auto-detection of system files."""

    def test_finds_prmtop(self, tmp_path):
        (tmp_path / "system.prmtop").write_text("stub")
        (tmp_path / "system.inpcrd").write_text("stub")
        (tmp_path / "system.pdb").write_text("stub")
        manager = OpenMMEquilibrationManager(working_dir=tmp_path)
        files = manager.find_system_files()
        assert files is not None
        assert "prmtop" in files

    def test_finds_inpcrd_extension(self, tmp_path):
        (tmp_path / "system.prmtop").write_text("stub")
        (tmp_path / "system.inpcrd").write_text("stub")
        manager = OpenMMEquilibrationManager(working_dir=tmp_path)
        files = manager.find_system_files()
        assert files is not None
        assert files["inpcrd"].endswith(".inpcrd")

    def test_finds_crd_fallback(self, tmp_path):
        (tmp_path / "system.prmtop").write_text("stub")
        (tmp_path / "system.crd").write_text("stub")
        manager = OpenMMEquilibrationManager(working_dir=tmp_path)
        files = manager.find_system_files()
        assert files is not None
        assert files["inpcrd"].endswith(".crd")

    def test_returns_none_if_no_prmtop(self, tmp_path):
        manager = OpenMMEquilibrationManager(working_dir=tmp_path)
        files = manager.find_system_files()
        assert files is None


# ============================================================================
# SECTION 7: OpenMMLogAnalyzer
# ============================================================================

# Minimal synthetic OpenMM log content (StateDataReporter tab-separated format)
_LOG_HEADER = '#"Step"\t"Time (ps)"\t"Potential Energy (kJ/mole)"\t"Kinetic Energy (kJ/mole)"\t"Total Energy (kJ/mole)"\t"Temperature (K)"\t"Box Volume (nm^3)"\t"Density (g/mL)"\t"Progress (%)"\t"Remaining Time"\t"Speed (ns/day)"'
_LOG_ROWS = [
    "1000\t2.0\t-500000.0\t100000.0\t-400000.0\t303.0\t750.0\t0.95\t1.0\t1:00\t10.0",
    "2000\t4.0\t-501000.0\t101000.0\t-400000.0\t305.0\t748.0\t0.96\t2.0\t0:59\t10.1",
    "3000\t6.0\t-502000.0\t102000.0\t-400000.0\t307.0\t746.0\t0.97\t3.0\t0:58\t10.2",
]
_LOG_CONTENT = "\n".join([_LOG_HEADER] + _LOG_ROWS) + "\n"


@pytest.fixture
def sample_log(tmp_path):
    """Write a synthetic OpenMM log file and return its path."""
    log = tmp_path / "step1_equilibration.log"
    log.write_text(_LOG_CONTENT)
    return log


@pytest.fixture
def two_log_files(tmp_path):
    """Two synthetic log files for multi-file parsing tests."""
    log1 = tmp_path / "step1_equilibration.log"
    log2 = tmp_path / "step2_equilibration.log"
    log1.write_text(_LOG_CONTENT)
    # Second log: shift steps and time
    rows2 = [
        "4000\t8.0\t-503000.0\t103000.0\t-400000.0\t309.0\t744.0\t0.98\t4.0\t0:57\t10.3",
        "5000\t10.0\t-504000.0\t104000.0\t-400000.0\t311.0\t742.0\t0.99\t5.0\t0:56\t10.4",
    ]
    log2.write_text("\n".join([_LOG_HEADER] + rows2) + "\n")
    return log1, log2


class TestOpenMMLogAnalyzerParsing:
    """Test log file parsing."""

    def test_parses_potential(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert len(ana.data["potential"]) == 3
        assert ana.data["potential"][0] == pytest.approx(-500000.0)

    def test_parses_kinetic(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert len(ana.data["kinetic"]) == 3
        assert ana.data["kinetic"][1] == pytest.approx(101000.0)

    def test_parses_total(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert all(abs(v - (-400000.0)) < 1.0 for v in ana.data["total"])

    def test_parses_temperature(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert ana.data["temp"][0] == pytest.approx(303.0)

    def test_parses_volume(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert ana.data["volume"][0] == pytest.approx(750.0)

    def test_parses_density(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert ana.data["density"][0] == pytest.approx(0.95)

    def test_parses_time_ps(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert ana.data["time_ps"] == pytest.approx([2.0, 4.0, 6.0])

    def test_parses_step(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        assert ana.data["step"] == pytest.approx([1000.0, 2000.0, 3000.0])

    def test_missing_file_returns_empty(self, tmp_path):
        ana = OpenMMLogAnalyzer(tmp_path / "nonexistent.log")
        assert len(ana.data["step"]) == 0

    def test_multi_file_concatenation(self, two_log_files):
        log1, log2 = two_log_files
        ana = OpenMMLogAnalyzer([log1, log2])
        assert len(ana.data["step"]) == 5  # 3 from file1 + 2 from file2


class TestOpenMMLogAnalyzerStatistics:
    """Test get_statistics()."""

    def test_statistics_keys(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        stats = ana.get_statistics()
        for key in ("potential", "kinetic", "total", "temp", "volume", "density"):
            assert key in stats

    def test_statistics_mean_potential(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        stats = ana.get_statistics()
        expected_mean = (-500000.0 + -501000.0 + -502000.0) / 3
        assert stats["potential"]["mean"] == pytest.approx(expected_mean, rel=1e-6)

    def test_statistics_initial_final(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        stats = ana.get_statistics()
        assert stats["temp"]["initial"] == pytest.approx(303.0)
        assert stats["temp"]["final"] == pytest.approx(307.0)

    def test_statistics_min_max(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        stats = ana.get_statistics()
        assert stats["volume"]["min"] == pytest.approx(746.0)
        assert stats["volume"]["max"] == pytest.approx(750.0)

    def test_step_not_in_statistics(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        stats = ana.get_statistics()
        assert "step" not in stats


class TestOpenMMLogAnalyzerTimeArray:
    """Test _calculate_time_array()."""

    def test_time_from_time_ps_column(self, sample_log):
        import numpy as np

        ana = OpenMMLogAnalyzer(sample_log)
        t = ana._calculate_time_array()
        assert t == pytest.approx(np.array([2.0, 4.0, 6.0]) / 1000.0, rel=1e-6)

    def test_time_in_nanoseconds(self, sample_log):
        ana = OpenMMLogAnalyzer(sample_log)
        t = ana._calculate_time_array()
        assert t[0] == pytest.approx(0.002)  # 2 ps = 0.002 ns


class TestOpenMMLogAnalyzerPlot:
    """Test plot methods (output only; no matplotlib display)."""

    def test_plot_energy_creates_file(self, sample_log, tmp_path):
        ana = OpenMMLogAnalyzer(sample_log)
        out = str(tmp_path / "energy.png")
        ana.plot_energy(save=out, show=False)
        assert Path(out).exists()

    def test_plot_properties_creates_file(self, sample_log, tmp_path):
        ana = OpenMMLogAnalyzer(sample_log)
        out = str(tmp_path / "props.png")
        ana.plot_properties(properties=["potential", "temp"], save=out, show=False)
        assert Path(out).exists()

    def test_plot_energy_kcal(self, sample_log, tmp_path):
        ana = OpenMMLogAnalyzer(sample_log)
        out = str(tmp_path / "energy_kcal.png")
        ana.plot_energy(energy_units="kcal/mol", save=out, show=False)
        assert Path(out).exists()

    def test_plot_empty_log_no_crash(self, tmp_path):
        log = tmp_path / "empty.log"
        log.write_text("")
        ana = OpenMMLogAnalyzer(log)
        # Should not raise; just logs a warning
        ana.plot_energy(show=False)

    def test_column_map_has_modern_names(self):
        from gatewizard.utils.openmm_analysis import _COLUMN_MAP

        assert "Potential Energy (kJ/mole)" in _COLUMN_MAP
        assert "Temperature (K)" in _COLUMN_MAP
        assert "Box Volume (nm^3)" in _COLUMN_MAP

    def test_column_map_legacy_names(self):
        from gatewizard.utils.openmm_analysis import _COLUMN_MAP

        assert "Potential Energy (kJ/mol)" in _COLUMN_MAP


# ============================================================================
# STANDALONE RUNNER
# ============================================================================


if __name__ == "__main__":
    import subprocess

    result = subprocess.run(
        ["python", "-m", "pytest", __file__, "-v"],
        cwd=str(Path(__file__).parent.parent),
    )
    sys.exit(result.returncode)
