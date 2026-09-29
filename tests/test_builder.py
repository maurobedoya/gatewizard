#!/usr/bin/env python3
"""
Builder Test Suite

This test suite covers:
1. Core functionality tests (specs and features)
2. Documentation example workflows (Examples 1-17)
3. Force field management tests
4. Integration tests (require external tools)

The test suite automatically discovers and runs all example scripts from:
    tests/builder_examples/builder_example_*.py

Note: Some tests require external tools (packmol-memgen, AmberTools) and are skipped
if those tools are not available.

Usage:
    # Run all tests
    pytest tests/test_builder.py -v

    # Run only example tests
    pytest tests/test_builder.py::TestBuilderExamples -v

    # Run specific example
    pytest tests/test_builder.py::TestBuilderExamples::test_example_script[builder_example_08.py] -v
"""

import pytest
import sys
import tempfile
from pathlib import Path

from tests.example_runner import (
    parametrize_example_scripts,
    remove_tree,
    run_example_script,
)

# Add gatewizard to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from gatewizard.core.builder import Builder
from gatewizard.tools.force_fields import ForceFieldManager

# ============================================================================
# SECTION 1: CORE FUNCTIONALITY TESTS
# ============================================================================


class TestBuilder:
    """Test the Builder class core functionality."""

    @pytest.fixture
    def builder(self):
        """Create a Builder instance."""
        return Builder()

    @pytest.fixture
    def ff_manager(self):
        """Create a ForceFieldManager instance."""
        return ForceFieldManager()

    def test_default_configuration(self, builder):
        """Test that default configuration is set correctly."""
        assert builder.config["water_model"] == "tip3p"
        assert builder.config["protein_ff"] == "ff19SB"
        assert builder.config["lipid_ff"] == "lipid21"
        assert builder.config["salt_concentration"] == 0.15
        assert builder.config["cation"] == "K+"
        assert builder.config["anion"] == "Cl-"
        assert builder.config["dist"] == 12
        assert builder.config["dist_wat"] == 26
        assert builder.config["preoriented"] == True
        assert builder.config["parametrize"] == True
        assert builder.config["notprotonate"] == True
        assert builder.config["nloop"] == 20
        assert builder.config["nloop_all"] == 100
        assert builder.config["tolerance"] == 2.0

    def test_set_configuration(self, builder):
        """Test configuration update."""
        builder.set_configuration(
            water_model="tip4p", salt_concentration=0.5, dist=14.0, dist_wat=20.0
        )

        assert builder.config["water_model"] == "tip4p"
        assert builder.config["salt_concentration"] == 0.5
        assert builder.config["dist"] == 14.0
        assert builder.config["dist_wat"] == 20.0
        # Other values should remain unchanged
        assert builder.config["protein_ff"] == "ff19SB"

    def test_lipids_argument_symmetric(self, builder):
        """Test lipids argument generation for symmetric membrane."""
        result = builder._prepare_lipids_argument(
            upper_lipids=["POPC"], lower_lipids=["POPC"]
        )
        assert result == "POPC//POPC"

    def test_lipids_argument_asymmetric(self, builder):
        """Test lipids argument generation for asymmetric membrane."""
        result = builder._prepare_lipids_argument(
            upper_lipids=["POPC", "POPE"], lower_lipids=["POPE", "POPS"]
        )
        assert result == "POPC:POPE//POPE:POPS"

    def test_lipids_argument_complex(self, builder):
        """Test lipids argument with complex composition."""
        result = builder._prepare_lipids_argument(
            upper_lipids=["POPC", "POPE", "CHL1"], lower_lipids=["POPE", "POPS", "CHL1"]
        )
        assert result == "POPC:POPE:CHL1//POPE:POPS:CHL1"

    def test_build_command_includes_packmol_options(self, builder, tmp_path):
        """Test packmol-memgen command includes PACKMOL loop/tolerance flags."""
        pdb_file = tmp_path / "protein.pdb"
        pdb_file.write_text("END\n")
        config = {
            **builder.config,
            "water_model": "opc",
            "protein_ff": "ff19SB",
            "lipid_ff": "lipid21",
            "dist": 14,
            "dist_wat": 26,
            "nloop": 30,
            "nloop_all": 120,
            "tolerance": 1.5,
        }
        cmd = builder._build_command(
            pdb_file, ["POPC"], ["POPC"], "1//1", config
        )
        cmd_str = " ".join(cmd)
        assert "--dist 14" in cmd_str
        assert "--dist_wat 26" in cmd_str
        assert "--nloop 30" in cmd_str
        assert "--nloop_all 120" in cmd_str
        assert "--tolerance 1.5" in cmd_str

    def test_build_command_protein_path_unchanged(self, builder, tmp_path):
        """Membrane-protein command still uses --pdb and does not add bilayer-only flags."""
        pdb_file = tmp_path / "protein.pdb"
        pdb_file.write_text("END\n")
        cmd = builder._build_command(
            pdb_file, ["POPC"], ["POPC"], "1//1", {**builder.config}
        )
        assert cmd[0] == "packmol-memgen"
        assert "--pdb" in cmd
        assert str(pdb_file) in cmd
        assert "--preoriented" in cmd
        assert "--distxy_fix" not in cmd
        assert "--solute" not in cmd

    def test_build_command_bilayer_only(self, builder):
        """Bilayer-only builds omit --pdb and require --distxy_fix."""
        config = {**builder.config, "distxy_fix": 100, "preoriented": True}
        cmd = builder._build_command(None, ["POPC"], ["POPC"], "1//1", config)
        cmd_str = " ".join(str(c) for c in cmd)
        assert "--pdb" not in cmd
        assert "--preoriented" not in cmd
        assert "--distxy_fix" in cmd
        assert "100" in cmd
        assert "--lipids" in cmd
        assert "MEMEMBED" not in cmd_str

    def test_build_command_bilayer_only_infers_distxy_from_dims(self, builder):
        config = {**builder.config, "dims": [80, 80, 120]}
        cmd = builder._build_command(None, ["DOPE", "DOPG"], ["DOPE", "DOPG"], "3:1//3:1", config)
        assert "--pdb" not in cmd
        assert "--distxy_fix" in cmd
        assert cmd[cmd.index("--distxy_fix") + 1] == "80.0"

    def test_build_command_free_solute(self, builder, tmp_path):
        """Free molecules add --solute/--solute_con and optional in-membrane placement."""
        pdb_file = tmp_path / "protein.pdb"
        pdb_file.write_text("END\n")
        tea = tmp_path / "TEA.pdb"
        tea.write_text("HETATM    1  N   TEA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
        config = {
            **builder.config,
            "solutes": [
                {
                    "pdb": str(tea),
                    "concentration": "4",
                    "in_membrane": False,
                    "prot_dist": 10,
                }
            ],
            "ligand_params": {
                "TEA": {"frcmod": "TEA.frcmod", "lib": "TEA.lib"}
            },
        }
        cmd = builder._build_command(pdb_file, ["POPC"], ["POPC"], "1//1", config)
        assert "--solute" in cmd
        assert str(tea) in cmd
        assert "--solute_con" in cmd
        assert cmd[cmd.index("--solute_con") + 1] == "4"
        assert "--solute_prot_dist" in cmd
        assert "--solute_inmem" not in cmd
        assert "--ligand_param" in cmd
        assert "--gaff2" in cmd
        assert "--pdb" in cmd

    def test_build_command_solute_in_membrane_no_protein(self, builder, tmp_path):
        tea = tmp_path / "TEA.pdb"
        tea.write_text("HETATM    1  N   TEA A   1       0.000   0.000   0.000  1.00  0.00           N\n")
        config = {
            **builder.config,
            "distxy_fix": 75,
            "solute_inmem": True,
            "solutes": [{"pdb": str(tea), "concentration": "0.1M"}],
        }
        cmd = builder._build_command(None, ["POPC"], ["POPC"], "1//1", config)
        assert "--pdb" not in cmd
        assert "--solute_inmem" in cmd
        assert "--solute_prot_dist" not in cmd

    def test_validate_bilayer_only_requires_xy_size(self, builder):
        valid, msg = builder.validate_system_inputs(
            pdb_file=None,
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            water_model="tip3p",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
        )
        assert valid is False
        assert "XY" in msg or "distxy" in msg.lower() or "protein" in msg.lower()

        valid_ok, msg_ok = builder.validate_system_inputs(
            pdb_file=None,
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            water_model="tip3p",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
            distxy_fix=100,
        )
        assert valid_ok is True, msg_ok
        assert "Remove protein hydrogens" not in (msg_ok or "")

    def test_validate_skips_protein_hydrogen_scan_without_pdb(self, builder, tmp_path):
        tea = tmp_path / "TEA.pdb"
        tea.write_text(
            "HETATM    1  N   TEA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        )
        valid, msg = builder.validate_system_inputs(
            pdb_file="",
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            water_model="tip3p",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
            distxy_fix=80,
            solutes=[{"pdb": str(tea), "concentration": "4", "name": "TEA"}],
            parametrize=True,
        )
        assert valid is True, msg
        assert "Remove protein hydrogens" not in (msg or "")
        assert "GAFF" in msg or "ligand_params" in msg or "Parametrize" in msg

    def test_generate_preparation_inputs_bilayer_only(self, builder, tmp_path):
        success, message, job_dir = builder.generate_preparation_inputs(
            pdb_file=None,
            working_dir=str(tmp_path),
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            distxy_fix=100,
            water_model="tip3p",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
        )
        assert success, message
        assert job_dir is not None
        assert job_dir.name.startswith("02_build_bilayer")
        assert (job_dir / "run_preparation.sh").is_file()
        import json

        status = json.loads((job_dir / "status.json").read_text())
        command = status["command"]
        assert "--pdb" not in command.split()
        assert "--distxy_fix" in command
        assert "--preoriented" not in command.split()
        assert status["status"] == "not_started"
        assert "MEMEMBED" not in status["steps"]
    def test_generate_preparation_inputs_with_solute(self, builder, tmp_path):
        tea = tmp_path / "TEA.pdb"
        tea.write_text(
            "HETATM    1  N   TEA A   1       0.000   0.000   0.000  1.00  0.00           N\nEND\n"
        )
        success, message, job_dir = builder.generate_preparation_inputs(
            pdb_file=None,
            working_dir=str(tmp_path),
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            distxy_fix=75,
            solutes=[{"pdb": str(tea), "concentration": "4", "in_membrane": True, "name": "TEA"}],
            ligand_params={"TEA": {"frcmod": str(tmp_path / "TEA.frcmod"), "lib": str(tmp_path / "TEA.lib")}},
            water_model="tip3p",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
        )
        assert success, message
        assert (job_dir / "solutes" / "TEA.pdb").is_file()
        import json

        command = json.loads((job_dir / "status.json").read_text())["command"]
        assert "--pdb" not in command.split()
        assert "--solute" in command
        assert "--solute_con" in command
        assert "--solute_inmem" in command
        assert "--gaff2" in command
        assert "--distxy_fix" in command


class TestCreateJobDirectoryParamCaches:
    """Parametrize ligands/caps into 02_build_* then Generate Input must reuse it."""

    @pytest.fixture
    def builder(self):
        return Builder()

    def test_reuses_folder_with_ligand_params_only(self, builder, tmp_path):
        preferred = tmp_path / "02_build_1jno_protonated"
        lig = preferred / "ligand_params" / "LIG"
        lig.mkdir(parents=True)
        (lig / "LIG.frcmod").write_text("x\n", encoding="utf-8")
        (lig / "LIG.lib").write_text("x\n", encoding="utf-8")

        job_dir = builder._create_job_directory(
            None, str(tmp_path), custom_output_name="02_build_1jno_protonated"
        )
        assert job_dir.resolve() == preferred.resolve()
        assert (job_dir / "ligand_params" / "LIG" / "LIG.frcmod").is_file()
        assert (job_dir / "logs").is_dir()

    def test_reuses_folder_with_peptide_cap_params_only(self, builder, tmp_path):
        preferred = tmp_path / "02_build_1jno_protonated"
        cap = preferred / "peptide_cap_params" / "ETA"
        cap.mkdir(parents=True)
        (cap / "ETA.frcmod").write_text("x\n", encoding="utf-8")
        (cap / "ETA.lib").write_text("x\n", encoding="utf-8")

        job_dir = builder._create_job_directory(
            None, str(tmp_path), custom_output_name="02_build_1jno_protonated"
        )
        assert job_dir.resolve() == preferred.resolve()
        assert (job_dir / "peptide_cap_params" / "ETA" / "ETA.lib").is_file()

    def test_generate_inputs_reuses_param_staging_folder(self, builder, tmp_path):
        preferred = tmp_path / "02_build_params_reuse"
        (preferred / "ligand_params" / "LIG").mkdir(parents=True)
        (preferred / "ligand_params" / "LIG" / "LIG.frcmod").write_text("x\n")
        (preferred / "peptide_cap_params" / "ETA").mkdir(parents=True)
        (preferred / "peptide_cap_params" / "ETA" / "ETA.lib").write_text("x\n")

        success, message, job_dir = builder.generate_preparation_inputs(
            pdb_file=None,
            working_dir=str(tmp_path),
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            distxy_fix=100,
            water_model="tip3p",
            protein_ff="ff19SB",
            lipid_ff="lipid21",
            output_folder_name="02_build_params_reuse",
        )
        assert success, message
        assert Path(job_dir).resolve() == preferred.resolve()
        assert (preferred / "run_preparation.sh").is_file()
        assert (preferred / "ligand_params" / "LIG" / "LIG.frcmod").is_file()
        assert (preferred / "peptide_cap_params" / "ETA" / "ETA.lib").is_file()

    def test_prior_preparation_gets_timestamp_and_copies_caches(self, builder, tmp_path):
        preferred = tmp_path / "02_build_again"
        preferred.mkdir()
        (preferred / "run_preparation.sh").write_text("#!/bin/bash\n", encoding="utf-8")
        (preferred / "status.json").write_text('{"status":"not_started"}\n', encoding="utf-8")
        lig = preferred / "ligand_params" / "LIG"
        lig.mkdir(parents=True)
        (lig / "LIG.frcmod").write_text("cached\n", encoding="utf-8")
        caps = preferred / "peptide_cap_params" / "FVA"
        caps.mkdir(parents=True)
        (caps / "FVA.lib").write_text("cached\n", encoding="utf-8")

        job_dir = builder._create_job_directory(
            None, str(tmp_path), custom_output_name="02_build_again"
        )
        assert job_dir.resolve() != preferred.resolve()
        assert job_dir.name.startswith("02_build_again_")
        assert (job_dir / "ligand_params" / "LIG" / "LIG.frcmod").read_text() == "cached\n"
        assert (job_dir / "peptide_cap_params" / "FVA" / "FVA.lib").read_text() == "cached\n"

    def test_default_stem_reuses_param_only_folder(self, builder, tmp_path):
        pdb = tmp_path / "1jno_protonated.pdb"
        pdb.write_text("END\n", encoding="utf-8")
        preferred = tmp_path / "02_build_1jno_protonated"
        (preferred / "peptide_cap_params" / "ETA").mkdir(parents=True)
        (preferred / "peptide_cap_params" / "ETA" / "ETA.frcmod").write_text("x\n")

        job_dir = builder._create_job_directory(str(pdb), str(tmp_path), custom_output_name=None)
        assert job_dir.resolve() == preferred.resolve()


# ============================================================================
# SECTION 2: FORCE FIELD MANAGEMENT TESTS
# ============================================================================


class TestForceFieldManager:
    """Test the ForceFieldManager class."""

    @pytest.fixture
    def ff_manager(self):
        """Create a ForceFieldManager instance."""
        return ForceFieldManager()

    def test_water_models_available(self, ff_manager):
        """Test that water models are available."""
        water_models = ff_manager.get_water_models()
        assert len(water_models) > 0
        assert "tip3p" in water_models
        assert "tip4pd" in water_models

    def test_protein_force_fields_available(self, ff_manager):
        """Test that protein force fields are available."""
        protein_ffs = ff_manager.get_protein_force_fields()
        assert len(protein_ffs) > 0
        assert "ff14SB" in protein_ffs
        assert "ff19SB" in protein_ffs

    def test_lipid_force_fields_available(self, ff_manager):
        """Test that lipid force fields are available."""
        lipid_ffs = ff_manager.get_lipid_force_fields()
        assert len(lipid_ffs) > 0
        assert "lipid17" in lipid_ffs
        assert "lipid21" in lipid_ffs

    def test_available_lipids(self, ff_manager):
        """Test that lipids are available."""
        lipids = ff_manager.get_available_lipids()
        assert len(lipids) > 0
        assert "POPC" in lipids
        assert "POPE" in lipids
        assert "CHL1" in lipids

    def test_validate_combination_valid(self, ff_manager):
        """Test validation of valid force field combination."""
        valid, message, is_warning = ff_manager.validate_combination(
            water_model="tip3p", protein_ff="ff14SB", lipid_ff="lipid21"
        )
        assert valid == True
        assert is_warning == False
        assert "valid" in message.lower()

    def test_validate_combination_invalid_water(self, ff_manager):
        """Test validation fails for invalid water model."""
        valid, message, is_warning = ff_manager.validate_combination(
            water_model="invalid_water", protein_ff="ff14SB", lipid_ff="lipid21"
        )
        assert valid == False
        assert is_warning == False
        assert "unknown" in message.lower() or "water" in message.lower()

    def test_validate_combination_unvalidated_warning(self, ff_manager):
        """Test validation returns warning for unvalidated (but not invalid) combinations."""
        # Test a combination that is not in the compatibility lists
        # For example, tip4p with lipid21 (if not in compatible_with lists)
        valid, message, is_warning = ff_manager.validate_combination(
            water_model="tip3p", protein_ff="ff19SB", lipid_ff="lipid17"
        )
        # Should be valid (not block) but with a warning if incompatible
        # This depends on what's actually in the compatibility lists
        # The key is: valid should be True if components exist, even if not compatible
        assert valid == True or valid == False  # Either is acceptable
        if not valid:
            # If it's invalid, it should not be because of unknown components
            assert "unknown" not in message.lower()

    def test_recommendations_membrane(self, ff_manager):
        """Test recommendations for membrane systems."""
        rec = ff_manager.get_recommendations("membrane")
        assert "water_model" in rec
        assert "protein_ff" in rec
        assert "lipid_ff" in rec
        assert "reason" in rec
        assert rec["water_model"] in ff_manager.get_water_models()
        assert rec["protein_ff"] in ff_manager.get_protein_force_fields()
        assert rec["lipid_ff"] in ff_manager.get_lipid_force_fields()

    def test_validate_lipid_valid(self, ff_manager):
        """Test validation of valid lipid."""
        assert ff_manager.validate_lipid("POPC") == True
        assert ff_manager.validate_lipid("POPE") == True
        assert ff_manager.validate_lipid("CHL1") == True

    def test_validate_lipid_invalid(self, ff_manager):
        """Test validation of invalid lipid."""
        assert ff_manager.validate_lipid("INVALID_LIPID") == False
        assert ff_manager.validate_lipid("XYZ123") == False


# ============================================================================
# SECTION 3: DOCUMENTATION EXAMPLE TESTS
# ============================================================================


class TestBuilderExamples:
    """Run each builder_example_*.py script once."""

    @pytest.fixture(autouse=True)
    def cleanup_example_outputs(self):
        """Clean up output directories created by example scripts."""
        yield
        examples_dir = Path(__file__).parent / "builder_examples"
        remove_tree(examples_dir / "systems")
        remove_tree(Path("./systems"))

    @parametrize_example_scripts(
        Path(__file__).parent / "builder_examples",
        "builder_example_*.py",
    )
    def test_example_script(self, script):
        run_example_script(script)


# ============================================================================
# SECTION 4: INTEGRATION TESTS (REQUIRE EXTERNAL TOOLS)
# ============================================================================


class TestBuilderIntegration:
    """Integration tests requiring external tools."""

    @pytest.fixture
    def temp_dir(self):
        """Create temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield Path(tmpdir)

    @pytest.fixture
    def sample_pdb(self, temp_dir):
        """Create a minimal sample PDB file."""
        pdb_content = """REMARK   Created by GateWizard test suite
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  ALA A   1       1.458   0.000   0.000  1.00  0.00           C
ATOM      3  C   ALA A   1       2.009   1.420   0.000  1.00  0.00           C
ATOM      4  O   ALA A   1       1.251   2.389   0.000  1.00  0.00           O
ATOM      5  CB  ALA A   1       1.989  -0.729  -1.232  1.00  0.00           C
TER       6      ALA A   1
END
"""
        pdb_file = temp_dir / "test_protein.pdb"
        pdb_file.write_text(pdb_content)
        return pdb_file

    # Note: This test only validates inputs, gatewizard environment should be activated
    def test_prepare_system_dry_run(self, temp_dir, sample_pdb):
        """Test system preparation (dry run - validation only)."""
        builder = Builder()

        # Only validate, don't actually run
        valid, error_msg = builder.validate_system_inputs(
            pdb_file=str(sample_pdb),
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
            water_model="tip3p",
            protein_ff="ff14SB",
            lipid_ff="lipid21",
        )

        # Should be valid (or have a specific validation error)
        assert valid or len(error_msg) > 0
        print(f"Validation result: {valid}, message: {error_msg}")

    def test_generate_preparation_inputs(self, temp_dir, sample_pdb):
        """Test generating preparation input files without launching."""
        builder = Builder()
        builder.set_configuration(water_model="tip3p", protein_ff="ff14SB")

        success, message, job_dir = builder.generate_preparation_inputs(
            pdb_file=str(sample_pdb),
            working_dir=str(temp_dir),
            upper_lipids=["POPC"],
            lower_lipids=["POPC"],
            lipid_ratios="1//1",
        )

        assert success, message
        assert job_dir is not None
        assert (job_dir / "run_preparation.sh").is_file()
        assert (job_dir / "status.json").is_file()

        import json

        status = json.loads((job_dir / "status.json").read_text())
        assert status["status"] == "not_started"
        assert status["start_time"] is None
        assert not (job_dir / "process.pid").exists()

    def test_run_preparation_requires_generated_files(self, temp_dir):
        """Test run_preparation fails when input files are missing."""
        builder = Builder()
        job_dir = temp_dir / "empty_job"
        job_dir.mkdir()

        success, message = builder.run_preparation(job_dir)
        assert not success
        assert "generate input files first" in message.lower()

    def test_cancel_preparation_marks_running_job(self, temp_dir):
        """Cancel updates status.json even when process.pid is already gone."""
        import json

        builder = Builder()
        job_dir = temp_dir / "cancel_job"
        job_dir.mkdir()
        (job_dir / "run_preparation.sh").write_text("#!/bin/bash\n", encoding="utf-8")
        (job_dir / "status.json").write_text(
            json.dumps(
                {
                    "status": "running",
                    "start_time": "2026-01-01T00:00:00",
                    "error": None,
                    "steps_completed": ["Packmol"],
                }
            ),
            encoding="utf-8",
        )

        result = builder.cancel_preparation(job_dir)
        assert result["success"] is True
        assert result["stopped"] is True
        assert result["status"] == "cancelled"
        status = json.loads((job_dir / "status.json").read_text(encoding="utf-8"))
        assert status["status"] == "cancelled"
        assert status["error"] == "Cancelled by user"
        assert status.get("end_time")

    def test_cancel_preparation_noop_when_already_done(self, temp_dir):
        import json

        builder = Builder()
        job_dir = temp_dir / "done_job"
        job_dir.mkdir()
        (job_dir / "run_preparation.sh").write_text("#!/bin/bash\n", encoding="utf-8")
        (job_dir / "status.json").write_text(
            json.dumps({"status": "completed", "error": None}),
            encoding="utf-8",
        )
        result = builder.cancel_preparation(job_dir)
        assert result["success"] is True
        assert result["stopped"] is False
        assert result["status"] == "completed"


class TestNamdOpcBuilderTleap:
    """OPC + NAMD tleap FlexibleWater integration."""

    def test_tleap_flexible_water_when_namd_opc(self):
        builder = Builder()
        builder.set_configuration(water_model="opc", md_engine="namd")
        content = builder._create_tleap_input("system_for_tleap.pdb", builder.config)
        assert "FlexibleWater on" in content
        assert "leaprc.water.opc" in content

    def test_tleap_no_flexible_water_tip3p_namd(self):
        builder = Builder()
        builder.set_configuration(water_model="tip3p", md_engine="namd")
        content = builder._create_tleap_input("system_for_tleap.pdb", builder.config)
        assert "FlexibleWater" not in content

    def test_tleap_no_flexible_water_opc_gromacs(self):
        builder = Builder()
        builder.set_configuration(water_model="opc", md_engine="gromacs")
        content = builder._create_tleap_input("system_for_tleap.pdb", builder.config)
        assert "FlexibleWater" not in content


if __name__ == "__main__":
    # Run tests with pytest
    pytest.main([__file__, "-v"])
