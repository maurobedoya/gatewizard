"""Tests for NAMD GPU-resident template customization."""

from __future__ import annotations

from pathlib import Path

from gatewizard.tools.equilibration import NAMDEquilibrationManager


def _manager() -> NAMDEquilibrationManager:
    return NAMDEquilibrationManager(Path("."), namd_executable="namd3")


def test_gpu_resident_block_empty_when_disabled() -> None:
    mgr = _manager()
    mgr.gpu_resident = False
    assert mgr._gpu_resident_block("Equilibration 1", 0) == ""
    assert mgr._gpu_resident_block("Production", 6) == ""


def test_gpu_resident_block_skipped_on_equilibration_stages() -> None:
    mgr = _manager()
    mgr.gpu_resident = True
    for index, name in enumerate(
        [
            "Equilibration 1",
            "Equilibration 2",
            "Equilibration 3",
            "Equilibration 4",
            "Equilibration 5",
            "Equilibration 6",
        ]
    ):
        block = mgr._gpu_resident_block(name, index)
        assert "GPUresident             on" not in block
        assert "production" in block.lower()


def test_gpu_resident_block_on_production_only() -> None:
    mgr = _manager()
    mgr.gpu_resident = True
    block = mgr._gpu_resident_block("Production", 6)
    assert "GPUresident             on" in block
    assert "CUDASOAIntegrate" not in block


def test_customize_keeps_reassign_on_equilibration() -> None:
    mgr = _manager()
    mgr.gpu_resident = True
    mgr.water_model = "tip3p"
    template = """
{WATER_MODEL_BLOCK}
{GPU_RESIDENT_BLOCK}
reassignFreq            500;
reassignTemp            $temp;
{INITIAL_TEMPERATURE_DIRECTIVE}
{CELL_BASIS_VECTORS}
{PME_SETTINGS}
{RESTRAINTS_BLOCK}
margin                  {MARGIN};
set time {TIME_NS};
set tstep {TIMESTEP};
timestep           {TIMESTEP}
minimize                {MINIMIZE_STEPS}
run                     [expr int($time * 1e6 / $tstep)]
"""
    stage0 = mgr._customize_charmm_gui_template(
        template,
        stage_name="Equilibration 1",
        stage_params={"temperature": 303.15, "time_ns": 0.125, "timestep": 1.0},
        stage_index=0,
        system_files={},
    )
    assert "GPUresident             on" not in stage0
    assert "reassignFreq            500;" in stage0
    assert "reassignTemp            $temp;" in stage0

    stage1 = mgr._customize_charmm_gui_template(
        template,
        stage_name="Equilibration 2",
        stage_params={"temperature": 303.15, "time_ns": 0.125, "timestep": 1.0},
        stage_index=1,
        system_files={},
    )
    assert "GPUresident             on" not in stage1
    assert "reassignFreq            500;" in stage1
    assert "reassignTemp            $temp;" in stage1

    prod = mgr._customize_charmm_gui_template(
        template,
        stage_name="Production",
        stage_params={"temperature": 303.15, "time_ns": 200.0, "timestep": 2.0},
        stage_index=6,
        system_files={},
    )
    assert "GPUresident             on" in prod


def test_write_job_metadata_persists_gpu_resident(tmp_path: Path) -> None:
    from gatewizard.utils.equilibration_job_metadata import (
        infer_equilibration_job_metadata,
        write_equilibration_job_metadata,
    )

    eq = tmp_path / "03_equilibration_demo"
    eq.mkdir()
    write_equilibration_job_metadata(
        eq,
        input_dir=str(tmp_path / "02_build"),
        ensemble="NPT",
        protocol={"name": "Demo", "description": "", "stages": [{"name": "Eq1"}]},
        engine="namd",
        gpu_resident=True,
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    assert meta["gpu_resident"] is True
