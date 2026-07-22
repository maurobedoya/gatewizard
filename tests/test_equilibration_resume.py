"""Tests for equilibration stage-level resume detection."""

from pathlib import Path

from gatewizard.utils.equilibration_resume import (
    equilibration_script_supports_resume,
    get_equilibration_resume_point,
    refresh_equilibration_run_script,
)


def test_get_equilibration_resume_point_openmm_partial(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    (eq / "run_equilibration.sh").write_text(
        'RESUME="${RESUME:-0}"\n_gw_openmm_stage_done() { true; }\n'
    )
    stem = "step1_equilibration"
    (eq / f"{stem}.inp").write_text("nstep = 125000\ndt = 0.001\n")
    (eq / f"{stem}.rst").write_text("rst")
    (eq / f"{stem}.log").write_text(
        '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
        "100.0%\t125000\t125.0\t40.0\n"
    )
    (eq / "step2_equilibration.inp").write_text("nstep = 125000\ndt = 0.001\n")

    point = get_equilibration_resume_point(eq, "openmm")

    assert point.can_resume is True
    assert point.stage_index == 1
    assert point.stage_name == "Equilibration 2"
    assert point.stage_stem == "step2_equilibration"
    assert point.completed_stages == 1


def test_get_equilibration_resume_point_openmm_stage2_interrupted(tmp_path: Path) -> None:
    """Stage 1 finished in log; kill mid stage 2 → resume Equilibration 2."""
    eq = tmp_path / "job"
    eq.mkdir()
    (eq / "run_equilibration.sh").write_text(
        'RESUME="${RESUME:-0}"\n_gw_openmm_stage_done() { true; }\n'
    )
    (eq / "equilibration_start_time.txt").write_text("2026-01-01T00:00:00+00:00")
    for stem, pct in (("step1_equilibration", "100.0"), ("step2_equilibration", "20.0")):
        (eq / f"{stem}.inp").write_text("nstep = 125000\ndt = 0.001\n")
        (eq / f"{stem}.log").write_text(
            '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
            f"{pct}%\t1000\t1.0\t40.0\n"
        )
    for stem in ("step3_equilibration", "step4_equilibration"):
        (eq / f"{stem}.inp").write_text("nstep = 1\n")

    point = get_equilibration_resume_point(eq, "openmm")

    assert point.can_resume is True
    assert point.stage_index == 1
    assert point.stage_name == "Equilibration 2"
    assert point.stage_stem == "step2_equilibration"
    assert point.completed_stages == 1


def test_get_equilibration_resume_point_interrupted_first_stage(tmp_path: Path) -> None:
    """Kill/crash mid-stage-1 should still offer Continue (restart stage 1)."""
    eq = tmp_path / "job"
    eq.mkdir()
    (eq / "run_equilibration.sh").write_text(
        'RESUME="${RESUME:-0}"\n_gw_openmm_stage_done() { true; }\n'
    )
    (eq / "equilibration_start_time.txt").write_text("2026-01-01T00:00:00+00:00")
    stem = "step1_equilibration"
    (eq / f"{stem}.inp").write_text("nstep = 125000\ndt = 0.001\n")
    (eq / f"{stem}.log").write_text(
        '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
        "35.0%\t43750\t43.75\t38.0\n"
    )
    (eq / "step2_equilibration.inp").write_text("nstep = 125000\ndt = 0.001\n")

    point = get_equilibration_resume_point(eq, "openmm")

    assert point.can_resume is True
    assert point.stage_index == 0
    assert point.stage_name == "Equilibration 1"
    assert point.stage_stem == stem
    assert point.completed_stages == 0


def test_get_equilibration_resume_point_none_completed(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    (eq / "run_equilibration.sh").write_text("#!/bin/bash\n")
    (eq / "step1_equilibration.inp").write_text("nstep = 1\n")

    point = get_equilibration_resume_point(eq, "openmm")

    assert point.can_resume is False
    assert "No completed stages" in point.reason


def test_equilibration_script_supports_resume(tmp_path: Path) -> None:
    script = tmp_path / "run_equilibration.sh"
    script.write_text("#!/bin/bash\necho hello\n")
    assert equilibration_script_supports_resume(script) is False

    script.write_text('RESUME="${RESUME:-0}"\n_gw_openmm_stage_done() { true; }\n')
    assert equilibration_script_supports_resume(script) is True


def test_refresh_equilibration_run_script_openmm(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    (eq / "step1_equilibration.inp").write_text("nstep = 1\n")
    (eq / "run_equilibration.sh").write_text(
        "#!/bin/bash\n"
        'PRMTOP="system.prmtop"\n'
        'INPCRD="system.inpcrd"\n'
        "echo legacy\n"
    )

    ok = refresh_equilibration_run_script(eq, "openmm")

    assert ok is True
    text = (eq / "run_equilibration.sh").read_text(encoding="utf-8")
    assert equilibration_script_supports_resume(eq / "run_equilibration.sh")
    assert "_gw_openmm_stage_done" in text
