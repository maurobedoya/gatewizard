"""Tests for OpenMM equilibration log progress parsing."""

from pathlib import Path

from gatewizard.utils.openmm_analysis import parse_openmm_log


def test_parse_openmm_log_uses_stage_local_progress_percent(tmp_path: Path) -> None:
    """Cumulative Step values must not inflate per-stage ns (chained -irst stages)."""
    inp = tmp_path / "step2_equilibration.inp"
    log = tmp_path / "step2_equilibration.log"
    inp.write_text("nstep = 125000\ndt = 0.001\n")
    log.write_text(
        '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
        "50.0%\t187500\t187.5\t38.0\n"
    )

    info = parse_openmm_log(log, inp)

    assert info.total_steps == 125000
    assert info.steps_completed == 62_500


def test_parse_openmm_log_completed_stage(tmp_path: Path) -> None:
    inp = tmp_path / "step1_equilibration.inp"
    log = tmp_path / "step1_equilibration.log"
    inp.write_text("nstep = 125000\ndt = 0.001\n")
    log.write_text(
        '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
        "100.0%\t125000\t125.0\t40.0\n"
    )

    info = parse_openmm_log(log, inp)

    assert info.completed is True
    assert info.steps_completed == 125_000
