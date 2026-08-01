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


def test_parse_openmm_log_completed_at_99_percent(tmp_path: Path) -> None:
    """OpenMM often reports 99.2% on the final row instead of 100%."""
    inp = tmp_path / "step1_equilibration.inp"
    log = tmp_path / "step1_equilibration.log"
    inp.write_text("nstep = 125000\ndt = 0.001\n")
    log.write_text(
        '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
        "98.4%\t124000\t124.0\t82.0\n"
        "99.2%\t125000\t125.0\t81.8\n"
    )

    info = parse_openmm_log(log, inp)

    assert info.completed is True
    assert info.steps_completed == 125_000


def test_parse_openmm_log_production_100_percent(tmp_path: Path) -> None:
    """Production at 100% with cumulative Step column still completes the stage."""
    inp = tmp_path / "step7_production.inp"
    log = tmp_path / "step7_production.log"
    inp.write_text("nstep = 100000000\ndt = 0.002\n")
    log.write_text(
        '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
        "100.0%\t101375000\t201875.0\t159.0\n"
    )

    info = parse_openmm_log(log, inp)

    assert info.completed is True
    assert info.steps_completed == 100_000_000
    assert abs(info.steps_completed * info.timestep_fs * 1e-6 - 200.0) < 0.01


def test_parse_openmm_log_reads_tail_of_large_file(tmp_path: Path) -> None:
    """Progress parsing must not read multi-MB production logs line-by-line."""
    inp = tmp_path / "step7_production.inp"
    log = tmp_path / "step7_production.log"
    inp.write_text("nstep = 100000000\ndt = 0.002\n")
    header = '#"Progress (%)"\t"Step"\t"Time (ps)"\t"Speed (ns/day)"\n'
    filler = "50.0%\t50000000\t100000.0\t159.0\n" * 5000
    tail = "100.0%\t100000000\t200000.0\t159.0\n"
    log.write_text(header + filler + tail)

    info = parse_openmm_log(log, inp)

    assert info.completed is True
    assert info.steps_completed == 100_000_000
