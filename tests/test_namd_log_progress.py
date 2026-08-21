"""Regression tests for NAMD equilibration log progress / performance parsing."""

from pathlib import Path

from gatewizard.utils.namd_analysis import parse_namd_log

EXAMPLES = Path(__file__).parent / "analysis_examples" / "equilibration_folder"


def test_parse_namd_production_example_has_performance():
    timing = parse_namd_log(EXAMPLES / "step7_production.log")
    assert timing.steps_completed > 0
    assert timing.ns_per_day > 0
    assert timing.wall_elapsed_seconds > 0


def test_parse_namd_production_without_final_timing_uses_wallclock(tmp_path: Path):
    """Production often ends after the last TIMING print; WallClock must still yield ns/day.

    This matched a GUI bug: completed production showed ``0.0 ns/day`` and no elapsed time.
    """
    src = (EXAMPLES / "step7_production.log").read_text(encoding="utf-8", errors="ignore")
    # Drop every TIMING line whose step equals the final OUTPUT write step
    lines = []
    for line in src.splitlines(keepends=True):
        if line.startswith("TIMING:") and " 360000 " in line:
            continue
        lines.append(line)
    log = tmp_path / "step7_production.log"
    log.write_text("".join(lines), encoding="utf-8")

    timing = parse_namd_log(log)
    assert timing.steps_completed == 50000
    assert timing.simulated_time_ns == 0.1
    assert timing.wall_elapsed_seconds > 100  # WallClock ~172 s in example
    assert timing.ns_per_day > 1.0
    # ~0.1 ns in ~172 s ≈ 50 ns/day
    assert 20.0 < timing.ns_per_day < 80.0


def test_parse_namd_energy_only_with_wallclock(tmp_path: Path):
    log = tmp_path / "step7_production.log"
    log.write_text(
        "\n".join(
            [
                "Info: TIMESTEP               2",
                "Info: FIRST TIMESTEP         0",
                "TCL: Running for 1000 steps",
                "ENERGY:       500 ...",
                "ENERGY:      1000 ...",
                "WRITING COORDINATES TO OUTPUT FILE AT STEP 1000",
                "WallClock: 100.0  CPUTime: 99.0  Memory: 0.000000 MB",
                "",
            ]
        ),
        encoding="utf-8",
    )
    timing = parse_namd_log(log)
    assert timing.steps_completed == 1000
    assert timing.wall_elapsed_seconds == 100.0
    # 1000 * 2 fs = 0.002 ns in 100 s → 1.728 ns/day
    assert abs(timing.ns_per_day - 1.728) < 1e-6


def test_namd_midrun_restart_write_is_running_not_complete(tmp_path: Path):
    """Periodic NAMD restart/DCD writes must not mark the stage complete or failed."""
    from gatewizard.utils.namd_analysis import get_equilibration_progress

    log = tmp_path / "step2_equilibration.log"
    log.write_text(
        "\n".join(
            [
                "Info: TIMESTEP               2",
                "TCL: Running for 500000 steps",
                "ENERGY:    240000  ...",
                "WRITING EXTENDED SYSTEM TO RESTART FILE AT STEP 240000",
                "WRITING COORDINATES TO DCD FILE step2_equilibration.dcd AT STEP 240000",
                "WRITING COORDINATES TO RESTART FILE AT STEP 240000",
                "FINISHED WRITING RESTART COORDINATES",
                "",
            ]
        ),
        encoding="utf-8",
    )
    progress = get_equilibration_progress(tmp_path)
    assert progress["equilibration_2"].status == "running"


def test_namd_end_of_program_marks_stage_complete(tmp_path: Path):
    from gatewizard.utils.namd_analysis import get_equilibration_progress

    log = tmp_path / "step2_equilibration.log"
    log.write_text(
        "\n".join(
            [
                "Info: TIMESTEP               2",
                "TCL: Running for 1000 steps",
                "ENERGY:      1000  ...",
                "WRITING COORDINATES TO OUTPUT FILE AT STEP 1000",
                "End of program",
                "",
            ]
        ),
        encoding="utf-8",
    )
    progress = get_equilibration_progress(tmp_path)
    assert progress["equilibration_2"].status == "completed"


def test_parse_namd_large_log_reads_head_and_tail(tmp_path: Path):
    """Multi-MB ENERGY dumps must not hide header TIMESTEP or tail TIMING/ns/day."""
    from gatewizard.utils.namd_analysis import get_equilibration_progress

    header = "\n".join(
        [
            "Info: TIMESTEP               2",
            "Info: FIRST TIMESTEP         0",
            "TCL: Running for 500000 steps",
            "Info: 100000 ATOMS",
            "",
        ]
    )
    footer = "\n".join(
        [
            "TIMING: 360000  CPU: 169.927, 0.00336/step  Wall: 170.557, 0.00332/step",
            "ENERGY:      360000 ...",
            "WRITING COORDINATES TO OUTPUT FILE AT STEP 360000",
            "WallClock: 172.0  CPUTime: 170.0  Memory: 0.000000 MB",
            "End of program",
            "",
        ]
    )
    junk = ("ENERGY:      12345  " + ("x" * 80) + "\n").encode("utf-8")
    log = tmp_path / "step2_equilibration.log"
    with log.open("wb") as handle:
        handle.write(header.encode("utf-8"))
        written = handle.tell()
        target = 6 * 1024 * 1024
        while written < target:
            handle.write(junk)
            written += len(junk)
        handle.write(b"\n")
        handle.write(footer.encode("utf-8"))

    timing = parse_namd_log(log)
    assert timing.timestep_fs == 2.0
    assert timing.total_steps == 500000
    assert timing.steps_completed == 360000
    assert timing.wall_elapsed_seconds == 172.0
    assert timing.ns_per_day > 0

    progress = get_equilibration_progress(tmp_path)
    assert progress["equilibration_2"].status == "completed"
    assert progress["equilibration_2"].timing is not None
    assert progress["equilibration_2"].timing.ns_per_day > 0
