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
