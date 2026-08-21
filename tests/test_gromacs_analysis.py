"""Tests for GROMACS log parsing, including energy minimization stages."""

from __future__ import annotations

from pathlib import Path

from gatewizard.utils.gromacs_analysis import parse_gromacs_log


_MINIMIZATION_LOG = """\
integrator              = steep
nsteps                  = 5000
Started mdrun on rank 0 Wed Jul 22 12:00:00 2026

           Step           Time
              0      0.00000
           2500      0.00000

Steepest Descents converged to Fmax < 1000 in 2345 steps
Finished mdrun on rank 0 Wed Jul 22 12:00:42 2026
"""

_MD_LOG = """\
integrator              = md
nsteps                  = 125000
dt                      = 0.001
Started mdrun on rank 0 Wed Jul 22 12:01:00 2026

           Step           Time
              0      0.00000
           62500      62.50000
          125000     125.00000

Performance:            6.100
Finished mdrun on rank 0 Wed Jul 22 12:01:30 2026
"""

_MD_KILLED_LOG = """\
integrator              = md
nsteps                  = 50000000
dt                      = 0.002
Started mdrun on rank 0 Sat Jul 25 16:17:24 2026

           Step           Time
              0        0.00000
        5671650    11343.30000

Performance:       64.001        0.375        2.700           34.133
Finished mdrun on rank 0 Sat Jul 25 20:32:38 2026
"""


def test_parse_gromacs_minimization_log(tmp_path: Path) -> None:
    log = tmp_path / "step0_minimization.log"
    log.write_text(_MINIMIZATION_LOG, encoding="utf-8")

    info = parse_gromacs_log(log, is_minimization=True)

    assert info.is_minimization is True
    assert info.total_steps == 5000
    assert info.steps_completed == 2345
    assert info.completed is True
    assert info.ns_per_day == 0.0
    assert info.timestep_fs == 0.0
    assert info.wall_elapsed_seconds == 42.0
    assert info.converged_early is True


def test_parse_gromacs_minimization_completed_all_steps(tmp_path: Path) -> None:
    log = tmp_path / "step0_minimization.log"
    log.write_text(
        """\
integrator              = steep
nsteps                  = 5000
Started mdrun on rank 0 Wed Jul 22 12:00:00 2026
           Step           Time
           5000      0.00000
Finished mdrun on rank 0 Wed Jul 22 12:00:42 2026
""",
        encoding="utf-8",
    )

    info = parse_gromacs_log(log, is_minimization=True)

    assert info.steps_completed == 5000
    assert info.converged_early is False


def test_parse_gromacs_2026_steepest_descents_start_banner(tmp_path: Path) -> None:
    """GROMACS 2026 EM logs use 'Started Steepest Descents', not 'Started mdrun'."""
    log = tmp_path / "step0_minimization.log"
    log.write_text(
        """\
integrator              = steep
nsteps                  = 10000
Started Steepest Descents on rank 0 Sat Jul 25 14:51:01 2026
           Step           Time
           109      109.00000
Finished mdrun on rank 0 Sat Jul 25 14:53:43 2026
""",
        encoding="utf-8",
    )

    info = parse_gromacs_log(log, is_minimization=True)

    assert info.is_minimization is True
    assert info.wall_elapsed_seconds == 162.0  # 14:53:43 - 14:51:01
    assert info.completed is True


def test_parse_gromacs_md_log(tmp_path: Path) -> None:
    log = tmp_path / "step1_equilibration.log"
    log.write_text(_MD_LOG, encoding="utf-8")

    info = parse_gromacs_log(log)

    assert info.is_minimization is False
    assert info.total_steps == 125000
    assert info.steps_completed == 125000
    assert info.timestep_fs == 1.0
    assert info.ns_per_day == 6.1
    assert info.completed is True
    assert info.interrupted is False
    assert abs(info.wall_elapsed_seconds - 30.0) < 0.01


def test_parse_gromacs_md_killed_before_nsteps(tmp_path: Path) -> None:
    """Kill MD still prints Performance/Finished — must not report 100% complete."""
    log = tmp_path / "step7_production.log"
    log.write_text(_MD_KILLED_LOG, encoding="utf-8")

    info = parse_gromacs_log(log)

    assert info.completed is False
    assert info.interrupted is True
    assert info.steps_completed == 5671650
    assert info.total_steps == 50000000
    assert info.ns_per_day == 64.001


def test_parse_gromacs_large_log_finds_started_after_topology(tmp_path: Path) -> None:
    """Topology dump sits between MDP echo and Started mdrun; progress uses head+tail."""
    log = tmp_path / "step7_production.log"
    header = (
        "integrator              = md\n"
        "nsteps                  = 125000\n"
        "dt                      = 0.001\n"
    )
    started = "Started mdrun on rank 0 Wed Jul 22 12:01:00 2026\n"
    footer = (
        "           Step           Time\n"
        "          125000     125.00000\n"
        "\n"
        "Performance:            6.100\n"
        "Finished mdrun on rank 0 Wed Jul 22 12:01:30 2026\n"
    )
    with log.open("wb") as handle:
        handle.write(header.encode("utf-8"))
        handle.write(b"x" * (200 * 1024))
        handle.write(b"\n")
        handle.write(started.encode("utf-8"))
        row = b"           Step           Time\n          1000       1.00000\n"
        while handle.tell() < 6 * 1024 * 1024:
            handle.write(row)
        handle.write(footer.encode("utf-8"))

    info = parse_gromacs_log(log)

    assert info.total_steps == 125000
    assert info.steps_completed == 125000
    assert info.ns_per_day == 6.1
    assert info.completed is True
    assert abs(info.wall_elapsed_seconds - 30.0) < 0.01
