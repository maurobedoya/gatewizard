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

Performance:            6.100
Finished mdrun on rank 0 Wed Jul 22 12:01:30 2026
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
    assert abs(info.wall_elapsed_seconds - 30.0) < 0.01
