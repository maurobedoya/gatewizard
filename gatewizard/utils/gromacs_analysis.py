"""
GROMACS log analysis utilities for extracting timing and performance information.

Mirrors the namd_analysis interface so the two are interchangeable in the
get-equilibration-status endpoint.
"""

import re
import time
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from .logger import get_logger

logger = get_logger(__name__)


@dataclass
class GROMACSTimingInfo:
    """Timing/progress container — matches the NAMDTiming field names used by app.py."""

    steps_completed: int = 0
    total_steps: int = 0
    timestep_fs: float = 0.0  # femtoseconds per step
    ns_per_day: float = 0.0
    # Internal flags (not used by app.py but useful for status logic)
    completed: bool = False
    has_error: bool = False


@dataclass
class GROMACSStageProgress:
    """Stage progress container — matches the NAMDProgress field names used by app.py."""

    stage_name: str = ""
    status: str = "not_started"  # not_started | running | completed | error
    timing: Optional[GROMACSTimingInfo] = None
    log_file: Optional[Path] = None


def parse_gromacs_log(log_file: Path) -> GROMACSTimingInfo:
    """
    Parse a GROMACS .log file and return timing/progress information.

    While running:
      - ``nsteps`` / ``dt`` give total steps and timestep.
      - The ``Step  Time`` block rows give the current step.
      - ``Started mdrun on rank 0 <timestamp>`` compared to the current clock
        gives elapsed wall time → live ns/day estimate.

    At completion:
      - ``Performance:`` gives the official ns/day.
      - ``Finished mdrun`` marks completion.
    """
    info = GROMACSTimingInfo()

    if not log_file.exists():
        return info

    try:
        with open(log_file, "r", encoding="utf-8", errors="ignore") as fh:
            content = fh.read()

        # ── MDP parameters (echoed near the top of every GROMACS log) ───────────
        nsteps_m = re.search(r"\bnsteps\s*=\s*(\d+)", content)
        if nsteps_m:
            info.total_steps = int(nsteps_m.group(1))

        # Timestep dt is in ps → convert to fs
        dt_m = re.search(r"\bdt\s*=\s*([\d.]+)", content)
        if dt_m:
            info.timestep_fs = float(dt_m.group(1)) * 1000.0

        # ── Current step from periodic energy-output blocks ──────────────────────
        # GROMACS writes at each nstlog interval:
        #            Step           Time
        #          100000      200.00000
        step_matches = re.findall(
            r"^\s{5,15}(\d+)\s+\d+\.\d+\s*$", content, re.MULTILINE
        )
        if step_matches:
            info.steps_completed = int(step_matches[-1])

        # ── Performance (only present after stage finishes) ──────────────────────
        # "Performance:    487.654          0.049"  (ns/day  hours/ns ...)
        perf_matches = re.findall(r"^Performance:\s+([\d.]+)", content, re.MULTILINE)
        if perf_matches:
            info.ns_per_day = float(perf_matches[-1])

        # ── Live estimate: parse start timestamp, compare to now ─────────────────
        # "Started mdrun on rank 0 Thu May 28 15:40:23 2026"
        # Use this only when Performance hasn't been written yet (stage still running).
        if info.ns_per_day == 0.0 and info.steps_completed > 0 and info.timestep_fs > 0:
            start_m = re.search(
                r"Started mdrun on rank \d+ \w+ (\w+ +\d+ +\d+:\d+:\d+ \d+)",
                content,
            )
            if start_m:
                try:
                    start_dt = datetime.strptime(
                        start_m.group(1).strip(), "%b %d %H:%M:%S %Y"
                    )
                    wall_elapsed_s = time.time() - start_dt.timestamp()
                    if wall_elapsed_s > 1.0:
                        simulated_ns = info.steps_completed * info.timestep_fs * 1e-6
                        info.ns_per_day = simulated_ns / wall_elapsed_s * 86400
                except ValueError:
                    pass

        # ── Completion / error markers ───────────────────────────────────────────
        if "Finished mdrun" in content:
            info.completed = True
            if info.total_steps > 0:
                info.steps_completed = info.total_steps

        if "Fatal error:" in content or "Error in user input:" in content:
            info.has_error = True

    except Exception as exc:  # pragma: no cover
        logger.debug(f"Error parsing GROMACS log {log_file}: {exc}")

    return info


def get_equilibration_progress(
    equilibration_dir: Path,
) -> Dict[str, GROMACSStageProgress]:
    """
    Return a progress dict for all standard GROMACS equilibration stages.

    Stages without a log file are included with ``status='not_started'`` so
    the GUI can display the full stage list.

    Args:
        equilibration_dir: The working directory that was passed to the run
            script (i.e. the directory that contains ``run_equilibration.sh``
            and where GROMACS writes ``step*.log`` files).

    Returns:
        Mapping of stage-name → :class:`GROMACSStageProgress`.
    """
    # Ordered mapping: display name → expected log filename
    stage_log_map: Dict[str, str] = {
        "minimization": "step0_minimization.log",
        "equilibration_1": "step1_equilibration.log",
        "equilibration_2": "step2_equilibration.log",
        "equilibration_3": "step3_equilibration.log",
        "equilibration_4": "step4_equilibration.log",
        "equilibration_5": "step5_equilibration.log",
        "equilibration_6": "step6_equilibration.log",
        "production": "step7_production.log",
    }

    progress: Dict[str, GROMACSStageProgress] = {}

    for stage_name, log_name in stage_log_map.items():
        stage = GROMACSStageProgress(stage_name=stage_name)
        log_file = equilibration_dir / log_name

        if log_file.exists():
            stage.log_file = log_file
            timing = parse_gromacs_log(log_file)
            stage.timing = timing

            if timing.has_error:
                stage.status = "error"
            elif timing.completed:
                stage.status = "completed"
            elif timing.steps_completed > 0:
                stage.status = "running"
            else:
                # Log file exists but no step output yet — job just started
                stage.status = "running"

        progress[stage_name] = stage

    # Trim trailing not_started stages so we don't clutter the UI when the
    # user hasn't configured all 8 stages (e.g. no production stage yet).
    keys = list(progress.keys())
    while keys and progress[keys[-1]].status == "not_started":
        del progress[keys[-1]]
        keys.pop()

    return progress
