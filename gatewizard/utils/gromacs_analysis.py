"""
GROMACS log analysis utilities for extracting timing and performance information.

Mirrors the namd_analysis interface so the two are interchangeable in the
get-equilibration-status endpoint.
"""

import math
import re
import time
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .logger import get_logger

logger = get_logger(__name__)


@dataclass
class GROMACSTimingInfo:
    """Timing/progress container — matches the NAMDTiming field names used by app.py."""

    steps_completed: int = 0
    total_steps: int = 0
    timestep_fs: float = 0.0  # femtoseconds per step
    ns_per_day: float = 0.0
    is_minimization: bool = False
    wall_elapsed_seconds: float = 0.0
    converged_early: bool = False
    # Internal flags (not used by app.py but useful for status logic)
    completed: bool = False
    has_error: bool = False
    # True when the log ends with ``Finished mdrun`` but MD did not reach nsteps
    # (typical after Kill MD / SIGTERM — GROMACS still prints Performance/Finished).
    interrupted: bool = False


@dataclass
class GROMACSStageProgress:
    """Stage progress container — matches the NAMDProgress field names used by app.py."""

    stage_name: str = ""
    status: str = "not_started"  # not_started | running | completed | error
    timing: Optional[GROMACSTimingInfo] = None
    log_file: Optional[Path] = None


def _parse_mdrun_timestamp(content: str, marker: str) -> Optional[datetime]:
    """Parse ``Started …`` / ``Finished mdrun`` timestamps from a GROMACS log."""
    pattern = rf"{marker} on rank \d+ \w+ (\w+ +\d+ +\d+:\d+:\d+ \d+)"
    match = re.search(pattern, content)
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1).strip(), "%b %d %H:%M:%S %Y")
    except ValueError:
        return None


def _parse_started_timestamp(content: str) -> Optional[datetime]:
    """Parse the wall-clock start time from a GROMACS mdrun log.

    Dynamics logs use ``Started mdrun …``. Energy minimisation in GROMACS
    2026+ uses ``Started Steepest Descents …`` / ``Started Conjugate Gradients …``.
    """
    for marker in (
        "Started mdrun",
        "Started Steepest Descents",
        "Started Conjugate Gradients",
    ):
        dt = _parse_mdrun_timestamp(content, marker)
        if dt is not None:
            return dt
    # Fallback for future integrator-specific banners.
    match = re.search(
        r"Started .+ on rank \d+ \w+ (\w+ +\d+ +\d+:\d+:\d+ \d+)",
        content,
    )
    if not match:
        return None
    try:
        return datetime.strptime(match.group(1).strip(), "%b %d %H:%M:%S %Y")
    except ValueError:
        return None


def _wall_elapsed_seconds(content: str, *, running: bool) -> float:
    start_dt = _parse_started_timestamp(content)
    if not start_dt:
        return 0.0
    end_dt = _parse_mdrun_timestamp(content, "Finished mdrun")
    if end_dt:
        return max(0.0, end_dt.timestamp() - start_dt.timestamp())
    if running:
        return max(0.0, time.time() - start_dt.timestamp())
    return 0.0


def parse_gromacs_log(
    log_file: Path, *, is_minimization: bool = False
) -> GROMACSTimingInfo:
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

        integrator_m = re.search(r"integrator\s*=\s*(\S+)", content)
        if integrator_m:
            integrator = integrator_m.group(1).strip().lower()
            if integrator in {"steep", "cg"}:
                is_minimization = True
        info.is_minimization = is_minimization

        # ── MDP parameters (echoed near the top of every GROMACS log) ───────────
        nsteps_m = re.search(r"\bnsteps\s*=\s*(\d+)", content)
        if nsteps_m:
            info.total_steps = int(nsteps_m.group(1))

        if is_minimization:
            requested_steps = info.total_steps
            step_matches = re.findall(
                r"^\s{5,15}(\d+)\s+\d+\.\d+\s*$", content, re.MULTILINE
            )
            if step_matches:
                info.steps_completed = int(step_matches[-1])

            converged_m = re.search(
                r"converged.*?in\s+(\d+)\s+steps", content, re.IGNORECASE
            )
            if converged_m:
                info.steps_completed = int(converged_m.group(1))

            if "Finished mdrun" in content:
                info.completed = True
                if info.total_steps > 0 and info.steps_completed == 0:
                    info.steps_completed = info.total_steps

            if (
                info.completed
                and requested_steps > 0
                and 0 < info.steps_completed < requested_steps
                and converged_m
            ):
                info.converged_early = True
                info.total_steps = requested_steps

            info.wall_elapsed_seconds = _wall_elapsed_seconds(
                content, running=not info.completed
            )
        else:
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
            if (
                info.ns_per_day == 0.0
                and info.steps_completed > 0
                and info.timestep_fs > 0
            ):
                start_dt = _parse_started_timestamp(content)
                if start_dt:
                    wall_elapsed_s = max(0.0, time.time() - start_dt.timestamp())
                    if wall_elapsed_s > 1.0:
                        simulated_ns = info.steps_completed * info.timestep_fs * 1e-6
                        info.ns_per_day = simulated_ns / wall_elapsed_s * 86400
                        info.wall_elapsed_seconds = wall_elapsed_s

            # ── Completion / error markers ───────────────────────────────────────────
            # GROMACS often still writes Performance + "Finished mdrun" after Kill MD.
            # Only treat the stage as complete when the last logged step reached nsteps.
            if "Finished mdrun" in content:
                if info.total_steps > 0 and info.steps_completed >= info.total_steps:
                    info.completed = True
                    info.steps_completed = info.total_steps
                elif info.total_steps > 0 and info.steps_completed > 0:
                    info.completed = False
                    info.interrupted = True
                elif info.total_steps > 0:
                    # Finished banner with no step rows — treat as incomplete.
                    info.completed = False
                    info.interrupted = True
                else:
                    info.completed = True
                if info.wall_elapsed_seconds == 0.0:
                    info.wall_elapsed_seconds = _wall_elapsed_seconds(
                        content, running=False
                    )

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
            timing = parse_gromacs_log(
                log_file, is_minimization=(stage_name == "minimization")
            )
            stage.timing = timing

            if timing.has_error or timing.interrupted:
                stage.status = "error"
            elif timing.completed:
                stage.status = "completed"
            elif timing.steps_completed > 0:
                stage.status = "running"
            else:
                # Log file exists but no step output yet — job just started
                stage.status = "running"

        progress[stage_name] = stage

    return progress


# ---------------------------------------------------------------------------
# GROMACS log energy analysis
# ---------------------------------------------------------------------------

# Labels that represent energy quantities (kJ/mol in GROMACS log)
_GROMACS_ENERGY_LABELS: frozenset = frozenset(
    {
        "Bond",
        "Angle",
        "Proper Dih.",
        "Improper Dih.",
        "LJ-14",
        "Coulomb-14",
        "LJ (SR)",
        "LJ (LR)",
        "Disper. corr.",
        "Coulomb (SR)",
        "Coulomb (LR)",
        "Coul. recip.",
        "Potential",
        "Kinetic En.",
        "Total Energy",
        "Conserved En.",
        "G96Bond",
        "G96Angle",
    }
)
_GROMACS_TEMP_LABELS: frozenset = frozenset({"Temperature"})
_GROMACS_PRESSURE_LABELS: frozenset = frozenset({"Pressure (bar)", "Pres. DC (bar)"})

# Field width used in GROMACS log energies blocks
_GROMACS_FIELD_WIDTH = 15


def _gromacs_label_unit_type(label: str) -> str:
    """Return 'energy' | 'temperature' | 'pressure' | 'other' for a GROMACS label."""
    if label in _GROMACS_ENERGY_LABELS:
        return "energy"
    if label in _GROMACS_TEMP_LABELS:
        return "temperature"
    if label in _GROMACS_PRESSURE_LABELS:
        return "pressure"
    return "other"


def _parse_gromacs_energies_lines(block_lines: List[str]) -> Dict[str, float]:
    """Parse the alternating label/value line pairs in a GROMACS Energies block.

    GROMACS writes 15-char fixed-width fields for both labels (right-aligned
    text) and values (scientific notation floats).  Lines come in pairs:
    one header row followed by one data row, repeating.
    """
    result: Dict[str, float] = {}
    FW = _GROMACS_FIELD_WIDTH
    i = 0
    while i + 1 < len(block_lines):
        label_line = block_lines[i].rstrip("\n")
        value_line = block_lines[i + 1].rstrip("\n")

        # Skip blank lines
        if not label_line.strip() or not value_line.strip():
            i += 1
            continue

        # Pad both lines to the same length so field slicing is safe
        max_len = max(len(label_line), len(value_line))
        label_line = label_line.ljust(max_len)
        value_line = value_line.ljust(max_len)

        n_fields = (max_len + FW - 1) // FW

        # Check that value_line is all-numeric (reject label/label pairs)
        value_fields_ok = True
        parsed_pairs: List[tuple] = []
        for j in range(n_fields):
            s = j * FW
            e = s + FW
            lbl = label_line[s:e].strip()
            val_str = value_line[s:e].strip()
            if lbl and val_str:
                try:
                    parsed_pairs.append((lbl, float(val_str)))
                except ValueError:
                    # value slot contains text → this is a label/label pair; skip both
                    value_fields_ok = False
                    break

        if value_fields_ok:
            for lbl, val in parsed_pairs:
                result[lbl] = val
            i += 2
        else:
            i += 1

    return result


class GROMACSLogEnergyAnalyzer:
    """Parse energy data from GROMACS .log text files.

    GROMACS writes periodic energy blocks in the log file:

        Step           Time
      100000      200.00000

      Energies (kJ/mol)
              Bond          Angle    Proper Dih.  ...
       1.66553e+02    2.66869e+02    3.64803e+02  ...

    This class collects all such blocks and provides the same interface as
    :class:`~gatewizard.utils.namd_analysis.EnergyAnalyzer` and
    :class:`~gatewizard.utils.openmm_analysis.OpenMMLogAnalyzer`.

    Args:
        log_files: Path(s) to GROMACS .log file(s).
        file_times: Optional ``{filename: duration_ns}`` override for the time
            axis (same semantics as in the NAMD/OpenMM analyzers).
    """

    def __init__(
        self,
        log_files: Union[Path, str, List[Union[Path, str]]],
        file_times: Optional[Dict[str, float]] = None,
    ) -> None:
        if isinstance(log_files, (str, Path)):
            self.log_files = [Path(log_files)]
        else:
            self.log_files = [Path(f) for f in log_files]
        self.file_times = file_times or {}
        self._file_ranges: Dict[str, tuple] = {}
        self.data = self._parse_all()

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def _parse_file(self, log_file: Path) -> Dict[str, List[float]]:
        """Parse a single GROMACS log file; return per-label time-series."""
        data: Dict[str, List[float]] = {"time_ns": []}
        try:
            with open(log_file, "r", encoding="utf-8", errors="ignore") as fh:
                lines = fh.readlines()
        except Exception as exc:
            logger.warning(f"Cannot read GROMACS log {log_file}: {exc}")
            return data

        i = 0
        current_time_ns: Optional[float] = None

        while i < len(lines):
            line = lines[i]

            # ── "Step    Time" header ────────────────────────────────
            if re.match(r"^\s+Step\s+Time\s*$", line):
                i += 1
                if i < len(lines):
                    m = re.match(r"^\s+(\d+)\s+([\d.]+)\s*$", lines[i])
                    if m:
                        current_time_ns = float(m.group(2)) / 1000.0  # ps → ns
                i += 1
                continue

            # ── "Energies (kJ/mol)" block ───────────────────────────
            if "Energies (kJ/mol)" in line and current_time_ns is not None:
                i += 1
                block_lines: List[str] = []
                while i < len(lines) and lines[i].strip():
                    block_lines.append(lines[i])
                    i += 1

                energies = _parse_gromacs_energies_lines(block_lines)
                if not energies:
                    continue

                data["time_ns"].append(current_time_ns)
                n_prev = len(data["time_ns"]) - 1

                # Merge new labels, backfill with NaN if not seen before
                for lbl in list(data.keys()):
                    if lbl == "time_ns":
                        continue
                    if lbl in energies:
                        data[lbl].append(energies[lbl])
                    else:
                        data[lbl].append(math.nan)

                for lbl, val in energies.items():
                    if lbl not in data:
                        data[lbl] = [math.nan] * n_prev + [val]

                continue

            i += 1

        return data

    def _parse_all(self) -> Dict[str, List[float]]:
        """Parse all log files, concatenating time-series with proper offsets."""
        combined: Dict[str, List[float]] = {"time_ns": []}
        cumulative_ns = 0.0

        for log_file in self.log_files:
            if not log_file.exists():
                logger.warning(f"GROMACS log not found: {log_file}")
                continue

            start_idx = len(combined["time_ns"])
            file_data = self._parse_file(log_file)
            n_points = len(file_data.get("time_ns", []))

            if n_points == 0:
                continue

            # Apply file_times override
            fname = log_file.name
            if fname in self.file_times and self.file_times[fname] > 0:
                import numpy as np

                duration = self.file_times[fname]
                file_data["time_ns"] = np.linspace(
                    0, duration, n_points, endpoint=False
                ).tolist()

            # Add cumulative offset
            file_data["time_ns"] = [t + cumulative_ns for t in file_data["time_ns"]]
            last_t = file_data["time_ns"][-1]
            dt = (
                (file_data["time_ns"][-1] - file_data["time_ns"][-2])
                if n_points > 1
                else 0.0
            )
            cumulative_ns = last_t + dt

            # Merge into combined
            for key, vals in file_data.items():
                if key not in combined:
                    n_so_far = len(combined["time_ns"]) - n_points
                    combined[key] = [math.nan] * n_so_far
                combined[key].extend(vals)

            # Backfill missing keys in this file
            total_now = len(combined["time_ns"])
            for key in list(combined.keys()):
                if len(combined[key]) < total_now:
                    combined[key].extend([math.nan] * (total_now - len(combined[key])))

            end_idx = len(combined["time_ns"])
            self._file_ranges[str(log_file)] = (start_idx, end_idx)

        return combined

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def _calculate_time_array(self) -> "np.ndarray":
        """Return time in nanoseconds as a numpy array."""
        import numpy as np

        return np.array(self.data.get("time_ns", []), dtype=float)

    def get_available_properties(self) -> List[str]:
        """Return labels that have at least one non-NaN value, sorted."""
        props = []
        for key, vals in self.data.items():
            if key == "time_ns":
                continue
            if vals and any(not math.isnan(v) for v in vals):
                props.append(key)
        return sorted(props)

    def get_statistics(self) -> Dict[str, Dict[str, float]]:
        """Return mean/std/min/max/initial/final for every non-empty property."""
        import numpy as np

        stats: Dict[str, Dict[str, float]] = {}
        for key, values in self.data.items():
            if key == "time_ns" or not values:
                continue
            arr = np.array(values, dtype=float)
            valid = arr[~np.isnan(arr)]
            if len(valid) == 0:
                continue
            stats[key] = {
                "mean": float(np.mean(valid)),
                "std": float(np.std(valid)),
                "min": float(np.min(valid)),
                "max": float(np.max(valid)),
                "initial": float(valid[0]),
                "final": float(valid[-1]),
            }
        return stats


def _gromacs_convert_units(
    arr: "np.ndarray",
    label: str,
    energy_units: str,
    pressure_units: str,
    temperature_units: str,
) -> tuple:
    """Convert a GROMACS property array from native units.

    GROMACS log files always use kJ/mol (energy), K (temperature), bar (pressure).
    Returns ``(converted_array, unit_label)``.
    """
    utype = _gromacs_label_unit_type(label)

    if utype == "energy":
        if energy_units == "kcal/mol":
            return arr / 4.184, "kcal/mol"
        return arr.copy(), "kJ/mol"

    if utype == "temperature":
        if temperature_units == "°C":
            return arr - 273.15, "°C"
        if temperature_units == "°F":
            return (arr - 273.15) * 9 / 5 + 32, "°F"
        return arr.copy(), "K"

    if utype == "pressure":
        # GROMACS writes bar; convert if needed
        _bar_factors = {"atm": 1 / 1.01325, "bar": 1.0, "kPa": 100.0, "MPa": 0.1}
        factor = _bar_factors.get(pressure_units, 1.0)
        unit_str = pressure_units if pressure_units in _bar_factors else "bar"
        return arr * factor, unit_str

    return arr.copy(), ""


def list_gromacs_energy_properties(
    log_files: List[Union[str, Path]],
    file_times: Optional[Dict[str, float]] = None,
) -> List[str]:
    """Return available energy property labels from GROMACS log file(s).

    Mirrors :func:`~gatewizard.utils.namd_analysis.list_namd_energy_properties`.
    """
    logs = [Path(f) for f in log_files]
    analyzer = GROMACSLogEnergyAnalyzer(logs, file_times=file_times)
    return analyzer.get_available_properties()


def run_gromacs_energetic_analysis(
    log_files: List[Union[str, Path]],
    properties: Optional[List[str]] = None,
    file_times: Optional[Dict[str, float]] = None,
    time_units: str = "ns",
    energy_units: str = "kcal/mol",
    pressure_units: str = "atm",
    temperature_units: str = "K",
    volume_units: str = "Å³",
) -> Dict[str, Any]:
    """Run GROMACS log energetic analysis and return JSON-serializable arrays.

    Returns the same dict structure as
    :func:`~gatewizard.utils.namd_analysis.run_energetic_analysis`:
    ``{x, x_label, series: [{name, key, unit, y}], statistics}``.
    """
    import numpy as np

    logs = [Path(f) for f in log_files]
    analyzer = GROMACSLogEnergyAnalyzer(logs, file_times=file_times)

    # Time axis
    x = analyzer._calculate_time_array()
    x_label = "Time (ns)"
    if time_units == "ps":
        x = x * 1000.0
        x_label = "Time (ps)"
    elif time_units in {"us", "µs"}:
        x = x / 1000.0
        x_label = "Time (µs)"

    selected = (
        properties if properties is not None else analyzer.get_available_properties()
    )

    series = []
    statistics: Dict[str, Any] = {}

    for label in selected:
        raw = analyzer.data.get(label, [])
        if not raw:
            continue
        arr = np.array(raw, dtype=float)
        if not (~np.isnan(arr)).any():
            continue

        converted, unit_label = _gromacs_convert_units(
            arr, label, energy_units, pressure_units, temperature_units
        )

        series.append(
            {
                "name": label,
                "key": label,
                "unit": unit_label,
                "y": converted.tolist(),
            }
        )

        valid_c = converted[~np.isnan(converted)]
        if len(valid_c) > 0:
            statistics[label] = {
                "mean": float(np.nanmean(valid_c)),
                "std": float(np.nanstd(valid_c)),
                "min": float(np.nanmin(valid_c)),
                "max": float(np.nanmax(valid_c)),
                "initial": float(valid_c[0]),
                "final": float(valid_c[-1]),
            }

    return {
        "x": x.tolist(),
        "x_label": x_label,
        "series": series,
        "statistics": statistics,
    }
