"""
NAMD log analysis utilities for extracting timing, performance, and ENERGY data.

This module parses NAMD log files for simulation progress, performance metrics,
and ENERGY series. Engine-agnostic trajectory / structural analysis lives in
``gatewizard.utils.trajectory_analysis``; bilayer analysis lives in
``gatewizard.utils.lipid_bilayer_analysis``.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass
from .energy_stride import lookup_file_map as _lookup_file_map
from .logger import get_logger
from .log_io import read_text_head_tail

logger = get_logger(__name__)


@dataclass
class NAMDTiming:
    """Container for NAMD timing information."""

    steps_completed: int = 0
    total_steps: int = 0
    simulated_time_ns: float = 0.0
    real_time_hours: float = 0.0
    ns_per_day: float = 0.0
    sec_per_step: float = 0.0
    processors: int = 0
    gpus: int = 0
    atoms: int = 0
    timestep_fs: float = 0.0
    first_timestep: int = 0
    hostname: str = ""
    wall_elapsed_seconds: float = 0.0


@dataclass
class NAMDProgress:
    """Container for NAMD progress information."""

    stage_name: str = ""
    status: str = "not_started"  # not_started, running, completed, error
    progress_percent: float = 0.0
    timing: Optional[NAMDTiming] = None
    log_file: Optional[Path] = None
    last_updated: float = 0.0


def _parse_namd_timing_wall_samples(content: str) -> list[tuple[int, float]]:
    """Return ``(step, wall_seconds)`` samples from NAMD TIMING lines.

    Modern NAMD 3 logs look like::

        TIMING: 360000  CPU: 169.927, 0.00336/step  Wall: 170.557, 0.00332/step, ...

    Older column-style TIMING lines are also accepted.
    """
    samples: list[tuple[int, float]] = []

    # Preferred: explicit Wall: field (NAMD 2.12+ / 3.x)
    for step_s, wall_s in re.findall(
        r"^TIMING:\s+(\d+)\s+.*?Wall:\s*([\d.eE+-]+)",
        content,
        re.MULTILINE,
    ):
        try:
            samples.append((int(step_s), float(wall_s)))
        except ValueError:
            continue
    if samples:
        return samples

    # Legacy whitespace columns: TIMING: step ... wall_time
    for step_s, wall_s in re.findall(
        r"^TIMING:\s+(\d+)\s+[\d.\-+eE]+\s+[\d.\-+eE]+\s+[\d.\-+eE]+\s+"
        r"[\d.\-+eE]+\s+[\d.\-+eE]+\s+[\d.\-+eE]+\s+([\d.\-+eE]+)",
        content,
        re.MULTILINE,
    ):
        try:
            samples.append((int(step_s), float(wall_s)))
        except ValueError:
            continue
    return samples


def _namd_wallclock_seconds(content: str) -> float:
    """Final ``WallClock:`` summary line, if present."""
    matches = re.findall(r"^WallClock:\s*([\d.eE+-]+)", content, re.MULTILINE)
    if not matches:
        return 0.0
    try:
        return float(matches[-1])
    except ValueError:
        return 0.0


def _namd_benchmark_ns_per_day(content: str) -> float:
    """Invert Benchmark ``days/ns`` to ns/day (average of samples)."""
    days_per_ns = re.findall(
        r"Benchmark time:.*?([\d.eE+-]+)\s+days/ns",
        content,
    )
    vals = []
    for d in days_per_ns:
        try:
            days = float(d)
        except ValueError:
            continue
        if days > 0:
            vals.append(1.0 / days)
    if not vals:
        return 0.0
    return sum(vals) / len(vals)


def _read_namd_log_text(log_file_path: Path) -> str:
    """Read a NAMD log; sample head+tail when the file is large."""
    return read_text_head_tail(log_file_path)


def parse_namd_log(log_file_path: Path) -> NAMDTiming:
    """
    Parse a NAMD log file to extract timing and performance information.

    Based on the namd_timing script functionality.
    """
    timing = NAMDTiming()

    if not log_file_path.exists():
        logger.debug(f"Log file does not exist: {log_file_path}")
        return timing

    try:
        content = _read_namd_log_text(log_file_path)

        logger.debug(
            f"Parsing log file {log_file_path.name}, size: {len(content)} chars"
        )

        # Extract basic system information
        proc_match = re.search(r"Running on (\d+) processors", content)
        if proc_match:
            timing.processors = int(proc_match.group(1))

        atoms_match = re.search(r"Info: (\d+) ATOMS", content)
        if atoms_match:
            timing.atoms = int(atoms_match.group(1))

        timestep_match = re.search(r"Info: TIMESTEP\s+(\d+(?:\.\d+)?)", content)
        if timestep_match:
            timing.timestep_fs = float(timestep_match.group(1))

        first_ts_match = re.search(r"FIRST TIMESTEP\s+(\d+)", content)
        if first_ts_match:
            timing.first_timestep = int(first_ts_match.group(1))

        # Extract hostname
        host_match = re.search(r"Info: \d+ NAMD.*?(\S+)", content)
        if host_match:
            timing.hostname = host_match.group(1)

        # Count CUDA devices
        cuda_matches = re.findall(r"CUDA device \d+", content)
        timing.gpus = len(cuda_matches)

        timing_samples = _parse_namd_timing_wall_samples(content)
        wall_clock_s = _namd_wallclock_seconds(content)

        # ENERGY lines give step progress when TIMING is sparse
        energy_steps = [
            int(s) for s in re.findall(r"^ENERGY:\s+(\d+)", content, re.MULTILINE)
        ]

        # Final coordinate write marks the true last step (may be after last TIMING)
        final_step_matches = [
            int(s)
            for s in re.findall(r"WRITING.*?TO OUTPUT FILE AT STEP (\d+)", content)
        ]

        last_step = 0
        last_wall_s = 0.0
        if timing_samples:
            last_step, last_wall_s = timing_samples[-1]
            logger.debug(
                f" Found {len(timing_samples)} TIMING lines; last step={last_step} wall={last_wall_s}"
            )
        elif energy_steps:
            last_step = energy_steps[-1]
            logger.debug(f" Found {len(energy_steps)} ENERGY lines (no TIMING wall)")

        if final_step_matches:
            final_step = final_step_matches[-1]
            if final_step > last_step:
                logger.debug(
                    f" Using final output step {final_step} (last TIMING/ENERGY was {last_step})"
                )
                last_step = final_step
            # Prefer WallClock for completed runs; else keep last TIMING wall
            if wall_clock_s > 0:
                last_wall_s = wall_clock_s
            elif last_wall_s <= 0 and timing_samples:
                last_wall_s = timing_samples[-1][1]
        elif wall_clock_s > last_wall_s:
            last_wall_s = wall_clock_s

        if last_step > 0:
            timing.steps_completed = max(0, last_step - timing.first_timestep)
            logger.debug(
                f" Last step: {last_step}, first: {timing.first_timestep}, "
                f"completed: {timing.steps_completed}"
            )

        if last_wall_s > 0:
            timing.wall_elapsed_seconds = last_wall_s
            timing.real_time_hours = last_wall_s / 3600.0
            if timing.steps_completed > 0:
                timing.sec_per_step = last_wall_s / timing.steps_completed

        if timing.timestep_fs > 0 and timing.steps_completed > 0:
            timing.simulated_time_ns = (
                timing.steps_completed * timing.timestep_fs
            ) / 1_000_000.0

            if timing.wall_elapsed_seconds > 0:
                timing.ns_per_day = timing.simulated_time_ns / (
                    timing.wall_elapsed_seconds / 86400.0
                )
            else:
                # Early in a run before any TIMING/WallClock: Benchmark estimate
                bench = _namd_benchmark_ns_per_day(content)
                if bench > 0:
                    timing.ns_per_day = bench

        # Try to get total steps from the configuration
        # Primary pattern: TCL commands (NAMD 3.0 format)
        tcl_run_match = re.search(
            r"TCL: Running for (\d+) steps", content, re.IGNORECASE
        )
        tcl_minimize_match = re.search(
            r"TCL: Minimizing for (\d+) steps", content, re.IGNORECASE
        )

        # Fallback patterns for older NAMD versions
        run_match = re.search(r"run\s+(\d+)", content)
        minimize_match = re.search(r"minimize\s+(\d+)", content)
        numsteps_match = re.search(r"numsteps\s+(\d+)", content)

        # Initialize with 0
        run_steps = 0
        minimize_steps = 0

        # Get run steps
        if tcl_run_match:
            run_steps = int(tcl_run_match.group(1))
            logger.debug(f" Found 'TCL: Running' pattern: {run_steps} steps")
        elif run_match:
            run_steps = int(run_match.group(1))
            logger.debug(f" Found 'run' pattern: {run_steps} steps")
        elif numsteps_match:
            run_steps = int(numsteps_match.group(1))
            logger.debug(f" Found 'numsteps' pattern: {run_steps} steps")

        # Get minimize steps
        if tcl_minimize_match:
            minimize_steps = int(tcl_minimize_match.group(1))
            logger.debug(f" Found 'TCL: Minimizing' pattern: {minimize_steps} steps")
        elif minimize_match:
            minimize_steps = int(minimize_match.group(1))
            logger.debug(f" Found 'minimize' pattern: {minimize_steps} steps")

        # Total steps is the sum of minimization and run steps
        if run_steps > 0 or minimize_steps > 0:
            timing.total_steps = run_steps + minimize_steps
            logger.debug(
                f" Total steps calculation: {minimize_steps} (minimize) + {run_steps} (run) = {timing.total_steps}"
            )

        # If we only found minimize steps but no run steps yet (simulation in progress),
        # try to read the input file to get the expected total steps
        if minimize_steps > 0 and run_steps == 0:
            inp_file_steps = _get_expected_steps_from_inp_file(log_file_path)
            if inp_file_steps > minimize_steps:
                timing.total_steps = inp_file_steps
                logger.debug(f" Using expected steps from input file: {inp_file_steps}")

        # Also try to find steps from NAMD output messages if we couldn't find run/minimize commands
        if timing.total_steps == 0:
            # Look for patterns like "Info: STEPS <number>"
            steps_info_match = re.search(
                r"Info:.*?STEPS?\s+(\d+)", content, re.IGNORECASE
            )
            if steps_info_match:
                timing.total_steps = int(steps_info_match.group(1))

            # Look for patterns in lines like "ETITLE:" or similar
            etitle_match = re.search(r"ETITLE:.*?(\d+)", content)
            if etitle_match and not timing.total_steps:
                # This might be less reliable, use as last resort
                pass

    except Exception as e:
        # Log error but return partial results
        pass

    return timing


def _get_expected_steps_from_inp_file(log_file_path: Path) -> int:
    """
    Try to read the corresponding input file to get expected total steps.

    Args:
        log_file_path: Path to the log file (e.g., eq1_equilibration.log)

    Returns:
        Expected total steps from input file, or 0 if not found
    """
    try:
        # Try to find corresponding .inp file
        inp_file = log_file_path.with_suffix(".inp")
        if not inp_file.exists():
            # Try common naming patterns
            base_name = log_file_path.stem.replace("_equilibration", "").replace(
                ".log", ""
            )
            possible_inp_files = [
                log_file_path.parent / f"{base_name}_equilibration.inp",
                log_file_path.parent / f"{base_name}.inp",
                log_file_path.parent
                / f"step6.{base_name.replace('eq', '')}_equilibration.inp",
            ]

            for possible_inp in possible_inp_files:
                if possible_inp.exists():
                    inp_file = possible_inp
                    break
            else:
                return 0

        with open(inp_file, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()

        # Look for minimize and run commands
        minimize_match = re.search(r"minimize\s+(\d+)", content)
        run_match = re.search(r"run\s+(\d+)", content)

        minimize_steps = int(minimize_match.group(1)) if minimize_match else 0
        run_steps = int(run_match.group(1)) if run_match else 0

        total_expected = minimize_steps + run_steps
        logger.debug(
            f" Expected steps from {inp_file.name}: {minimize_steps} (minimize) + {run_steps} (run) = {total_expected}"
        )

        return total_expected

    except Exception as e:
        logger.debug(f" Error reading input file for {log_file_path}: {e}")
        return 0


def get_equilibration_progress(equilibration_dir: Path) -> Dict[str, NAMDProgress]:
    """
    Get progress information for all equilibration stages.

    Looks for log files in the equilibration/namd directory structure.
    """
    progress = {}

    # Standard stage names and their expected log file patterns
    # Updated to match current equilibration structure (6 equilibration stages + production)
    # step1* is for equilibration_1, step2* for equilibration_2, etc.
    stage_patterns = {
        "equilibration_1": ["step1_equilibration*.log", "step1_*.log"],
        "equilibration_2": ["step2_equilibration*.log", "step2_*.log"],
        "equilibration_3": ["step3_equilibration*.log", "step3_*.log"],
        "equilibration_4": ["step4_equilibration*.log", "step4_*.log"],
        "equilibration_5": ["step5_equilibration*.log", "step5_*.log"],
        "equilibration_6": ["step6_equilibration*.log", "step6_*.log"],
        "production": ["step7_production*.log", "production*.log", "prod*.log"],
    }

    # Look for equilibration directory
    eq_namd_dir = equilibration_dir / "namd"
    logger.debug(f"Checking primary path: {eq_namd_dir}")
    logger.debug(f"Primary path exists: {eq_namd_dir.exists()}")

    if not eq_namd_dir.exists():
        # Try alternative locations
        eq_namd_dir = equilibration_dir / "equilibration" / "namd"
        logger.debug(f"Trying fallback path: {eq_namd_dir}")
        logger.debug(f"Fallback path exists: {eq_namd_dir.exists()}")

        if not eq_namd_dir.exists():
            eq_namd_dir = equilibration_dir
            logger.debug(f"Using base directory: {eq_namd_dir}")

    logger.debug(f"Final search directory: {eq_namd_dir}")

    if eq_namd_dir.exists():
        all_files = list(eq_namd_dir.glob("*"))
        logger.debug(f" All files in directory: {[f.name for f in all_files]}")
        log_files = list(eq_namd_dir.glob("*.log"))

    for stage_name, patterns in stage_patterns.items():
        stage_progress = NAMDProgress(stage_name=stage_name)

        # Look for log files matching the patterns
        log_file = None
        for pattern in patterns:
            matches = list(eq_namd_dir.glob(pattern))
            if matches:
                # Use the most recent log file
                log_file = max(matches, key=lambda f: f.stat().st_mtime)
                break

        if log_file and log_file.exists():
            stage_progress.log_file = log_file
            stage_progress.last_updated = log_file.stat().st_mtime

            # Parse the log file
            timing = parse_namd_log(log_file)
            stage_progress.timing = timing

            log_text = ""
            try:
                from gatewizard.utils.equilibration_failure import failure_line_from_text

                log_text = _read_namd_log_text(log_file)
                fatal = failure_line_from_text(log_text)
            except Exception:
                fatal = None

            # Determine status and progress
            has_end = "End of program" in log_text
            if fatal:
                stage_progress.status = "error"
                if timing.total_steps > 0 and timing.steps_completed > 0:
                    stage_progress.progress_percent = (
                        timing.steps_completed / timing.total_steps
                    ) * 100.0
            elif has_end:
                # Mid-run restart writes ("WRITING … RESTART/DCD") are not completion.
                stage_progress.status = "completed"
                stage_progress.progress_percent = 100.0
            elif timing.steps_completed > 0:
                stage_progress.status = "running"
                if timing.total_steps > 0:
                    stage_progress.progress_percent = (
                        timing.steps_completed / timing.total_steps
                    ) * 100.0
                else:
                    stage_progress.progress_percent = 50.0
            else:
                stage_progress.status = "not_started"
        else:
            pass  # No log file found for stage

        progress[stage_name] = stage_progress

    return progress


def format_timing_info(timing: NAMDTiming) -> str:
    """Format timing information for display."""
    if not timing or timing.steps_completed == 0:
        return "No timing data available"

    lines = []

    # Basic info
    if timing.hostname:
        lines.append(f"Host: {timing.hostname}")

    if timing.processors > 0:
        gpu_info = f", {timing.gpus} GPUs" if timing.gpus > 0 else ""
        lines.append(f"Resources: {timing.processors} processors{gpu_info}")

    if timing.atoms > 0:
        lines.append(f"System: {timing.atoms:,} atoms")

    # Progress info
    if timing.total_steps > 0:
        lines.append(f"Steps: {timing.steps_completed:,} / {timing.total_steps:,}")
    else:
        lines.append(f"Steps completed: {timing.steps_completed:,}")

    # Performance info
    if timing.simulated_time_ns > 0:
        lines.append(f"Simulated time: {timing.simulated_time_ns:.3f} ns")

    if timing.real_time_hours > 0:
        lines.append(f"Real time: {timing.real_time_hours:.2f} hours")

    if timing.ns_per_day > 0:
        lines.append(f"Performance: {timing.ns_per_day:.4f} ns/day")

    if timing.sec_per_step > 0:
        lines.append(f"Time per step: {timing.sec_per_step:.3f} sec")

    return "\n".join(lines)


def format_progress_summary(progress_dict: Dict[str, NAMDProgress]) -> str:
    """Format a summary of all stage progress."""
    lines = ["Equilibration Progress Summary:"]

    for stage_name, progress in progress_dict.items():
        status_icon = {
            "not_started": "⏸️",
            "running": "🏃",
            "completed": "✅",
            "error": "❌",
        }.get(progress.status, "❓")

        stage_display = stage_name.replace("_", " ").title()
        line = f"  {status_icon} {stage_display}: {progress.status}"

        if progress.progress_percent > 0:
            line += f" ({progress.progress_percent:.1f}%)"

        if progress.timing and progress.timing.ns_per_day > 0:
            line += f" - {progress.timing.ns_per_day:.4f} ns/day"

        lines.append(line)

    return "\n".join(lines)


# ============================================================================
# High-level Wrapper Classes for Easy Analysis
# ============================================================================


class EnergyAnalyzer:
    """
    Easy-to-use wrapper for NAMD energy analysis with built-in plotting and full GUI capabilities.

    Supports:
    - Single or multiple log files with custom time scaling
    - Selective energy property plotting
    - Multiple plots (same figure or separate)
    - Full customization (colors, grid, units, target values, etc.)

    Example:
        >>> # Single file - plot
        >>> analyzer = EnergyAnalyzer("step1_equilibration.log")
        >>> analyzer.plot_energy(save="energy.png")

        >>> # With custom target temperature and pressure
        >>> analyzer.plot_energy(
        ...     target_temperature=310.0,  # 310 K
        ...     target_pressure=1.01325,   # 1.01325 atm (1 bar)
        ...     save="energy.png"
        ... )

        >>> # Multiple files with custom time
        >>> analyzer = EnergyAnalyzer(
        ...     ["step1.log", "step2.log", "step3.log"],
        ...     file_times={"step1.log": 0.05, "step2.log": 0.05, "step3.log": 0.05}
        ... )

        >>> # Plot specific properties
        >>> analyzer.plot_properties(
        ...     properties=["Temperature", "Total Energy", "Pressure"],
        ...     save="specific_energies.png"
        ... )

        >>> # Plot each property separately
        >>> analyzer.plot_properties(
        ...     properties=["Temperature", "Pressure"],
        ...     separate_plots=True,
        ...     save_prefix="plot_"
        ... )
    """

    def __init__(
        self,
        log_file: Union[Path, str, List[Union[Path, str]]],
        file_times: Optional[Dict[str, float]] = None,
        file_strides: Optional[Dict[str, int]] = None,
    ):
        """
        Initialize energy analyzer with NAMD log file(s).

        Args:
            log_file: Path to NAMD log file, or list of paths for multi-file analysis
            file_times: Dict mapping filename (just the name, not full path) to duration in ns
                       Example: {"step1.log": 0.05, "step2.log": 0.05}
                       Time is the DURATION of each file, not cumulative
            file_strides: Dict mapping filename to keep-every-N ENERGY samples (≥1)
        """
        # Handle single file or list
        if isinstance(log_file, (str, Path)):
            self.log_files = [Path(log_file)]
        else:
            self.log_files = [Path(f) for f in log_file]

        # Store file times (duration of each file in ns)
        self.file_times = file_times or {}
        self.file_strides = file_strides or {}

        # Track file ranges for time calculation (must be before parsing!)
        self._file_ranges = {}  # {filepath: (start_idx, end_idx, min_ts, max_ts)}

        # Parse all log files
        self.data = self._parse_energy_data()
        self.timing = parse_namd_log(self.log_files[0]) if self.log_files else None

    def _parse_energy_data(self) -> Dict[str, List[float]]:
        """Parse ENERGY lines from NAMD log file(s)."""
        data = {
            "timestep": [],
            "bond": [],
            "angle": [],
            "dihedral": [],
            "improper": [],
            "elect": [],
            "vdw": [],
            "boundary": [],
            "misc": [],
            "kinetic": [],
            "total": [],
            "temp": [],
            "potential": [],
            "total3": [],
            "tempavg": [],
            "pressure": [],
            "gpressure": [],
            "volume": [],
            "pressavg": [],
            "gpressavg": [],
        }

        for log_file in self.log_files:
            if not log_file.exists():
                logger.warning(f"Log file not found: {log_file}")
                continue

            start_idx = len(data["timestep"])

            try:
                with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                    content = f.read()

                # Find ETITLE line to get column labels
                etitle_match = re.search(r"ETITLE:\s+(.+)", content)
                if not etitle_match:
                    logger.warning(f"No ETITLE line found in {log_file}")
                    continue

                # Parse ENERGY lines
                energy_lines = re.findall(r"^ENERGY:\s+(.+)$", content, re.MULTILINE)
                stride = max(1, int(_lookup_file_map(self.file_strides, log_file) or 1))

                min_ts, max_ts = None, None

                for line_i, line in enumerate(energy_lines):
                    if line_i % stride != 0:
                        continue
                    values = line.split()
                    if len(values) >= 14:  # Minimum expected columns
                        try:
                            ts = int(values[0])
                            if min_ts is None:
                                min_ts = ts
                            max_ts = ts

                            data["timestep"].append(ts)
                            data["bond"].append(float(values[1]))
                            data["angle"].append(float(values[2]))
                            data["dihedral"].append(float(values[3]))
                            data["improper"].append(float(values[4]))
                            data["elect"].append(float(values[5]))
                            data["vdw"].append(float(values[6]))
                            data["boundary"].append(float(values[7]))
                            data["misc"].append(float(values[8]))
                            data["kinetic"].append(float(values[9]))
                            data["total"].append(float(values[10]))
                            data["temp"].append(float(values[11]))
                            data["potential"].append(float(values[12]))
                            data["total3"].append(float(values[13]))

                            if len(values) >= 15:
                                data["tempavg"].append(float(values[14]))
                            if len(values) >= 16:
                                data["pressure"].append(float(values[15]))
                            if len(values) >= 17:
                                data["gpressure"].append(float(values[16]))
                            if len(values) >= 18:
                                data["volume"].append(float(values[17]))
                            if len(values) >= 19:
                                data["pressavg"].append(float(values[18]))
                            if len(values) >= 20:
                                data["gpressavg"].append(float(values[19]))
                        except (ValueError, IndexError) as e:
                            continue

                # Store file range information
                end_idx = len(data["timestep"])
                self._file_ranges[str(log_file)] = (
                    start_idx,
                    end_idx,
                    min_ts or 0,
                    max_ts or 0,
                )

            except Exception as e:
                logger.error(f"Error parsing {log_file}: {e}")

        return data

    def get_statistics(self) -> Dict[str, Dict[str, float]]:
        """
        Get statistical summary of energy data.

        Returns:
            Dictionary with statistics for each energy component
        """
        import numpy as np

        stats = {}
        for key, values in self.data.items():
            if values and key != "timestep":
                arr = np.array(values)
                stats[key] = {
                    "mean": float(np.mean(arr)),
                    "std": float(np.std(arr)),
                    "min": float(np.min(arr)),
                    "max": float(np.max(arr)),
                    "initial": float(arr[0]) if len(arr) > 0 else 0.0,
                    "final": float(arr[-1]) if len(arr) > 0 else 0.0,
                }

        return stats

    def _calculate_time_array(self) -> "np.ndarray":
        """
        Calculate time array for all data points with proper file time scaling.

        Returns:
            Time array in nanoseconds
        """
        import numpy as np

        time_array_ns = []

        # Check if we have custom time assignments
        has_custom_times = bool(self.file_times) and any(
            t > 0 for t in self.file_times.values()
        )

        if has_custom_times and self._file_ranges:
            # Use custom time assignments
            cumulative_time_ns = 0.0

            for log_file in self.log_files:
                filepath_str = str(log_file)
                if filepath_str not in self._file_ranges:
                    continue

                start_idx, end_idx, min_ts, max_ts = self._file_ranges[filepath_str]
                num_points = end_idx - start_idx

                if num_points == 0:
                    continue

                # Get assigned time for this file (duration in ns)
                filename = log_file.name
                assigned_time_ns = self.file_times.get(filename, 0.0)

                if assigned_time_ns <= 0:
                    # No time assigned, use timestep-based calculation
                    timestep_fs = (
                        self.timing.timestep_fs
                        if self.timing and self.timing.timestep_fs > 0
                        else 2.0
                    )
                    for i in range(num_points):
                        time_array_ns.append(
                            cumulative_time_ns + i * timestep_fs / 1_000_000.0
                        )
                    cumulative_time_ns += num_points * timestep_fs / 1_000_000.0
                else:
                    # Distribute points evenly across assigned time
                    if num_points == 1:
                        time_array_ns.append(cumulative_time_ns)
                    else:
                        file_times = np.linspace(
                            cumulative_time_ns,
                            cumulative_time_ns + assigned_time_ns,
                            num_points,
                        )
                        time_array_ns.extend(file_times.tolist())

                    cumulative_time_ns += assigned_time_ns

            return np.array(time_array_ns)
        else:
            # Use timestep-based calculation
            timestep_fs = (
                self.timing.timestep_fs
                if self.timing and self.timing.timestep_fs > 0
                else 2.0
            )
            timesteps = np.array(self.data["timestep"])
            return timesteps * timestep_fs / 1_000_000.0  # Convert fs to ns

    def _normalize_property_name(self, prop_name: str) -> Optional[str]:
        """
        Normalize property name to internal data key (case-insensitive).

        Handles various input formats:
        - Display names: "Total Energy", "Temperature", etc.
        - NAMD column names: "TOTAL", "TEMP", "PRESSURE", etc.
        - Short names: "total", "temp", "pressure", etc.
        - Any case variation: "Total", "TOTAL", "total", "ToTaL", etc.

        Args:
            prop_name: Property name in any format/case

        Returns:
            Normalized internal key (lowercase) or None if not recognized
        """
        # Master mapping of all possible property names to internal keys
        # Using lowercase keys for everything
        property_mappings = {
            # Full display names
            "total energy": "total",
            "potential energy": "potential",
            "kinetic energy": "kinetic",
            "electrostatic energy": "elect",
            "van der waals energy": "vdw",
            "bond energy": "bond",
            "angle energy": "angle",
            "dihedral energy": "dihedral",
            "improper energy": "improper",
            "temperature": "temp",
            "pressure": "pressure",
            "volume": "volume",
            # Short forms (NAMD column names and common usage)
            "total": "total",
            "potential": "potential",
            "kinetic": "kinetic",
            "elect": "elect",
            "electrostatic": "elect",
            "vdw": "vdw",
            "bond": "bond",
            "angle": "angle",
            "dihedral": "dihedral",
            "improper": "improper",
            "temp": "temp",
            "pressure": "pressure",
            "volume": "volume",
            # Additional aliases
            "pot": "potential",
            "kin": "kinetic",
            "elec": "elect",
            "press": "pressure",
            "vol": "volume",
        }

        # Normalize input to lowercase for case-insensitive matching
        normalized_input = prop_name.lower().strip()

        # Direct lookup
        if normalized_input in property_mappings:
            return property_mappings[normalized_input]

        return None

    def get_available_properties(self) -> List[str]:
        """
        Get list of available energy properties that can be plotted.

        Returns:
            List of property names with units
        """
        available = []
        property_map = {
            "total": "Total Energy",
            "potential": "Potential Energy",
            "kinetic": "Kinetic Energy",
            "elect": "Electrostatic Energy",
            "vdw": "Van der Waals Energy",
            "bond": "Bond Energy",
            "angle": "Angle Energy",
            "dihedral": "Dihedral Energy",
            "improper": "Improper Energy",
            "temp": "Temperature",
            "pressure": "Pressure",
            "volume": "Volume",
        }

        for key, name in property_map.items():
            if key in self.data and self.data[key]:
                available.append(name)

        return available

    def plot_properties(
        self,
        properties: Optional[List[str]] = None,
        separate_plots: bool = False,
        energy_units: str = "kcal/mol",
        time_units: str = "ns",
        pressure_units: str = "atm",
        temperature_units: str = "K",
        volume_units: str = "Å³",
        line_colors: Optional[List[str]] = None,
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        grid_color: Optional[str] = None,
        show_grid: bool = True,
        xlim: Optional[tuple] = None,
        ylim: Optional[tuple] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        save: Optional[str] = None,
        save_prefix: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (10, 6),
        dpi: int = 300,
    ):
        """
        Plot selected energy properties with full GUI-level customization.

        Property names are case-insensitive and support multiple formats:
        - Full names: "Temperature", "Total Energy", "Pressure" (any case)
        - Short names: "TEMP", "TOTAL", "PRESSURE" (any case)
        - Aliases: "pot" (potential), "kin" (kinetic), "elec" (electrostatic)

        Args:
            properties: List of property names to plot. If None, plots 4-panel view.
                       Available: "Total Energy", "Potential Energy", "Kinetic Energy",
                                 "Electrostatic Energy", "Van der Waals Energy", "Bond Energy",
                                 "Angle Energy", "Dihedral Energy", "Improper Energy",
                                 "Temperature", "Pressure", "Volume"
                       Note: Property names are case-insensitive, so "TEMP", "temp",
                             "Temperature" all work.
            separate_plots: If True, create separate plot files for each property
            energy_units: "kcal/mol" or "kJ/mol"
            time_units: "ps", "ns", or "µs"
            pressure_units: "atm", "bar", "Pa", "kPa", "MPa"
            temperature_units: "K", "°C", or "°F"
            volume_units: "Å³", "nm³", "mL", "L"
            line_colors: List of colors for each property line
            bg_color: Plot area background color (hex or "none" for transparent)
            fig_bg_color: Figure border background color (hex or "none")
            text_color: Text/axes color ("Auto", color name, or hex)
            grid_color: Grid line color (None to match text_color)
            show_grid: Show grid lines
            xlim: X-axis limits (min, max)
            ylim: Y-axis limits (min, max)
            title: Plot title (auto-generated if None)
            xlabel: X-axis label (auto-generated if None)
            ylabel: Y-axis label (auto-generated if None)
            save: Filename to save plot (only for single plot or combined)
            save_prefix: Prefix for filenames when separate_plots=True
            show: Display plot interactively
            figsize: Figure size (width, height) in inches
            dpi: Resolution for saved figure

        Example:
            >>> # Plot specific properties on same figure (case-insensitive)
            >>> analyzer.plot_properties(
            ...     properties=["TEMP", "total energy", "Pressure"],  # Any case works!
            ...     line_colors=["red", "blue", "green"],
            ...     save="combined.png"
            ... )

            >>> # Plot each property separately
            >>> analyzer.plot_properties(
            ...     properties=["Temperature", "Pressure"],
            ...     separate_plots=True,
            ...     save_prefix="energy_"
            ... )
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return

        if not self.data["timestep"]:
            logger.warning("No energy data to plot")
            return

        # If no properties specified, use the 4-panel plot
        if properties is None:
            return self.plot_energy(
                energy_units=energy_units,
                time_units=time_units,
                bg_color=bg_color,
                fig_bg_color=fig_bg_color,
                text_color=text_color,
                show_grid=show_grid,
                title=title,
                save=save,
                show=show,
                figsize=figsize,
                dpi=dpi,
            )

        # Calculate time array with proper file scaling
        time_ns = self._calculate_time_array()

        # Convert time units
        if time_units == "ps":
            plot_time = time_ns * 1000.0
        elif time_units == "µs":
            plot_time = time_ns / 1000.0
        else:  # ns
            plot_time = time_ns

        # Default line colors
        if line_colors is None:
            line_colors = [
                "blue",
                "red",
                "green",
                "orange",
                "purple",
                "brown",
                "pink",
                "cyan",
            ]

        from gatewizard.utils.plot_spec import plot_spec_from_plot_properties_kwargs
        from gatewizard.utils import matplotlib_renderer

        # Build series payload for renderer (units already applied via analyzer data)
        series_payload: List[Dict[str, Any]] = []
        for prop_name in properties:
            data_key = self._normalize_property_name(prop_name)
            if not data_key or data_key not in self.data or not self.data[data_key]:
                logger.warning(f"Property '{prop_name}' not available or not recognized")
                continue
            import numpy as np

            y_data = np.array(self.data[data_key])
            y_data, unit_label = self._convert_property_units(
                data_key,
                y_data,
                energy_units,
                pressure_units,
                temperature_units,
                volume_units,
            )
            series_payload.append(
                {
                    "key": data_key,
                    "name": prop_name,
                    "unit": unit_label,
                    "y": y_data.tolist(),
                    "x": plot_time.tolist() if hasattr(plot_time, "tolist") else list(plot_time),
                }
            )

        if not series_payload:
            logger.warning("No plottable properties resolved")
            return

        plot_spec = plot_spec_from_plot_properties_kwargs(
            [s["name"] for s in series_payload],
            separate_plots=separate_plots,
            line_colors=line_colors,
            energy_units=energy_units,
            time_units=time_units,
            pressure_units=pressure_units,
            temperature_units=temperature_units,
            volume_units=volume_units,
            bg_color=bg_color,
            fig_bg_color=fig_bg_color,
            text_color=text_color,
            grid_color=grid_color,
            show_grid=show_grid,
            xlim=xlim,
            ylim=ylim,
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            figsize=figsize,
            dpi=dpi,
        )
        # Attach resolved keys/colors to panels
        for i, panel in enumerate(plot_spec["panels"]):
            if i < len(series_payload):
                panel["key"] = series_payload[i]["key"]
                panel["ylabel"] = ylabel or f"{series_payload[i]['name']} ({series_payload[i]['unit']})"

        data = {"x": series_payload[0]["x"], "series": series_payload}

        if separate_plots:
            prefix = save_prefix or "plot_"
            for i, panel in enumerate(plot_spec["panels"]):
                single_spec = plot_spec_from_plot_properties_kwargs(
                    [panel.get("name") or panel["key"]],
                    separate_plots=False,
                    line_colors=[panel.get("line_color")],
                    energy_units=energy_units,
                    time_units=time_units,
                    pressure_units=pressure_units,
                    temperature_units=temperature_units,
                    volume_units=volume_units,
                    bg_color=bg_color,
                    fig_bg_color=fig_bg_color,
                    text_color=text_color,
                    grid_color=grid_color,
                    show_grid=show_grid,
                    xlim=xlim,
                    ylim=ylim,
                    title=title or panel.get("name"),
                    xlabel=xlabel,
                    ylabel=panel.get("ylabel"),
                    figsize=figsize,
                    dpi=dpi,
                )
                single_spec["panels"][0]["key"] = panel["key"]
                fig = matplotlib_renderer.render_energetic(data, single_spec)
                import matplotlib.pyplot as plt

                safe_name = str(panel.get("name") or panel["key"]).lower().replace(" ", "_")
                filename = f"{prefix}{safe_name}.png"
                fig.savefig(filename, dpi=dpi, bbox_inches="tight")
                plt.close(fig)
                logger.info(f"Plot saved: {filename}")
                if show:
                    plt.show()
        else:
            fig = matplotlib_renderer.render_energetic(data, plot_spec)
            import matplotlib.pyplot as plt

            if save:
                fig.savefig(save, dpi=dpi, bbox_inches="tight")
                logger.info(f"Plot saved: {save}")
            if show:
                plt.show()
            plt.close(fig)

    def _auto_text_color(self, bg_color: str) -> str:
        """Auto-determine text color based on background luminance."""
        if bg_color == "none":
            return "black"
        try:
            hex_color = bg_color.lstrip("#")
            r, g, b = (
                int(hex_color[0:2], 16),
                int(hex_color[2:4], 16),
                int(hex_color[4:6], 16),
            )
            luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
            return "black" if luminance > 0.5 else "white"
        except:
            return "white"

    def _convert_property_units(
        self,
        data_key: str,
        values: "np.ndarray",
        energy_units: str,
        pressure_units: str,
        temperature_units: str,
        volume_units: str,
    ) -> tuple:
        """Convert property values to requested units and return (values, unit_label)."""
        import numpy as np

        # Energy properties
        if data_key in [
            "total",
            "potential",
            "kinetic",
            "elect",
            "vdw",
            "bond",
            "angle",
            "dihedral",
            "improper",
        ]:
            if energy_units == "kJ/mol":
                return values * 4.184, "kJ/mol"
            return values, "kcal/mol"

        # Temperature
        elif data_key == "temp":
            if temperature_units == "°C":
                return values - 273.15, "°C"
            elif temperature_units == "°F":
                return (values - 273.15) * 9 / 5 + 32, "°F"
            return values, "K"

        # Pressure
        elif data_key == "pressure":
            conversions = {
                "atm": 1.0,
                "bar": 1.01325,
                "Pa": 101325.0,
                "kPa": 101.325,
                "MPa": 0.101325,
            }
            factor = conversions.get(pressure_units, 1.0)
            return values * factor, pressure_units

        # Volume
        elif data_key == "volume":
            conversions = {"Å³": 1.0, "nm³": 0.001, "mL": 1.66054e-24, "L": 1.66054e-27}
            factor = conversions.get(volume_units, 1.0)
            return values * factor, volume_units

        return values, ""

    def _plot_single_property(
        self,
        prop_name,
        plot_time,
        time_units,
        energy_units,
        pressure_units,
        temperature_units,
        volume_units,
        line_color,
        bg_color,
        fig_bg_color,
        text_color,
        grid_color,
        show_grid,
        xlim,
        ylim,
        title,
        xlabel,
        ylabel,
        save_prefix,
        show,
        figsize,
        dpi,
    ):
        """Plot a single property in its own figure."""
        import matplotlib.pyplot as plt
        import numpy as np

        data_key = self._normalize_property_name(prop_name)
        if not data_key or data_key not in self.data or not self.data[data_key]:
            logger.warning(f"Property '{prop_name}' not available")
            return

        fig, ax = plt.subplots(figsize=figsize)

        # Set backgrounds
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        if bg_color != "none":
            ax.set_facecolor(bg_color)

        # Get data and convert units
        y_data = np.array(self.data[data_key])
        y_data, unit_label = self._convert_property_units(
            data_key,
            y_data,
            energy_units,
            pressure_units,
            temperature_units,
            volume_units,
        )

        # Plot
        ax.plot(plot_time, y_data, color=line_color, linewidth=1.5)

        # Styling
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)

        if show_grid:
            ax.grid(True, alpha=0.3, color=grid_color)

        # Labels
        ax.set_xlabel(xlabel or f"Time ({time_units})", color=text_color)
        ax.set_ylabel(ylabel or f"{prop_name} ({unit_label})", color=text_color)
        ax.set_title(title or prop_name, color=text_color, fontweight="bold")

        # Limits
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        plt.tight_layout()

        # Save
        if save_prefix:
            safe_name = prop_name.lower().replace(" ", "_")
            filename = f"{save_prefix}{safe_name}.png"
            plt.savefig(filename, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {filename}")

        if show:
            plt.show()
        else:
            plt.close()

    def _plot_combined_properties(
        self,
        properties,
        plot_time,
        time_units,
        energy_units,
        pressure_units,
        temperature_units,
        volume_units,
        line_colors,
        bg_color,
        fig_bg_color,
        text_color,
        grid_color,
        show_grid,
        xlim,
        ylim,
        title,
        xlabel,
        ylabel,
        save,
        show,
        figsize,
        dpi,
    ):
        """Plot multiple properties on the same figure."""
        import matplotlib.pyplot as plt
        import numpy as np

        fig, ax = plt.subplots(figsize=figsize)

        # Set backgrounds
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        if bg_color != "none":
            ax.set_facecolor(bg_color)

        # Plot each property
        lines = []
        labels = []

        for i, prop_name in enumerate(properties):
            data_key = self._normalize_property_name(prop_name)
            if not data_key or data_key not in self.data or not self.data[data_key]:
                logger.warning(
                    f"Property '{prop_name}' not available or not recognized"
                )
                continue

            y_data = np.array(self.data[data_key])
            y_data, unit_label = self._convert_property_units(
                data_key,
                y_data,
                energy_units,
                pressure_units,
                temperature_units,
                volume_units,
            )

            color = line_colors[i % len(line_colors)]
            line = ax.plot(
                plot_time,
                y_data,
                color=color,
                linewidth=1.5,
                marker="o",
                markersize=2,
                label=f"{prop_name} ({unit_label})",
            )
            lines.extend(line)
            labels.append(f"{prop_name} ({unit_label})")

        # Styling
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)

        if show_grid:
            ax.grid(True, alpha=0.3, color=grid_color)

        # Labels
        ax.set_xlabel(xlabel or f"Time ({time_units})", color=text_color)
        ax.set_ylabel(ylabel or "Multiple Properties", color=text_color)
        ax.set_title(
            title or f"NAMD Analysis - {len(properties)} Properties",
            color=text_color,
            fontweight="bold",
        )

        # Legend
        if len(properties) > 1:
            legend = ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.setp(legend.get_texts(), color=text_color)

        # Limits
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        plt.tight_layout()

        # Save
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()

    def plot_energy(
        self,
        energy_units: str = "kcal/mol",
        time_units: str = "ns",
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        show_grid: bool = True,
        title: Optional[str] = None,
        target_temperature: Optional[float] = None,
        target_pressure: Optional[float] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (12, 10),
        dpi: int = 300,
    ):
        """
        Create energy analysis plot with full GUI customization.

        Args:
            energy_units: 'kcal/mol' or 'kJ/mol'
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            bg_color: Background color for plot area
            fig_bg_color: Background color for figure border
            text_color: Text/axes color ('Auto' or specific color)
            show_grid: Show grid lines on plots
            title: Main title for the figure (default: auto-generated)
            target_temperature: Target temperature in Kelvin (default: auto-calculated from last 50% of trajectory)
            target_pressure: Target pressure in atm (default: 1.0 atm)
            save: Filename to save plot (e.g., "energy.png")
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return

        if not self.data["timestep"]:
            logger.warning("No energy data to plot")
            return

        # Calculate time array with proper file scaling
        time_ns = self._calculate_time_array()

        # Convert time units
        if time_units == "ps":
            plot_time = time_ns * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = time_ns / 1000.0  # Convert ns to µs
        else:  # ns
            plot_time = time_ns

        # Energy unit conversion factor
        if energy_units == "kJ/mol":
            energy_factor = 4.184  # kcal/mol to kJ/mol
        else:  # kcal/mol
            energy_factor = 1.0

        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=figsize)

        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)

        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = "black"
            else:
                try:
                    hex_color = bg_color.lstrip("#")
                    r, g, b = (
                        int(hex_color[0:2], 16),
                        int(hex_color[2:4], 16),
                        int(hex_color[4:6], 16),
                    )
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = "black" if luminance > 0.5 else "white"
                except:
                    text_color = "white"

        # Set title
        if len(self.log_files) == 1:
            title_text = title or f"Energy Analysis - {self.log_files[0].name}"
        else:
            title_text = title or f"Energy Analysis - {len(self.log_files)} Files"
        fig.suptitle(title_text, fontsize=14, fontweight="bold", color=text_color)

        # Configure all subplots with common settings
        for ax in axes.flat:
            if bg_color != "none":
                ax.set_facecolor(bg_color)
            ax.tick_params(colors=text_color)
            for spine in ax.spines.values():
                spine.set_color(text_color)

        # Panel 1: Total Energy
        if self.data["total"]:
            energy_data = np.array(self.data["total"]) * energy_factor / 1000
            axes[0, 0].plot(plot_time, energy_data, "b-", linewidth=0.8)
            axes[0, 0].set_xlabel(f"Time ({time_units})", color=text_color)
            axes[0, 0].set_ylabel(
                f"Total Energy (×10³ {energy_units})", color=text_color
            )
            axes[0, 0].set_title("Total Energy Convergence", color=text_color)
            if show_grid:
                axes[0, 0].grid(True, alpha=0.3, color=text_color)

        # Panel 2: Potential and Kinetic Energy
        if self.data["potential"] and self.data["kinetic"]:
            pot_energy = np.array(self.data["potential"]) * energy_factor / 1000
            kin_energy = np.array(self.data["kinetic"]) * energy_factor / 1000
            axes[0, 1].plot(
                plot_time, pot_energy, "r-", linewidth=0.8, label="Potential", alpha=0.8
            )
            axes[0, 1].plot(
                plot_time, kin_energy, "g-", linewidth=0.8, label="Kinetic", alpha=0.8
            )
            axes[0, 1].set_xlabel(f"Time ({time_units})", color=text_color)
            axes[0, 1].set_ylabel(f"Energy (×10³ {energy_units})", color=text_color)
            axes[0, 1].set_title("Potential and Kinetic Energy", color=text_color)
            legend = axes[0, 1].legend()
            plt.setp(legend.get_texts(), color=text_color)
            if show_grid:
                axes[0, 1].grid(True, alpha=0.3, color=text_color)

        # Panel 3: Temperature
        if self.data["temp"]:
            temp_array = np.array(self.data["temp"])
            # Use user-provided target or auto-calculate from last 50% of trajectory
            target_temp = (
                target_temperature if target_temperature is not None else 300.0
            )
            axes[1, 0].plot(plot_time, temp_array, "orange", linewidth=0.8)
            axes[1, 0].axhline(
                y=target_temp,
                color=text_color,
                linestyle="--",
                linewidth=1,
                label=f"Target: {target_temp:.1f} K",
                alpha=0.7,
            )
            axes[1, 0].set_xlabel(f"Time ({time_units})", color=text_color)
            axes[1, 0].set_ylabel("Temperature (K)", color=text_color)
            axes[1, 0].set_title("Temperature Stability", color=text_color)
            legend = axes[1, 0].legend()
            plt.setp(legend.get_texts(), color=text_color)
            if show_grid:
                axes[1, 0].grid(True, alpha=0.3, color=text_color)

        # Panel 4: Pressure (if available)
        if self.data["pressure"]:
            # Use user-provided target or default to 1.0 atm
            target_press = target_pressure if target_pressure is not None else 1.0
            axes[1, 1].plot(plot_time, self.data["pressure"], "purple", linewidth=0.8)
            axes[1, 1].axhline(
                y=target_press,
                color=text_color,
                linestyle="--",
                linewidth=1,
                label=f"Target: {target_press:.1f} atm",
                alpha=0.7,
            )
            axes[1, 1].set_xlabel(f"Time ({time_units})", color=text_color)
            axes[1, 1].set_ylabel("Pressure (atm)", color=text_color)
            axes[1, 1].set_title("Pressure Fluctuations", color=text_color)
            legend = axes[1, 1].legend()
            plt.setp(legend.get_texts(), color=text_color)
            if show_grid:
                axes[1, 1].grid(True, alpha=0.3, color=text_color)
        else:
            # If no pressure, show energy components
            if self.data["elect"] and self.data["vdw"]:
                elect_energy = np.array(self.data["elect"]) * energy_factor / 1000
                vdw_energy = np.array(self.data["vdw"]) * energy_factor / 1000
                axes[1, 1].plot(
                    plot_time,
                    elect_energy,
                    "b-",
                    linewidth=0.8,
                    label="Electrostatic",
                    alpha=0.7,
                )
                axes[1, 1].plot(
                    plot_time,
                    vdw_energy,
                    "r-",
                    linewidth=0.8,
                    label="van der Waals",
                    alpha=0.7,
                )
                axes[1, 1].set_xlabel(f"Time ({time_units})", color=text_color)
                axes[1, 1].set_ylabel(f"Energy (×10³ {energy_units})", color=text_color)
                axes[1, 1].set_title("Energy Components", color=text_color)
                legend = axes[1, 1].legend()
                plt.setp(legend.get_texts(), color=text_color)
                if show_grid:
                    axes[1, 1].grid(True, alpha=0.3, color=text_color)

        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()


def _to_path_list(paths: List[Union[str, Path]]) -> List[Path]:
    """Normalize a list of filesystem paths to resolved Path objects."""
    return [Path(p).expanduser().resolve() for p in paths]


def list_namd_energy_properties(
    log_files: List[Union[str, Path]],
    file_times: Optional[Dict[str, float]] = None,
) -> List[str]:
    """Return available NAMD ENERGY properties detected from log files.

    Fast path: scan for ``ETITLE`` / first ``ENERGY`` line instead of parsing
    every ENERGY sample (full parse is reserved for run_energetic_analysis).
    """
    del file_times
    logs = _to_path_list(log_files)
    property_map = {
        "TOTAL": "Total Energy",
        "POTENTIAL": "Potential Energy",
        "KINETIC": "Kinetic Energy",
        "ELECT": "Electrostatic Energy",
        "VDW": "Van der Waals Energy",
        "BOND": "Bond Energy",
        "ANGLE": "Angle Energy",
        "DIHED": "Dihedral Energy",
        "IMPRP": "Improper Energy",
        "TEMP": "Temperature",
        "PRESSURE": "Pressure",
        "VOLUME": "Volume",
    }
    # Fallback order matching EnergyAnalyzer.get_available_properties
    default_order = [
        "Total Energy",
        "Potential Energy",
        "Kinetic Energy",
        "Electrostatic Energy",
        "Van der Waals Energy",
        "Bond Energy",
        "Angle Energy",
        "Dihedral Energy",
        "Improper Energy",
        "Temperature",
        "Pressure",
        "Volume",
    ]

    for log_file in logs:
        if not log_file.is_file():
            continue
        try:
            with open(log_file, "r", encoding="utf-8", errors="ignore") as fh:
                etitle_cols: Optional[List[str]] = None
                for _ in range(5000):
                    line = fh.readline()
                    if not line:
                        break
                    if line.startswith("ETITLE:"):
                        etitle_cols = line.split()[1:]  # drop ETITLE:
                        break
                    if line.startswith("ENERGY:"):
                        # No ETITLE yet — assume standard NAMD columns present
                        return [
                            "Total Energy",
                            "Potential Energy",
                            "Kinetic Energy",
                            "Temperature",
                            "Pressure",
                            "Volume",
                        ]
                if etitle_cols:
                    found = []
                    for col in etitle_cols:
                        name = property_map.get(col.upper())
                        if name and name not in found:
                            found.append(name)
                    if found:
                        # Stable order
                        return [n for n in default_order if n in found]
        except OSError as exc:
            logger.warning(f"Cannot peek NAMD log {log_file}: {exc}")
            continue

    # Last resort: full parse of first file
    if logs:
        analyzer = EnergyAnalyzer([logs[0]])
        return analyzer.get_available_properties()
    return []


def run_energetic_analysis(
    log_files: List[Union[str, Path]],
    properties: Optional[List[str]] = None,
    file_times: Optional[Dict[str, float]] = None,
    file_strides: Optional[Dict[str, int]] = None,
    time_units: str = "ns",
    energy_units: str = "kcal/mol",
    pressure_units: str = "atm",
    temperature_units: str = "K",
    volume_units: str = "Å³",
) -> Dict[str, Any]:
    """
    Run NAMD log energetic analysis and return JSON-serializable arrays.

    Returns:
        Dict with `x`, `x_label`, and `series` entries, where each series has
        `name`, `key`, `unit`, and `y`.
    """
    import numpy as np

    logs = _to_path_list(log_files)
    analyzer = EnergyAnalyzer(logs, file_times=file_times, file_strides=file_strides)

    # Time in ns from analyzer, then convert for display
    x = analyzer._calculate_time_array()
    x_label = "Time (ns)"
    if time_units == "ps":
        x = x * 1000.0
        x_label = "Time (ps)"
    elif time_units in {"us", "µs"}:
        x = x / 1000.0
        x_label = "Time (µs)"

    selected = properties or analyzer.get_available_properties()
    series = []
    for prop in selected:
        key = analyzer._normalize_property_name(prop)
        if key is None:
            continue
        raw = analyzer.data.get(key, [])
        if not raw:
            continue
        values, unit = analyzer._convert_property_units(
            key,
            np.asarray(raw, dtype=float),
            energy_units,
            pressure_units,
            temperature_units,
            volume_units,
        )
        series.append(
            {
                "name": prop,
                "key": key,
                "unit": unit,
                "y": values.tolist(),
            }
        )

    return {
        "x": x.tolist(),
        "x_label": x_label,
        "series": series,
        "statistics": analyzer.get_statistics(),
    }
