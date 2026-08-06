"""
Amber mdout analysis utilities for equilibration progress and energetic plots.

Mirrors the GROMACS / OpenMM analysis interfaces used by the GUI status and
Analysis energetic endpoints.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from .energy_stride import lookup_file_map
from .logger import get_logger

logger = get_logger(__name__)


@dataclass
class AmberTimingInfo:
    """Timing/progress container — field names aligned with app.py expectations."""

    steps_completed: int = 0
    total_steps: int = 0
    timestep_fs: float = 0.0
    ns_per_day: float = 0.0
    is_minimization: bool = False
    wall_elapsed_seconds: float = 0.0
    converged_early: bool = False
    completed: bool = False
    has_error: bool = False
    interrupted: bool = False


@dataclass
class AmberStageProgress:
    stage_name: str = ""
    status: str = "not_started"  # not_started | running | completed | error
    timing: Optional[AmberTimingInfo] = None
    log_file: Optional[Path] = None


_NSTEP_RE = re.compile(
    r"NSTEP\s*=\s*(\d+)\s+TIME\(PS\)\s*=\s*([\d.Ee+-]+)"
    r"(?:\s+TEMP\(K\)\s*=\s*([\d.Ee+-]+))?"
    r"(?:\s+PRESS\s*=\s*([\d.Ee+-]+))?",
    re.IGNORECASE,
)
_ENERGY_LINE_RE = re.compile(
    r"((?:1-4\s+)?[A-Za-z][A-Za-z0-9]*)\s*=\s*"
    r"([-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?)"
)
_MIN_HEADER_RE = re.compile(
    r"^\s*NSTEP\s+ENERGY\s+RMS\s+GMAX\s+NAME\s+NUMBER\s*$",
    re.MULTILINE | re.IGNORECASE,
)
_MIN_ROW_RE = re.compile(
    r"^\s*(\d+)\s+"
    r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?|\*{3,})\s+"
)
_ERROR_MARKERS = (
    "FATAL ERROR",
    "Error terminating",
    "INFINITE ENERGY",
    "NaN",
    "******",
)

# mdinfo "Current Timing Info" — prefer cumulative "all steps" over last-window sample.
_MDINFO_ALL_STEPS_TIMING_RE = re.compile(
    r"Average timings for all steps:.*?"
    r"Elapsed\(s\)\s*=\s*([\d.Ee+-]+).*?"
    r"ns/day\s*=\s*([\d.Ee+-]+)",
    re.I | re.DOTALL,
)
_MDINFO_LAST_STEPS_TIMING_RE = re.compile(
    r"Average timings for last\s+\d+\s+steps:.*?"
    r"Elapsed\(s\)\s*=\s*([\d.Ee+-]+).*?"
    r"ns/day\s*=\s*([\d.Ee+-]+)",
    re.I | re.DOTALL,
)
_MDINFO_COMPLETED_STEPS_RE = re.compile(
    r"Completed\s*:\s*(\d+)",
    re.I,
)


def _parse_mdinfo_timings(mdinfo_content: str) -> tuple[float, float]:
    """Return (elapsed_wall_s, ns_per_day) from Amber mdinfo timing blocks."""
    if not mdinfo_content:
        return 0.0, 0.0
    for pattern in (_MDINFO_ALL_STEPS_TIMING_RE, _MDINFO_LAST_STEPS_TIMING_RE):
        match = pattern.search(mdinfo_content)
        if match:
            return float(match.group(1)), float(match.group(2))
    return 0.0, 0.0


def _parse_mdin_totals(mdin_path: Optional[Path]) -> tuple[int, float, bool]:
    """Return (total_steps, timestep_fs, is_minimization) from an mdin file."""
    if mdin_path is None or not mdin_path.is_file():
        return 0, 0.0, False
    try:
        text = mdin_path.read_text(encoding="utf-8", errors="ignore")
    except OSError:
        return 0, 0.0, False
    is_min = bool(re.search(r"\bimin\s*=\s*1\b", text, re.I))
    if is_min:
        m = re.search(r"\bmaxcyc\s*=\s*(\d+)", text, re.I)
        return (int(m.group(1)) if m else 0), 0.0, True
    nstlim_m = re.search(r"\bnstlim\s*=\s*(\d+)", text, re.I)
    dt_m = re.search(r"\bdt\s*=\s*([\d.Ee+-]+)", text, re.I)
    nstlim = int(nstlim_m.group(1)) if nstlim_m else 0
    dt_ps = float(dt_m.group(1)) if dt_m else 0.0
    return nstlim, dt_ps * 1000.0, False


_MIN_STEP_RE = re.compile(
    r"^\s{0,6}(\d+)\s+[-+]?\d+\.\d+E[+-]\d+",
    re.MULTILINE,
)


def _steps_from_amber_text(content: str, *, is_minimization: bool) -> int:
    """Extract the latest NSTEP from mdout/mdinfo text."""
    if is_minimization:
        min_steps = _MIN_STEP_RE.findall(content)
        if min_steps:
            return int(min_steps[-1])
        # Some builds also echo NSTEP = N in mdinfo during minimization
    nstep_matches = list(_NSTEP_RE.finditer(content))
    if nstep_matches:
        return int(nstep_matches[-1].group(1))
    return 0


def parse_amber_mdout(
    mdout_file: Path,
    *,
    is_minimization: bool = False,
    mdin_file: Optional[Path] = None,
    mdinfo_file: Optional[Path] = None,
) -> AmberTimingInfo:
    """Parse an Amber ``.mdout`` (and optional ``.mdinfo``) for progress.

    Amber's ``.mdinfo`` is a short status file rewritten every ``ntpr`` steps —
    especially useful during minimization, when ``.mdout`` may lag due to
    Fortran / pmemd buffering. Step counts take the maximum of mdout and mdinfo.

    Completion requires a normal end marker in ``.mdout`` **and** reaching the
    requested step count when known. Partial kills must not report
    ``completed=True`` just because a restart was written.
    """
    info = AmberTimingInfo(is_minimization=is_minimization)
    mdout_file = Path(mdout_file)
    if mdinfo_file is None:
        mdinfo_file = mdout_file.with_suffix(".mdinfo")
    else:
        mdinfo_file = Path(mdinfo_file)

    if mdin_file is None:
        cand = mdout_file.with_suffix(".mdin")
        if cand.is_file():
            mdin_file = cand

    total_from_mdin, dt_from_mdin, min_from_mdin = _parse_mdin_totals(mdin_file)
    if min_from_mdin:
        info.is_minimization = True
    if total_from_mdin:
        info.total_steps = total_from_mdin
    if dt_from_mdin:
        info.timestep_fs = dt_from_mdin

    content = ""
    if mdout_file.is_file():
        try:
            content = mdout_file.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            content = ""

    mdinfo_content = ""
    if mdinfo_file.is_file():
        try:
            mdinfo_content = mdinfo_file.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            mdinfo_content = ""

    if not content and not mdinfo_content:
        return info

    combined_for_errors = content + "\n" + mdinfo_content
    if re.search(r"FATAL ERROR|Error terminating|INFINITE ENERGY", combined_for_errors, re.I):
        info.has_error = True

    steps_mdout = _steps_from_amber_text(content, is_minimization=info.is_minimization)
    steps_mdinfo = _steps_from_amber_text(
        mdinfo_content, is_minimization=info.is_minimization
    )
    info.steps_completed = max(steps_mdout, steps_mdinfo)
    completed_m = _MDINFO_COMPLETED_STEPS_RE.search(mdinfo_content)
    if completed_m:
        info.steps_completed = max(info.steps_completed, int(completed_m.group(1)))

    if not info.is_minimization:
        # Prefer dt / nstlim echoed in mdout control section when present
        dt_m = re.search(r"\bdt\s*=\s*([\d.Ee+-]+)", content, re.I)
        if dt_m and info.timestep_fs <= 0:
            info.timestep_fs = float(dt_m.group(1)) * 1000.0
        nstlim_m = re.search(r"\bnstlim\s*=\s*(\d+)", content, re.I)
        if nstlim_m and info.total_steps <= 0:
            info.total_steps = int(nstlim_m.group(1))

    # Live cumulative timings from mdinfo (pmemd flushes ~every ntpr / 60s).
    mdinfo_elapsed, mdinfo_nsday = _parse_mdinfo_timings(mdinfo_content)
    if mdinfo_nsday > 0:
        info.ns_per_day = mdinfo_nsday
    if mdinfo_elapsed > 0:
        info.wall_elapsed_seconds = mdinfo_elapsed

    # Wall time / ns/day from mdout TIMINGS when present (usually at stage end).
    nsday_m = re.search(
        r"(?:ns/day|ns per day)\s*[:=]\s*([\d.Ee+-]+)", content, re.I
    )
    if nsday_m:
        info.ns_per_day = float(nsday_m.group(1))

    wall_m = re.search(
        r"(?:Elapsed\(wallclock\)|Elapsed time|Total wall time)\s*[:=]?\s*"
        r"([\d.]+)\s*(?:seconds|s)?",
        content,
        re.I,
    )
    if wall_m:
        info.wall_elapsed_seconds = float(wall_m.group(1))

    has_final = bool(
        re.search(r"Final Results|TIMINGS|Final Performance", content, re.I)
    )
    reached = (
        info.total_steps > 0 and info.steps_completed >= info.total_steps
    ) or (
        info.is_minimization
        and has_final
        and info.steps_completed > 0
        and (info.total_steps == 0 or info.steps_completed >= info.total_steps)
    )

    if has_final and (reached or (info.is_minimization and info.steps_completed > 0)):
        if info.total_steps > 0 and info.steps_completed < info.total_steps:
            # Finished banner without reaching nstlim → interrupted / early stop
            info.interrupted = True
            info.completed = False
        else:
            info.completed = True
            if info.total_steps > 0 and info.steps_completed == 0:
                info.steps_completed = info.total_steps
    elif has_final and info.total_steps > 0 and info.steps_completed < info.total_steps:
        info.interrupted = True

    # ns/day from wall when neither mdinfo nor TIMINGS reported it.
    if (
        info.ns_per_day == 0.0
        and info.steps_completed > 0
        and info.timestep_fs > 0
        and info.wall_elapsed_seconds > 1.0
    ):
        sim_ns = info.steps_completed * info.timestep_fs * 1e-6
        days = info.wall_elapsed_seconds / 86400.0
        if days > 0:
            info.ns_per_day = sim_ns / days
    elif (
        info.wall_elapsed_seconds == 0.0
        and info.ns_per_day > 0
        and info.steps_completed > 0
        and info.timestep_fs > 0
    ):
        sim_ns = info.steps_completed * info.timestep_fs * 1e-6
        info.wall_elapsed_seconds = (sim_ns / info.ns_per_day) * 86400.0

    return info


def get_equilibration_progress(
    equilibration_dir: Path,
) -> Dict[str, AmberStageProgress]:
    """Return progress for standard Amber equilibration stages.

    While a stage is running, ``log_file`` points at ``.mdinfo`` when present
    (live snapshot rewritten each ``ntpr``). After completion it points at
    ``.mdout`` (full history, subject to pmemd flush buffering).
    """
    stage_log_map: Dict[str, str] = {
        "minimization": "step0_minimization.mdout",
        "equilibration_1": "step1_equilibration.mdout",
        "equilibration_2": "step2_equilibration.mdout",
        "equilibration_3": "step3_equilibration.mdout",
        "equilibration_4": "step4_equilibration.mdout",
        "equilibration_5": "step5_equilibration.mdout",
        "equilibration_6": "step6_equilibration.mdout",
        "production": "step7_production.mdout",
    }

    progress: Dict[str, AmberStageProgress] = {}
    equilibration_dir = Path(equilibration_dir)

    for stage_name, log_name in stage_log_map.items():
        stage = AmberStageProgress(stage_name=stage_name)
        log_file = equilibration_dir / log_name
        mdinfo = equilibration_dir / log_name.replace(".mdout", ".mdinfo")
        if log_file.exists() or mdinfo.exists():
            mdin = equilibration_dir / log_name.replace(".mdout", ".mdin")
            timing = parse_amber_mdout(
                log_file if log_file.exists() else mdinfo,
                is_minimization=(stage_name == "minimization"),
                mdin_file=mdin if mdin.is_file() else None,
                mdinfo_file=mdinfo if mdinfo.is_file() else None,
            )
            stage.timing = timing
            if timing.has_error or timing.interrupted:
                stage.status = "error"
            elif timing.completed:
                stage.status = "completed"
            else:
                stage.status = "running"
            # Live UI: mdinfo is the current ntpr snapshot; mdout is history.
            if stage.status == "running" and mdinfo.is_file():
                stage.log_file = mdinfo
            elif log_file.exists():
                stage.log_file = log_file
            elif mdinfo.is_file():
                stage.log_file = mdinfo
        progress[stage_name] = stage

    return progress


# ---------------------------------------------------------------------------
# Energetic analysis
# ---------------------------------------------------------------------------

_AMBER_ENERGY_KEYS = frozenset(
    {
        "Etot",
        "EKtot",
        "EPtot",
        "BOND",
        "ANGLE",
        "DIHED",
        "UB",
        "IMP",
        "VDWAALS",
        "EELEC",
        "EGB",
        "EHBOND",
        "RESTRAINT",
        "1-4 VDW",
        "1-4 EEL",
        "ENERGY",
    }
)
_AMBER_TEMP_KEYS = frozenset({"TEMP", "TEMP(K)"})
_AMBER_PRESS_KEYS = frozenset({"PRESS", "PRESSURE"})


def _amber_unit_type(key: str) -> str:
    if key in _AMBER_TEMP_KEYS or key.upper().startswith("TEMP"):
        return "temperature"
    if key in _AMBER_PRESS_KEYS or key.upper().startswith("PRESS"):
        return "pressure"
    if key in _AMBER_ENERGY_KEYS or key in {"Density", "VOLUME", "Volume"}:
        if key.upper() in {"DENSITY", "VOLUME"}:
            return "other"
        return "energy"
    return "other"


def _parse_float_token(token: str) -> Optional[float]:
    if not token or "*" in token:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def _parse_amber_md_energy_frames(content: str) -> List[Dict[str, float]]:
    """Extract MD energy prints (NSTEP = … TIME(PS) = …)."""
    frames: List[Dict[str, float]] = []
    parts = re.split(r"(?=^\s*NSTEP\s*=)", content, flags=re.MULTILINE)
    for part in parts:
        if not re.search(r"NSTEP\s*=", part):
            continue
        frame: Dict[str, float] = {}
        m = _NSTEP_RE.search(part)
        if not m:
            continue
        frame["NSTEP"] = float(m.group(1))
        frame["TIME"] = float(m.group(2))
        if m.group(3) is not None:
            frame["TEMP"] = float(m.group(3))
        if m.group(4) is not None:
            frame["PRESS"] = float(m.group(4))
        for line in part.splitlines()[1:]:
            if re.match(r"^-{5,}", line) or line.strip().startswith("|"):
                break
            for key, val in _ENERGY_LINE_RE.findall(line):
                parsed = _parse_float_token(val)
                if parsed is not None:
                    frame[key] = parsed
        if "EPtot" in frame or "Etot" in frame or "ENERGY" in frame:
            frames.append(frame)
    return frames


def _parse_amber_minimization_frames(content: str) -> List[Dict[str, float]]:
    """Extract minimization ENERGY table prints (no TIME(PS) column).

    Amber min mdouts look like::

        NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
           100      -1.0420E+05     4.0344E+00     ...

         BOND    = ...  ANGLE   = ...  DIHED      = ...
    """
    if not _MIN_HEADER_RE.search(content):
        return []

    frames: List[Dict[str, float]] = []
    parts = re.split(
        r"(?=^\s*NSTEP\s+ENERGY\s+RMS\s+GMAX)",
        content,
        flags=re.MULTILINE | re.IGNORECASE,
    )
    for part in parts:
        if not _MIN_HEADER_RE.search(part):
            continue
        lines = part.splitlines()
        data_line = None
        data_idx = -1
        for i, line in enumerate(lines):
            if i == 0:
                continue
            if not line.strip() or line.strip().startswith("|"):
                continue
            if _MIN_HEADER_RE.match(line):
                continue
            data_line = line
            data_idx = i
            break
        if data_line is None:
            continue
        m = _MIN_ROW_RE.match(data_line)
        if not m:
            continue
        energy = _parse_float_token(m.group(2))
        if energy is None:
            continue
        frame: Dict[str, float] = {
            "NSTEP": float(m.group(1)),
            "ENERGY": energy,
            "Etot": energy,  # alias so min shows up with MD Etot plots
        }
        for line in lines[data_idx + 1 :]:
            if not line.strip():
                continue
            if _MIN_HEADER_RE.match(line) or re.match(r"^-{5,}", line):
                break
            if line.strip().startswith("|"):
                break
            if re.match(r"^\s*\d+\s+", line) and "ENERGY" not in line:
                # Next numeric row without repeating header — stop
                if _MIN_ROW_RE.match(line):
                    break
            for key, val in _ENERGY_LINE_RE.findall(line):
                parsed = _parse_float_token(val)
                if parsed is not None:
                    frame[key] = parsed
        frames.append(frame)
    return frames


def _parse_amber_energy_frames(content: str) -> List[Dict[str, float]]:
    """Extract per-print energy dictionaries from Amber mdout text.

    Prefers MD ``TIME(PS)`` blocks; falls back to minimization ENERGY tables.
    """
    md_frames = _parse_amber_md_energy_frames(content)
    if md_frames:
        return md_frames
    return _parse_amber_minimization_frames(content)


class AmberLogEnergyAnalyzer:
    """Parse Amber mdout files into column-oriented energy series.

    ``file_times`` maps basename → duration in **ns** (same semantics as
    NAMD/GROMACS). When set, points from that file are spaced evenly across
    the assigned duration instead of using the mdout ``TIME(PS)`` column.
    """

    def __init__(
        self,
        log_files: List[Union[str, Path]],
        file_times: Optional[Dict[str, float]] = None,
        file_strides: Optional[Dict[str, int]] = None,
    ):
        self.log_files = [Path(f) for f in log_files]
        self.file_times = file_times or {}
        self.file_strides = file_strides or {}
        self.data: Dict[str, List[float]] = {}
        self._time_ns: List[float] = []
        self._file_ranges: Dict[str, tuple] = {}
        self._parse()

    def _lookup_file_time_ns(self, path: Path) -> Optional[float]:
        for key in (path.name, str(path)):
            if key in self.file_times:
                try:
                    value = float(self.file_times[key])
                except (TypeError, ValueError):
                    return None
                return value if value > 0 else None
        return None

    def _parse(self) -> None:
        import numpy as np

        cumulative_ns = 0.0
        for path in self.log_files:
            if not path.is_file():
                continue
            try:
                content = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            frames = _parse_amber_energy_frames(content)
            if not frames:
                continue

            stride = max(1, int(lookup_file_map(self.file_strides, path) or 1))
            if stride > 1:
                frames = frames[::stride]
            if not frames:
                continue

            start_idx = len(self._time_ns)
            n_points = len(frames)
            assigned_ns = self._lookup_file_time_ns(path)

            if assigned_ns is not None:
                # Match GROMACS/NAMD: redistribute points across user duration
                if n_points == 1:
                    local_times = [0.0]
                else:
                    local_times = np.linspace(
                        0.0, assigned_ns, n_points, endpoint=True
                    ).tolist()
                file_span_ns = assigned_ns
            elif "TIME" in frames[0]:
                local_t0 = float(frames[0].get("TIME", 0.0))
                local_times = [
                    (float(fr.get("TIME", 0.0)) - local_t0) / 1000.0  # ps → ns
                    for fr in frames
                ]
                file_span_ns = (
                    max(local_times[-1] - local_times[0], 0.0) if local_times else 0.0
                )
            else:
                # Minimization: no TIME(PS) — space by relative NSTEP on [0, 1] ns
                n0 = float(frames[0].get("NSTEP", 0.0))
                n1 = float(frames[-1].get("NSTEP", n0))
                span = max(n1 - n0, 1.0)
                local_times = [
                    (float(fr.get("NSTEP", n0)) - n0) / span for fr in frames
                ]
                file_span_ns = 1.0

            for fr, t_local in zip(frames, local_times):
                self._time_ns.append(float(t_local) + cumulative_ns)
                for key, val in fr.items():
                    if key in {"NSTEP", "TIME"}:
                        continue
                    self.data.setdefault(key, []).append(val)
                n = len(self._time_ns)
                for key, series in self.data.items():
                    while len(series) < n:
                        series.append(float("nan"))

            end_idx = len(self._time_ns)
            self._file_ranges[str(path)] = (start_idx, end_idx)
            cumulative_ns += file_span_ns

    def get_available_properties(self) -> List[str]:
        return sorted(k for k, v in self.data.items() if any(x == x for x in v))

    def _calculate_time_array(self):
        import numpy as np

        if not self._time_ns:
            return np.array([], dtype=float)
        return np.array(self._time_ns, dtype=float)


def list_amber_energy_properties(
    log_files: List[Union[str, Path]],
    file_times: Optional[Dict[str, float]] = None,
) -> List[str]:
    """Return available Amber mdout properties (fast header/frame peek)."""
    del file_times
    logs = [Path(f) for f in log_files]
    for path in logs:
        if not path.is_file():
            continue
        try:
            # Read a bounded prefix — enough for the first energy frame(s)
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                content = fh.read(512_000)
        except OSError:
            continue
        frames = _parse_amber_energy_frames(content)
        if frames:
            keys = sorted(
                k
                for k in frames[0].keys()
                if k not in {"NSTEP", "TIME"} and frames[0].get(k) == frames[0].get(k)
            )
            if keys:
                return keys
    if logs:
        analyzer = AmberLogEnergyAnalyzer([logs[0]])
        return analyzer.get_available_properties()
    return []


def run_amber_energetic_analysis(
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
    """Return JSON-serializable energetic series from Amber mdout file(s)."""
    import numpy as np

    analyzer = AmberLogEnergyAnalyzer(
        log_files, file_times=file_times, file_strides=file_strides
    )
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

        unit_type = _amber_unit_type(label)
        converted = arr.copy()
        unit_label = ""
        if unit_type == "energy":
            # Amber reports kcal/mol natively
            if energy_units.lower().startswith("kj"):
                converted = arr * 4.184
                unit_label = "kJ/mol"
            else:
                unit_label = "kcal/mol"
        elif unit_type == "temperature":
            unit_label = temperature_units or "K"
        elif unit_type == "pressure":
            # Amber PRESS is in bar; convert to atm if requested
            if pressure_units.lower().startswith("atm"):
                converted = arr / 1.01325
                unit_label = "atm"
            else:
                unit_label = "bar"
        else:
            unit_label = volume_units if "vol" in label.lower() else ""

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


def plot_amber_properties(
    log_files: List[Union[str, Path]],
    properties: Optional[List[str]] = None,
    file_times: Optional[Dict[str, float]] = None,
    *,
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
) -> None:
    """Plot Amber mdout energetic series using the shared PlotSpec renderer."""
    from gatewizard.utils import matplotlib_renderer
    from gatewizard.utils.plot_spec import plot_spec_from_plot_properties_kwargs

    data = run_amber_energetic_analysis(
        log_files,
        properties=properties,
        file_times=file_times,
        time_units=time_units,
        energy_units=energy_units,
        pressure_units=pressure_units,
        temperature_units=temperature_units,
        volume_units=volume_units,
    )
    if not data.get("series"):
        logger.warning("No Amber energetic data to plot")
        return

    plot_spec = plot_spec_from_plot_properties_kwargs(
        [s["name"] for s in data["series"]],
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
    for i, panel in enumerate(plot_spec["panels"]):
        if i < len(data["series"]):
            panel["key"] = data["series"][i]["key"]

    if separate_plots:
        prefix = save_prefix or "amber_"
        for i, panel in enumerate(plot_spec["panels"]):
            single = plot_spec_from_plot_properties_kwargs(
                [panel.get("name") or panel["key"]],
                separate_plots=False,
                line_colors=[panel.get("line_color")],
                energy_units=energy_units,
                time_units=time_units,
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
            single["panels"][0]["key"] = panel["key"]
            import matplotlib.pyplot as plt

            fig = matplotlib_renderer.render_energetic(data, single)
            safe = str(panel.get("name") or panel["key"]).lower().replace(" ", "_")
            out = f"{prefix}{safe}.png"
            fig.savefig(out, dpi=dpi, bbox_inches="tight")
            plt.close(fig)
            logger.info("Plot saved: %s", out)
            if show:
                plt.show()
    else:
        if len(data["series"]) > 1:
            plot_spec["layout"] = "grid"
            plot_spec["cols"] = min(len(data["series"]), 2)
        import matplotlib.pyplot as plt

        fig = matplotlib_renderer.render_energetic(data, plot_spec)
        if save:
            fig.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info("Plot saved: %s", save)
        if show:
            plt.show()
        plt.close(fig)
