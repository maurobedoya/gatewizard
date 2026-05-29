# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
OpenMM log file analysis — mirrors the EnergyAnalyzer (NAMD) interface.

OpenMM writes energy data via StateDataReporter to stdout in a tab-separated
format with a '#'-prefixed header line.  This module parses those log files
and provides the same plotting API as the NAMD EnergyAnalyzer so analysis
workflows are interchangeable.

Energy values are stored internally in **kJ/mol** (OpenMM native units).
All plot methods accept ``energy_units='kJ/mol'`` (default) or
``energy_units='kcal/mol'`` for conversion.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, List, Optional, Union

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Column-name → internal key mapping
# ---------------------------------------------------------------------------
_COLUMN_MAP: Dict[str, str] = {
    "Step": "step",
    "Time (ps)": "time_ps",
    "Potential Energy (kJ/mole)": "potential",
    "Kinetic Energy (kJ/mole)": "kinetic",
    "Total Energy (kJ/mole)": "total",
    "Temperature (K)": "temp",
    "Box Volume (nm^3)": "volume",
    "Volume (nm^3)": "volume",
    "Density (g/mL)": "density",
    "Density (g/cm^3)": "density",
    # Older OpenMM column names
    "Potential Energy (kJ/mol)": "potential",
    "Kinetic Energy (kJ/mol)": "kinetic",
    "Total Energy (kJ/mol)": "total",
}

_NUMERIC_KEYS = {
    "step",
    "time_ps",
    "potential",
    "kinetic",
    "total",
    "temp",
    "volume",
    "density",
}


class OpenMMLogAnalyzer:
    """Parse and analyze OpenMM StateDataReporter log files.

    The interface mirrors :class:`~gatewizard.utils.namd_analysis.EnergyAnalyzer`
    so the two can be used interchangeably in analysis workflows.

    Energy values are stored in **kJ/mol**.  Pass ``energy_units='kcal/mol'``
    to any plot method to convert automatically.

    Args:
        log_file: Path to an OpenMM log file, or a list of paths for
            multi-stage analysis (times are concatenated).
        file_times: Optional dict ``{filename: duration_ns}`` that overrides
            the time axis derived from the log's "Time (ps)" column.  Useful
            when log files lack the time column.

    Example:
        >>> from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer
        >>> ana = OpenMMLogAnalyzer("step1_equilibration.log")
        >>> stats = ana.get_statistics()
        >>> print(stats["potential"]["mean"])
        >>> ana.plot_energy(save="energy.png")
    """

    def __init__(
        self,
        log_file: Union[Path, str, List[Union[Path, str]]],
        file_times: Optional[Dict[str, float]] = None,
    ):
        if isinstance(log_file, (str, Path)):
            self.log_files = [Path(log_file)]
        else:
            self.log_files = [Path(f) for f in log_file]

        self.file_times = file_times or {}
        self._file_ranges: Dict[str, tuple] = {}
        self.data = self._parse_log_data()

    # ------------------------------------------------------------------
    # Parsing
    # ------------------------------------------------------------------

    def _parse_log_data(self) -> Dict[str, List[float]]:
        """Parse tab-separated StateDataReporter output from all log files."""
        data: Dict[str, List] = {k: [] for k in _NUMERIC_KEYS}

        for log_file in self.log_files:
            if not log_file.exists():
                logger.warning(f"Log file not found: {log_file}")
                continue

            start_idx = len(data["step"])
            col_indices: Dict[str, int] = {}  # key → column index

            try:
                with open(log_file, "r", encoding="utf-8", errors="ignore") as fh:
                    for raw_line in fh:
                        line = raw_line.strip()
                        if not line:
                            continue

                        # Header line starts with '#'
                        if line.startswith("#"):
                            header = line.lstrip("#").strip()
                            parts = [p.strip().strip('"') for p in header.split("\t")]
                            col_indices = {}
                            for idx, col_name in enumerate(parts):
                                key = _COLUMN_MAP.get(col_name)
                                if key in _NUMERIC_KEYS:
                                    col_indices[key] = idx
                            continue

                        if not col_indices:
                            continue

                        parts = line.split("\t")
                        row: Dict[str, float] = {}
                        valid = True
                        for key, idx in col_indices.items():
                            if idx >= len(parts):
                                valid = False
                                break
                            try:
                                row[key] = float(parts[idx])
                            except ValueError:
                                valid = False
                                break
                        if not valid or "step" not in row:
                            continue

                        for key in _NUMERIC_KEYS:
                            if key in row:
                                data[key].append(row[key])
                            elif data[key] or key == "step":
                                # Keep arrays aligned: pad with NaN for missing cols
                                import math

                                data[key].append(math.nan)

            except Exception as exc:
                logger.error(f"Error parsing {log_file}: {exc}")
                continue

            end_idx = len(data["step"])
            first_step = int(data["step"][start_idx]) if start_idx < end_idx else 0
            last_step = int(data["step"][end_idx - 1]) if start_idx < end_idx else 0
            self._file_ranges[str(log_file)] = (
                start_idx,
                end_idx,
                first_step,
                last_step,
            )

        return data

    # ------------------------------------------------------------------
    # Statistics
    # ------------------------------------------------------------------

    def get_statistics(self) -> Dict[str, Dict[str, float]]:
        """Return mean/std/min/max/initial/final for every numeric column.

        Returns:
            ``{key: {"mean": …, "std": …, "min": …, "max": …,
                      "initial": …, "final": …}}``
        """
        import numpy as np

        stats: Dict[str, Dict[str, float]] = {}
        for key, values in self.data.items():
            if key == "step" or not values:
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

    # ------------------------------------------------------------------
    # Time axis
    # ------------------------------------------------------------------

    def _calculate_time_array(self) -> "np.ndarray":
        """Return time in nanoseconds for all data points."""
        import numpy as np

        if self.data["time_ps"] and not all(
            __import__("math").isnan(v) for v in self.data["time_ps"]
        ):
            # Direct from log: "Time (ps)" column
            time_ns = np.array(self.data["time_ps"], dtype=float) / 1000.0

            # If file_times provided, override with per-file durations
            if self.file_times and self._file_ranges:
                time_ns_adj = np.empty_like(time_ns)
                cumulative = 0.0
                for log_file in self.log_files:
                    key = str(log_file)
                    if key not in self._file_ranges:
                        continue
                    start_idx, end_idx, _, _ = self._file_ranges[key]
                    fname = log_file.name
                    duration = self.file_times.get(fname, None)
                    n_points = end_idx - start_idx
                    if duration is not None and n_points > 0:
                        time_ns_adj[start_idx:end_idx] = (
                            np.linspace(0, duration, n_points, endpoint=False)
                            + cumulative
                        )
                        cumulative += duration
                    else:
                        segment = time_ns[start_idx:end_idx]
                        if len(segment):
                            offset = cumulative - (
                                segment[0] if start_idx == 0 else segment[0]
                            )
                            time_ns_adj[start_idx:end_idx] = segment + offset
                            cumulative = float(time_ns_adj[end_idx - 1]) + (
                                (segment[-1] - segment[-2]) / 1000.0
                                if len(segment) > 1
                                else 0
                            )
                return time_ns_adj
            return time_ns

        # Fallback: use step numbers with a 2 fs timestep assumption
        steps = np.array(self.data["step"], dtype=float)
        return steps * 2e-6  # 2 fs per step → ns

    # ------------------------------------------------------------------
    # plot_properties  (mirrors EnergyAnalyzer.plot_properties)
    # ------------------------------------------------------------------

    def plot_properties(
        self,
        properties: Optional[List[str]] = None,
        energy_units: str = "kJ/mol",
        time_units: str = "ns",
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        show_grid: bool = True,
        separate_plots: bool = False,
        save: Optional[str] = None,
        show: bool = False,
        figsize: Optional[tuple] = None,
        dpi: int = 300,
        xlim: Optional[tuple] = None,
        ylim: Optional[tuple] = None,
        colors: Optional[List[str]] = None,
    ):
        """Plot selected properties vs time.

        Args:
            properties: List of property keys to plot, e.g.
                ``["potential", "kinetic", "temp", "volume"]``.
                Defaults to all non-empty numeric columns.
            energy_units: ``'kJ/mol'`` (default) or ``'kcal/mol'``.
            time_units: ``'ps'``, ``'ns'`` (default), or ``'µs'``.
            bg_color: Plot area background color.
            fig_bg_color: Figure border background color.
            text_color: Axis text color (``'Auto'`` or explicit color).
            show_grid: Draw grid lines.
            separate_plots: Save/show each property as its own figure.
            save: Filename prefix for saved plots.
            show: Display plots interactively.
            figsize: Figure size tuple.
            dpi: Dots per inch for saved figures.
            xlim: x-axis limits as ``(min, max)``.
            ylim: y-axis limits as ``(min, max)``.
            colors: Line color list (one per property).
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy required for plotting")
            return

        if not self.data["step"]:
            logger.warning("No data to plot")
            return

        if properties is None:
            properties = [
                k
                for k in _NUMERIC_KEYS
                if k not in ("step", "time_ps") and self.data[k]
            ]

        time_ns = self._calculate_time_array()
        plot_time, time_label = _convert_time(time_ns, time_units)
        energy_factor, energy_label = _energy_conversion(energy_units)
        energy_keys = {"potential", "kinetic", "total"}

        _tc = _auto_text_color(text_color, bg_color)
        default_colors = [
            "#61afef",
            "#98c379",
            "#e06c75",
            "#e5c07b",
            "#c678dd",
            "#56b6c2",
            "#d19a66",
            "#abb2bf",
        ]

        def _plot_one(ax, key, color):
            raw = np.array(self.data.get(key, []), dtype=float)
            if len(raw) == 0:
                return
            y = raw * energy_factor if key in energy_keys else raw
            y_label = _property_label(key, energy_label)
            ax.plot(plot_time[: len(y)], y, color=color, linewidth=0.8)
            ax.set_xlabel(f"Time ({time_label})", color=_tc)
            ax.set_ylabel(y_label, color=_tc)
            ax.set_title(key.replace("_", " ").title(), color=_tc)
            ax.tick_params(colors=_tc)
            for spine in ax.spines.values():
                spine.set_edgecolor(_tc)
            if bg_color != "none":
                ax.set_facecolor(bg_color)
            if show_grid:
                ax.grid(True, alpha=0.3, color=_tc, linewidth=0.5)
            if xlim:
                ax.set_xlim(xlim)
            if ylim:
                ax.set_ylim(ylim)

        if separate_plots:
            for prop_i, key in enumerate(properties):
                fig, ax = plt.subplots(figsize=figsize or (10, 4))
                if fig_bg_color != "none":
                    fig.patch.set_facecolor(fig_bg_color)
                color = (
                    colors[prop_i]
                    if colors and prop_i < len(colors)
                    else default_colors[prop_i % len(default_colors)]
                )
                _plot_one(ax, key, color)
                plt.tight_layout()
                if save:
                    fname = (
                        f"{save}_{key}.png"
                        if not save.endswith(".png")
                        else f"{key}_{save}"
                    )
                    plt.savefig(fname, dpi=dpi, bbox_inches="tight")
                if show:
                    plt.show()
                plt.close(fig)
        else:
            n = len(properties)
            cols = min(n, 2)
            rows = (n + cols - 1) // cols
            fig, axes_arr = plt.subplots(
                rows, cols, figsize=figsize or (12, 4 * rows), squeeze=False
            )
            if fig_bg_color != "none":
                fig.patch.set_facecolor(fig_bg_color)
            for prop_i, key in enumerate(properties):
                r, c = divmod(prop_i, cols)
                color = (
                    colors[prop_i]
                    if colors and prop_i < len(colors)
                    else default_colors[prop_i % len(default_colors)]
                )
                _plot_one(axes_arr[r][c], key, color)
            # Hide unused subplots
            for prop_i in range(len(properties), rows * cols):
                r, c = divmod(prop_i, cols)
                axes_arr[r][c].set_visible(False)
            plt.tight_layout()
            if save:
                plt.savefig(save, dpi=dpi, bbox_inches="tight")
            if show:
                plt.show()
            plt.close(fig)

    # ------------------------------------------------------------------
    # plot_energy  (mirrors EnergyAnalyzer.plot_energy — 2×2 summary)
    # ------------------------------------------------------------------

    def plot_energy(
        self,
        energy_units: str = "kJ/mol",
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
        """2×2 energy summary plot (total energy, potential+kinetic,
        temperature, volume/density).

        Args:
            energy_units: ``'kJ/mol'`` (default) or ``'kcal/mol'``.
            time_units: ``'ps'``, ``'ns'`` (default), or ``'µs'``.
            bg_color: Plot area background color.
            fig_bg_color: Figure border background color.
            text_color: Axis text color (``'Auto'`` or explicit color).
            show_grid: Draw grid lines.
            title: Figure title (auto-generated if None).
            target_temperature: Reference temperature in K.
            target_pressure: Reference pressure in bar (for annotation).
            save: Filename for saved figure.
            show: Display interactively.
            figsize: Figure dimensions.
            dpi: Resolution for saved figure.
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy required for plotting")
            return

        if not self.data["step"]:
            logger.warning("No energy data to plot")
            return

        time_ns = self._calculate_time_array()
        plot_time, time_label = _convert_time(time_ns, time_units)
        energy_factor, energy_label = _energy_conversion(energy_units)
        _tc = _auto_text_color(text_color, bg_color)

        fig, axes = plt.subplots(2, 2, figsize=figsize)
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)

        def _setup(ax, xlabel, ylabel, subplot_title):
            ax.set_xlabel(xlabel, color=_tc)
            ax.set_ylabel(ylabel, color=_tc)
            ax.set_title(subplot_title, color=_tc)
            ax.tick_params(colors=_tc)
            for spine in ax.spines.values():
                spine.set_edgecolor(_tc)
            if bg_color != "none":
                ax.set_facecolor(bg_color)
            if show_grid:
                ax.grid(True, alpha=0.3, color=_tc, linewidth=0.5)

        tl = f"Time ({time_label})"

        # ── top-left: Total energy ──────────────────────────────────────
        ax = axes[0, 0]
        if self.data["total"]:
            y = np.array(self.data["total"], dtype=float) * energy_factor
            ax.plot(
                plot_time[: len(y)], y, color="#61afef", linewidth=0.8, label="Total"
            )
        _setup(ax, tl, f"Total Energy ({energy_label})", "Total Energy")

        # ── top-right: Potential + Kinetic ──────────────────────────────
        ax = axes[0, 1]
        if self.data["potential"]:
            y = np.array(self.data["potential"], dtype=float) * energy_factor
            ax.plot(
                plot_time[: len(y)],
                y,
                color="#e06c75",
                linewidth=0.8,
                label="Potential",
            )
        if self.data["kinetic"]:
            y = np.array(self.data["kinetic"], dtype=float) * energy_factor
            ax.plot(
                plot_time[: len(y)], y, color="#98c379", linewidth=0.8, label="Kinetic"
            )
        ax.legend(
            fontsize=8,
            facecolor=bg_color if bg_color != "none" else "white",
            labelcolor=_tc,
        )
        _setup(ax, tl, f"Energy ({energy_label})", "Potential & Kinetic Energy")

        # ── bottom-left: Temperature ────────────────────────────────────
        ax = axes[1, 0]
        if self.data["temp"]:
            y = np.array(self.data["temp"], dtype=float)
            ax.plot(plot_time[: len(y)], y, color="#e5c07b", linewidth=0.8)
            if target_temperature is not None:
                ax.axhline(
                    target_temperature,
                    color="#abb2bf",
                    linestyle="--",
                    linewidth=0.8,
                    label=f"Target {target_temperature} K",
                )
            elif len(y) > 0:
                n50 = max(1, len(y) // 2)
                avg_t = float(np.nanmean(y[-n50:]))
                ax.axhline(
                    avg_t,
                    color="#abb2bf",
                    linestyle="--",
                    linewidth=0.8,
                    label=f"Avg {avg_t:.1f} K",
                )
            ax.legend(
                fontsize=8,
                facecolor=bg_color if bg_color != "none" else "white",
                labelcolor=_tc,
            )
        _setup(ax, tl, "Temperature (K)", "Temperature")

        # ── bottom-right: Volume or Density ────────────────────────────
        ax = axes[1, 1]
        if self.data["volume"] and not all(
            __import__("math").isnan(v) for v in self.data["volume"]
        ):
            y = np.array(self.data["volume"], dtype=float)
            ax.plot(plot_time[: len(y)], y, color="#c678dd", linewidth=0.8)
            _setup(ax, tl, "Volume (nm³)", "Box Volume")
        elif self.data["density"] and not all(
            __import__("math").isnan(v) for v in self.data["density"]
        ):
            y = np.array(self.data["density"], dtype=float)
            ax.plot(plot_time[: len(y)], y, color="#c678dd", linewidth=0.8)
            _setup(ax, tl, "Density (g/mL)", "Density")
        else:
            ax.text(
                0.5,
                0.5,
                "No volume/density data\n(NVT ensemble)",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color=_tc,
            )
            _setup(ax, tl, "", "Box Volume / Density")

        main_title = title or f"OpenMM Energy Analysis ({energy_units})"
        fig.suptitle(main_title, color=_tc, fontsize=13, fontweight="bold")
        plt.tight_layout()

        if save:
            plt.savefig(
                save,
                dpi=dpi,
                bbox_inches="tight",
                facecolor=fig_bg_color if fig_bg_color != "none" else "white",
            )
        if show:
            plt.show()
        plt.close(fig)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _convert_time(time_ns: "np.ndarray", units: str):
    if units == "ps":
        return time_ns * 1000.0, "ps"
    if units == "µs":
        return time_ns / 1000.0, "µs"
    return time_ns, "ns"


def _energy_conversion(energy_units: str):
    if energy_units == "kcal/mol":
        return 1.0 / 4.184, "kcal/mol"
    return 1.0, "kJ/mol"


def _auto_text_color(text_color: str, bg_color: str) -> str:
    if text_color != "Auto":
        return text_color
    if bg_color in ("none", "white", "#ffffff"):
        return "black"
    try:
        hx = bg_color.lstrip("#")
        r, g, b = int(hx[0:2], 16), int(hx[2:4], 16), int(hx[4:6], 16)
        luminance = 0.299 * r + 0.587 * g + 0.114 * b
        return "white" if luminance < 128 else "black"
    except Exception:
        return "white"


def _property_label(key: str, energy_label: str) -> str:
    energy_keys = {"potential", "kinetic", "total"}
    labels = {
        "potential": f"Potential Energy ({energy_label})",
        "kinetic": f"Kinetic Energy ({energy_label})",
        "total": f"Total Energy ({energy_label})",
        "temp": "Temperature (K)",
        "volume": "Volume (nm³)",
        "density": "Density (g/mL)",
        "time_ps": "Time (ps)",
    }
    return labels.get(key, key.replace("_", " ").title())


# ---------------------------------------------------------------------------
# Equilibration progress tracking (mirrors gromacs_analysis interface)
# ---------------------------------------------------------------------------

from dataclasses import dataclass


@dataclass
class OpenMMTimingInfo:
    """Progress/timing container — field names match NAMDTiming used by app.py."""

    steps_completed: int = 0
    total_steps: int = 0
    timestep_fs: float = 0.0  # femtoseconds per step
    ns_per_day: float = 0.0
    completed: bool = False
    has_error: bool = False


@dataclass
class OpenMMStageProgress:
    """Stage progress container — field names match NAMDProgress used by app.py."""

    stage_name: str = ""
    status: str = "not_started"  # not_started | running | completed | error
    timing: Optional[OpenMMTimingInfo] = None
    log_file: Optional[Path] = None


def _parse_openmm_inp(inp_file: Path) -> tuple[int, float]:
    """Return (nstep, dt_fs) from an OpenMM .inp parameter file."""
    nstep = 0
    dt_fs = 0.0
    if not inp_file.exists():
        return nstep, dt_fs
    try:
        content = inp_file.read_text(encoding="utf-8", errors="ignore")
        m = re.search(r"\bnstep\s*=\s*(\d+)", content)
        if m:
            nstep = int(m.group(1))
        m = re.search(r"\bdt\s*=\s*([\d.]+)", content)
        if m:
            dt_fs = float(m.group(1)) * 1000.0  # ps → fs
    except Exception:
        pass
    return nstep, dt_fs


def parse_openmm_log(
    log_file: Path, inp_file: Optional[Path] = None
) -> OpenMMTimingInfo:
    """
    Parse an OpenMM StateDataReporter log file and return timing/progress info.

    The log format (written by openmm_run.py via StateDataReporter) is:
        #"Progress (%)"\\t"Step"\\t"Time (ps)"\\t...\\t"Speed (ns/day)"\\t"Time Remaining"
        0.8%\\t1000\\t1.0\\t...\\t1.41\\t2:06:13

    ``nstep`` and ``dt`` are read from the companion ``.inp`` file when available.
    """
    info = OpenMMTimingInfo()

    if not log_file.exists():
        return info

    try:
        content = log_file.read_text(encoding="utf-8", errors="ignore")

        # ── Total steps & timestep from .inp file ────────────────────────────────
        if inp_file is None:
            inp_file = log_file.with_suffix(".inp")
        nstep, dt_fs = _parse_openmm_inp(inp_file)
        info.total_steps = nstep
        info.timestep_fs = dt_fs

        # ── Parse header to find column positions ────────────────────────────────
        # Header line starts with '#"Progress'
        header_m = re.search(r'^#"Progress[^\n]*', content, re.MULTILINE)
        if not header_m:
            return info

        header_cols = [
            c.strip().strip('"') for c in header_m.group().lstrip("#").split("\t")
        ]
        try:
            step_idx = header_cols.index("Step")
        except ValueError:
            step_idx = 1  # fallback: second column
        try:
            speed_idx = header_cols.index("Speed (ns/day)")
        except ValueError:
            speed_idx = len(header_cols) - 2  # second to last

        # ── Parse data rows ──────────────────────────────────────────────────────
        # Data rows start with a percentage like "0.8%"
        data_rows = re.findall(r"^\d+\.?\d*%\t[^\n]+", content, re.MULTILINE)
        if not data_rows:
            return info

        last_row = data_rows[-1].split("\t")

        try:
            info.steps_completed = int(last_row[step_idx])
        except (IndexError, ValueError):
            pass

        try:
            speed_val = float(last_row[speed_idx])
            if speed_val > 0:
                info.ns_per_day = speed_val
        except (IndexError, ValueError):
            pass

        # ── Completion / error markers ───────────────────────────────────────────
        if "Equilibration complete" in content or (
            info.total_steps > 0 and info.steps_completed >= info.total_steps
        ):
            info.completed = True

        if re.search(r"(Error|Traceback|failed)", content, re.IGNORECASE):
            # Avoid false positives from log messages that mention 'error' casually
            if re.search(r"^(Error|Traceback)", content, re.MULTILINE | re.IGNORECASE):
                info.has_error = True

    except Exception as exc:
        logger.debug(f"Error parsing OpenMM log {log_file}: {exc}")

    return info


def get_equilibration_progress(
    equilibration_dir: Path,
) -> Dict[str, OpenMMStageProgress]:
    """
    Return a progress dict for all standard OpenMM equilibration stages.

    Looks for ``step1_equilibration.log`` … ``step6_equilibration.log`` and
    ``step7_production.log`` directly in *equilibration_dir*.

    Returns:
        Ordered mapping of stage-name → :class:`OpenMMStageProgress`.
        Trailing ``not_started`` stages are trimmed so the GUI stays clean.
    """
    stage_log_map: Dict[str, str] = {
        "equilibration_1": "step1_equilibration",
        "equilibration_2": "step2_equilibration",
        "equilibration_3": "step3_equilibration",
        "equilibration_4": "step4_equilibration",
        "equilibration_5": "step5_equilibration",
        "equilibration_6": "step6_equilibration",
        "production": "step7_production",
    }

    progress: Dict[str, OpenMMStageProgress] = {}

    for stage_name, stem in stage_log_map.items():
        stage = OpenMMStageProgress(stage_name=stage_name)
        log_file = equilibration_dir / f"{stem}.log"
        inp_file = equilibration_dir / f"{stem}.inp"

        if log_file.exists():
            stage.log_file = log_file
            timing = parse_openmm_log(log_file, inp_file)
            stage.timing = timing

            if timing.has_error:
                stage.status = "error"
            elif timing.completed:
                stage.status = "completed"
            elif timing.steps_completed > 0:
                stage.status = "running"
            else:
                stage.status = "running"  # log exists but no output yet

        progress[stage_name] = stage

    # Trim trailing not_started stages
    keys = list(progress.keys())
    while keys and progress[keys[-1]].status == "not_started":
        del progress[keys[-1]]
        keys.pop()

    return progress
