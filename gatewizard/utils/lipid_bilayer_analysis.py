"""
Lipid bilayer trajectory analysis using lipyphilic.

Provides area-per-lipid and membrane-thickness analysis for MD trajectories
via Voronoi tessellation and interleaflet headgroup distances (lipyphilic).
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

_LIPYPHILIC_INSTALL = "pip install lipyphilic  # or: pip install -e . from the gatewizard repo"


def _require_lipyphilic():
    """Import lipyphilic or raise a clear installation error."""
    try:
        import lipyphilic  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "lipyphilic is required for lipid bilayer analysis. "
            f"Install with: {_LIPYPHILIC_INSTALL}"
        ) from exc


def _to_path_list(paths: List[Union[str, Path]]) -> List[Path]:
    return [Path(p).expanduser().resolve() for p in paths]


def _analysis_result(analysis: Any, attribute: str) -> Any:
    """Read lipyphilic analysis output from ``results`` or legacy attributes."""
    if hasattr(analysis, "results") and hasattr(analysis.results, attribute):
        return getattr(analysis.results, attribute)
    return getattr(analysis, attribute)


def _leaflet_means_per_frame(
    areas: "np.ndarray", leaflets: "np.ndarray"
) -> Dict[str, "np.ndarray"]:
    """Compute mean area per lipid for each leaflet at each frame."""
    import numpy as np

    n_lipids, n_frames = areas.shape
    upper = np.full(n_frames, np.nan)
    lower = np.full(n_frames, np.nan)

    if leaflets.ndim == 1:
        upper_mask = leaflets == 1
        lower_mask = leaflets == -1
        if np.any(upper_mask):
            upper[:] = np.nanmean(areas[upper_mask, :], axis=0)
        if np.any(lower_mask):
            lower[:] = np.nanmean(areas[lower_mask, :], axis=0)
        return {"upper": upper, "lower": lower}

    for frame in range(n_frames):
        frame_leaflets = leaflets[:, frame]
        upper_mask = frame_leaflets == 1
        lower_mask = frame_leaflets == -1
        if np.any(upper_mask):
            upper[frame] = np.nanmean(areas[upper_mask, frame])
        if np.any(lower_mask):
            lower[frame] = np.nanmean(areas[lower_mask, frame])

    return {"upper": upper, "lower": lower}


def _stats_from_series(values: "np.ndarray") -> Dict[str, float]:
    import numpy as np

    clean = np.asarray(values, dtype=float)
    clean = clean[np.isfinite(clean)]
    if clean.size == 0:
        return {"mean": float("nan"), "std": float("nan"), "min": float("nan"), "max": float("nan")}
    return {
        "mean": float(np.mean(clean)),
        "std": float(np.std(clean)),
        "min": float(np.min(clean)),
        "max": float(np.max(clean)),
    }


def _correct_pbc_straddling_thickness(
    thickness: "np.ndarray", box_z: float
) -> "np.ndarray":
    """If thickness is the long periodic path through water, fold to the bilayer gap.

    ``lipyphilic.MembThickness`` wraps coordinates then does
    ``z_upper - z_lower``. When the bilayer straddles the periodic z boundary,
    that difference becomes ``L_z - d`` (water thickness) instead of ``d``.
    """
    import numpy as np

    values = np.abs(np.asarray(thickness, dtype=float))
    if not np.isfinite(box_z) or box_z <= 0:
        return values
    straddling = values > (0.5 * box_z)
    if np.any(straddling):
        values = values.copy()
        values[straddling] = box_z - values[straddling]
    return values


class BilayerTrajectoryAnalyzer:
    """
    Lipid bilayer analysis wrapper built on lipyphilic and MDAnalysis.

    Supports:
    - Area per lipid (Voronoi tessellation)
    - Membrane thickness (interleaflet headgroup distance)

    Example:
        >>> analyzer = BilayerTrajectoryAnalyzer("bilayer.pdb", "traj.xtc")
        >>> data = analyzer.calculate_area_per_lipid(lipid_sel="name GL1 GL2")
        >>> thickness = analyzer.calculate_membrane_thickness(lipid_sel="name PO4")
    """

    def __init__(
        self,
        topology: Union[str, Path],
        trajectory: Union[str, Path, List[Union[str, Path]]],
        file_times: Optional[Dict[str, float]] = None,
    ):
        from gatewizard.utils.namd_analysis import TrajectoryAnalyzer

        self._trajectory = TrajectoryAnalyzer(topology, trajectory, file_times=file_times)
        self._z_centered_for: Optional[str] = None

    @property
    def universe(self):
        return self._trajectory.universe

    def _calculate_time_array(self):
        return self._trajectory._calculate_time_array()

    def _ensure_membrane_centered_in_z(self, lipid_sel: str) -> None:
        """Shift the bilayer so it does not straddle the periodic z boundary.

        ``lipyphilic.MembThickness`` wraps atoms into the primary cell before
        subtracting leaflet heights. If the membrane sits across z=0, wrapping
        moves one leaflet to the top of the box and the reported "thickness"
        becomes the water gap (``L_z - d``). Centering the lipid COM at the box
        midplane without wrapping keeps both leaflets contiguous.
        """
        if self._z_centered_for == lipid_sel:
            return
        if self._z_centered_for is not None and self._z_centered_for != lipid_sel:
            logger.debug(
                "Membrane already z-centered for %r; skipping re-center for %r",
                self._z_centered_for,
                lipid_sel,
            )
            return

        try:
            from MDAnalysis.transformations import center_in_box
        except ImportError:
            return

        ag = self.universe.select_atoms(lipid_sel)
        if len(ag) == 0:
            return

        try:
            self.universe.trajectory.add_transformations(
                center_in_box(ag, center="mass", wrap=False)
            )
            self._z_centered_for = lipid_sel
            logger.debug(
                "Centered bilayer in z using %d atoms from %r", len(ag), lipid_sel
            )
        except ValueError as exc:
            # Transformations already locked (trajectory previously iterated).
            logger.warning(
                "Could not add membrane z-centering transformation (%s). "
                "Thickness may be wrong if the bilayer straddles the periodic boundary.",
                exc,
            )

    def _assign_leaflets(
        self,
        lipid_sel: str,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        verbose: bool = False,
    ):
        _require_lipyphilic()
        from lipyphilic.leaflets.assign_leaflets import AssignLeaflets

        leaflets = AssignLeaflets(universe=self.universe, lipid_sel=lipid_sel)
        leaflets.run(start=start, stop=stop, step=step, verbose=verbose)
        return leaflets

    def _lipid_residue_metadata(self, lipid_sel: str) -> Dict[str, List]:
        """Return residue IDs/names for lipids matching the selection."""
        atoms = self.universe.select_atoms(lipid_sel)
        residues = atoms.residues
        return {
            "resids": residues.resids.tolist(),
            "resnames": list(residues.resnames),
        }

    def calculate_area_per_lipid(
        self,
        lipid_sel: str = "name PO4",
        leaflet_lipid_sel: Optional[str] = None,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """
        Calculate area per lipid via 2D Voronoi tessellation (lipyphilic).

        Args:
            lipid_sel: Atom selection for Voronoi tessellation (e.g. MARTINI
                ``name GL1 GL2 ROH`` or all-atom ``name PO4``).
            leaflet_lipid_sel: Selection for leaflet assignment. Defaults to
                ``lipid_sel``.
            start, stop, step: Trajectory frame range passed to lipyphilic.
            verbose: Show lipyphilic progress bars.

        Returns:
            Dict with time (ns), per-lipid areas, leaflet means, and statistics.

        Note:
            Only lipid sites enter the Voronoi tessellation. In protein–membrane
            systems the protein footprint is not excluded, so mean APL can be
            close to ``L_x * L_y / (n_lipids / 2)`` (upper bound).
        """
        import numpy as np

        _require_lipyphilic()
        from lipyphilic.analysis.area_per_lipid import AreaPerLipid

        leaflet_sel = leaflet_lipid_sel or lipid_sel
        self._ensure_membrane_centered_in_z(leaflet_sel)
        leaflets = self._assign_leaflets(
            leaflet_sel, start=start, stop=stop, step=step, verbose=verbose
        )

        areas = AreaPerLipid(
            universe=self.universe,
            lipid_sel=lipid_sel,
            leaflets=leaflets.leaflets,
        )
        areas.run(start=start, stop=stop, step=step, verbose=verbose)

        area_array = np.asarray(_analysis_result(areas, "areas"), dtype=float)
        leaflet_data = _analysis_result(leaflets, "leaflets")
        leaflet_means = _leaflet_means_per_frame(area_array, leaflet_data)
        mean_per_frame = np.nanmean(area_array, axis=0)

        n_frames = area_array.shape[1]
        time_ns = self._calculate_time_array()
        if len(time_ns) != n_frames:
            time_ns = time_ns[:n_frames] if len(time_ns) > n_frames else np.linspace(
                0.0, max(len(time_ns) - 1, 0) * 0.01, n_frames
            )

        metadata = self._lipid_residue_metadata(lipid_sel)

        return {
            "time": time_ns,
            "areas": area_array,
            "mean_area_per_lipid": mean_per_frame,
            "mean_upper_leaflet": leaflet_means["upper"],
            "mean_lower_leaflet": leaflet_means["lower"],
            "resids": metadata["resids"],
            "resnames": metadata["resnames"],
        }

    def calculate_membrane_thickness(
        self,
        lipid_sel: str = "name PO4",
        leaflet_lipid_sel: Optional[str] = None,
        leaflet_filter_sel: Optional[str] = None,
        n_bins: int = 1,
        interpolate: bool = False,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """
        Calculate bilayer thickness from interleaflet headgroup distances.

        Args:
            lipid_sel: Headgroup atom selection for the thickness calculation.
            leaflet_lipid_sel: Selection for leaflet assignment. Defaults to
                ``lipid_sel``.
            leaflet_filter_sel: Optional selection passed to
                ``AssignLeaflets.filter_leaflets()`` to exclude species (e.g.
                cholesterol) from the thickness calculation.
            n_bins: Grid resolution for intrinsic surface construction.
            interpolate: Interpolate missing grid values (slower).
            start, stop, step: Trajectory frame range passed to lipyphilic.
            verbose: Show lipyphilic progress bars.

        Returns:
            Dict with time (ns), thickness (Å), and statistics.
        """
        import numpy as np

        _require_lipyphilic()
        from lipyphilic.analysis.memb_thickness import MembThickness

        leaflet_sel = leaflet_lipid_sel or lipid_sel
        self._ensure_membrane_centered_in_z(leaflet_sel)
        leaflets = self._assign_leaflets(
            leaflet_sel, start=start, stop=stop, step=step, verbose=verbose
        )

        if leaflet_filter_sel:
            leaflet_data = leaflets.filter_leaflets(leaflet_filter_sel)
        else:
            leaflet_data = _analysis_result(leaflets, "leaflets")

        memb_thickness = MembThickness(
            universe=self.universe,
            lipid_sel=lipid_sel,
            leaflets=leaflet_data,
            n_bins=n_bins,
            interpolate=interpolate,
        )
        memb_thickness.run(start=start, stop=stop, step=step, verbose=verbose)

        thickness = np.asarray(
            _analysis_result(memb_thickness, "memb_thickness"), dtype=float
        ).ravel()
        box_z = float(self.universe.dimensions[2]) if self.universe.dimensions is not None else 0.0
        thickness = _correct_pbc_straddling_thickness(thickness, box_z)

        n_frames = thickness.size
        time_ns = self._calculate_time_array()
        if len(time_ns) != n_frames:
            time_ns = time_ns[:n_frames] if len(time_ns) > n_frames else np.linspace(
                0.0, max(len(time_ns) - 1, 0) * 0.01, n_frames
            )

        return {
            "time": time_ns,
            "thickness": thickness,
        }

    def plot_area_per_lipid(
        self,
        lipid_sel: str = "name PO4",
        leaflet_lipid_sel: Optional[str] = None,
        series: str = "mean",
        time_units: str = "ns",
        area_units: str = "Å²",
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        verbose: bool = False,
        line_color: str = "#1f77b4",
        bg_color: str = "#ffffff",
        fig_bg_color: str = "#ffffff",
        text_color: str = "black",
        show_grid: bool = True,
        title: Optional[str] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (10, 6),
        dpi: int = 300,
    ):
        """Plot area-per-lipid time series."""
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return

        data = self.calculate_area_per_lipid(
            lipid_sel=lipid_sel,
            leaflet_lipid_sel=leaflet_lipid_sel,
            start=start,
            stop=stop,
            step=step,
            verbose=verbose,
        )

        series_key = {
            "mean": "mean_area_per_lipid",
            "upper": "mean_upper_leaflet",
            "lower": "mean_lower_leaflet",
        }.get(series, "mean_area_per_lipid")

        plot_time = np.asarray(data["time"], dtype=float)
        if time_units == "ps":
            plot_time = plot_time * 1000.0
            xlabel = "Time (ps)"
        elif time_units in {"us", "µs"}:
            plot_time = plot_time / 1000.0
            xlabel = "Time (µs)"
        else:
            xlabel = "Time (ns)"

        y = np.asarray(data[series_key], dtype=float)

        fig, ax = plt.subplots(figsize=figsize)
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        if bg_color != "none":
            ax.set_facecolor(bg_color)

        ax.plot(plot_time, y, color=line_color, linewidth=1.5)
        ax.set_xlabel(xlabel, color=text_color)
        ax.set_ylabel(f"Area per lipid ({area_units})", color=text_color)
        ax.set_title(
            title or f"Area per lipid ({series} leaflet)" if series != "mean" else "Area per lipid",
            color=text_color,
        )
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        if show_grid:
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")
        if show:
            plt.show()
        else:
            plt.close()

    def plot_membrane_thickness(
        self,
        lipid_sel: str = "name PO4",
        leaflet_lipid_sel: Optional[str] = None,
        leaflet_filter_sel: Optional[str] = None,
        n_bins: int = 1,
        interpolate: bool = False,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        verbose: bool = False,
        time_units: str = "ns",
        thickness_units: str = "Å",
        line_color: str = "#2ca02c",
        bg_color: str = "#ffffff",
        fig_bg_color: str = "#ffffff",
        text_color: str = "black",
        show_grid: bool = True,
        title: Optional[str] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (10, 6),
        dpi: int = 300,
    ):
        """Plot membrane thickness time series."""
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return

        data = self.calculate_membrane_thickness(
            lipid_sel=lipid_sel,
            leaflet_lipid_sel=leaflet_lipid_sel,
            leaflet_filter_sel=leaflet_filter_sel,
            n_bins=n_bins,
            interpolate=interpolate,
            start=start,
            stop=stop,
            step=step,
            verbose=verbose,
        )

        plot_time = np.asarray(data["time"], dtype=float)
        if time_units == "ps":
            plot_time = plot_time * 1000.0
            xlabel = "Time (ps)"
        elif time_units in {"us", "µs"}:
            plot_time = plot_time / 1000.0
            xlabel = "Time (µs)"
        else:
            xlabel = "Time (ns)"

        y = np.asarray(data["thickness"], dtype=float)

        fig, ax = plt.subplots(figsize=figsize)
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        if bg_color != "none":
            ax.set_facecolor(bg_color)

        ax.plot(plot_time, y, color=line_color, linewidth=1.5)
        ax.set_xlabel(xlabel, color=text_color)
        ax.set_ylabel(f"Membrane thickness ({thickness_units})", color=text_color)
        ax.set_title(title or "Membrane thickness", color=text_color)
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        if show_grid:
            ax.grid(True, alpha=0.3)

        plt.tight_layout()
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")
        if show:
            plt.show()
        else:
            plt.close()


def run_bilayer_analysis(
    topology_file: Union[str, Path],
    trajectory_files: List[Union[str, Path]],
    analysis_type: str,
    lipid_sel: str = "name PO4",
    leaflet_lipid_sel: Optional[str] = None,
    leaflet_filter_sel: Optional[str] = None,
    n_bins: int = 1,
    interpolate: bool = False,
    file_times: Optional[Dict[str, float]] = None,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Run lipid bilayer analysis and return JSON-serializable arrays.

    Supported analysis types: ``area_per_lipid``, ``membrane_thickness``.
    """
    import numpy as np

    top = Path(topology_file).expanduser().resolve()
    trajs = _to_path_list(trajectory_files)
    analyzer = BilayerTrajectoryAnalyzer(top, trajs, file_times=file_times)

    atype = analysis_type.strip().lower().replace(" ", "_").replace("-", "_")
    if atype in {"area_per_lipid", "apl"}:
        data = analyzer.calculate_area_per_lipid(
            lipid_sel=lipid_sel,
            leaflet_lipid_sel=leaflet_lipid_sel,
            start=start,
            stop=stop,
            step=step,
            verbose=verbose,
        )
        mean_y = np.asarray(data["mean_area_per_lipid"], dtype=float)
        return {
            "analysis_type": "area_per_lipid",
            "x": np.asarray(data["time"], dtype=float).tolist(),
            "y": mean_y.tolist(),
            "x_label": "Time (ns)",
            "y_label": "Area per lipid (Å²)",
            "series_name": "Mean area per lipid",
            "mean_upper_leaflet": np.asarray(
                data["mean_upper_leaflet"], dtype=float
            ).tolist(),
            "mean_lower_leaflet": np.asarray(
                data["mean_lower_leaflet"], dtype=float
            ).tolist(),
            "lipid_resids": data["resids"],
            "lipid_resnames": data["resnames"],
            "per_lipid_areas": np.asarray(data["areas"], dtype=float).tolist(),
            "stats": _stats_from_series(mean_y),
        }

    if atype in {"membrane_thickness", "memb_thickness", "thickness"}:
        data = analyzer.calculate_membrane_thickness(
            lipid_sel=lipid_sel,
            leaflet_lipid_sel=leaflet_lipid_sel,
            leaflet_filter_sel=leaflet_filter_sel,
            n_bins=n_bins,
            interpolate=interpolate,
            start=start,
            stop=stop,
            step=step,
            verbose=verbose,
        )
        y = np.asarray(data["thickness"], dtype=float)
        return {
            "analysis_type": "membrane_thickness",
            "x": np.asarray(data["time"], dtype=float).tolist(),
            "y": y.tolist(),
            "x_label": "Time (ns)",
            "y_label": "Membrane thickness (Å)",
            "series_name": "Membrane thickness",
            "n_bins": n_bins,
            "stats": _stats_from_series(y),
        }

    raise ValueError(
        f"Unsupported bilayer analysis type: {analysis_type}. "
        "Supported: area_per_lipid, membrane_thickness"
    )
