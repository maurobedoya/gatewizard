"""
Engine-agnostic MD trajectory / structural analysis utilities.

Uses MDAnalysis for RMSD, RMSF, distance, radius of gyration, and related
helpers. Works with trajectories from NAMD, GROMACS, Amber, OpenMM, etc.
"""

from pathlib import Path
from typing import Dict, List, Optional, Any, TYPE_CHECKING, Union

from .logger import get_logger

if TYPE_CHECKING:
    import numpy as np

logger = get_logger(__name__)


def _to_path_list(paths: List[Union[str, Path]]) -> List[Path]:
    """Normalize a list of filesystem paths to resolved Path objects."""
    return [Path(p).expanduser().resolve() for p in paths]


_STRUCTURE_SNAPSHOT_SUFFIXES = frozenset({".pdb", ".ent", ".gro"})


def is_structure_snapshot(path: Union[str, Path]) -> bool:
    """True for static coordinate files (PDB/GRO), not DCD/XTC trajectories."""
    return Path(path).suffix.lower() in _STRUCTURE_SNAPSHOT_SUFFIXES


def split_analysis_trajectories(
    paths: List[Union[str, Path]],
) -> tuple[List[Path], List[Path]]:
    """Split ``paths`` into ``(snapshots, trajectories)``."""
    snapshots: List[Path] = []
    trajectories: List[Path] = []
    for raw in paths or []:
        p = Path(raw).expanduser()
        try:
            p = p.resolve()
        except OSError:
            pass
        if is_structure_snapshot(p):
            snapshots.append(p)
        else:
            trajectories.append(p)
    return snapshots, trajectories


def prepare_structural_inputs(
    trajectory_files: List[Union[str, Path]],
    *,
    analysis_type: str = "",
    reference_structure: Optional[Union[str, Path]] = None,
) -> tuple[List[Path], Optional[Path]]:
    """Drop static PDB/GRO snapshots from the traj chain when real trajs exist.

    Extra PDBs as frame 0 have no periodic box and crash membrane thickness / APL
    (``Box is None`` / ``NoneType`` subscript). For RMSD, a leftover snapshot or
    ``reference_structure`` is returned as the reference, not as a time-series frame.
    """
    snapshots, trajectories = split_analysis_trajectories(trajectory_files)
    coord = list(trajectories) if trajectories else list(snapshots)
    atype = str(analysis_type or "").strip().lower().replace(" ", "_").replace("-", "_")
    ref: Optional[Path] = None
    if reference_structure:
        ref = Path(reference_structure).expanduser()
        try:
            ref = ref.resolve()
        except OSError:
            pass
    elif atype in {"rmsd"} and snapshots and trajectories:
        ref = snapshots[0]
    if snapshots and trajectories:
        logger.warning(
            "Ignoring structure file(s) in the trajectory list (%s). "
            "PDB/GRO frames usually have no periodic box and break membrane "
            "thickness / APL. For RMSD vs a starting structure, set "
            "reference_structure instead of listing the PDB as a trajectory.",
            ", ".join(p.name for p in snapshots),
        )
    return coord, ref


def _fill_missing_box_dimensions(dimensions) -> None:
    """Forward/back-fill zero unit cells so MemoryReader frames keep a box."""
    n = int(getattr(dimensions, "shape", [0])[0] or 0)
    if n == 0:
        return
    last = None
    for i in range(n):
        if float(dimensions[i, 0]) > 0:
            last = dimensions[i].copy()
        elif last is not None:
            dimensions[i] = last
    nxt = None
    for i in range(n - 1, -1, -1):
        if float(dimensions[i, 0]) > 0:
            nxt = dimensions[i].copy()
        elif nxt is not None:
            dimensions[i] = nxt


def _lookup_file_map(file_map: Optional[Dict[str, Any]], path: Path) -> Any:
    """Look up a per-file value using basename, with case-insensitive fallback."""
    if not file_map:
        return None
    name = path.name
    if name in file_map:
        return file_map[name]
    name_lower = name.lower()
    for key, val in file_map.items():
        if key.lower() == name_lower:
            return val
    return None


def _align_time_to_frame_count(
    full_time: "np.ndarray",
    n_frames: int,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
) -> "np.ndarray":
    """
    Slice/resample a full time axis to match an analyzed frame count.

    Never substitutes a fake 0.01 ns/frame axis when ``full_time`` spans a
    positive duration — resamples within the assigned range instead.
    """
    import numpy as np

    if n_frames <= 0:
        return np.asarray([], dtype=float)

    full = np.asarray(full_time, dtype=float)
    if len(full) == 0:
        return np.linspace(0.0, max(n_frames - 1, 0) * 0.002, n_frames)

    s = 0 if start is None else max(0, int(start))
    st = 1 if step is None else max(1, int(step))
    stop_idx = len(full) if stop is None else min(int(stop), len(full))
    sliced = full[s:stop_idx:st]

    if len(sliced) == n_frames:
        return sliced.astype(float)
    if len(sliced) > n_frames:
        return sliced[:n_frames].astype(float)

    t_start = float(full[0])
    t_end = float(full[-1])
    if t_end > t_start:
        return np.linspace(t_start, t_end, n_frames)
    if len(sliced) >= 1:
        return np.linspace(float(sliced[0]), float(sliced[-1]), n_frames)
    return np.linspace(t_start, t_end if t_end > t_start else t_start, n_frames)


def run_structural_analysis(
    topology_file: Union[str, Path],
    trajectory_files: List[Union[str, Path]],
    analysis_type: str,
    selection: str = "protein and backbone",
    selection2: str = "",
    reference_frame: int = 0,
    align: bool = True,
    file_times: Optional[Dict[str, float]] = None,
    file_strides: Optional[Dict[str, int]] = None,
    rmsf_xaxis_type: str = "residue_number",
    reference_structure: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """
    Run trajectory structural analysis and return JSON-serializable arrays.

    Supported analysis types: `rmsd`, `rmsf`, `distance`, `radius_of_gyration`.

    ``reference_structure`` is an optional PDB/GRO used as the RMSD reference
    instead of ``reference_frame``. Static PDB files in ``trajectory_files`` are
    dropped when DCD/XTC files are also present (they have no periodic box).
    """
    import gc

    top = Path(topology_file).expanduser().resolve()
    trajs, ref_struct = prepare_structural_inputs(
        trajectory_files,
        analysis_type=analysis_type,
        reference_structure=reference_structure,
    )
    analyzer = TrajectoryAnalyzer(
        top, trajs, file_times=file_times, file_strides=file_strides
    )

    try:
        return _run_structural_analysis_body(
            analyzer,
            analysis_type,
            selection,
            selection2,
            reference_frame,
            align,
            rmsf_xaxis_type,
            reference_structure=ref_struct,
        )
    finally:
        analyzer.clear_analysis_cache()
        del analyzer
        gc.collect()


def _run_structural_analysis_body(
    analyzer: "TrajectoryAnalyzer",
    analysis_type: str,
    selection: str,
    selection2: str,
    reference_frame: int,
    align: bool,
    rmsf_xaxis_type: str,
    reference_structure: Optional[Union[str, Path]] = None,
) -> Dict[str, Any]:
    """Compute structural analysis arrays (caller owns analyzer lifecycle)."""
    import numpy as np

    atype = analysis_type.strip().lower().replace(" ", "_")
    if atype in {"rmsd"}:
        data = analyzer.calculate_rmsd(
            selection=selection,
            reference_frame=reference_frame,
            align=align,
            reference_structure=reference_structure,
        )
        y = np.asarray(data["rmsd"], dtype=float)
        return {
            "analysis_type": "rmsd",
            "x": np.asarray(data["time"], dtype=float).tolist(),
            "y": y.tolist(),
            "x_label": "Time (ns)",
            "y_label": "RMSD (Å)",
            "series_name": "RMSD",
            "stats": {
                "mean": float(np.mean(y)),
                "std": float(np.std(y)),
                "min": float(np.min(y)),
                "max": float(np.max(y)),
            },
        }

    if atype in {"rmsf"}:
        data = analyzer.calculate_rmsf(selection=selection)
        resids = np.asarray(data["resids"]).tolist()
        resnames = np.asarray(data.get("resnames", [])).tolist()
        atom_indices = list(range(len(resids)))

        if rmsf_xaxis_type == "residue_type_number":
            x_values = list(range(len(resids)))
            labels = []
            for i, rid in enumerate(resids):
                rname = resnames[i] if i < len(resnames) else ""
                labels.append(f"{rname}{rid}")
        elif rmsf_xaxis_type == "atom_index":
            x_values = atom_indices
            labels = [str(v) for v in atom_indices]
        else:  # residue_number (default)
            x_values = resids
            labels = [str(rid) for rid in resids]

        y = np.asarray(data["rmsf"], dtype=float)
        return {
            "analysis_type": "rmsf",
            "x": x_values,
            "x_labels": labels,
            "y": y.tolist(),
            "x_label": "Residue" if rmsf_xaxis_type != "atom_index" else "Atom index",
            "y_label": "RMSF (Å)",
            "series_name": "RMSF",
            "stats": {
                "mean": float(np.mean(y)),
                "std": float(np.std(y)),
                "min": float(np.min(y)),
                "max": float(np.max(y)),
            },
        }

    if atype in {"distance", "distances"}:
        if not selection2.strip():
            raise ValueError("selection2 is required for distance analysis")
        data = analyzer.calculate_distances({"distance": (selection, selection2)})[
            "distance"
        ]
        y = np.asarray(data["distance"], dtype=float)
        return {
            "analysis_type": "distance",
            "x": np.asarray(data["time"], dtype=float).tolist(),
            "y": y.tolist(),
            "x_label": "Time (ns)",
            "y_label": "Distance (Å)",
            "series_name": "Distance",
            "stats": {
                "mean": float(np.mean(y)),
                "std": float(np.std(y)),
                "min": float(np.min(y)),
                "max": float(np.max(y)),
            },
        }

    if atype in {"radius_of_gyration", "rg", "radius"}:
        data = analyzer.calculate_radius_of_gyration(selection=selection)
        y = np.asarray(data["rg"], dtype=float)
        return {
            "analysis_type": "radius_of_gyration",
            "x": np.asarray(data["time"], dtype=float).tolist(),
            "y": y.tolist(),
            "x_label": "Time (ns)",
            "y_label": "Radius of Gyration (Å)",
            "series_name": "Radius of Gyration",
            "stats": {
                "mean": float(np.mean(y)),
                "std": float(np.std(y)),
                "min": float(np.min(y)),
                "max": float(np.max(y)),
            },
        }

    raise ValueError(f"Unsupported structural analysis type: {analysis_type}")


class TrajectoryAnalyzer:
    """
    Easy-to-use wrapper for MD trajectory analysis with built-in plotting.

    Requires MDAnalysis to be installed.

    Supports all GUI options including:
    - Multiple trajectory files with time scaling
    - Full plot customization (colors, units, labels, limits)
    - RMSF-specific X-axis formatting with residue labels
    - Unit conversions (Å/nm, ps/ns/µs, kcal/kJ)

    Example:
        >>> analyzer = TrajectoryAnalyzer("system.psf", "trajectory.dcd")
        >>> analyzer.plot_rmsd(selection="protein and backbone", save="rmsd.png")
        >>>
        >>> # Multiple trajectories with time scaling
        >>> analyzer = TrajectoryAnalyzer(
        ...     "system.psf",
        ...     ["eq1.dcd", "eq2.dcd", "prod.dcd"],
        ...     file_times={"eq1.dcd": 1.0, "eq2.dcd": 2.0, "prod.dcd": 10.0}  # durations in ns
        ... )
        >>> analyzer.plot_rmsd(selection="protein", time_units="ns", distance_units="Å")
        >>>
        >>> # Customized plotting
        >>> analyzer.plot_rmsf(
        ...     selection="protein and name CA",
        ...     xaxis_type="residue_type_number",  # Show "ALA123" style labels
        ...     residue_name_format="triple",      # Use 3-letter codes
        ...     label_frequency="every_5",         # Label every 5th residue
        ...     line_color="#1f77b4",
        ...     bg_color="#2b2b2b",
        ...     show_grid=True,
        ...     save="rmsf.png"
        ... )
    """

    def __init__(
        self,
        topology: Path,
        trajectory: Path | List[Path],
        file_times: Optional[Dict[str, float]] = None,
        file_strides: Optional[Dict[str, int]] = None,
    ):
        """
        Initialize trajectory analyzer.

        Args:
            topology: Path to topology file (PSF, PDB, PRMTOP, etc.)
            trajectory: Path(s) to trajectory file(s) (DCD, XTC, TRR, etc.)
                       Can be a single path or list of paths for concatenated analysis
            file_times: Optional dictionary mapping trajectory filenames to their
                       simulation durations in nanoseconds. Used for proper time scaling.
                       Example: {"eq1.dcd": 1.0, "eq2.dcd": 2.0, "prod.dcd": 10.0}
            file_strides: Optional per-file frame stride (≥1). Skips frames when iterating.
        """
        try:
            import MDAnalysis as mda
        except ImportError:
            raise ImportError(
                "MDAnalysis is required for trajectory analysis. "
            )

        self.topology = Path(topology)

        # Handle single or multiple trajectories
        if isinstance(trajectory, (str, Path)):
            self.trajectories = [Path(trajectory)]
        else:
            self.trajectories = [Path(t) for t in trajectory]

        # Store file times for proper time scaling
        self.file_times = file_times or {}
        self.file_strides = {
            k: max(1, int(v)) for k, v in (file_strides or {}).items()
        }
        self._file_frame_counts: Dict[str, int] = {}

        # Load trajectories into MDAnalysis
        if len(self.trajectories) == 1:
            self.universe = mda.Universe(str(self.topology), str(self.trajectories[0]))
        else:
            # Concatenate multiple trajectories
            self.universe = mda.Universe(
                str(self.topology), [str(t) for t in self.trajectories]
            )

        logger.info(
            f"Loaded trajectory: {len(self.universe.trajectory)} frames "
            f"from {len(self.trajectories)} file(s)"
        )

    def _ensure_frame_counts(self) -> None:
        if self._file_frame_counts:
            return
        try:
            import MDAnalysis as mda
        except ImportError:
            raise ImportError("MDAnalysis is required")

        for traj_path in self.trajectories:
            temp_universe = mda.Universe(str(self.topology), str(traj_path))
            self._file_frame_counts[traj_path.name] = len(temp_universe.trajectory)

    def _uses_stride(self) -> bool:
        if not self.file_strides:
            return False
        return any(
            max(1, int(_lookup_file_map(self.file_strides, p) or 1)) > 1
            for p in self.trajectories
        )

    def _uniform_stride(self) -> Optional[int]:
        if not self.file_strides:
            return None
        strides = [
            max(1, int(_lookup_file_map(self.file_strides, p) or 1))
            for p in self.trajectories
        ]
        if len(set(strides)) == 1:
            return strides[0]
        return None

    def _analysis_frame_indices(self) -> List[int]:
        """Global frame indices kept when per-file strides are applied."""
        self._ensure_frame_counts()
        indices: List[int] = []
        offset = 0
        for traj_path in self.trajectories:
            n_frames = self._file_frame_counts.get(traj_path.name, 0)
            stride = max(1, int(_lookup_file_map(self.file_strides, traj_path) or 1))
            for local_i in range(0, n_frames, stride):
                indices.append(offset + local_i)
            offset += n_frames
        return indices

    def _analysis_universe_and_ref(
        self, reference_frame: int = 0
    ) -> tuple["Any", int]:
        """
        Return the universe used for analysis and a valid reference frame index.

        When per-file strides are set (>1), builds an in-memory universe containing
        only the kept frames so alignment and statistics never touch skipped frames.
        """
        import MDAnalysis as mda
        import numpy as np

        n_total = len(self.universe.trajectory)
        if n_total == 0:
            return self.universe, 0

        ref_global = max(0, min(int(reference_frame), n_total - 1))

        if not self._uses_stride():
            return self.universe, ref_global

        cache = getattr(self, "_analysis_u_cache", None)
        indices = self._analysis_frame_indices()
        if not indices:
            return self.universe, ref_global

        if cache is None:
            from MDAnalysis.coordinates.memory import MemoryReader

            n_atoms = self.universe.atoms.n_atoms
            coordinates = np.empty((len(indices), n_atoms, 3), dtype=np.float32)
            dimensions = np.zeros((len(indices), 6), dtype=np.float32)
            for i, fi in enumerate(indices):
                self.universe.trajectory[fi]
                coordinates[i] = self.universe.atoms.positions.astype(
                    np.float32, copy=False
                )
                dims = self.universe.trajectory.ts.dimensions
                if dims is not None and len(dims) >= 6:
                    dimensions[i] = dims[:6]
            _fill_missing_box_dimensions(dimensions)
            dt = float(getattr(self.universe.trajectory.ts, "dt", 1.0))
            # Pass the array + format=MemoryReader (order fac = frames, atoms, xyz).
            # Do not pass a pre-built MemoryReader — MDAnalysis wraps it in ChainReader
            # and raises "tuple index out of range".
            self._analysis_u_cache = mda.Universe(
                str(self.topology),
                coordinates,
                format=MemoryReader,
                order="fac",
                dt=dt,
                dimensions=dimensions,
            )
            self._analysis_index_map = indices
            logger.info(
                "Using strided in-memory trajectory: %d / %d frames for analysis",
                len(indices),
                n_total,
            )

        ref_local = (
            self._analysis_index_map.index(ref_global)
            if ref_global in self._analysis_index_map
            else 0
        )
        return self._analysis_u_cache, ref_local

    def clear_analysis_cache(self) -> None:
        """Release in-memory strided trajectory copies to reduce RAM after analysis."""
        self._analysis_u_cache = None
        if hasattr(self, "_analysis_index_map"):
            del self._analysis_index_map

    def time_array_for_analysis(
        self,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
    ) -> "np.ndarray":
        """Time axis aligned to analyzed frames (stride / start / stop / step)."""
        import numpy as np

        full = self._calculate_time_array()
        if start is not None or stop is not None or step is not None:
            s = 0 if start is None else int(start)
            st = 1 if step is None else max(1, int(step))
            stop_idx = len(full) if stop is None else int(stop)
            return np.asarray(full[s:stop_idx:st], dtype=float)
        if self._uses_stride():
            indices = self._analysis_frame_indices()
            return np.asarray([full[i] for i in indices if i < len(full)], dtype=float)
        return np.asarray(full, dtype=float)

    def _calculate_time_array(self) -> "np.ndarray":
        """
        Calculate proper time array based on file_times dict.

        Returns:
            Array of time values in nanoseconds (one per trajectory frame).
        """
        try:
            import numpy as np
            import MDAnalysis as mda
        except ImportError:
            raise ImportError("numpy and MDAnalysis are required")

        if not self.file_times:
            # Fallback: use frame indices with default timestep (2 fs)
            timestep_ps = 0.002  # 2 fs in ps
            n_frames = len(self.universe.trajectory)
            return np.arange(n_frames) * timestep_ps / 1000.0  # Convert to ns

        has_positive_times = any(
            float(v) > 0 for v in self.file_times.values()
        )

        # Calculate time array based on user-specified file durations
        time_array = []
        cumulative_time_ns = 0.0

        for traj_path in self.trajectories:
            temp_universe = mda.Universe(str(self.topology), str(traj_path))
            n_frames = len(temp_universe.trajectory)
            self._file_frame_counts[traj_path.name] = n_frames

            duration_ns = float(_lookup_file_map(self.file_times, traj_path) or 0.0)

            if duration_ns > 0:
                if n_frames == 1:
                    file_times = np.array([cumulative_time_ns], dtype=float)
                else:
                    file_times = np.linspace(
                        cumulative_time_ns,
                        cumulative_time_ns + duration_ns,
                        n_frames,
                    )
                time_array.extend(file_times.tolist())
                cumulative_time_ns += duration_ns
            elif has_positive_times:
                # Other files have assigned times; keep index spacing minimal here
                if n_frames == 1:
                    file_times = np.array([cumulative_time_ns], dtype=float)
                else:
                    file_times = np.linspace(
                        cumulative_time_ns,
                        cumulative_time_ns + max(n_frames - 1, 0) * 0.002,
                        n_frames,
                    )
                time_array.extend(file_times.tolist())
                cumulative_time_ns = float(file_times[-1])
            else:
                file_times = cumulative_time_ns + np.arange(n_frames) * 0.01
                time_array.extend(file_times.tolist())
                cumulative_time_ns += n_frames * 0.01

        return np.array(time_array)

    def calculate_rmsd(
        self,
        selection: str = "protein and backbone",
        reference_frame: int = 0,
        align: bool = True,
        reference_structure: Optional[Union[str, Path]] = None,
    ) -> Dict[str, "np.ndarray"]:
        """
        Calculate RMSD for selected atoms.

        Args:
            selection: MDAnalysis selection string
            reference_frame: Frame to use as reference (0 = first trajectory frame)
            align: If True, perform alignment (rotation + translation) before RMSD.
                  If False, calculate raw coordinate RMSD without alignment.
            reference_structure: Optional PDB/GRO used as the reference instead of
                ``reference_frame``. Topology atom order must match the trajectories.

        Returns:
            Dictionary with 'time' (ns) and 'rmsd' (Angstroms) arrays
        """
        try:
            from MDAnalysis.analysis import rms
            import MDAnalysis as mda
            import numpy as np
        except ImportError:
            raise ImportError("MDAnalysis and numpy are required")

        # Select atoms on the analysis universe (strided in-memory when stride > 1)
        u, ref_local = self._analysis_universe_and_ref(reference_frame)
        ref = u
        if reference_structure:
            ref_path = Path(reference_structure).expanduser()
            if not ref_path.is_file():
                raise ValueError(f"RMSD reference structure not found: {ref_path}")
            ref = mda.Universe(str(self.topology), str(ref_path))
            ref_local = 0
            logger.info("RMSD reference structure: %s (not a trajectory frame)", ref_path.name)

        if align:
            # Single-pass aligned RMSD (QCP superposition per frame). Faster than
            # AlignTraj(in_memory=True) plus a redundant per-frame rms.rmsd loop.
            rmsd_analysis = rms.RMSD(
                u,
                ref,
                select=selection,
                ref_frame=ref_local,
            )
            rmsd_analysis.run()
            rmsd_array = np.asarray(rmsd_analysis.results.rmsd[:, 2], dtype=float)
            logger.info(
                "Aligned RMSD computed for %d analysis frame(s)",
                len(rmsd_array),
            )
        else:
            atoms = u.select_atoms(selection)
            ref_atoms = ref.select_atoms(selection)
            if len(atoms) == 0 or len(ref_atoms) == 0:
                raise ValueError(f"RMSD selection {selection!r} matched 0 atoms.")
            if len(atoms) != len(ref_atoms):
                raise ValueError(
                    f"RMSD selection {selection!r} has {len(atoms)} atoms in the "
                    f"trajectory but {len(ref_atoms)} in the reference structure."
                )
            ref.trajectory[ref_local]
            ref_coords = ref_atoms.positions.astype(np.float64, copy=True)
            n_frames = len(u.trajectory)
            n_atoms = atoms.n_atoms
            coords = np.empty((n_frames, n_atoms, 3), dtype=np.float64)
            for i, _ts in enumerate(u.trajectory):
                coords[i] = atoms.positions
            diff = coords - ref_coords
            rmsd_array = np.sqrt(np.mean(np.sum(diff * diff, axis=2), axis=1))
            logger.info(
                "Unaligned RMSD computed for %d analysis frame(s)",
                len(rmsd_array),
            )

        # Get proper time array (in nanoseconds)
        time_ns = self.time_array_for_analysis()

        return {"time": time_ns, "rmsd": rmsd_array}  # RMSD values in Angstroms

    def calculate_rmsf(
        self, selection: str = "protein and name CA"
    ) -> Dict[str, "np.ndarray"]:
        """
        Calculate RMSF for selected atoms.

        Args:
            selection: MDAnalysis selection string

        Returns:
            Dictionary with 'resids', 'rmsf' (Angstroms), 'resnames', and 'atom_indices' arrays
        """
        try:
            import MDAnalysis as mda
            from MDAnalysis.analysis import rms
            import numpy as np
        except ImportError:
            raise ImportError("MDAnalysis and numpy are required")

        u, _ = self._analysis_universe_and_ref(0)
        atoms = u.select_atoms(selection)

        rmsf_analysis = rms.RMSF(atoms).run()
        rmsf_vals = rmsf_analysis.results.rmsf

        return {
            "resids": atoms.resids,
            "rmsf": rmsf_vals,  # RMSF in Angstroms
            "resnames": atoms.resnames,  # Residue names (e.g., ALA, GLY)
            "atom_indices": atoms.indices,  # Atom indices
        }

    def calculate_distances(
        self, selections: Dict[str, tuple]
    ) -> Dict[str, Dict[str, "np.ndarray"]]:
        """
        Calculate distances between atom selections.

        Args:
            selections: Dictionary of {name: (selection1, selection2)}

        Example:
            selections = {
                "gate": ("resid 50-70 and name CA", "resid 150-170 and name CA"),
                "salt_bridge": ("resid 125 and name NH1", "resid 200 and name OD1 OD2")
            }

        Returns:
            Dictionary with distance data for each named selection.
            Each entry contains 'time' (ns) and 'distance' (Angstroms) arrays.
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("numpy is required")

        u, _ = self._analysis_universe_and_ref(0)
        results = {}

        for name, (sel1, sel2) in selections.items():
            atoms1 = u.select_atoms(sel1)
            atoms2 = u.select_atoms(sel2)

            distances = []

            for fi in range(len(u.trajectory)):
                u.trajectory[fi]
                # Calculate center of mass distance in Angstroms
                dist = np.linalg.norm(atoms1.center_of_mass() - atoms2.center_of_mass())
                distances.append(dist)

            # Get proper time array (in nanoseconds)
            time_ns = self.time_array_for_analysis()

            results[name] = {
                "time": time_ns,
                "distance": np.array(distances),  # Distance in Angstroms
            }

        return results

    def calculate_radius_of_gyration(
        self, selection: str = "protein"
    ) -> Dict[str, "np.ndarray"]:
        """
        Calculate radius of gyration over trajectory.

        Args:
            selection: MDAnalysis selection string

        Returns:
            Dictionary with 'time' (ns) and 'rg' (Angstroms) arrays
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("numpy is required")

        u, _ = self._analysis_universe_and_ref(0)
        atoms = u.select_atoms(selection)

        rg_values = []

        for fi in range(len(u.trajectory)):
            u.trajectory[fi]
            rg_values.append(atoms.radius_of_gyration())  # In Angstroms

        # Get proper time array (in nanoseconds)
        time_ns = self.time_array_for_analysis()

        return {"time": time_ns, "rg": np.array(rg_values)}  # Rg in Angstroms

    def plot_rmsd(
        self,
        selection: str = "protein and backbone",
        reference_frame: int = 0,
        align: bool = True,
        distance_units: str = "Å",
        time_units: str = "ns",
        line_color: str = "blue",
        line_width: float = 1.2,
        line_style: str = "-",
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        show_grid: bool = True,
        xlim: Optional[tuple] = None,
        ylim: Optional[tuple] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        highlight_threshold: Optional[float] = None,
        highlight_color: str = "orange",
        highlight_alpha: float = 0.2,
        show_convergence: bool = True,
        convergence_color: str = "red",
        convergence_style: str = "--",
        convergence_width: float = 1.5,
        hlines: Optional[List[float]] = None,
        hline_colors: Optional[List[str]] = None,
        hline_styles: Optional[List[str]] = None,
        hline_widths: Optional[List[float]] = None,
        vlines: Optional[List[float]] = None,
        vline_colors: Optional[List[str]] = None,
        vline_styles: Optional[List[str]] = None,
        vline_widths: Optional[List[float]] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (10, 6),
        dpi: int = 300,
    ):
        """
        Plot RMSD with full GUI customization options.

        Args:
            selection: MDAnalysis selection string
            reference_frame: Frame to use as reference
            align: Perform alignment before RMSD calculation
            distance_units: 'Å' (angstrom) or 'nm' (nanometer)
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            line_color: Color for the plot line (matplotlib color string or hex)
            line_width: Width of the plot line (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area (hex or 'none' for transparent)
            fig_bg_color: Background color for figure border (hex or 'none')
            text_color: Text/axes color ('Auto', matplotlib color, or hex)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: auto-generated)
            xlabel: X-axis label (default: auto-generated with units)
            ylabel: Y-axis label (default: auto-generated with units)
            highlight_threshold: If set, highlight regions above this RMSD value
            highlight_color: Color for highlight region and line (default: 'orange')
            highlight_alpha: Alpha transparency for highlight fill (default: 0.2)
            show_convergence: Show convergence line (mean of last 20% of trajectory)
            convergence_color: Color for convergence line (default: 'red')
            convergence_style: Line style for convergence line (default: '--')
            convergence_width: Width of convergence line (default: 1.5)
            hlines: List of Y values for horizontal reference lines
            hline_colors: List of colors for horizontal lines (default: cycle through standard colors)
            hline_styles: List of line styles for horizontal lines (default: '--')
            hline_widths: List of line widths for horizontal lines (default: 1.0)
            vlines: List of X values for vertical reference lines
            vline_colors: List of colors for vertical lines (default: cycle through standard colors)
            vline_styles: List of line styles for vertical lines (default: '--')
            vline_widths: List of line widths for vertical lines (default: 1.0)
            save: Filename to save plot
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

        data = self.calculate_rmsd(selection, reference_frame, align)

        # Convert units
        plot_time = data["time"].copy()  # Time is in ns
        plot_rmsd = data["rmsd"].copy()  # RMSD is in Angstroms

        # Convert distance units
        if distance_units == "nm":
            plot_rmsd = plot_rmsd / 10.0  # Convert Å to nm

        # Convert time units
        if time_units == "ps":
            plot_time = plot_time * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = plot_time / 1000.0  # Convert ns to µs

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)

        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)

        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = "black"
            else:
                # Calculate luminance
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

        # Plot data
        ax.plot(
            plot_time,
            plot_rmsd,
            color=line_color,
            linewidth=line_width,
            linestyle=line_style,
            alpha=0.7,
        )

        # Add highlight threshold if specified
        if highlight_threshold is not None:
            threshold_display = (
                highlight_threshold
                if distance_units == "Å"
                else highlight_threshold / 10.0
            )
            ax.axhline(
                y=threshold_display,
                color=highlight_color,
                linestyle=":",
                linewidth=1.5,
                alpha=0.8,
                label=f"Threshold: {highlight_threshold} {distance_units}",
            )
            # Fill region above threshold
            # Convert boolean array to list for type compatibility
            where_condition = (plot_rmsd >= threshold_display).tolist()  # type: ignore
            ax.fill_between(
                plot_time,
                threshold_display,
                plot_rmsd,
                where=where_condition,
                alpha=highlight_alpha,
                color=highlight_color,
                label="Above threshold",
            )

        # Add convergence line (last 20% of trajectory)
        if show_convergence:
            cutoff_idx = int(len(plot_rmsd) * 0.8)
            converged_value = float(np.mean(plot_rmsd[cutoff_idx:]))
            ax.axhline(
                y=converged_value,
                color=convergence_color,
                linestyle=convergence_style,
                linewidth=convergence_width,
                label=f"Converged: {converged_value:.2f} {distance_units}",
            )

        # Add custom horizontal reference lines
        if hlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, yval in enumerate(hlines):
                color = (
                    hline_colors[i]
                    if hline_colors and i < len(hline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    hline_styles[i] if hline_styles and i < len(hline_styles) else "--"
                )
                width = (
                    hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                )
                ax.axhline(
                    y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Add custom vertical reference lines
        if vlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, xval in enumerate(vlines):
                color = (
                    vline_colors[i]
                    if vline_colors and i < len(vline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    vline_styles[i] if vline_styles and i < len(vline_styles) else "--"
                )
                width = (
                    vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                )
                ax.axvline(
                    x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Set labels with appropriate color
        xlabel_text = xlabel or f"Time ({time_units})"
        ylabel_text = ylabel or f"RMSD ({distance_units})"
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)

        # Set title
        title_text = title or f"RMSD - {selection}"
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight="bold")

        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)

        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color)

        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()

    def plot_rmsf(
        self,
        selection: str = "protein and name CA",
        xaxis_type: str = "residue_number",
        show_residue_labels: bool = True,
        residue_name_format: str = "single",
        label_frequency: str = "auto",
        distance_units: str = "Å",
        line_color: str = "blue",
        line_width: float = 1.2,
        line_style: str = "-",
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        show_grid: bool = True,
        xlim: Optional[tuple] = None,
        ylim: Optional[tuple] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        highlight_threshold: Optional[float] = None,
        highlight_color: str = "orange",
        highlight_alpha: float = 0.2,
        hlines: Optional[List[float]] = None,
        hline_colors: Optional[List[str]] = None,
        hline_styles: Optional[List[str]] = None,
        hline_widths: Optional[List[float]] = None,
        vlines: Optional[List[float]] = None,
        vline_colors: Optional[List[str]] = None,
        vline_styles: Optional[List[str]] = None,
        vline_widths: Optional[List[float]] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (12, 6),
        dpi: int = 300,
    ):
        """
        Plot RMSF with full GUI customization including residue labeling.

        Args:
            selection: MDAnalysis selection string
            xaxis_type: X-axis type - 'residue_number', 'residue_type_number', or 'atom_index'
            show_residue_labels: Show residue labels on X-axis
            residue_name_format: 'single' (A, G, V) or 'triple' (ALA, GLY, VAL)
            label_frequency: 'all', 'auto', 'every_2', 'every_5', 'every_10', 'every_20'
            distance_units: 'Å' (angstrom) or 'nm' (nanometer) for Y-axis
            line_color: Color for the plot line
            line_width: Width of the plot line (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area
            fig_bg_color: Background color for figure border
            text_color: Text/axes color ('Auto' or specific color)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: auto-generated)
            xlabel: X-axis label (default: auto-generated)
            ylabel: Y-axis label (default: auto-generated with units)
            highlight_threshold: If set, highlight residues above this RMSF value
            highlight_color: Color for highlight region and line (default: 'orange')
            highlight_alpha: Alpha transparency for highlight fill (default: 0.2)
            hlines: List of Y values for horizontal reference lines
            hline_colors: List of colors for horizontal lines
            hline_styles: List of line styles for horizontal lines (default: '--')
            hline_widths: List of line widths for horizontal lines (default: 1.0)
            vlines: List of X values for vertical reference lines
            vline_colors: List of colors for vertical lines
            vline_styles: List of line styles for vertical lines (default: '--')
            vline_widths: List of line widths for vertical lines (default: 1.0)
            save: Filename to save plot
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

        data = self.calculate_rmsf(selection)

        # Convert distance units
        plot_rmsf = data["rmsf"].copy()  # RMSF is in Angstroms
        if distance_units == "nm":
            plot_rmsf = plot_rmsf / 10.0  # Convert Å to nm

        # Prepare X-axis data and labels
        if xaxis_type == "residue_number":
            x_data = data["resids"]
            xlabel_default = "Residue Number"
        elif xaxis_type == "residue_type_number":
            x_data = data["resids"]
            xlabel_default = "Residue"
        elif xaxis_type == "atom_index":
            x_data = data["atom_indices"]
            xlabel_default = "Atom Index"
        else:
            x_data = data["resids"]
            xlabel_default = "Residue Number"

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)

        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)

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

        # Plot data
        ax.plot(
            x_data,
            plot_rmsf,
            color=line_color,
            linewidth=line_width,
            linestyle=line_style,
        )
        ax.fill_between(x_data, 0, plot_rmsf, alpha=0.3, color=line_color)

        # Highlight flexible regions if threshold specified
        if highlight_threshold is not None:
            threshold_display = (
                highlight_threshold
                if distance_units == "Å"
                else highlight_threshold / 10.0
            )
            ax.axhline(
                y=threshold_display,
                color=highlight_color,
                linestyle=":",
                linewidth=1.5,
                alpha=0.8,
                label=f"Threshold: {highlight_threshold} {distance_units}",
            )
            # Fill region above threshold
            # Convert boolean array to list for type compatibility
            where_condition = (plot_rmsf >= threshold_display).tolist()  # type: ignore
            ax.fill_between(
                x_data,
                threshold_display,
                plot_rmsf,
                where=where_condition,
                alpha=highlight_alpha,
                color=highlight_color,
                label="Above threshold",
            )

        # Add custom horizontal reference lines
        if hlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, yval in enumerate(hlines):
                color = (
                    hline_colors[i]
                    if hline_colors and i < len(hline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    hline_styles[i] if hline_styles and i < len(hline_styles) else "--"
                )
                width = (
                    hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                )
                ax.axhline(
                    y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Add custom vertical reference lines
        if vlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, xval in enumerate(vlines):
                color = (
                    vline_colors[i]
                    if vline_colors and i < len(vline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    vline_styles[i] if vline_styles and i < len(vline_styles) else "--"
                )
                width = (
                    vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                )
                ax.axvline(
                    x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Set labels
        xlabel_text = xlabel or xlabel_default
        ylabel_text = ylabel or f"RMSF ({distance_units})"
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)

        # Set title
        title_text = title or f"RMSF - {selection}"
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight="bold")

        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)

        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color, axis="y")

        # Handle residue labels on X-axis
        if show_residue_labels and xaxis_type == "residue_type_number":
            # Create residue labels with names
            resnames = data["resnames"]
            resids = data["resids"]

            # Convert residue names if needed
            if residue_name_format == "single":
                # Use 1-letter amino acid codes
                aa_codes = {
                    "ALA": "A",
                    "ARG": "R",
                    "ASN": "N",
                    "ASP": "D",
                    "CYS": "C",
                    "GLN": "Q",
                    "GLU": "E",
                    "GLY": "G",
                    "HIS": "H",
                    "ILE": "I",
                    "LEU": "L",
                    "LYS": "K",
                    "MET": "M",
                    "PHE": "F",
                    "PRO": "P",
                    "SER": "S",
                    "THR": "T",
                    "TRP": "W",
                    "TYR": "Y",
                    "VAL": "V",
                }
                labels = [
                    f"{aa_codes.get(name, name)}{resid}"
                    for name, resid in zip(resnames, resids)
                ]
            else:  # triple
                labels = [f"{name}{resid}" for name, resid in zip(resnames, resids)]

            # Determine label frequency
            n_residues = len(resids)
            if label_frequency == "all":
                step = 1
            elif label_frequency == "auto":
                # Auto-determine based on number of residues
                if n_residues < 20:
                    step = 1
                elif n_residues < 50:
                    step = 2
                elif n_residues < 100:
                    step = 5
                elif n_residues < 200:
                    step = 10
                else:
                    step = 20
            elif label_frequency == "every_2":
                step = 2
            elif label_frequency == "every_5":
                step = 5
            elif label_frequency == "every_10":
                step = 10
            elif label_frequency == "every_20":
                step = 20
            else:
                step = 1

            # Set tick positions and labels
            tick_positions = resids[::step]
            tick_labels = [labels[i] for i in range(0, len(labels), step)]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, rotation=45, ha="right")

        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if highlight_threshold is not None:
            ax.legend()

        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()

    def plot_distances(
        self,
        selections: Dict[str, tuple],
        distance_units: str = "Å",
        time_units: str = "ns",
        line_colors: Optional[List[str]] = None,
        line_width: float = 1.2,
        line_style: str = "-",
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        show_grid: bool = True,
        xlim: Optional[tuple] = None,
        ylim: Optional[tuple] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        show_mean_lines: bool = True,
        hlines: Optional[List[float]] = None,
        hline_colors: Optional[List[str]] = None,
        hline_styles: Optional[List[str]] = None,
        hline_widths: Optional[List[float]] = None,
        vlines: Optional[List[float]] = None,
        vline_colors: Optional[List[str]] = None,
        vline_styles: Optional[List[str]] = None,
        vline_widths: Optional[List[float]] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (10, 6),
        dpi: int = 300,
    ):
        """
        Plot distances with full GUI customization.

        Args:
            selections: Dictionary of {name: (selection1, selection2)}
            distance_units: 'Å' (angstrom) or 'nm' (nanometer)
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            line_colors: List of colors for each distance pair (default: auto-cycle)
            line_width: Width of the plot lines (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area
            fig_bg_color: Background color for figure border
            text_color: Text/axes color ('Auto' or specific color)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: "Distance Analysis")
            xlabel: X-axis label (default: auto-generated with units)
            ylabel: Y-axis label (default: auto-generated with units)
            show_mean_lines: Show mean distance lines
            save: Filename to save plot
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

        results = self.calculate_distances(selections)

        # Default colors if not provided
        if line_colors is None:
            line_colors = [
                "#1f77b4",
                "#ff7f0e",
                "#2ca02c",
                "#d62728",
                "#9467bd",
                "#8c564b",
            ]

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)

        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)

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

        # Plot each distance pair
        for i, (name, data) in enumerate(results.items()):
            color = line_colors[i % len(line_colors)]

            # Convert units
            plot_time = data["time"].copy()  # Time is in ns
            plot_distance = data["distance"].copy()  # Distance is in Å

            # Convert distance units
            if distance_units == "nm":
                plot_distance = plot_distance / 10.0  # Convert Å to nm

            # Convert time units
            if time_units == "ps":
                plot_time = plot_time * 1000.0  # Convert ns to ps
            elif time_units == "µs":
                plot_time = plot_time / 1000.0  # Convert ns to µs

            ax.plot(
                plot_time,
                plot_distance,
                color=color,
                linewidth=line_width,
                linestyle=line_style,
                label=name,
                alpha=0.7,
            )

            # Add mean line (calculated from second half of trajectory)
            if show_mean_lines:
                mean_dist = float(
                    np.mean(plot_distance[int(len(plot_distance) * 0.5) :])
                )
                ax.axhline(
                    y=mean_dist, color=color, linestyle="--", linewidth=1.0, alpha=0.5
                )

        # Set labels
        xlabel_text = xlabel or f"Time ({time_units})"
        ylabel_text = ylabel or f"Distance ({distance_units})"
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)

        # Set title
        title_text = title or "Distance Analysis"
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight="bold")

        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)

        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color)

        # Add custom horizontal reference lines
        if hlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, yval in enumerate(hlines):
                color = (
                    hline_colors[i]
                    if hline_colors and i < len(hline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    hline_styles[i] if hline_styles and i < len(hline_styles) else "--"
                )
                width = (
                    hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                )
                ax.axhline(
                    y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Add custom vertical reference lines
        if vlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, xval in enumerate(vlines):
                color = (
                    vline_colors[i]
                    if vline_colors and i < len(vline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    vline_styles[i] if vline_styles and i < len(vline_styles) else "--"
                )
                width = (
                    vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                )
                ax.axvline(
                    x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend()
        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()

    def plot_radius_of_gyration(
        self,
        selection: str = "protein",
        distance_units: str = "Å",
        time_units: str = "ns",
        line_color: str = "purple",
        line_width: float = 1.2,
        line_style: str = "-",
        bg_color: str = "#2b2b2b",
        fig_bg_color: str = "#212121",
        text_color: str = "Auto",
        show_grid: bool = True,
        xlim: Optional[tuple] = None,
        ylim: Optional[tuple] = None,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        show_convergence: bool = True,
        convergence_color: str = "red",
        convergence_style: str = "--",
        convergence_width: float = 1.5,
        hlines: Optional[List[float]] = None,
        hline_colors: Optional[List[str]] = None,
        hline_styles: Optional[List[str]] = None,
        hline_widths: Optional[List[float]] = None,
        vlines: Optional[List[float]] = None,
        vline_colors: Optional[List[str]] = None,
        vline_styles: Optional[List[str]] = None,
        vline_widths: Optional[List[float]] = None,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (10, 6),
        dpi: int = 300,
    ):
        """
        Plot radius of gyration with full GUI customization options.

        Args:
            selection: MDAnalysis selection string (default: "protein")
            distance_units: 'Å' (angstrom) or 'nm' (nanometer)
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            line_color: Color for the plot line (matplotlib color string or hex)
            line_width: Width of the plot line (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area (hex or 'none' for transparent)
            fig_bg_color: Background color for figure border (hex or 'none')
            text_color: Text/axes color ('Auto', matplotlib color, or hex)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: auto-generated)
            xlabel: X-axis label (default: auto-generated with units)
            ylabel: Y-axis label (default: auto-generated with units)
            show_convergence: Show convergence line (mean of last 20% of trajectory)
            convergence_color: Color for convergence line (default: 'red')
            convergence_style: Line style for convergence line (default: '--')
            convergence_width: Width of convergence line (default: 1.5)
            save: Filename to save plot
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

        data = self.calculate_radius_of_gyration(selection)

        # Convert units
        plot_time = data["time"].copy()  # Time is in ns
        plot_rg = data["rg"].copy()  # Rg is in Angstroms

        # Convert distance units
        if distance_units == "nm":
            plot_rg = plot_rg / 10.0  # Convert Å to nm

        # Convert time units
        if time_units == "ps":
            plot_time = plot_time * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = plot_time / 1000.0  # Convert ns to µs

        # Create figure
        fig, ax = plt.subplots(figsize=figsize)

        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)

        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)

        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = "black"
            else:
                # Calculate luminance
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

        # Plot data
        ax.plot(
            plot_time,
            plot_rg,
            color=line_color,
            linewidth=line_width,
            linestyle=line_style,
            alpha=0.7,
        )

        # Add convergence line (last 20% of trajectory)
        if show_convergence:
            cutoff_idx = int(len(plot_rg) * 0.8)
            converged_value = float(np.mean(plot_rg[cutoff_idx:]))
            ax.axhline(
                y=converged_value,
                color=convergence_color,
                linestyle=convergence_style,
                linewidth=convergence_width,
                label=f"Converged: {converged_value:.2f} {distance_units}",
            )

        # Add custom horizontal reference lines
        if hlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, yval in enumerate(hlines):
                color = (
                    hline_colors[i]
                    if hline_colors and i < len(hline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    hline_styles[i] if hline_styles and i < len(hline_styles) else "--"
                )
                width = (
                    hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                )
                ax.axhline(
                    y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Add custom vertical reference lines
        if vlines:
            default_colors = ["gray", "darkgray", "lightgray", "silver"]
            for i, xval in enumerate(vlines):
                color = (
                    vline_colors[i]
                    if vline_colors and i < len(vline_colors)
                    else default_colors[i % len(default_colors)]
                )
                style = (
                    vline_styles[i] if vline_styles and i < len(vline_styles) else "--"
                )
                width = (
                    vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                )
                ax.axvline(
                    x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7
                )

        # Set labels with appropriate color
        xlabel_text = xlabel or f"Time ({time_units})"
        ylabel_text = ylabel or f"Radius of Gyration ({distance_units})"
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)

        # Set title
        title_text = title or f"Radius of Gyration - {selection}"
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight="bold")

        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)

        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color)

        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)

        if show_convergence:
            ax.legend()

        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()

    def plot_summary(
        self,
        save: Optional[str] = None,
        show: bool = False,
        figsize: tuple = (14, 10),
        dpi: int = 300,
    ):
        """
        Create summary analysis plot with RMSD, RMSF, and Rg.

        Args:
            save: Filename to save plot
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

        # Calculate all metrics
        rmsd_data = self.calculate_rmsd("protein and backbone")
        rmsf_data = self.calculate_rmsf("protein and name CA")
        rg_data = self.calculate_radius_of_gyration("protein")

        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle("Trajectory Analysis Summary", fontsize=14, fontweight="bold")

        # RMSD
        axes[0, 0].plot(
            rmsd_data["time"], rmsd_data["rmsd"], "b-", linewidth=1.2, alpha=0.7
        )
        cutoff = int(len(rmsd_data["rmsd"]) * 0.8)
        converged = np.mean(rmsd_data["rmsd"][cutoff:])
        axes[0, 0].axhline(
            y=converged,
            color="r",
            linestyle="--",
            linewidth=1.5,
            label=f"Converged: {converged:.2f} Å",
        )
        axes[0, 0].set_xlabel("Time (ns)")
        axes[0, 0].set_ylabel("RMSD (Å)")
        axes[0, 0].set_title("RMSD - Protein Backbone")
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)

        # RMSF
        axes[0, 1].plot(rmsf_data["resids"], rmsf_data["rmsf"], "g-", linewidth=1.2)
        axes[0, 1].fill_between(
            rmsf_data["resids"], 0, rmsf_data["rmsf"], alpha=0.3, color="g"
        )
        axes[0, 1].set_xlabel("Residue Number")
        axes[0, 1].set_ylabel("RMSF (Å)")
        axes[0, 1].set_title("RMSF - C-alpha Atoms")
        axes[0, 1].grid(True, alpha=0.3, axis="y")

        # Radius of Gyration
        axes[1, 0].plot(
            rg_data["time"], rg_data["rg"], "purple", linewidth=1.2, alpha=0.7
        )
        cutoff = int(len(rg_data["rg"]) * 0.8)
        converged_rg = np.mean(rg_data["rg"][cutoff:])
        axes[1, 0].axhline(
            y=converged_rg,
            color="r",
            linestyle="--",
            linewidth=1.5,
            label=f"Converged: {converged_rg:.2f} Å",
        )
        axes[1, 0].set_xlabel("Time (ns)")
        axes[1, 0].set_ylabel("Rg (Å)")
        axes[1, 0].set_title("Radius of Gyration")
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)

        # Statistics summary
        axes[1, 1].axis("off")
        stats_text = f"""
Statistics Summary:

RMSD:
  Mean: {np.mean(rmsd_data['rmsd']):.2f} ± {np.std(rmsd_data['rmsd']):.2f} Å
  Converged: {converged:.2f} Å

RMSF:
  Mean: {np.mean(rmsf_data['rmsf']):.2f} ± {np.std(rmsf_data['rmsf']):.2f} Å
  Max: {np.max(rmsf_data['rmsf']):.2f} Å (Residue {rmsf_data['resids'][np.argmax(rmsf_data['rmsf'])]})

Radius of Gyration:
  Mean: {np.mean(rg_data['rg']):.2f} ± {np.std(rg_data['rg']):.2f} Å
  Converged: {converged_rg:.2f} Å
  
Trajectory Info:
  Frames: {len(self.universe.trajectory)}
  Time: {rmsd_data['time'][-1]:.2f} ns
"""
        axes[1, 1].text(
            0.1,
            0.5,
            stats_text,
            fontsize=10,
            family="monospace",
            verticalalignment="center",
        )

        plt.tight_layout()

        if save:
            plt.savefig(save, dpi=dpi, bbox_inches="tight")
            logger.info(f"Plot saved: {save}")

        if show:
            plt.show()
        else:
            plt.close()
