"""
Lipid bilayer trajectory analysis.

Area per lipid methods:
- **FATSLiM** (default): external archived ``fatslim apl`` CLI (GPLv3 companion
  env; see ``gatewizard.utils.fatslim_apl`` / ``scripts/install_fatslim_env.sh``).
- **EVAPL** (experimental / not yet validated): freud Voronoi + half-plane clips
  for occupants whose Z lies in the leaflet headgroup Z-range.
- **lipyphilic**: official ``lipyphilic.analysis.area_per_lipid.AreaPerLipid``
  (exclude_sel when the installed lipyphilic supports it — see
  ``requirements-lipyphilic-git.txt`` until a PyPI release includes PR #164).
- **gridmat** (experimental): GateWizard in-process GridMAT-style grid
  assignment, aligned toward GridMAT-MD.pl (not bit-identical).
- **gridmat_md**: external ``GridMAT-MD.pl`` via subprocess (literature-faithful;
  see ``gatewizard.utils.gridmat_md_apl`` / ``scripts/install_gridmat_md.sh``).
- **vtmc**: Voronoi + Monte Carlo reference.

Leaflet assignment and membrane thickness always use lipyphilic.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Union

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

_LIPYPHILIC_INSTALL = "pip install lipyphilic  # or: pip install -e . from the gatewizard repo"
_LIPYPHILIC_GIT_INSTALL = (
    "pip install -r requirements-lipyphilic-git.txt  "
    "# AreaPerLipid exclude_sel (PR #164) until PyPI ships it"
)
_FREUD_INSTALL = "pip install freud-analysis  # or: pip install -e . from the gatewizard repo"


def _require_lipyphilic():
    """Import lipyphilic or raise a clear installation error."""
    try:
        import lipyphilic  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "lipyphilic is required for lipid bilayer analysis. "
            f"Install with: {_LIPYPHILIC_INSTALL}"
        ) from exc


def _require_freud():
    """Import freud or raise a clear installation error."""
    try:
        import freud  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "freud is required for area-per-lipid analysis. "
            f"Install with: {_FREUD_INSTALL}"
        ) from exc


def _lipyphilic_areaperlipid_supports_exclude() -> bool:
    """True when installed lipyphilic AreaPerLipid accepts ``exclude_sel`` (PR #164+)."""
    import inspect

    from lipyphilic.analysis.area_per_lipid import AreaPerLipid

    return "exclude_sel" in inspect.signature(AreaPerLipid.__init__).parameters


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


def _remove_overlapping_xy(positions: "np.ndarray") -> None:
    """Nudge duplicate XY coordinates so freud Voronoi does not fail."""
    import numpy as np

    _, indices, counts = np.unique(
        positions, return_index=True, return_counts=True, axis=0
    )
    if np.max(counts) > 1:
        for duplicate_index in indices[counts > 1]:
            positions[duplicate_index, 0] += 0.001


def _normalize_xy_positions(positions: "np.ndarray") -> "np.ndarray":
    """Return N×3 XY positions with z=0 for freud 2D Voronoi."""
    import numpy as np

    pos = np.asarray(positions, dtype=float).copy()
    if pos.ndim != 2 or pos.shape[1] < 2:
        return np.zeros((0, 3), dtype=float)
    if pos.shape[1] == 2:
        pos3 = np.zeros((pos.shape[0], 3), dtype=float)
        pos3[:, :2] = pos
        pos = pos3
    else:
        pos[:, 2] = 0.0
    _remove_overlapping_xy(pos)
    return pos


def _voronoi_compute(positions: "np.ndarray", lx: float, ly: float):
    """Run periodic 2D freud Voronoi and return the result object."""
    import freud

    pos = _normalize_xy_positions(positions)
    if pos.shape[0] == 0:
        return None
    voro = freud.locality.Voronoi()
    return voro.compute(
        system=({"Lx": float(lx), "Ly": float(ly), "dimensions": 2}, pos)
    )


def _voronoi_atom_areas(positions: "np.ndarray", lx: float, ly: float) -> "np.ndarray":
    """Periodic 2D Voronoi cell areas for XY seed positions (z ignored)."""
    import numpy as np

    result = _voronoi_compute(positions, lx, ly)
    if result is None:
        return np.asarray([], dtype=float)
    return np.asarray(result.volumes, dtype=float)


def _unwrap_xy_to_ref(
    points: "np.ndarray", ref: "np.ndarray", lx: float, ly: float
) -> "np.ndarray":
    """Minimum-image unwrap of XY points around ``ref``."""
    import numpy as np

    pts = np.asarray(points, dtype=float)
    if pts.size == 0:
        return np.empty((0, 2), dtype=float)
    xy = np.atleast_2d(pts)[:, :2].copy()
    ref_xy = np.asarray(ref, dtype=float)[:2]
    d = xy - ref_xy
    if lx > 0:
        d[:, 0] -= float(lx) * np.round(d[:, 0] / float(lx))
    if ly > 0:
        d[:, 1] -= float(ly) * np.round(d[:, 1] / float(ly))
    return d + ref_xy


def _polygon_area_xy(vertices: "np.ndarray") -> float:
    """Shoelace area of a 2D polygon (absolute value)."""
    import numpy as np

    verts = np.asarray(vertices, dtype=float)
    if verts.ndim != 2 or verts.shape[0] < 3:
        return 0.0
    x = verts[:, 0]
    y = verts[:, 1]
    return float(0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))))


def _clip_polygon_halfplane(
    vertices: "np.ndarray",
    ref_xy: "np.ndarray",
    clip_xy: "np.ndarray",
) -> "np.ndarray":
    """Clip a polygon with the perpendicular bisector of ``ref`` and ``clip``.

    Keeps the half-plane containing ``ref``. EVAPL treats each in-cell exclude
    atom as an extra Voronoi site *for this lipid only* (neighbors are not
    re-tessellated).
    """
    import numpy as np

    verts = np.asarray(vertices, dtype=float)
    if verts.ndim != 2 or verts.shape[0] < 3:
        return np.empty((0, 2), dtype=float)
    verts = verts[:, :2]
    ref = np.asarray(ref_xy, dtype=float)[:2]
    clip = np.asarray(clip_xy, dtype=float)[:2]
    mid = 0.5 * (ref + clip)
    normal = ref - clip
    if float(np.dot(normal, normal)) < 1e-16:
        return verts.copy()

    def _inside(point: "np.ndarray") -> bool:
        return float(np.dot(point - mid, normal)) >= -1e-12

    def _intersect(p1: "np.ndarray", p2: "np.ndarray") -> "np.ndarray":
        d1 = float(np.dot(p1 - mid, normal))
        d2 = float(np.dot(p2 - mid, normal))
        denom = d1 - d2
        if abs(denom) < 1e-16:
            return p2.copy()
        t = min(1.0, max(0.0, d1 / denom))
        return p1 + t * (p2 - p1)

    out: list = []
    prev = verts[-1]
    prev_in = _inside(prev)
    for curr in verts:
        curr_in = _inside(curr)
        if curr_in:
            if not prev_in:
                out.append(_intersect(prev, curr))
            out.append(curr.copy())
        elif prev_in:
            out.append(_intersect(prev, curr))
        prev = curr
        prev_in = curr_in
    if len(out) < 3:
        return np.empty((0, 2), dtype=float)
    return np.asarray(out, dtype=float)


def _evapl_clip_areas(
    lipid_positions: "np.ndarray",
    exclude_xyz: "np.ndarray",
    lx: float,
    ly: float,
) -> "np.ndarray":
    """EVAPL multi-atom clip on a global periodic Voronoi.

    Exclude atoms (already gated to the leaflet headgroup Z-range) are assigned
    to the nearest lipid (Voronoi membership, PBC). Each occupied cell is
    successively clipped by the perpendicular bisector of the lipid and each
    in-cell exclude atom (points closer to the lipid than to any such atom).
    Neighbors are not re-tessellated.

    Unclipped cells tile ``Lx × Ly`` exactly; the occupant footprint is
    ``box − sum(clipped areas)`` per leaflet.
    """
    import numpy as np
    from scipy.spatial import cKDTree

    pos = _normalize_xy_positions(lipid_positions)
    n_lipid = int(pos.shape[0])
    if n_lipid == 0:
        return np.asarray([], dtype=float)

    base = _voronoi_compute(pos, lx, ly)
    if base is None:
        return np.full(n_lipid, np.nan, dtype=float)
    areas = np.asarray(base.volumes, dtype=float).copy()

    exclude = np.asarray(exclude_xyz, dtype=float)
    if exclude.ndim != 2 or exclude.shape[0] == 0 or exclude.shape[1] < 2:
        return areas

    lipid_xy = pos[:, :2]
    exclude_xy = exclude[:, :2]
    lipid_pbc = _gridmat_expand_xy_periodic(lipid_xy, lx, ly)
    tree = cKDTree(lipid_pbc)
    _, nn = tree.query(exclude_xy, k=1)
    owners = np.asarray(nn, dtype=int) % n_lipid

    polytopes = base.polytopes
    for lipid_index in np.unique(owners):
        mask = owners == lipid_index
        pts = exclude[mask]
        ref = lipid_xy[int(lipid_index)]
        unwrapped = _unwrap_xy_to_ref(pts[:, :2], ref, lx, ly)

        verts = np.asarray(polytopes[int(lipid_index)], dtype=float)
        if verts.ndim != 2 or verts.shape[0] < 3:
            continue
        verts_xy = _unwrap_xy_to_ref(verts[:, :2], ref, lx, ly)
        old_area = float(areas[lipid_index])

        for atom_xy in unwrapped:
            clipped = _clip_polygon_halfplane(verts_xy, ref, atom_xy)
            if clipped.shape[0] < 3:
                areas[lipid_index] = 0.0
                verts_xy = clipped
                break
            verts_xy = clipped
        else:
            new_area = _polygon_area_xy(verts_xy)
            if np.isfinite(new_area) and 0.0 <= new_area <= old_area + 1e-6:
                areas[lipid_index] = new_area
    return areas


def _clip_lipid_area_evapl(
    lipid_positions: "np.ndarray",
    ref_index: int,
    exclude_xy: "np.ndarray",
    lx: float,
    ly: float,
) -> float:
    """EVAPL area for one lipid (see ``_evapl_clip_areas``)."""
    pos = _normalize_xy_positions(lipid_positions)
    if pos.shape[0] == 0 or ref_index < 0 or ref_index >= pos.shape[0]:
        return float("nan")
    areas = _evapl_clip_areas(pos, exclude_xy, lx, ly)
    if areas.size <= ref_index:
        return float("nan")
    return float(areas[ref_index])


def _leaflet_headgroup_z_range(leaflet_atoms) -> tuple:
    """Return ``(z_min, z_max)`` of leaflet headgroup atoms."""
    import numpy as np

    z = np.asarray(leaflet_atoms.positions[:, 2], dtype=float)
    if z.size == 0:
        return (0.0, 0.0)
    return (float(np.min(z)), float(np.max(z)))


def _exclude_atoms_in_z_range(
    universe,
    exclude_sel: str,
    z_min: float,
    z_max: float,
    keep_z: bool = True,
) -> "np.ndarray":
    """Exclude-atom positions with Z in ``[z_min, z_max]`` (heavy atoms preferred)."""
    import numpy as np

    sel = (exclude_sel or "").strip()
    if not sel:
        return np.empty((0, 3), dtype=float)
    try:
        exclude_ag = universe.select_atoms(sel)
    except Exception as exc:
        logger.warning("Invalid exclude selection %r (%s); skipping exclusion", sel, exc)
        return np.empty((0, 3), dtype=float)
    if len(exclude_ag) == 0:
        return np.empty((0, 3), dtype=float)
    try:
        heavy = exclude_ag.select_atoms("not name H*")
        if len(heavy) > 0:
            exclude_ag = heavy
    except Exception:
        pass
    z = exclude_ag.positions[:, 2]
    mask = (z >= float(z_min)) & (z <= float(z_max))
    if not np.any(mask):
        return np.empty((0, 3), dtype=float)
    exclude_ag = exclude_ag[mask]
    exclude_ag.wrap(inplace=True)
    pos = exclude_ag.positions.copy()
    if not keep_z:
        pos[:, 2] = 0.0
    return pos


def _evapl_leaflet_exclude_xyz(universe, leaflet_atoms, exclude_sel: str) -> "np.ndarray":
    """Protein/exclude XYZ in the leaflet headgroup Z-range for EVAPL."""
    z_min, z_max = _leaflet_headgroup_z_range(leaflet_atoms)
    return _exclude_atoms_in_z_range(universe, exclude_sel, z_min, z_max, keep_z=True)



def _gridmat_round_coord(value: float) -> float:
    """Round Å coords to Perl nm-at-3-decimals precision (0.01 Å)."""
    return round(float(value) / 0.01) * 0.01


def _gridmat_build_xy_grid(
    lx: float,
    ly: float,
    grid_n: int,
    conserve_ratio: bool,
    *,
    x_max: Optional[float] = None,
    y_max: Optional[float] = None,
) -> "np.ndarray":
    """Return N×2 grid points in Å (GridMAT-MD.pl-style frame when maxima given).

    Perl anchors the conserve-ratio grid on lipid-reference XY maxima so the
    span is ``[x_max−Lx, x_max] × [y_max−Ly, y_max]``. When maxima are omitted,
    defaults to ``x_max=Lx``, ``y_max=Ly`` (origin at zero).
    """
    import numpy as np

    n = max(2, int(grid_n))
    xmax = float(lx) if x_max is None else float(x_max)
    ymax = float(ly) if y_max is None else float(y_max)
    xmin = xmax - float(lx)
    ymin = ymax - float(ly)
    if conserve_ratio:
        interval_x = float(lx) / (n - 1)
        grid_y = max(2, int(round(float(ly) / interval_x) + 1))
        xs = np.linspace(xmin, xmax, n)
        ys = np.linspace(ymin, ymax, grid_y)
    else:
        xs = np.linspace(xmin, xmax, n)
        ys = np.linspace(ymin, ymax, n)
    xs = np.asarray([_gridmat_round_coord(v) for v in xs], dtype=float)
    ys = np.asarray([_gridmat_round_coord(v) for v in ys], dtype=float)
    gx, gy = np.meshgrid(xs, ys)
    return np.column_stack([gx.ravel(), gy.ravel()])


def _gridmat_lipid_xy_site(positions: "np.ndarray") -> "np.ndarray":
    """Reference XY for one lipid: single atom or unweighted mean (Perl)."""
    import numpy as np

    pos = np.asarray(positions, dtype=float)
    if pos.ndim != 2 or pos.shape[0] == 0:
        return np.array([np.nan, np.nan], dtype=float)
    if pos.shape[0] == 1:
        return pos[0, :2].copy()
    return np.mean(pos[:, :2], axis=0)


def _gridmat_midplane_leaflets(z_values: "np.ndarray") -> "np.ndarray":
    """Leaflet signs from global Z midplane: ``z > mid → +1`` else ``-1``."""
    import numpy as np

    z = np.asarray(z_values, dtype=float)
    if z.size == 0:
        return np.empty(0, dtype=int)
    mid = float(np.mean(z))
    return np.where(z > mid, 1, -1).astype(int)


def _gridmat_protein_sites(
    universe,
    leaflet_atoms,
    exclude_sel: str,
    precision: float,
    z_min: float,
    z_max: float,
    lx: float,
    ly: float,
) -> "np.ndarray":
    """Protein XY sites embedded in a leaflet (GridMAT-MD 'offending atoms' logic)."""
    import numpy as np

    try:
        protein = universe.select_atoms(exclude_sel)
    except Exception:
        return np.empty((0, 3), dtype=float)
    if len(protein) == 0 or len(leaflet_atoms) == 0:
        return np.empty((0, 3), dtype=float)

    prec = float(precision)
    hg = leaflet_atoms.positions
    hg_xy = _gridmat_expand_xy_periodic(hg[:, :2], lx, ly)
    hg_z = np.tile(hg[:, 2], 9)
    sites = []
    for pos in protein.positions:
        px, py, pz = pos
        if pz < z_min or pz > z_max:
            continue
        dists = np.hypot(hg_xy[:, 0] - px, hg_xy[:, 1] - py)
        near = dists <= prec
        if not np.any(near):
            continue
        z_lo = float(np.min(hg_z[near]))
        z_hi = float(np.max(hg_z[near]))
        if z_lo <= pz <= z_hi:
            sites.append([px, py, pz])
    if not sites:
        return np.empty((0, 3), dtype=float)
    return np.asarray(sites, dtype=float)


def _gridmat_expand_xy_periodic(xy: "np.ndarray", lx: float, ly: float) -> "np.ndarray":
    """Replicate XY sites across 3×3 periodic images (GridMAT-MD style)."""
    import numpy as np

    pts = np.asarray(xy, dtype=float)
    if pts.size == 0:
        return pts.reshape(0, 2)
    shifts = [(0.0, 0.0), (-lx, ly), (0.0, ly), (lx, ly), (-lx, 0.0), (lx, 0.0), (-lx, -ly), (0.0, -ly), (lx, -ly)]
    out = []
    for dx, dy in shifts:
        dup = pts.copy()
        dup[:, 0] += dx
        dup[:, 1] += dy
        out.append(dup)
    return np.vstack(out)


def _gridmat_assign_areas(
    grid_xy: "np.ndarray",
    lipid_xy: "np.ndarray",
    lipid_resindices: "np.ndarray",
    protein_xy: "np.ndarray",
    lx: float,
    ly: float,
) -> Dict[int, float]:
    """Assign grid cells to nearest site; return area (Å²) per residue index."""
    import numpy as np
    from collections import Counter

    n_grid = grid_xy.shape[0]
    if n_grid == 0:
        return {}
    box_area = float(lx) * float(ly)
    cell_frac = 1.0 / n_grid

    lipid_xy_pbc = _gridmat_expand_xy_periodic(lipid_xy, lx, ly)
    lipid_res_pbc = np.tile(lipid_resindices, 9) if lipid_xy.shape[0] else lipid_resindices
    protein_xy_pbc = _gridmat_expand_xy_periodic(protein_xy, lx, ly)

    owners: list = []
    for gxy in grid_xy:
        best_dist = np.inf
        owner = None
        if lipid_xy_pbc.shape[0] > 0:
            dists = np.hypot(lipid_xy_pbc[:, 0] - gxy[0], lipid_xy_pbc[:, 1] - gxy[1])
            idx = int(np.argmin(dists))
            best_dist = float(dists[idx])
            owner = ("lipid", int(lipid_res_pbc[idx]))
        if protein_xy_pbc.shape[0] > 0:
            pdists = np.hypot(protein_xy_pbc[:, 0] - gxy[0], protein_xy_pbc[:, 1] - gxy[1])
            pidx = int(np.argmin(pdists))
            if float(pdists[pidx]) < best_dist:
                owner = ("protein", -1)
        if owner is not None:
            owners.append(owner)

    counts = Counter(owners)
    out: Dict[int, float] = {}
    for resindex in np.unique(lipid_resindices):
        n_owned = counts.get(("lipid", int(resindex)), 0)
        out[int(resindex)] = n_owned * cell_frac * box_area
    return out


def _area_per_lipid_frame_gridmat(
    membrane_atoms,
    frame_leaflets: "np.ndarray",
    box_xy: tuple,
    exclude_sel: Optional[str],
    gridmat_n: int,
    gridmat_precision: float,
    out_areas: "np.ndarray",
    frame_index: int,
) -> None:
    """Experimental in-process GridMAT (aligned toward GridMAT-MD.pl).

    Uses global lipid-reference Z midplane for leaflet ownership (not
    ``AssignLeaflets``), Perl-style grid anchoring on XY maxima, and
    unweighted multi-atom XY sites. Not bit-identical to ``GridMAT-MD.pl``;
    use ``apl_method='gridmat_md'`` for the original tool.
    """
    import numpy as np

    del frame_leaflets  # GridMAT cell ownership uses midplane, not AssignLeaflets
    lx, ly = box_xy
    universe = membrane_atoms.universe
    exclude_sel = (exclude_sel or "").strip()

    # Optional centering (matches analysis_03 / GridMAT-MD CLI prep).
    dims = universe.dimensions
    if len(membrane_atoms) > 0 and dims is not None and len(dims) >= 3:
        com = membrane_atoms.center_of_mass()
        box_center = np.asarray(dims[:3], dtype=float) * 0.5
        universe.atoms.translate(box_center - com)
        universe.atoms.wrap(inplace=True)

    membrane_atoms.wrap(inplace=True)
    resindices_all = membrane_atoms.residues.resindices

    # One reference site per residue (unweighted XY / Z for midplane).
    lipid_xy_all = []
    lipid_z_all = []
    lipid_resindices_all = []
    for res in membrane_atoms.residues:
        ra = res.atoms.intersection(membrane_atoms)
        if len(ra) == 0:
            continue
        pos = ra.positions
        lipid_xy_all.append(_gridmat_lipid_xy_site(pos))
        lipid_z_all.append(float(np.mean(pos[:, 2])))
        lipid_resindices_all.append(int(res.resindex))

    if not lipid_resindices_all:
        return

    lipid_xy_all = np.asarray(lipid_xy_all, dtype=float)
    lipid_z_all = np.asarray(lipid_z_all, dtype=float)
    lipid_resindices_all = np.asarray(lipid_resindices_all, dtype=int)
    midplane_signs = _gridmat_midplane_leaflets(lipid_z_all)

    x_max = float(np.max(lipid_xy_all[:, 0]))
    y_max = float(np.max(lipid_xy_all[:, 1]))
    grid_xy = _gridmat_build_xy_grid(
        lx, ly, gridmat_n, conserve_ratio=True, x_max=x_max, y_max=y_max
    )

    for leaflet_sign in (-1, 1):
        mask = midplane_signs == leaflet_sign
        if not np.any(mask):
            continue
        leaflet_xy = lipid_xy_all[mask]
        leaflet_res = lipid_resindices_all[mask]
        leaflet_z = lipid_z_all[mask]
        z_min = float(np.min(leaflet_z))
        z_max = float(np.max(leaflet_z))

        leaf_atom_mask = np.isin(membrane_atoms.resindices, leaflet_res)
        leaflet_atoms = membrane_atoms[leaf_atom_mask]

        protein_xy = np.empty((0, 2), dtype=float)
        if exclude_sel and len(leaflet_atoms) > 0:
            prot = _gridmat_protein_sites(
                universe,
                leaflet_atoms,
                exclude_sel,
                gridmat_precision,
                z_min,
                z_max,
                lx,
                ly,
            )
            if prot.shape[0] > 0:
                protein_xy = prot[:, :2]

        areas = _gridmat_assign_areas(
            grid_xy, leaflet_xy, leaflet_res, protein_xy, lx, ly
        )
        for resindex, area in areas.items():
            m = resindices_all == resindex
            if np.any(m):
                out_areas[m, frame_index] = area


def _vtmc_leaflet_protein_xy(
    universe,
    leaflet_atoms,
    exclude_sel: str,
    z_min: float,
    z_max: float,
) -> "np.ndarray":
    """Protein XY sites in the leaflet Z slab for VTMC disk sampling."""
    import numpy as np

    del leaflet_atoms  # Z range is passed explicitly (same as EVAPL gate)
    pos = _exclude_atoms_in_z_range(universe, exclude_sel, z_min, z_max, keep_z=False)
    if pos.shape[0] == 0:
        return np.empty((0, 2), dtype=float)
    return np.asarray(pos[:, :2], dtype=float)


def _vtmc_assign_areas(
    lipid_xy: "np.ndarray",
    lipid_resindices: "np.ndarray",
    protein_xy: "np.ndarray",
    lx: float,
    ly: float,
    n_samples: int,
    protein_radius: float,
    rng: "np.random.Generator",
) -> Dict[int, float]:
    """Monte Carlo Voronoi APL (Mori, Ogushi, Sugita 2012 VTMC-style).

    Random XY points are assigned to the nearest lipid (periodic Voronoi).
    Points that fall inside a protein atom disk of radius ``protein_radius``
    are excluded from lipid area (Monte Carlo protein footprint).
    """
    import numpy as np
    from scipy.spatial import cKDTree

    n_lip = int(lipid_xy.shape[0])
    if n_lip == 0:
        return {}
    n_mc = max(1000, int(n_samples))
    box_area = float(lx) * float(ly)
    samples = np.column_stack(
        [
            rng.uniform(0.0, float(lx), size=n_mc),
            rng.uniform(0.0, float(ly), size=n_mc),
        ]
    )

    lipid_pbc = _gridmat_expand_xy_periodic(lipid_xy, lx, ly)
    lipid_res_pbc = np.tile(lipid_resindices, 9)
    lipid_tree = cKDTree(lipid_pbc)
    _, lipid_nn = lipid_tree.query(samples, k=1)
    owners = lipid_res_pbc[np.asarray(lipid_nn, dtype=int)]

    is_protein = np.zeros(n_mc, dtype=bool)
    if protein_xy.shape[0] > 0 and float(protein_radius) > 0:
        prot_pbc = _gridmat_expand_xy_periodic(protein_xy, lx, ly)
        prot_tree = cKDTree(prot_pbc)
        prot_dist, _ = prot_tree.query(samples, k=1)
        is_protein = np.asarray(prot_dist, dtype=float) <= float(protein_radius)

    lipid_mask = ~is_protein
    out: Dict[int, float] = {}
    for resindex in np.unique(lipid_resindices):
        n_owned = int(np.sum((owners == int(resindex)) & lipid_mask))
        out[int(resindex)] = n_owned / n_mc * box_area
    return out


def _area_per_lipid_frame_vtmc(
    membrane_atoms,
    frame_leaflets: "np.ndarray",
    box_xy: tuple,
    exclude_sel: Optional[str],
    exclude_cutoff: float,
    exclude_dim: int,
    vtmc_n_samples: int,
    vtmc_protein_radius: float,
    out_areas: "np.ndarray",
    frame_index: int,
) -> None:
    """VTMC-style APL: Voronoi ownership + Monte Carlo protein disks.

    Reference: Mori, Ogushi, Sugita, J. Comput. Chem. (2012) doi:10.1002/jcc.21973.

    ``exclude_cutoff`` / ``exclude_dim`` are accepted for API symmetry with other
    APL methods; VTMC selects protein heavy atoms by leaflet headgroup Z-range.
    """
    import numpy as np

    del exclude_cutoff, exclude_dim  # unused; Z-slab selection instead
    lx, ly = box_xy
    universe = membrane_atoms.universe
    exclude_sel = (exclude_sel or "").strip()
    # Deterministic per-frame RNG so repeated runs match.
    rng = np.random.default_rng(10_000 + int(frame_index))
    resindices_all = membrane_atoms.residues.resindices

    for leaflet_sign in (-1, 1):
        leaflet_res = membrane_atoms.residues[frame_leaflets == leaflet_sign]
        if len(leaflet_res) == 0:
            continue
        leaflet_atoms = leaflet_res.atoms.intersection(membrane_atoms)
        if len(leaflet_atoms) == 0:
            continue
        leaflet_atoms.wrap(inplace=True)
        pos = leaflet_atoms.positions
        z_min = float(np.min(pos[:, 2]))
        z_max = float(np.max(pos[:, 2]))

        lipid_xy = []
        lipid_resindices = []
        for res in leaflet_res:
            ra = res.atoms.intersection(membrane_atoms)
            pos_hg = ra.positions
            if pos_hg.shape[0] == 1:
                xy = pos_hg[0, :2]
            else:
                xy = ra.center_of_mass()[:2]
            lipid_xy.append(xy)
            lipid_resindices.append(int(res.resindex))
        lipid_xy = np.asarray(lipid_xy, dtype=float)
        lipid_resindices = np.asarray(lipid_resindices, dtype=int)

        protein_xy = np.empty((0, 2), dtype=float)
        if exclude_sel:
            protein_xy = _vtmc_leaflet_protein_xy(
                universe, leaflet_atoms, exclude_sel, z_min, z_max
            )

        areas = _vtmc_assign_areas(
            lipid_xy,
            lipid_resindices,
            protein_xy,
            lx,
            ly,
            vtmc_n_samples,
            vtmc_protein_radius,
            rng,
        )
        for resindex, area in areas.items():
            mask = resindices_all == resindex
            if np.any(mask):
                out_areas[mask, frame_index] = area


def _resolve_apl_method(apl_method: Optional[str], exclude_sel: Optional[str] = None) -> str:
    """Normalize APL backend. Default (including ``auto``) is ``fatslim``.

    ``exclude_sel`` is accepted for call-site compatibility; it does not change
    resolution (unlike the former EVAPL-on-exclude auto rule).
    """
    _ = exclude_sel
    method = (apl_method or "auto").strip().lower().replace("-", "_")
    if method in {"", "auto"}:
        return "fatslim"
    if method in {"fatslim", "fatslim_cli"}:
        return "fatslim"
    if method in {"lipyphilic", "voronoi", "standard"}:
        return "lipyphilic"
    if method == "evapl":
        return "evapl"
    if method in {"gridmat", "grid"}:
        return "gridmat"
    if method in {"gridmat_md", "gridmat_pl", "gridmat_cli"}:
        return "gridmat_md"
    if method in {"vtmc", "voronoi_mc", "mori"}:
        return "vtmc"
    raise ValueError(
        f"Unsupported apl_method {apl_method!r}. "
        "Supported: auto, fatslim, evapl, lipyphilic, gridmat, gridmat_md, vtmc"
    )


def _exclude_atom_positions(
    universe,
    leaflet_atoms,
    exclude_sel: Optional[str],
    exclude_cutoff: float,
    exclude_dim: int,
    keep_z: bool = False,
) -> "np.ndarray":
    """Positions of exclude atoms near a leaflet (empty if none).

    ``keep_z=False`` (default) zeros Z for 2D Voronoi sites. EVAPL does not
    use this helper for clipping; it selects exclude atoms by leaflet
    headgroup Z-range via ``_evapl_leaflet_exclude_xyz``.
    """
    import numpy as np
    from MDAnalysis.lib.distances import distance_array

    sel = (exclude_sel or "").strip()
    if not sel:
        return np.empty((0, 3), dtype=float)
    try:
        exclude_ag = universe.select_atoms(sel)
    except Exception as exc:
        logger.warning("Invalid exclude selection %r (%s); skipping exclusion", sel, exc)
        return np.empty((0, 3), dtype=float)
    if len(exclude_ag) == 0 or len(leaflet_atoms) == 0:
        return np.empty((0, 3), dtype=float)

    cutoff = float(exclude_cutoff or 0.0)
    if cutoff > 0:
        if int(exclude_dim) == 1:
            mid_z = float(np.mean(leaflet_atoms.positions[:, 2]))
            mask = np.abs(exclude_ag.positions[:, 2] - mid_z) <= cutoff
        else:
            dists = distance_array(
                exclude_ag.positions,
                leaflet_atoms.positions,
                box=universe.dimensions,
            )
            mask = np.min(dists, axis=1) <= cutoff
        if not np.any(mask):
            return np.empty((0, 3), dtype=float)
        exclude_ag = exclude_ag[mask]

    exclude_ag.wrap(inplace=True)
    pos = exclude_ag.positions.copy()
    if not keep_z:
        pos[:, 2] = 0.0
    return pos


def _area_per_lipid_frame(
    membrane_atoms,
    frame_leaflets: "np.ndarray",
    box_xy: tuple,
    exclude_sel: Optional[str],
    exclude_cutoff: float,
    exclude_dim: int,
    lipid_species: "np.ndarray",
    num_seeds: Dict[str, int],
    out_areas: "np.ndarray",
    frame_index: int,
    apl_method: str = "evapl",
    gridmat_n: int = 20,
    gridmat_precision: float = 13.0,
    vtmc_n_samples: int = 50_000,
    vtmc_protein_radius: float = 1.7,
) -> None:
    """Fill one frame of per-residue APL (upper and lower leaflets)."""
    if apl_method == "gridmat":
        _area_per_lipid_frame_gridmat(
            membrane_atoms,
            frame_leaflets,
            box_xy,
            exclude_sel,
            gridmat_n,
            gridmat_precision,
            out_areas,
            frame_index,
        )
        return
    if apl_method == "vtmc":
        _area_per_lipid_frame_vtmc(
            membrane_atoms,
            frame_leaflets,
            box_xy,
            exclude_sel,
            exclude_cutoff,
            exclude_dim,
            vtmc_n_samples,
            vtmc_protein_radius,
            out_areas,
            frame_index,
        )
        return
    if apl_method == "evapl":
        _area_per_lipid_frame_evapl(
            membrane_atoms,
            frame_leaflets,
            box_xy,
            exclude_sel,
            exclude_cutoff,
            exclude_dim,
            lipid_species,
            num_seeds,
            out_areas,
            frame_index,
        )
        return

    raise ValueError(
        f"Unknown apl_method for freud frame loop: {apl_method!r}. "
        "Use evapl, gridmat, or vtmc here; lipyphilic uses AreaPerLipid."
    )


def _area_per_lipid_frame_evapl(
    membrane_atoms,
    frame_leaflets: "np.ndarray",
    box_xy: tuple,
    exclude_sel: Optional[str],
    exclude_cutoff: float,
    exclude_dim: int,
    lipid_species: "np.ndarray",
    num_seeds: Dict[str, int],
    out_areas: "np.ndarray",
    frame_index: int,
) -> None:
    """EVAPL per-lipid APL with multi-atom headgroup-Z clips.

    ``exclude_cutoff`` / ``exclude_dim`` are accepted for API symmetry with other
    APL methods; EVAPL selects exclude atoms by leaflet headgroup Z-range only.
    """
    import numpy as np

    del exclude_cutoff, exclude_dim  # unused; headgroup Z-range instead
    lx, ly = box_xy
    universe = membrane_atoms.universe
    exclude_sel = (exclude_sel or "").strip()

    for leaflet_sign in (-1, 1):
        leaflet_res = membrane_atoms.residues[frame_leaflets == leaflet_sign]
        if len(leaflet_res) == 0:
            continue
        leaflet_atoms = leaflet_res.atoms.intersection(membrane_atoms)
        if len(leaflet_atoms) == 0:
            continue
        leaflet_atoms.wrap(inplace=True)
        raw_pos = leaflet_atoms.positions
        lipid_pos = _normalize_xy_positions(raw_pos)

        if not exclude_sel:
            atom_areas = _voronoi_atom_areas(lipid_pos, lx, ly)
        else:
            leaflet_exclude = _evapl_leaflet_exclude_xyz(
                universe, leaflet_atoms, exclude_sel
            )
            atom_areas = _evapl_clip_areas(lipid_pos, leaflet_exclude, lx, ly)

        for species in lipid_species:
            species_indices = leaflet_atoms.resnames == species
            if not np.any(species_indices):
                continue
            species_apl = atom_areas[species_indices]
            species_atoms = leaflet_atoms[species_indices]
            seeds = int(num_seeds[species])
            species_apl = np.sum(
                species_apl.reshape(species_atoms.n_residues, seeds),
                axis=1,
            )
            species_resindices = np.isin(
                membrane_atoms.residues.resindices,
                species_atoms.residues.resindices,
                assume_unique=True,
            )
            out_areas[species_resindices, frame_index] = species_apl


class BilayerTrajectoryAnalyzer:
    """
    Lipid bilayer analysis wrapper built on MDAnalysis, freud, and lipyphilic.

    Supports:
    - Area per lipid (periodic Voronoi; optional protein footprint exclusion)
    - Membrane thickness (lipyphilic interleaflet headgroup distance)

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
        file_strides: Optional[Dict[str, int]] = None,
    ):
        from gatewizard.utils.trajectory_analysis import TrajectoryAnalyzer

        self._trajectory = TrajectoryAnalyzer(
            topology,
            trajectory,
            file_times=file_times,
            file_strides=file_strides,
        )
        self._z_centered_for: Optional[str] = None

    @property
    def universe(self):
        u, _ = self._trajectory._analysis_universe_and_ref(0)
        return u

    def _calculate_time_array(self):
        return self._trajectory._calculate_time_array()

    def _align_time_ns(
        self,
        n_frames: int,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
    ):
        if self._trajectory._uses_stride():
            return self._trajectory.time_array_for_analysis()
        from gatewizard.utils.trajectory_analysis import _align_time_to_frame_count

        full = self._calculate_time_array()
        return _align_time_to_frame_count(full, n_frames, start=start, stop=stop, step=step)

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
        exclude_sel: Optional[str] = "protein",
        exclude_cutoff: float = 30.0,
        exclude_dim: int = 3,
        apl_method: Optional[str] = "fatslim",
        gridmat_n: int = 20,
        gridmat_precision: float = 13.0,
        vtmc_n_samples: int = 50_000,
        vtmc_protein_radius: float = 1.7,
        fatslim_nthreads: int = 1,
        fatslim_jobs: int = 1,
        gridmat_md_jobs: Optional[int] = None,
        start: Optional[int] = None,
        stop: Optional[int] = None,
        step: Optional[int] = None,
        verbose: bool = False,
    ) -> Dict[str, Any]:
        """
        Calculate area per lipid (default: external FATSLiM CLI).

        Args:
            lipid_sel: Atom selection for Voronoi tessellation (e.g. MARTINI
                ``name GL1 GL2 ROH`` or all-atom ``name PO4``).
            leaflet_lipid_sel: Selection for leaflet assignment. Defaults to
                ``lipid_sel``.
            exclude_sel: Atoms treated as membrane occupants / intercalators
                (e.g. ``protein``, a peptide, DNA, ligands). Empty/None
                disables exclusion. With ``apl_method='fatslim'`` (default),
                mapped to FATSLiM ``--interacting-group``. With
                ``apl_method='evapl'`` (experimental / not yet validated),
                exclude atoms whose Z lies in the leaflet headgroup Z-range
                shrink each owning lipid's Voronoi cell via successive
                half-plane clips. With ``apl_method='lipyphilic'``, passed to
                upstream ``AreaPerLipid`` when the installed lipyphilic supports
                ``exclude_sel`` (PR #164+); otherwise raises.
            exclude_cutoff: Å cutoff for exclude atoms (lipyphilic AreaPerLipid,
                VTMC). **Ignored by FATSLiM, EVAPL, and both GridMAT paths**.
                Default 30 Å.
            exclude_dim: ``3`` = full 3D distance; ``1`` = |z − leaflet
                midplane|. Used by lipyphilic AreaPerLipid and VTMC.
                **Ignored by FATSLiM, EVAPL, and GridMAT.**
            apl_method: ``fatslim`` / ``auto`` (default; external FATSLiM CLI),
                ``evapl`` (experimental Exclusion-aware Voronoi),
                ``lipyphilic`` (official ``lipyphilic.AreaPerLipid``),
                ``gridmat`` (experimental GW GridMAT; ``gridmat_n`` /
                ``gridmat_precision``), ``gridmat_md`` (external GridMAT-MD.pl),
                or ``vtmc`` (experimental GW VTMC; ``vtmc_n_samples`` /
                ``vtmc_protein_radius``).
            fatslim_nthreads: FATSLiM ``--nthreads`` per process (``-1`` = all
                CPUs). Prefer 1–4 when using ``fatslim_jobs`` > 1.
            fatslim_jobs: Parallel FATSLiM frame chunks (``--begin-frame`` /
                ``--end-frame``). Prefer several jobs with small nthreads over
                one job with all CPUs.
            gridmat_md_jobs: Parallel Perl workers for ``gridmat_md`` (default
                ``min(8, cpu_count)``).
            start, stop, step: Trajectory frame range.
            verbose: Show a progress bar while iterating frames.

        Returns:
            Dict with time (ns), per-lipid areas, leaflet means, and statistics.

        Note:
            Pure bilayers (no matching exclude atoms) behave like classic box
            Voronoi for in-process methods. FATSLiM requires the companion
            ``fatslim`` binary (see ``scripts/install_fatslim_env.sh``).
            ``gridmat_md`` requires ``GridMAT-MD.pl`` (see
            ``scripts/install_gridmat_md.sh``).
        """
        import numpy as np
        from gatewizard.utils.fatslim_apl import run_fatslim_apl

        resolved_method = _resolve_apl_method(apl_method, exclude_sel)
        _require_lipyphilic()
        leaflet_sel = leaflet_lipid_sel or lipid_sel
        self._ensure_membrane_centered_in_z(leaflet_sel)
        leaflets = self._assign_leaflets(
            leaflet_sel, start=start, stop=stop, step=step, verbose=verbose
        )
        leaflet_data = np.asarray(_analysis_result(leaflets, "leaflets"))

        u = self.universe
        dims = u.dimensions
        if dims is not None and not np.allclose(dims[3:], 90.0):
            raise ValueError(
                "Area per lipid requires an orthorhombic box — triclinic systems "
                "are not supported."
            )

        membrane = u.select_atoms(lipid_sel, updating=False)
        if membrane.n_residues == 0:
            raise ValueError(f"No residues match lipid selection {lipid_sel!r}")
        if leaflet_data.shape[0] != membrane.n_residues:
            raise ValueError(
                f"'leaflets' has {leaflet_data.shape[0]} residues but lipid_sel "
                f"matches {membrane.n_residues}"
            )

        traj = u.trajectory
        n_traj = len(traj)
        frame_indices = list(range(*slice(start, stop, step).indices(n_traj)))
        n_frames = len(frame_indices)
        if leaflet_data.ndim == 2 and leaflet_data.shape[1] != n_frames:
            raise ValueError(
                "The frames to analyse must be identical to those used in assigning "
                "lipids to leaflets."
            )

        if resolved_method == "fatslim":
            import tempfile

            work_dir = Path(tempfile.mkdtemp(prefix="gatewizard_fatslim_"))
            logger.info("FATSLiM export directory: %s", work_dir)
            fatslim_out = run_fatslim_apl(
                u,
                lipid_sel=lipid_sel,
                exclude_sel=exclude_sel,
                frame_indices=frame_indices,
                work_dir=work_dir,
                fatslim_nthreads=fatslim_nthreads,
                fatslim_jobs=fatslim_jobs,
            )
            area_array = np.asarray(fatslim_out["areas"], dtype=float)
            if area_array.shape[1] != n_frames:
                n_frames = int(area_array.shape[1])
                if leaflet_data.ndim == 2:
                    leaflet_data = leaflet_data[:, :n_frames]
        elif resolved_method == "gridmat_md":
            import tempfile

            from gatewizard.utils.gridmat_md_apl import run_gridmat_md_apl

            work_dir = Path(tempfile.mkdtemp(prefix="gatewizard_gridmat_md_"))
            logger.info("GridMAT-MD.pl export directory: %s", work_dir)
            gridmat_out = run_gridmat_md_apl(
                u,
                lipid_sel=lipid_sel,
                exclude_sel=exclude_sel,
                frame_indices=frame_indices,
                work_dir=work_dir,
                gridmat_n=gridmat_n,
                gridmat_precision=gridmat_precision,
                gridmat_md_jobs=gridmat_md_jobs,
                leaflet_data=leaflet_data,
            )
            area_array = np.asarray(gridmat_out["areas"], dtype=float)
            if area_array.shape[1] != n_frames:
                n_frames = int(area_array.shape[1])
                if leaflet_data.ndim == 2:
                    leaflet_data = leaflet_data[:, :n_frames]
        elif resolved_method == "lipyphilic":
            area_array = self._area_per_lipid_via_areaperlipid(
                lipid_sel=lipid_sel,
                leaflet_data=leaflet_data,
                exclude_sel=exclude_sel,
                exclude_cutoff=exclude_cutoff,
                exclude_dim=exclude_dim,
                start=start,
                stop=stop,
                step=step,
                verbose=verbose,
            )
        else:
            if resolved_method == "evapl":
                _require_freud()
            lipid_species = np.unique(membrane.resnames)
            num_lipids = {
                lipid: int(np.sum(membrane.residues.resnames == lipid))
                for lipid in lipid_species
            }
            num_seeds = {
                lipid: int(np.sum(membrane.resnames == lipid) // num_lipids[lipid])
                for lipid in lipid_species
            }
            area_array = np.full((membrane.n_residues, n_frames), np.nan, dtype=float)
            from MDAnalysis.lib.log import ProgressBar

            frame_slice = traj[start:stop:step]
            iterator = ProgressBar(frame_slice) if verbose else frame_slice
            for frame_index, ts in enumerate(iterator):
                frame_leaflets = (
                    leaflet_data[:, frame_index] if leaflet_data.ndim == 2 else leaflet_data
                )
                lx = float(ts.dimensions[0])
                ly = float(ts.dimensions[1])
                _area_per_lipid_frame(
                    membrane,
                    frame_leaflets,
                    (lx, ly),
                    exclude_sel,
                    exclude_cutoff,
                    exclude_dim,
                    lipid_species,
                    num_seeds,
                    area_array,
                    frame_index,
                    apl_method=resolved_method,
                    gridmat_n=gridmat_n,
                    gridmat_precision=gridmat_precision,
                    vtmc_n_samples=vtmc_n_samples,
                    vtmc_protein_radius=vtmc_protein_radius,
                )

        return self._finalize_area_per_lipid_result(
            area_array=area_array,
            leaflet_data=leaflet_data,
            lipid_sel=lipid_sel,
            exclude_sel=exclude_sel,
            exclude_cutoff=exclude_cutoff,
            resolved_method=resolved_method,
            n_frames=n_frames,
            start=start,
            stop=stop,
            step=step,
        )

    def _area_per_lipid_via_areaperlipid(
        self,
        *,
        lipid_sel: str,
        leaflet_data: "np.ndarray",
        exclude_sel: Optional[str],
        exclude_cutoff: float,
        exclude_dim: int,
        start: Optional[int],
        stop: Optional[int],
        step: Optional[int],
        verbose: bool,
    ) -> "np.ndarray":
        """Run official ``lipyphilic.AreaPerLipid`` and return (n_lipids, n_frames)."""
        import numpy as np
        from lipyphilic.analysis.area_per_lipid import AreaPerLipid

        excl = (exclude_sel or "").strip()
        kwargs: Dict[str, Any] = {
            "universe": self.universe,
            "lipid_sel": lipid_sel,
            "leaflets": leaflet_data,
        }
        if excl:
            if not _lipyphilic_areaperlipid_supports_exclude():
                raise ValueError(
                    f"apl_method='lipyphilic' with exclude_sel={exclude_sel!r} requires "
                    "lipyphilic AreaPerLipid exclude support (PR #164). "
                    f"Install with: {_LIPYPHILIC_GIT_INSTALL} — "
                    "or use apl_method='fatslim' (default) or 'evapl' (experimental)."
                )
            kwargs["exclude_sel"] = excl
            kwargs["exclude_cutoff"] = float(exclude_cutoff)
            kwargs["exclude_dim"] = int(exclude_dim)
        elif _lipyphilic_areaperlipid_supports_exclude():
            kwargs["exclude_sel"] = ""
            kwargs["exclude_cutoff"] = float(exclude_cutoff)
            kwargs["exclude_dim"] = int(exclude_dim)

        apl = AreaPerLipid(**kwargs)
        apl.run(start=start, stop=stop, step=step, verbose=verbose)
        areas = np.asarray(_analysis_result(apl, "areas"), dtype=float)
        if areas.ndim != 2:
            raise RuntimeError(
                f"lipyphilic AreaPerLipid returned unexpected areas shape {areas.shape}"
            )
        return areas

    def _finalize_area_per_lipid_result(
        self,
        *,
        area_array: "np.ndarray",
        leaflet_data: "np.ndarray",
        lipid_sel: str,
        exclude_sel: Optional[str],
        exclude_cutoff: float,
        resolved_method: str,
        n_frames: int,
        start: Optional[int],
        stop: Optional[int],
        step: Optional[int],
    ) -> Dict[str, Any]:
        """Shared stats / logging / return dict for all APL backends."""
        import numpy as np

        leaflet_means = _leaflet_means_per_frame(area_array, leaflet_data)
        # GridMAT (GW + .pl): headline mean matches Perl Ave APL aggregation.
        if resolved_method in {"gridmat", "gridmat_md"}:
            mean_per_frame = 0.5 * (
                np.asarray(leaflet_means["upper"], dtype=float)
                + np.asarray(leaflet_means["lower"], dtype=float)
            )
        else:
            mean_per_frame = np.nanmean(area_array, axis=0)
        time_ns = self._align_time_ns(n_frames, start=start, stop=stop, step=step)
        metadata = self._lipid_residue_metadata(lipid_sel)
        n_lipids = int(area_array.shape[0])
        u = self.universe

        box_areas = []
        sample_idx = list(dict.fromkeys([0, max(0, n_frames // 2), max(0, n_frames - 1)]))
        for fi in sample_idx:
            try:
                ts_i = (start or 0) + fi * (step or 1)
                u.trajectory[ts_i]
                sample_dims = u.dimensions
                if sample_dims is not None and len(sample_dims) >= 2:
                    box_areas.append(float(sample_dims[0]) * float(sample_dims[1]))
            except Exception:
                continue
        mean_apl = float(np.nanmean(mean_per_frame)) if n_frames else float("nan")
        std_apl = float(np.nanstd(mean_per_frame)) if n_frames else float("nan")
        min_apl = float(np.nanmin(mean_per_frame)) if n_frames else float("nan")
        max_apl = float(np.nanmax(mean_per_frame)) if n_frames else float("nan")
        box_mean = float(np.mean(box_areas)) if box_areas else float("nan")
        expected = (
            (2.0 * box_mean / n_lipids) if n_lipids and np.isfinite(box_mean) else float("nan")
        )
        nearly_flat = bool(
            np.isfinite(std_apl) and std_apl < max(1e-6, 0.01 * abs(mean_apl))
        )
        excl_note = (
            f" method={resolved_method} exclude={exclude_sel!r}"
            + (
                f" cutoff={exclude_cutoff}"
                if resolved_method not in {"evapl", "vtmc", "fatslim"}
                else ""
            )
            if (exclude_sel or "").strip()
            else f" method={resolved_method} exclude=off"
        )
        logger.info(
            "Area per lipid: selection=%r n_lipids=%d n_frames=%d "
            "mean=%.4f std=%.4f min=%.4f max=%.4f Å² | "
            "sample_box_XY=%.2f Å² expected_mean≈box/n=%.4f Å² | "
            "flat_mean=%s%s",
            lipid_sel,
            n_lipids,
            n_frames,
            mean_apl,
            std_apl,
            min_apl,
            max_apl,
            box_mean,
            expected,
            nearly_flat,
            excl_note,
        )
        if nearly_flat and not (exclude_sel or "").strip():
            logger.info(
                "Area per lipid mean is nearly constant across frames. "
                "Without exclude sites, Voronoi lipid areas tile the periodic XY box, so "
                "mean ≈ Lx·Ly / n_lipids. Fixed-area ensembles (NVT, some NPAT) "
                "yield a flat mean while upper/lower leaflet means can still fluctuate."
            )

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
        if box_z <= 1.0:
            raise ValueError(
                "Trajectory has no periodic box (Box is None). Membrane thickness "
                "needs unit-cell dimensions from DCD/XTC/TRR. Remove starting PDB/GRO "
                "files from the trajectory list; to RMSD against a starting structure, "
                "use the RMSD reference PDB field instead."
            )
        thickness = _correct_pbc_straddling_thickness(thickness, box_z)

        n_frames = thickness.size
        time_ns = self._align_time_ns(n_frames, start=start, stop=stop, step=step)

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
    exclude_sel: Optional[str] = "protein",
    exclude_cutoff: float = 30.0,
    exclude_dim: int = 3,
    apl_method: Optional[str] = "fatslim",
    gridmat_n: int = 20,
    gridmat_precision: float = 13.0,
    vtmc_n_samples: int = 50_000,
    vtmc_protein_radius: float = 1.7,
    fatslim_nthreads: int = 1,
    fatslim_jobs: int = 1,
    gridmat_md_jobs: Optional[int] = None,
    file_times: Optional[Dict[str, float]] = None,
    file_strides: Optional[Dict[str, int]] = None,
    start: Optional[int] = None,
    stop: Optional[int] = None,
    step: Optional[int] = None,
    verbose: bool = False,
) -> Dict[str, Any]:
    """
    Run lipid bilayer analysis and return JSON-serializable arrays.

    Supported analysis types: ``area_per_lipid``, ``membrane_thickness``.
    """
    import gc
    import numpy as np
    from gatewizard.utils.trajectory_analysis import prepare_structural_inputs, _lookup_file_map

    top = Path(topology_file).expanduser().resolve()
    trajs, _ = prepare_structural_inputs(
        trajectory_files, analysis_type=analysis_type
    )
    effective_step = step
    if file_strides and effective_step is None:
        strides = [
            max(1, int(_lookup_file_map(file_strides, p) or 1)) for p in trajs
        ]
        # Per-file stride is applied via in-memory subsampling on TrajectoryAnalyzer;
        # only use a global step when there is a single file and no stride map yet.
        if len(trajs) == 1 and len(set(strides)) == 1 and strides[0] > 1:
            effective_step = strides[0]
    analyzer = BilayerTrajectoryAnalyzer(
        top, trajs, file_times=file_times, file_strides=file_strides
    )
    if analyzer._trajectory._uses_stride():
        effective_step = 1

    try:
        atype = analysis_type.strip().lower().replace(" ", "_").replace("-", "_")
        if atype in {"area_per_lipid", "apl"}:
            data = analyzer.calculate_area_per_lipid(
                lipid_sel=lipid_sel,
                leaflet_lipid_sel=leaflet_lipid_sel,
                exclude_sel=exclude_sel,
                exclude_cutoff=exclude_cutoff,
                exclude_dim=exclude_dim,
                apl_method=apl_method,
                gridmat_n=gridmat_n,
                gridmat_precision=gridmat_precision,
                vtmc_n_samples=vtmc_n_samples,
                vtmc_protein_radius=vtmc_protein_radius,
                fatslim_nthreads=fatslim_nthreads,
                fatslim_jobs=fatslim_jobs,
                gridmat_md_jobs=gridmat_md_jobs,
                start=start,
                stop=stop,
                step=effective_step,
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
                step=effective_step,
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
    finally:
        analyzer._trajectory.clear_analysis_cache()
        del analyzer
        gc.collect()
