"""
Cavity hydration via standalone PACKMOL (AmberTools).

Estimate free volume inside a user-defined box, build PACKMOL input, and fill
internal cavities with TIP3P waters.  Designed for heavy-atom-only structures
from Visualize using inflated exclusion radii.
"""

from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from gatewizard.core.structure_manager import Atom, parse_pdb
from gatewizard.utils.logger import get_logger
from gatewizard.utils.optional_deps import (
    get_external_tool_versions,
    probe_executable_version,
    resolve_executable,
)

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

WATER_NUMBER_DENSITY_A3 = 1.0 / 30.0  # ~1 molecule per 30 Å³
MAX_GRID_POINTS = 2_000_000
DEFAULT_GRID_SPACING = 0.5
DEFAULT_TOLERANCE = 2.0
DEFAULT_NLOOP = 20
DEFAULT_SOLUTE_RADIUS_HEAVY = 2.5
DEFAULT_SOLUTE_RADIUS_EXPLICIT = 1.5
PACKING_EFFICIENCY_HEAVY = 0.85
PACKING_EFFICIENCY_EXPLICIT = 1.0

# Bondi-like VdW radii (Å)
_VDW_EXPLICIT: Dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "P": 1.80,
    "F": 1.47,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "FE": 1.80,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 2.00,
    "NA": 2.27,
    "K": 2.75,
}

# Implicit-H envelope added in heavy-atom-safe mode
_HEAVY_INFLATION: Dict[str, float] = {
    "C": 0.5,
    "N": 1.0,
    "O": 1.0,
    "S": 1.0,
    "P": 1.0,
}

PROTEIN_RESIDUES = frozenset(
    {
        "ALA",
        "ARG",
        "ASN",
        "ASP",
        "CYS",
        "GLN",
        "GLU",
        "GLY",
        "HIS",
        "ILE",
        "LEU",
        "LYS",
        "MET",
        "PHE",
        "PRO",
        "SER",
        "THR",
        "TRP",
        "TYR",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "CYX",
        "ASH",
        "GLH",
        "LYN",
        "ARN",
        "TYM",
        "CYM",
        "HSD",
        "HSE",
        "HSP",
    }
)

HydrogenStatus = str  # "full" | "partial" | "none"
ExclusionMode = str  # "heavy_atom_safe" | "explicit"

BoxTuple = Tuple[float, float, float]


class PackmolHydrationError(Exception):
    """Raised when hydration preparation or execution fails."""


@dataclass
class VolumeEstimate:
    box_volume_A3: float
    occupied_volume_A3: float
    free_volume_A3: float
    suggested_waters: int
    hydrogen_status: HydrogenStatus
    exclusion_mode: ExclusionMode
    packing_efficiency: float
    grid_spacing_used: float

    def as_dict(self) -> Dict[str, Any]:
        return {
            "box_volume_A3": self.box_volume_A3,
            "occupied_volume_A3": self.occupied_volume_A3,
            "free_volume_A3": self.free_volume_A3,
            "suggested_waters": self.suggested_waters,
            "hydrogen_status": self.hydrogen_status,
            "exclusion_mode": self.exclusion_mode,
            "packing_efficiency": self.packing_efficiency,
            "grid_spacing_used": self.grid_spacing_used,
        }


@dataclass
class HydrationJobResult:
    job_dir: str
    output_pdb: str
    packmol_log: str
    packmol_inp_path: str
    volumes: VolumeEstimate
    hydrogen_status: HydrogenStatus
    exclusion_mode: ExclusionMode
    success: bool
    message: str

    def as_dict(self) -> Dict[str, Any]:
        return {
            "job_dir": self.job_dir,
            "output_pdb": self.output_pdb,
            "packmol_log": self.packmol_log,
            "packmol_inp_path": self.packmol_inp_path,
            "volumes": self.volumes.as_dict(),
            "hydrogen_status": self.hydrogen_status,
            "exclusion_mode": self.exclusion_mode,
            "success": self.success,
            "message": self.message,
        }


# ---------------------------------------------------------------------------
# Resource paths
# ---------------------------------------------------------------------------


def _bundled_tip3p_path() -> Path:
    return Path(__file__).resolve().parent.parent / "resources" / "water" / "TIP3P.pdb"


def _resolve_tip3p_template() -> Path:
    bundled = _bundled_tip3p_path()
    if bundled.is_file():
        return bundled
    for env_var in ("AMBERHOME", "CONDA_PREFIX"):
        root = os.environ.get(env_var)
        if not root:
            continue
        for rel in (
            "dat/leap/pdb/TIP3P.pdb",
            "dat/leap/pdb/tip3p.pdb",
        ):
            candidate = Path(root) / rel
            if candidate.is_file():
                return candidate
    raise PackmolHydrationError(
        "TIP3P template not found. Expected bundled file or AmberTools dat/leap/pdb/TIP3P.pdb"
    )


# ---------------------------------------------------------------------------
# Public API — availability
# ---------------------------------------------------------------------------


def check_packmol_available() -> Dict[str, Any]:
    """Return PACKMOL availability, version, and resolved executable path."""
    path = resolve_executable(("packmol",))
    version = None
    if path:
        version = probe_executable_version(path, "packmol")
    return {
        "available": path is not None,
        "version": version,
        "resolved_path": path,
    }


# ---------------------------------------------------------------------------
# Hydrogen detection
# ---------------------------------------------------------------------------


def _is_hydrogen(atom: Atom) -> bool:
    el = (atom.element or "").upper()
    if el == "H":
        return True
    name = atom.name.strip().upper()
    return name.startswith("H") and el in ("", "H")


def _is_protein_heavy(atom: Atom) -> bool:
    res = (atom.res_name or "").upper()
    if res in PROTEIN_RESIDUES:
        return not _is_hydrogen(atom)
    return False


def detect_hydrogen_status(
    pdb_file: Optional[str] = None,
    atoms: Optional[Sequence[Atom]] = None,
) -> HydrogenStatus:
    """
    Classify hydrogen content on protein atoms.

    Returns ``full``, ``partial``, or ``none``.
    """
    if atoms is None:
        if not pdb_file:
            raise ValueError("Either pdb_file or atoms must be provided")
        atoms = parse_pdb(pdb_file).atoms

    protein_heavy = 0
    protein_h = 0
    for atom in atoms:
        res = (atom.res_name or "").upper()
        if res not in PROTEIN_RESIDUES:
            continue
        if _is_hydrogen(atom):
            protein_h += 1
        else:
            protein_heavy += 1

    if protein_heavy == 0:
        return "none"
    ratio = protein_h / protein_heavy
    if ratio >= 0.8:
        return "full"
    if ratio >= 0.05:
        return "partial"
    return "none"


def _default_exclusion_mode(
    hydrogen_status: HydrogenStatus,
    exclusion_mode: Optional[ExclusionMode],
) -> ExclusionMode:
    if exclusion_mode in ("heavy_atom_safe", "explicit"):
        return exclusion_mode  # type: ignore[return-value]
    if hydrogen_status == "full":
        return "explicit"
    return "heavy_atom_safe"


def _default_solute_radius(
    exclusion_mode: ExclusionMode,
    solute_radius: Optional[float],
) -> float:
    if solute_radius is not None and solute_radius > 0:
        return float(solute_radius)
    if exclusion_mode == "explicit":
        return DEFAULT_SOLUTE_RADIUS_EXPLICIT
    return DEFAULT_SOLUTE_RADIUS_HEAVY


def _packing_efficiency(exclusion_mode: ExclusionMode) -> float:
    if exclusion_mode == "explicit":
        return PACKING_EFFICIENCY_EXPLICIT
    return PACKING_EFFICIENCY_HEAVY


# ---------------------------------------------------------------------------
# Radii helpers
# ---------------------------------------------------------------------------


def _normalize_element(atom: Atom) -> str:
    el = (atom.element or "").strip().upper()
    if el:
        return el
    name = atom.name.strip()
    if name:
        letter = name[0].upper()
        if letter.isalpha():
            return letter
    return "C"


def vdw_radius(element: str, exclusion_mode: ExclusionMode = "explicit") -> float:
    """Return VdW radius for an element, optionally inflated for heavy-atom-safe mode."""
    el = element.upper()
    base = _VDW_EXPLICIT.get(el, 1.70)
    if exclusion_mode == "heavy_atom_safe":
        base += _HEAVY_INFLATION.get(el, 0.3)
    return base


def exclusion_radius(
    atom: Atom,
    exclusion_mode: ExclusionMode,
    solute_radius: float,
) -> float:
    """Effective exclusion sphere radius for grid occupancy tests."""
    return vdw_radius(_normalize_element(atom), exclusion_mode) + solute_radius


# ---------------------------------------------------------------------------
# Volume estimation
# ---------------------------------------------------------------------------


def _parse_box(box_min: Sequence[float], box_max: Sequence[float]) -> Tuple[BoxTuple, BoxTuple]:
    bmin = (float(box_min[0]), float(box_min[1]), float(box_min[2]))
    bmax = (float(box_max[0]), float(box_max[1]), float(box_max[2]))
    for i in range(3):
        if bmin[i] >= bmax[i]:
            raise ValueError(
                f"Invalid box: min[{i}]={bmin[i]} must be less than max[{i}]={bmax[i]}"
            )
    return bmin, bmax


def _filter_atoms(
    atoms: Sequence[Atom],
    atom_indices: Optional[Sequence[int]],
) -> List[Atom]:
    if atom_indices is None:
        return list(atoms)
    index_set = set(int(i) for i in atom_indices)
    return [a for i, a in enumerate(atoms) if i in index_set]


def _choose_grid_spacing(
    bmin: BoxTuple,
    bmax: BoxTuple,
    grid_spacing: Optional[float],
) -> float:
    spacing = float(grid_spacing or DEFAULT_GRID_SPACING)
    if spacing <= 0:
        raise ValueError("grid_spacing must be positive")

    nx = max(1, int(math.ceil((bmax[0] - bmin[0]) / spacing)))
    ny = max(1, int(math.ceil((bmax[1] - bmin[1]) / spacing)))
    nz = max(1, int(math.ceil((bmax[2] - bmin[2]) / spacing)))
    points = nx * ny * nz
    if points <= MAX_GRID_POINTS:
        return spacing

    scale = (points / MAX_GRID_POINTS) ** (1.0 / 3.0)
    coarsened = spacing * scale
    logger.warning(
        "Grid would exceed %d points; coarsening spacing %.3f → %.3f Å",
        MAX_GRID_POINTS,
        spacing,
        coarsened,
    )
    return coarsened


def _atoms_near_box(
    atoms: Sequence[Atom],
    atom_radii: Sequence[float],
    bmin: BoxTuple,
    bmax: BoxTuple,
) -> Tuple[List[Tuple[float, float, float]], List[float]]:
    """Keep atoms whose exclusion spheres may intersect the box."""
    max_r = max(atom_radii) if atom_radii else 0.0
    kept_coords: List[Tuple[float, float, float]] = []
    kept_radii: List[float] = []
    for atom, ar in zip(atoms, atom_radii):
        x, y, z = atom.coord
        if (
            x < bmin[0] - max_r
            or x > bmax[0] + max_r
            or y < bmin[1] - max_r
            or y > bmax[1] + max_r
            or z < bmin[2] - max_r
            or z > bmax[2] + max_r
        ):
            continue
        kept_coords.append((x, y, z))
        kept_radii.append(ar)
    return kept_coords, kept_radii


def _count_free_cells_numpy(
    bmin: BoxTuple,
    bmax: BoxTuple,
    spacing: float,
    atom_coords: Sequence[Tuple[float, float, float]],
    atom_radii: Sequence[float],
) -> Tuple[int, int]:
    import numpy as np

    nx = max(1, int(math.ceil((bmax[0] - bmin[0]) / spacing)))
    ny = max(1, int(math.ceil((bmax[1] - bmin[1]) / spacing)))
    nz = max(1, int(math.ceil((bmax[2] - bmin[2]) / spacing)))
    total_cells = nx * ny * nz

    if not atom_coords:
        return total_cells, total_cells

    xs = bmin[0] + (np.arange(nx) + 0.5) * spacing
    ys = bmin[1] + (np.arange(ny) + 0.5) * spacing
    zs = bmin[2] + (np.arange(nz) + 0.5) * spacing
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    points = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)

    occupied = np.zeros(points.shape[0], dtype=bool)
    coords = np.array(atom_coords, dtype=float)
    radii = np.array(atom_radii, dtype=float)
    for i in range(len(coords)):
        d = points - coords[i]
        dist2 = np.einsum("ij,ij->i", d, d)
        occupied |= dist2 <= radii[i] * radii[i]

    free_cells = int(total_cells - occupied.sum())
    return free_cells, total_cells


def _count_free_cells_python(
    bmin: BoxTuple,
    bmax: BoxTuple,
    spacing: float,
    atom_coords: Sequence[Tuple[float, float, float]],
    atom_radii: Sequence[float],
) -> Tuple[int, int]:
    nx = max(1, int(math.ceil((bmax[0] - bmin[0]) / spacing)))
    ny = max(1, int(math.ceil((bmax[1] - bmin[1]) / spacing)))
    nz = max(1, int(math.ceil((bmax[2] - bmin[2]) / spacing)))
    total_cells = nx * ny * nz
    free_cells = 0

    for ix in range(nx):
        x = bmin[0] + (ix + 0.5) * spacing
        for iy in range(ny):
            y = bmin[1] + (iy + 0.5) * spacing
            for iz in range(nz):
                z = bmin[2] + (iz + 0.5) * spacing
                occupied = False
                for (ax, ay, az), ar in zip(atom_coords, atom_radii):
                    dx, dy, dz = x - ax, y - ay, z - az
                    if dx * dx + dy * dy + dz * dz <= ar * ar:
                        occupied = True
                        break
                if not occupied:
                    free_cells += 1
    return free_cells, total_cells


def estimate_cavity_volume(
    pdb_file: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    solute_radius: Optional[float] = None,
    exclusion_mode: Optional[ExclusionMode] = None,
    grid_spacing: Optional[float] = None,
    atom_indices: Optional[Sequence[int]] = None,
) -> VolumeEstimate:
    """
    Estimate free volume inside a box using a regular grid and atom exclusion spheres.
    """
    bmin, bmax = _parse_box(box_min, box_max)
    structure = parse_pdb(pdb_file)
    atoms = _filter_atoms(structure.atoms, atom_indices)

    hydrogen_status = detect_hydrogen_status(atoms=atoms)
    mode = _default_exclusion_mode(hydrogen_status, exclusion_mode)
    radius = _default_solute_radius(mode, solute_radius)
    spacing = _choose_grid_spacing(bmin, bmax, grid_spacing)
    efficiency = _packing_efficiency(mode)

    box_volume = (
        (bmax[0] - bmin[0]) * (bmax[1] - bmin[1]) * (bmax[2] - bmin[2])
    )

    # Precompute atom exclusion radii
    atom_radii = [exclusion_radius(a, mode, radius) for a in atoms]
    near_coords, near_radii = _atoms_near_box(atoms, atom_radii, bmin, bmax)

    cell_volume = spacing ** 3
    try:
        free_cells, _total = _count_free_cells_numpy(
            bmin, bmax, spacing, near_coords, near_radii
        )
    except Exception:
        free_cells, _total = _count_free_cells_python(
            bmin, bmax, spacing, near_coords, near_radii
        )

    free_volume = free_cells * cell_volume
    occupied_volume = box_volume - free_volume
    suggested = max(
        0,
        int(free_volume * WATER_NUMBER_DENSITY_A3 * efficiency),
    )

    return VolumeEstimate(
        box_volume_A3=box_volume,
        occupied_volume_A3=occupied_volume,
        free_volume_A3=free_volume,
        suggested_waters=suggested,
        hydrogen_status=hydrogen_status,
        exclusion_mode=mode,
        packing_efficiency=efficiency,
        grid_spacing_used=spacing,
    )


# ---------------------------------------------------------------------------
# PACKMOL input
# ---------------------------------------------------------------------------


def build_hydrate_inp_text(
    protein_path: str,
    tip3p_path: str,
    output_pdb: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: int,
    solute_radius: float = DEFAULT_SOLUTE_RADIUS_HEAVY,
    tolerance: float = DEFAULT_TOLERANCE,
    nloop: int = DEFAULT_NLOOP,
) -> str:
    """Build PACKMOL input text for cavity hydration."""
    bmin, bmax = _parse_box(box_min, box_max)
    if n_waters < 1:
        raise ValueError("n_waters must be at least 1")

    lines = [
        f"tolerance {tolerance}",
        "filetype pdb",
        f"output {output_pdb}",
        "",
        f"structure {protein_path}",
        "  number 1",
        f"  radius {solute_radius}",
        "end structure",
        "",
        f"structure {tip3p_path}",
        f"  nloop {nloop}",
        f"  number {n_waters}",
        "  inside box "
        f"{bmin[0]:.3f} {bmin[1]:.3f} {bmin[2]:.3f} "
        f"{bmax[0]:.3f} {bmax[1]:.3f} {bmax[2]:.3f}",
        "end structure",
        "",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Job preparation
# ---------------------------------------------------------------------------


def prepare_hydration_job(
    pdb_file: str,
    job_dir: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: int,
    solute_radius: Optional[float] = None,
    exclusion_mode: Optional[ExclusionMode] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    nloop: int = DEFAULT_NLOOP,
    grid_spacing: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Create job directory, copy PDB and TIP3P template, write packmol.inp.

    Returns paths and volume estimate metadata.
    """
    job_path = Path(job_dir)
    job_path.mkdir(parents=True, exist_ok=True)

    pdb_src = Path(pdb_file)
    if not pdb_src.is_file():
        raise PackmolHydrationError(f"PDB file not found: {pdb_file}")

    protein_name = pdb_src.name
    protein_dest = job_path / protein_name
    shutil.copy2(pdb_src, protein_dest)

    tip3p_src = _resolve_tip3p_template()
    tip3p_dest = job_path / "TIP3P.pdb"
    shutil.copy2(tip3p_src, tip3p_dest)

    volumes = estimate_cavity_volume(
        pdb_file=str(protein_dest),
        box_min=box_min,
        box_max=box_max,
        solute_radius=solute_radius,
        exclusion_mode=exclusion_mode,
        grid_spacing=grid_spacing,
    )
    radius = _default_solute_radius(volumes.exclusion_mode, solute_radius)
    stem = pdb_src.stem
    output_name = f"{stem}_hydrated.pdb"
    output_path = job_path / output_name

    inp_text = build_hydrate_inp_text(
        protein_path=str(protein_dest.resolve()),
        tip3p_path=str(tip3p_dest.resolve()),
        output_pdb=str(output_path.resolve()),
        box_min=box_min,
        box_max=box_max,
        n_waters=n_waters,
        solute_radius=radius,
        tolerance=tolerance,
        nloop=nloop,
    )
    inp_path = job_path / "packmol.inp"
    inp_path.write_text(inp_text, encoding="utf-8")

    return {
        "job_dir": str(job_path.resolve()),
        "protein_path": str(protein_dest.resolve()),
        "tip3p_path": str(tip3p_dest.resolve()),
        "packmol_inp_path": str(inp_path.resolve()),
        "output_pdb_name": output_name,
        "output_pdb_path": str(output_path.resolve()),
        "inp_text": inp_text,
        "volumes": volumes.as_dict(),
        "hydrogen_status": volumes.hydrogen_status,
        "exclusion_mode": volumes.exclusion_mode,
    }


def preview_hydrate_inp(
    pdb_file: str,
    job_dir: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: int,
    solute_radius: Optional[float] = None,
    exclusion_mode: Optional[ExclusionMode] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    nloop: int = DEFAULT_NLOOP,
    grid_spacing: Optional[float] = None,
) -> Dict[str, Any]:
    """Prepare job files without running PACKMOL (alias for prepare_hydration_job output shape)."""
    result = prepare_hydration_job(
        pdb_file=pdb_file,
        job_dir=job_dir,
        box_min=box_min,
        box_max=box_max,
        n_waters=n_waters,
        solute_radius=solute_radius,
        exclusion_mode=exclusion_mode,
        tolerance=tolerance,
        nloop=nloop,
        grid_spacing=grid_spacing,
    )
    return {
        "inp_text": result["inp_text"],
        "job_dir_relative_paths": {
            "protein": Path(result["protein_path"]).name,
            "tip3p": "TIP3P.pdb",
            "inp": "packmol.inp",
            "output": result["output_pdb_name"],
        },
        "output_pdb_name": result["output_pdb_name"],
        "exclusion_mode": result["exclusion_mode"],
        "volumes": result["volumes"],
        "hydrogen_status": result["hydrogen_status"],
        "packmol_inp_path": result["packmol_inp_path"],
    }


# ---------------------------------------------------------------------------
# PACKMOL execution
# ---------------------------------------------------------------------------


def run_packmol(
    inp_path: str,
    cwd: Optional[str] = None,
    timeout: int = 600,
) -> Tuple[bool, str]:
    """Run PACKMOL on an input file. Returns (success, log_text)."""
    packmol = resolve_executable(("packmol",))
    if not packmol:
        raise PackmolHydrationError(
            "PACKMOL executable not found. Install AmberTools (conda-forge)."
        )

    work_dir = cwd or str(Path(inp_path).parent)
    cmd = [packmol, "-i", str(Path(inp_path).name)]
    try:
        proc = subprocess.run(
            cmd,
            cwd=work_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        log = (exc.stdout or "") + (exc.stderr or "")
        return False, log + "\nPACKMOL timed out."

    log = (proc.stdout or "") + (proc.stderr or "")
    success = proc.returncode == 0 and "SUCCESS" in log.upper()
    return success, log


def hydrate_cavity(
    pdb_file: str,
    working_dir: str,
    output_folder_name: str,
    box_min: Sequence[float],
    box_max: Sequence[float],
    n_waters: Optional[int] = None,
    solute_radius: Optional[float] = None,
    exclusion_mode: Optional[ExclusionMode] = None,
    tolerance: float = DEFAULT_TOLERANCE,
    nloop: int = DEFAULT_NLOOP,
    grid_spacing: Optional[float] = None,
) -> HydrationJobResult:
    """Prepare job, run PACKMOL, and return result metadata."""
    if not working_dir:
        raise PackmolHydrationError("working_dir is required")
    if not output_folder_name.strip():
        raise PackmolHydrationError("output_folder_name is required")

    job_dir = str(Path(working_dir) / output_folder_name.strip())

    volumes = estimate_cavity_volume(
        pdb_file=pdb_file,
        box_min=box_min,
        box_max=box_max,
        solute_radius=solute_radius,
        exclusion_mode=exclusion_mode,
        grid_spacing=grid_spacing,
    )
    count = n_waters if n_waters is not None and n_waters > 0 else volumes.suggested_waters
    if count < 1:
        raise PackmolHydrationError(
            "No waters to place (free volume too small or n_waters=0)."
        )

    prep = prepare_hydration_job(
        pdb_file=pdb_file,
        job_dir=job_dir,
        box_min=box_min,
        box_max=box_max,
        n_waters=count,
        solute_radius=solute_radius,
        exclusion_mode=exclusion_mode,
        tolerance=tolerance,
        nloop=nloop,
        grid_spacing=grid_spacing,
    )
    volumes = VolumeEstimate(**prep["volumes"])

    inp_path = prep["packmol_inp_path"]
    success, log = run_packmol(inp_path)
    log_path = Path(job_dir) / "packmol.log"
    log_path.write_text(log, encoding="utf-8")

    output_pdb = prep["output_pdb_path"]
    message = "PACKMOL completed successfully." if success else "PACKMOL failed; see log."
    if success and not Path(output_pdb).is_file():
        success = False
        message = f"PACKMOL reported success but output not found: {output_pdb}"

    return HydrationJobResult(
        job_dir=job_dir,
        output_pdb=output_pdb,
        packmol_log=str(log_path.resolve()),
        packmol_inp_path=inp_path,
        volumes=volumes,
        hydrogen_status=prep["hydrogen_status"],
        exclusion_mode=prep["exclusion_mode"],
        success=success,
        message=message,
    )


def run_custom_packmol(
    inp_text: str,
    working_dir: str,
    output_folder_name: str,
    inp_filename: str = "packmol_custom.inp",
    timeout: int = 600,
) -> Dict[str, Any]:
    """Write user PACKMOL input and execute (Phase 2 custom tab)."""
    if not working_dir:
        raise PackmolHydrationError("working_dir is required")

    job_dir = Path(working_dir) / output_folder_name.strip()
    job_dir.mkdir(parents=True, exist_ok=True)
    inp_path = job_dir / inp_filename
    inp_path.write_text(inp_text, encoding="utf-8")

    success, log = run_packmol(str(inp_path), cwd=str(job_dir), timeout=timeout)
    log_path = job_dir / "packmol.log"
    log_path.write_text(log, encoding="utf-8")

    output_match = re.search(r"^output\s+(\S+)", inp_text, re.MULTILINE | re.IGNORECASE)
    output_pdb = str((job_dir / output_match.group(1)).resolve()) if output_match else ""

    return {
        "job_dir": str(job_dir.resolve()),
        "packmol_inp_path": str(inp_path.resolve()),
        "packmol_log": str(log_path.resolve()),
        "output_pdb": output_pdb,
        "success": success,
        "message": "PACKMOL completed successfully." if success else "PACKMOL failed; see log.",
    }
