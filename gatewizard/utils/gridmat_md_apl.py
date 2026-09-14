"""
Optional GridMAT-MD.pl CLI bridge for area-per-lipid analysis.

GridMAT-MD is GPLv3; GateWizard (MIT) never vendors or imports the ``.pl``
script. This module exports per-frame GRO files, writes a param file, runs
``perl GridMAT-MD.pl``, and parses leaflet ``Ave APL`` values (Å²).

Discovery: ``GATEWIZARD_GRIDMAT_MD`` → ``PATH`` → common clone locations.
Perl: ``GATEWIZARD_PERL`` → ``PATH``.
See ``scripts/install_gridmat_md.sh`` and ``docs/analysis.md``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

GRIDMAT_MD_ENV_VAR = "GATEWIZARD_GRIDMAT_MD"
PERL_ENV_VAR = "GATEWIZARD_PERL"
GRIDMAT_MD_SCRIPT_NAME = "GridMAT-MD.pl"
GRIDMAT_MD_INSTALL_HINT = (
    "GridMAT-MD.pl is required for apl_method='gridmat_md'. "
    "Clone https://github.com/jalemkul/gridmat-md and set "
    "GATEWIZARD_GRIDMAT_MD to the GridMAT-MD.pl path (perl must be on PATH "
    "or set GATEWIZARD_PERL). See scripts/install_gridmat_md.sh and "
    "docs/analysis.md. Example: "
    "bash scripts/install_gridmat_md.sh && "
    "export GATEWIZARD_GRIDMAT_MD=\"$HOME/gridmat-md/GridMAT-MD.pl\". "
    "Or use apl_method='gridmat' (GateWizard in-process, experimental) "
    "or apl_method='fatslim' (default)."
)

_AVE_APL_RE = re.compile(r"Ave APL\s*=\s*([0-9.]+)")

_COMMON_SOLVENT_RESNAMES = (
    "WAT",
    "HOH",
    "SOL",
    "TIP3",
    "TIP3P",
    "TIP4",
    "W",
    "PW",
    "WS",
)
_COMMON_ION_RESNAMES = (
    "K+",
    "Cl-",
    "Na+",
    "NA+",
    "CL-",
    "SOD",
    "CLA",
    "POT",
    "MG",
    "CA",
    "ZN",
)


class GridmatMdError(RuntimeError):
    """Raised when GridMAT-MD.pl is missing or a run fails."""


def _is_readable_file(path: Path) -> bool:
    try:
        return path.is_file()
    except OSError:
        return False


def _common_clone_candidates() -> List[Path]:
    """Likely ``…/gridmat-md/GridMAT-MD.pl`` locations."""
    home = Path.home()
    roots: List[Path] = [
        home / "gridmat-md",
        home / "src" / "gridmat-md",
        home / "opt" / "gridmat-md",
        home / "software" / "gridmat-md",
        home / "vendor" / "gridmat-md",
        Path.cwd() / "vendor" / "gridmat-md",
        Path.cwd() / "gridmat-md",
        Path("/opt/gridmat-md"),
        Path("/usr/local/gridmat-md"),
    ]
    # Sibling of a GateWizard checkout
    try:
        pkg_root = Path(__file__).resolve().parents[2]
        roots.append(pkg_root / "vendor" / "gridmat-md")
        roots.append(pkg_root.parent / "gridmat-md")
    except (IndexError, OSError):
        pass

    out: List[Path] = []
    seen: set[str] = set()
    for root in roots:
        cand = root / GRIDMAT_MD_SCRIPT_NAME
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        out.append(cand)
    return out


def resolve_gridmat_md_script() -> Optional[str]:
    """Return path to ``GridMAT-MD.pl`` or ``None`` if not found."""
    env_path = (os.environ.get(GRIDMAT_MD_ENV_VAR) or "").strip()
    if env_path:
        path = Path(env_path).expanduser()
        if _is_readable_file(path):
            return str(path.resolve())
        # Allow pointing at the clone directory
        if path.is_dir():
            nested = path / GRIDMAT_MD_SCRIPT_NAME
            if _is_readable_file(nested):
                return str(nested.resolve())
        which_env = shutil.which(env_path)
        if which_env and _is_readable_file(Path(which_env)):
            return which_env

    which = shutil.which(GRIDMAT_MD_SCRIPT_NAME)
    if which and _is_readable_file(Path(which)):
        return which

    for cand in _common_clone_candidates():
        if _is_readable_file(cand):
            return str(cand.resolve())
    return None


def resolve_perl_executable() -> Optional[str]:
    """Return path to ``perl`` or ``None`` if not found."""
    env_path = (os.environ.get(PERL_ENV_VAR) or "").strip()
    if env_path:
        path = Path(env_path).expanduser()
        if _is_readable_file(path):
            return str(path.resolve())
        which_env = shutil.which(env_path)
        if which_env:
            return which_env

    return shutil.which("perl")


def is_available() -> bool:
    """True when both GridMAT-MD.pl and perl can be resolved."""
    return resolve_gridmat_md_script() is not None and resolve_perl_executable() is not None


def require_gridmat_md() -> Tuple[str, str]:
    """Return ``(script, perl)`` or raise :class:`GridmatMdError` with install hint."""
    script = resolve_gridmat_md_script()
    perl = resolve_perl_executable()
    missing = []
    if not script:
        missing.append("GridMAT-MD.pl")
    if not perl:
        missing.append("perl")
    if missing:
        raise GridmatMdError(
            f"Missing {', '.join(missing)}. {GRIDMAT_MD_INSTALL_HINT}"
        )
    return script, perl


def normalize_gridmat_md_jobs(jobs: Optional[int]) -> int:
    """Frame workers. ``None`` / ``<=0`` → ``min(8, cpu_count)``."""
    if jobs is None or int(jobs) <= 0:
        return max(1, min(8, os.cpu_count() or 1))
    return int(jobs)


def parse_gridmat_ave_apl(path: Union[str, Path]) -> Optional[float]:
    """Parse ``Ave APL = …`` from a GridMAT-MD areas ``.dat`` file."""
    p = Path(path)
    if not p.is_file():
        return None
    match = _AVE_APL_RE.search(p.read_text(encoding="utf-8", errors="replace"))
    return float(match.group(1)) if match else None


def write_gridmat_param(
    path: Union[str, Path],
    *,
    gro_name: str,
    output_prefix: str,
    lipid_types: Sequence[Tuple[str, str]],
    solvent: str,
    ions: str,
    grid_n: int,
    precision_nm: float,
    protein: bool = True,
) -> None:
    """Write a one-frame GridMAT-MD parameter file (GRO, area only)."""
    types = list(lipid_types)[:3]
    while len(types) < 3:
        types.append(("XX", "__none__"))

    lines = [
        f"coord_file {gro_name}",
        "file_type gro",
        "num_frames 1",
        f"num_lipid_types {len(types)}",
    ]
    for i, (resname, atomname) in enumerate(types, start=1):
        lines.append(f"resname{i} {resname}")
        lines.append(f"atomname{i} {atomname}")
    lines.extend(
        [
            f"solvent {solvent}",
            f"ions {ions}",
            "box_size vectors",
            f"grid {int(grid_n)}",
            "conserve_ratio yes",
            f"protein {'yes' if protein else 'no'}",
            f"precision {float(precision_nm):.4g}",
            "P_value 5.0",
            f"output_prefix {output_prefix}",
            "output_format vector",
            "thickness no",
            "area yes",
            "",
        ]
    )
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def renumber_gro_resids_for_gridmat(
    path: Union[str, Path],
    lipid_types: Sequence[Tuple[str, str]],
) -> int:
    """Give each lipid headgroup a unique GRO residue number.

    GridMAT-MD groups lipids by resid, then does APL = area / n_leaflet.
    Some force-field GRO exports reuse the same resid for every lipid of a
    type, collapsing all headgroups into one reference site.
    """
    p = Path(path)
    lines = p.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        return 0
    try:
        n_atoms = int(lines[1].split()[0])
    except (ValueError, IndexError):
        return 0

    match_set = {(r.strip(), a.strip()) for r, a in lipid_types if a.strip() != "__none__"}
    body = lines[2 : 2 + n_atoms]
    tail = lines[2 + n_atoms :]
    next_resid = 1
    n_renumbered = 0
    new_body: List[str] = []
    for line in body:
        padded = line.ljust(15)
        rname = padded[5:10].strip()
        aname = padded[10:15].strip()
        if (rname, aname) in match_set:
            line = f"{next_resid:5d}{line[5:]}"
            next_resid += 1
            n_renumbered += 1
        new_body.append(line)
    p.write_text("\n".join([lines[0], lines[1], *new_body, *tail, ""]), encoding="utf-8")
    return n_renumbered


def _center_bilayer_in_box(universe, lipid_sel: str) -> None:
    """Translate lipid COM to box center and wrap (GridMAT leaflet ID)."""
    import numpy as np

    lipids = universe.select_atoms(lipid_sel)
    if len(lipids) == 0:
        return
    dims = universe.dimensions
    if dims is None or len(dims) < 3:
        return
    com = lipids.center_of_mass()
    box_center = np.asarray(dims[:3], dtype=float) * 0.5
    shift = box_center - com
    universe.atoms.translate(shift)
    universe.atoms.wrap(inplace=True)


def _infer_lipid_types(universe, lipid_sel: str) -> List[Tuple[str, str]]:
    """Up to three (resname, atomname) pairs from ``lipid_sel``."""
    from collections import Counter

    ag = universe.select_atoms(lipid_sel)
    if len(ag) == 0:
        raise GridmatMdError(f"No atoms match lipid selection {lipid_sel!r}")
    pairs = Counter(zip(ag.resnames.tolist(), ag.names.tolist()))
    ranked = [pair for pair, _ in pairs.most_common(3)]
    return [(str(r), str(a)) for r, a in ranked]


def _infer_solvent_ions(universe) -> Tuple[str, str]:
    """Pick solvent / ion tokens present in the topology (comma-separated)."""
    try:
        resnames = {str(r).strip() for r in universe.residues.resnames}
    except Exception:
        resnames = set()

    solvents = [n for n in _COMMON_SOLVENT_RESNAMES if n in resnames]
    ions = [n for n in _COMMON_ION_RESNAMES if n in resnames]
    solvent = ",".join(solvents) if solvents else "WAT"
    ion_str = ",".join(ions) if ions else "K+,Cl-"
    return solvent, ion_str


def _run_gridmat_perl_job(
    perl: str,
    script: str,
    frame_dir: str,
    frame_i: int,
) -> Tuple[int, Optional[float], Optional[float], str]:
    """Worker: run ``perl GridMAT-MD.pl param`` in a prepared frame directory."""
    result = subprocess.run(
        [perl, script, "param"],
        capture_output=True,
        text=True,
        cwd=frame_dir,
        check=False,
    )
    if result.returncode != 0:
        err = (result.stderr or result.stdout or "").strip()[-400:]
        return frame_i, None, None, f"GridMAT-MD.pl failed frame {frame_i}: {err}"

    work = Path(frame_dir)
    top_apl = parse_gridmat_ave_apl(work / "out.top_areas.dat")
    bot_apl = parse_gridmat_ave_apl(work / "out.bottom_areas.dat")
    if top_apl is None or bot_apl is None:
        tops = list(work.glob("*top_areas.dat"))
        bots = list(work.glob("*bottom_areas.dat"))
        if tops and bots:
            top_apl = parse_gridmat_ave_apl(tops[0])
            bot_apl = parse_gridmat_ave_apl(bots[0])
    if top_apl is None or bot_apl is None:
        return frame_i, None, None, f"GridMAT-MD.pl frame {frame_i}: parse Ave APL failed"
    return frame_i, float(top_apl), float(bot_apl), ""


def run_gridmat_md_apl(
    universe,
    *,
    lipid_sel: str,
    exclude_sel: Optional[str],
    frame_indices: Sequence[int],
    work_dir: Union[str, Path],
    gridmat_n: int = 20,
    gridmat_precision: float = 13.0,
    gridmat_md_jobs: Optional[int] = None,
    leaflet_data: Optional[Any] = None,
    gridmat_script: Optional[str] = None,
    perl_exe: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Export strided frames, run GridMAT-MD.pl, return per-lipid areas.

    Frames are prepared sequentially (GRO + resid renumber + param). Only the
    ``perl`` invocations run in a :class:`ProcessPoolExecutor`.

    ``gridmat_precision`` is in Å (GateWizard convention); written to the param
    file in nm (÷10), matching GridMAT-MD.pl.

    Returns:
        Dict with ``areas`` (n_lipids, n_frames) in Å² aligned to
        ``universe.select_atoms(lipid_sel).residues``, plus ``meta``.
        Per-lipid values are leaflet Ave APL (or frame mean when leaflets are
        unavailable) — GridMAT-MD.pl does not expose per-lipid areas in the
        Ave APL summary files.
    """
    import numpy as np

    script = gridmat_script
    perl = perl_exe
    if not script or not perl:
        resolved_script, resolved_perl = require_gridmat_md()
        script = script or resolved_script
        perl = perl or resolved_perl

    indices = [int(i) for i in frame_indices]
    if not indices:
        raise GridmatMdError("No frames selected for GridMAT-MD APL.")

    membrane = universe.select_atoms(lipid_sel, updating=False)
    if membrane.n_residues == 0:
        raise GridmatMdError(f"No residues match lipid selection {lipid_sel!r}")

    lipid_types = _infer_lipid_types(universe, lipid_sel)
    solvent, ions = _infer_solvent_ions(universe)
    protein = bool((exclude_sel or "").strip())
    precision_nm = float(gridmat_precision) / 10.0

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)

    import MDAnalysis as mda

    n_frames = len(indices)
    frame_dirs: List[Path] = []
    for i, idx in enumerate(indices):
        universe.trajectory[idx]
        _center_bilayer_in_box(universe, lipid_sel)
        frame_dir = work / f"frame_{i:05d}"
        frame_dir.mkdir(parents=True, exist_ok=True)
        gro = frame_dir / "frame.gro"
        param = frame_dir / "param"
        with mda.Writer(str(gro), universe.atoms.n_atoms) as w:
            w.write(universe.atoms)
        n_renum = renumber_gro_resids_for_gridmat(gro, lipid_types)
        if i == 0 and n_renum:
            logger.info(
                "GridMAT-MD: renumbered %d lipid headgroup GRO resids (%s)",
                n_renum,
                ", ".join(f"{r}/{a}" for r, a in lipid_types),
            )
        write_gridmat_param(
            param,
            gro_name="frame.gro",
            output_prefix="out",
            lipid_types=lipid_types,
            solvent=solvent,
            ions=ions,
            grid_n=gridmat_n,
            precision_nm=precision_nm,
            protein=protein,
        )
        frame_dirs.append(frame_dir)

    n_jobs = normalize_gridmat_md_jobs(gridmat_md_jobs)
    n_workers = max(1, min(n_jobs, n_frames))
    tops = np.full(n_frames, np.nan, dtype=float)
    bots = np.full(n_frames, np.nan, dtype=float)
    errors: List[str] = []

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        futs = [
            pool.submit(
                _run_gridmat_perl_job,
                perl,
                script,
                str(frame_dirs[i]),
                i,
            )
            for i in range(n_frames)
        ]
        for fut in as_completed(futs):
            frame_i, top_apl, bot_apl, err = fut.result()
            if err:
                errors.append(err)
                continue
            tops[frame_i] = float(top_apl)
            bots[frame_i] = float(bot_apl)

    if errors:
        raise GridmatMdError(
            "; ".join(errors[:3])
            + (f" (+{len(errors) - 3} more)" if len(errors) > 3 else "")
            + f"\n{GRIDMAT_MD_INSTALL_HINT}"
        )

    if np.all(np.isnan(tops)) or np.all(np.isnan(bots)):
        raise GridmatMdError(
            "GridMAT-MD.pl produced no usable Ave APL values. "
            "Check lipid_sel / GRO resid renumbering / solvent+ions tokens."
        )

    means = 0.5 * (tops + bots)
    area_array = np.full((membrane.n_residues, n_frames), np.nan, dtype=float)

    leaflets = None
    if leaflet_data is not None:
        leaflets = np.asarray(leaflet_data)
        if leaflets.ndim == 1:
            leaflets = np.broadcast_to(leaflets[:, None], (leaflets.shape[0], n_frames))

    if leaflets is not None and leaflets.shape[0] == membrane.n_residues:
        for fi in range(n_frames):
            signs = leaflets[:, fi]
            upper = signs > 0
            lower = signs < 0
            if np.any(upper):
                area_array[upper, fi] = tops[fi]
            if np.any(lower):
                area_array[lower, fi] = bots[fi]
            neither = ~(upper | lower)
            if np.any(neither):
                area_array[neither, fi] = means[fi]
    else:
        area_array[:, :] = means[None, :]

    return {
        "areas": area_array,
        "meta": {
            "source": "gridmat_md_cli",
            "script": script,
            "perl": perl,
            "n_frames_requested": n_frames,
            "n_frames_parsed": int(np.sum(np.isfinite(means))),
            "centered_bilayer": True,
            "grid": int(gridmat_n),
            "precision_angstrom": float(gridmat_precision),
            "precision_nm": precision_nm,
            "lipid_types": lipid_types,
            "solvent": solvent,
            "ions": ions,
            "protein": protein,
            "exclude_sel": (exclude_sel or "").strip() or None,
            "work_dir": str(work),
            "jobs": n_workers,
            "upper_ave_apl": tops.tolist(),
            "lower_ave_apl": bots.tolist(),
            "mean_ave_apl": means.tolist(),
        },
    }
