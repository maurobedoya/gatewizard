"""
Optional FATSLiM CLI bridge for area-per-lipid analysis.

FATSLiM is GPLv3 and unmaintained; GateWizard (MIT) never imports it.
This module exports GRO/XTC/NDX with GateWizard's MDAnalysis, runs the
external ``fatslim apl`` binary, and parses ``apl_raw*.csv`` (nm² → Å²).

Discovery order: ``GATEWIZARD_FATSLIM`` → ``PATH`` → conda env ``fatslim-py38``.
See ``scripts/install_fatslim_env.sh`` and ``docs/analysis.md``.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

FATSLIM_ENV_VAR = "GATEWIZARD_FATSLIM"
FATSLIM_CONDA_ENV = "fatslim-py38"
FATSLIM_INSTALL_HINT = (
    "FATSLiM is required for apl_method='fatslim' (GateWizard default). "
    "Install a companion conda env (Python ≤3.8) — do not install into "
    "GateWizard's Python — then set GATEWIZARD_FATSLIM to the fatslim binary. "
    "See scripts/install_fatslim_env.sh and docs/analysis.md "
    "(API / GUI / WSL). Example: "
    "conda create -n fatslim-py38 python=3.8 -y && "
    "conda activate fatslim-py38 && "
    "pip install 'numpy<1.24' 'cython<3' fatslim 'gsd<3.0' 'MDAnalysis==2.4.3' && "
    "export GATEWIZARD_FATSLIM=\"$(which fatslim)\". "
    "Or choose apl_method='evapl' (experimental / not yet validated), "
    "'lipyphilic', 'gridmat', or 'vtmc'."
)


class FatslimError(RuntimeError):
    """Raised when the FATSLiM CLI is missing or a run fails."""


def _is_executable(path: Path) -> bool:
    if not path.is_file():
        return False
    if os.name == "nt":
        return True
    return os.access(path, os.X_OK)


def _conda_env_bin_candidates(env_name: str) -> List[Path]:
    """Likely ``…/envs/<env_name>/bin/fatslim`` locations."""
    roots: List[Path] = []
    for key in ("MAMBA_ROOT_PREFIX", "CONDA_ROOT", "CONDA_PREFIX"):
        val = os.environ.get(key)
        if not val:
            continue
        p = Path(val)
        roots.append(p)
        if p.name != "envs":
            roots.append(p / "envs")
        # When CONDA_PREFIX is an env itself, sibling envs live under ../envs
        if p.parent.name == "envs":
            roots.append(p.parent)

    for base in (
        Path.home() / "micromamba",
        Path.home() / "miniconda3",
        Path.home() / "anaconda3",
        Path.home() / "miniforge3",
        Path.home() / "mambaforge",
        Path("/opt/conda"),
    ):
        roots.append(base)
        roots.append(base / "envs")

    seen: set[str] = set()
    out: List[Path] = []
    for root in roots:
        try:
            if root.name == "envs" and root.is_dir():
                cand = root / env_name / "bin" / "fatslim"
            elif (root / "envs").is_dir():
                cand = root / "envs" / env_name / "bin" / "fatslim"
            else:
                continue
        except OSError:
            continue
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        out.append(cand)
    return out


def resolve_fatslim_executable() -> Optional[str]:
    """Return path to ``fatslim`` or ``None`` if not found."""
    env_path = (os.environ.get(FATSLIM_ENV_VAR) or "").strip()
    if env_path:
        path = Path(env_path).expanduser()
        if _is_executable(path):
            return str(path.resolve())
        which_env = shutil.which(env_path)
        if which_env:
            return which_env

    which = shutil.which("fatslim")
    if which:
        return which

    for cand in _conda_env_bin_candidates(FATSLIM_CONDA_ENV):
        if _is_executable(cand):
            return str(cand.resolve())
    return None


def is_available() -> bool:
    """True when a FATSLiM executable can be resolved."""
    return resolve_fatslim_executable() is not None


def require_fatslim() -> str:
    """Return executable path or raise :class:`FatslimError` with install hint."""
    exe = resolve_fatslim_executable()
    if not exe:
        raise FatslimError(FATSLIM_INSTALL_HINT)
    return exe


def _center_bilayer_in_box(universe, lipid_sel: str) -> None:
    """Translate lipid COM to box center and wrap (FATSLiM leaflet ID)."""
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


def _write_ndx(
    path: Path,
    headgroup_indices_0based: Sequence[int],
    protein_indices_0based: Sequence[int],
    *,
    include_protein: bool,
) -> None:
    def _write_group(fh, title: str, indices: Sequence[int]) -> None:
        fh.write(f"[ {title} ]\n")
        ids = [str(int(i) + 1) for i in indices]
        for i in range(0, len(ids), 15):
            fh.write(" ".join(ids[i : i + 15]) + "\n")

    with path.open("w", encoding="utf-8") as fh:
        _write_group(fh, "headgroups", headgroup_indices_0based)
        if include_protein:
            _write_group(fh, "protein", protein_indices_0based)


def parse_apl_raw_csv(path: Union[str, Path]) -> Dict[int, float]:
    """Parse one FATSLiM ``apl_raw`` CSV → resid → APL in Å²."""
    resid_to_apl: Dict[int, float] = {}
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.lower().startswith("resid"):
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) < 6:
                continue
            try:
                resid = int(float(parts[0]))
                apl_nm2 = float(parts[-1])
            except ValueError:
                continue
            resid_to_apl[resid] = apl_nm2 * 100.0  # nm² → Å²
    return resid_to_apl


def _sorted_frame_csvs(work_dir: Path, pattern: str = "apl_raw*.csv") -> List[Path]:
    files = list(work_dir.glob(pattern))
    if not files:
        return []

    def _key(p: Path) -> Tuple[int, str]:
        m = re.search(r"(\d+)", p.stem)
        return (int(m.group(1)) if m else 0, p.name)

    return sorted(files, key=_key)


def normalize_fatslim_jobs(jobs: Optional[int]) -> int:
    """Frame-chunk workers. ``None`` / ``<=0`` → 1 (single process)."""
    if jobs is None or int(jobs) <= 0:
        return 1
    return int(jobs)


def frame_index_chunks(n_frames: int, n_jobs: int) -> List[Tuple[int, int]]:
    """Split ``[0, n_frames)`` into half-open ``(begin, end)`` for FATSLiM.

    FATSLiM ``--begin-frame`` / ``--end-frame`` use a half-open range
    (``end`` exclusive).
    """
    if n_frames <= 0:
        return []
    workers = max(1, min(int(n_jobs), n_frames))
    base, rem = divmod(n_frames, workers)
    chunks: List[Tuple[int, int]] = []
    start = 0
    for i in range(workers):
        size = base + (1 if i < rem else 0)
        end = start + size
        chunks.append((start, end))
        start = end
    return chunks


def _fill_frame_column(
    area_array,
    frame_i: int,
    mapping: Dict[int, float],
    resid_to_row: Dict[int, int],
) -> None:
    import numpy as np

    matched = 0
    for resid, apl in mapping.items():
        row = resid_to_row.get(int(resid))
        if row is not None:
            area_array[row, frame_i] = apl
            matched += 1
    if matched == 0:
        mean_apl = float(np.mean(list(mapping.values())))
        area_array[:, frame_i] = mean_apl
        logger.warning(
            "FATSLiM resids did not match lipid_sel resids on frame %d; "
            "using frame mean APL=%.2f Å² for all lipids",
            frame_i,
            mean_apl,
        )


def _run_fatslim_chunk(
    exe: str,
    gro: Path,
    ndx: Path,
    xtc: Path,
    begin: int,
    end: int,
    nthreads: int,
    export_prefix: Path,
    *,
    include_protein: bool,
) -> Tuple[int, int, List[Path], str]:
    """Run ``fatslim apl`` on half-open ``[begin, end)`` frames of ``xtc``."""
    export_prefix.parent.mkdir(parents=True, exist_ok=True)
    for old in export_prefix.parent.glob(export_prefix.name + "*.csv"):
        old.unlink()
    cmd = [
        exe,
        "apl",
        "-c",
        str(gro),
        "-n",
        str(ndx),
        "-t",
        str(xtc),
        "--idfreq",
        "0",
        "--begin-frame",
        str(begin),
        "--end-frame",
        str(end),
        "--nthreads",
        str(int(nthreads)),
        "--export-apl-raw",
        str(export_prefix),
    ]
    if include_protein:
        cmd.extend(["--interacting-group", "protein"])

    logger.info("Running FATSLiM: %s", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        tail = (result.stderr or result.stdout or "").strip()[-800:]
        return (
            begin,
            end,
            [],
            f"FATSLiM apl failed frames [{begin},{end}) "
            f"(exit {result.returncode}). {tail}",
        )

    csv_files = _sorted_frame_csvs(export_prefix.parent, f"{export_prefix.name}*.csv")
    expected = end - begin
    if len(csv_files) != expected:
        return (
            begin,
            end,
            csv_files,
            f"FATSLiM frames [{begin},{end}): got {len(csv_files)} CSV files, "
            f"expected {expected}",
        )
    return begin, end, csv_files, ""


def run_fatslim_apl(
    universe,
    *,
    lipid_sel: str,
    exclude_sel: Optional[str],
    frame_indices: Sequence[int],
    work_dir: Union[str, Path],
    fatslim_exe: Optional[str] = None,
    fatslim_nthreads: int = 1,
    fatslim_jobs: int = 1,
) -> Dict[str, Any]:
    """
    Export strided frames, run ``fatslim apl``, return per-lipid areas.

    ``fatslim_nthreads`` is passed to each ``fatslim apl --nthreads`` process
    (``-1`` = all CPUs). ``fatslim_jobs`` splits the exported trajectory into
    parallel frame ranges via ``--begin-frame`` / ``--end-frame``. Prefer a
    small nthreads (1–4) with several jobs over one job with all CPUs.

    Returns:
        Dict with ``areas`` (n_lipids, n_frames) in Å² aligned to
        ``universe.select_atoms(lipid_sel).residues``, plus ``meta``.
    """
    import numpy as np

    exe = fatslim_exe or require_fatslim()
    indices = [int(i) for i in frame_indices]
    if not indices:
        raise FatslimError("No frames selected for FATSLiM APL.")

    membrane = universe.select_atoms(lipid_sel, updating=False)
    if membrane.n_residues == 0:
        raise FatslimError(f"No residues match lipid selection {lipid_sel!r}")
    resids = np.asarray(membrane.residues.resids, dtype=int)

    excl = (exclude_sel or "").strip()
    include_protein = bool(excl)
    prot_indices: List[int] = []
    if include_protein:
        try:
            prot = universe.select_atoms(excl)
        except Exception as exc:
            raise FatslimError(f"Invalid exclude_sel {excl!r}: {exc}") from exc
        if len(prot) == 0:
            raise FatslimError(
                f"exclude_sel {excl!r} matched no atoms for FATSLiM interacting group."
            )
        heavy = prot.select_atoms("not name H*")
        if len(heavy) > 0:
            prot = heavy
        prot_indices = [int(i) for i in prot.indices]

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    for old in work.glob("apl_raw*.csv"):
        old.unlink()
    for chunk_dir in work.glob("chunk_*"):
        if chunk_dir.is_dir():
            shutil.rmtree(chunk_dir)

    gro = work / "frame0.gro"
    xtc = work / "strided.xtc"
    ndx = work / "system.ndx"

    import MDAnalysis as mda

    universe.trajectory[indices[0]]
    _center_bilayer_in_box(universe, lipid_sel)
    with mda.Writer(str(gro), universe.atoms.n_atoms) as w:
        w.write(universe.atoms)

    with mda.Writer(str(xtc), universe.atoms.n_atoms) as w:
        for idx in indices:
            universe.trajectory[idx]
            _center_bilayer_in_box(universe, lipid_sel)
            w.write(universe.atoms)

    # Re-select after writes (indices stable on topology)
    hg = universe.select_atoms(lipid_sel)
    _write_ndx(
        ndx,
        hg.indices,
        prot_indices,
        include_protein=include_protein,
    )

    n_frames = len(indices)
    n_jobs = normalize_fatslim_jobs(fatslim_jobs)
    nthreads = int(fatslim_nthreads)
    chunks = frame_index_chunks(n_frames, n_jobs)
    resid_to_row = {int(r): i for i, r in enumerate(resids.tolist())}
    area_array = np.full((len(resids), n_frames), np.nan, dtype=float)

    errors: List[str] = []
    with ThreadPoolExecutor(max_workers=len(chunks)) as pool:
        futs = []
        for i, (begin, end) in enumerate(chunks):
            chunk_dir = work / f"chunk_{i:03d}"
            chunk_dir.mkdir(parents=True, exist_ok=True)
            export_prefix = chunk_dir / "apl_raw"
            futs.append(
                pool.submit(
                    _run_fatslim_chunk,
                    exe,
                    gro,
                    ndx,
                    xtc,
                    begin,
                    end,
                    nthreads,
                    export_prefix,
                    include_protein=include_protein,
                )
            )
        for fut in as_completed(futs):
            begin, end, csv_files, err = fut.result()
            if err:
                errors.append(err)
                continue
            for local_i, csv_path in enumerate(csv_files):
                frame_i = begin + local_i
                mapping = parse_apl_raw_csv(csv_path)
                if not mapping:
                    errors.append(f"No APL values in {csv_path.name}")
                    continue
                _fill_frame_column(area_array, frame_i, mapping, resid_to_row)

    if errors:
        raise FatslimError(
            "; ".join(errors[:3])
            + (f" (+{len(errors) - 3} more)" if len(errors) > 3 else "")
            + f"\n{FATSLIM_INSTALL_HINT}"
        )

    if np.all(np.isnan(area_array)):
        raise FatslimError(
            "FATSLiM produced no usable APL values. "
            "Check headgroups NDX / bilayer centering / lipid_sel."
        )

    return {
        "areas": area_array,
        "meta": {
            "source": "fatslim_cli",
            "executable": exe,
            "n_frames_requested": n_frames,
            "n_frames_parsed": n_frames,
            "centered_bilayer": True,
            "idfreq": 0,
            "exclude_sel": excl or None,
            "work_dir": str(work),
            "nthreads": nthreads,
            "jobs": len(chunks),
        },
    }
