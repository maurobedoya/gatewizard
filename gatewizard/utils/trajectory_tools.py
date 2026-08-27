"""
Post-processing utilities for MD trajectories.

Engine-aware PBC fixing:
  - GROMACS → ``gmx trjconv`` (whole → nojump → cluster → mol/center)
  - Amber / NAMD / OpenMM → ``cpptraj`` ``autoimage``
  - Fallback → MDAnalysis molecule-aware make_whole + wrap

Jobs write ``status.json`` + ``logs/fix_pbc.log`` for GUI polling
(same shape as preparation jobs).
"""

from __future__ import annotations

import json
import os
import re
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Optional

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

_VALID_ENGINES = frozenset({"auto", "gromacs", "amber", "namd", "openmm", "mdanalysis"})
# Amber lipid21/lipid17 split residues (packmol-memgen POPC → PA/PC/OL) plus
# common CHARMM-style whole-lipid names. Used to default Fix PBC centering to
# protein + membrane for Amber/cpptraj.
_AMBER_LIPID_RESNAMES = frozenset(
    {
        "PA",
        "PC",
        "PE",
        "PS",
        "PG",
        "PI",
        "OL",
        "LA",
        "MY",
        "ST",
        "AR",
        "CHL",
        "CHOL",
        "CHL1",
        "POPC",
        "POPE",
        "POPS",
        "POPG",
        "POPA",
        "POPI",
        "DOPC",
        "DOPE",
        "DOPS",
        "DOPG",
        "DPPC",
        "DPPE",
        "DMPC",
        "DLPC",
        "DSPC",
        "SM",
        "BSM",
        "SPM",
        "PSM",
        "NSM",
        "OLE",
        "PAL",
        "STE",
        "LIN",
        "CHOLEST",
        "PALM",
        "OLEO",
        "STEROL",
    }
)
_AMBER_LIPID_RESNAME_ORDER = (
    "PA",
    "PC",
    "PE",
    "PS",
    "PG",
    "PI",
    "OL",
    "LA",
    "MY",
    "ST",
    "AR",
    "CHL",
    "CHOL",
    "CHL1",
    "POPC",
    "POPE",
    "POPS",
    "POPG",
    "POPA",
    "POPI",
    "DOPC",
    "DOPE",
    "DOPS",
    "DOPG",
    "DPPC",
    "DPPE",
    "DMPC",
    "DLPC",
    "DSPC",
    "SM",
    "BSM",
    "SPM",
    "PSM",
    "NSM",
    "OLE",
    "PAL",
    "STE",
    "LIN",
    "CHOLEST",
    "PALM",
    "OLEO",
    "STEROL",
)
# Compact protein + bilayer default when the topology has not been parsed yet
# (or has no recognised lipids). Extra names match nothing and are harmless.
_DEFAULT_BILAYER_CENTER_RESNAMES = (
    "PA",
    "PC",
    "OL",
    "PE",
    "PS",
    "PG",
    "PI",
    "POPC",
    "POPE",
    "POPS",
    "DPPC",
    "DMPC",
    "DOPC",
    "DSPC",
    "CHL",
    "CHOL",
    "CHL1",
    "SM",
)
_GROMACS_TRAJ = {".xtc", ".trr"}
_AMBER_TRAJ = {".nc", ".ncdf", ".mdcrd", ".crd"}
_DCD_TRAJ = {".dcd"}
_COORD_TRAJ = _GROMACS_TRAJ | _AMBER_TRAJ | _DCD_TRAJ | {".pdb"}


_INTERRUPT_ERROR = "Interrupted — worker process exited unexpectedly"
# Scan can run after status.json exists but before process.pid is written.
_START_GRACE_S = 120
# gmx trjconv may be quiet between frame log lines on large files.
_LOG_ACTIVE_S = 15 * 60


class JobCancelled(Exception):
    """Raised when a Tools job is cancelled by the user."""


def _job_key(job_dir: Path | str) -> str:
    return str(Path(job_dir).expanduser().resolve())


def _cancel_flag_path(job_dir: Path) -> Path:
    return Path(job_dir) / "cancel.requested"


def _normalize_stride(stride: Optional[int]) -> int:
    try:
        n = int(stride) if stride is not None else 1
    except (TypeError, ValueError):
        n = 1
    return max(1, n)


def _resolve_file_strides(
    trajectory_files: list[str],
    file_strides: Optional[dict[str, int]] = None,
    default_stride: int = 1,
) -> dict[str, int]:
    """Map resolved trajectory paths → stride (≥1)."""
    default = _normalize_stride(default_stride)
    raw = file_strides or {}
    # Normalize lookup keys once
    keyed: dict[str, int] = {}
    for k, v in raw.items():
        try:
            keyed[str(Path(k).expanduser().resolve())] = _normalize_stride(v)
        except OSError:
            keyed[str(k)] = _normalize_stride(v)
        keyed[str(k)] = _normalize_stride(v)

    out: dict[str, int] = {}
    for traj_raw in trajectory_files:
        traj = Path(traj_raw).expanduser().resolve()
        key = str(traj)
        out[key] = (
            keyed.get(key)
            or keyed.get(str(traj_raw))
            or keyed.get(traj.name)
            or default
        )
    return out


def _is_cancel_requested(job_dir: Optional[Path]) -> bool:
    if job_dir is None:
        return False
    return _cancel_flag_path(job_dir).is_file()


def _check_cancel(job_dir: Optional[Path]) -> None:
    if _is_cancel_requested(job_dir):
        raise JobCancelled("Job cancelled by user")


def _read_pid_file(job_dir: Path) -> Optional[int]:
    pid_path = Path(job_dir) / "process.pid"
    if not pid_path.is_file():
        return None
    try:
        text = pid_path.read_text(encoding="utf-8").strip()
        return int(text) if text else None
    except (OSError, ValueError):
        return None


def _pid_alive(pid: Optional[int]) -> bool:
    if pid is None or pid <= 0:
        return False
    try:
        os.kill(pid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True  # exists but not owned by us
    except OSError:
        return False


def _parse_iso(value: Any) -> Optional[datetime]:
    raw = str(value or "").strip()
    if not raw:
        return None
    try:
        return datetime.fromisoformat(raw)
    except ValueError:
        return None


def _looks_like_false_interrupt(status: dict) -> bool:
    if str(status.get("status") or "").lower() != "error":
        return False
    return _INTERRUPT_ERROR in str(status.get("error") or "")


def _job_within_start_grace(status: dict, *, now: Optional[float] = None) -> bool:
    started = _parse_iso(status.get("start_time"))
    if started is None:
        return False
    age = (now if now is not None else time.time()) - started.timestamp()
    return 0 <= age < _START_GRACE_S


def _job_log_recently_active(job_dir: Path, *, now: Optional[float] = None) -> bool:
    """True when the worker is still appending logs (even if PID is briefly missing)."""
    latest = 0.0
    for rel in ("logs/fix_pbc.log", "logs/worker_stderr.log"):
        path = Path(job_dir) / rel
        try:
            latest = max(latest, path.stat().st_mtime)
        except OSError:
            continue
    if latest <= 0:
        return False
    age = (now if now is not None else time.time()) - latest
    return age < _LOG_ACTIVE_S


def _worker_appears_alive(job_dir: Path, status: dict, *, now: Optional[float] = None) -> bool:
    if _pid_alive(_read_pid_file(job_dir)):
        return True
    if _job_log_recently_active(job_dir, now=now):
        return True
    if _job_within_start_grace(status, now=now):
        return True
    return False


def _kill_process_group(pid: int) -> bool:
    """SIGTERM then SIGKILL the process group for a detached worker."""
    killed = False
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(os.getpgid(pid), sig)
            killed = True
        except ProcessLookupError:
            return True
        except (OSError, AttributeError):
            try:
                os.kill(pid, sig)
                killed = True
            except ProcessLookupError:
                return True
            except OSError:
                pass
        if sig == signal.SIGTERM:
            deadline = time.time() + 2.0
            while time.time() < deadline:
                if not _pid_alive(pid):
                    return True
                time.sleep(0.1)
    return killed or not _pid_alive(pid)


def _run_cancellable(
    cmd: list[str],
    *,
    stdin_text: Optional[str],
    log_fn: Callable[[str], None],
    job_dir: Optional[Path],
) -> subprocess.CompletedProcess:
    """Run a subprocess; keep it in the worker session so cancel can killpg."""
    _check_cancel(job_dir)
    # When running under the detached worker (job_dir set), children must stay in
    # the worker process group so Cancel can kill the whole tree. Outside a job,
    # isolate with a new session so a stray gmx does not take down the caller.
    proc = subprocess.Popen(
        cmd,
        stdin=subprocess.PIPE if stdin_text is not None else subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=job_dir is None,
    )
    try:
        try:
            stdout, stderr = proc.communicate(input=stdin_text)
        except Exception:
            proc.kill()
            raise
        return subprocess.CompletedProcess(cmd, proc.returncode, stdout, stderr)
    finally:
        _check_cancel(job_dir)


def cancel_fix_pbc_job(job_dir: str) -> dict:
    """
    Cancel a running Fix PBC job.

    Writes ``cancel.requested``, then kills the detached worker process group
    (same pattern as equilibration stop). Works after the GUI/backend restarts.
    """
    path = Path(job_dir).expanduser().resolve()
    if not path.is_dir():
        raise FileNotFoundError(f"Job directory not found: {path}")
    if not (path / "tools_job.json").is_file() and not (path / "status.json").is_file():
        raise FileNotFoundError(f"Not a Tools job directory: {path}")

    status_path = path / "status.json"
    status = "unknown"
    if status_path.is_file():
        try:
            status = json.loads(status_path.read_text(encoding="utf-8")).get(
                "status", "unknown"
            )
        except (json.JSONDecodeError, OSError):
            pass

    if status in {"completed", "error", "cancelled"}:
        return {
            "success": True,
            "job_dir": str(path),
            "stopped": False,
            "status": status,
            "message": f"Job already finished ({status})",
        }

    _cancel_flag_path(path).write_text("1\n", encoding="utf-8")
    pid = _read_pid_file(path)
    killed = False
    if pid is not None and _pid_alive(pid):
        killed = _kill_process_group(pid)

    _update_status(
        path,
        status="cancelled",
        error="Cancelled by user",
        end_time=_now_iso(),
    )
    _append_log(path, "CANCELLED by user")
    (path / "process.pid").unlink(missing_ok=True)

    return {
        "success": True,
        "job_dir": str(path),
        "stopped": True,
        "killed_process": killed,
        "status": "cancelled",
        "message": "Cancel requested",
    }


# ── status helpers ────────────────────────────────────────────────────────────


def _now_iso() -> str:
    return datetime.now().isoformat()


def _write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def _update_status(job_dir: Path, **fields: Any) -> None:
    status_path = job_dir / "status.json"
    data: dict[str, Any] = {}
    if status_path.is_file():
        try:
            data = json.loads(status_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            data = {}
    data.update(fields)
    _write_json(status_path, data)


def _append_log(job_dir: Path, message: str) -> None:
    log_path = job_dir / "logs" / "fix_pbc.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(message.rstrip() + "\n")


def _mark_step(job_dir: Path, step: str, *, completed: bool = False) -> None:
    status_path = job_dir / "status.json"
    data: dict[str, Any] = {}
    if status_path.is_file():
        try:
            data = json.loads(status_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            data = {}
    steps = list(data.get("steps") or [])
    done = list(data.get("steps_completed") or [])
    if step not in steps:
        steps.append(step)
    if completed and step not in done:
        done.append(step)
    try:
        current = steps.index(step) + (1 if completed else 0)
    except ValueError:
        current = len(done)
    data["steps"] = steps
    data["steps_completed"] = done
    data["current_step"] = current
    # A live worker must not stay stuck on a false "Interrupted" scan result.
    if _looks_like_false_interrupt(data) or str(data.get("status") or "").lower() == "running":
        data["status"] = "running"
        data["error"] = None
        data["end_time"] = None
    _write_json(status_path, data)
    _append_log(job_dir, f"[{'done' if completed else 'step'}] {step}")


# ── executable resolution ─────────────────────────────────────────────────────


def resolve_gmx_executable(explicit: Optional[str] = None) -> Optional[str]:
    """Resolve a ``gmx`` binary path."""
    if explicit and Path(explicit).is_file():
        return str(Path(explicit).resolve())
    try:
        from gatewizard.utils.optional_deps import list_md_engine_candidates

        for cand in list_md_engine_candidates("gromacs"):
            exe = cand.get("executable")
            if exe and Path(exe).is_file():
                return str(Path(exe).resolve())
    except Exception:
        pass
    which = shutil.which("gmx") or shutil.which("gmx_mpi")
    return str(Path(which).resolve()) if which else None


def resolve_cpptraj_executable(explicit: Optional[str] = None) -> Optional[str]:
    """Resolve a ``cpptraj`` binary (AmberTools)."""
    if explicit and Path(explicit).is_file():
        return str(Path(explicit).resolve())
    which = shutil.which("cpptraj")
    if which:
        return str(Path(which).resolve())
    amberhome = os.environ.get("AMBERHOME")
    if amberhome:
        cand = Path(amberhome) / "bin" / "cpptraj"
        if cand.is_file():
            return str(cand.resolve())
    conda = os.environ.get("CONDA_PREFIX")
    if conda:
        cand = Path(conda) / "bin" / "cpptraj"
        if cand.is_file():
            return str(cand.resolve())
    return None


# ── engine detection ──────────────────────────────────────────────────────────


def _sibling_engine_hint(path: Path) -> Optional[str]:
    """Infer engine from files near a topology/trajectory path."""
    search_dirs = [path.parent]
    # Walk up one level (engine subfolder under equilibration/)
    if path.parent.parent != path.parent:
        search_dirs.append(path.parent.parent)

    for d in search_dirs:
        meta = d / "equilibration_job.json"
        if meta.is_file():
            try:
                raw = json.loads(meta.read_text(encoding="utf-8"))
                eng = (raw.get("engine") or "").strip().lower()
                if eng in {"namd", "gromacs", "openmm", "amber"}:
                    return eng
            except Exception:
                pass
        if (d / "openmm_run.py").is_file() or list(d.glob("step*_equilibration.inp")):
            return "openmm"
        if list(d.glob("*.mdp")) or (d / "index.ndx").is_file() or list(d.glob("*.tpr")):
            return "gromacs"
        if list(d.glob("step*.mdin")) or list(d.glob("*.mdout")):
            return "amber"
        if list(d.glob("step*_equilibration.conf")) or list(
            d.glob("step*_minimization.conf")
        ):
            return "namd"
    return None


def _find_companion_tpr(
    trajectory: Optional[Path] = None,
    topology: Optional[Path] = None,
    explicit: Optional[Path] = None,
) -> Optional[Path]:
    """Locate a GROMACS ``.tpr`` near the trajectory and/or topology.

    Search order:
    1. Explicit path (user Browse for GROMACS TPR)
    2. Topology itself if it is a ``.tpr``
    3. Same stem / folder as the trajectory
    4. Same stem / folder as the topology (e.g. ``system.pdb`` → ``step*.tpr``)
    """
    if explicit is not None:
        exp = Path(explicit).expanduser()
        if exp.is_file() and exp.suffix.lower() == ".tpr":
            return exp.resolve()
    if topology and topology.suffix.lower() == ".tpr" and topology.is_file():
        return topology.resolve()

    candidates: list[Path] = []
    if trajectory is not None:
        candidates.extend(
            [
                trajectory.with_suffix(".tpr"),
                trajectory.parent / f"{trajectory.stem}.tpr",
            ]
        )
        if topology:
            candidates.append(topology.parent / f"{trajectory.stem}.tpr")
    if topology is not None:
        candidates.append(topology.with_suffix(".tpr"))
        # Prefer CHARMM-GUI / GateWizard stage TPRs next to the topology.
        parent = topology.parent
        for pattern in ("step7_production.tpr", "step6_*.tpr", "step0_*.tpr", "*.tpr"):
            for hit in sorted(parent.glob(pattern)):
                if hit.is_file():
                    candidates.append(hit)

    seen: set[Path] = set()
    for c in candidates:
        try:
            resolved = c.resolve()
        except OSError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        if resolved.is_file():
            return resolved

    if trajectory is not None:
        tprs = sorted(p for p in trajectory.parent.glob("*.tpr") if p.is_file())
        if tprs:
            return tprs[0].resolve()
    return None


def _find_companion_ndx(
    trajectory: Optional[Path] = None,
    tpr: Optional[Path] = None,
    topology: Optional[Path] = None,
    explicit: Optional[Path] = None,
) -> Optional[Path]:
    if explicit is not None:
        exp = Path(explicit).expanduser()
        if exp.is_file():
            return exp.resolve()
    search: list[Path] = []
    if trajectory is not None:
        search.append(trajectory.parent)
    if tpr is not None:
        search.append(tpr.parent)
    if topology is not None:
        search.append(topology.parent)
    seen: set[Path] = set()
    for d in search:
        try:
            resolved = d.resolve()
        except OSError:
            continue
        if resolved in seen:
            continue
        seen.add(resolved)
        ndx = resolved / "index.ndx"
        if ndx.is_file():
            return ndx.resolve()
    return None


def _prmtop_residue_names(topology: Path) -> list[str]:
    """Unique residue labels from an Amber ``prmtop``/``parm7`` (``RESIDUE_LABEL``)."""
    try:
        text = topology.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    match = re.search(
        r"%FLAG\s+RESIDUE_LABEL\s*\n%FORMAT\([^)]*\)\s*\n(.*?)(?=\n%FLAG|\Z)",
        text,
        flags=re.DOTALL | re.IGNORECASE,
    )
    if not match:
        return []
    blob = match.group(1).replace("\r", "").replace("\n", "")
    names: list[str] = []
    seen: set[str] = set()
    for i in range(0, len(blob), 4):
        name = blob[i : i + 4].strip()
        if not name or name in seen:
            continue
        seen.add(name)
        names.append(name)
    return names


def _pdb_residue_names(topology: Path) -> list[str]:
    """Unique residue names from PDB ATOM/HETATM records (columns 18–21)."""
    names: list[str] = []
    seen: set[str] = set()
    try:
        with topology.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if not line.startswith(("ATOM", "HETATM")):
                    continue
                name = line[17:21].strip().upper()
                if not name or name in seen:
                    continue
                seen.add(name)
                names.append(name)
    except OSError:
        return []
    return names


def _psf_residue_names(topology: Path) -> list[str]:
    """Unique residue names from a CHARMM/NAMD PSF atom section (4th field)."""
    names: list[str] = []
    seen: set[str] = set()
    try:
        with topology.open(encoding="utf-8", errors="replace") as fh:
            in_atoms = False
            remaining = 0
            for line in fh:
                upper = line.upper()
                if not in_atoms:
                    if "!NATOM" not in upper:
                        continue
                    in_atoms = True
                    try:
                        remaining = int(line.split("!")[0].split()[0])
                    except (IndexError, ValueError):
                        remaining = 10**9
                    continue
                s = line.strip()
                if not s:
                    continue
                if s.startswith("!") or "!N" in upper:
                    break
                parts = s.split()
                if len(parts) < 4:
                    continue
                remaining -= 1
                name = parts[3].strip().upper()
                if name and name not in seen:
                    seen.add(name)
                    names.append(name)
                if remaining <= 0:
                    break
    except OSError:
        return []
    return names


def _gro_residue_names(topology: Path) -> list[str]:
    """Unique residue names from a GROMACS GRO (columns 6–10)."""
    names: list[str] = []
    seen: set[str] = set()
    try:
        with topology.open(encoding="utf-8", errors="replace") as fh:
            for i, line in enumerate(fh):
                if i < 2:
                    continue
                if len(line) < 10:
                    continue
                name = line[5:10].strip().upper()
                if not name or name in seen:
                    continue
                seen.add(name)
                names.append(name)
    except OSError:
        return []
    return names


def _topology_residue_names(topology: Path) -> list[str]:
    """Residue names from prmtop, PSF, PDB, or GRO — used to pick lipid defaults."""
    ext = topology.suffix.lower()
    if ext in {".prmtop", ".parm7"}:
        return _prmtop_residue_names(topology)
    if ext == ".psf":
        return _psf_residue_names(topology)
    if ext in {".pdb", ".ent"}:
        return _pdb_residue_names(topology)
    if ext == ".gro":
        return _gro_residue_names(topology)
    return []


def ordered_lipid_resnames(residue_names: list[str] | set[str]) -> list[str]:
    """Lipid residue names from a topology, in a stable display order."""
    found = {str(n).strip().upper() for n in residue_names if str(n).strip()}
    lipids = found & _AMBER_LIPID_RESNAMES
    if not lipids:
        return []
    ordered = [n for n in _AMBER_LIPID_RESNAME_ORDER if n in lipids]
    extra = sorted(lipids - set(ordered))
    return ordered + extra


def default_cpptraj_membrane_center_selection() -> str:
    """Protein + common bilayer residue names (Amber lipid21 and CHARMM lipids)."""
    return "protein or resname " + " ".join(_DEFAULT_BILAYER_CENTER_RESNAMES)


def recommend_cpptraj_center_selection(
    topology_file: str | Path | None = None,
    residue_names: Optional[list[str]] = None,
) -> tuple[str, list[str]]:
    """
    Default cpptraj/MDA center selection for Amber/NAMD/OpenMM membranes.

    Protein-only autoimage splits the bilayer. When lipid residues are present
    in the topology (prmtop, PSF, PDB, GRO), return ``protein or resname …``.
    Otherwise return protein + a compact common-lipid list so PSF/PDB loads
    do not fall back to protein-only.
    """
    names = list(residue_names or [])
    if not names and topology_file:
        path = Path(topology_file).expanduser()
        if path.is_file():
            names = _topology_residue_names(path)
    lipids = ordered_lipid_resnames(names)
    if not lipids:
        return default_cpptraj_membrane_center_selection(), []
    return f"protein or resname {' '.join(lipids)}", lipids


def detect_pbc_engine(
    topology_file: str,
    trajectory_files: list[str],
    engine_hint: Optional[str] = None,
    tpr_path: Optional[str] = None,
    ndx_path: Optional[str] = None,
) -> dict:
    """
    Detect which MD engine / tool should fix PBC for the given inputs.

    Returns dict with ``engine``, ``method``, ``reason``, ``tpr``, ``ndx``,
    ``warnings``, and for Amber/NAMD/OpenMM ``lipid_resnames`` plus
    ``recommended_center_selection`` (protein + detected or common lipids).

    ``tpr_path`` / ``ndx_path`` are optional explicit overrides from the Tools UI
    (GROMACS TPR / Index browse fields).
    """
    hint = (engine_hint or "auto").strip().lower()
    if hint not in _VALID_ENGINES:
        hint = "auto"

    top = Path(topology_file).expanduser().resolve() if topology_file else None
    trajs = [Path(p).expanduser().resolve() for p in trajectory_files]
    warnings: list[str] = []

    top_ext = top.suffix.lower() if top else ""
    traj_exts = {t.suffix.lower() for t in trajs}

    sibling = None
    for p in ([top] if top else []) + trajs:
        if p:
            sibling = _sibling_engine_hint(p)
            if sibling:
                break

    engine = hint if hint != "auto" else None
    reason = ""

    if engine is None:
        if top_ext == ".tpr" or traj_exts & _GROMACS_TRAJ or sibling == "gromacs":
            engine = "gromacs"
            reason = "GROMACS inputs (.tpr/.xtc/.trr or job folder)"
        elif traj_exts & _AMBER_TRAJ or sibling == "amber":
            engine = "amber"
            reason = "Amber NetCDF trajectory or job folder"
        elif sibling == "openmm":
            engine = "openmm"
            reason = "OpenMM job folder markers"
        elif sibling == "namd" or top_ext == ".psf":
            engine = "namd"
            reason = "NAMD job folder or PSF topology"
        elif top_ext in {".prmtop", ".parm7"} and traj_exts & _DCD_TRAJ:
            engine = "namd"
            reason = "Amber topology + DCD (typical GateWizard NAMD/OpenMM)"
            warnings.append(
                "Could not distinguish NAMD vs OpenMM; using cpptraj (works for both)."
            )
        else:
            engine = "mdanalysis"
            reason = "No engine-specific inputs detected; using MDAnalysis fallback"

    if hint != "auto":
        reason = f"User selected {engine}"

    method = {
        "gromacs": "gmx trjconv",
        "amber": "cpptraj autoimage",
        "namd": "cpptraj autoimage",
        "openmm": "cpptraj autoimage",
        "mdanalysis": "MDAnalysis make_whole + wrap",
    }.get(engine, "unknown")

    tpr = _find_companion_tpr(
        trajs[0] if trajs else None,
        top,
        explicit=Path(tpr_path).expanduser() if tpr_path else None,
    )
    ndx = _find_companion_ndx(
        trajs[0] if trajs else None,
        tpr,
        topology=top,
        explicit=Path(ndx_path).expanduser() if ndx_path else None,
    )

    if engine == "gromacs" and tpr is None and trajs:
        warnings.append(
            "No .tpr found next to the trajectory or topology. "
            "Browse for the matching step*.tpr under GROMACS TPR."
        )
    if engine in {"amber", "namd", "openmm"} and top_ext not in {
        ".prmtop",
        ".parm7",
        ".top",
        ".psf",
    }:
        warnings.append(
            f"Topology {top_ext or '(missing)'} may not work with cpptraj; prefer system.prmtop."
        )

    center_groups: list[dict] = []
    recommended_center = None
    recommended_center_groups: list[str] = []
    recommended_output = "System"
    if engine == "gromacs" and (tpr or ndx):
        group_info = list_gromacs_index_groups(
            ndx_path=str(ndx) if ndx else None,
            tpr_path=str(tpr) if tpr else None,
        )
        center_groups = group_info.get("groups") or []
        recommended_center = group_info.get("recommended")
        recommended_center_groups = list(group_info.get("recommended_groups") or [])
        warnings.extend(group_info.get("warnings") or [])
        if any(g.get("name") == "System" for g in center_groups):
            recommended_output = "System"
        elif center_groups:
            recommended_output = center_groups[0]["name"]

    lipid_resnames: list[str] = []
    recommended_center_selection = "protein"
    if engine != "gromacs":
        recommended_center_selection, lipid_resnames = recommend_cpptraj_center_selection(
            top
        )

    return {
        "engine": engine,
        "method": method,
        "reason": reason,
        "tpr": str(tpr) if tpr else None,
        "ndx": str(ndx) if ndx else None,
        "topology": str(top) if top else None,
        "warnings": warnings,
        "center_groups": center_groups,
        "recommended_center": recommended_center,
        "recommended_center_groups": recommended_center_groups,
        "recommended_output": recommended_output,
        "lipid_resnames": lipid_resnames,
        "recommended_center_selection": recommended_center_selection,
        "supported_output_formats": (
            ["xtc", "trr", "same"]
            if engine == "gromacs"
            else ["dcd", "xtc", "nc", "same"]
        ),
    }


# ── GROMACS trjconv ───────────────────────────────────────────────────────────

# Common lipid / membrane index names (CHARMM-GUI / GateWizard / custom).
_LIPID_LIKE_GROUP_NAMES = frozenset(
    {
        "PA",
        "PC",
        "OL",
        "PE",
        "PS",
        "PG",
        "PI",
        "SM",
        "CHL",
        "CHOL",
        "POPC",
        "POPE",
        "POPS",
        "POPG",
        "DOPC",
        "DOPE",
        "DPPC",
        "DPPE",
        "MEMB",
        "Membrane",
        "membrane",
        "Lipid",
        "lipid",
        "LIPID",
        "Lipids",
        "lipids",
    }
)
_PROTEIN_LIKE_GROUP_NAMES = ("SOLU_MEMB", "Protein", "SOLU", "Protein-H")


def _ndx_has_group(ndx_path: Path, name: str) -> bool:
    try:
        text = ndx_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return re.search(rf"\[\s*{re.escape(name)}\s*\]", text) is not None


def _is_lipid_like_group(name: str) -> bool:
    n = (name or "").strip()
    if not n:
        return False
    if n in _LIPID_LIKE_GROUP_NAMES:
        return True
    # Heuristic: short all-caps lipid codes (2–4 letters) excluding water/ions
    if re.fullmatch(r"[A-Z]{2,4}", n) and n not in {
        "SOL",
        "WAT",
        "HOH",
        "ION",
        "NA",
        "CL",
        "POT",
        "SOD",
        "CLA",
    }:
        return True
    return False


def _normalize_group_name_list(
    names: Optional[list[str]] = None,
    singular: Optional[str] = None,
) -> list[str]:
    """Unique non-empty group names; ``names`` wins over singular when non-empty."""
    raw_list = list(names or [])
    if any(str(x or "").strip() for x in raw_list):
        seq = raw_list
    else:
        seq = [singular] if singular else []
    out: list[str] = []
    seen: set[str] = set()
    for raw in seq:
        name = str(raw or "").strip()
        if not name or name in seen:
            continue
        seen.add(name)
        out.append(name)
    return out


def list_ndx_groups(ndx_path: str | Path) -> list[dict]:
    """Parse ``[ GroupName ]`` entries from a GROMACS index file."""
    path = Path(ndx_path).expanduser().resolve()
    if not path.is_file():
        return []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []

    groups: list[dict] = []
    current: Optional[str] = None
    atoms: list[str] = []
    header_re = re.compile(r"\[\s*([^\]]+?)\s*\]")

    def _flush() -> None:
        nonlocal current, atoms
        if current is None:
            return
        n = 0
        for tok in atoms:
            n += len(tok.split())
        groups.append(
            {
                "name": current,
                "index": len(groups),
                "n_atoms": n,
                "recommended": current in {"SOLU_MEMB", "Protein", "SOLU"},
                "lipid_like": _is_lipid_like_group(current),
            }
        )
        current = None
        atoms = []

    for line in text.splitlines():
        m = header_re.match(line.strip())
        if m:
            _flush()
            current = m.group(1).strip()
            atoms = []
            continue
        if current is not None:
            atoms.append(line)
    _flush()
    return groups


def parse_ndx_group_atoms(ndx_path: str | Path, group_name: str) -> list[int]:
    """Return atom indices (1-based GROMACS) for ``group_name`` in an ndx file."""
    path = Path(ndx_path).expanduser().resolve()
    want = (group_name or "").strip()
    if not path.is_file() or not want:
        return []
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return []
    header_re = re.compile(r"\[\s*([^\]]+?)\s*\]")
    current: Optional[str] = None
    atoms: list[int] = []
    for line in text.splitlines():
        m = header_re.match(line.strip())
        if m:
            if current == want:
                break
            current = m.group(1).strip()
            atoms = []
            continue
        if current == want:
            for tok in line.split():
                try:
                    atoms.append(int(tok))
                except ValueError:
                    continue
    return atoms


def build_compound_ndx(
    source_ndx: Path,
    group_names: list[str],
    compound_name: str,
    dest_ndx: Path,
) -> dict:
    """
    Copy ``source_ndx`` and append a compound group that is the union of
    ``group_names`` atom indices.

    Returns ``{compound_name, source_groups, n_atoms, path}``.
    """
    src = Path(source_ndx).expanduser().resolve()
    if not src.is_file():
        raise FileNotFoundError(f"Index file not found: {src}")
    names = _normalize_group_name_list(group_names)
    if not names:
        raise ValueError("No index groups selected for compound ndx")
    compound = (compound_name or "GW_COMPOUND").strip() or "GW_COMPOUND"

    union: list[int] = []
    seen: set[int] = set()
    resolved: list[str] = []
    for name in names:
        atoms = parse_ndx_group_atoms(src, name)
        if not atoms:
            raise ValueError(f"Index group {name!r} not found (or empty) in {src.name}")
        resolved.append(name)
        for a in atoms:
            if a not in seen:
                seen.add(a)
                union.append(a)

    try:
        base_text = src.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        raise RuntimeError(f"Could not read index file: {exc}") from exc

    # Drop an existing compound section with the same name, then append fresh.
    header_re = re.compile(rf"^\[\s*{re.escape(compound)}\s*\]\s*$", re.MULTILINE)
    if header_re.search(base_text):
        parts = re.split(r"(?=^\[)", base_text, flags=re.MULTILINE)
        kept = []
        for part in parts:
            if not part.strip():
                continue
            if re.match(rf"\[\s*{re.escape(compound)}\s*\]", part.strip()):
                continue
            kept.append(part)
        base_text = "".join(kept)
    if base_text and not base_text.endswith("\n"):
        base_text += "\n"

    lines = [f"[ {compound} ]"]
    row: list[str] = []
    for i, atom in enumerate(union, start=1):
        row.append(f"{atom:4d}")
        if i % 15 == 0:
            lines.append(" ".join(row))
            row = []
    if row:
        lines.append(" ".join(row))
    lines.append("")

    dest = Path(dest_ndx).expanduser().resolve()
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(base_text + "\n".join(lines), encoding="utf-8")
    return {
        "compound_name": compound,
        "source_groups": resolved,
        "n_atoms": len(union),
        "path": str(dest),
    }


def _ensure_source_ndx(
    *,
    ndx: Optional[Path],
    tpr: Optional[Path],
    gmx: str,
    work_dir: Path,
) -> Path:
    """Return an ndx path suitable for merging (use existing or make_ndx from TPR)."""
    if ndx and ndx.is_file():
        return ndx
    if tpr is None or not tpr.is_file():
        raise FileNotFoundError(
            "Multi-group Fix PBC requires index.ndx or a .tpr to build one."
        )
    work_dir.mkdir(parents=True, exist_ok=True)
    out = work_dir / "gw_base_groups.ndx"
    proc = subprocess.run(
        [gmx, "make_ndx", "-f", str(tpr), "-o", str(out)],
        input="q\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=120,
    )
    if not out.is_file():
        raise RuntimeError(
            "gmx make_ndx failed to create a temporary index for multi-group merge: "
            f"{(proc.stderr or proc.stdout or '')[-400:]}"
        )
    return out


def _prepare_gromacs_fix_groups(
    *,
    ndx: Optional[Path],
    tpr: Optional[Path],
    gmx: str,
    work_dir: Path,
    center_group: Optional[str] = None,
    center_groups: Optional[list[str]] = None,
    output_group: Optional[str] = None,
    output_groups: Optional[list[str]] = None,
) -> dict:
    """
    Resolve center/output names and an ndx path for the trjconv pipeline.

    Multi-select merges selected groups into temporary ``GW_CENTER`` /
    ``GW_OUTPUT`` entries on a copied index file.
    """
    center_sel = _normalize_group_name_list(center_groups, center_group)
    output_sel = _normalize_group_name_list(output_groups, output_group)
    ndx_file = ndx if ndx and ndx.is_file() else None

    need_compound = len(center_sel) > 1 or len(output_sel) > 1
    working_ndx = ndx_file
    if need_compound:
        source = _ensure_source_ndx(ndx=ndx_file, tpr=tpr, gmx=gmx, work_dir=work_dir)
        working_ndx = work_dir / "gw_fix_pbc_groups.ndx"
        # Start from a fresh copy of the source index
        shutil.copy2(source, working_ndx)

    if len(center_sel) > 1:
        assert working_ndx is not None
        info = build_compound_ndx(working_ndx, center_sel, "GW_CENTER", working_ndx)
        cg = info["compound_name"]
        center_sources = list(info["source_groups"])
        center_label = f"{cg} = {'+'.join(center_sources)}"
    else:
        cg = (
            center_sel[0]
            if center_sel
            else _gromacs_center_group(ndx_file, None)
        )
        center_sources = [cg] if cg else []
        center_label = cg

    if len(output_sel) > 1:
        assert working_ndx is not None
        info = build_compound_ndx(working_ndx, output_sel, "GW_OUTPUT", working_ndx)
        og = info["compound_name"]
        output_sources = list(info["source_groups"])
        output_label = f"{og} = {'+'.join(output_sources)}"
    else:
        og = (
            output_sel[0]
            if output_sel
            else _gromacs_output_group(ndx_file, None)
        )
        output_sources = [og] if og else []
        output_label = og

    return {
        "center_group": cg,
        "output_group": og,
        "center_sources": center_sources,
        "output_sources": output_sources,
        "center_label": center_label,
        "output_label": output_label,
        "ndx": working_ndx,
    }


def _recommend_center_groups(groups: list[dict]) -> tuple[str, list[str]]:
    """Pick recommended center name + multi-select list for the GUI."""
    names = [str(g.get("name") or "") for g in groups if g.get("name")]
    name_set = set(names)
    for pref in _PROTEIN_LIKE_GROUP_NAMES:
        if pref == "SOLU_MEMB" and pref in name_set:
            return pref, [pref]
    lipids = [n for n in names if _is_lipid_like_group(n)]
    # Prefer named MEMB/Membrane as a single pick when no split lipids
    if "SOLU_MEMB" not in name_set and lipids:
        protein = next((p for p in ("Protein", "Protein-H", "SOLU") if p in name_set), None)
        # Split lipid codes (PA/PC/OL) → multi-select; keep order from ndx
        splitish = [n for n in lipids if n.upper() not in {"MEMB", "MEMBRANE", "LIPID", "LIPIDS"}]
        if len(splitish) >= 2:
            rec = ([protein] if protein else []) + splitish
            return rec[0], rec
        if protein and lipids:
            rec = [protein] + lipids
            return rec[0], rec
        return lipids[0], lipids
    for pref in ("Protein", "SOLU", "Protein-H"):
        if pref in name_set:
            return pref, [pref]
    if names:
        return names[0], [names[0]]
    return "Protein", ["Protein"]


def list_gromacs_index_groups(
    *,
    ndx_path: Optional[str] = None,
    tpr_path: Optional[str] = None,
    gmx_executable: Optional[str] = None,
) -> dict:
    """
    List centering/output groups for GROMACS Fix PBC.

    Prefers parsing ``index.ndx``. If missing, runs ``gmx make_ndx -f tpr``
    (quit immediately) and parses the printed group table.
    """
    warnings: list[str] = []
    groups: list[dict] = []
    source = None

    ndx = Path(ndx_path).expanduser().resolve() if ndx_path else None
    if ndx and ndx.is_file():
        groups = list_ndx_groups(ndx)
        source = str(ndx)
    elif tpr_path:
        tpr = Path(tpr_path).expanduser().resolve()
        gmx = resolve_gmx_executable(gmx_executable)
        if tpr.is_file() and gmx:
            try:
                import tempfile

                with tempfile.TemporaryDirectory() as tmp:
                    tmp_ndx = Path(tmp) / "tmp.ndx"
                    proc = subprocess.run(
                        [gmx, "make_ndx", "-f", str(tpr), "-o", str(tmp_ndx)],
                        input="q\n",
                        text=True,
                        capture_output=True,
                        check=False,
                        timeout=60,
                    )
                blob = (proc.stdout or "") + "\n" + (proc.stderr or "")
                # Lines like: "  0 System              : 12345 atoms"
                for m in re.finditer(
                    r"^\s*(\d+)\s+(.+?)\s*:\s*(\d+)\s+atoms",
                    blob,
                    re.MULTILINE,
                ):
                    name = m.group(2).strip()
                    groups.append(
                        {
                            "name": name,
                            "index": int(m.group(1)),
                            "n_atoms": int(m.group(3)),
                            "recommended": name
                            in {"Protein", "SOLU_MEMB", "SOLU", "Protein-H"},
                            "lipid_like": _is_lipid_like_group(name),
                        }
                    )
                source = f"gmx make_ndx:{tpr}"
            except Exception as exc:
                warnings.append(f"Could not list groups from tpr: {exc}")
        elif not gmx:
            warnings.append("gmx not found; cannot list default TPR groups.")
        else:
            warnings.append(f"TPR not found: {tpr_path}")
    else:
        warnings.append("Provide index.ndx or a .tpr to list centering groups.")

    # Deduplicate by name keeping first
    seen: set[str] = set()
    unique: list[dict] = []
    for g in groups:
        if g["name"] in seen:
            continue
        seen.add(g["name"])
        unique.append(g)
    groups = unique

    recommended, recommended_groups = _recommend_center_groups(groups)
    for g in groups:
        g["recommended"] = g["name"] in set(recommended_groups) or g["name"] == recommended

    return {
        "groups": groups,
        "recommended": recommended,
        "recommended_groups": recommended_groups,
        "source": source,
        "warnings": warnings,
    }


def _gromacs_center_group(ndx: Optional[Path], explicit: Optional[str] = None) -> str:
    if explicit and str(explicit).strip():
        return str(explicit).strip()
    if ndx and _ndx_has_group(ndx, "SOLU_MEMB"):
        return "SOLU_MEMB"
    if ndx and _ndx_has_group(ndx, "Protein"):
        return "Protein"
    return "Protein"


def _gromacs_output_group(ndx: Optional[Path], explicit: Optional[str] = None) -> str:
    if explicit and str(explicit).strip():
        return str(explicit).strip()
    if ndx and _ndx_has_group(ndx, "System"):
        return "System"
    return "System"


def _normalize_gromacs_output_format(output_format: str) -> tuple[str, Optional[str]]:
    """
    GROMACS trjconv cannot write DCD/NetCDF.

    Returns (normalized_format, warning_or_None).
    """
    fmt = (output_format or "same").strip().lower()
    if fmt in {"xtc", "trr", "same"}:
        return fmt, None
    if fmt in {"dcd", "nc", "ncdf", "mdcrd"}:
        return "xtc", (
            f"GROMACS trjconv cannot write {fmt.upper()}; using XTC instead. "
            "Select XTC (or Same) for GROMACS jobs."
        )
    return "xtc", f"Unsupported GROMACS output format {fmt!r}; using XTC."


def _run_trjconv(
    gmx: str,
    *,
    tpr: Path,
    traj_in: Path,
    traj_out: Path,
    ndx: Optional[Path],
    pbc: str,
    center: bool,
    center_group: str,
    output_group: str,
    ur: Optional[str],
    log_fn: Callable[[str], None],
    job_dir: Optional[Path] = None,
    skip: int = 1,
) -> None:
    cmd = [gmx, "trjconv", "-s", str(tpr), "-f", str(traj_in), "-o", str(traj_out), "-pbc", pbc]
    if ndx:
        cmd.extend(["-n", str(ndx)])
    if center:
        cmd.append("-center")
    if ur:
        cmd.extend(["-ur", ur])
    skip_n = _normalize_stride(skip)
    if skip_n > 1:
        cmd.extend(["-skip", str(skip_n)])

    if center:
        stdin = f"{center_group}\n{output_group}\n"
    elif pbc == "cluster":
        stdin = f"{center_group}\n{output_group}\n"
    else:
        stdin = f"{output_group}\n"

    log_fn(f"$ {' '.join(cmd)}")
    log_fn(f"  stdin groups: {stdin.strip().replace(chr(10), ' | ')}")
    proc = _run_cancellable(cmd, stdin_text=stdin, log_fn=log_fn, job_dir=job_dir)
    if proc.stdout:
        log_fn(proc.stdout[-2000:])
    if proc.stderr:
        log_fn(proc.stderr[-2000:])
    if proc.returncode != 0:
        _check_cancel(job_dir)
        raise RuntimeError(
            f"gmx trjconv -pbc {pbc} failed (exit {proc.returncode}): "
            f"{(proc.stderr or proc.stdout or '')[-500:]}"
        )
    if not traj_out.is_file():
        # Some gmx builds ignore unsupported extensions (e.g. .dcd) and write .xtc.
        alts = [
            traj_out.with_suffix(ext)
            for ext in (".xtc", ".trr")
            if ext != traj_out.suffix.lower()
        ]
        found = next((p for p in alts if p.is_file()), None)
        if found is not None:
            raise RuntimeError(
                f"gmx trjconv wrote {found.name} instead of {traj_out.name}. "
                "GROMACS cannot write DCD — use output format XTC or TRR."
            )
        raise RuntimeError(
            f"gmx trjconv did not write {traj_out}. "
            "GROMACS cannot write DCD — use output format XTC or TRR."
        )


def _fix_pbc_gromacs_one(
    *,
    tpr: Path,
    trajectory: Path,
    output_path: Path,
    ndx: Optional[Path],
    gmx: str,
    center_group: Optional[str],
    output_group: Optional[str] = None,
    center_groups: Optional[list[str]] = None,
    output_groups: Optional[list[str]] = None,
    skip_cluster: bool = False,
    stride: int = 1,
    log_fn: Callable[[str], None],
    job_dir: Optional[Path] = None,
) -> dict:
    """Recommended membrane-protein sequence: whole → nojump → cluster → mol+center."""
    work = output_path.parent / f".tmp_{trajectory.stem}_pbc"
    work.mkdir(parents=True, exist_ok=True)
    try:
        resolved = _prepare_gromacs_fix_groups(
            ndx=ndx,
            tpr=tpr,
            gmx=gmx,
            work_dir=work,
            center_group=center_group,
            center_groups=center_groups,
            output_group=output_group,
            output_groups=output_groups,
        )
        cg = resolved["center_group"]
        out_group = resolved["output_group"]
        use_ndx = resolved["ndx"]
        log_fn(f"Center group: {resolved['center_label']}")
        log_fn(f"Output group: {resolved['output_label']}")
        if len(resolved["center_sources"]) > 1:
            log_fn(
                "Merged center sources: "
                + "+".join(resolved["center_sources"])
            )
        if len(resolved["output_sources"]) > 1:
            log_fn(
                "Merged output sources: "
                + "+".join(resolved["output_sources"])
            )

        skip = _normalize_stride(stride)
        whole = work / "01_whole.xtc"
        nojump = work / "02_nojump.xtc"
        cluster = work / "03_cluster.xtc"

        # Apply stride only on the first read of the original trajectory.
        _run_trjconv(
            gmx,
            tpr=tpr,
            traj_in=trajectory,
            traj_out=whole,
            ndx=use_ndx,
            pbc="whole",
            center=False,
            center_group=cg,
            output_group=out_group,
            ur=None,
            log_fn=log_fn,
            job_dir=job_dir,
            skip=skip,
        )
        _run_trjconv(
            gmx,
            tpr=tpr,
            traj_in=whole,
            traj_out=nojump,
            ndx=use_ndx,
            pbc="nojump",
            center=False,
            center_group=cg,
            output_group=out_group,
            ur=None,
            log_fn=log_fn,
            job_dir=job_dir,
        )
        # cluster can fail for some systems; fall back to mol+center on nojump
        clustered_ok = False
        if skip_cluster:
            log_fn("cluster step skipped (user option)")
        else:
            try:
                _run_trjconv(
                    gmx,
                    tpr=tpr,
                    traj_in=nojump,
                    traj_out=cluster,
                    ndx=use_ndx,
                    pbc="cluster",
                    center=False,
                    center_group=cg,
                    output_group=out_group,
                    ur=None,
                    log_fn=log_fn,
                    job_dir=job_dir,
                )
                clustered_ok = True
            except JobCancelled:
                raise
            except RuntimeError as exc:
                log_fn(f"cluster step skipped: {exc}")

        src = cluster if clustered_ok else nojump
        _run_trjconv(
            gmx,
            tpr=tpr,
            traj_in=src,
            traj_out=output_path,
            ndx=use_ndx,
            pbc="mol",
            center=True,
            center_group=cg,
            output_group=out_group,
            ur="compact",
            log_fn=log_fn,
            job_dir=job_dir,
        )

        n_frames = _count_frames_mda(str(tpr), str(output_path))
        return {
            "input": str(trajectory),
            "output": str(output_path),
            "n_frames": n_frames,
            "format": output_path.suffix.lstrip(".").lower(),
            "ok": True,
            "error": None,
            "method": "gmx trjconv",
            "center_group": cg,
            "output_group": out_group,
            "center_groups": resolved["center_sources"],
            "output_groups": resolved["output_sources"],
            "center_label": resolved["center_label"],
            "output_label": resolved["output_label"],
            "skip_cluster": bool(skip_cluster),
            "stride": skip,
        }
    finally:
        # Best-effort cleanup of intermediates
        try:
            for p in work.glob("*"):
                p.unlink(missing_ok=True)
            work.rmdir()
        except OSError:
            pass


# ── cpptraj ───────────────────────────────────────────────────────────────────


_PROTEIN_CPPTRAJ_MASK = (
    ":ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,HID,HIE,HIP,ILE,LEU,LYS,"
    "MET,PHE,PRO,SER,THR,TRP,TYR,VAL,ACE,NME"
)

_LIPID_HINT_RE = re.compile(
    r"(?i)\b(membrane|lipid|lipids|popc|pope|pops|popg|dppc|dppe|chol|cholesterol|"
    r"dlpc|dlpe|dopc|dope|dspc|palmitoyl|oleoyl)\b|@?p31\b|"
    r":(?:[A-Za-z0-9_]*,)*(?:PA|PC|OL|PHOS|OLEO|PALM)(?:,[A-Za-z0-9_]*)*"
)


def _mask_looks_like_membrane(text: str) -> bool:
    """Heuristic: selection targets bilayer lipids / phosphates."""
    return bool(_LIPID_HINT_RE.search(text or ""))


def _split_cpptraj_or_masks(mask: str) -> list[str]:
    """Split ``(A)|(B)`` / ``A|B`` style cpptraj OR masks into parts."""
    raw = (mask or "").strip()
    if not raw:
        return []
    if re.search(r"\)\s*\|\s*\(", raw):
        # "(A)|(B)" — do not strip outer parens; each side already has them.
        return [
            p.strip().strip("()")
            for p in re.split(r"\)\s*\|\s*\(", raw)
            if p.strip()
        ]
    if "|" in raw and not raw.lstrip().startswith(":"):
        return [p.strip().strip("()") for p in raw.split("|") if p.strip()]
    return [raw]


def _is_protein_cpptraj_mask(mask: str) -> bool:
    u = (mask or "").upper()
    return "ALA" in u and ("ARG" in u or "LYS" in u or "GLY" in u)


def compact_autoimage_anchor(amber_mask: str) -> Optional[str]:
    """
    Choose a cpptraj ``autoimage`` anchor that is preferably one molecule.

    Multi-molecule anchors (all lipids, phosphate atoms, protein|membrane) make
    ``autoimage`` fail and fall back to the first molecule — often wrong for
    bilayers. Prefer the protein half of an OR mask; omit lipid-only multi-res
    anchors so the caller can use ``mode byvec moveanchor`` instead.
    """
    mask = (amber_mask or "").strip()
    if not mask:
        return None
    parts = _split_cpptraj_or_masks(mask)
    if len(parts) > 1:
        for part in parts:
            if _is_protein_cpptraj_mask(part):
                return part
        # Ambiguous multi-region — omit custom anchor
        return None
    # Atom selections spanning many lipids (e.g. @P31)
    if mask.startswith("@"):
        return None
    # Many short lipid residue names → many molecules
    if re.fullmatch(r":[A-Za-z0-9_,]+", mask):
        names = [n for n in mask[1:].split(",") if n]
        if names and not _is_protein_cpptraj_mask(mask):
            if len(names) >= 1 and all(1 <= len(n) <= 4 for n in names):
                return None
    return mask


def build_cpptraj_autoimage_attempts(
    amber_mask: str, *, membrane_like: Optional[bool] = None
) -> list[str]:
    """
    Ordered cpptraj action blocks (no parm/trajin/trajout/run).

    Never uses ``origin`` (that parks the system on the box corner). For
    membranes prefer ``mode byvec moveanchor``, then ``center`` + ``image``.
    """
    mask = (amber_mask or "").strip() or _PROTEIN_CPPTRAJ_MASK
    is_memb = (
        bool(membrane_like)
        if membrane_like is not None
        else _mask_looks_like_membrane(mask)
    )
    anchor = compact_autoimage_anchor(mask)
    attempts: list[str] = []

    def _add(block: str) -> None:
        text = "\n".join(line for line in block.strip().splitlines() if line.strip())
        if text and text not in attempts:
            attempts.append(text)

    if is_memb:
        if anchor:
            _add(
                f"autoimage anchor {anchor} mode byvec moveanchor\n"
                f"center {mask} mass\n"
                f"image center"
            )
            _add(f"autoimage anchor {anchor} mode byvec moveanchor")
        _add(
            f"autoimage mode byvec moveanchor\n"
            f"center {mask} mass\n"
            f"image center"
        )
        _add("autoimage mode byvec moveanchor")
        _add("autoimage")
    else:
        if anchor:
            _add(f"autoimage anchor {anchor}")
        _add(f"autoimage anchor {mask}")
        _add("autoimage")
    return attempts


def mda_selection_to_cpptraj_mask(selection: str) -> str:
    """
    Convert a simple MDAnalysis-style selection to an Amber/cpptraj atom mask.

    cpptraj does **not** accept space-separated ``resname PA PC OL`` (MDA syntax).
    That becomes ``anchor resname`` with leftover tokens and fails with
    ``Not all arguments handled: [ PA PC OL ]``. Use ``:PA,PC,OL`` instead.

    Supported conversions (best-effort):
      - ``protein`` → standard amino-acid residue mask
      - ``resname PA PC OL`` / ``resname PA,PC,OL`` → ``:PA,PC,OL``
      - ``name P31`` / ``name P31 P32`` → ``@P31`` / ``@P31,P32``
      - ``resid 1-100`` → ``:1-100``
      - ``A or B`` → ``(maskA)|(maskB)``
      - already Amber-like (``:…``, ``@…``) → unchanged
    """
    raw = (selection or "").strip()
    if not raw:
        return _PROTEIN_CPPTRAJ_MASK

    # Already Amber/cpptraj mask syntax
    if raw[0] in {":", "@", "*", "(", "/"}:
        return raw

    low = raw.lower()
    if low in {"protein", "protein and backbone", "backbone"}:
        return _PROTEIN_CPPTRAJ_MASK

    # Split mild boolean OR (MDA "or" / Amber "|")
    if re.search(r"\s+or\s+", raw, flags=re.IGNORECASE) or "||" in raw:
        parts = re.split(r"\s+or\s+|\|\|", raw, flags=re.IGNORECASE)
        masks = [mda_selection_to_cpptraj_mask(p.strip()) for p in parts if p.strip()]
        if len(masks) == 1:
            return masks[0]
        return "(" + ")|(".join(masks) + ")"

    m = re.match(r"^resname\s+(.+)$", raw, flags=re.IGNORECASE)
    if m:
        names = [n for n in re.split(r"[\s,]+", m.group(1).strip()) if n]
        if not names:
            raise ValueError(f"Empty resname list in selection: {selection!r}")
        return ":" + ",".join(names)

    m = re.match(r"^name\s+(.+)$", raw, flags=re.IGNORECASE)
    if m:
        names = [n for n in re.split(r"[\s,]+", m.group(1).strip()) if n]
        if not names:
            raise ValueError(f"Empty atom name list in selection: {selection!r}")
        return "@" + ",".join(names)

    m = re.match(r"^resid(?:ue)?\s+(.+)$", raw, flags=re.IGNORECASE)
    if m:
        return ":" + re.sub(r"\s+", "", m.group(1).strip())

    # Bare comma/space-separated residue names (common user shorthand)
    if re.fullmatch(r"[A-Za-z0-9_,\s]+", raw) and not re.search(
        r"\b(and|not|around|byres|same)\b", low
    ):
        names = [n for n in re.split(r"[\s,]+", raw) if n]
        if names and all(re.fullmatch(r"[A-Za-z0-9_]+", n) for n in names):
            return ":" + ",".join(names)

    # Last resort: pass through (may still fail in cpptraj)
    return raw


def _fix_pbc_cpptraj_one(
    *,
    topology: Path,
    trajectory: Path,
    output_path: Path,
    cpptraj: str,
    center_mask: str,
    stride: int = 1,
    log_fn: Callable[[str], None],
    job_dir: Optional[Path] = None,
) -> dict:
    """AmberTools cpptraj autoimage (preferred for Amber/NAMD/OpenMM Amber-FF)."""
    skip = _normalize_stride(stride)
    # trajin <file> [<start> [<stop> [<offset>]]]
    trajin_line = (
        f"trajin {trajectory}" if skip <= 1 else f"trajin {trajectory} 1 last {skip}"
    )
    amber_mask = mda_selection_to_cpptraj_mask(center_mask)
    if amber_mask != (center_mask or "").strip():
        log_fn(
            f"Converted center selection {center_mask!r} → cpptraj mask {amber_mask!r}"
        )

    membrane_like = _mask_looks_like_membrane(center_mask) or _mask_looks_like_membrane(
        amber_mask
    )
    anchor = compact_autoimage_anchor(amber_mask)
    if membrane_like:
        log_fn(
            "Membrane-like center selection: using autoimage mode byvec moveanchor "
            "(no 'origin' — that centers on the box corner)."
        )
        if anchor:
            log_fn(f"Compact autoimage anchor: {anchor}")
        else:
            log_fn(
                "No single-molecule anchor derived from the selection "
                "(lipid masks span many molecules); relying on byvec/moveanchor."
            )

    attempts = build_cpptraj_autoimage_attempts(
        amber_mask, membrane_like=membrane_like
    )
    last_err = ""
    for i, actions in enumerate(attempts):
        if output_path.is_file():
            try:
                output_path.unlink()
            except OSError:
                pass
        script = f"""parm {topology}
{trajin_line}
{actions}
trajout {output_path}
run
"""
        label = f"attempt {i + 1}/{len(attempts)}"
        log_fn(f"$ cpptraj ({label}) <<EOF\n{script}EOF")
        proc = _run_cancellable(
            [cpptraj], stdin_text=script, log_fn=log_fn, job_dir=job_dir
        )
        if proc.stdout:
            log_fn(proc.stdout[-2000:])
        if proc.stderr:
            log_fn(proc.stderr[-2000:])
        if proc.returncode == 0 and output_path.is_file():
            n_frames = _count_frames_mda(str(topology), str(output_path))
            return {
                "input": str(trajectory),
                "output": str(output_path),
                "n_frames": n_frames,
                "format": output_path.suffix.lstrip(".").lower(),
                "ok": True,
                "error": None,
                "method": "cpptraj autoimage",
                "center_mask": amber_mask,
                "center_selection_input": center_mask,
                "autoimage_actions": actions,
                "stride": skip,
            }
        _check_cancel(job_dir)
        last_err = (proc.stderr or proc.stdout or f"exit {proc.returncode}")[-500:]
        log_fn(f"WARNING: {label} failed — trying next autoimage strategy.")

    raise RuntimeError(
        f"cpptraj autoimage failed after {len(attempts)} strategies: {last_err}"
    )


# ── MDAnalysis fallback (molecule-aware) ──────────────────────────────────────


def _count_frames_mda(topology: str, trajectory: str) -> int:
    try:
        import MDAnalysis as mda

        u = mda.Universe(topology, trajectory)
        return len(u.trajectory)
    except Exception:
        return 0


def _fix_pbc_mda_one(
    *,
    topology: Path,
    trajectory: Path,
    output_path: Path,
    center_selection: str,
    stride: int = 1,
    log_fn: Callable[[str], None],
    job_dir: Optional[Path] = None,
) -> dict:
    """
    Improved MDAnalysis fallback: make molecules whole, center selection, wrap molecules.

    Still inferior to gmx/cpptraj for membrane systems, but better than residue-only wrap.
    """
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.lib.mdamath import make_whole

    _check_cancel(job_dir)
    skip = _normalize_stride(stride)
    log_fn(f"MDAnalysis fallback: {trajectory.name} → {output_path.name} (stride={skip})")
    try:
        u = mda.Universe(str(topology), str(trajectory), to_guess=("bonds",))
    except Exception:
        try:
            u = mda.Universe(str(topology), str(trajectory), guess_bonds=True)
        except Exception:
            u = mda.Universe(str(topology), str(trajectory))

    sel = (center_selection or "protein").strip() or "protein"
    center_ag = u.select_atoms(sel)
    if len(center_ag) == 0:
        center_ag = u.select_atoms("not resname WAT HOH TIP3 TIP4 SOL H2O NA CL K")
    if len(center_ag) == 0:
        raise ValueError(f"Center selection {sel!r} matched no atoms")

    n_frames = 0
    with mda.Writer(str(output_path), n_atoms=u.atoms.n_atoms) as writer:
        for frame_i, _ts in enumerate(u.trajectory):
            if skip > 1 and (frame_i % skip) != 0:
                continue
            if frame_i % 25 == 0:
                _check_cancel(job_dir)
            box = u.dimensions
            if box is None or np.any(np.asarray(box)[:3] <= 0):
                raise ValueError("Missing or invalid box dimensions")

            # Make each fragment/residue contiguous across PBC
            try:
                for frag in u.atoms.fragments:
                    if len(frag) > 1:
                        make_whole(frag)
            except Exception:
                for residue in u.residues:
                    if len(residue.atoms) > 1:
                        try:
                            make_whole(residue.atoms)
                        except Exception:
                            continue

            try:
                com = center_ag.center_of_mass()
            except Exception:
                com = center_ag.positions.mean(axis=0)
            u.atoms.positions = u.atoms.positions + (box[:3] / 2.0 - com)
            wrapped = False
            for compound in ("fragments", "residues", "atoms"):
                try:
                    u.atoms.wrap(compound=compound, center="geometry")
                    wrapped = True
                    break
                except Exception:
                    try:
                        u.atoms.wrap(compound=compound)
                        wrapped = True
                        break
                    except Exception:
                        continue
            if not wrapped:
                # Manual residue-rigid wrap into [0, L)
                L = np.asarray(box[:3], dtype=float)
                pos = u.atoms.positions.copy()
                for residue in u.residues:
                    idx = residue.atoms.indices
                    res_com = pos[idx].mean(axis=0)
                    shift = -np.floor(res_com / L) * L
                    pos[idx] += shift
                u.atoms.positions = pos
            writer.write(u.atoms)
            n_frames += 1

    return {
        "input": str(trajectory),
        "output": str(output_path),
        "n_frames": n_frames,
        "format": output_path.suffix.lstrip(".").lower(),
        "ok": True,
        "error": None,
        "method": "MDAnalysis make_whole + wrap",
        "stride": skip,
    }


# ── output naming ─────────────────────────────────────────────────────────────


def _output_path_for(
    traj: Path,
    out_dir: Path,
    *,
    engine: str,
    output_format: str,
) -> Path:
    fmt = (output_format or "same").strip().lower()
    if engine == "gromacs":
        fmt, _warn = _normalize_gromacs_output_format(fmt)
        if fmt == "same":
            ext = traj.suffix.lstrip(".").lower()
            if ext not in {"xtc", "trr"}:
                ext = "xtc"
        else:
            ext = fmt  # xtc | trr
        return out_dir / f"{traj.stem}_pbcfix.{ext}"

    if fmt == "same":
        ext = traj.suffix.lstrip(".").lower() or "dcd"
    elif fmt in {"dcd", "xtc", "nc", "trr"}:
        ext = fmt
    else:
        ext = "dcd"
    if engine == "amber" and fmt == "same":
        ext = "nc"
    return out_dir / f"{traj.stem}_pbcfix.{ext}"


# ── public sync API ───────────────────────────────────────────────────────────


def fix_pbc_trajectories(
    topology_file: str,
    trajectory_files: list[str],
    output_dir: str,
    *,
    engine: str = "auto",
    center_selection: str = "protein",
    center_group: Optional[str] = None,
    output_group: Optional[str] = None,
    center_groups: Optional[list[str]] = None,
    output_groups: Optional[list[str]] = None,
    skip_cluster: bool = False,
    tpr_path: Optional[str] = None,
    ndx_path: Optional[str] = None,
    gmx_executable: Optional[str] = None,
    cpptraj_executable: Optional[str] = None,
    output_format: str = "same",
    stride: int = 1,
    file_strides: Optional[dict[str, int]] = None,
    # legacy kwargs (ignored for engine-native paths; used by MDA fallback)
    wrap_mode: str = "center+wrap",
    residue_as_rigid: bool = True,
    log_fn: Optional[Callable[[str], None]] = None,
    progress_fn: Optional[Callable[[str, bool], None]] = None,
    job_dir: Optional[str] = None,
) -> dict:
    """
    Fix PBC using the engine-recommended tool.

    ``engine``: auto | gromacs | amber | namd | openmm | mdanalysis
    ``file_strides``: optional per-trajectory stride map (path → N); falls back to ``stride``.
    """
    del wrap_mode, residue_as_rigid  # legacy; engine-native paths supersede them
    job_path = Path(job_dir).expanduser().resolve() if job_dir else None
    stride_map = _resolve_file_strides(trajectory_files, file_strides, stride)

    def _log(msg: str) -> None:
        logger.info(msg)
        if log_fn:
            log_fn(msg)

    def _progress(step: str, done: bool = False) -> None:
        if progress_fn:
            progress_fn(step, done)

    if not trajectory_files:
        raise ValueError("No trajectory files provided")

    _check_cancel(job_path)
    detected = detect_pbc_engine(
        topology_file,
        trajectory_files,
        engine_hint=engine,
        tpr_path=tpr_path,
        ndx_path=ndx_path,
    )
    eng = detected["engine"]
    _log(f"Engine: {eng} ({detected['method']}) — {detected['reason']}")
    sel = (center_selection or "protein").strip() or "protein"
    rec_sel = str(detected.get("recommended_center_selection") or "").strip()
    if (
        eng != "gromacs"
        and rec_sel
        and rec_sel != "protein"
        and sel.lower() in {"protein", ""}
    ):
        _log(f"Default center selection (protein + detected lipids): {rec_sel}")
        center_selection = rec_sel
    if any(n > 1 for n in stride_map.values()):
        for p, n in stride_map.items():
            if n > 1:
                _log(f"Stride {Path(p).name}: every {n} frame(s)")
    for w in detected.get("warnings") or []:
        _log(f"Warning: {w}")

    top = Path(topology_file).expanduser().resolve()
    if not top.is_file():
        raise FileNotFoundError(f"Topology file not found: {top}")

    out_dir = Path(output_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    tpr = Path(tpr_path).expanduser().resolve() if tpr_path else None
    if tpr is None and detected.get("tpr"):
        tpr = Path(detected["tpr"])
    ndx = Path(ndx_path).expanduser().resolve() if ndx_path else None
    if ndx is None and detected.get("ndx"):
        ndx = Path(detected["ndx"])

    effective_format = output_format
    if eng == "gromacs":
        effective_format, fmt_warn = _normalize_gromacs_output_format(output_format)
        if fmt_warn:
            _log(f"Warning: {fmt_warn}")

    outputs: list[dict] = []

    if eng == "gromacs":
        gmx = resolve_gmx_executable(gmx_executable)
        if not gmx:
            raise RuntimeError(
                "gmx not found. Install GROMACS or set the executable path."
            )
        _log(f"Using gmx: {gmx}")
        center_sel = _normalize_group_name_list(center_groups, center_group)
        output_sel = _normalize_group_name_list(output_groups, output_group)
        ndx_for_groups = ndx if ndx and ndx.is_file() else None
        if center_sel:
            _log(
                "Center group(s): "
                + ("+".join(center_sel) if len(center_sel) > 1 else center_sel[0])
            )
        else:
            _log(f"Center group: {_gromacs_center_group(ndx_for_groups, None)}")
        if output_sel:
            _log(
                "Output group(s): "
                + ("+".join(output_sel) if len(output_sel) > 1 else output_sel[0])
            )
        else:
            _log(f"Output group: {_gromacs_output_group(ndx_for_groups, None)}")
        if skip_cluster:
            _log("trjconv cluster step: disabled")
        if tpr is None or not tpr.is_file():
            raise FileNotFoundError(
                "GROMACS Fix PBC requires a .tpr (structure/topology). "
                "Select the matching step*.tpr next to your .xtc."
            )

        for traj_raw in trajectory_files:
            traj = Path(traj_raw).expanduser().resolve()
            traj_stride = stride_map.get(str(traj), 1)
            step = f"Fix {traj.name}"
            _progress(step, False)
            entry: dict[str, Any] = {
                "input": str(traj),
                "output": None,
                "n_frames": 0,
                "format": None,
                "ok": False,
                "error": None,
                "method": "gmx trjconv",
            }
            try:
                if not traj.is_file():
                    raise FileNotFoundError(f"Trajectory not found: {traj}")
                out_path = _output_path_for(
                    traj, out_dir, engine=eng, output_format=effective_format
                )
                entry["format"] = out_path.suffix.lstrip(".").lower()
                entry["output"] = str(out_path)
                result = _fix_pbc_gromacs_one(
                    tpr=tpr,
                    trajectory=traj,
                    output_path=out_path,
                    ndx=ndx if ndx and ndx.is_file() else None,
                    gmx=gmx,
                    center_group=center_group,
                    output_group=output_group,
                    center_groups=center_groups,
                    output_groups=output_groups,
                    skip_cluster=skip_cluster,
                    stride=traj_stride,
                    log_fn=_log,
                    job_dir=job_path,
                )
                entry.update(result)
            except JobCancelled:
                raise
            except Exception as exc:
                entry["error"] = str(exc)
                _log(f"ERROR {traj.name}: {exc}")
            outputs.append(entry)
            _progress(step, True)

    elif eng in {"amber", "namd", "openmm"}:
        cpptraj = resolve_cpptraj_executable(cpptraj_executable)
        if not cpptraj:
            _log("cpptraj not found; falling back to MDAnalysis")
            eng = "mdanalysis"
        else:
            _log(f"Using cpptraj: {cpptraj}")
            # Convert MDA-style selections (e.g. "resname PA PC OL") to Amber masks
            # (":PA,PC,OL"). Passing MDA syntax makes cpptraj fail the anchor and
            # fall back to plain autoimage on the first molecule.
            try:
                mask = mda_selection_to_cpptraj_mask(center_selection)
            except ValueError as exc:
                raise RuntimeError(str(exc)) from exc
            _log(f"cpptraj autoimage anchor mask: {mask} (from {center_selection!r})")

            for traj_raw in trajectory_files:
                _check_cancel(job_path)
                traj = Path(traj_raw).expanduser().resolve()
                traj_stride = stride_map.get(str(traj), 1)
                step = f"Fix {traj.name}"
                _progress(step, False)
                entry = {
                    "input": str(traj),
                    "output": None,
                    "n_frames": 0,
                    "format": None,
                    "ok": False,
                    "error": None,
                    "method": "cpptraj autoimage",
                }
                try:
                    if not traj.is_file():
                        raise FileNotFoundError(f"Trajectory not found: {traj}")
                    out_path = _output_path_for(
                        traj, out_dir, engine=eng, output_format=output_format
                    )
                    entry["format"] = out_path.suffix.lstrip(".").lower()
                    entry["output"] = str(out_path)
                    result = _fix_pbc_cpptraj_one(
                        topology=top,
                        trajectory=traj,
                        output_path=out_path,
                        cpptraj=cpptraj,
                        center_mask=mask,
                        stride=traj_stride,
                        log_fn=_log,
                        job_dir=job_path,
                    )
                    entry.update(result)
                except JobCancelled:
                    raise
                except Exception as exc:
                    entry["error"] = str(exc)
                    _log(f"ERROR {traj.name}: {exc}")
                outputs.append(entry)
                _progress(step, True)

    if eng == "mdanalysis":
        for traj_raw in trajectory_files:
            _check_cancel(job_path)
            traj = Path(traj_raw).expanduser().resolve()
            traj_stride = stride_map.get(str(traj), 1)
            step = f"Fix {traj.name}"
            _progress(step, False)
            entry = {
                "input": str(traj),
                "output": None,
                "n_frames": 0,
                "format": None,
                "ok": False,
                "error": None,
                "method": "MDAnalysis",
            }
            try:
                if not traj.is_file():
                    raise FileNotFoundError(f"Trajectory not found: {traj}")
                out_path = _output_path_for(
                    traj, out_dir, engine="mdanalysis", output_format=output_format
                )
                entry["format"] = out_path.suffix.lstrip(".").lower()
                entry["output"] = str(out_path)
                result = _fix_pbc_mda_one(
                    topology=top,
                    trajectory=traj,
                    output_path=out_path,
                    center_selection=center_selection,
                    stride=traj_stride,
                    log_fn=_log,
                    job_dir=job_path,
                )
                entry.update(result)
            except JobCancelled:
                raise
            except Exception as exc:
                entry["error"] = str(exc)
                _log(f"ERROR {traj.name}: {exc}")
            outputs.append(entry)
            _progress(step, True)

    return {
        "outputs": outputs,
        "engine": eng if eng != "mdanalysis" or detected["engine"] == "mdanalysis" else "mdanalysis",
        "method": detected["method"] if eng != "mdanalysis" else "MDAnalysis make_whole + wrap",
        "center_selection": center_selection,
        "center_group": center_group,
        "output_group": output_group,
        "center_groups": _normalize_group_name_list(center_groups, center_group),
        "output_groups": _normalize_group_name_list(output_groups, output_group),
        "skip_cluster": bool(skip_cluster),
        "stride": _normalize_stride(stride),
        "file_strides": stride_map,
        "tpr": str(tpr) if tpr else None,
        "ndx": str(ndx) if ndx else None,
        "output_dir": str(out_dir),
        "warnings": detected.get("warnings") or [],
        "detect": detected,
    }


def run_fix_pbc(*args: Any, **kwargs: Any) -> dict:
    """Alias for :func:`fix_pbc_trajectories`."""
    return fix_pbc_trajectories(*args, **kwargs)


# ── async job API (GUI) ───────────────────────────────────────────────────────


def execute_fix_pbc_job(job_dir: Path | str) -> int:
    """
    Run a Fix PBC job from its directory (called by the detached worker).

    Returns a process exit code (0 success, 1 failure, 130 cancelled).
    """
    path = Path(job_dir).expanduser().resolve()
    meta_path = path / "tools_job.json"
    if not meta_path.is_file():
        raise FileNotFoundError(f"Missing tools_job.json in {path}")

    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        raise RuntimeError(f"Cannot read tools_job.json: {exc}") from exc

    try:
        (path / "process.pid").write_text(str(os.getpid()), encoding="utf-8")
    except OSError:
        pass

    exit_code = 0
    try:
        _check_cancel(path)
        _mark_step(path, "Detect engine", completed=True)

        def log_fn(msg: str) -> None:
            _append_log(path, msg)

        def progress_fn(step: str, done: bool) -> None:
            _mark_step(path, step, completed=done)

        topology = meta.get("topology")
        trajectories = list(meta.get("trajectories") or [])
        if not topology or not trajectories:
            raise RuntimeError("tools_job.json missing topology or trajectories")

        result = fix_pbc_trajectories(
            topology,
            trajectories,
            str(path / "outputs"),
            engine=meta.get("engine_hint") or meta.get("engine") or "auto",
            center_selection=meta.get("center_selection") or "protein",
            center_group=meta.get("center_group"),
            output_group=meta.get("output_group"),
            center_groups=meta.get("center_groups"),
            output_groups=meta.get("output_groups"),
            skip_cluster=bool(meta.get("skip_cluster")),
            tpr_path=meta.get("tpr"),
            ndx_path=meta.get("ndx"),
            gmx_executable=meta.get("gmx_executable"),
            cpptraj_executable=meta.get("cpptraj_executable"),
            output_format=meta.get("output_format") or "same",
            stride=_normalize_stride(meta.get("stride")),
            file_strides=meta.get("file_strides"),
            log_fn=log_fn,
            progress_fn=progress_fn,
            job_dir=str(path),
        )
        if _is_cancel_requested(path):
            raise JobCancelled("Job cancelled by user")
        _mark_step(path, "Finalize", completed=True)
        ok_any = any(o.get("ok") for o in result.get("outputs") or [])
        if not ok_any:
            err = next(
                (o.get("error") for o in result.get("outputs") or [] if o.get("error")),
                "All trajectories failed",
            )
            _update_status(
                path,
                status="error",
                error=err,
                end_time=_now_iso(),
                outputs=result.get("outputs") or [],
                engine=result.get("engine"),
                method=result.get("method"),
            )
            _append_log(path, f"FAILED: {err}")
            exit_code = 1
        else:
            partial = [o for o in result.get("outputs") or [] if not o.get("ok")]
            _update_status(
                path,
                status="completed",
                error=(f"{len(partial)} trajectory(ies) failed" if partial else None),
                end_time=_now_iso(),
                outputs=result.get("outputs") or [],
                engine=result.get("engine"),
                method=result.get("method"),
            )
            _append_log(path, "Job completed")
            exit_code = 0
    except JobCancelled:
        _update_status(
            path,
            status="cancelled",
            error="Cancelled by user",
            end_time=_now_iso(),
        )
        _append_log(path, "CANCELLED by user")
        exit_code = 130
    except Exception as exc:
        if _is_cancel_requested(path):
            _update_status(
                path,
                status="cancelled",
                error="Cancelled by user",
                end_time=_now_iso(),
            )
            _append_log(path, "CANCELLED by user")
            exit_code = 130
        else:
            _update_status(
                path,
                status="error",
                error=str(exc),
                end_time=_now_iso(),
            )
            _append_log(path, f"FAILED: {exc}")
            logger.exception("Fix PBC job failed: %s", path)
            exit_code = 1
    finally:
        (path / "process.pid").unlink(missing_ok=True)
    return exit_code


def start_fix_pbc_job(
    topology_file: str,
    trajectory_files: list[str],
    output_parent: str,
    *,
    engine: str = "auto",
    center_selection: str = "protein",
    center_group: Optional[str] = None,
    output_group: Optional[str] = None,
    center_groups: Optional[list[str]] = None,
    output_groups: Optional[list[str]] = None,
    skip_cluster: bool = False,
    tpr_path: Optional[str] = None,
    ndx_path: Optional[str] = None,
    gmx_executable: Optional[str] = None,
    cpptraj_executable: Optional[str] = None,
    output_format: str = "same",
    stride: int = 1,
    file_strides: Optional[dict[str, int]] = None,
    job_name: Optional[str] = None,
) -> dict:
    """
    Create a job directory and run Fix PBC in a detached worker process.

    The worker survives GUI/backend close (Builder/MemPro-style). Cancel kills
    the worker process group via ``process.pid``.

    ``output_parent`` is typically the project working directory.
    ``job_name`` is the folder name created directly under it (e.g. ``05_tools``
    or a user-chosen name). If that folder already exists, a numeric suffix is
    added (``name_2``, …) so concurrent jobs do not collide.

    Returns immediately with ``job_dir``. Progress is written to
    ``status.json`` and ``logs/fix_pbc.log`` for polling via the GUI.
    """
    parent = Path(output_parent).expanduser().resolve()
    parent.mkdir(parents=True, exist_ok=True)

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    raw_name = (job_name or "").strip() or f"fix_pbc_{stamp}"
    # Keep a single path segment — no nested "parent/child" from the name field.
    base = Path(raw_name).name.replace("\\", "_").replace("/", "_") or f"fix_pbc_{stamp}"
    job_dir = parent / base
    n = 2
    while job_dir.exists():
        job_dir = parent / f"{base}_{n}"
        n += 1
    job_dir.mkdir(parents=True, exist_ok=False)
    (job_dir / "outputs").mkdir()
    (job_dir / "logs").mkdir()

    stride_n = _normalize_stride(stride)
    stride_map = _resolve_file_strides(trajectory_files, file_strides, stride_n)
    detected = detect_pbc_engine(
        topology_file,
        trajectory_files,
        engine_hint=engine,
        tpr_path=tpr_path,
        ndx_path=ndx_path,
    )
    steps = ["Detect engine"] + [f"Fix {Path(p).name}" for p in trajectory_files] + [
        "Finalize"
    ]

    center_sel = _normalize_group_name_list(center_groups, center_group)
    output_sel = _normalize_group_name_list(output_groups, output_group)
    meta = {
        "type": "fix_pbc",
        "engine": detected["engine"],
        "engine_hint": engine,
        "method": detected["method"],
        "topology": str(Path(topology_file).expanduser().resolve()),
        "trajectories": [str(Path(p).expanduser().resolve()) for p in trajectory_files],
        "center_selection": center_selection,
        "center_group": center_sel[0] if len(center_sel) == 1 else center_group,
        "output_group": output_sel[0] if len(output_sel) == 1 else output_group,
        "center_groups": center_sel,
        "output_groups": output_sel,
        "skip_cluster": bool(skip_cluster),
        "tpr": tpr_path or detected.get("tpr"),
        "ndx": ndx_path or detected.get("ndx"),
        "gmx_executable": gmx_executable,
        "cpptraj_executable": cpptraj_executable,
        "output_format": output_format,
        "stride": stride_n,
        "file_strides": stride_map,
    }
    _write_json(job_dir / "tools_job.json", meta)
    status_payload = {
        "status": "running",
        "current_step": 0,
        "steps": steps,
        "steps_completed": [],
        "error": None,
        "start_time": _now_iso(),
        "end_time": None,
        "outputs": [],
        "engine": detected["engine"],
        "method": detected["method"],
        "config": meta,
    }
    _append_log(job_dir, f"Job started: {job_dir.name}")
    _append_log(job_dir, f"Engine detect: {detected}")
    if center_sel:
        label = "+".join(center_sel) if len(center_sel) > 1 else center_sel[0]
        _append_log(
            job_dir,
            f"Center group(s): {label}"
            + (f" → GW_CENTER" if len(center_sel) > 1 else ""),
        )
    if output_sel:
        label = "+".join(output_sel) if len(output_sel) > 1 else output_sel[0]
        _append_log(
            job_dir,
            f"Output group(s): {label}"
            + (f" → GW_OUTPUT" if len(output_sel) > 1 else ""),
        )
    if skip_cluster:
        _append_log(job_dir, "trjconv cluster step: disabled")
    if any(n > 1 for n in stride_map.values()):
        for p, n in stride_map.items():
            if n > 1:
                _append_log(job_dir, f"Stride {Path(p).name}: every {n} frame(s)")

    # Detached worker (new session) so the job survives Electron/backend exit.
    worker_log = job_dir / "logs" / "worker_stderr.log"
    try:
        err_fh = open(worker_log, "a", encoding="utf-8")
    except OSError:
        err_fh = subprocess.DEVNULL
    try:
        proc = subprocess.Popen(
            [sys.executable, "-m", "gatewizard.utils.fix_pbc_worker", str(job_dir)],
            cwd=str(job_dir),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=err_fh,
            start_new_session=True,
            env=os.environ.copy(),
        )
    finally:
        if err_fh is not subprocess.DEVNULL:
            try:
                err_fh.close()
            except Exception:
                pass

    (job_dir / "process.pid").write_text(str(proc.pid), encoding="utf-8")
    _write_json(job_dir / "status.json", status_payload)
    _append_log(job_dir, f"Detached worker pid={proc.pid}")

    return {
        "success": True,
        "job_dir": str(job_dir),
        "engine": detected["engine"],
        "method": detected["method"],
        "pid": proc.pid,
        "message": f"Fix PBC started ({detected['method']})",
        "detect": detected,
    }


def _reconcile_stale_running_job(job_dir: Path, status: dict) -> dict:
    """Keep running jobs alive unless the worker is gone and logs have gone quiet.

    A scan often races the PID file (nested output folders under the working
    directory are found immediately). Do not mark interrupted in that window,
    and revive a false interrupt if the log is still growing.
    """
    status = dict(status)
    state = str(status.get("status") or "").lower()
    alive = _worker_appears_alive(job_dir, status)

    if _looks_like_false_interrupt(status) and alive:
        status["status"] = "running"
        status["error"] = None
        status["end_time"] = None
        try:
            _write_json(job_dir / "status.json", status)
            _append_log(job_dir, "Restored running (worker still active)")
        except OSError:
            pass
        return status

    if state != "running":
        return status

    if alive:
        return status

    if _is_cancel_requested(job_dir):
        status["status"] = "cancelled"
        status["error"] = status.get("error") or "Cancelled by user"
    else:
        status["status"] = "error"
        status["error"] = status.get("error") or _INTERRUPT_ERROR
    status["end_time"] = status.get("end_time") or _now_iso()
    try:
        _write_json(job_dir / "status.json", status)
        _append_log(
            job_dir,
            f"Marked {status['status']} (worker pid={_read_pid_file(job_dir)} not running)",
        )
    except OSError:
        pass
    return status


def scan_tools_jobs(directory: str) -> list[dict]:
    """Scan for Tools job directories containing ``tools_job.json`` + ``status.json``.

    Searches the working directory (and one nesting level) so jobs are found
    after app restart regardless of the Output folder field value.
    """
    base = Path(directory).expanduser().resolve()
    if not base.is_dir():
        return []
    found: list[dict] = []
    # Job dir = workingDir/<name>/  (current)
    # Legacy: workingDir/<parent>/fix_pbc_*/ 
    candidates = list(base.glob("*/tools_job.json")) + list(
        base.glob("*/*/tools_job.json")
    )
    # Also accept tools_job.json directly in *base* if someone pointed WD at a job
    if (base / "tools_job.json").is_file():
        candidates.insert(0, base / "tools_job.json")

    seen: set[str] = set()
    for meta_path in candidates:
        job_dir = meta_path.parent
        key = str(job_dir.resolve())
        if key in seen:
            continue
        seen.add(key)
        status_path = job_dir / "status.json"
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            status = (
                json.loads(status_path.read_text(encoding="utf-8"))
                if status_path.is_file()
                else {}
            )
        except Exception:
            continue
        status = _reconcile_stale_running_job(job_dir, status)
        center_sel = _normalize_group_name_list(
            meta.get("center_groups"), meta.get("center_group")
        )
        output_sel = _normalize_group_name_list(
            meta.get("output_groups"), meta.get("output_group")
        )
        center_label = (
            (f"GW_CENTER = {'+'.join(center_sel)}" if len(center_sel) > 1 else (center_sel[0] if center_sel else None))
        )
        output_label = (
            (f"GW_OUTPUT = {'+'.join(output_sel)}" if len(output_sel) > 1 else (output_sel[0] if output_sel else None))
        )
        # Prefer labels from completed traj outputs when present
        for out in status.get("outputs") or []:
            if out.get("center_label"):
                center_label = out["center_label"]
            if out.get("output_label"):
                output_label = out["output_label"]
            break
        found.append(
            {
                "job_dir": str(job_dir),
                "name": job_dir.name,
                "type": meta.get("type", "fix_pbc"),
                "engine": status.get("engine") or meta.get("engine"),
                "method": status.get("method") or meta.get("method"),
                "status": status.get("status", "unknown"),
                "current_step": status.get("current_step", 0),
                "steps": status.get("steps", []),
                "steps_completed": status.get("steps_completed", []),
                "error": status.get("error"),
                "start_time": status.get("start_time"),
                "end_time": status.get("end_time"),
                "outputs": status.get("outputs", []),
                "center_groups": center_sel,
                "output_groups": output_sel,
                "center_label": center_label,
                "output_label": output_label,
                "skip_cluster": bool(meta.get("skip_cluster")),
            }
        )
    found.sort(key=lambda j: j.get("start_time") or "", reverse=True)
    return found
