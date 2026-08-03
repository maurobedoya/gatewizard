"""Helpers for dual local / cluster equilibration run scripts."""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Optional

LOCAL_RUN_SCRIPT = "run_equilibration.sh"
CLUSTER_RUN_SCRIPT = "run_equilibration_cluster.sh"

_CLUSTER_HEADER = (
    "# Cluster runner — use after module load; called by run_equilibration.slurm\n"
)


def cluster_engine_executable(
    engine: str,
    local_executable: Optional[str] = None,
    *,
    use_gpu: Optional[bool] = None,
) -> str:
    """Return a module-friendly command name for remote Slurm jobs.

    Absolute / WSL paths from the GUI Executable picker are replaced with the
    usual Environment Modules binaries. Bare command names are kept as-is.

    For Amber, ``use_gpu=True`` forces ``pmemd.cuda`` (cluster GPU submits),
    even when the local picker was plain ``pmemd``.
    """
    eng = (engine or "").strip().lower()
    local = (local_executable or "").strip()
    name = Path(local).name if local else ""
    bare = bool(
        local
        and ("/" not in local.replace("\\", "/"))
        and not (len(local) >= 2 and local[1] == ":")
    )

    if eng == "amber":
        return _cluster_amber_executable(local=local, name=name, use_gpu=use_gpu)

    # Already a bare command (no directory separators) → keep it.
    if bare and eng in {"namd", "gromacs", "openmm"}:
        return local

    if eng == "namd":
        if name.lower().startswith("namd"):
            return name  # namd3 / namd2
        return "namd3"
    if eng == "gromacs":
        return "gmx"
    if eng == "openmm":
        if name.lower().startswith("python"):
            return name
        return "python3"
    return name or local or "true"


def _cluster_amber_executable(
    *, local: str, name: str, use_gpu: Optional[bool]
) -> str:
    lower = (name or local or "").lower()
    wants_mpi = "mpi" in lower

    if use_gpu is True:
        return "pmemd.cuda.MPI" if wants_mpi else "pmemd.cuda"
    if use_gpu is False:
        if wants_mpi:
            return "pmemd.MPI"
        if name and "cuda" not in name.lower():
            return name  # sander / pmemd / …
        return "pmemd"

    # Infer from local executable name when use_gpu is unspecified.
    if "cuda" in lower:
        return "pmemd.cuda.MPI" if wants_mpi else "pmemd.cuda"
    if name:
        return name
    return "pmemd"


def stamp_cluster_run_script_header(content: str) -> str:
    """Insert the cluster-runner comment after the shebang line."""
    text = content or ""
    if "Cluster runner" in text:
        return text
    lines = text.splitlines(keepends=True)
    if not lines:
        return _CLUSTER_HEADER
    if lines[0].startswith("#!"):
        return lines[0] + _CLUSTER_HEADER + "".join(lines[1:])
    return _CLUSTER_HEADER + text


def write_cluster_run_script(eq_dir: Path, content: str) -> Path:
    """Write ``run_equilibration_cluster.sh`` with the standard header."""
    eq_dir = Path(eq_dir)
    path = eq_dir / CLUSTER_RUN_SCRIPT
    path.write_text(stamp_cluster_run_script_header(content), encoding="utf-8")
    try:
        path.chmod(path.stat().st_mode | 0o111)
    except OSError:
        os.chmod(path, 0o755)
    return path


def resolve_cluster_launch_script(eq_dir: Path) -> Path:
    """Prefer cluster runner; fall back to local for older job folders."""
    eq_dir = Path(eq_dir)
    cluster = eq_dir / CLUSTER_RUN_SCRIPT
    if cluster.is_file():
        return cluster
    return eq_dir / LOCAL_RUN_SCRIPT


def script_has_wsl_or_windows_path(text: str) -> bool:
    """True if a run script embeds a Windows/WSL absolute path."""
    if not text:
        return False
    lower = text.lower()
    if "/mnt/c/" in lower or "/mnt/d/" in lower:
        return True
    return bool(re.search(r"[A-Za-z]:\\", text))


def _parse_script_assignment(text: str, var: str) -> str:
    match = re.search(
        rf'^{re.escape(var)}="([^"]*)"', text, flags=re.MULTILINE
    )
    return match.group(1).strip() if match else ""


def ensure_amber_cluster_runner_for_gpus(eq_dir: Path, *, gpus: int) -> bool:
    """Rewrite Amber ``run_equilibration_cluster.sh`` for the submit GPU count.

    When ``gpus > 0``, both minimization and dynamics use ``pmemd.cuda``.
    When ``gpus == 0``, both use CPU ``pmemd``.

    Returns True if this looks like an Amber job folder and the cluster runner
    was updated (or already correct after rewrite).
    """
    eq_dir = Path(eq_dir)
    stage_stems = [p.stem for p in sorted(eq_dir.glob("step*.mdin"))]
    if not stage_stems:
        return False

    src = eq_dir / CLUSTER_RUN_SCRIPT
    if not src.is_file():
        src = eq_dir / LOCAL_RUN_SCRIPT
    if not src.is_file():
        return False

    text = src.read_text(encoding="utf-8", errors="replace")
    lower = text.lower()
    if "amber=" not in lower and "pmemd" not in lower and "sander" not in lower:
        return False

    from gatewizard.tools.equilibration import AmberEquilibrationManager
    from gatewizard.utils.equilibration_resources import (
        resolve_compute_resources_from_eq_dir,
    )

    want_gpu = int(gpus or 0) > 0
    compute = resolve_compute_resources_from_eq_dir(eq_dir)
    local_amber = _parse_script_assignment(text, "AMBER") or "pmemd"
    cluster_exe = cluster_engine_executable(
        "amber", local_amber, use_gpu=want_gpu
    )
    num_gpus = max(1, int(gpus)) if want_gpu else max(1, int(compute.get("num_gpus") or 1))

    manager = AmberEquilibrationManager(eq_dir)
    path = manager.generate_run_script(
        amber_dir=eq_dir,
        prmtop_name=_parse_script_assignment(text, "PRMTOP") or "system.prmtop",
        inpcrd_name=_parse_script_assignment(text, "INPCRD") or "system.inpcrd",
        stage_stems=stage_stems,
        amber_executable=cluster_exe,
        cpu_cores=compute.get("cpu_cores"),
        use_gpu=want_gpu,
        gpu_id=int(compute.get("gpu_id") or 0),
        num_gpus=num_gpus,
        script_filename=CLUSTER_RUN_SCRIPT,
    )
    path.write_text(
        stamp_cluster_run_script_header(path.read_text(encoding="utf-8")),
        encoding="utf-8",
    )
    return True
