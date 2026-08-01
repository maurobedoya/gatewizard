"""Mid-run progress sync from node-local scratch back to the submit directory."""

from __future__ import annotations

import re
import shlex
from pathlib import Path
from typing import Optional, Tuple

from gatewizard.utils.cluster.ssh import ClusterSSHError, run_remote

# Lightweight progress artifacts copied during mid-run sync / scratch pull.
_SCRATCH_PROGRESS_INCLUDES = (
    "step*.log step*.xst step*.mdinfo step*.mdout step*.rst7 step*.rst "
    "step*_minimization.log equilibration_background.log *.out *.err "
    "run_equilibration.sh run_equilibration_cluster.sh openmm_nvt.*.out"
)

# rsync --include list for Watching (logs only — not trajectories / full Pull).
PROGRESS_RSYNC_FILTERS = [
    "*/",
    "step*.log",
    "step*_minimization.log",
    "step*.mdout",
    "step*.mdinfo",
    "step*.xst",
    "equilibration_background.log",
    "*.out",
    "*.err",
]


def resolve_compute_node(session_id: str, job_id: str) -> str:
    """Return the first allocated node name for a Slurm job, or ``\"\"``."""
    jid = shlex.quote(str(job_id))
    _rc, out, err = run_remote(
        session_id,
        f"squeue -j {jid} -h -o '%N' 2>/dev/null || true",
        timeout=30,
    )
    text = (out or err or "").strip()
    if not text:
        return ""
    first = text.split(",")[0].strip()
    if "[" in first:
        _rc2, out2, _ = run_remote(
            session_id,
            f"scontrol show hostnames {shlex.quote(first)} 2>/dev/null | head -1 || true",
            timeout=30,
        )
        host = (out2 or "").strip().splitlines()
        return host[0].strip() if host else ""
    return first


def sync_scratch_progress_to_submit(
    session_id: str,
    *,
    job_id: str,
    node: str,
    scratch_root: str,
    remote_submit_dir: str,
) -> Tuple[bool, str]:
    """Copy lightweight progress files from node scratch → submit directory.

    On many clusters ``/scratch/$USER/$SLURM_JOB_ID`` exists only on the
    compute node. We SSH from the login node to that host and stream a tar of
    log/progress files into the submit directory (which Pull then downloads).
    """
    node = (node or "").strip()
    if not node or not job_id:
        return False, "no compute node / job id"
    scratch_root = (scratch_root or "").rstrip("/")
    if not scratch_root:
        return False, "scratch_root empty"
    remote_submit_dir = (remote_submit_dir or "").rstrip("/")
    if not remote_submit_dir:
        return False, "remote submit dir empty"

    scratch_job = f"{scratch_root}/{job_id}"
    node_q = shlex.quote(node)
    scratch_q = shlex.quote(scratch_job)
    submit_q = shlex.quote(remote_submit_dir)

    # Stream selected files login ← compute. Exclude heavy trajectories.
    cmd = (
        f"ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new {node_q} "
        f"'if [ ! -d {scratch_q} ]; then echo GW_MIDRUN_MISSING; exit 0; fi; "
        f"cd {scratch_q} && tar czf - "
        f"--ignore-failed-read "
        f"{_SCRATCH_PROGRESS_INCLUDES} "
        f"2>/dev/null' "
        f"| tar xzf - -C {submit_q} 2>&1; "
        f"echo GW_MIDRUN_DONE"
    )
    try:
        rc, out, err = run_remote(session_id, cmd, timeout=180)
    except ClusterSSHError as ex:
        return False, str(ex)
    text = "\n".join(p for p in (out, err) if p and p.strip()).strip()
    if "GW_MIDRUN_MISSING" in text:
        return False, f"scratch not found on {node}:{scratch_job}"
    if rc != 0 and "GW_MIDRUN_DONE" not in text:
        return False, text or f"mid-run sync failed (exit {rc})"
    return True, f"synced scratch logs from {node}:{scratch_job}"


def expand_scratch_job_dir(scratch_root: str, job_id: str) -> str:
    root = (scratch_root or "").rstrip("/")
    return f"{root}/{job_id}" if root else ""


def remote_path_is_dir(session_id: str, path: str) -> bool:
    """Return True when *path* exists as a directory on the login node."""
    path = (path or "").strip().rstrip("/")
    if not path:
        return False
    path_q = shlex.quote(path)
    _rc, out, _err = run_remote(
        session_id,
        f"test -d {path_q} && echo 1 || echo 0",
        timeout=20,
    )
    return (out or "").strip().endswith("1")


def resolve_slurm_workdir(session_id: str, job_id: str) -> str:
    """Return Slurm WorkDir for *job_id* from ``sacct``, or ``\"\"``."""
    jid = shlex.quote(str(job_id))
    _rc, out, _err = run_remote(
        session_id,
        f"sacct -j {jid} -n -X -o WorkDir -P 2>/dev/null | head -1",
        timeout=30,
    )
    line = (out or "").strip().splitlines()[0] if (out or "").strip() else ""
    if not line or line.lower() == "workdir":
        return ""
    return line.strip().rstrip("/")


def build_remote_submit_path(
    submit_root: str,
    job_folder: str,
    *,
    username: str = "",
) -> str:
    """Build ``submit_root/job_folder`` with ``$USER``-style expansion."""
    from gatewizard.utils.cluster.paths import expand_remote_path, join_remote

    root = expand_remote_path(submit_root or "", username=username).rstrip("/")
    folder = (job_folder or "").strip().strip("/")
    if not root or not folder:
        return ""
    return join_remote(root, folder)


def resolve_remote_job_dir(
    session_id: str,
    *,
    stored_path: str = "",
    job_id: str = "",
    submit_root: str = "",
    username: str = "",
    job_folder: str = "",
) -> Tuple[str, str, list[str]]:
    """Return ``(path, source, tried)`` — first candidate that exists on the login node."""
    tried: list[str] = []
    seen: set[str] = set()
    candidates: list[tuple[str, str]] = []

    def add(path: str, label: str) -> None:
        p = (path or "").strip().rstrip("/")
        if not p or p in seen:
            return
        seen.add(p)
        candidates.append((p, label))

    add(stored_path, "stored")
    if job_id:
        add(resolve_slurm_workdir(session_id, str(job_id)), "sacct WorkDir")
    add(
        build_remote_submit_path(submit_root, job_folder, username=username),
        "profile submit_root",
    )

    for path, source in candidates:
        tried.append(path)
        if remote_path_is_dir(session_id, path):
            return path, source, tried

    if candidates:
        path, source = candidates[0]
        return path, f"unverified ({source})", tried
    return "", "", tried


def parse_scratch_workdir_from_slurm(local_dir: Path) -> Tuple[str, str]:
    """Parse ``(scratch_root, workdir_strategy)`` from a local ``run_equilibration.slurm``."""
    slurm = Path(local_dir) / "run_equilibration.slurm"
    if not slurm.is_file():
        return "", ""
    try:
        text = slurm.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return "", ""
    m = re.search(r'workdir="([^"$]+)/\$SLURM_JOB_ID"', text)
    if m:
        return m.group(1).rstrip("/"), "scratch_job_id"
    m = re.search(r'workdir="([^"]+)/([^"/\$]+)"', text)
    if m and "SLURM" not in m.group(2):
        return m.group(1).rstrip("/"), "scratch_named"
    return "", ""


def stage_scratch_to_login(
    session_id: str,
    *,
    node: str,
    scratch_job_dir: str,
    staging_dir: str,
    full: bool = False,
) -> Tuple[bool, str]:
    """Copy a scratch job directory to a login-node staging folder via compute SSH."""
    node = (node or "").strip()
    scratch_job_dir = (scratch_job_dir or "").rstrip("/")
    staging_dir = (staging_dir or "").rstrip("/")
    if not node or not scratch_job_dir or not staging_dir:
        return False, "node / scratch / staging path missing"

    node_q = shlex.quote(node)
    scratch_q = shlex.quote(scratch_job_dir)
    staging_q = shlex.quote(staging_dir)
    if full:
        tar_cmd = (
            f"cd {scratch_q} && tar czf - --ignore-failed-read "
            f"--exclude='*.dcd' --exclude='*.xtc' --exclude='*.nc' "
            f"--exclude='*.cpt' --exclude='*.chk' . 2>/dev/null"
        )
    else:
        tar_cmd = (
            f"cd {scratch_q} && tar czf - --ignore-failed-read "
            f"{_SCRATCH_PROGRESS_INCLUDES} 2>/dev/null"
        )

    cmd = (
        f"rm -rf {staging_q} && mkdir -p {staging_q} && "
        f"ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new {node_q} "
        f"'if [ ! -d {scratch_q} ]; then echo GW_PULL_MISSING; exit 0; fi; "
        f"{tar_cmd}' "
        f"| tar xzf - -C {staging_q} 2>&1; "
        f"echo GW_PULL_DONE"
    )
    try:
        rc, out, err = run_remote(session_id, cmd, timeout=600)
    except ClusterSSHError as ex:
        return False, str(ex)
    text = "\n".join(p for p in (out, err) if p and p.strip()).strip()
    if "GW_PULL_MISSING" in text:
        return False, f"scratch not found on {node}:{scratch_job_dir}"
    if rc != 0 and "GW_PULL_DONE" not in text:
        return False, text or f"scratch staging failed (exit {rc})"
    return True, f"staged from {node}:{scratch_job_dir}"


def resolve_scratch_job_dir(
    *,
    job_id: str,
    execution: dict,
    profile_scratch_root: str = "",
    local_dir: Optional[Path] = None,
) -> str:
    """Best-effort scratch job path for fallback pulls."""
    execution = execution or {}
    strategy = str(execution.get("workdir_strategy") or "").strip()
    scratch_root = str(execution.get("scratch_root") or profile_scratch_root or "").strip()
    if not scratch_root and local_dir is not None:
        parsed_root, parsed_strategy = parse_scratch_workdir_from_slurm(local_dir)
        if parsed_root:
            scratch_root = parsed_root
        if parsed_strategy and not strategy:
            strategy = parsed_strategy
    if not scratch_root or scratch_root.startswith("$"):
        return ""
    scratch_root = scratch_root.rstrip("/")
    if strategy == "scratch_named":
        folder = (
            str(execution.get("job_folder_name") or "").strip()
            or (Path(local_dir).name if local_dir else "")
        )
        return f"{scratch_root}/{folder}" if folder else ""
    return expand_scratch_job_dir(scratch_root, job_id)
