"""Detect MD failures from stage logs and Slurm batch outputs."""

from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional, Tuple

# Real failure markers — NOT bare "error" (GROMACS warnings say "spelling error").
EQ_FAILURE_RE = re.compile(
    r"(?i)"
    r"("
    r"fatal error"
    r"|error in user input:"
    r"|error in stage\b"
    r"|error: namd executable"
    r"|cuda driver is a stub library"
    r"|cuda initialization error"
    r"|minimisation failed"
    r"|minimization failed"
    r"|production failed"
    r"|stage\s+\d+.*failed"
    r"|equilibration failed"
    r"|command not found"
    r"|segmentation fault"
    r")"
)

# Local metadata that must not be overwritten by a remote pull.
PULL_PRESERVE_EXCLUDES = [
    "equilibration_job.json",
    "equilibration.pid",
]


def failure_line_from_text(text: str) -> Optional[str]:
    """Return the first line that indicates a real MD failure, else None."""
    if not text:
        return None
    for line in text.splitlines():
        if EQ_FAILURE_RE.search(line):
            return line.strip() or line
    return None


def _tail_text(path: Path, *, max_bytes: int = 64_000) -> str:
    try:
        raw = path.read_bytes()
    except OSError:
        return ""
    if len(raw) > max_bytes:
        raw = raw[-max_bytes:]
    return raw.decode("utf-8", errors="replace")


def iter_failure_log_candidates(eq_dir: Path) -> List[Path]:
    """Prefer stage engine logs, then Slurm outs, then the local background log."""
    eq_dir = Path(eq_dir)
    found: List[Path] = []
    seen = set()

    def _add(paths) -> None:
        for p in paths:
            try:
                key = p.resolve()
            except OSError:
                key = p
            if key in seen or not p.is_file():
                continue
            seen.add(key)
            found.append(p)

    _add(sorted(eq_dir.glob("step*_equilibration*.log")))
    _add(sorted(eq_dir.glob("step*_production*.log")))
    _add(sorted(eq_dir.glob("step*.log")))
    _add(sorted(eq_dir.glob("step*.mdout")))
    _add(sorted(eq_dir.glob("*.out")))  # Slurm stdout (often empty on early fail)
    _add(sorted(eq_dir.glob("*.err")))
    bg = eq_dir / "equilibration_background.log"
    if bg.is_file():
        _add([bg])
    return found


def find_equilibration_failure(
    eq_dir: Path,
    *,
    newer_than: Optional[float] = None,
) -> Tuple[Optional[str], Optional[str]]:
    """Scan job folder for a failure message.

    Returns ``(message, source_filename)``. Empty Slurm ``.out`` files are
    ignored so stage logs (where NAMD writes FATAL ERROR) win.

    If ``newer_than`` is a UNIX mtime, logs older than that timestamp are
    ignored (used after cluster resubmit so prior FATAL lines do not stick).
    """
    eq_dir = Path(eq_dir)
    for path in iter_failure_log_candidates(eq_dir):
        # Empty Slurm outs are common when the batch wrapper barely started.
        try:
            st = path.stat()
            if st.st_size == 0:
                continue
            if newer_than is not None and st.st_mtime < float(newer_than):
                continue
        except OSError:
            continue
        msg = failure_line_from_text(_tail_text(path))
        if msg:
            return msg, path.name
    return None, None


def parse_submitted_at(execution: Optional[dict]) -> Optional[float]:
    """Return UNIX timestamp from ``execution.submitted_at`` ISO string, if any."""
    if not isinstance(execution, dict):
        return None
    raw = execution.get("submitted_at")
    if not raw:
        return None
    try:
        from datetime import datetime

        text = str(raw).strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        return datetime.fromisoformat(text).timestamp()
    except Exception:
        return None


def remote_job_is_active(execution: Optional[dict]) -> bool:
    """True when Slurm still has the job in a live queue state."""
    if not isinstance(execution, dict) or execution.get("mode") != "remote":
        return False
    state = str(execution.get("last_remote_state") or "").upper()
    return state in {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING", "REQUEUED"}


# Outputs left by a previous MD attempt — moved aside on cluster resubmit.
_ARCHIVE_GLOBS = (
    "step*.log",
    "step*.mdout",
    "step*.mdinfo",
    "step*.edr",
    "step*.xtc",
    "step*.trr",
    "step*.cpt",
    "step*.gro",
    "step*.coor",
    "step*.vel",
    "step*.xsc",
    "step*.xst",
    "step*.dcd",
    "equilibration_background.log",
    "equilibration_start_time.txt",
    "*.out",
    "*.err",
)


def archive_previous_run_outputs(eq_dir: Path, *, engine: Optional[str] = None) -> int:
    """Move prior MD outputs into ``_previous_cluster_run/<stamp>/``.

    Keeps inputs (``.conf`` / ``.mdp`` / ``.tpr`` / scripts) so a resubmit
    starts clean for status UI and remote upload.

    When ``engine`` is set (or inferred), checkpoint files for completed
    stages are kept so ``RESUME=1`` can skip finished stages on resubmit.
    """
    from datetime import datetime

    from gatewizard.utils.equilibration_resume import resume_checkpoint_paths

    eq_dir = Path(eq_dir)
    if engine is None:
        try:
            from gatewizard.utils.equilibration_job_metadata import (
                infer_equilibration_job_metadata,
            )

            meta = infer_equilibration_job_metadata(eq_dir, heal=False)
            engine = meta.get("engine") if isinstance(meta, dict) else None
        except Exception:
            engine = None
    preserve = resume_checkpoint_paths(eq_dir, engine or "")
    moved = 0
    to_move: List[Path] = []
    for pattern in _ARCHIVE_GLOBS:
        for path in eq_dir.glob(pattern):
            if not path.is_file():
                continue
            # Never archive the archive tree itself
            if "_previous_cluster_run" in path.parts:
                continue
            if preserve:
                try:
                    if path.resolve() in preserve:
                        continue
                except OSError:
                    if path in preserve:
                        continue
            to_move.append(path)
    if not to_move:
        return 0
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dest_root = eq_dir / "_previous_cluster_run" / stamp
    dest_root.mkdir(parents=True, exist_ok=True)
    for path in to_move:
        target = dest_root / path.name
        try:
            if target.exists():
                target = dest_root / f"{path.stem}_{moved}{path.suffix}"
            path.replace(target)
            moved += 1
        except OSError:
            continue
    return moved


def summarize_slurm_outputs(eq_dir: Path) -> Optional[str]:
    """Note empty/missing Slurm stdout when stage logs already show a failure."""
    eq_dir = Path(eq_dir)
    outs = sorted(eq_dir.glob("*.out"))
    if not outs:
        return None
    nonempty = [p for p in outs if p.is_file() and p.stat().st_size > 0]
    if nonempty:
        return None
    return (
        f"Slurm output empty ({outs[-1].name}); check stage logs "
        "(e.g. step1_equilibration.log) for the engine error"
    )
