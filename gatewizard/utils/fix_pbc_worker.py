"""Detached Fix PBC worker — survives GUI/backend close.

Usage:
    python -m gatewizard.utils.fix_pbc_worker <job_dir>
"""

from __future__ import annotations

import sys
from pathlib import Path


def main(argv: list[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    if not args:
        print("Usage: python -m gatewizard.utils.fix_pbc_worker <job_dir>", file=sys.stderr)
        return 2

    job_dir = Path(args[0]).expanduser().resolve()
    if not job_dir.is_dir():
        print(f"Job directory not found: {job_dir}", file=sys.stderr)
        return 1

    # Ensure PID file points at this worker (parent also writes it on spawn).
    try:
        (job_dir / "process.pid").write_text(str(__import__("os").getpid()), encoding="utf-8")
    except OSError:
        pass

    from gatewizard.utils.trajectory_tools import execute_fix_pbc_job

    return execute_fix_pbc_job(job_dir)


if __name__ == "__main__":
    raise SystemExit(main())
