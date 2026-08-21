"""Fix PBC job status: do not mark a live worker as interrupted."""

from __future__ import annotations

import json
from datetime import datetime, timedelta
from pathlib import Path

from gatewizard.utils import trajectory_tools as tt


def _write_status(job_dir: Path, **fields) -> dict:
    data = {
        "status": "running",
        "start_time": datetime.now().isoformat(),
        "steps": ["Detect engine", "Fix a.xtc"],
        "steps_completed": ["Detect engine"],
        "error": None,
    }
    data.update(fields)
    (job_dir / "status.json").write_text(json.dumps(data), encoding="utf-8")
    return data


def test_reconcile_keeps_running_without_pid_during_start_grace(tmp_path: Path):
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    status = _write_status(job_dir)
    out = tt._reconcile_stale_running_job(job_dir, status)
    assert out["status"] == "running"
    assert out.get("error") in (None, "")


def test_reconcile_keeps_running_when_log_is_fresh(tmp_path: Path):
    job_dir = tmp_path / "job"
    (job_dir / "logs").mkdir(parents=True)
    (job_dir / "logs" / "fix_pbc.log").write_text("[step] Fix a.xtc\n", encoding="utf-8")
    status = _write_status(
        job_dir,
        start_time=(datetime.now() - timedelta(minutes=10)).isoformat(),
    )
    out = tt._reconcile_stale_running_job(job_dir, status)
    assert out["status"] == "running"


def test_reconcile_marks_interrupted_when_worker_is_gone(tmp_path: Path):
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    status = _write_status(
        job_dir,
        start_time=(datetime.now() - timedelta(hours=2)).isoformat(),
    )
    out = tt._reconcile_stale_running_job(job_dir, status)
    assert out["status"] == "error"
    assert "exited unexpectedly" in (out.get("error") or "")
    assert (job_dir / "process.pid").exists() is False


def test_reconcile_revives_false_interrupt_when_log_grows(tmp_path: Path):
    job_dir = tmp_path / "job"
    (job_dir / "logs").mkdir(parents=True)
    (job_dir / "logs" / "fix_pbc.log").write_text("[done] Fix a.xtc\n", encoding="utf-8")
    status = _write_status(
        job_dir,
        status="error",
        error=tt._INTERRUPT_ERROR,
        start_time=(datetime.now() - timedelta(minutes=5)).isoformat(),
    )
    out = tt._reconcile_stale_running_job(job_dir, status)
    assert out["status"] == "running"
    assert out.get("error") is None
    saved = json.loads((job_dir / "status.json").read_text(encoding="utf-8"))
    assert saved["status"] == "running"


def test_mark_step_clears_false_interrupt(tmp_path: Path):
    job_dir = tmp_path / "job"
    job_dir.mkdir()
    _write_status(
        job_dir,
        status="error",
        error=tt._INTERRUPT_ERROR,
        steps=["Detect engine", "Fix a.xtc"],
        steps_completed=["Detect engine"],
    )
    tt._mark_step(job_dir, "Fix a.xtc", completed=True)
    saved = json.loads((job_dir / "status.json").read_text(encoding="utf-8"))
    assert saved["status"] == "running"
    assert saved.get("error") is None
    assert "Fix a.xtc" in saved["steps_completed"]
