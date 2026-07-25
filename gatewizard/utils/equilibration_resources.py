"""
Read and aggregate CPU / GPU settings for equilibration job folders.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from gatewizard.utils.equilibration_resume import _parse_namd_protocols_from_script


def _aggregate_stage_resources(stage_items: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Summarize cpu_cores / gpu settings across protocol stages."""
    cpus: List[int] = []
    gpu_ids: List[int] = []
    num_gpus_vals: List[int] = []
    use_gpu = False

    for stage in stage_items:
        if not isinstance(stage, dict):
            continue
        cpus.append(int(stage.get("cpu_cores") or 1))
        if stage.get("use_gpu"):
            use_gpu = True
            gpu_ids.append(int(stage.get("gpu_id") or 0))
            num_gpus_vals.append(int(stage.get("num_gpus") or 1))

    return {
        "use_gpu": use_gpu,
        "cpu_cores_min": min(cpus) if cpus else None,
        "cpu_cores_max": max(cpus) if cpus else None,
        "gpu_id_min": min(gpu_ids) if gpu_ids else None,
        "gpu_id_max": max(gpu_ids) if gpu_ids else None,
        "num_gpus": max(num_gpus_vals) if num_gpus_vals else (1 if use_gpu else 0),
        "platform": None,
    }


def resolve_compute_resources_from_stages(
    stage_items: Optional[List[Dict[str, Any]]],
) -> Dict[str, Any]:
    """Flat CPU/GPU settings for run-script generation (max cores / first GPU range)."""
    summary = _aggregate_stage_resources(stage_items or [])
    return {
        "cpu_cores": int(summary["cpu_cores_max"] or summary["cpu_cores_min"] or 1),
        "use_gpu": bool(summary["use_gpu"]),
        "gpu_id": int(summary["gpu_id_min"] if summary["gpu_id_min"] is not None else 0),
        "num_gpus": int(summary["num_gpus"] or (1 if summary["use_gpu"] else 0)),
        "platform": summary.get("platform"),
    }


def resolve_compute_resources_from_eq_dir(eq_dir: Path) -> Dict[str, Any]:
    """Load flat CPU/GPU settings from ``equilibration_resources.json`` if present."""
    eq_dir = Path(eq_dir)
    path = eq_dir / "equilibration_resources.json"
    if not path.is_file():
        return {
            "cpu_cores": 1,
            "use_gpu": False,
            "gpu_id": 0,
            "num_gpus": 0,
            "platform": None,
        }
    try:
        data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except (json.JSONDecodeError, OSError):
        data = {}
    if not isinstance(data, dict):
        data = {}
    cpu = data.get("cpu_cores_max")
    if cpu is None:
        cpu = data.get("cpu_cores_min")
    if cpu is None:
        cpu = data.get("cpu_cores")
    gpu_id = data.get("gpu_id_min")
    if gpu_id is None:
        gpu_id = data.get("gpu_id")
    use_gpu = bool(data.get("use_gpu"))
    num_gpus = data.get("num_gpus")
    if num_gpus is None:
        num_gpus = 1 if use_gpu else 0
    return {
        "cpu_cores": int(cpu or 1),
        "use_gpu": use_gpu,
        "gpu_id": int(gpu_id if gpu_id is not None else 0),
        "num_gpus": int(num_gpus or 0),
        "platform": data.get("platform"),
    }


def write_equilibration_resources(
    eq_dir: Path,
    engine: str,
    stages: List[Dict[str, Any]],
    *,
    openmm_platform: Optional[str] = None,
) -> Path:
    """Persist resource summary next to run_equilibration.sh (on input generation)."""
    eq_dir = Path(eq_dir)
    summary = _aggregate_stage_resources(stages)
    summary["engine"] = engine.lower().strip()
    if openmm_platform:
        summary["platform"] = openmm_platform
        summary["use_gpu"] = openmm_platform.upper() != "CPU"
    path = eq_dir / "equilibration_resources.json"
    path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return path


def _openmm_platform_from_script(script_text: str) -> str:
    if re.search(r'^PLATFORM="[^"]+"', script_text, re.MULTILINE):
        match = re.search(r'^PLATFORM="([^"]+)"', script_text, re.MULTILINE)
        if match and match.group(1):
            return match.group(1)
    return "auto"


def infer_equilibration_resources(eq_dir: Path, engine: str) -> Dict[str, Any]:
    """
    Return CPU / GPU resource summary for a job directory.

    Prefers ``equilibration_resources.json``, then engine-specific fallbacks.
    """
    eq_dir = Path(eq_dir)
    engine = (engine or "").strip().lower()
    resources_file = eq_dir / "equilibration_resources.json"
    if resources_file.is_file():
        try:
            data = json.loads(resources_file.read_text(encoding="utf-8", errors="replace"))
            if isinstance(data, dict):
                return data
        except (json.JSONDecodeError, OSError):
            pass

    if engine == "namd":
        summary_file = eq_dir / "protocol_summary.json"
        if summary_file.is_file():
            try:
                summary = json.loads(summary_file.read_text(encoding="utf-8", errors="replace"))
                stages = summary.get("stages") or {}
                if isinstance(stages, dict):
                    items = list(stages.values())
                else:
                    items = list(stages)
                if items:
                    out = _aggregate_stage_resources(items)
                    out["engine"] = "namd"
                    return out
            except (json.JSONDecodeError, OSError, TypeError, ValueError):
                pass
        script = eq_dir / "run_equilibration.sh"
        if script.is_file():
            try:
                text = script.read_text(encoding="utf-8", errors="replace")
                _, protocols = _parse_namd_protocols_from_script(text)
                if protocols:
                    out = _aggregate_stage_resources(list(protocols.values()))
                    out["engine"] = "namd"
                    return out
            except OSError:
                pass

    if engine == "openmm":
        script = eq_dir / "run_equilibration.sh"
        platform = "auto"
        if script.is_file():
            try:
                platform = _openmm_platform_from_script(
                    script.read_text(encoding="utf-8", errors="replace")
                )
            except OSError:
                pass
        return {
            "engine": "openmm",
            "use_gpu": platform.upper() not in {"", "AUTO", "CPU", "REFERENCE"},
            "platform": platform,
            "cpu_cores_min": None,
            "cpu_cores_max": None,
            "gpu_id_min": None,
            "gpu_id_max": None,
            "num_gpus": 0,
        }

    return {
        "engine": engine,
        "use_gpu": None,
        "platform": None,
        "cpu_cores_min": None,
        "cpu_cores_max": None,
        "gpu_id_min": None,
        "gpu_id_max": None,
        "num_gpus": None,
    }
