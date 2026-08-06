"""
Read, resolve, and persist per-stage CPU / GPU settings for equilibration jobs.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

from gatewizard.utils.equilibration_resume import _parse_namd_protocols_from_script

RESOURCES_VERSION = 2

DEFAULT_COMPUTE_DEFAULTS: Dict[str, Any] = {
    "cpu_cores": 1,
    "gpu_id": 0,
    "num_gpus": 1,
    "use_gpu": True,
    "compute_target": "auto",
}

DEFAULT_MINIMIZATION_RESOURCES: Dict[str, Any] = {
    "cpu_cores": 6,
    "gpu_id": 0,
    "num_gpus": 0,
    "use_gpu": False,
}

DEFAULT_PRODUCTION_RESOURCES: Dict[str, Any] = {
    "cpu_cores": 1,
    "gpu_id": 0,
    "num_gpus": 1,
    "use_gpu": True,
}


def engine_resource_profile(engine: str) -> Dict[str, Any]:
    """Per-engine default compute settings for minimization, MD, and production.

    GROMACS: CPU minimization; equilibration and production use CPU×6 + GPU×1.
    Amber: entire equilibration on CPU×6; production uses CPU×1 + GPU×1 (pmemd.cuda).
    """
    engine = (engine or "").strip().lower()
    cpu_md = {
        "cpu_cores": 6,
        "gpu_id": 0,
        "num_gpus": 1,
        "use_gpu": True,
    }
    cpu_only_md = {
        "cpu_cores": 6,
        "gpu_id": 0,
        "num_gpus": 0,
        "use_gpu": False,
    }
    mini = dict(DEFAULT_MINIMIZATION_RESOURCES)
    prod_gpu = dict(DEFAULT_PRODUCTION_RESOURCES)

    if engine == "amber":
        return {
            "compute_defaults": {**cpu_only_md, "compute_target": "auto"},
            "minimization": mini,
            "equilibration": cpu_only_md,
            "production": prod_gpu,
        }
    if engine == "gromacs":
        return {
            "compute_defaults": {**cpu_md, "compute_target": "auto"},
            "minimization": mini,
            "equilibration": cpu_md,
            "production": dict(cpu_md),
        }
    # NAMD, OpenMM, and unknown engines: GPU for MD stages, lighter host CPU for prod
    return {
        "compute_defaults": {**cpu_md, "compute_target": "auto"},
        "minimization": mini,
        "equilibration": cpu_md,
        "production": prod_gpu,
    }


def infer_stage_kind(stage: Dict[str, Any]) -> str:
    """Return ``minimization``, ``production``, or ``equilibration``."""
    explicit = str(stage.get("stage_kind") or "").strip().lower()
    if explicit in {"minimization", "production", "equilibration"}:
        return explicit
    name = str(stage.get("name") or "").strip().lower()
    if name in {"minimization", "energy minimization", "energy_minimization"}:
        return "minimization"
    if name == "production":
        return "production"
    if int(stage.get("minimize_steps") or 0) > 0 and stage.get("time_ns", 0) in (
        0,
        0.0,
        None,
    ):
        return "minimization"
    return "equilibration"


def stage_stem_for_index(engine: str, index: int, stage_kind: str) -> str:
    """Map protocol stage index to on-disk stem (``step0_minimization``, …)."""
    engine = (engine or "").strip().lower()
    if stage_kind == "minimization":
        return "step0_minimization"
    if stage_kind == "production":
        return "step7_production"
    # equilibration stages after optional minimization occupy step1..step6
    return f"step{index}_equilibration" if engine == "namd" else f"step{index}_equilibration"


def _stems_for_stage_list(
    engine: str, stages: List[Dict[str, Any]]
) -> List[str]:
    """Return filesystem stems aligned with resolved protocol stages."""
    engine = (engine or "").strip().lower()
    stems: List[str] = []
    eq_num = 0
    for stage in stages:
        kind = infer_stage_kind(stage)
        if kind == "minimization":
            stems.append("step0_minimization")
        elif kind == "production":
            stems.append("step7_production")
        else:
            eq_num += 1
            stems.append(f"step{eq_num}_equilibration")
    return stems


def _merge_defaults(
    base: Dict[str, Any], override: Optional[Dict[str, Any]]
) -> Dict[str, Any]:
    out = dict(base)
    if override:
        for key, val in override.items():
            if val is not None:
                out[key] = val
    return out


def resolve_stage_resources(
    stage: Dict[str, Any],
    compute_defaults: Optional[Dict[str, Any]] = None,
    *,
    engine: str = "",
) -> Dict[str, Any]:
    """Resolve one stage's compute settings from inherit flag and kind defaults."""
    profile = engine_resource_profile(engine)
    defaults = _merge_defaults(DEFAULT_COMPUTE_DEFAULTS, profile["compute_defaults"])
    defaults = _merge_defaults(defaults, compute_defaults)
    kind = infer_stage_kind(stage)
    inherit = stage.get("resources_inherit")
    if inherit is None:
        has_explicit = any(
            stage.get(k) is not None for k in ("cpu_cores", "gpu_id", "num_gpus", "use_gpu")
        )
        inherit = kind == "equilibration" and not has_explicit

    if kind == "minimization":
        base = dict(profile["minimization"])
        inherit = False
    elif kind == "production":
        base = dict(profile["production"])
        if inherit:
            base.update(
                {
                    "cpu_cores": int(defaults.get("cpu_cores") or 1),
                    "gpu_id": int(defaults.get("gpu_id") or 0),
                    "num_gpus": int(defaults.get("num_gpus") or 1),
                    "use_gpu": bool(defaults.get("use_gpu")),
                }
            )
            inherit = False
    elif inherit:
        eq_base = profile["equilibration"]
        base = {
            "cpu_cores": int(defaults.get("cpu_cores") or eq_base["cpu_cores"]),
            "gpu_id": int(defaults.get("gpu_id") if defaults.get("gpu_id") is not None else eq_base["gpu_id"]),
            "num_gpus": int(defaults.get("num_gpus") if defaults.get("num_gpus") is not None else eq_base["num_gpus"]),
            "use_gpu": bool(defaults.get("use_gpu") if "use_gpu" in defaults else eq_base["use_gpu"]),
        }
    else:
        base = dict(defaults)

    resolved = dict(base)
    if not inherit:
        for key in ("cpu_cores", "gpu_id", "num_gpus", "use_gpu", "compute_target"):
            if key in stage and stage[key] is not None:
                resolved[key] = stage[key]

    if kind == "minimization":
        resolved["use_gpu"] = False
        resolved["num_gpus"] = 0

    resolved["cpu_cores"] = max(1, int(resolved.get("cpu_cores") or 1))
    resolved["gpu_id"] = int(resolved.get("gpu_id") or 0)
    if resolved.get("use_gpu"):
        resolved["num_gpus"] = max(1, int(resolved.get("num_gpus") or 1))
    else:
        resolved["num_gpus"] = 0

    resolved["stage_kind"] = kind
    resolved["name"] = stage.get("name") or kind.title()
    return resolved


def resolve_all_stage_resources(
    stages: List[Dict[str, Any]],
    compute_defaults: Optional[Dict[str, Any]] = None,
    *,
    engine: str = "",
    stems: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """Resolve compute settings for every protocol stage."""
    if stems is None:
        stems = _stems_for_stage_list(engine, stages)
    resolved: List[Dict[str, Any]] = []
    for idx, stage in enumerate(stages):
        if not isinstance(stage, dict):
            continue
        item = resolve_stage_resources(stage, compute_defaults, engine=engine)
        if idx < len(stems):
            item["stem"] = stems[idx]
        resolved.append(item)
    return resolved


def aggregate_slurm_resources(stage_items: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Slurm allocation: max CPU cores and max GPU count across stages."""
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
        "cpu_cores": max(cpus) if cpus else 1,
        "num_gpus": max(num_gpus_vals) if num_gpus_vals else (1 if use_gpu else 0),
        "use_gpu": use_gpu,
        "gpu_id_min": min(gpu_ids) if gpu_ids else None,
        "gpu_id_max": max(gpu_ids) if gpu_ids else None,
    }


def enrich_resources_display(data: Dict[str, Any]) -> Dict[str, Any]:
    """Add cpu_cores_min/max, gpu_id_min/max, and summary for job cards / form restore."""
    stages = data.get("stages") or []
    slurm = data.get("slurm") or {}
    cpus = [
        int(s.get("cpu_cores") or 1)
        for s in stages
        if isinstance(s, dict)
    ]
    if not cpus and slurm.get("cpu_cores"):
        cpus = [int(slurm["cpu_cores"])]
    if cpus:
        data["cpu_cores_min"] = min(cpus)
        data["cpu_cores_max"] = max(cpus)
    gpu_ids = [
        int(s.get("gpu_id") or 0)
        for s in stages
        if isinstance(s, dict) and s.get("use_gpu")
    ]
    if gpu_ids:
        data["gpu_id_min"] = min(gpu_ids)
        data["gpu_id_max"] = max(gpu_ids)
    elif slurm.get("gpu_id_min") is not None:
        data["gpu_id_min"] = slurm.get("gpu_id_min")
        data["gpu_id_max"] = slurm.get("gpu_id_max")
    if slurm.get("num_gpus") is not None:
        data["num_gpus"] = int(slurm["num_gpus"])
    if slurm.get("use_gpu") is not None:
        data["use_gpu"] = bool(slurm["use_gpu"])
    data.setdefault("summary", resources_summary_for_display(data))
    return data


def resources_summary_for_display(
    data: Dict[str, Any],
) -> str:
    """Compact label for progress cards, e.g. ``Min CPU×4 · MD GPU×1 · Prod GPU×1``."""
    stages = data.get("stages") or []
    if not isinstance(stages, list) or not stages:
        slurm = data.get("slurm") or data
        cpu = slurm.get("cpu_cores") or slurm.get("cpu_cores_max") or 1
        ngpu = slurm.get("num_gpus") or 0
        return f"CPU×{cpu}" + (f" · GPU×{ngpu}" if ngpu else "")

    parts: List[str] = []
    for stage in stages:
        if not isinstance(stage, dict):
            continue
        kind = stage.get("stage_kind") or infer_stage_kind(stage)
        cpu = int(stage.get("cpu_cores") or 1)
        if kind == "minimization":
            parts.append(f"Min CPU×{cpu}")
        elif kind == "production":
            if stage.get("use_gpu"):
                parts.append(f"Prod GPU×{int(stage.get('num_gpus') or 1)}")
            else:
                parts.append(f"Prod CPU×{cpu}")
        elif stage.get("use_gpu"):
            parts.append(f"MD GPU×{int(stage.get('num_gpus') or 1)}")
        else:
            parts.append(f"MD CPU×{cpu}")
    # Deduplicate consecutive identical MD lines
    compact: List[str] = []
    md_seen = False
    for part in parts:
        if part.startswith("MD "):
            if md_seen:
                continue
            md_seen = True
        compact.append(part)
    slurm = data.get("slurm") or {}
    if slurm:
        sc = slurm.get("cpu_cores")
        sg = slurm.get("num_gpus")
        if sc is not None:
            compact.append(f"Slurm CPU×{sc}/GPU×{sg or 0}")
    return " · ".join(compact)


def write_equilibration_resources(
    eq_dir: Path,
    engine: str,
    stages: List[Dict[str, Any]],
    *,
    openmm_platform: Optional[str] = None,
    compute_defaults: Optional[Dict[str, Any]] = None,
    stems: Optional[List[str]] = None,
) -> Path:
    """Persist v2 per-stage resource summary next to run_equilibration.sh."""
    eq_dir = Path(eq_dir)
    engine = engine.lower().strip()
    resolved = resolve_all_stage_resources(
        stages,
        compute_defaults,
        engine=engine,
        stems=stems,
    )
    if openmm_platform:
        plat = openmm_platform.upper()
        for item in resolved:
            if item.get("stage_kind") != "minimization":
                item["platform"] = plat
                item["use_gpu"] = plat not in {"", "AUTO", "CPU", "REFERENCE"}

    payload: Dict[str, Any] = {
        "version": RESOURCES_VERSION,
        "engine": engine,
        "compute_defaults": _merge_defaults(
            DEFAULT_COMPUTE_DEFAULTS, compute_defaults
        ),
        "stages": resolved,
        "slurm": aggregate_slurm_resources(resolved),
    }
    if openmm_platform:
        payload["platform"] = openmm_platform

    path = eq_dir / "equilibration_resources.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def load_equilibration_resources(eq_dir: Path) -> Dict[str, Any]:
    """Load ``equilibration_resources.json`` or return empty dict."""
    path = Path(eq_dir) / "equilibration_resources.json"
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    except (json.JSONDecodeError, OSError):
        return {}
    return data if isinstance(data, dict) else {}


def resolve_compute_resources_from_stages(
    stage_items: Optional[List[Dict[str, Any]]],
    compute_defaults: Optional[Dict[str, Any]] = None,
    *,
    engine: str = "",
) -> Dict[str, Any]:
    """Flat CPU/GPU for legacy single-profile callers (max across stages)."""
    resolved = resolve_all_stage_resources(
        stage_items or [], compute_defaults, engine=engine
    )
    slurm = aggregate_slurm_resources(resolved)
    platform = None
    for item in resolved:
        if item.get("platform"):
            platform = item["platform"]
            break
    return {
        "cpu_cores": int(slurm["cpu_cores"]),
        "use_gpu": bool(slurm["use_gpu"]),
        "gpu_id": int(slurm.get("gpu_id_min") or 0),
        "num_gpus": int(slurm.get("num_gpus") or 0),
        "platform": platform,
        "stages": resolved,
    }


def resolve_compute_resources_from_eq_dir(eq_dir: Path) -> Dict[str, Any]:
    """Load flat + per-stage settings from ``equilibration_resources.json``."""
    eq_dir = Path(eq_dir)
    data = load_equilibration_resources(eq_dir)
    if not data:
        return {
            "cpu_cores": 1,
            "use_gpu": False,
            "gpu_id": 0,
            "num_gpus": 0,
            "platform": None,
            "stages": [],
            "slurm": {},
        }

    if int(data.get("version") or 0) >= RESOURCES_VERSION and data.get("stages"):
        slurm = data.get("slurm") or aggregate_slurm_resources(data["stages"])
        platform = data.get("platform")
        return {
            "cpu_cores": int(slurm.get("cpu_cores") or 1),
            "use_gpu": bool(slurm.get("use_gpu")),
            "gpu_id": int(slurm.get("gpu_id_min") or 0),
            "num_gpus": int(slurm.get("num_gpus") or 0),
            "platform": platform,
            "stages": list(data.get("stages") or []),
            "slurm": slurm,
            "compute_defaults": data.get("compute_defaults"),
        }

    cpu = data.get("cpu_cores_max") or data.get("cpu_cores_min") or data.get("cpu_cores")
    gpu_id = data.get("gpu_id_min") if data.get("gpu_id_min") is not None else data.get("gpu_id")
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
        "stages": list(data.get("stages") or []),
        "slurm": data.get("slurm") or {},
    }


def slurm_resources_from_eq_dir(eq_dir: Path) -> Dict[str, Any]:
    """Return the Slurm aggregate block for cluster submission."""
    flat = resolve_compute_resources_from_eq_dir(eq_dir)
    slurm = flat.get("slurm") or {}
    if slurm:
        return slurm
    return {
        "cpu_cores": flat["cpu_cores"],
        "num_gpus": flat["num_gpus"],
        "use_gpu": flat["use_gpu"],
        "gpu_id_min": flat["gpu_id"],
        "gpu_id_max": flat["gpu_id"],
    }


def _openmm_platform_from_script(script_text: str) -> str:
    if re.search(r'^PLATFORM="[^"]+"', script_text, re.MULTILINE):
        match = re.search(r'^PLATFORM="([^"]+)"', script_text, re.MULTILINE)
        if match and match.group(1):
            return match.group(1)
    return "auto"


def infer_equilibration_resources(eq_dir: Path, engine: str) -> Dict[str, Any]:
    """Return resource summary for job cards and form restore."""
    eq_dir = Path(eq_dir)
    engine = (engine or "").strip().lower()
    data = load_equilibration_resources(eq_dir)
    if data:
        if int(data.get("version") or 0) >= RESOURCES_VERSION:
            data.setdefault("engine", engine)
            return enrich_resources_display(data)
        # Legacy flat file — wrap as v2-like for display
        slurm = {
            "cpu_cores": int(
                data.get("cpu_cores_max") or data.get("cpu_cores_min") or 1
            ),
            "num_gpus": int(data.get("num_gpus") or 0),
            "use_gpu": bool(data.get("use_gpu")),
            "gpu_id_min": data.get("gpu_id_min"),
            "gpu_id_max": data.get("gpu_id_max"),
        }
        wrapped = {
            "version": RESOURCES_VERSION,
            "engine": data.get("engine") or engine,
            "stages": data.get("stages") or [],
            "slurm": slurm,
            "cpu_cores_min": data.get("cpu_cores_min"),
            "cpu_cores_max": data.get("cpu_cores_max"),
            "use_gpu": data.get("use_gpu"),
            "num_gpus": data.get("num_gpus"),
            "gpu_id_min": data.get("gpu_id_min"),
            "gpu_id_max": data.get("gpu_id_max"),
            "platform": data.get("platform"),
        }
        return enrich_resources_display(wrapped)

    if engine == "namd":
        summary_file = eq_dir / "protocol_summary.json"
        if summary_file.is_file():
            try:
                summary = json.loads(summary_file.read_text(encoding="utf-8", errors="replace"))
                stage_map = summary.get("stages") or {}
                items = list(stage_map.values()) if isinstance(stage_map, dict) else list(stage_map)
                if items:
                    resolved = resolve_all_stage_resources(items, engine="namd")
                    payload = {
                        "version": RESOURCES_VERSION,
                        "engine": "namd",
                        "stages": resolved,
                        "slurm": aggregate_slurm_resources(resolved),
                    }
                    return enrich_resources_display(payload)
            except (json.JSONDecodeError, OSError, TypeError, ValueError):
                pass
        script = eq_dir / "run_equilibration.sh"
        if script.is_file():
            try:
                text = script.read_text(encoding="utf-8", errors="replace")
                _, protocols = _parse_namd_protocols_from_script(text)
                if protocols:
                    items = list(protocols.values())
                    resolved = resolve_all_stage_resources(items, engine="namd")
                    payload = {
                        "version": RESOURCES_VERSION,
                        "engine": "namd",
                        "stages": resolved,
                        "slurm": aggregate_slurm_resources(resolved),
                    }
                    return enrich_resources_display(payload)
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
            "version": RESOURCES_VERSION,
            "engine": "openmm",
            "use_gpu": platform.upper() not in {"", "AUTO", "CPU", "REFERENCE"},
            "platform": platform,
            "stages": [],
            "slurm": {"cpu_cores": 1, "num_gpus": 0, "use_gpu": False},
            "summary": f"Platform {platform}",
        }

    return {
        "version": RESOURCES_VERSION,
        "engine": engine,
        "use_gpu": None,
        "platform": None,
        "stages": [],
        "slurm": {},
        "summary": "",
    }
