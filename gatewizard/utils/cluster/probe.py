"""Probe a remote cluster and update equilibration job execution metadata."""

from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from gatewizard.utils.cluster.modules import group_engine_modules, parse_module_avail
from gatewizard.utils.cluster.paths import expand_remote_path, suggest_scratch_root, suggest_submit_root
from gatewizard.utils.cluster.resources import parse_sinfo, parse_sinfo_nodes
from gatewizard.utils.cluster.ssh import run_remote
from gatewizard.utils.cluster.types import ClusterProfile, ProbeResult
from gatewizard.utils.equilibration_job_metadata import JOB_METADATA_FILE


def probe_cluster(
    session_id: str,
    profile: Optional[ClusterProfile] = None,
) -> ProbeResult:
    """Collect modules, paths, and Slurm inventory from an open SSH session."""
    errors: List[str] = []
    now = datetime.now(timezone.utc).isoformat()

    def _cmd(command: str) -> str:
        rc, out, err = run_remote(session_id, command, timeout=90)
        # module avail / sinfo often write to stderr; merge both streams.
        text = "\n".join(part for part in (out or "", err or "") if part.strip())
        if rc != 0 and not text.strip():
            errors.append(f"{command}: {err.strip() or f'exit {rc}'}")
            return ""
        return text

    hostname = _cmd("hostname").strip()
    home = _cmd("echo $HOME").strip()
    data_dir = _cmd("echo ${DATA_DIR:-}").strip()
    scratch_dir = _cmd("echo ${SCRATCH_DIR:-}").strip()
    if not data_dir:
        # common fallbacks
        user = _cmd("echo $USER").strip()
        candidate = f"/data/{user}" if user else ""
        check = _cmd(f"test -d {candidate} && echo {candidate} || true").strip()
        data_dir = check
    if not scratch_dir:
        user = _cmd("echo $USER").strip()
        candidate = f"/scratch/{user}" if user else ""
        check = _cmd(f"test -d {candidate} && echo {candidate} || true").strip()
        scratch_dir = check

    # Prefer full ``module avail`` first. Some Lmod versions return exit 0 from
    # ``module avail -t`` with "No module(s) found" help text, which would
    # short-circuit ``|| module avail`` if -t were tried first.
    raw_modules = _cmd(
        "for f in "
        "/etc/profile "
        "/etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /etc/profile.d/modules.sh "
        "/opt/ohpc/admin/lmod/lmod/init/bash "
        "/usr/share/lmod/lmod/init/bash "
        "$HOME/.bash_profile $HOME/.bashrc; do "
        "[ -r \"$f\" ] && . \"$f\" >/dev/null 2>&1 || true; "
        "done; "
        "command -v module >/dev/null 2>&1 || type module >/dev/null 2>&1 || "
        "echo 'GW_PROBE: module command not found after sourcing profiles'; "
        "out=$(module avail 2>&1 || true); "
        "if ! echo \"$out\" | grep -qE '/|md/|cuda/'; then "
        "  out=$(module --default avail 2>&1 || module avail -t 2>&1 || ml av 2>&1 || true); "
        "fi; "
        "printf '%s\\n' \"$out\""
    )
    modules = parse_module_avail(raw_modules)
    if not modules:
        raw_lower = (raw_modules or "").lower()
        if "gw_probe: module command not found" in raw_lower:
            errors.append(
                "module command not found on the login node after sourcing profiles. "
                "Ask your admin how Environment Modules / Lmod is initialized, "
                "or type module paths manually in the dialog."
            )
        elif not raw_modules.strip():
            errors.append(
                "module avail returned no output. Load the modules environment "
                "(e.g. source /etc/profile or module init) and retry Connect & probe."
            )
        elif any(
            token in raw_lower
            for token in ("no module", "not found", "command not found", "unknown command")
        ):
            errors.append(
                "module avail failed or found nothing. Check that Environment Modules / Lmod "
                "is available on the login node."
            )
        else:
            errors.append(
                "Could not parse any MD software modules from module avail. "
                "Raw output may be help text or an unsupported format — set modules manually."
            )
    hints = profile.module_hints if profile else None
    engine_modules = group_engine_modules(modules, hints=hints)

    raw_sinfo = _cmd(
        "command -v sinfo >/dev/null 2>&1 || echo 'GW_PROBE: sinfo not found'; "
        "sinfo -o '%P %a %D %c %G %m %l' 2>&1 || sinfo 2>&1 || true"
    )
    partitions = parse_sinfo(raw_sinfo)
    # Deduplicate partition names (sinfo may repeat one name per node feature set).
    if partitions:
        seen_names: set = set()
        deduped = []
        for part in partitions:
            key = (part.name or "").rstrip("*")
            if key in seen_names:
                for i, existing in enumerate(deduped):
                    if (existing.name or "").rstrip("*") == key:
                        if (part.max_gpus or 0) > (existing.max_gpus or 0):
                            deduped[i] = part
                        break
                continue
            seen_names.add(key)
            part.name = key
            deduped.append(part)
        partitions = deduped

    raw_sinfo_nodes = _cmd(
        "command -v sinfo >/dev/null 2>&1 || true; "
        "sinfo -N -h -o '%N|%P|%T|%c|%G|%f' 2>&1 || "
        "sinfo -N -o '%N %P %T %c %G %f' 2>&1 || true"
    )
    nodes = parse_sinfo_nodes(raw_sinfo_nodes)

    if "gw_probe: sinfo not found" in (raw_sinfo or "").lower():
        errors.append(
            "sinfo not found on this host (Slurm client missing on the login node?). "
            "Enter the partition name manually."
        )
    elif not partitions and raw_sinfo.strip():
        errors.append(
            "Could not parse Slurm partitions from sinfo. Enter the partition name manually."
        )
    if partitions and not nodes and raw_sinfo_nodes.strip():
        errors.append(
            "Could not parse node names from sinfo -N. You can still type a nodelist manually."
        )

    return ProbeResult(
        hostname=hostname,
        home=home,
        data_dir=data_dir,
        scratch_dir=scratch_dir,
        modules=modules,
        engine_modules=engine_modules,
        partitions=partitions,
        nodes=nodes,
        raw_module_avail=raw_modules,
        raw_sinfo=raw_sinfo,
        raw_sinfo_nodes=raw_sinfo_nodes,
        probed_at=now,
        errors=errors,
    )


def apply_probe_defaults(profile: ClusterProfile, probe: ProbeResult) -> ClusterProfile:
    """Fill empty / templated submit/scratch roots from probe results.

    Never replace an explicit absolute path the user configured in Settings.
    """
    username = profile.username
    submit = (profile.submit_root or "").strip()
    if not submit or "$" in submit:
        if "$" in submit:
            profile.submit_root = expand_remote_path(
                submit,
                username=username,
                home=probe.home,
                data_dir=probe.data_dir,
                scratch_dir=probe.scratch_dir,
            )
        if not profile.submit_root or "$" in (profile.submit_root or ""):
            suggested = suggest_submit_root(
                data_dir=probe.data_dir, home=probe.home, username=username
            )
            if suggested:
                profile.submit_root = suggested
    # else: keep user's absolute submit_root unchanged

    scratch = (profile.scratch_root or "").strip()
    if not scratch or scratch in {"$SCRATCH_DIR", "${SCRATCH_DIR}"} or "$" in scratch:
        if scratch and "$" in scratch and scratch not in {"$SCRATCH_DIR", "${SCRATCH_DIR}"}:
            profile.scratch_root = expand_remote_path(
                scratch,
                username=username,
                home=probe.home,
                data_dir=probe.data_dir,
                scratch_dir=probe.scratch_dir,
            )
        if (
            not profile.scratch_root
            or profile.scratch_root in {"$SCRATCH_DIR", "${SCRATCH_DIR}"}
            or "$" in (profile.scratch_root or "")
        ):
            profile.scratch_root = suggest_scratch_root(
                scratch_dir=probe.scratch_dir or "$SCRATCH_DIR", username=username
            )
    # else: keep user's absolute scratch_root unchanged

    profile.last_probe = probe.to_dict()
    return profile


def read_job_metadata(eq_dir: Path) -> Dict[str, Any]:
    path = Path(eq_dir) / JOB_METADATA_FILE
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def write_execution_metadata(eq_dir: Path, execution: Dict[str, Any]) -> Path:
    """Merge ``execution`` into equilibration_job.json (create minimal file if needed)."""
    eq_dir = Path(eq_dir)
    eq_dir.mkdir(parents=True, exist_ok=True)
    path = eq_dir / JOB_METADATA_FILE
    payload = read_job_metadata(eq_dir)
    payload["execution"] = execution
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def update_execution_fields(eq_dir: Path, **fields: Any) -> Dict[str, Any]:
    payload = read_job_metadata(eq_dir)
    execution = dict(payload.get("execution") or {})
    execution.update({k: v for k, v in fields.items() if v is not None})
    payload["execution"] = execution
    path = Path(eq_dir) / JOB_METADATA_FILE
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return execution


def read_batch_script_resources(eq_dir: Path) -> Dict[str, Any]:
    """Parse ``#SBATCH -c`` / ``--gpus`` from the local batch script."""
    eq_dir = Path(eq_dir)
    for name in ("run_equilibration.slurm", "run_equilibration.sbatch"):
        path = eq_dir / name
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        cpus = 0
        gpus = 0
        for line in text.splitlines():
            s = line.strip()
            if not s.startswith("#SBATCH"):
                continue
            m = re.search(r"(?:-c|--cpus-per-task=)\s*(\d+)", s)
            if m:
                cpus = int(m.group(1))
            m = re.search(r"(?:--gpus(?:-per-node)?|--gres=gpu:)=?\s*(\d+)", s, re.I)
            if m:
                gpus = int(m.group(1))
            else:
                m = re.search(r"--gres=gpu(?::[^\s:]+)?:(\d+)", s, re.I)
                if m:
                    gpus = int(m.group(1))
                elif re.search(r"--gres=gpu\b", s, re.I) and gpus == 0:
                    gpus = 1
        return {"cpus": cpus, "gpus": gpus, "batch_script": name}
    return {}


def enrich_execution_resources(
    eq_dir: Path,
    execution: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Fill missing ``allocated_cpus`` / ``resources`` from the batch script."""
    eq_dir = Path(eq_dir)
    if execution is None:
        execution = dict(read_job_metadata(eq_dir).get("execution") or {})
    else:
        execution = dict(execution)
    batch = read_batch_script_resources(eq_dir)
    res = dict(execution.get("resources") or {}) if isinstance(execution.get("resources"), dict) else {}
    if batch.get("cpus") and not res.get("cpus"):
        res["cpus"] = int(batch["cpus"])
    if batch.get("gpus") and not res.get("gpus"):
        res["gpus"] = int(batch["gpus"])
    if res:
        execution["resources"] = res
    if not execution.get("allocated_cpus") and res.get("cpus"):
        execution["allocated_cpus"] = int(res["cpus"])
    return execution
