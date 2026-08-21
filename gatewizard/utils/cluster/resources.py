"""Parse Slurm inventory (``sinfo``) into a common partition shape."""

from __future__ import annotations

import re
import shlex
from typing import Any, Dict, List, Optional

from gatewizard.utils.cluster.types import NodeInfo, PartitionInfo


def parse_sinfo(text: str) -> List[PartitionInfo]:
    """Parse ``sinfo -o '%P %a %D %c %G %m %l'`` or default ``sinfo`` tables."""
    lines = [ln.rstrip() for ln in text.splitlines() if ln.strip()]
    if not lines:
        return []

    # Prefer machine-readable custom format: PARTITION AVAIL NODES CPUS GRES MEMORY TIMELIMIT
    custom = _parse_custom_sinfo(lines)
    if custom:
        return custom
    return _parse_default_sinfo(lines)


def parse_sinfo_nodes(text: str) -> List[NodeInfo]:
    """Parse ``sinfo -N -h -o '%N|%P|%T|%c|%G|%f'`` (pipe-separated) or space-separated."""
    nodes: List[NodeInfo] = []
    seen_keys = set()
    for line in text.splitlines():
        raw = line.strip()
        if not raw or raw.lower().startswith(("nodelist", "node", "sinfo", "bash")):
            continue
        parts = [p.strip() for p in raw.split("|")] if "|" in raw else raw.split()
        if len(parts) < 2:
            continue
        name = parts[0]
        partition = parts[1].rstrip("*") if len(parts) > 1 else ""
        key = (name, partition)
        if key in seen_keys:
            continue
        if ":" in name or name.lower() in {"bash", "error", "command"}:
            continue
        state = parts[2] if len(parts) > 2 else ""
        cpus = _to_int(parts[3]) if len(parts) > 3 else 0
        gres = parts[4] if len(parts) > 4 else ""
        features = parts[5] if len(parts) > 5 else ""
        if gres in {"(null)", "N/A"}:
            gres = ""
        if features in {"(null)", "N/A"}:
            features = ""
        seen_keys.add(key)
        nodes.append(
            NodeInfo(
                name=name,
                partition=partition,
                state=state,
                cpus=cpus,
                gpus=_gpus_from_gres(gres),
                gres=gres,
                features=features,
            )
        )
    return nodes


def prefer_nodes(
    nodes: List[NodeInfo],
    *,
    want_gpu: bool,
    partition: str = "",
) -> List[NodeInfo]:
    """Sort / filter nodes for the dialog: prefer GPU + idle when requesting GPUs."""
    part = (partition or "").rstrip("*").lower()
    filtered = []
    for node in nodes:
        state = (node.state or "").lower()
        if any(tok in state for tok in ("down", "drain", "fail", "unk", "not")):
            continue
        if part and (node.partition or "").rstrip("*").lower() != part:
            continue
        filtered.append(node)

    def score(node: NodeInfo) -> tuple:
        state = (node.state or "").lower()
        idle = 0 if "idle" in state else (1 if "mix" in state else 2)
        has_gpu = (node.gpus or 0) > 0 or "gpu" in (node.gres or "").lower()
        if want_gpu:
            return (0 if has_gpu else 1, idle, -(node.gpus or 0), node.name)
        return (0 if not has_gpu else 1, idle, node.name)

    return sorted(filtered, key=score)


def _parse_custom_sinfo(lines: List[str]) -> List[PartitionInfo]:
    header = lines[0].lower()
    if "partition" not in header and "%p" not in header:
        # Try without header: space-separated fields from -o
        pass
    results: List[PartitionInfo] = []
    start = 1 if any(h in header for h in ("partition", "avail", "nodes")) else 0
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 2:
            continue
        name = parts[0].rstrip("*")
        # Skip shell noise / non-partition tokens (e.g. "bash:" from stderr merge)
        if not name or ":" in name or name.lower() in {"bash", "sinfo", "error"}:
            continue
        avail = parts[1] if len(parts) > 1 else ""
        if avail.lower() not in {"up", "down", "drain", "inq", "mix*", "mix", "alloc", "idle", "unk"} and not avail.lower().startswith(
            ("up", "down", "drain", "mix", "alloc", "idle")
        ):
            # Custom format expects AVAIL as second field; skip junk rows
            if avail.lower() in {"cannot", "no", "command"}:
                continue
        nodes = _to_int(parts[2]) if len(parts) > 2 else 0
        cpus = _to_int(parts[3]) if len(parts) > 3 else 0
        gres = parts[4] if len(parts) > 4 else ""
        mem = parts[5] if len(parts) > 5 else ""
        tmax = parts[6] if len(parts) > 6 else ""
        results.append(
            PartitionInfo(
                name=name,
                nodes=nodes,
                cpus_per_node=cpus,
                max_gpus=_gpus_from_gres(gres),
                max_mem=mem if mem not in {"(null)", "N/A"} else "",
                max_time=tmax if tmax not in {"infinite", "UNLIMITED"} else tmax,
                avail=avail,
            )
        )
    return results


def _parse_default_sinfo(lines: List[str]) -> List[PartitionInfo]:
    results: List[PartitionInfo] = []
    for line in lines:
        if line.lower().startswith("partition"):
            continue
        parts = line.split()
        if not parts:
            continue
        name = parts[0].rstrip("*")
        avail = parts[1] if len(parts) > 1 else ""
        nodes = _to_int(parts[3]) if len(parts) > 3 else 0
        results.append(PartitionInfo(name=name, avail=avail, nodes=nodes))
    return results


def _gpus_from_gres(gres: str) -> int:
    """Count GPUs in Slurm GRES (sum multi-type lines like ``gpu:2080ti:2,gpu:3090:1``)."""
    if not gres or gres in {"(null)", "N/A"}:
        return 0
    counts = [int(x) for x in re.findall(r"gpu(?::\w+)?:(\d+)", gres, re.I)]
    if counts:
        return sum(counts)
    if "gpu" in gres.lower():
        return 1
    return 0


def parse_gpu_types_from_gres(gres: str) -> List[Dict[str, Any]]:
    """Named GPU types from a Slurm GRES string.

    Examples::

        gpu:2080ti:2,gpu:3090:1 → [{"type": "2080ti", "count": 2}, {"type": "3090", "count": 1}]
        gpu:2                   → []  (untyped; UI shows only Any)
        gpu:l4:2                → [{"type": "l4", "count": 2}]
    """
    if not gres or gres in {"(null)", "N/A", "n/a"}:
        return []
    out: List[Dict[str, Any]] = []
    for chunk in gres.split(","):
        chunk = chunk.strip()
        if not chunk.lower().startswith("gpu"):
            continue
        bits = [b.strip() for b in chunk.split(":") if b.strip()]
        # Named types need gpu:TYPE:N (len 3). Untyped gpu:N has len 2 — skip.
        if len(bits) < 3:
            continue
        type_name = bits[1]
        count_s = bits[2]
        # TYPE may be numeric (3090) — that is still a named GRES type.
        if not type_name or not re.fullmatch(r"[A-Za-z0-9_+\-.]+", type_name):
            continue
        try:
            count = int(count_s)
        except ValueError:
            continue
        if count <= 0:
            continue
        out.append({"type": type_name, "count": count})
    return out


def gpu_types_from_nodes(
    nodes: List[NodeInfo],
    *,
    partition: str = "",
    nodelist: str = "",
) -> List[Dict[str, Any]]:
    """Union of named GPU types across nodes (optionally filtered).

    When ``nodelist`` is set, only that node’s types are returned.
    Counts are the max seen for each type across matching nodes.
    """
    part = (partition or "").rstrip("*").lower()
    want_node = (nodelist or "").split(",")[0].strip().lower()
    merged: Dict[str, int] = {}
    for node in nodes or []:
        if want_node and (node.name or "").lower() != want_node:
            continue
        if part and (node.partition or "").rstrip("*").lower() != part:
            continue
        for item in parse_gpu_types_from_gres(node.gres or ""):
            t = str(item.get("type") or "")
            c = int(item.get("count") or 0)
            if not t or c <= 0:
                continue
            merged[t] = max(merged.get(t, 0), c)
    return [{"type": t, "count": merged[t]} for t in sorted(merged.keys())]


def normalize_gpu_type(value: Optional[str]) -> str:
    """Sanitize a user/Slurm GPU type token (empty → any)."""
    t = (value or "").strip()
    if not t or t.lower() in {"any", "auto", "none", "*"}:
        return ""
    # Slurm GRES types are typically [A-Za-z0-9_+.-]
    if not re.fullmatch(r"[A-Za-z0-9_+\-.]+", t):
        return ""
    return t


def canonicalize_slurm_state(state: str) -> str:
    """Normalize Slurm/sacct states (e.g. ``CANCELLED by 1002`` → ``CANCELLED``)."""
    s = (state or "").strip()
    if not s:
        return ""
    u = s.upper()
    # Longest-first so NODE_FAIL wins over FAIL*, OUT_OF_MEMORY over MEMORY, etc.
    bases = (
        "OUT_OF_MEMORY",
        "NODE_FAIL",
        "BOOT_FAIL",
        "CONFIGURING",
        "COMPLETING",
        "CANCELLED",
        "COMPLETED",
        "PREEMPTED",
        "SUSPENDED",
        "REQUEUED",
        "TIMEOUT",
        "RUNNING",
        "PENDING",
        "FAILED",
        "COMPLETE",
    )
    for base in bases:
        if u == base or u.startswith(base + " ") or u.startswith(base + "+"):
            return "COMPLETED" if base == "COMPLETE" else base
    token = u.split()[0].split("+")[0]
    return token


def _to_int(value: str) -> int:
    try:
        return int(re.sub(r"[^0-9]", "", value) or "0")
    except ValueError:
        return 0


def prefer_partitions(
    partitions: List[PartitionInfo], *, want_gpu: bool
) -> List[PartitionInfo]:
    """Sort partitions: prefer GPU-capable when ``want_gpu``, else default/normal."""

    def score(part: PartitionInfo) -> tuple:
        name = (part.name or "").lower()
        has_gpu = (part.max_gpus or 0) > 0 or "gpu" in name
        is_defaultish = name in {"normal", "batch", "default", "compute", "main"} or name.endswith(
            "default"
        )
        if want_gpu:
            return (0 if has_gpu else 1, 0 if is_defaultish else 1, name)
        return (0 if not has_gpu else 1, 0 if is_defaultish else 1, name)

    return sorted(partitions, key=score)


def parse_squeue_me(text: str) -> List[dict]:
    """Parse ``squeue`` rows into job dicts.

    Prefers pipe-separated ``%i|%T|%j|%P|%R|%N|%M|%C`` (robust when ``%R`` is
    empty). Falls back to space-separated for older fixtures.
    """
    jobs: List[dict] = []
    lines = [ln for ln in text.splitlines() if ln.strip()]
    start = 0
    if lines and lines[0].lower().startswith("jobid"):
        start = 1
    for line in lines[start:]:
        raw = line.strip()
        if "|" in raw:
            parts = [p.strip() for p in raw.split("|")]
            if len(parts) < 2:
                continue
            reason = parts[4] if len(parts) > 4 else ""
            node_list = parts[5] if len(parts) > 5 else ""
            # Prefer %N; if empty and %R is not a pending reason, use %R.
            if not node_list and reason and not reason.startswith("("):
                node_list = reason
            cpus_raw = parts[7] if len(parts) > 7 else ""
            try:
                cpus = int(cpus_raw) if cpus_raw else 0
            except ValueError:
                cpus = 0
            jobs.append(
                {
                    "job_id": parts[0],
                    "state": parts[1] if len(parts) > 1 else "",
                    "name": parts[2] if len(parts) > 2 else "",
                    "partition": parts[3] if len(parts) > 3 else "",
                    "reason": reason,
                    "node_list": node_list,
                    "elapsed": parts[6] if len(parts) > 6 else "",
                    "cpus": cpus,
                    "raw": line,
                }
            )
            continue
        parts = line.split(None, 6)
        if len(parts) < 2:
            continue
        jobs.append(
            {
                "job_id": parts[0],
                "state": parts[1] if len(parts) > 1 else "",
                "name": parts[2] if len(parts) > 2 else "",
                "partition": parts[3] if len(parts) > 3 else "",
                "reason": parts[4] if len(parts) > 4 else "",
                "node_list": parts[5] if len(parts) > 5 else "",
                "elapsed": parts[6] if len(parts) > 6 else "",
                "cpus": 0,
                "raw": line,
            }
        )
    return jobs


def summarize_node_gpu_label(gres: str) -> str:
    """Human GPU type(s) from ``sinfo`` GRES (no counts).

    Example: ``gpu:l4:2,gpu:2080ti:1`` → ``L4 / 2080 Ti``.
    """
    if not gres or gres in {"(null)", "N/A", "n/a"}:
        return ""
    labels: List[str] = []
    for chunk in gres.split(","):
        chunk = chunk.strip()
        if not chunk.lower().startswith("gpu"):
            continue
        # gpu[:type[:count]]
        bits = chunk.split(":")
        if len(bits) >= 2 and bits[1] and not bits[1].isdigit():
            pretty = bits[1].replace("_", " ")
            # 2080ti → 2080 Ti, l4 → L4, a100 → A100
            if re.fullmatch(r"\d+ti", pretty, re.I):
                pretty = re.sub(r"ti$", " Ti", pretty, flags=re.I)
            elif re.fullmatch(r"\d+super", pretty, re.I):
                pretty = re.sub(r"super$", " Super", pretty, flags=re.I)
            else:
                pretty = pretty.upper() if len(pretty) <= 4 else pretty
            labels.append(pretty)
        else:
            labels.append("GPU")
    seen = set()
    out: List[str] = []
    for lab in labels:
        key = lab.lower()
        if key in seen:
            continue
        seen.add(key)
        out.append(lab)
    return " / ".join(out)


def query_node_gres(session_id: str, node: str) -> str:
    """Return GRES string for a compute node via ``sinfo -N``."""
    from gatewizard.utils.cluster.ssh import run_remote

    node = (node or "").split(",")[0].strip()
    if not node:
        return ""
    cmd = f"sinfo -N -n {shlex.quote(node)} -h -o '%G' 2>/dev/null | head -1 || true"
    try:
        _rc, out, _err = run_remote(session_id, cmd, timeout=30)
    except Exception:
        return ""
    gres = (out or "").strip().splitlines()
    val = gres[0].strip() if gres else ""
    if val in {"(null)", "N/A", "n/a"}:
        return ""
    return val


def parse_sbatch_output(text: str) -> str:
    """Extract job id from ``Submitted batch job 238``."""
    match = re.search(r"Submitted batch job\s+(\d+)", text)
    if match:
        return match.group(1)
    match = re.search(r"\b(\d+)\b", text.strip())
    return match.group(1) if match else ""
