"""Render scheduler batch wrappers for equilibration jobs."""

from __future__ import annotations

from typing import Dict, List, Optional

from gatewizard.utils.cluster.types import BatchScriptRequest, WORKDIR_STRATEGIES

DEFAULT_SCRATCH_JOB_ID_TEMPLATE = """#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -o {{job_name}}.%j.out
#SBATCH -c {{cpus}}
#SBATCH -t {{time_limit}}
{{extra_sbatch}}
{{mail_block}}

{{purge_block}}
{{module_loads}}

SUBMIT_DIR="$SLURM_SUBMIT_DIR"
workdir="{{scratch_root}}/$SLURM_JOB_ID"
mkdir -p "$workdir" || exit 1
rsync -a --exclude='*.dcd' --exclude='*.xtc' --exclude='*.nc' "$SUBMIT_DIR"/ "$workdir"/ || exit 1
cd "$workdir" || exit 1

# Periodically copy lightweight logs back so Watching can show mid-run progress
# (scratch is often node-local and invisible from the login node until job end).
(
  while true; do
    sleep 60
    rsync -a --include='step*.log' --include='step*.xst' --include='step*.mdinfo' \
      --include='step*.mdout' --include='step*.rst7' --include='step*.rst' \
      --include='step*.coor' --include='step*.gro' --exclude='*' \
      "$workdir"/ "$SUBMIT_DIR"/ 2>/dev/null || true
  done
) &
_GW_SYNC_PID=$!
trap 'kill $_GW_SYNC_PID 2>/dev/null || true' EXIT

{{run_command}}
status=$?

kill $_GW_SYNC_PID 2>/dev/null || true
wait $_GW_SYNC_PID 2>/dev/null || true
rsync -a "$workdir"/ "$SUBMIT_DIR"/ || exit 1
rm -rf "$workdir"
exit $status
"""

DEFAULT_RUN_IN_PLACE_TEMPLATE = """#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -o {{job_name}}.%j.out
#SBATCH -c {{cpus}}
#SBATCH -t {{time_limit}}
{{extra_sbatch}}
{{mail_block}}

{{purge_block}}
{{module_loads}}

cd "$SLURM_SUBMIT_DIR" || exit 1
{{run_command}}
"""

DEFAULT_SCRATCH_NAMED_TEMPLATE = """#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -o {{job_name}}.%j.out
#SBATCH -c {{cpus}}
#SBATCH -t {{time_limit}}
{{extra_sbatch}}
{{mail_block}}

{{purge_block}}
{{module_loads}}

SUBMIT_DIR="$SLURM_SUBMIT_DIR"
workdir="{{scratch_root}}/{{job_folder_name}}"
mkdir -p "$workdir" || exit 1
rsync -a --exclude='*.dcd' --exclude='*.xtc' --exclude='*.nc' "$SUBMIT_DIR"/ "$workdir"/ || exit 1
cd "$workdir" || exit 1

(
  while true; do
    sleep 60
    rsync -a --include='step*.log' --include='step*.xst' --include='step*.mdinfo' \
      --include='step*.mdout' --include='step*.rst7' --include='step*.rst' \
      --include='step*.coor' --include='step*.gro' --exclude='*' \
      "$workdir"/ "$SUBMIT_DIR"/ 2>/dev/null || true
  done
) &
_GW_SYNC_PID=$!
trap 'kill $_GW_SYNC_PID 2>/dev/null || true' EXIT

{{run_command}}
status=$?

kill $_GW_SYNC_PID 2>/dev/null || true
wait $_GW_SYNC_PID 2>/dev/null || true
rsync -a "$workdir"/ "$SUBMIT_DIR"/ || exit 1
rm -rf "$workdir"
exit $status
"""

DEFAULT_TMPDIR_TEMPLATE = """#!/bin/bash
#SBATCH -J {{job_name}}
#SBATCH -o {{job_name}}.%j.out
#SBATCH -c {{cpus}}
#SBATCH -t {{time_limit}}
{{extra_sbatch}}
{{mail_block}}

{{purge_block}}
{{module_loads}}

SUBMIT_DIR="$SLURM_SUBMIT_DIR"
workdir="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/$SLURM_JOB_ID"
mkdir -p "$workdir" || exit 1
rsync -a --exclude='*.dcd' --exclude='*.xtc' --exclude='*.nc' "$SUBMIT_DIR"/ "$workdir"/ || exit 1
cd "$workdir" || exit 1

(
  while true; do
    sleep 60
    rsync -a --include='step*.log' --include='step*.xst' --include='step*.mdinfo' \
      --include='step*.mdout' --include='step*.rst7' --include='step*.rst' \
      --include='step*.coor' --include='step*.gro' --exclude='*' \
      "$workdir"/ "$SUBMIT_DIR"/ 2>/dev/null || true
  done
) &
_GW_SYNC_PID=$!
trap 'kill $_GW_SYNC_PID 2>/dev/null || true' EXIT

{{run_command}}
status=$?

kill $_GW_SYNC_PID 2>/dev/null || true
wait $_GW_SYNC_PID 2>/dev/null || true
rsync -a "$workdir"/ "$SUBMIT_DIR"/ || exit 1
rm -rf "$workdir"
exit $status
"""

_STRATEGY_TEMPLATES = {
    "run_in_place": DEFAULT_RUN_IN_PLACE_TEMPLATE,
    "scratch_job_id": DEFAULT_SCRATCH_JOB_ID_TEMPLATE,
    "scratch_named": DEFAULT_SCRATCH_NAMED_TEMPLATE,
    "tmpdir": DEFAULT_TMPDIR_TEMPLATE,
    "custom_template": DEFAULT_SCRATCH_JOB_ID_TEMPLATE,
}


def _render_simple(template: str, values: Dict[str, str]) -> str:
    out = template
    for key, value in values.items():
        out = out.replace("{{" + key + "}}", value)
    # Drop leftover unknown placeholders to empty
    while "{{" in out and "}}" in out:
        start = out.index("{{")
        end = out.index("}}", start) + 2
        out = out[:start] + out[end:]
    # Clean excessive blank lines from empty optional blocks
    lines = out.splitlines()
    cleaned: List[str] = []
    blank = 0
    for line in lines:
        if line.strip() == "":
            blank += 1
            if blank <= 1:
                cleaned.append(line)
        else:
            blank = 0
            cleaned.append(line)
    return "\n".join(cleaned).rstrip() + "\n"


def default_template_for_strategy(strategy: str) -> str:
    if strategy not in WORKDIR_STRATEGIES:
        strategy = "scratch_job_id"
    return _STRATEGY_TEMPLATES.get(strategy, DEFAULT_SCRATCH_JOB_ID_TEMPLATE)


def render_batch_script(req: BatchScriptRequest) -> str:
    """Render a Slurm batch script from structured options and/or custom template."""
    strategy = req.workdir_strategy if req.workdir_strategy in WORKDIR_STRATEGIES else "scratch_job_id"
    template = req.template or default_template_for_strategy(strategy)

    extra_lines = list(req.extra_sbatch_lines or [])
    if req.gpus and not any("gpu" in ln.lower() for ln in extra_lines):
        extra_lines.append(f"#SBATCH --gpus={int(req.gpus)}")
    if req.partition and not any(
        "--partition" in ln or ln.strip().startswith("#SBATCH -p") for ln in extra_lines
    ):
        extra_lines.insert(0, f"#SBATCH --partition={req.partition}")
    if req.nodelist and not any("--nodelist" in ln or "-w " in ln for ln in extra_lines):
        extra_lines.append(f"#SBATCH --nodelist={req.nodelist.strip()}")
    if req.constraint and not any("--constraint" in ln or "-C " in ln for ln in extra_lines):
        extra_lines.append(f"#SBATCH --constraint={req.constraint.strip()}")
    if req.mem and not any("--mem" in ln for ln in extra_lines):
        extra_lines.append(f"#SBATCH --mem={req.mem}")

    mail_block = ""
    if req.mail_user and req.mail_type and req.mail_type.upper() != "NONE":
        mail_block = (
            f"#SBATCH --mail-user={req.mail_user}\n"
            f"#SBATCH --mail-type={req.mail_type}"
        )

    purge_block = "ml purge" if req.purge_modules else ""
    module_loads = "\n".join(f"module load {m}" for m in req.modules if m)

    values = {
        "job_name": _safe_job_name(req.job_name),
        "cpus": str(max(1, int(req.cpus or 1))),
        "gpus": str(max(0, int(req.gpus or 0))),
        "time_limit": req.time_limit or "24:00:00",
        "partition": req.partition or "",
        "mail_user": req.mail_user or "",
        "mail_type": req.mail_type or "NONE",
        "extra_sbatch": "\n".join(extra_lines),
        "mail_block": mail_block,
        "purge_block": purge_block,
        "purge_modules": purge_block,
        "module_loads": module_loads,
        "modules": module_loads,
        "scratch_root": req.scratch_root or "$SCRATCH_DIR",
        "submit_dir": "$SLURM_SUBMIT_DIR",
        "workdir_strategy": strategy,
        "job_folder_name": req.job_folder_name or _safe_job_name(req.job_name),
        "run_command": req.run_command or "bash run_equilibration_cluster.sh",
    }
    return _render_simple(template, values)


def _safe_job_name(name: str) -> str:
    cleaned = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in (name or "gatewizard"))
    return (cleaned[:64] or "gatewizard").strip("._")
