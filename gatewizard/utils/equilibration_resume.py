"""
Equilibration resume helpers — stage-level continue (skip completed stages).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

from gatewizard.utils import gromacs_analysis, namd_analysis, openmm_analysis
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)


@dataclass
class EquilibrationResumePoint:
    can_resume: bool
    stage_index: int = -1
    stage_name: str = ""
    stage_stem: str = ""
    mode: str = "stage"
    reason: str = ""
    completed_stages: int = 0
    total_stages: int = 0


OPENMM_RESUME_SHELL = """
# Stage-level resume: skip stages that already finished (not mid-stage checkpoints).
# RESUME=1 bash run_equilibration.sh
RESUME="${RESUME:-0}"

_gw_openmm_stage_done() {
  local stem="$1"
  [ -f "${stem}.rst" ] && grep -qE '100(\\.0)?%' "${stem}.log" 2>/dev/null
}
""".strip()

GROMACS_RESUME_SHELL = """
# Stage-level resume: skip stages that already finished (not mid-stage checkpoints).
# RESUME=1 bash run_equilibration.sh
RESUME="${RESUME:-0}"

_gw_gromacs_stage_done() {
  local prefix="$1"
  [ -f "${prefix}.gro" ] && grep -q "Finished mdrun" "${prefix}.log" 2>/dev/null
}
""".strip()

NAMD_RESUME_SHELL = """
# Stage-level resume: skip stages that already finished (not mid-stage checkpoints).
# RESUME=1 bash run_equilibration.sh
RESUME="${RESUME:-0}"

_gw_namd_stage_done() {
  local stem="$1"
  [ -f "${stem}.coor" ] && [ -f "${stem}.log" ] && ! grep -qi "Error in Stage" "${stem}.log" 2>/dev/null
}
""".strip()


OPENMM_STAGE_ORDER: List[tuple[str, str]] = [
    ("step1_equilibration", "Equilibration 1"),
    ("step2_equilibration", "Equilibration 2"),
    ("step3_equilibration", "Equilibration 3"),
    ("step4_equilibration", "Equilibration 4"),
    ("step5_equilibration", "Equilibration 5"),
    ("step6_equilibration", "Equilibration 6"),
    ("step7_production", "Production"),
]


def _openmm_stage_complete(eq_dir: Path, stem: str) -> bool:
    """Match get_equilibration_progress / GUI: completion from the stage log."""
    log = eq_dir / f"{stem}.log"
    if not log.is_file():
        return False
    try:
        timing = openmm_analysis.parse_openmm_log(log, eq_dir / f"{stem}.inp")
        return timing.completed
    except Exception:
        return False


def _gromacs_stage_complete(eq_dir: Path, prefix: str) -> bool:
    gro = eq_dir / f"{prefix}.gro"
    log = eq_dir / f"{prefix}.log"
    if not gro.is_file() or not log.is_file():
        return False
    try:
        return "Finished mdrun" in log.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return False


def _namd_stage_complete(eq_dir: Path, stem: str) -> bool:
    coor = eq_dir / f"{stem}.coor"
    log = eq_dir / f"{stem}.log"
    if not coor.is_file() or not log.is_file():
        return False
    try:
        text = log.read_text(encoding="utf-8", errors="replace")
        if "Error in Stage" in text:
            return False
        timing = namd_analysis.parse_namd_log(log)
        if timing.total_steps > 0 and timing.steps_completed >= timing.total_steps:
            return True
        return "WRITING" in text or "End of program" in text
    except Exception:
        return False


def _stage_stems_on_disk(eq_dir: Path, engine: str) -> List[tuple[str, str, str]]:
    """Return (key, display_name, filesystem_stem) in protocol order."""
    engine = engine.lower().strip()
    entries: List[tuple[str, str, str]] = []
    if engine == "openmm":
        for stem, label in OPENMM_STAGE_ORDER:
            if (eq_dir / f"{stem}.inp").is_file():
                entries.append((stem, label, stem))
    elif engine == "gromacs":
        for mdp in sorted(eq_dir.glob("step*.mdp")):
            stem = mdp.stem
            entries.append((stem, stem.replace("_", " ").title(), stem))
    else:
        for conf in sorted(eq_dir.glob("step*.conf")):
            if "_restraints" in conf.name:
                continue
            stem = conf.stem
            if stem.endswith("_equilibration") or stem == "step7_production":
                entries.append((stem, stem.replace("_", " ").title(), stem))
    return entries


def _is_stage_complete(eq_dir: Path, engine: str, stem: str) -> bool:
    engine = engine.lower().strip()
    if engine == "openmm":
        return _openmm_stage_complete(eq_dir, stem)
    if engine == "gromacs":
        return _gromacs_stage_complete(eq_dir, stem)
    return _namd_stage_complete(eq_dir, stem)


def _stage_has_partial_output(eq_dir: Path, engine: str, stem: str) -> bool:
    """True when a stage has run output but is not fully complete."""
    if _is_stage_complete(eq_dir, engine, stem):
        return False
    if engine == "openmm":
        log = eq_dir / f"{stem}.log"
        return log.is_file() and log.stat().st_size > 0
    if engine == "gromacs":
        log = eq_dir / f"{stem}.log"
        return log.is_file() and log.stat().st_size > 0
    log = eq_dir / f"{stem}.log"
    return log.is_file() and log.stat().st_size > 0


def protocol_was_interrupted(eq_dir: Path, engine: str) -> bool:
    """True when a protocol was started but has not finished all stages."""
    eq_dir = Path(eq_dir)
    if (eq_dir / "equilibration_start_time.txt").is_file():
        return True
    bg_log = eq_dir / "equilibration_background.log"
    if bg_log.is_file() and bg_log.stat().st_size > 0:
        return True
    for _, _, stem in _stage_stems_on_disk(eq_dir, engine):
        if _stage_has_partial_output(eq_dir, engine, stem):
            return True
    return False


def equilibration_script_supports_resume(script_path: Path) -> bool:
    """True when run_equilibration.sh includes stage-level RESUME helpers."""
    try:
        text = script_path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return False
    return "_gw_" in text and 'RESUME="${RESUME' in text


def _parse_script_var(text: str, name: str) -> Optional[str]:
    match = re.search(rf'^{re.escape(name)}="([^"]*)"', text, re.MULTILINE)
    return match.group(1) if match else None


def _parse_bash_default_var(text: str, name: str) -> Optional[str]:
    """Parse ``NAME="${NAME:-default}"`` and return the default (may be empty)."""
    match = re.search(
        rf'^{re.escape(name)}="\$\{{{re.escape(name)}:-([^}}]*)\}}"',
        text,
        re.MULTILINE,
    )
    if match:
        return match.group(1)
    # Plain assignment fallback: NAME="value"
    return _parse_script_var(text, name)


def _parse_gmxrc_path(text: str) -> Optional[str]:
    match = re.search(r'^source "([^"]+)"', text, re.MULTILINE)
    return match.group(1) if match else None


def _parse_namd_protocols_from_script(text: str) -> Tuple[Optional[str], dict]:
    """Rebuild a minimal NAMD protocols dict from an existing run script."""
    namd_exe = _parse_script_var(text, "NAMD")
    protocols: dict = {}
    pattern = re.compile(
        r"\$NAMD(\s+\+p\d+)?(\s+\+devices\s+[\d,]+)?\s+(\S+\.conf)\s*>\s*(\S+\.log)",
        re.MULTILINE,
    )
    for index, match in enumerate(pattern.finditer(text)):
        command = match.group(0)
        cpu_match = re.search(r"\+p(\d+)", command)
        gpu_match = re.search(r"\+devices\s+([\d,]+)", command)
        conf_stem = Path(match.group(3)).stem
        if conf_stem == "step7_production":
            stage_key = "step7_production"
        else:
            stage_key = conf_stem.removesuffix("_equilibration")
        protocols[stage_key] = {
            "name": stage_key.replace("_", " ").title(),
            "steps": "N/A",
            "timestep": "N/A",
            "use_gpu": gpu_match is not None,
            "cpu_cores": int(cpu_match.group(1)) if cpu_match else 1,
            "gpu_id": int(gpu_match.group(1).split(",")[0]) if gpu_match else 0,
            "num_gpus": len(gpu_match.group(1).split(",")) if gpu_match else 1,
        }
    return namd_exe, protocols


def refresh_equilibration_run_script(eq_dir: Path, engine: str) -> bool:
    """
    Rewrite ``run_equilibration.sh`` from on-disk stage files (configs unchanged).

    Used when an older script predates stage-level RESUME support.
    """
    eq_dir = Path(eq_dir)
    engine = (engine or "").strip().lower()
    script = eq_dir / "run_equilibration.sh"
    if not script.is_file():
        return False

    try:
        if engine == "openmm":
            from gatewizard.tools.equilibration import OpenMMEquilibrationManager
            from gatewizard.utils.equilibration_resources import (
                resolve_compute_resources_from_eq_dir,
            )

            text = script.read_text(encoding="utf-8", errors="replace")
            stage_config_names = [p.stem for p in sorted(eq_dir.glob("step*.inp"))]
            if not stage_config_names:
                return False
            compute = resolve_compute_resources_from_eq_dir(eq_dir)
            # Prefer PLATFORM default already written into the script (backend may patch it).
            script_platform = _parse_bash_default_var(text, "PLATFORM")
            manager = OpenMMEquilibrationManager(eq_dir)
            manager.generate_run_script(
                stage_config_names=stage_config_names,
                openmm_dir=eq_dir,
                prmtop_name=_parse_script_var(text, "PRMTOP") or "system.prmtop",
                inpcrd_name=_parse_script_var(text, "INPCRD") or "system.inpcrd",
                bilayer_pdb_name=_parse_script_var(text, "BILAYER_PDB"),
                cpu_cores=compute["cpu_cores"],
                use_gpu=compute["use_gpu"],
                gpu_id=compute["gpu_id"],
                num_gpus=compute["num_gpus"] or 1,
                platform=script_platform or compute.get("platform"),
            )
            return True

        if engine == "gromacs":
            from gatewizard.tools.equilibration import GROMACSEquilibrationManager
            from gatewizard.utils.equilibration_resources import (
                resolve_compute_resources_from_eq_dir,
            )

            text = script.read_text(encoding="utf-8", errors="replace")
            n_stages = len(list(eq_dir.glob("step[1-9]_equilibration.mdp")))
            if n_stages == 0:
                return False
            compute = resolve_compute_resources_from_eq_dir(eq_dir)
            manager = GROMACSEquilibrationManager(eq_dir)
            ndx_name = _parse_script_var(text, "NDX")
            manager.generate_run_script(
                gromacs_dir=eq_dir,
                gro_name=_parse_script_var(text, "GRO") or "step5_input.gro",
                top_name=_parse_script_var(text, "TOP") or "topol_posres.top",
                ndx_name=ndx_name,
                n_stages=n_stages,
                gmx_executable=_parse_script_var(text, "GMX") or "gmx",
                gmxrc_path=_parse_gmxrc_path(text),
                cpu_cores=compute["cpu_cores"],
                use_gpu=compute["use_gpu"],
                gpu_id=compute["gpu_id"],
                num_gpus=compute["num_gpus"] or 1,
            )
            return True

        if engine == "namd":
            from gatewizard.tools.equilibration import NAMDEquilibrationManager

            summary_file = eq_dir / "protocol_summary.json"
            namd_exe: Optional[str] = None
            protocols: dict = {}
            if summary_file.is_file():
                try:
                    summary = json.loads(
                        summary_file.read_text(encoding="utf-8", errors="replace")
                    )
                    protocols = summary.get("stages") or {}
                    namd_exe = summary.get("namd_executable")
                except (json.JSONDecodeError, OSError):
                    protocols = {}
            if not protocols:
                text = script.read_text(encoding="utf-8", errors="replace")
                namd_exe, protocols = _parse_namd_protocols_from_script(text)
            if not protocols:
                return False
            manager = NAMDEquilibrationManager(eq_dir, namd_executable=namd_exe or "namd3")
            script.write_text(manager.generate_run_script(protocols, namd_exe))
            script.chmod(0o755)
            return True
    except Exception as exc:
        logger.warning("Failed to refresh equilibration run script in %s: %s", eq_dir, exc)
        return False

    return False


def get_equilibration_resume_point(workdir: Path, engine: str) -> EquilibrationResumePoint:
    """
    Detect whether a stage-level continue is possible.

    Resume is allowed when at least one stage finished, or when a run was
    interrupted (including kill/crash mid-stage).
    """
    eq_dir = Path(workdir)
    engine = (engine or "").strip().lower()
    if not (eq_dir / "run_equilibration.sh").is_file():
        return EquilibrationResumePoint(
            can_resume=False, reason="run_equilibration.sh not found"
        )

    stems = _stage_stems_on_disk(eq_dir, engine)
    if not stems:
        return EquilibrationResumePoint(
            can_resume=False, reason="No stage input files found"
        )

    completed = 0
    first_incomplete_idx = -1
    first_incomplete_name = ""
    first_incomplete_stem = ""

    for idx, (_key, display_name, stem) in enumerate(stems):
        if _is_stage_complete(eq_dir, engine, stem):
            completed += 1
        elif first_incomplete_idx < 0:
            first_incomplete_idx = idx
            first_incomplete_name = display_name
            first_incomplete_stem = stem

    total = len(stems)
    if completed == 0:
        if protocol_was_interrupted(eq_dir, engine):
            resume_idx = first_incomplete_idx if first_incomplete_idx >= 0 else 0
            _, display_name, stem = stems[resume_idx]
            return EquilibrationResumePoint(
                can_resume=True,
                stage_index=resume_idx,
                stage_name=display_name,
                stage_stem=stem,
                mode="stage",
                reason="",
                completed_stages=0,
                total_stages=total,
            )
        return EquilibrationResumePoint(
            can_resume=False,
            reason="No completed stages yet — use Run Equilibration to start",
            completed_stages=0,
            total_stages=total,
        )
    if first_incomplete_idx < 0:
        return EquilibrationResumePoint(
            can_resume=False,
            reason="All stages already completed",
            completed_stages=completed,
            total_stages=total,
        )

    return EquilibrationResumePoint(
        can_resume=True,
        stage_index=first_incomplete_idx,
        stage_name=first_incomplete_name,
        stage_stem=first_incomplete_stem,
        mode="stage",
        reason="",
        completed_stages=completed,
        total_stages=total,
    )
