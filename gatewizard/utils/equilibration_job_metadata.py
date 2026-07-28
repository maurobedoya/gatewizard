"""
Persist and infer equilibration job form metadata (input dir, ensemble, protocol).
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

JOB_METADATA_FILE = "equilibration_job.json"
_EQ_FOLDER_PREFIX = "03_equilibration_"
_ENGINE_SUFFIXES = ("_openmm", "_gromacs", "_namd")

_CONSTRAINT_DISPLAY_NAMES = {
    "protein_backbone": "Protein backbone",
    "protein_sidechain": "Protein sidechain",
    "lipid_head": "Lipid head",
    "lipid_tail": "Lipid tail",
    "water": "Water",
    "ions": "Ions",
    "ion": "Ions",
    "other": "Other",
}


def write_equilibration_job_metadata(
    eq_dir: Path,
    *,
    input_dir: str,
    ensemble: str,
    protocol: Dict[str, Any],
    engine: str,
    openmm_platform: Optional[str] = None,
) -> Path:
    """Persist GUI form state used to generate this equilibration job."""
    eq_dir = Path(eq_dir)
    payload: Dict[str, Any] = {
        "input_dir": str(Path(input_dir).resolve()),
        "ensemble": ensemble.strip().upper(),
        "protocol": protocol,
        "engine": engine.lower().strip(),
    }
    if openmm_platform:
        payload["openmm_platform"] = openmm_platform
    path = eq_dir / JOB_METADATA_FILE
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path


def _constraints_dict_to_list(constraints: Any) -> List[Dict[str, Any]]:
    if isinstance(constraints, list):
        return [dict(item) for item in constraints if isinstance(item, dict)]
    if not isinstance(constraints, dict):
        return []
    items: List[Dict[str, Any]] = []
    for key, force in constraints.items():
        key_str = str(key)
        items.append(
            {
                "name": _CONSTRAINT_DISPLAY_NAMES.get(
                    key_str, key_str.replace("_", " ").title()
                ),
                "force_constant": float(force),
                "selection": key_str,
            }
        )
    return items


def _stage_dict_to_gui_stage(stage: Dict[str, Any]) -> Dict[str, Any]:
    constraints = stage.get("constraints")
    gui_stage = {
        key: value
        for key, value in stage.items()
        if key != "constraints"
    }
    gui_stage["constraints"] = _constraints_dict_to_list(constraints)
    return gui_stage


def _standard_gui_selections() -> Dict[str, str]:
    """Default MDAnalysis aliases used by the GUI protocol editor."""
    return {
        "protein_backbone": "protein and backbone",
        "protein_sidechain": "protein and not backbone",
        "lipid_head": (
            "(resname POPC POPE POPS DPPC DMPC DOPC DSPC PC PE PS PA PG PI SM "
            "OL LA MY ST AR OLE PAL STE LIN CHOL CHL CHOLEST PALM OLEO STEROL) "
            "and (name P O11 O12 O13 O14 O21 O22 O31 O32 O33 O34 O1P O2P O3P O4P "
            "OP1 OP2 OP3 OP4 N C11 C12 C13 C14 N31 C32 C33 C34 C35 C1 C2 C3 "
            "HN1 HN2 HN3 HO2 HO3 HS)"
        ),
        "lipid_tail": (
            "(resname POPC POPE POPS DPPC DMPC DOPC DSPC PC PE PS PA PG PI SM "
            "OL LA MY ST AR OLE PAL STE LIN CHOL CHL CHOLEST PALM OLEO STEROL) "
            "and not (name P O11 O12 O13 O14 O21 O22 O31 O32 O33 O34 O1P O2P O3P "
            "O4P OP1 OP2 OP3 OP4 N C11 C12 C13 C14 N31 C32 C33 C34 C35 C1 C2 C3 "
            "HN1 HN2 HN3 HO2 HO3 HS)"
        ),
        "water": "resname TIP3 HOH WAT SOL TIP4 SPC T3P T4P",
        "ions": (
            "resname NA CL K CA MG ZN FE CU SOD CLA POT CAL MAG ZIN IRN COP "
            "Na+ Cl- K+ Ca2+ Mg2+ Zn2+ Fe2+ Fe3+ Cu2+ NA+ CL- LIT RUB CES BAR"
        ),
        "other": (
            "not (protein or (resname POPC POPE POPS DPPC DMPC DOPC DSPC PC PE "
            "PS PA PG PI SM OL LA MY ST AR OLE PAL STE LIN CHOL CHL CHOLEST "
            "PALM OLEO STEROL) or (resname TIP3 HOH WAT SOL TIP4 SPC T3P T4P) "
            "or (resname NA CL K CA MG ZN FE CU SOD CLA POT CAL MAG ZIN IRN "
            "COP Na+ Cl- K+ Ca2+ Mg2+ Zn2+ Fe2+ Fe3+ Cu2+ NA+ CL- LIT RUB CES BAR))"
        ),
    }


def _protocol_from_namd_summary(summary: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    stages_raw = summary.get("stages")
    if not stages_raw:
        return None

    if isinstance(stages_raw, dict):
        stage_items = list(stages_raw.values())
    elif isinstance(stages_raw, list):
        stage_items = stages_raw
    else:
        return None

    stages: List[Dict[str, Any]] = []
    for item in stage_items:
        if isinstance(item, dict):
            stages.append(_stage_dict_to_gui_stage(item))

    if not stages:
        return None

    return {
        "name": summary.get("protocol_name") or "Equilibration Protocol",
        "description": summary.get("description")
        or f"{summary.get('scheme_type', 'Equilibration')} protocol recovered from job folder",
        "selections": _standard_gui_selections(),
        "stages": stages,
    }


def _normalize_gui_protocol(protocol: Any) -> Optional[Dict[str, Any]]:
    """Ensure protocol is GUI-shaped (list constraints + selections map)."""
    if not isinstance(protocol, dict):
        return None
    stages_raw = protocol.get("stages")
    if not isinstance(stages_raw, list) or not stages_raw:
        return None
    stages = [
        _stage_dict_to_gui_stage(stage) if isinstance(stage, dict) else stage
        for stage in stages_raw
        if isinstance(stage, dict)
    ]
    if not stages:
        return None
    selections = protocol.get("selections")
    if not isinstance(selections, dict) or not selections:
        selections = _standard_gui_selections()
    return {
        "name": protocol.get("name") or "Equilibration Protocol",
        "description": protocol.get("description") or "",
        "selections": selections,
        "stages": stages,
    }


def _candidate_input_basenames(job_name: str) -> List[str]:
    if not job_name.startswith(_EQ_FOLDER_PREFIX):
        return []
    base = job_name[len(_EQ_FOLDER_PREFIX) :]
    if not base:
        return []
    candidates = [base]
    for suffix in _ENGINE_SUFFIXES:
        if base.endswith(suffix):
            stripped = base[: -len(suffix)]
            if stripped and stripped not in candidates:
                candidates.append(stripped)
    return candidates


def _looks_like_build_dir(path: Path) -> bool:
    if not path.is_dir():
        return False
    markers = (
        "system.prmtop",
        "system.inpcrd",
        "system.rst",
        "status.json",
        "builder_status.json",
    )
    return any((path / name).is_file() for name in markers)


def _infer_input_dir(eq_dir: Path, working_dir: Optional[Path]) -> Optional[str]:
    search_roots: List[Path] = []
    if working_dir is not None:
        search_roots.append(Path(working_dir))
    parent = eq_dir.parent
    if parent not in search_roots:
        search_roots.append(parent)

    for basename in _candidate_input_basenames(eq_dir.name):
        for root in search_roots:
            candidate = root / basename
            if _looks_like_build_dir(candidate):
                return str(candidate.resolve())

    return None


def _infer_ensemble(eq_dir: Path) -> Optional[str]:
    summary_file = eq_dir / "protocol_summary.json"
    if summary_file.is_file():
        try:
            summary = json.loads(summary_file.read_text(encoding="utf-8", errors="replace"))
            scheme = summary.get("scheme_type")
            if isinstance(scheme, str) and scheme.strip():
                return scheme.strip().upper()
            stages = summary.get("stages")
            if isinstance(stages, dict):
                for stage in stages.values():
                    if isinstance(stage, dict):
                        ens = stage.get("ensemble")
                        if isinstance(ens, str) and ens.strip():
                            return ens.strip().upper()
            elif isinstance(stages, list):
                for stage in stages:
                    if isinstance(stage, dict):
                        ens = stage.get("ensemble")
                        if isinstance(ens, str) and ens.strip():
                            return ens.strip().upper()
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            pass

    run_script = eq_dir / "run_equilibration.sh"
    if run_script.is_file():
        try:
            text = run_script.read_text(encoding="utf-8", errors="replace")
            for scheme in ("NPgT", "NPAT", "NPT", "NVT"):
                if re.search(rf"\b{scheme}\b", text):
                    return scheme
        except OSError:
            pass

    return None


def _infer_protocol(eq_dir: Path) -> Optional[Dict[str, Any]]:
    summary_file = eq_dir / "protocol_summary.json"
    if summary_file.is_file():
        try:
            summary = json.loads(summary_file.read_text(encoding="utf-8", errors="replace"))
            protocol = _protocol_from_namd_summary(summary)
            if protocol:
                return protocol
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            pass
    return None


def infer_equilibration_job_metadata(
    eq_dir: Path,
    working_dir: Optional[Path] = None,
) -> Dict[str, Any]:
    """Best-effort recovery of form metadata for an existing equilibration job."""
    eq_dir = Path(eq_dir)
    metadata_file = eq_dir / JOB_METADATA_FILE
    if metadata_file.is_file():
        try:
            data = json.loads(metadata_file.read_text(encoding="utf-8", errors="replace"))
            if isinstance(data, dict):
                input_dir = data.get("input_dir")
                ensemble = data.get("ensemble")
                protocol = data.get("protocol")
                resolved_input = (
                    str(Path(input_dir).resolve())
                    if isinstance(input_dir, str) and input_dir.strip()
                    else None
                )
                return {
                    "input_dir": resolved_input,
                    "ensemble": ensemble.strip().upper()
                    if isinstance(ensemble, str) and ensemble.strip()
                    else None,
                    "protocol": _normalize_gui_protocol(protocol),
                }
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            pass

    return {
        "input_dir": _infer_input_dir(eq_dir, Path(working_dir) if working_dir else None),
        "ensemble": _infer_ensemble(eq_dir),
        "protocol": _infer_protocol(eq_dir),
    }
