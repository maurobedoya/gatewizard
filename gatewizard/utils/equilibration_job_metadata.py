"""
Persist and infer equilibration job form metadata (input dir, ensemble, protocol).
"""

from __future__ import annotations

import hashlib
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
    existing: Dict[str, Any] = {}
    path = eq_dir / JOB_METADATA_FILE
    if path.is_file():
        try:
            loaded = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(loaded, dict):
                existing = loaded
        except Exception:
            existing = {}
    payload: Dict[str, Any] = {
        "input_dir": str(Path(input_dir).resolve()),
        "ensemble": ensemble.strip().upper(),
        "protocol": protocol,
        "engine": engine.lower().strip(),
    }
    if openmm_platform:
        payload["openmm_platform"] = openmm_platform
    # Preserve remote execution block across local re-generates
    if isinstance(existing.get("execution"), dict):
        payload["execution"] = existing["execution"]
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


def _is_equilibration_job_dir(path: Path) -> bool:
    return (path / "run_equilibration.sh").is_file()


def _is_builder_dir(path: Path) -> bool:
    """True for GateWizard builder output folders (not equilibration job dirs)."""
    if not path.is_dir() or _is_equilibration_job_dir(path):
        return False
    if (path / "status.json").is_file() or (path / "builder_status.json").is_file():
        return True
    if re.match(r"^\d{2}_build", path.name) and (path / "system.prmtop").is_file():
        return True
    return False


def _looks_like_build_dir(path: Path) -> bool:
    return _is_builder_dir(path)


def _list_builder_dirs(search_roots: List[Path]) -> List[Path]:
    seen: set[str] = set()
    builders: List[Path] = []
    for root in search_roots:
        if not root.is_dir():
            continue
        try:
            children = list(root.iterdir())
        except OSError:
            continue
        for child in children:
            if not child.is_dir():
                continue
            key = str(child.resolve())
            if key in seen or not _is_builder_dir(child):
                continue
            seen.add(key)
            builders.append(child)
    return sorted(builders, key=lambda p: p.name)


def _file_digest(path: Path) -> Optional[str]:
    try:
        return hashlib.md5(path.read_bytes()).hexdigest()
    except OSError:
        return None


def _inpcrd_atom_count(path: Path) -> Optional[int]:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        if len(lines) < 2:
            return None
        return int(lines[1].split()[0])
    except (OSError, ValueError, IndexError):
        return None


def _gro_atom_count(path: Path) -> Optional[int]:
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        if len(lines) < 2:
            return None
        return int(lines[1].strip())
    except (OSError, ValueError, IndexError):
        return None


def _match_builder_by_topology(eq_dir: Path, builders: List[Path]) -> Optional[Path]:
    """Pick the builder folder whose topology matches files copied into the eq job."""
    eq_prmtop = eq_dir / "system.prmtop"
    eq_inpcrd = eq_dir / "system.inpcrd"
    if eq_prmtop.is_file():
        prmtop_hash = _file_digest(eq_prmtop)
        inpcrd_hash = _file_digest(eq_inpcrd) if eq_inpcrd.is_file() else None
        for builder in builders:
            b_prmtop = builder / "system.prmtop"
            if not b_prmtop.is_file():
                continue
            if prmtop_hash and _file_digest(b_prmtop) == prmtop_hash:
                b_inpcrd = builder / "system.inpcrd"
                if inpcrd_hash and b_inpcrd.is_file():
                    if _file_digest(b_inpcrd) == inpcrd_hash:
                        return builder
                else:
                    return builder

    eq_gro = eq_dir / "system.gro"
    if eq_gro.is_file():
        gro_atoms = _gro_atom_count(eq_gro)
        if gro_atoms is not None:
            for builder in builders:
                b_inpcrd = builder / "system.inpcrd"
                if b_inpcrd.is_file() and _inpcrd_atom_count(b_inpcrd) == gro_atoms:
                    return builder
    return None


def _match_builder_by_bilayer(
    eq_dir: Path, builders: List[Path], *, script_text: str = ""
) -> Optional[Path]:
    bilayer_names: List[str] = []
    for path in sorted(eq_dir.glob("bilayer*.pdb")):
        bilayer_names.append(path.name)
    if script_text:
        for match in re.finditer(r'BILAYER_PDB="([^"]+)"', script_text):
            name = match.group(1).strip()
            if name and name not in bilayer_names:
                bilayer_names.append(name)
    if not bilayer_names:
        return None
    for builder in builders:
        if any((builder / name).is_file() for name in bilayer_names):
            return builder
    return None


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
            if _is_builder_dir(candidate):
                return str(candidate.resolve())

    builders = _list_builder_dirs(search_roots)
    if not builders:
        return None

    script_text = ""
    run_script = eq_dir / "run_equilibration.sh"
    if run_script.is_file():
        try:
            script_text = run_script.read_text(encoding="utf-8", errors="replace")
        except OSError:
            pass

    matched = _match_builder_by_topology(eq_dir, builders)
    if matched is not None:
        return str(matched.resolve())
    matched = _match_builder_by_bilayer(eq_dir, builders, script_text=script_text)
    if matched is not None:
        return str(matched.resolve())
    if len(builders) == 1:
        return str(builders[0].resolve())
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


def _mdp_get(text: str, key: str) -> Optional[str]:
    match = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", text, flags=re.M | re.I)
    return match.group(1) if match else None


def _stage_name_from_mdp(path: Path) -> str:
    stem = path.stem  # step0_minimization / step1_equilibration / step7_production
    lower = stem.lower()
    if "minimiz" in lower:
        return "Minimization"
    if "production" in lower:
        return "Production"
    match = re.search(r"step(\d+)", lower)
    if match and "equilibr" in lower:
        # step1 → Equilibration 1
        return f"Equilibration {int(match.group(1))}"
    return stem.replace("_", " ").title()


# GROMACS/OpenMM write kJ/mol/nm²; GUI stores kcal/mol/Å² (same factor as equilibration managers).
_KJ_NM2_TO_KCAL_A2 = 1.0 / 418.4

_STANDARD_CONSTRAINT_KEYS = (
    "protein_backbone",
    "protein_sidechain",
    "lipid_head",
    "lipid_tail",
    "water",
    "ions",
    "other",
)

_AMBER_RESTRAINT_GROUP_NAMES = {
    "protein backbone": "protein_backbone",
    "protein sidechain": "protein_sidechain",
    "lipid head": "lipid_head",
    "lipid tail": "lipid_tail",
    "water": "water",
    "ions": "ions",
    "ion": "ions",
    "other": "other",
}


def _kj_nm2_to_kcal_a2(value: float) -> float:
    return round(float(value) * _KJ_NM2_TO_KCAL_A2, 6)


def _standard_constraint_dict(**overrides: float) -> Dict[str, float]:
    base = {key: 0.0 for key in _STANDARD_CONSTRAINT_KEYS}
    base.update(overrides)
    return base


def _constraints_from_gromacs_mdp(text: str) -> List[Dict[str, Any]]:
    """Positional restraints from GROMACS ``define = -DPOSRES_FC_*`` macros."""
    if "POSRES" not in text:
        return _constraints_dict_to_list(_standard_constraint_dict())
    forces = _standard_constraint_dict()
    macro_map = {
        "POSRES_FC_BB": "protein_backbone",
        "POSRES_FC_SC": "protein_sidechain",
        "POSRES_FC_LIPID": "lipid_head",
        "POSRES_FC_WATER": "water",
        "POSRES_FC_ION": "ions",
        "POSRES_FC_OTHER": "other",
    }
    for macro, key in macro_map.items():
        match = re.search(rf"{macro}=([0-9.]+)", text)
        if not match:
            continue
        kcal = _kj_nm2_to_kcal_a2(float(match.group(1)))
        forces[key] = kcal
        if macro == "POSRES_FC_LIPID":
            forces["lipid_tail"] = kcal
    for match in re.finditer(r"POSRES_FC_([A-Z0-9_]+)=([0-9.]+)", text):
        macro_tail = match.group(1)
        if macro_tail in {"BB", "SC", "LIPID", "WATER", "ION", "OTHER"}:
            continue
        key = macro_tail.lower()
        forces[key] = _kj_nm2_to_kcal_a2(float(match.group(2)))
    return _constraints_dict_to_list(forces)


def _constraints_from_openmm_inp(text: str) -> List[Dict[str, Any]]:
    """Positional restraints from OpenMM ``fc_bb`` / ``fc_sc`` / ``fc_lpos`` (kJ/mol/nm²)."""

    def _fc(name: str) -> float:
        match = re.search(rf"^\s*{re.escape(name)}\s*=\s*([0-9.]+)", text, flags=re.M | re.I)
        return _kj_nm2_to_kcal_a2(float(match.group(1))) if match else 0.0

    lipid = _fc("fc_lpos")
    return _constraints_dict_to_list(
        _standard_constraint_dict(
            protein_backbone=_fc("fc_bb"),
            protein_sidechain=_fc("fc_sc"),
            lipid_head=lipid,
            lipid_tail=lipid,
        )
    )


def _constraints_from_amber_mdin(text: str) -> List[Dict[str, Any]]:
    """Positional restraints from Amber ``ntr=1`` groups after the ``&cntrl`` block."""
    if not re.search(r"\bntr\s*=\s*1\b", text, flags=re.I):
        return _constraints_dict_to_list(_standard_constraint_dict())
    parts = re.split(r"^/\s*$", text, maxsplit=1, flags=re.M)
    if len(parts) < 2:
        return _constraints_dict_to_list(_standard_constraint_dict())
    forces = _standard_constraint_dict()
    lines = parts[1].splitlines()
    idx = 0
    while idx < len(lines):
        label = lines[idx].strip().lower()
        if label in _AMBER_RESTRAINT_GROUP_NAMES:
            key = _AMBER_RESTRAINT_GROUP_NAMES[label]
            idx += 1
            if idx >= len(lines):
                break
            try:
                forces[key] = float(lines[idx].strip())
            except ValueError:
                idx += 1
                continue
            idx += 1
            while idx < len(lines):
                token = lines[idx].strip()
                if not token:
                    idx += 1
                    continue
                if token.lower() in _AMBER_RESTRAINT_GROUP_NAMES or token == "END":
                    break
                idx += 1
            continue
        idx += 1
    return _constraints_dict_to_list(forces)


def _stage_constraints_missing(stage: Dict[str, Any]) -> bool:
    constraints = stage.get("constraints")
    if constraints is None:
        return True
    if isinstance(constraints, dict):
        return len(constraints) == 0
    if isinstance(constraints, list):
        return len(constraints) == 0
    return True


def _patch_protocol_constraints_from_inputs(
    eq_dir: Path, protocol: Dict[str, Any]
) -> Dict[str, Any]:
    """Fill empty stage constraint lists from on-disk MDP/conf/mdin/inp files."""
    if not isinstance(protocol.get("stages"), list):
        return protocol
    if not any(
        isinstance(stage, dict) and _stage_constraints_missing(stage)
        for stage in protocol["stages"]
    ):
        return protocol
    inferred = _infer_protocol_from_engine_inputs(eq_dir)
    if not inferred:
        return protocol
    by_name = {
        str(s.get("name") or "").strip().lower(): s
        for s in inferred.get("stages") or []
        if isinstance(s, dict)
    }
    for stage in protocol["stages"]:
        if not isinstance(stage, dict) or not _stage_constraints_missing(stage):
            continue
        key = str(stage.get("name") or "").strip().lower()
        src = by_name.get(key)
        if src is None:
            continue
        src_constraints = src.get("constraints")
        if isinstance(src_constraints, list) and src_constraints:
            stage["constraints"] = [dict(item) for item in src_constraints if isinstance(item, dict)]
    return protocol


def _stage_dict(
    *,
    name: str,
    ensemble: str,
    time_ns: float,
    steps: int,
    timestep_fs: float,
    temperature: float = 303.15,
    constraints: Optional[List[Dict[str, Any]]] = None,
) -> Dict[str, Any]:
    return {
        "name": name,
        "ensemble": ensemble,
        "time_ns": 0.0 if "minimiz" in name.lower() else round(float(time_ns), 6),
        "steps": int(steps or 0),
        "timestep": float(timestep_fs or 0.0),
        "temperature": float(temperature or 303.15),
        "constraints": constraints if constraints is not None else [],
    }


def _infer_protocol_from_gromacs_mdps(eq_dir: Path) -> Optional[Dict[str, Any]]:
    """Rebuild a minimal GUI protocol from ``step*.mdp`` when summaries are missing."""
    mdps = sorted(eq_dir.glob("step*.mdp"))
    if not mdps:
        return None
    stages: List[Dict[str, Any]] = []
    scheme = "NVT"
    for mdp in mdps:
        try:
            text = mdp.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        nsteps_s = _mdp_get(text, "nsteps")
        dt_s = _mdp_get(text, "dt")  # ps
        temp_s = _mdp_get(text, "ref_t") or _mdp_get(text, "ref-t")
        pcoupl = (_mdp_get(text, "pcoupl") or "no").lower()
        try:
            nsteps = int(float(nsteps_s)) if nsteps_s else 0
        except ValueError:
            nsteps = 0
        try:
            dt_ps = float(dt_s) if dt_s else 0.0
        except ValueError:
            dt_ps = 0.0
        timestep_fs = dt_ps * 1000.0 if dt_ps > 0 else 0.0
        time_ns = (nsteps * dt_ps / 1000.0) if nsteps and dt_ps else 0.0
        try:
            temperature = float(temp_s.split()[0]) if temp_s else 303.15
        except ValueError:
            temperature = 303.15
        if pcoupl not in {"", "no", "none"}:
            scheme = "NPT"
        is_min = "minimiz" in mdp.stem.lower() or _mdp_get(text, "integrator") in {
            "steep",
            "cg",
            "l-bfgs",
        }
        name = _stage_name_from_mdp(mdp)
        stages.append(
            _stage_dict(
                name=name,
                ensemble="minimization" if is_min else scheme,
                time_ns=0.0 if is_min else time_ns,
                steps=nsteps,
                timestep_fs=timestep_fs,
                temperature=temperature,
                constraints=_constraints_from_gromacs_mdp(text),
            )
        )
    if not stages:
        return None
    return {
        "name": f"{scheme} Equilibration Protocol",
        "description": f"{scheme} protocol recovered from GROMACS MDP files",
        "selections": _standard_gui_selections(),
        "stages": stages,
    }


def _infer_protocol_from_namd_confs(eq_dir: Path) -> Optional[Dict[str, Any]]:
    """Rebuild protocol times from ``step*.conf`` (``set time`` / ``set tstep``)."""
    confs = sorted(eq_dir.glob("step*.conf"))
    if not confs:
        return None
    stages: List[Dict[str, Any]] = []
    scheme = "NVT"
    for conf in confs:
        try:
            text = conf.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        lower = conf.stem.lower()
        is_min = "minimiz" in lower
        time_m = re.search(r"set\s+time\s+([0-9.]+)", text, flags=re.I)
        tstep_m = re.search(r"set\s+tstep\s+([0-9.]+)", text, flags=re.I)
        try:
            time_ns = float(time_m.group(1)) if time_m else 0.0
        except ValueError:
            time_ns = 0.0
        try:
            timestep_fs = float(tstep_m.group(1)) if tstep_m else 2.0
        except ValueError:
            timestep_fs = 2.0
        steps = int(round(time_ns * 1_000_000 / timestep_fs)) if time_ns > 0 and timestep_fs > 0 else 0
        if re.search(r"langevinPiston\s+on", text, flags=re.I) or re.search(
            r"useFlexibleCell\s+yes", text, flags=re.I
        ):
            scheme = "NPT"
        name = _stage_name_from_mdp(conf)
        stages.append(
            _stage_dict(
                name=name,
                ensemble="minimization" if is_min else scheme,
                time_ns=0.0 if is_min else time_ns,
                steps=0 if is_min else steps,
                timestep_fs=timestep_fs,
            )
        )
    if not stages:
        return None
    return {
        "name": f"{scheme} Equilibration Protocol",
        "description": f"{scheme} protocol recovered from NAMD conf files",
        "selections": _standard_gui_selections(),
        "stages": stages,
    }


def _amber_mdin_get(text: str, key: str) -> Optional[str]:
    match = re.search(rf"\b{re.escape(key)}\s*=\s*([^,!\n]+)", text, flags=re.I)
    return match.group(1).strip() if match else None


def _infer_protocol_from_amber_mdins(eq_dir: Path) -> Optional[Dict[str, Any]]:
    mdins = sorted(eq_dir.glob("step*.mdin"))
    if not mdins:
        return None
    stages: List[Dict[str, Any]] = []
    scheme = "NVT"
    for mdin in mdins:
        try:
            text = mdin.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        lower = mdin.stem.lower()
        is_min = "minimiz" in lower or (_amber_mdin_get(text, "imin") or "0") != "0"
        nstlim_s = _amber_mdin_get(text, "nstlim") or _amber_mdin_get(text, "maxcyc")
        dt_s = _amber_mdin_get(text, "dt")
        temp_s = _amber_mdin_get(text, "temp0")
        try:
            nsteps = int(float(nstlim_s)) if nstlim_s else 0
        except ValueError:
            nsteps = 0
        try:
            dt_ps = float(dt_s) if dt_s else 0.0
        except ValueError:
            dt_ps = 0.0
        timestep_fs = dt_ps * 1000.0 if dt_ps > 0 else 0.0
        time_ns = (nsteps * dt_ps / 1000.0) if nsteps and dt_ps else 0.0
        try:
            temperature = float(temp_s) if temp_s else 303.15
        except ValueError:
            temperature = 303.15
        if re.search(r"\bntp\s*=\s*[1-9]", text, flags=re.I):
            scheme = "NPT"
        name = _stage_name_from_mdp(mdin)
        stages.append(
            _stage_dict(
                name=name,
                ensemble="minimization" if is_min else scheme,
                time_ns=0.0 if is_min else time_ns,
                steps=nsteps,
                timestep_fs=timestep_fs or 2.0,
                temperature=temperature,
                constraints=_constraints_from_amber_mdin(text),
            )
        )
    if not stages:
        return None
    return {
        "name": f"{scheme} Equilibration Protocol",
        "description": f"{scheme} protocol recovered from Amber mdin files",
        "selections": _standard_gui_selections(),
        "stages": stages,
    }


def _openmm_inp_get(text: str, key: str) -> Optional[str]:
    match = re.search(rf"^\s*{re.escape(key)}\s*=\s*(\S+)", text, flags=re.M | re.I)
    return match.group(1) if match else None


def _infer_protocol_from_openmm_inps(eq_dir: Path) -> Optional[Dict[str, Any]]:
    """Rebuild protocol times from ``step*.inp`` (``nstep`` / ``dt``)."""
    inps = sorted(eq_dir.glob("step*.inp"))
    if not inps:
        return None
    stages: List[Dict[str, Any]] = []
    scheme = "NVT"
    for inp in inps:
        try:
            text = inp.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        lower = inp.stem.lower()
        is_min = "minimiz" in lower
        nstep_s = _openmm_inp_get(text, "nstep")
        dt_s = _openmm_inp_get(text, "dt")
        temp_s = _openmm_inp_get(text, "temp") or _openmm_inp_get(text, "gen_temp")
        try:
            nsteps = int(float(nstep_s)) if nstep_s else 0
        except ValueError:
            nsteps = 0
        try:
            dt_ps = float(dt_s) if dt_s else 0.0
        except ValueError:
            dt_ps = 0.0
        timestep_fs = dt_ps * 1000.0 if dt_ps > 0 else 0.0
        time_ns = (nsteps * dt_ps / 1000.0) if nsteps and dt_ps else 0.0
        try:
            temperature = float(temp_s) if temp_s else 303.15
        except ValueError:
            temperature = 303.15
        pcouple = (_openmm_inp_get(text, "pcouple") or "no").lower()
        if pcouple not in {"", "no", "none"}:
            scheme = "NPT"
        name = _stage_name_from_mdp(inp)
        stages.append(
            _stage_dict(
                name=name,
                ensemble="minimization" if is_min else scheme,
                time_ns=0.0 if is_min else time_ns,
                steps=nsteps,
                timestep_fs=timestep_fs or 2.0,
                temperature=temperature,
                constraints=_constraints_from_openmm_inp(text),
            )
        )
    if not stages:
        return None
    return {
        "name": f"{scheme} Equilibration Protocol",
        "description": f"{scheme} protocol recovered from OpenMM inp files",
        "selections": _standard_gui_selections(),
        "stages": stages,
    }


def _infer_protocol_from_engine_inputs(eq_dir: Path) -> Optional[Dict[str, Any]]:
    """Best-effort protocol from on-disk engine inputs (MDP / conf / mdin / inp)."""
    return (
        _infer_protocol_from_gromacs_mdps(eq_dir)
        or _infer_protocol_from_namd_confs(eq_dir)
        or _infer_protocol_from_amber_mdins(eq_dir)
        or _infer_protocol_from_openmm_inps(eq_dir)
    )


def _production_stage_from_file(eq_dir: Path) -> Optional[Dict[str, Any]]:
    """Read only the production input file — fast path for time_ns correction."""
    mdp = eq_dir / "step7_production.mdp"
    if mdp.is_file():
        try:
            text = mdp.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        if text:
            nsteps_s = _mdp_get(text, "nsteps")
            dt_s = _mdp_get(text, "dt")
            try:
                nsteps = int(float(nsteps_s)) if nsteps_s else 0
            except ValueError:
                nsteps = 0
            try:
                dt_ps = float(dt_s) if dt_s else 0.0
            except ValueError:
                dt_ps = 0.0
            time_ns = (nsteps * dt_ps / 1000.0) if nsteps and dt_ps else 0.0
            return _stage_dict(
                name="Production",
                ensemble="NVT",
                time_ns=time_ns,
                steps=nsteps,
                timestep_fs=dt_ps * 1000.0 if dt_ps > 0 else 2.0,
            )

    conf = eq_dir / "step7_production.conf"
    if conf.is_file():
        try:
            text = conf.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        if text:
            time_m = re.search(r"set\s+time\s+([0-9.]+)", text, flags=re.I)
            tstep_m = re.search(r"set\s+tstep\s+([0-9.]+)", text, flags=re.I)
            try:
                time_ns = float(time_m.group(1)) if time_m else 0.0
            except ValueError:
                time_ns = 0.0
            try:
                timestep_fs = float(tstep_m.group(1)) if tstep_m else 2.0
            except ValueError:
                timestep_fs = 2.0
            steps = (
                int(round(time_ns * 1_000_000 / timestep_fs))
                if time_ns > 0 and timestep_fs > 0
                else 0
            )
            return _stage_dict(
                name="Production",
                ensemble="NVT",
                time_ns=time_ns,
                steps=steps,
                timestep_fs=timestep_fs,
            )

    mdin = eq_dir / "step7_production.mdin"
    if mdin.is_file():
        try:
            text = mdin.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        if text:
            nstlim_s = _amber_mdin_get(text, "nstlim") or _amber_mdin_get(text, "maxcyc")
            dt_s = _amber_mdin_get(text, "dt")
            try:
                nsteps = int(float(nstlim_s)) if nstlim_s else 0
            except ValueError:
                nsteps = 0
            try:
                dt_ps = float(dt_s) if dt_s else 0.0
            except ValueError:
                dt_ps = 0.0
            time_ns = (nsteps * dt_ps / 1000.0) if nsteps and dt_ps else 0.0
            return _stage_dict(
                name="Production",
                ensemble="NVT",
                time_ns=time_ns,
                steps=nsteps,
                timestep_fs=dt_ps * 1000.0 if dt_ps > 0 else 2.0,
            )

    inp = eq_dir / "step7_production.inp"
    if inp.is_file():
        try:
            text = inp.read_text(encoding="utf-8", errors="replace")
        except OSError:
            text = ""
        if text:
            nstep_s = _openmm_inp_get(text, "nstep")
            dt_s = _openmm_inp_get(text, "dt")
            try:
                nsteps = int(float(nstep_s)) if nstep_s else 0
            except ValueError:
                nsteps = 0
            try:
                dt_ps = float(dt_s) if dt_s else 0.0
            except ValueError:
                dt_ps = 0.0
            time_ns = (nsteps * dt_ps / 1000.0) if nsteps and dt_ps else 0.0
            return _stage_dict(
                name="Production",
                ensemble="NVT",
                time_ns=time_ns,
                steps=nsteps,
                timestep_fs=dt_ps * 1000.0 if dt_ps > 0 else 2.0,
            )
    return None


def _infer_engine_from_dir(eq_dir: Path) -> Optional[str]:
    resources_file = eq_dir / "equilibration_resources.json"
    if resources_file.is_file():
        try:
            data = json.loads(resources_file.read_text(encoding="utf-8", errors="replace"))
            eng = data.get("engine") if isinstance(data, dict) else None
            if isinstance(eng, str) and eng.strip():
                return eng.strip().lower()
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            pass
    if (eq_dir / "openmm_run.py").is_file() or list(eq_dir.glob("step*_equilibration.inp")):
        return "openmm"
    if list(eq_dir.glob("*.mdp")) or (eq_dir / "index.ndx").is_file():
        return "gromacs"
    if list(eq_dir.glob("step*.mdin")) or list(eq_dir.glob("*.mdout")):
        return "amber"
    if list(eq_dir.glob("step*_equilibration.conf")) or list(
        eq_dir.glob("step*_minimization.conf")
    ):
        return "namd"
    run_script = eq_dir / "run_equilibration.sh"
    if run_script.is_file():
        try:
            script = run_script.read_text(encoding="utf-8", errors="replace").lower()
            if "openmm_run.py" in script or "from openmm" in script:
                return "openmm"
            if "gmx" in script or "grompp" in script or "mdrun" in script:
                return "gromacs"
            if "pmemd" in script or "sander" in script or "amber=" in script:
                return "amber"
            if "namd" in script:
                return "namd"
        except OSError:
            pass
    return None


def _stage_time_missing(stage: Dict[str, Any]) -> bool:
    try:
        return float(stage.get("time_ns") or 0) <= 0 and int(stage.get("steps") or 0) <= 0
    except (TypeError, ValueError):
        return True


def _protocol_needs_input_time_patch(protocol: Dict[str, Any]) -> bool:
    """True when we still need to open MDP/conf/mdin (avoid OneDrive scans when complete)."""
    stages = protocol.get("stages")
    if not isinstance(stages, list) or not stages:
        return True
    have_prod = False
    for stage in stages:
        if not isinstance(stage, dict):
            continue
        name = str(stage.get("name") or "").lower()
        if "production" in name:
            have_prod = True
            try:
                if float(stage.get("time_ns") or 0) <= 0 and int(stage.get("steps") or 0) <= 0:
                    return True
            except (TypeError, ValueError):
                return True
        elif _stage_time_missing(stage) and "minimiz" not in name:
            # Equilibration stages with no time — worth a patch pass.
            return True
    return not have_prod


def _sync_production_stage_from_disk(
    eq_dir: Path, protocol: Dict[str, Any]
) -> bool:
    """Overwrite Production time/steps from step7_* when present. Returns True if updated."""
    if not isinstance(protocol.get("stages"), list):
        return False
    prod_src = _production_stage_from_file(eq_dir)
    if not prod_src:
        return False
    updated = False
    for stage in protocol["stages"]:
        if not isinstance(stage, dict):
            continue
        if "production" not in str(stage.get("name") or "").lower():
            continue
        src_time = float(prod_src.get("time_ns") or 0)
        src_steps = int(prod_src.get("steps") or 0)
        cur_time = float(stage.get("time_ns") or 0)
        cur_steps = int(stage.get("steps") or 0)
        if src_time > 0 and abs(cur_time - src_time) > 1e-6:
            stage["time_ns"] = prod_src["time_ns"]
            updated = True
        if src_steps > 0 and cur_steps != src_steps:
            stage["steps"] = prod_src["steps"]
            updated = True
        if prod_src.get("timestep") and not stage.get("timestep"):
            stage["timestep"] = prod_src["timestep"]
            updated = True
        return updated or src_time > 0 or src_steps > 0
    protocol["stages"].append(dict(prod_src))
    return True


def _patch_protocol_times_from_inputs(
    eq_dir: Path, protocol: Dict[str, Any]
) -> Dict[str, Any]:
    """Fill missing production/eq times from engine inputs; append Production if absent."""
    if not isinstance(protocol.get("stages"), list):
        return protocol
    _sync_production_stage_from_disk(eq_dir, protocol)
    # Skip scanning all step*.mdp/conf/mdin when the protocol already has times
    # (Use in form on OneDrive folders was hanging on this redundant IO).
    if not _protocol_needs_input_time_patch(protocol):
        return _patch_protocol_constraints_from_inputs(eq_dir, protocol)
    inferred = _infer_protocol_from_engine_inputs(eq_dir)
    if not inferred:
        return _patch_protocol_constraints_from_inputs(eq_dir, protocol)
    by_name = {
        str(s.get("name") or "").strip().lower(): s
        for s in inferred.get("stages") or []
        if isinstance(s, dict)
    }
    for stage in protocol["stages"]:
        if not isinstance(stage, dict):
            continue
        key = str(stage.get("name") or "").strip().lower()
        src = by_name.get(key)
        if src is None and "production" in key:
            src = next((s for s in by_name.values() if "production" in s.get("name", "").lower()), None)
        if src is None:
            continue
        if _stage_time_missing(stage) or (
            "production" in key
            and (
                float(stage.get("time_ns") or 0) <= 0
                or abs(float(stage.get("time_ns") or 0) - float(src.get("time_ns") or 0))
                > 1e-6
            )
        ):
            if src.get("time_ns") is not None:
                stage["time_ns"] = src["time_ns"]
            if src.get("steps"):
                stage["steps"] = src["steps"]
            if src.get("timestep"):
                stage["timestep"] = src["timestep"]
    have_prod = any(
        "production" in str(s.get("name") or "").lower()
        for s in protocol["stages"]
        if isinstance(s, dict)
    )
    if not have_prod:
        prod = next(
            (
                s
                for s in inferred.get("stages") or []
                if isinstance(s, dict) and "production" in str(s.get("name") or "").lower()
            ),
            None,
        )
        if prod:
            protocol["stages"].append(dict(prod))
    return _patch_protocol_constraints_from_inputs(eq_dir, protocol)


def _infer_protocol(eq_dir: Path) -> Optional[Dict[str, Any]]:
    summary_file = eq_dir / "protocol_summary.json"
    if summary_file.is_file():
        try:
            summary = json.loads(summary_file.read_text(encoding="utf-8", errors="replace"))
            protocol = _protocol_from_namd_summary(summary)
            if protocol:
                protocol = _patch_protocol_times_from_inputs(eq_dir, protocol)
                return _patch_protocol_constraints_from_inputs(eq_dir, protocol)
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            pass
    inferred = _infer_protocol_from_engine_inputs(eq_dir)
    if inferred:
        return _patch_protocol_constraints_from_inputs(eq_dir, inferred)
    return None


def infer_equilibration_job_metadata(
    eq_dir: Path,
    working_dir: Optional[Path] = None,
    *,
    heal: bool = True,
) -> Dict[str, Any]:
    """Best-effort recovery of form metadata for an existing equilibration job.

    Cluster Watching may leave ``equilibration_job.json`` with only an
    ``execution`` block. In that case fall back to ``protocol_summary.json``
    (and folder heuristics) for protocol / ensemble / input_dir — and optionally
    write the recovered fields back so **Use in form** keeps working.
    """
    eq_dir = Path(eq_dir)
    work = Path(working_dir) if working_dir else None
    resolved_input: Optional[str] = None
    ensemble: Optional[str] = None
    protocol: Optional[Dict[str, Any]] = None
    engine: Optional[str] = None
    existing: Dict[str, Any] = {}

    metadata_file = eq_dir / JOB_METADATA_FILE
    if metadata_file.is_file():
        try:
            data = json.loads(metadata_file.read_text(encoding="utf-8", errors="replace"))
            if isinstance(data, dict):
                existing = data
                input_dir = data.get("input_dir")
                if isinstance(input_dir, str) and input_dir.strip():
                    resolved_input = str(Path(input_dir).resolve())
                ens = data.get("ensemble")
                if isinstance(ens, str) and ens.strip():
                    ensemble = ens.strip().upper()
                protocol = _normalize_gui_protocol(data.get("protocol"))
                eng = data.get("engine")
                if isinstance(eng, str) and eng.strip():
                    engine = eng.strip().lower()
        except (json.JSONDecodeError, OSError, TypeError, ValueError):
            existing = {}

    if not resolved_input:
        resolved_input = _infer_input_dir(eq_dir, work)
    if not ensemble:
        ensemble = _infer_ensemble(eq_dir)
    if not protocol:
        protocol = _infer_protocol(eq_dir)
    elif isinstance(protocol, dict):
        # Job JSON / summary may omit production time_ns — recover from inputs.
        protocol = _patch_protocol_times_from_inputs(eq_dir, protocol)
        protocol = _patch_protocol_constraints_from_inputs(eq_dir, protocol)
    if not engine:
        engine = _infer_engine_from_dir(eq_dir)
    if not ensemble and protocol and isinstance(protocol.get("stages"), list):
        for stage in protocol["stages"]:
            if not isinstance(stage, dict):
                continue
            ens = stage.get("ensemble")
            if isinstance(ens, str) and ens.strip() and ens.strip().lower() != "minimization":
                ensemble = ens.strip().upper()
                break
        if not ensemble and isinstance(protocol.get("name"), str):
            for scheme in ("NPGT", "NPAT", "NPT", "NVT"):
                if scheme in protocol["name"].upper():
                    ensemble = "NPgT" if scheme == "NPGT" else scheme
                    break

    # Persist recovered GUI fields when the job JSON was thinned to execution-only.
    if heal and metadata_file.is_file() and (protocol or ensemble or resolved_input or engine):
        missing_protocol = not _normalize_gui_protocol(existing.get("protocol"))
        missing_ensemble = not (
            isinstance(existing.get("ensemble"), str) and existing["ensemble"].strip()
        )
        missing_input = not (
            isinstance(existing.get("input_dir"), str) and existing["input_dir"].strip()
        )
        missing_engine = not (
            isinstance(existing.get("engine"), str) and existing["engine"].strip()
        )
        if missing_protocol or missing_ensemble or missing_input or missing_engine:
            payload = dict(existing)
            if protocol and missing_protocol:
                payload["protocol"] = protocol
            if ensemble and missing_ensemble:
                payload["ensemble"] = ensemble
            if resolved_input and missing_input:
                payload["input_dir"] = resolved_input
            if engine and missing_engine:
                payload["engine"] = engine
            try:
                metadata_file.write_text(
                    json.dumps(payload, indent=2), encoding="utf-8"
                )
            except OSError:
                pass

    return {
        "input_dir": resolved_input,
        "ensemble": ensemble,
        "protocol": protocol,
        "engine": engine,
    }
