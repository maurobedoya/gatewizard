# gatewizard/tools/peptide_cap_parametrization.py
# SPDX-License-Identifier: MIT
"""GAFF parametrization for nonstandard peptide polymer termini (FVA/FOR/ETA).

Unlike free ligands, these residues stay covalently linked in the chain and need
head/tail connect atoms in the Amber OFF library for tleap loadPDB.
"""

from __future__ import annotations

import json
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gatewizard.tools.ligand_parametrization import (
    DEFAULT_ATOM_TYPE,
    DEFAULT_CHARGE_METHOD,
    LigandParametrizationError,
    _write_status,
)
from gatewizard.utils.helpers import get_clean_env, resolve_conda_executable
from gatewizard.utils.logger import get_logger
from gatewizard.utils.peptide_residues import (
    AMBER_UNSUPPORTED_PEPTIDE_WARN,
    PEPTIDE_CAP_CONNECT,
    detect_peptide_caps_needing_gaff,
    normalize_resname,
)

logger = get_logger(__name__)


def extract_polymer_residue_pdb(
    pdb_file: str, residue_name: str, output_dir: str
) -> str:
    """Extract one copy of a polymer residue (ATOM or HETATM) to its own PDB."""
    pdb_path = Path(pdb_file).resolve()
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    output_pdb = out_dir / f"{residue_name}.pdb"
    target = normalize_resname(residue_name)

    copies: Dict[Tuple[str, str], List[str]] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 20:
                continue
            if normalize_resname(line[17:20]) != target:
                continue
            chain = line[21:22] if len(line) > 21 else " "
            res_seq = line[22:27].strip() if len(line) > 26 else ""
            copies.setdefault((chain, res_seq), []).append(line)

    if not copies:
        raise LigandParametrizationError(
            f"Polymer residue '{residue_name}' not found in {pdb_file}"
        )

    first_key = sorted(copies.keys())[0]
    lines = copies[first_key]
    if len(copies) > 1:
        logger.info(
            "Found %s copies of %s; extracting first (chain=%s, resSeq=%s, %s atoms)",
            len(copies),
            residue_name,
            first_key[0],
            first_key[1],
            len(lines),
        )

    with open(output_pdb, "w", encoding="utf-8") as handle:
        handle.writelines(lines)
        handle.write("END\n")
    return str(output_pdb)


def _pdb_has_hydrogen_atoms(pdb_file: str) -> bool:
    with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 20:
                continue
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            name = line[12:16].strip().upper()
            if element == "H" or name.startswith("H"):
                return True
    return False


def _normalize_cap_pdb_resname(pdb_file: str, residue_name: str, output_pdb: str) -> str:
    """Force a single residue name/id and contiguous serials for antechamber."""
    name = normalize_resname(residue_name)
    out_lines: list[str] = []
    serial = 1
    with open(pdb_file, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                atom_name = line[12:16]
                alt_loc = line[16:17] if len(line) > 16 else " "
                coords = line[30:54]
                rest = line[54:].rstrip("\n") if len(line) > 54 else ""
                if len(rest) < 22:
                    rest = rest + (" " * (22 - len(rest)))
                out_lines.append(
                    f"HETATM{serial:5d} {atom_name}{alt_loc}{name:>3} A   1    "
                    f"{coords}{rest}\n"
                )
                serial += 1
            elif line.startswith("END"):
                out_lines.append("END\n")
    if not any(line.startswith("END") for line in out_lines):
        out_lines.append("END\n")
    Path(output_pdb).write_text("".join(out_lines), encoding="utf-8")
    return output_pdb


def _add_hydrogens_rdkit(input_pdb: str, output_pdb: str, residue_name: str) -> bool:
    # Optional dependency — only needed for heavy-atom deposited caps.
    try:
        from rdkit import Chem
    except ImportError:
        return False
    mol = Chem.MolFromPDBFile(str(input_pdb), removeHs=False, sanitize=False)
    if mol is None or mol.GetNumAtoms() == 0:
        return False
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        # Still try AddHs; some exotic caps sanitize poorly
        pass
    try:
        mol_h = Chem.AddHs(mol, addCoords=True)
    except Exception as exc:
        logger.warning("RDKit AddHs failed for %s: %s", residue_name, exc)
        return False
    tmp = Path(output_pdb).with_suffix(".rdkit.pdb")
    Chem.MolToPDBFile(mol_h, str(tmp))
    _normalize_cap_pdb_resname(str(tmp), residue_name, output_pdb)
    try:
        tmp.unlink(missing_ok=True)
    except OSError:
        pass
    return Path(output_pdb).is_file() and _pdb_has_hydrogen_atoms(output_pdb)


def _add_hydrogens_reduce(input_pdb: str, output_pdb: str, residue_name: str) -> bool:
    try:
        reduce_exe = resolve_conda_executable("reduce")
    except Exception:
        reduce_exe = shutil.which("reduce") or ""
    if not reduce_exe:
        return False
    try:
        result = subprocess.run(
            [reduce_exe, "-NOFLIP", str(input_pdb)],
            capture_output=True,
            text=True,
            timeout=120,
            env=get_clean_env(),
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("reduce failed for %s: %s", residue_name, exc)
        return False
    # reduce often returns 1 with warnings but still prints PDB
    if not result.stdout.strip():
        return False
    tmp = Path(output_pdb).with_suffix(".reduce.pdb")
    tmp.write_text(result.stdout, encoding="utf-8")
    _normalize_cap_pdb_resname(str(tmp), residue_name, output_pdb)
    try:
        tmp.unlink(missing_ok=True)
    except OSError:
        pass
    return Path(output_pdb).is_file() and _pdb_has_hydrogen_atoms(output_pdb)


def ensure_cap_hydrogens(pdb_file: str, residue_name: str) -> str:
    """Ensure the extracted cap PDB has hydrogens for antechamber/SQM.

    Deposited peptide HETs (ETA/FVA) are often heavy-atom only. Without H,
    antechamber reports an odd electron count and SQM aborts.
    """
    pdb_path = Path(pdb_file).resolve()
    if _pdb_has_hydrogen_atoms(str(pdb_path)):
        normalized = pdb_path.parent / f"{normalize_resname(residue_name)}_norm.pdb"
        return _normalize_cap_pdb_resname(str(pdb_path), residue_name, str(normalized))

    out_pdb = pdb_path.parent / f"{normalize_resname(residue_name)}_h.pdb"
    if _add_hydrogens_rdkit(str(pdb_path), str(out_pdb), residue_name):
        logger.info("Added hydrogens to peptide cap %s with RDKit", residue_name)
        return str(out_pdb)
    if _add_hydrogens_reduce(str(pdb_path), str(out_pdb), residue_name):
        logger.info("Added hydrogens to peptide cap %s with reduce", residue_name)
        return str(out_pdb)
    raise LigandParametrizationError(
        f"Peptide cap {residue_name} has no hydrogens and GateWizard could not "
        "add them (need RDKit or Amber reduce). Add hydrogens to the PDB and retry."
    )


def _tleap_head_tail_lines(residue_name: str) -> str:
    meta = PEPTIDE_CAP_CONNECT.get(normalize_resname(residue_name), {})
    head = meta.get("head_atom")
    tail = meta.get("tail_atom")
    lines: List[str] = []
    if head:
        lines.append(f"set {residue_name} head {residue_name}.1.{head}")
    else:
        lines.append(f"set {residue_name} head null")
    if tail:
        lines.append(f"set {residue_name} tail {residue_name}.1.{tail}")
    else:
        lines.append(f"set {residue_name} tail null")
    return "\n".join(lines)


def parametrize_peptide_cap(
    residue_pdb: str,
    residue_name: str,
    output_dir: str,
    charge: int = 0,
    charge_method: str = DEFAULT_CHARGE_METHOD,
    multiplicity: int = 1,
    atom_type: str = DEFAULT_ATOM_TYPE,
    sqm_keywords: str = "",
) -> Dict[str, Any]:
    """Parametrize a polymer cap with antechamber/parmchk2/tleap (head/tail set)."""
    name = normalize_resname(residue_name)
    if name not in AMBER_UNSUPPORTED_PEPTIDE_WARN:
        raise LigandParametrizationError(
            f"Residue '{residue_name}' is not a supported GAFF peptide cap "
            f"({sorted(AMBER_UNSUPPORTED_PEPTIDE_WARN)})"
        )

    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    residue_pdb = str(Path(residue_pdb).resolve())

    log_dir = out_dir / "logs"
    log_dir.mkdir(exist_ok=True)

    mol2_file = out_dir / f"{name}.mol2"
    frcmod_file = out_dir / f"{name}.frcmod"
    lib_file = out_dir / f"{name}.lib"
    prmtop_file = out_dir / f"{name}.prmtop"
    inpcrd_file = out_dir / f"{name}.inpcrd"
    status_file = out_dir / "status.json"
    status: Dict[str, Any] = {
        "residue_name": name,
        "kind": "peptide_cap",
        "charge": charge,
        "charge_method": charge_method,
        "atom_type": atom_type,
        "status": "running",
        "current_step": "antechamber",
        "steps_completed": [],
        "start_time": datetime.now().isoformat(),
        "error": None,
    }
    _write_status(status_file, status)

    try:
        residue_pdb = ensure_cap_hydrogens(residue_pdb, name)
        status["input_pdb_prepared"] = residue_pdb
        _write_status(status_file, status)

        def _run_antechamber(method: str) -> subprocess.CompletedProcess:
            cmd = [
                "antechamber",
                "-i",
                residue_pdb,
                "-fi",
                "pdb",
                "-o",
                str(mol2_file),
                "-fo",
                "mol2",
                "-c",
                method,
                "-nc",
                str(charge),
                "-m",
                str(multiplicity),
                "-rn",
                name,
                "-s",
                "2",
                "-at",
                atom_type,
            ]
            ek = sqm_keywords.strip()
            if ek:
                cmd.extend(["-ek", ek])
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                cwd=str(out_dir),
                timeout=600,
                env=get_clean_env(),
            )
            with open(log_dir / "antechamber.log", "a", encoding="utf-8") as handle:
                handle.write(f"COMMAND: {' '.join(cmd)}\n\n")
                handle.write(f"STDOUT:\n{proc.stdout}\n\nSTDERR:\n{proc.stderr}\n\n")
            return proc

        (log_dir / "antechamber.log").write_text("", encoding="utf-8")
        result = _run_antechamber(charge_method)
        used_method = charge_method
        if (result.returncode != 0 or not mol2_file.exists()) and charge_method != "gas":
            # SQM-based methods often fail on awkward caps; Gasteiger is a fallback.
            logger.warning(
                "Antechamber (%s) failed for peptide cap %s; retrying with -c gas",
                charge_method,
                name,
            )
            if mol2_file.exists():
                mol2_file.unlink()
            result = _run_antechamber("gas")
            used_method = "gas"
        if result.returncode != 0 or not mol2_file.exists():
            sqm_hint = ""
            sqm_out = out_dir / "sqm.out"
            if sqm_out.is_file():
                sqm_hint = f" See {sqm_out} for SQM details."
            raise LigandParametrizationError(
                f"Antechamber failed for peptide cap {name}: {result.stderr}{sqm_hint}"
            )
        status["charge_method"] = used_method
        if used_method != charge_method:
            status["charge_method_fallback"] = True

        status["steps_completed"].append("antechamber")
        status["current_step"] = "parmchk2"
        _write_status(status_file, status)

        parmchk_cmd = [
            "parmchk2",
            "-i",
            str(mol2_file),
            "-f",
            "mol2",
            "-o",
            str(frcmod_file),
            "-s",
            atom_type,
        ]
        result = subprocess.run(
            parmchk_cmd,
            capture_output=True,
            text=True,
            cwd=str(out_dir),
            timeout=120,
            env=get_clean_env(),
        )
        with open(log_dir / "parmchk2.log", "w", encoding="utf-8") as handle:
            handle.write(f"COMMAND: {' '.join(parmchk_cmd)}\n\n")
            handle.write(f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n")
        if result.returncode != 0 or not frcmod_file.exists():
            raise LigandParametrizationError(
                f"Parmchk2 failed for peptide cap {name}: {result.stderr}"
            )

        status["steps_completed"].append("parmchk2")
        status["current_step"] = "tleap"
        _write_status(status_file, status)

        leaprc = f"leaprc.{atom_type}"
        head_tail = _tleap_head_tail_lines(name)
        tleap_content = f"""source {leaprc}
loadamberparams {name}.frcmod
{name} = loadmol2 {name}.mol2
{head_tail}
check {name}
saveoff {name} {name}.lib
saveamberparm {name} {name}.prmtop {name}.inpcrd
quit
"""
        tleap_input = out_dir / "tleap.in"
        tleap_input.write_text(tleap_content, encoding="utf-8")
        result = subprocess.run(
            ["tleap", "-f", "tleap.in"],
            capture_output=True,
            text=True,
            cwd=str(out_dir),
            timeout=120,
            env=get_clean_env(),
        )
        with open(log_dir / "tleap.log", "w", encoding="utf-8") as handle:
            handle.write("COMMAND: tleap -f tleap.in\n\n")
            handle.write(f"INPUT:\n{tleap_content}\n\n")
            handle.write(f"STDOUT:\n{result.stdout}\n\nSTDERR:\n{result.stderr}\n")
        if result.returncode != 0 or not lib_file.exists():
            raise LigandParametrizationError(
                f"tleap failed for peptide cap {name}: {result.stderr}"
            )

        status["steps_completed"].append("tleap")
        status["status"] = "completed"
        status["current_step"] = "completed"
        status["end_time"] = datetime.now().isoformat()
        _write_status(status_file, status)

        return {
            "mol2": str(mol2_file),
            "frcmod": str(frcmod_file),
            "lib": str(lib_file),
            "prmtop": str(prmtop_file),
            "inpcrd": str(inpcrd_file),
            "charge": charge,
            "atom_type": atom_type,
            "name": name,
            "role": PEPTIDE_CAP_CONNECT.get(name, {}).get("role", "unknown"),
        }
    except LigandParametrizationError as exc:
        status["status"] = "error"
        status["error"] = str(exc)
        status["end_time"] = datetime.now().isoformat()
        _write_status(status_file, status)
        raise
    except subprocess.TimeoutExpired as exc:
        msg = f"Timeout during peptide-cap parametrization of {name}"
        status["status"] = "error"
        status["error"] = msg
        status["end_time"] = datetime.now().isoformat()
        _write_status(status_file, status)
        raise LigandParametrizationError(msg) from exc
    except Exception as exc:
        msg = f"Unexpected error parametrizing peptide cap {name}: {exc}"
        status["status"] = "error"
        status["error"] = msg
        status["end_time"] = datetime.now().isoformat()
        _write_status(status_file, status)
        raise LigandParametrizationError(msg) from exc


def parametrize_peptide_cap_from_system_pdb(
    pdb_file: str,
    residue_name: str,
    output_dir: str,
    charge: int = 0,
    charge_method: str = DEFAULT_CHARGE_METHOD,
    multiplicity: int = 1,
    atom_type: str = DEFAULT_ATOM_TYPE,
    sqm_keywords: str = "",
) -> Dict[str, Any]:
    """Extract and parametrize one polymer cap under ``peptide_cap_params/{name}/``."""
    name = normalize_resname(residue_name)
    base_dir = Path(output_dir).resolve()
    cap_dir = base_dir / "peptide_cap_params" / name
    cap_dir.mkdir(parents=True, exist_ok=True)
    cap_pdb = extract_polymer_residue_pdb(pdb_file, name, str(cap_dir))
    return parametrize_peptide_cap(
        residue_pdb=cap_pdb,
        residue_name=name,
        output_dir=str(cap_dir),
        charge=charge,
        charge_method=charge_method,
        multiplicity=multiplicity,
        atom_type=atom_type,
        sqm_keywords=sqm_keywords,
    )


def parametrize_all_peptide_caps_from_system_pdb(
    pdb_file: str,
    output_dir: str,
    charges: Optional[Dict[str, int]] = None,
    charge_method: str = DEFAULT_CHARGE_METHOD,
    atom_type: str = DEFAULT_ATOM_TYPE,
    sqm_keywords: str = "",
) -> Dict[str, Dict[str, Any]]:
    """Detect and parametrize all GAFF-needed polymer caps in *pdb_file*."""
    charges = charges or {}
    caps = detect_peptide_caps_needing_gaff(pdb_file)
    results: Dict[str, Dict[str, Any]] = {}
    for cap in caps:
        name = cap["name"]
        results[name] = parametrize_peptide_cap_from_system_pdb(
            pdb_file=pdb_file,
            residue_name=name,
            output_dir=output_dir,
            charge=int(charges.get(name, 0)),
            charge_method=charge_method,
            atom_type=atom_type,
            sqm_keywords=sqm_keywords,
        )
    return results


def peptide_cap_params_charge_delta(
    peptide_cap_params: Dict[str, Dict[str, Any]],
) -> int:
    """Sum formal charges for packmol-memgen ``--charge_pdb_delta``."""
    total = 0
    for files in peptide_cap_params.values():
        try:
            total += int(files.get("charge", 0) or 0)
        except (TypeError, ValueError):
            continue
    return total


def check_peptide_cap_parametrization(
    output_dir: str, cap_names: List[str]
) -> Dict[str, Dict[str, Any]]:
    """Return cached completed peptide-cap params under ``peptide_cap_params/``."""
    base = Path(output_dir).resolve() / "peptide_cap_params"
    found: Dict[str, Dict[str, Any]] = {}
    for name in cap_names:
        name_u = normalize_resname(name)
        cap_dir = base / name_u
        frcmod = cap_dir / f"{name_u}.frcmod"
        lib = cap_dir / f"{name_u}.lib"
        status_path = cap_dir / "status.json"
        charge = 0
        if status_path.is_file():
            try:
                data = json.loads(status_path.read_text(encoding="utf-8"))
                if data.get("status") != "completed":
                    continue
                charge = int(data.get("charge", 0) or 0)
            except (OSError, ValueError, json.JSONDecodeError):
                continue
        if frcmod.is_file() and lib.is_file():
            found[name_u] = {
                "frcmod": str(frcmod),
                "lib": str(lib),
                "charge": charge,
                "name": name_u,
                "role": PEPTIDE_CAP_CONNECT.get(name_u, {}).get("role", "unknown"),
            }
    return found
