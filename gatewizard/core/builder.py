# gatewizard/core/builder.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
This module handles the preparation of membrane protein systems using
packmol-memgen and related tools.
"""

import json
import os
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from gatewizard.utils.logger import get_logger
from gatewizard.utils.helpers import get_clean_env, get_clean_env_shell_snippet, resolve_conda_executable, subprocess_argv_for_script
from gatewizard.core.preparation import (
    _resolve_pdb4amber_executable,
    count_protein_hydrogens,
    strip_protein_hydrogens,
)
from gatewizard.tools.force_fields import ForceFieldManager
from gatewizard.tools.validators import SystemValidator
from gatewizard.tools.ligand_parametrization import (
    build_ligand_param_args,
    build_tleap_ligand_lines,
)

logger = get_logger(__name__)


def optional_pdb_path(pdb_file: Optional[str]) -> Optional[str]:
    """Return a non-empty PDB path, or ``None`` when building without a protein."""
    if pdb_file is None:
        return None
    text = str(pdb_file).strip()
    if not text or text.lower() in {"none", "null"}:
        return None
    return text


def normalize_builder_solutes(config: Optional[Dict[str, Any]] = None) -> List[Dict[str, Any]]:
    """Normalize free-molecule entries for packmol-memgen ``--solute`` / ``--solute_con``.

    Accepts ``solutes`` as a list of dicts (``pdb``, ``concentration``, optional
    ``in_membrane`` / ``prot_dist`` / ``name``) or the shorthand ``solute`` plus
    ``solute_con`` strings/lists used in the Amber tutorial.
    """
    cfg = config or {}
    raw = cfg.get("solutes")
    out: List[Dict[str, Any]] = []
    if raw:
        items = raw if isinstance(raw, (list, tuple)) else [raw]
        for item in items:
            if isinstance(item, str):
                path = item.strip()
                if path:
                    out.append(
                        {
                            "pdb": path,
                            "concentration": "1",
                            "in_membrane": False,
                            "prot_dist": None,
                            "name": "",
                        }
                    )
                continue
            if not isinstance(item, dict):
                continue
            path = str(item.get("pdb") or item.get("path") or "").strip()
            if not path:
                continue
            prot_dist = item.get("prot_dist", item.get("solute_prot_dist"))
            out.append(
                {
                    "pdb": path,
                    "concentration": str(
                        item.get("concentration", item.get("solute_con", "1"))
                    ).strip()
                    or "1",
                    "in_membrane": bool(
                        item.get("in_membrane", item.get("solute_inmem", False))
                    ),
                    "prot_dist": prot_dist,
                    "name": str(item.get("name") or "").strip(),
                }
            )
        return out

    solute = cfg.get("solute")
    if not solute:
        return []
    paths = [solute] if isinstance(solute, str) else list(solute)
    cons = cfg.get("solute_con", "1")
    if isinstance(cons, (list, tuple)):
        cons_list = [str(c).strip() or "1" for c in cons]
    else:
        cons_list = [str(cons).strip() or "1"]
    in_mem = bool(cfg.get("solute_inmem", False))
    prot_dist = cfg.get("solute_prot_dist")
    for i, path in enumerate(paths):
        text = str(path).strip()
        if not text:
            continue
        con = cons_list[i] if i < len(cons_list) else cons_list[-1]
        out.append(
            {
                "pdb": text,
                "concentration": con,
                "in_membrane": in_mem,
                "prot_dist": prot_dist,
                "name": "",
            }
        )
    return out


def effective_distxy_fix(config: Optional[Dict[str, Any]] = None) -> Optional[float]:
    """Membrane XY size (Å) for packmol-memgen ``--distxy_fix``.

    Uses ``distxy_fix`` when set; otherwise the X length of ``dims``.
    Required by packmol-memgen when no protein PDB is given.
    """
    cfg = config or {}
    raw = cfg.get("distxy_fix")
    if raw is not None and str(raw).strip() != "":
        try:
            value = float(raw)
        except (TypeError, ValueError):
            return None
        return value if value > 0 else None
    dims = cfg.get("dims")
    if dims and len(dims) >= 1:
        try:
            value = float(dims[0])
        except (TypeError, ValueError):
            return None
        return value if value > 0 else None
    return None


class BuilderError(Exception):
    """Custom exception for system building errors."""

    pass


class Builder:
    """
    Main class for building membrane protein systems.

    This class handles the complete workflow of preparing membrane protein
    systems including lipid composition, force field selection, and system
    parameterization.
    """

    def __init__(self):
        """Initialize the Builder."""
        self.ff_manager = ForceFieldManager()
        self.validator = SystemValidator()

        # Default configuration
        self.config = {
            "water_model": "tip3p",
            "protein_ff": "ff19SB",
            "lipid_ff": "lipid21",
            "preoriented": True,
            "parametrize": True,
            "salt_concentration": 0.15,
            "cation": "K+",
            "anion": "Cl-",
            "dist": 12,  # Minimum solute-to-box-boundary distance in Angstroms
            "dist_wat": 26,  # Default water layer thickness in Angstroms
            "dims": None,  # Explicit box dimensions [X, Y, Z] in Angstroms
            "notprotonate": True,  # Preserve PropKa residue names (GLH/ASH/…); skip reduce re-protonation
            "remove_protein_h": False,  # Strip protein H before packmol-memgen (keep ligand/hetero H)
            "two_stage_process": False,  # Enable two-stage packing + parametrization
            "pack_only": False,  # Only perform packing stage
            "parametrize_only": False,  # Only perform parametrization stage
            "ligand_params": {},  # Dict of ligand name -> {frcmod, lib} file paths
            "nloop": 20,  # GENCAN loops for PACKMOL (packmol-memgen default)
            "nloop_all": 100,  # GENCAN loops for all-together packing
            "tolerance": 2.0,  # PACKMOL clash tolerance (radius1+radius2)
            "md_engine": None,  # Target MD engine: namd, gromacs, openmm, amber, or None
            "distxy_fix": None,  # Fixed membrane XY size (Å); required without a protein PDB
            "solutes": [],  # Free molecules: [{pdb, concentration, in_membrane, prot_dist, name}]
            "solute_inmem": False,  # Place free molecules in the membrane (else water)
            "solute_prot_dist": None,  # Keep free molecules this far from the protein (Å)
        }

    def set_configuration(self, **kwargs):
        """Update configuration parameters."""
        self.config.update(kwargs)
        logger.info(f"Configuration updated: {kwargs}")

    def validate_system_inputs(
        self,
        pdb_file: Optional[str] = None,
        upper_lipids: Optional[List[str]] = None,
        lower_lipids: Optional[List[str]] = None,
        lipid_ratios: str = "",
        **kwargs,
    ) -> Tuple[bool, str]:
        """
        Validate all inputs for system preparation.

        Args:
            pdb_file: Path to input PDB file, or empty/None for a bilayer-only build
            upper_lipids: List of lipids for upper leaflet
            lower_lipids: List of lipids for lower leaflet
            lipid_ratios: Ratio string (format: upper_ratios//lower_ratios)
            **kwargs: Additional parameters to validate

        Returns:
            Tuple of (is_valid, error_message). When valid with a non-empty
            message, the message is a warning (e.g. force-field note or
            protein hydrogens detected).
        """
        pdb_path = optional_pdb_path(pdb_file)
        valid, error_msg = self.validator.validate_system_inputs(
            pdb_path, upper_lipids or [], lower_lipids or [], lipid_ratios, **kwargs
        )
        if not valid:
            return False, error_msg

        warnings: List[str] = []
        if error_msg:
            warnings.append(error_msg)

        # Foreign / non-Amber protein H often break tleap after packmol-memgen.
        # Skip the warning when the user already opted to strip them, or when
        # there is no protein PDB (bilayer-only / free-molecule builds).
        if pdb_path and not kwargs.get("remove_protein_h", False):
            try:
                n_h = count_protein_hydrogens(pdb_path)
            except OSError as exc:
                logger.warning(f"Could not scan protein hydrogens: {exc}")
                n_h = 0
            if n_h > 0:
                warnings.append(
                    f"Protein has {n_h} hydrogen atom(s). Hydrogens from other "
                    "software (e.g. Schrödinger) can break tleap after "
                    "packmol-memgen. Open Advanced settings and enable "
                    "'Remove protein hydrogens' before generating input files. "
                    "Only protein hydrogens are removed; ligands and other "
                    "heteroatoms keep theirs."
                )

        solutes = normalize_builder_solutes(kwargs)
        if solutes and kwargs.get("parametrize", True):
            ligand_params = kwargs.get("ligand_params") or {}
            named = [s.get("name") or "" for s in solutes]
            missing = [n for n in named if n and n not in ligand_params]
            if not ligand_params or missing:
                warnings.append(
                    "Free molecules need GAFF parameters (frcmod/lib) when "
                    "parametrizing. Parametrize each molecule before generating "
                    "inputs, or pass ligand_params / --ligand_param."
                )

        return True, "\n\n".join(warnings)

    def prepare_system(
        self,
        pdb_file: Optional[str] = None,
        working_dir: str = "",
        upper_lipids: Optional[List[str]] = None,
        lower_lipids: Optional[List[str]] = None,
        lipid_ratios: str = "",
        **kwargs,
    ) -> Tuple[bool, str, Optional[Path]]:
        """
        Prepare a complete membrane protein system.

        Args:
            pdb_file: Path to input PDB file, or empty/None for a bilayer-only build
            working_dir: Working directory for outputs
            upper_lipids: List of lipids for upper leaflet
            lower_lipids: List of lipids for lower leaflet
            lipid_ratios: Lipid ratios string
            **kwargs: Additional configuration parameters.  Notable extras:

                * **wait** (*bool*) – Block until the job finishes or errors
                  (default ``False``).
                * **wait_timeout** (*float | None*) – Max seconds to wait
                  (default ``None`` = unlimited).
                * **wait_poll_interval** (*float*) – Seconds between status
                  checks when *wait=True* (default ``5.0``).
                * **wait_verbose** (*bool*) – Print progress while waiting
                  (default ``True``).
                * **distxy_fix** (*float*) – Membrane XY size in Å. Required
                  when *pdb_file* is omitted.
                * **solutes** (*list*) – Free molecules to pack (``pdb`` +
                  ``concentration``; see ``normalize_builder_solutes``).

        Returns:
            Tuple of (success, message, job_directory)
        """
        # Update configuration with provided parameters
        config = {**self.config, **kwargs}
        upper_lipids = upper_lipids or []
        lower_lipids = lower_lipids or []

        # Validate inputs
        valid, error_msg = self.validate_system_inputs(
            pdb_file, upper_lipids, lower_lipids, lipid_ratios, **config
        )
        if not valid:
            return False, error_msg, None

        try:
            job_dir, local_pdb, config = self._prepare_job_workspace(
                pdb_file, working_dir, config
            )

            # Build packmol-memgen command
            cmd = self._build_command(
                local_pdb, upper_lipids, lower_lipids, lipid_ratios, config
            )

            # Create and launch preparation script
            script_path = self._generate_preparation_files(cmd, job_dir, config)
            success = self._run_preparation_script(job_dir, script_path, config)

            if success:
                # If wait=True, block until the job finishes
                wait = config.get("wait", False)
                if wait:
                    wait_timeout = config.get("wait_timeout", None)
                    wait_poll = config.get("wait_poll_interval", 5.0)
                    wait_verbose = config.get("wait_verbose", True)
                    completed, wait_msg = self.wait_for_completion(
                        job_dir,
                        poll_interval=wait_poll,
                        timeout=wait_timeout,
                        verbose=wait_verbose,
                    )
                    return completed, wait_msg, job_dir

                return True, "System preparation started successfully", job_dir
            else:
                return False, "Failed to start system preparation", job_dir

        except Exception as e:
            logger.error(f"Error preparing system: {e}", exc_info=True)
            return False, f"Error preparing system: {str(e)}", None

    def generate_preparation_inputs(
        self,
        pdb_file: Optional[str] = None,
        working_dir: str = "",
        upper_lipids: Optional[List[str]] = None,
        lower_lipids: Optional[List[str]] = None,
        lipid_ratios: str = "",
        **kwargs,
    ) -> Tuple[bool, str, Optional[Path]]:
        """
        Generate preparation input files without launching the job.

        Creates a job directory with ``run_preparation.sh`` and ``status.json``
        (status ``not_started``).  Use :meth:`run_preparation` to execute the job.

        *pdb_file* may be empty/None for a bilayer-only build (set
        ``distxy_fix`` or ``dims``). Free molecules use ``solutes``.
        """
        config = {**self.config, **kwargs}
        upper_lipids = upper_lipids or []
        lower_lipids = lower_lipids or []

        valid, error_msg = self.validate_system_inputs(
            pdb_file, upper_lipids, lower_lipids, lipid_ratios, **config
        )
        if not valid:
            return False, error_msg, None

        try:
            job_dir, local_pdb, config = self._prepare_job_workspace(
                pdb_file, working_dir, config
            )
            cmd = self._build_command(
                local_pdb, upper_lipids, lower_lipids, lipid_ratios, config
            )
            self._generate_preparation_files(cmd, job_dir, config)
            return True, "Preparation input files generated successfully", job_dir

        except Exception as e:
            logger.error(f"Error generating preparation inputs: {e}", exc_info=True)
            return False, f"Error generating preparation inputs: {str(e)}", None

    def run_preparation(self, job_dir: str | Path) -> Tuple[bool, str]:
        """
        Launch a previously generated preparation job.

        Args:
            job_dir: Path to the job directory containing ``run_preparation.sh``.

        Returns:
            Tuple of (success, message)
        """
        job_path = Path(job_dir).resolve()
        script_path = job_path / "run_preparation.sh"
        if not script_path.is_file():
            return False, "run_preparation.sh not found — generate input files first"

        status_file = job_path / "status.json"
        config: Dict[str, Any] = {}
        if status_file.is_file():
            try:
                with open(status_file, "r", encoding="utf-8") as f:
                    status_data = json.load(f)
                if status_data.get("status") == "running":
                    return False, "Preparation is already running"
                config = status_data.get("config", {})
            except (json.JSONDecodeError, OSError) as e:
                logger.warning(f"Could not read status.json: {e}")

        success = self._run_preparation_script(job_path, script_path, config)
        if success:
            return True, "System preparation started successfully"
        return False, "Failed to start system preparation"

    def cancel_preparation(self, job_dir: str | Path) -> Dict[str, Any]:
        """
        Cancel a running preparation job from the GUI.

        Kills the detached ``run_preparation.sh`` process group (PID in
        ``process.pid``) and marks ``status.json`` as ``cancelled``.
        """
        job_path = Path(job_dir).expanduser().resolve()
        if not job_path.is_dir():
            raise FileNotFoundError(f"Job directory not found: {job_path}")
        status_file = job_path / "status.json"
        if not status_file.is_file() and not (job_path / "run_preparation.sh").is_file():
            raise FileNotFoundError(f"Not a preparation job directory: {job_path}")

        status = "unknown"
        status_data: Dict[str, Any] = {}
        if status_file.is_file():
            try:
                status_data = json.loads(status_file.read_text(encoding="utf-8"))
                status = status_data.get("status", "unknown")
            except (json.JSONDecodeError, OSError):
                status_data = {}

        if status in {"completed", "error", "cancelled", "not_started"}:
            return {
                "success": True,
                "job_dir": str(job_path),
                "stopped": False,
                "status": status,
                "message": f"Job already finished ({status})",
            }

        pid = self._read_preparation_pid(job_path)
        killed = False
        if pid is not None and self._pid_alive(pid):
            killed = self._kill_process_group(pid)

        status_data["status"] = "cancelled"
        status_data["error"] = "Cancelled by user"
        status_data["end_time"] = datetime.now().isoformat()
        try:
            status_file.write_text(
                json.dumps(status_data, indent=2), encoding="utf-8"
            )
        except OSError as exc:
            logger.warning("Could not update status.json after cancel: %s", exc)

        log_path = job_path / "logs" / "preparation.log"
        try:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with open(log_path, "a", encoding="utf-8") as fh:
                fh.write(
                    f"\n[{datetime.now().isoformat()}] CANCELLED by user"
                    f"{' (process killed)' if killed else ''}\n"
                )
        except OSError:
            pass

        (job_path / "process.pid").unlink(missing_ok=True)
        logger.info("Cancelled preparation job %s (killed=%s)", job_path.name, killed)
        return {
            "success": True,
            "job_dir": str(job_path),
            "stopped": True,
            "killed_process": killed,
            "status": "cancelled",
            "message": "Cancel requested",
        }

    @staticmethod
    def _read_preparation_pid(job_dir: Path) -> Optional[int]:
        pid_path = job_dir / "process.pid"
        if not pid_path.is_file():
            return None
        try:
            text = pid_path.read_text(encoding="utf-8").strip()
            return int(text) if text else None
        except (OSError, ValueError):
            return None

    @staticmethod
    def _pid_alive(pid: Optional[int]) -> bool:
        if pid is None or pid <= 0:
            return False
        try:
            os.kill(pid, 0)
            return True
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        except OSError:
            return False

    @staticmethod
    def _kill_process_group(pid: int) -> bool:
        """SIGTERM then SIGKILL the detached preparation process group."""
        killed = False
        for sig in (signal.SIGTERM, signal.SIGKILL):
            try:
                os.killpg(os.getpgid(pid), sig)
                killed = True
            except ProcessLookupError:
                return True
            except (OSError, AttributeError):
                try:
                    os.kill(pid, sig)
                    killed = True
                except ProcessLookupError:
                    return True
                except OSError:
                    pass
            if sig == signal.SIGTERM:
                deadline = time.time() + 2.0
                while time.time() < deadline:
                    if not Builder._pid_alive(pid):
                        return True
                    time.sleep(0.1)
        return killed or not Builder._pid_alive(pid)

    def wait_for_completion(
        self,
        job_dir,
        poll_interval: float = 5.0,
        timeout: Optional[float] = None,
        verbose: bool = True,
    ) -> Tuple[bool, str]:
        """
        Block until a preparation job completes or fails.

        This method polls the job's ``status.json`` file and optionally
        verifies the PID is still alive.  It is useful for scripting
        multiple sequential preparations.

        Args:
            job_dir: Path to the job directory (string or Path).
            poll_interval: Seconds between status checks (default 5).
            timeout: Maximum seconds to wait.  ``None`` means no limit.
            verbose: Print elapsed-time progress to stdout.

        Returns:
            Tuple of (success, message).
        """
        job_dir = Path(job_dir)
        status_file = job_dir / "status.json"
        pid_file = job_dir / "process.pid"
        start_wait = time.time()

        if verbose:
            print(f"Waiting for job to complete: {job_dir.name}")

        while True:
            # --- timeout check ---
            elapsed = time.time() - start_wait
            if timeout is not None and elapsed > timeout:
                msg = f"Timeout after {timeout:.0f}s waiting for {job_dir.name}"
                if verbose:
                    print(f"\n{msg}")
                return False, msg

            # --- read status.json ---
            step_label = ""
            if status_file.exists():
                try:
                    with open(status_file, "r", encoding="utf-8") as f:
                        content = f.read().strip()
                    if content:
                        status_data = json.loads(content)
                        job_status = status_data.get("status", "unknown")
                        step_label = status_data.get("current_step", "")

                        if job_status == "completed":
                            msg = (
                                f"Job completed successfully: {job_dir.name} "
                                f"({elapsed:.1f}s)"
                            )
                            if verbose:
                                print(f"\n{msg}")
                            return True, msg

                        if job_status == "error":
                            error = status_data.get("error", "Unknown error")
                            msg = f"Job failed: {error}"
                            if verbose:
                                print(f"\n{msg}")
                            return False, msg
                except (json.JSONDecodeError, IOError):
                    pass  # file may be mid-write

            # --- check PID still alive ---
            if pid_file.exists():
                try:
                    pid = int(pid_file.read_text().strip())
                    os.kill(pid, 0)  # signal 0 = existence check
                except (ProcessLookupError, PermissionError):
                    # Process gone — give status.json one last chance to update
                    time.sleep(2)
                    if status_file.exists():
                        try:
                            with open(status_file, "r", encoding="utf-8") as f:
                                status_data = json.loads(f.read().strip())
                            if status_data.get("status") == "completed":
                                msg = (
                                    f"Job completed successfully: {job_dir.name} "
                                    f"({elapsed:.1f}s)"
                                )
                                if verbose:
                                    print(f"\n{msg}")
                                return True, msg
                            if status_data.get("status") == "error":
                                error = status_data.get("error", "Unknown error")
                                msg = f"Job failed: {error}"
                                if verbose:
                                    print(f"\n{msg}")
                                return False, msg
                        except (json.JSONDecodeError, IOError):
                            pass
                    msg = (
                        f"Process terminated without updating status for "
                        f"{job_dir.name}"
                    )
                    if verbose:
                        print(f"\n{msg}")
                    return False, msg
                except (ValueError, IOError):
                    pass  # PID file not readable yet

            # --- progress line ---
            if verbose:
                mins, secs = divmod(int(elapsed), 60)
                step_str = f" [step {step_label}]" if step_label else ""
                print(
                    f"\r  Elapsed: {mins:02d}:{secs:02d}{step_str}   ",
                    end="",
                    flush=True,
                )

            time.sleep(poll_interval)

    def _prepare_job_workspace(
        self,
        pdb_file: Optional[str],
        working_dir: str,
        config: Dict[str, Any],
    ) -> Tuple[Path, Optional[Path], Dict[str, Any]]:
        """Create the job directory and copy protein / free-molecule PDBs."""
        pdb_path = optional_pdb_path(pdb_file)
        custom_output_name = config.get("output_folder_name", None)
        job_dir = self._create_job_directory(
            pdb_path, working_dir, custom_output_name
        )

        local_pdb: Optional[Path] = None
        if pdb_path:
            local_pdb = self._prepare_input_files(
                pdb_path, job_dir, remove_protein_h=bool(config.get("remove_protein_h"))
            )
            if config.get("notprotonate", False):
                logger.info(
                    "Using notprotonate mode - assuming PDB has correct protonation states"
                )
                config["_workflow_note"] = (
                    "IMPORTANT: When using --notprotonate, ensure your PDB file has correct protonation states and residue names (e.g., GLU->GLH for protonated glutamate)"
                )

        solutes = normalize_builder_solutes(config)
        if solutes:
            config["solutes"] = self._prepare_solute_files(job_dir, solutes)

        config["_has_protein"] = local_pdb is not None
        return job_dir, local_pdb, config

    def _create_job_directory(
        self,
        pdb_file: Optional[str],
        working_dir: str,
        custom_output_name: Optional[str] = None,
    ) -> Path:
        """Create unique job directory."""
        work_dir = Path(working_dir).resolve()
        work_dir.mkdir(parents=True, exist_ok=True)

        if custom_output_name:
            # Use custom output name as-is, only add timestamp if directory already exists
            job_dir = work_dir / custom_output_name
            if job_dir.exists():
                # Only add timestamp if directory already exists to ensure uniqueness
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                job_dir = work_dir / f"{custom_output_name}_{timestamp}"
        else:
            # Use default naming scheme
            pdb_path = optional_pdb_path(pdb_file)
            if pdb_path:
                pdb_name = os.path.splitext(os.path.basename(pdb_path))[0]
                stem = f"02_build_{pdb_name}"
            else:
                stem = "02_build_bilayer"
            job_dir = work_dir / stem
            if job_dir.exists():
                timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                job_dir = work_dir / f"{stem}_{timestamp}"

        job_dir.mkdir(parents=True, exist_ok=True)

        # Create logs directory
        (job_dir / "logs").mkdir(exist_ok=True)

        logger.info(f"Created job directory: {job_dir}")
        return job_dir

    def _prepare_input_files(
        self,
        pdb_file: str,
        job_dir: Path,
        remove_protein_h: bool = False,
    ) -> Path:
        """Copy input PDB into the job directory; optionally strip protein H."""
        pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
        local_pdb = job_dir / f"{pdb_name}.pdb"
        shutil.copy2(pdb_file, local_pdb)

        if remove_protein_h:
            result = strip_protein_hydrogens(str(local_pdb), str(local_pdb))
            logger.info(
                "Removed %s protein hydrogen atom(s) from Builder input PDB",
                result["removed"],
            )

        logger.info(f"Copied PDB file to: {local_pdb}")
        return local_pdb

    def _prepare_solute_files(
        self, job_dir: Path, solutes: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Copy free-molecule PDBs into the job directory."""
        dest_dir = job_dir / "solutes"
        dest_dir.mkdir(exist_ok=True)
        local: List[Dict[str, Any]] = []
        for index, solute in enumerate(solutes):
            src = Path(solute["pdb"]).expanduser().resolve()
            dest = dest_dir / src.name
            if dest.exists() and dest.resolve() != src:
                dest = dest_dir / f"{src.stem}_{index}{src.suffix}"
            shutil.copy2(src, dest)
            local.append({**solute, "pdb": str(dest)})
            logger.info(f"Copied solute PDB to: {dest}")
        return local

    def _build_command(
        self,
        pdb_file: Optional[Path],
        upper_lipids: List[str],
        lower_lipids: List[str],
        lipid_ratios: str,
        config: Dict[str, Any],
    ) -> List[str]:
        """Build packmol-memgen command."""
        cmd = ["packmol-memgen"]
        if pdb_file is not None:
            cmd.extend(["--pdb", str(pdb_file)])

        # Add lipids (only needed for packing stage or single-stage process)
        if not config.get("parametrize_only", False):
            lipids_arg = self._prepare_lipids_argument(upper_lipids, lower_lipids)
            if lipids_arg:
                cmd.extend(["--lipids", lipids_arg])

            # Add ratios (only for packing stage)
            if lipid_ratios:
                cmd.extend(["--ratio", lipid_ratios])

        # Add force fields
        cmd.extend(["--ffwat", config["water_model"]])
        cmd.extend(["--ffprot", config["protein_ff"]])
        cmd.extend(["--fflip", config["lipid_ff"]])

        # Add flags. --preoriented is only meaningful with a protein PDB
        # (otherwise packmol-memgen would run MEMEMBED on nothing).
        has_protein = pdb_file is not None
        if has_protein and config.get("preoriented", True):
            cmd.append("--preoriented")

        # Add --parametrize to generate bilayer*_lipid.pdb with CRYST1 information
        # This file is needed for reading box dimensions during equilibration
        if config.get("parametrize", True):
            cmd.append("--parametrize")

        # Preserve PropKa / preparation protonation (GLH, ASH, HIP, …).
        # Without this, packmol-memgen's reduce can rename atoms (e.g. HA→HCA) and break tleap.
        if config.get("notprotonate", True):
            cmd.append("--notprotonate")

        # Add salt options (only for packing stage or single-stage)
        if not config.get("parametrize_only", False) and config.get("add_salt", True):
            cmd.append("--salt")
            if config.get("salt_concentration"):
                cmd.extend(["--saltcon", str(config["salt_concentration"])])
            if config.get("cation") and config["cation"] != "K+":
                cmd.extend(["--salt_c", config["cation"]])
            if config.get("anion") and config["anion"] != "Cl-":
                cmd.extend(["--salt_a", config["anion"]])

        # Add box boundary distance (include even if default value for explicit control)
        if config.get("dist") is not None:
            cmd.extend(["--dist", str(config["dist"])])

        # Add water layer distance (include even if default value for explicit control)
        if config.get("dist_wat") is not None:
            cmd.extend(["--dist_wat", str(config["dist_wat"])])

        # PACKMOL options (packing stage only)
        if not config.get("parametrize_only", False):
            if config.get("nloop") is not None:
                cmd.extend(["--nloop", str(config["nloop"])])
            if config.get("nloop_all") is not None:
                cmd.extend(["--nloop_all", str(config["nloop_all"])])
            if config.get("tolerance") is not None:
                cmd.extend(["--tolerance", str(config["tolerance"])])

        # Add explicit box dimensions (mutually exclusive with dist_wat)
        dims = config.get("dims")
        if dims and len(dims) == 3:
            cmd.extend(["--dims", str(dims[0]), str(dims[1]), str(dims[2])])

        # Fixed membrane XY size. Required by packmol-memgen when --pdb is omitted.
        # Do not infer from --dims when a protein is present — that would change
        # the existing membrane-protein command.
        explicit_distxy = config.get("distxy_fix")
        if explicit_distxy is not None and str(explicit_distxy).strip() != "":
            cmd.extend(["--distxy_fix", str(explicit_distxy)])
        elif not has_protein:
            inferred = effective_distxy_fix(config)
            if inferred is not None:
                cmd.extend(["--distxy_fix", str(inferred)])

        # Free molecules packed into water (default) or the membrane.
        solutes = normalize_builder_solutes(config)
        if solutes and not config.get("parametrize_only", False):
            for solute in solutes:
                cmd.extend(["--solute", str(solute["pdb"])])
                cmd.extend(["--solute_con", str(solute["concentration"])])
            in_mem = bool(config.get("solute_inmem", False)) or any(
                bool(s.get("in_membrane")) for s in solutes
            )
            if in_mem:
                cmd.append("--solute_inmem")
            prot_dist = config.get("solute_prot_dist")
            if prot_dist is None:
                for solute in solutes:
                    if solute.get("prot_dist") not in (None, ""):
                        prot_dist = solute["prot_dist"]
                        break
            if has_protein and prot_dist not in (None, ""):
                cmd.extend(["--solute_prot_dist", str(prot_dist)])

        # Add ligand parameters (--ligand_param frcmod:lib for each ligand)
        ligand_params = config.get("ligand_params", {})
        if ligand_params:
            ligand_args = build_ligand_param_args(ligand_params)
            cmd.extend(ligand_args)
            # Also add --gaff2 flag when ligands are present
            cmd.append("--gaff2")
            logger.info(f"Added ligand parameters for: {list(ligand_params.keys())}")

        logger.info(f"Built command: {' '.join(cmd)}")
        return cmd

    def _prepare_lipids_argument(
        self, upper_lipids: List[str], lower_lipids: List[str]
    ) -> str:
        """Prepare lipids argument for packmol-memgen."""
        lipids_parts = []

        if upper_lipids:
            lipids_parts.append(":".join(upper_lipids))

        if lower_lipids:
            if not upper_lipids:
                lipids_parts.append("")  # Empty upper leaflet
            lipids_parts.append(":".join(lower_lipids))

        return "//".join(lipids_parts) if lipids_parts else ""

    def _launch_preparation(
        self, cmd: List[str], job_dir: Path, config: Dict[str, Any]
    ) -> bool:
        """Generate preparation files and launch the process."""
        try:
            script_path = self._generate_preparation_files(cmd, job_dir, config)
            return self._run_preparation_script(job_dir, script_path, config)
        except Exception as e:
            logger.error(f"Error launching preparation: {e}", exc_info=True)
            return False

    def _get_workflow_steps(self, config: Dict[str, Any]) -> List[str]:
        """Return workflow step names for the given configuration."""
        has_protein = config.get("_has_protein", True)
        preoriented = config.get("preoriented", True) or not has_protein
        parametrize = config.get("parametrize", True)
        if preoriented:
            return ["Packmol", "pdb4amber", "tleap"] if parametrize else ["Packmol"]
        return (
            ["MEMEMBED", "Packmol", "pdb4amber", "tleap"]
            if parametrize
            else ["MEMEMBED", "Packmol"]
        )

    def _write_status_file(
        self,
        cmd: List[str],
        job_dir: Path,
        config: Dict[str, Any],
        *,
        status: str,
        start_time: Optional[str],
    ) -> None:
        """Write status.json for job monitoring."""
        status_data = {
            "command": " ".join(cmd),
            "start_time": start_time,
            "end_time": None,
            "current_step": 0,
            "status": status,
            "error": None,
            "steps": self._get_workflow_steps(config),
            "steps_completed": [],
            "step_messages": [],
            "config": config,
        }

        with open(job_dir / "status.json", "w", encoding="utf-8") as f:
            json.dump(status_data, f, indent=2)

    def _generate_preparation_files(
        self, cmd: List[str], job_dir: Path, config: Dict[str, Any]
    ) -> Path:
        """Create status.json (not_started) and run_preparation.sh."""
        self._write_status_file(
            cmd, job_dir, config, status="not_started", start_time=None
        )
        script_path = self._create_execution_script(cmd, job_dir, config)
        logger.info(f"Generated preparation input files in: {job_dir}")
        return script_path

    def _mark_preparation_running(self, job_dir: Path) -> None:
        """Update status.json to running when the job is launched."""
        status_file = job_dir / "status.json"
        if not status_file.is_file():
            return
        try:
            with open(status_file, "r", encoding="utf-8") as f:
                status_data = json.load(f)
            status_data["status"] = "running"
            status_data["start_time"] = datetime.now().isoformat()
            with open(status_file, "w", encoding="utf-8") as f:
                json.dump(status_data, f, indent=2)
        except (json.JSONDecodeError, OSError) as e:
            logger.warning(f"Could not update status.json: {e}")

    def _run_preparation_script(
        self, job_dir: Path, script_path: Path, config: Dict[str, Any]
    ) -> bool:
        """Launch the preparation script asynchronously.

        Parametrization (pdb4amber / tleap) is performed inside
        ``run_preparation.sh`` after packmol-memgen finishes.  Do not call
        :meth:`_post_process_files` here — the job is still running and bilayer
        outputs do not exist yet, which produced false WARNING noise.
        """
        try:
            self._mark_preparation_running(job_dir)
            return self._execute_script(script_path)

        except Exception as e:
            logger.error(f"Error running preparation script: {e}", exc_info=True)
            return False

    def _create_initial_status(
        self, cmd: List[str], job_dir: Path, config: Dict[str, Any]
    ):
        """Create initial status file for job monitoring (running)."""
        self._write_status_file(
            cmd,
            job_dir,
            config,
            status="running",
            start_time=datetime.now().isoformat(),
        )

    def _create_execution_script(
        self, cmd: List[str], job_dir: Path, config: Dict[str, Any]
    ) -> Path:
        """Create bash script for execution."""
        script_path = job_dir / "run_preparation.sh"
        cmd_str = " ".join(f'"{c}"' if " " in c else c for c in cmd)
        prepared_pdb_for_leap = (job_dir / "system_for_tleap.pdb").resolve()

        # Add workflow-specific notes
        workflow_note = ""
        if config.get("notprotonate", False):
            if config.get("pack_only", False):
                workflow_note = """#
    # Stage 1 of Propka Workflow: PACKING ONLY
    # This stage packs the protein into the membrane without parametrization.
    # After this completes:
    # 1. Analyze the original PDB with Propka
    # 2. Modify residue names in the packed PDB file (bilayer_*.pdb)
    # 3. Run Stage 2 parametrization with --notprotonate
    #"""
            elif config.get("parametrize_only", False):
                workflow_note = """#
    # Stage 2 of Propka Workflow: PARAMETRIZATION ONLY
    # This stage assumes residue names have been modified based on Propka results.
    # Using --notprotonate to preserve the corrected protonation states.
    #"""
        elif config.get("_workflow_note"):
            workflow_note = f"""#
    # {config['_workflow_note']}
    #"""

        from gatewizard.tools.namd_water import tleap_flexible_water_lines

        tleap_flexible_water_block = tleap_flexible_water_lines(
            config.get("md_engine"), config.get("water_model", "tip3p")
        )

        script_content = f"""#!/bin/bash
    # Run preparation script for Gatewizard
    # Generated: {datetime.now().isoformat()}
    {workflow_note}
    # New Gatewizard workflow:
    # Stage 1: Pack protein into membrane WITH --parametrize (packmol-memgen) for CRYST1 info
    # Stage 2: Run pdb4amber to prepare PDB file for Amber compatibility  
    # Stage 3: Run tleap to parametrize the system using prepared PDB
    #
    # The bilayer*_lipid.pdb file from Stage 1 provides essential CRYST1 box information
    # This workflow produces more reliable .prmtop/.inpcrd files for MD simulations
{get_clean_env_shell_snippet()}
    # Save PID for tracking
    echo $$ > "{job_dir}/process.pid"

    # Change to job directory
    cd "{job_dir}"

    # Create logs directory
    mkdir -p logs

    # Define steps for monitoring - dynamic based on configuration
    preoriented={str(config.get('preoriented', True)).lower()}
    parametrize_requested={str(config.get('parametrize', True)).lower()}
    
    if [ "$preoriented" = "true" ]; then
        if [ "$parametrize_requested" = "true" ]; then
            declare -a STEPS=("Packmol" "pdb4amber" "tleap")
        else
            declare -a STEPS=("Packmol")
        fi
    else
        if [ "$parametrize_requested" = "true" ]; then
            declare -a STEPS=("MEMEMBED" "Packmol" "pdb4amber" "tleap")
        else
            declare -a STEPS=("MEMEMBED" "Packmol")
        fi
    fi
    
    echo "INFO: Workflow steps: ${{STEPS[@]}}" | tee -a logs/preparation.log

    # Create Python status utility script
    cat > status_utils.py << 'EOF'
import json
import sys
import os
from datetime import datetime

def update_status(step, msg):
    try:
        # Debug: Print current working directory
        print('DEBUG: Working in directory:', os.getcwd())
        
        # Read current status
        try:
            with open('status.json', 'r') as f:
                status = json.load(f)
            print('DEBUG: Read existing status file')
        except Exception as e:
            print('DEBUG: Creating new status file, error was:', str(e))
            status = {{
                'status': 'running',
                'start_time': datetime.now().isoformat(),
                'current_step': 0,
                'total_steps': 5,
                'steps_completed': [],
                'last_update': None
            }}

        # Update status
        status['current_step'] = step
        status['last_update'] = datetime.now().isoformat()
        
        # Update message for current step
        if 'step_messages' not in status:
            status['step_messages'] = []
        while len(status['step_messages']) < step:
            status['step_messages'].append("")
        if step > 0:
            status['step_messages'][step-1] = msg
        
        # Write updated status
        with open('status.json', 'w') as f:
            json.dump(status, f, indent=2)
        
        print(f'DEBUG: Updated status - Step {{step}}: {{msg}}')
        
    except Exception as e:
        print('ERROR updating status:', str(e))
        import traceback
        traceback.print_exc()

def mark_complete():
    try:
        print('DEBUG: Marking job complete in', os.getcwd())
        
        with open('status.json', 'r') as f:
            status = json.load(f)

        status['status'] = 'completed'
        status['end_time'] = datetime.now().isoformat()
        status['current_step'] = max(status.get('current_step', 0), 5)  # Ensure final step
        
        # Add completion to steps if not already there
        if 'Completed' not in status.get('steps_completed', []):
            status['steps_completed'].append('Completed')

        with open('status.json', 'w') as f:
            json.dump(status, f, indent=2)
            
        print('DEBUG: Job marked as complete')
    except Exception as e:
        print('ERROR marking complete:', str(e))

def handle_error(error_msg):
    try:
        with open('status.json', 'r') as f:
            status = json.load(f)

        status['status'] = 'error'
        status['error'] = error_msg
        status['end_time'] = datetime.now().isoformat()

        with open('status.json', 'w') as f:
            json.dump(status, f, indent=2)
            
        print('ERROR: Job marked as failed -', error_msg)
    except Exception as e:
        print('ERROR handling error:', str(e))

if __name__ == "__main__":
    action = sys.argv[1]
    if action == "update":
        update_status(int(sys.argv[2]), sys.argv[3])
    elif action == "complete":
        mark_complete()
    elif action == "error":
        handle_error(sys.argv[2])
EOF

    # Function to update status
    update_status() {{
        step=$1
        msg=$2
        python3 status_utils.py update "$step" "$msg" 2>&1 | tee -a logs/status_updates.log
    }}

    # Function to mark complete
    mark_complete() {{
        python3 status_utils.py complete 2>&1 | tee -a logs/status_updates.log
    }}

    # Function to handle errors
    handle_error() {{
        error_msg=$1
        python3 status_utils.py error "$error_msg" 2>&1 | tee -a logs/status_updates.log
        exit 1
    }}

    # Monitor output for different phases
    monitor_output() {{
        line_count=0
        last_update_time=$(date +%s)
        last_step=-1
        
        while IFS= read -r line; do
            echo "$line" | tee -a logs/preparation.log
            ((line_count++))

            # Convert line to lowercase for easier pattern matching
            lower_line=$(echo "$line" | tr '[:upper:]' '[:lower:]')
            
            # Detect step progression based on actual output patterns
            step_detected=0
            step_msg=""
            
            # Dynamic step detection based on configured workflow steps
            # Step indices depend on whether MEMEMBED is included
            
            # MEMEMBED/Orientation phase (only if not pre-oriented)
            if [[ "$lower_line" == *"memembed"* ]] || [[ "$lower_line" == *"orientation"* ]] || [[ "$lower_line" == *"embed"* ]]; then
                # Only count MEMEMBED if it's in our steps array
                for i in "${{!STEPS[@]}}"; do
                    if [ "${{STEPS[i]}}" = "MEMEMBED" ]; then
                        step_detected=$((i))
                        step_msg="MEMEMBED"
                        break
                    fi
                done
            # Packmol phase  
            elif [[ "$lower_line" == *"packmol"* ]] || [[ "$lower_line" == *"packing"* ]] || [[ "$lower_line" == *"pack"* ]]; then
                for i in "${{!STEPS[@]}}"; do
                    if [ "${{STEPS[i]}}" = "Packmol" ]; then
                        step_detected=$((i))
                        step_msg="Packmol"
                        break
                    fi
                done
            # pdb4amber phase
            elif [[ "$lower_line" == *"pdb4amber"* ]] || [[ "$lower_line" == *"preparing pdb"* ]]; then
                for i in "${{!STEPS[@]}}"; do
                    if [ "${{STEPS[i]}}" = "pdb4amber" ]; then
                        step_detected=$((i))
                        step_msg="pdb4amber"
                        break
                    fi
                done
            # tleap phase
            elif [[ "$lower_line" == *"tleap"* ]] || [[ "$lower_line" == *"leap"* ]] || [[ "$lower_line" == *"amber"* ]]; then
                for i in "${{!STEPS[@]}}"; do
                    if [ "${{STEPS[i]}}" = "tleap" ]; then
                        step_detected=$((i))
                        step_msg="tleap"
                        break
                    fi
                done
            fi
            
            # Update status if we detected a new step
            if (( step_detected > last_step )); then
                update_status $step_detected "$step_msg"
                last_step=$step_detected
                last_update_time=$(date +%s)
            fi
            
            # Fallback: Update based on line count and time (every 100 lines or 45 seconds)
            current_time=$(date +%s)
            elapsed=$((current_time - last_update_time))
            
            if (( line_count % 100 == 0 || elapsed > 45 )); then
                if (( last_step < 1 && line_count > 50 )); then
                    update_status 1 "Processing"
                    last_step=1
                elif (( last_step < 2 && line_count > 200 )); then
                    update_status 2 "Building"
                    last_step=2
                elif (( last_step < 3 && line_count > 400 )); then
                    update_status 3 "Finalizing"
                    last_step=3
                fi
                last_update_time=$current_time
            fi
        done
    }}

    # Main execution
    echo "Starting membrane preparation at $(date)" | tee -a logs/preparation.log
    echo "Command: {cmd_str}" | tee -a logs/preparation.log
    
    # Initial progress update - ensure it works
    echo "Updating initial status..." | tee -a logs/preparation.log
    update_status 0 "Starting"
    
    # Brief delay to ensure status file is written
    sleep 1

    # Execute command with monitoring and capture exit status
    echo "Launching main command..." | tee -a logs/preparation.log
    {cmd_str} 2>&1 | monitor_output &
    CMD_PID=$!
    
    # Enhanced completion checker - run in background
    {{
        echo "Background completion checker started" >> logs/preparation.log
        sleep 30  # Wait 30 seconds before first check
        
        while kill -0 $CMD_PID 2>/dev/null; do
            sleep 15  # Check every 15 seconds
            
            # Check if output files exist to determine completion
            output_files_found=0
            for pattern in "bilayer_*.pdb" "*.pdb" "*.prmtop" "*.inpcrd"; do
                if compgen -G "$pattern" > /dev/null 2>&1; then
                    output_files_found=1
                    echo "Found output files matching $pattern - job progressing" >> logs/preparation.log
                    break
                fi
            done
            
            if [ $output_files_found -eq 1 ]; then
                # Update to Packmol step when files are first detected (still in packing phase)
                update_status 2 "Packmol"
            fi
        done
        echo "Background completion checker finished" >> logs/preparation.log
    }} &
    COMPLETION_CHECKER_PID=$!
    
    # Wait for main command to complete
    wait $CMD_PID
    EXIT_STATUS=$?

    # Clean up background checker
    if kill -0 $COMPLETION_CHECKER_PID 2>/dev/null; then
        kill $COMPLETION_CHECKER_PID 2>/dev/null || true
    fi

    # Check exit status and handle completion
    echo "Command completed with exit status $EXIT_STATUS" | tee -a logs/preparation.log
    
    if [ $EXIT_STATUS -eq 0 ]; then
        echo "Membrane preparation completed successfully at $(date)" | tee -a logs/preparation.log
        
        # Check if parametrization was requested
        parametrize_requested={str(config.get('parametrize', True)).lower()}
        
        if [ "$parametrize_requested" = "true" ]; then
            echo "Starting parametrization workflow..." | tee -a logs/preparation.log
            
            # Find the bilayer PDB file generated by packmol-memgen
            bilayer_pdb=""
            for pdb_file in bilayer_*.pdb; do
                if [ -f "$pdb_file" ]; then
                    bilayer_pdb="$pdb_file"
                    echo "Found bilayer PDB file: $bilayer_pdb" | tee -a logs/preparation.log
                    break
                fi
            done
            
            if [ -z "$bilayer_pdb" ]; then
                echo "❌ No bilayer PDB file found for parametrization" | tee -a logs/preparation.log
                handle_error "Bilayer PDB file not found"
            fi
            
            # Step 1: Run pdb4amber
            echo "Running pdb4amber on $bilayer_pdb..." | tee -a logs/preparation.log
            
            # Find the correct step index for pdb4amber in the dynamic steps array
            pdb4amber_step_index=-1
            for i in "${{!STEPS[@]}}"; do
                if [ "${{STEPS[i]}}" = "pdb4amber" ]; then
                    pdb4amber_step_index=$i
                    break
                fi
            done
            
            if [ $pdb4amber_step_index -ge 0 ]; then
                update_status $pdb4amber_step_index "pdb4amber"
            else
                echo "Warning: pdb4amber step not found in steps array" | tee -a logs/preparation.log
            fi
            
            prepared_pdb="system_for_tleap.pdb"
            if pdb4amber -i "$bilayer_pdb" -o "$prepared_pdb" 2>&1 | tee -a logs/parametrization.log; then
                echo "✅ pdb4amber completed successfully" | tee -a logs/preparation.log
                
                # Step 2: Run tleap
                echo "Running tleap parametrization..." | tee -a logs/preparation.log
                
                # Find the correct step index for tleap in the dynamic steps array
                tleap_step_index=-1
                for i in "${{!STEPS[@]}}"; do
                    if [ "${{STEPS[i]}}" = "tleap" ]; then
                        tleap_step_index=$i
                        break
                    fi
                done
                
                if [ $tleap_step_index -ge 0 ]; then
                    update_status $tleap_step_index "tleap"
                else
                    echo "Warning: tleap step not found in steps array" | tee -a logs/preparation.log
                fi
                
                # Get force field configuration
                protein_ff="{config.get('protein_ff', 'ff14SB')}"
                lipid_ff="{config.get('lipid_ff', 'lipid21')}"
                water_model="{config.get('water_model', 'tip3p')}"
                
                # Map force fields to leaprc files
                case "$protein_ff" in
                    "ff14SB") protein_leaprc="leaprc.protein.ff14SB" ;;
                    "ff19SB") protein_leaprc="leaprc.protein.ff19SB" ;;
                    "ff99SB") protein_leaprc="leaprc.protein.ff99SB" ;;
                    "ff03") protein_leaprc="leaprc.protein.ff03" ;;
                    *) protein_leaprc="leaprc.protein.ff14SB" ;;
                esac
                
                case "$lipid_ff" in
                    "lipid21") lipid_leaprc="leaprc.lipid21" ;;
                    "lipid17") lipid_leaprc="leaprc.lipid17" ;;
                    "GAFF") lipid_leaprc="leaprc.gaff" ;;
                    *) lipid_leaprc="leaprc.lipid21" ;;
                esac
                
                case "$water_model" in
                    "tip3p") water_leaprc="leaprc.water.tip3p" ;;
                    "tip4pd") water_leaprc="leaprc.water.tip4pd" ;;
                    "tip4pew") water_leaprc="leaprc.water.tip4pew" ;;
                    "spce") water_leaprc="leaprc.water.spce" ;;
                    "opc") water_leaprc="leaprc.water.opc" ;;
                    "opc3") water_leaprc="leaprc.water.opc3" ;;
                    "spceb") water_leaprc="leaprc.water.spceb" ;;
                    "fb3") water_leaprc="leaprc.water.fb3" ;;
                    *) water_leaprc="leaprc.water.tip3p" ;;
                esac
                
                echo "Force field configuration:" | tee -a logs/preparation.log
                echo "  Protein FF: $protein_ff -> $protein_leaprc" | tee -a logs/preparation.log
                echo "  Lipid FF: $lipid_ff -> $lipid_leaprc" | tee -a logs/preparation.log
                echo "  Water model: $water_model -> $water_leaprc" | tee -a logs/preparation.log
                echo "  MD engine target: {config.get('md_engine') or '(not set)'}" | tee -a logs/preparation.log
                
                # Create dynamic tleap input file
                cat > leap_parametrize.in << EOF
# Load force field for proteins $protein_ff
source $protein_leaprc

# Load force field for lipids $lipid_ff
source $lipid_leaprc

# Load water model {config.get('water_model', 'tip3p').upper()}
source $water_leaprc

{tleap_flexible_water_block}
{self._generate_bash_ligand_tleap_lines(config)}

# Load PDB file prepared by pdb4amber (protein + membrane + water + neutralized)
system = loadPDB "{prepared_pdb_for_leap}"

# Check total system charge
charge system

# Neutralize total charge
addIonsRand system Na+ 0
addIonsRand system Cl- 0

# Save parameter and coordinate files
saveAmberParm system system.prmtop system.inpcrd

# Save the system processed by tleap as PDB
savePDB system system.pdb

# Exit
quit
EOF
                
                # Execute tleap
                if tleap -f leap_parametrize.in 2>&1 | tee -a logs/parametrization.log; then
                    echo "✅ tleap parametrization completed successfully" | tee -a logs/preparation.log
                    
                    # Check if output files were created
                    output_files=("system.prmtop" "system.inpcrd" "system.pdb")
                    all_created=true
                    
                    echo "Generated files:" | tee -a logs/preparation.log
                    for filename in "${{output_files[@]}}"; do
                        if [ -f "$filename" ]; then
                            size=$(stat -c%s "$filename" 2>/dev/null || stat -f%z "$filename" 2>/dev/null || echo "unknown")
                            echo "   $filename: $size bytes" | tee -a logs/preparation.log
                        else
                            echo "   ❌ $filename: NOT CREATED" | tee -a logs/preparation.log
                            all_created=false
                        fi
                    done
                    
                    if [ "$all_created" = true ]; then
                        echo "✅ All parametrization files created successfully" | tee -a logs/preparation.log
                    else
                        echo "❌ Some parametrization files were not created" | tee -a logs/preparation.log
                        handle_error "Incomplete parametrization output"
                    fi
                    
                else
                    echo "❌ tleap parametrization failed" | tee -a logs/preparation.log
                    handle_error "tleap parametrization failed"
                fi
                
            else
                echo "❌ pdb4amber failed" | tee -a logs/preparation.log
                handle_error "pdb4amber failed"
            fi
            
        else
            echo "Skipping parametrization - not requested" | tee -a logs/preparation.log
        fi
        
        # Final status update
        update_status 5 "Completed"
        
        # Mark as completed 
        mark_complete
        
    else
        echo "Command failed with exit status $EXIT_STATUS" | tee -a logs/preparation.log
        
        # Check if output files exist anyway (sometimes commands exit with error but produce output)
        output_found=0
        for pattern in "bilayer_*.pdb" "*.pdb" "*.prmtop" "*.inpcrd"; do
            if compgen -G "$pattern" > /dev/null 2>&1; then
                output_found=1
                echo "Found output files despite exit status - marking as completed" | tee -a logs/preparation.log
                update_status 5 "Completed"
                mark_complete
                break
            fi
        done
        
        if [ $output_found -eq 0 ]; then
            echo "No output files found - marking as failed" | tee -a logs/preparation.log
            handle_error "Preparation failed with exit code $EXIT_STATUS"
        fi
    fi
    """

        # Write and make executable
        with open(script_path, "w", encoding="utf-8", newline="\n") as f:
            f.write(script_content)
        script_path.chmod(0o755)

        logger.info(f"Created execution script: {script_path}")
        return script_path

    def _execute_script(self, script_path: Path) -> bool:
        """Execute the preparation script."""
        try:
            # Check if we're in WSL environment
            is_wsl = (
                os.path.exists("/proc/version")
                and "microsoft" in open("/proc/version").read().lower()
            )

            if sys.platform.startswith("win") and not is_wsl:
                # Native Windows execution - use Windows-specific constants
                DETACHED_PROCESS = 0x00000008
                CREATE_NEW_PROCESS_GROUP = 0x00000200

                git_bash = r"C:\Program Files\Git\bin\bash.exe"
                if os.path.exists(git_bash):
                    subprocess.Popen(
                        [git_bash, str(script_path)],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                        env=get_clean_env(),
                        creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
                    )
                else:
                    # Fallback to WSL from Windows
                    subprocess.Popen(
                        ["wsl", "bash", str(script_path)],
                        stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL,
                        env=get_clean_env(),
                        creationflags=DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP,
                    )
            else:
                # Unix-like systems (including WSL)
                subprocess.Popen(
                    ["bash", str(script_path)],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    env=get_clean_env(),
                    start_new_session=True,
                )

            logger.info("Preparation script launched successfully")
            return True

        except Exception as e:
            logger.error(f"Error executing script: {e}")
            return False

    def _post_process_files(self, job_dir: Path, config: Dict[str, Any]) -> bool:
        """Post-process generated files with pdb4amber / tleap (legacy helper).

        The generated ``run_preparation.sh`` already runs this workflow after
        packmol-memgen completes.  Prefer that path for async jobs; this method
        remains for synchronous/manual use when bilayer outputs already exist.
        """
        try:
            # Check if parametrization was requested
            if not config.get("parametrize", True):
                logger.info("Skipping parametrization - parametrize was False")
                return True

            # Find the generated bilayer PDB file (from packmol-memgen with --parametrize)
            bilayer_pdb_files = list(job_dir.glob("bilayer_*.pdb"))

            if not bilayer_pdb_files:
                logger.debug(
                    f"No bilayer PDB files found yet in {job_dir} "
                    "(expected until packmol-memgen finishes)"
                )
                return False

            bilayer_pdb = bilayer_pdb_files[0]  # Use first found
            logger.info(f"Found bilayer PDB file: {bilayer_pdb.name}")

            # Change to job directory for processing
            original_cwd = Path.cwd()
            os.chdir(job_dir)

            try:
                # Step 1: Run pdb4amber to prepare the PDB file
                prepared_pdb = self._run_pdb4amber(bilayer_pdb.name, config)
                if not prepared_pdb:
                    return False

                # Step 2: Run tleap parametrization
                success = self._run_tleap_parametrization(prepared_pdb, config)
                if success:
                    logger.info("Parametrization completed successfully")
                    return True
                else:
                    logger.error("Parametrization failed")
                    return False

            except Exception as e:
                logger.error(f"Error during parametrization: {e}")
                return False
            finally:
                os.chdir(original_cwd)

        except Exception as e:
            logger.error(f"Error in post-processing: {e}")
            return False

    def _run_pdb4amber(self, input_pdb: str, config: Dict[str, Any]) -> Optional[str]:
        """Run pdb4amber to prepare PDB file for Amber parametrization."""
        try:
            output_pdb = "system_for_tleap.pdb"
            logger.info(f"Running pdb4amber on {input_pdb}")

            # Build pdb4amber command
            pdb4amber_exe = _resolve_pdb4amber_executable()
            cmd = subprocess_argv_for_script(
                pdb4amber_exe, ["-i", input_pdb, "-o", output_pdb]
            )

            # Execute pdb4amber
            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                env=get_clean_env(),
            )

            if Path(output_pdb).exists():
                logger.info(f"pdb4amber completed successfully: {output_pdb}")
                logger.debug(f"pdb4amber output: {result.stdout}")
                return output_pdb
            else:
                logger.error("pdb4amber did not create output file")
                return None

        except subprocess.CalledProcessError as e:
            logger.error(f"pdb4amber failed: {e.stderr if e.stderr else str(e)}")
            return None
        except FileNotFoundError:
            logger.error("pdb4amber not found. Please install AmberTools.")
            return None
        except Exception as e:
            logger.error(f"Error running pdb4amber: {e}")
            return None

    def _run_tleap_parametrization(
        self, input_pdb: str, config: Dict[str, Any]
    ) -> bool:
        """Run tleap parametrization using the new workflow."""
        try:
            # Create tleap input file
            leap_input = self._create_tleap_input(input_pdb, config)

            # Write leap input file
            leap_input_file = "leap_parametrize.in"
            with open(leap_input_file, "w") as f:
                f.write(leap_input)

            logger.info(f"Created tleap input file: {leap_input_file}")

            # Execute tleap
            tleap_exe = resolve_conda_executable("tleap")
            cmd = subprocess_argv_for_script(tleap_exe, ["-f", leap_input_file])

            result = subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
                stdin=subprocess.DEVNULL,
                env=get_clean_env(),
            )

            # Check if output files were created
            expected_files = ["system.prmtop", "system.inpcrd", "system.pdb"]
            success = True

            for filename in expected_files:
                if Path(filename).exists():
                    logger.info(f"Created: {filename}")
                else:
                    logger.warning(f"Expected file not created: {filename}")
                    success = False

            if success:
                logger.info("tleap parametrization completed successfully")
                logger.debug(f"tleap output: {result.stdout}")

            return success

        except subprocess.CalledProcessError as e:
            logger.error(f"tleap failed: {e.stderr if e.stderr else str(e)}")
            return False
        except FileNotFoundError:
            logger.error("tleap not found. Please install AmberTools.")
            return False
        except Exception as e:
            logger.error(f"Error running tleap: {e}")
            return False

    def _create_tleap_input(self, input_pdb: str, config: Dict[str, Any]) -> str:
        """Create tleap input file content."""
        # Get absolute path for PDB file
        pdb_path = Path(input_pdb).absolute()

        # Get force field parameters from config
        protein_ff = config.get("protein_ff", "ff14SB")
        lipid_ff = config.get("lipid_ff", "lipid21")
        water_model = config.get("water_model", "tip3p")

        # Map protein force fields to leaprc files
        protein_leaprc_map = {
            "ff14SB": "leaprc.protein.ff14SB",
            "ff19SB": "leaprc.protein.ff19SB",
            "ff99SB": "leaprc.protein.ff99SB",
            "ff03": "leaprc.protein.ff03",
        }

        # Map lipid force fields to leaprc files
        lipid_leaprc_map = {
            "lipid21": "leaprc.lipid21",
            "lipid17": "leaprc.lipid17",
            "GAFF": "leaprc.gaff",
        }

        # Map water models to leaprc files
        water_leaprc_map = {
            "tip3p": "leaprc.water.tip3p",
            "tip4pd": "leaprc.water.tip4pd",
            "tip4pew": "leaprc.water.tip4pew",
            "spce": "leaprc.water.spce",
            "opc": "leaprc.water.opc",
            "opc3": "leaprc.water.opc3",
            "spceb": "leaprc.water.spceb",
            "fb3": "leaprc.water.fb3",
        }

        # Get the appropriate leaprc files
        protein_leaprc = protein_leaprc_map.get(protein_ff, "leaprc.protein.ff14SB")
        lipid_leaprc = lipid_leaprc_map.get(lipid_ff, "leaprc.lipid21")
        water_leaprc = water_leaprc_map.get(water_model, "leaprc.water.tip3p")

        md_engine = config.get("md_engine")
        from gatewizard.tools.namd_water import tleap_flexible_water_lines

        flexible_water = tleap_flexible_water_lines(md_engine, water_model)

        # Generate ligand parameter lines if ligands are present
        ligand_params = config.get("ligand_params", {})
        # Extract atom type from parametrized ligand info (all ligands share the same type)
        ligand_atom_type = "gaff2"
        for _lig_info in ligand_params.values():
            if isinstance(_lig_info, dict) and "atom_type" in _lig_info:
                ligand_atom_type = _lig_info["atom_type"]
                break
        ligand_lines = build_tleap_ligand_lines(
            ligand_params, atom_type=ligand_atom_type
        )

        # Generate tleap input content
        leap_content = f"""# Load force field for proteins {protein_ff}
source {protein_leaprc}

# Load force field for lipids {lipid_ff}
source {lipid_leaprc}

# Load water model {water_model.upper()}
source {water_leaprc}

{flexible_water}
{ligand_lines}

# Load PDB file (protein + membrane + water + neutralized)
system = loadPDB {pdb_path}

# Check total system charge
charge system

# Neutralize total charge
addIonsRand system Na+ 0
addIonsRand system Cl- 0

# Save parameter and coordinate files
saveAmberParm system system.prmtop system.inpcrd

# Save the system processed by tleap as PDB
savePDB system system.pdb

# Exit
quit
"""

        return leap_content

    def _generate_bash_ligand_tleap_lines(self, config: Dict[str, Any]) -> str:
        """Generate tleap ligand parameter lines for the bash execution script."""
        ligand_params = config.get("ligand_params", {})
        if not ligand_params:
            return ""
        lines = ["# Load GAFF2 and ligand parameters", "source leaprc.gaff2"]
        for name, files in ligand_params.items():
            frcmod = files.get("frcmod", "")
            lib = files.get("lib", "")
            if frcmod:
                lines.append(f"loadamberparams {frcmod}")
            if lib:
                lines.append(f"loadoff {lib}")
        return "\n".join(lines)

    def prepare_system_stage1_for_propka(
        self,
        pdb_file: str,
        working_dir: str,
        upper_lipids: List[str],
        lower_lipids: List[str],
        lipid_ratios: str = "",
        **kwargs,
    ) -> Tuple[bool, str, Optional[Path]]:
        """
        Prepare system - Stage 1 for Propka workflow: Packing only.

        This method implements the first stage of the two-stage workflow
        recommended in the packmol-memgen tutorial for Propka integration:
        1. Pack protein into membrane without parametrization
        2. Apply Propka protonation states to packed file
        3. Run parametrization on modified packed file

        Args:
            pdb_file: Path to input PDB file
            working_dir: Working directory for outputs
            upper_lipids: List of lipids for upper leaflet
            lower_lipids: List of lipids for lower leaflet
            lipid_ratios: Lipid ratios string
            **kwargs: Additional configuration parameters

        Returns:
            Tuple of (success, message, job_directory)
        """
        # Configure for stage 1: pack only, let packmol-memgen handle initial protonation
        stage1_config = {
            **self.config,
            **kwargs,
            "pack_only": True,
            "parametrize_only": False,
            "parametrize": False,  # Disable parametrization for stage 1
            "notprotonate": False,  # Let packmol-memgen handle initial protonation during packing
            "_workflow_note": (
                "Stage 1 of Propka workflow: Packing protein into membrane. "
                "After this completes, modify residue names in the packed PDB "
                "based on Propka results, then run Stage 2 parametrization."
            ),
        }

        return self.prepare_system(
            pdb_file,
            working_dir,
            upper_lipids,
            lower_lipids,
            lipid_ratios,
            **stage1_config,
        )

    def prepare_system_stage2_for_propka(
        self, packed_pdb_file: str, working_dir: str, **kwargs
    ) -> Tuple[bool, str, Optional[Path]]:
        """
        Prepare system - Stage 2 for Propka workflow: Parametrization only.

        This method implements the second stage of the two-stage workflow.
        It assumes the packed PDB file has been modified with correct
        residue names based on Propka analysis.

        Args:
            packed_pdb_file: Path to packed PDB file with correct residue names
            working_dir: Working directory for outputs
            **kwargs: Additional configuration parameters

        Returns:
            Tuple of (success, message, job_directory)
        """
        # Configure for stage 2: parametrize only with notprotonate
        stage2_config = {
            **self.config,
            **kwargs,
            "pack_only": False,
            "parametrize_only": True,
            "parametrize": True,  # Enable parametrization for stage 2
            "notprotonate": True,  # Skip protonation since residue names are already correct
            "_workflow_note": (
                "Stage 2 of Propka workflow: Parametrizing system with "
                "Propka-modified residue names using --notprotonate."
            ),
        }

        # For stage 2, we don't need to specify lipids since the system is already packed
        return self.prepare_system(
            packed_pdb_file, working_dir, [], [], "", **stage2_config
        )

    def modify_residue_names_for_propka(
        self,
        packed_pdb_file: str,
        propka_results: Dict[str, str],
        output_file: Optional[str] = None,
    ) -> Tuple[bool, str]:
        """
        Modify residue names in packed PDB file based on propka results.

        This function implements the workflow described in the packmol-memgen tutorial
        where residue names are changed according to propka protonation states
        between the packing and parametrization stages.

        Example from tutorial:
        for i in " 49" 146 243 340; do sed -i "/GLU . $i/s/GLU/GLH/g" bilayer_1BL8.pdb; done

        Args:
            packed_pdb_file: Path to packed PDB file
            propka_results: Dictionary mapping residue IDs to new residue names
                          e.g., {"A:71": "GLH", "B:71": "GLH"} or {"49": "GLH", "146": "GLH"}
            output_file: Output file path (default: same as input with _modified suffix)

        Returns:
            Tuple of (success, message)
        """
        try:
            if output_file is None:
                # Use basename with _propka suffix for Propka workflow
                base_name = os.path.splitext(packed_pdb_file)[0]
                output_file = f"{base_name}_propka.pdb"

            with open(packed_pdb_file, "r") as f:
                lines = f.readlines()

            modified_lines = []
            modifications_made = []

            for line in lines:
                if line.startswith(("ATOM", "HETATM")):
                    # Parse PDB line
                    chain = line[21:22].strip()
                    res_num = line[22:26].strip()
                    res_name = line[17:20].strip()

                    # Create residue identifier - try both chain:resnum and resnum formats
                    res_id_with_chain = f"{chain}:{res_num}" if chain else res_num
                    res_id_no_chain = res_num

                    # Check if this residue should be modified
                    new_res_name = None
                    matched_id = None

                    if res_id_with_chain in propka_results:
                        new_res_name = propka_results[res_id_with_chain]
                        matched_id = res_id_with_chain
                    elif res_id_no_chain in propka_results:
                        new_res_name = propka_results[res_id_no_chain]
                        matched_id = res_id_no_chain

                    if new_res_name and new_res_name != res_name:
                        # Replace residue name in the line
                        modified_line = line[:17] + f"{new_res_name:>3}" + line[20:]
                        modified_lines.append(modified_line)

                        # Track modifications (avoid duplicates)
                        mod_key = (matched_id, res_name, new_res_name)
                        if mod_key not in modifications_made:
                            modifications_made.append(mod_key)
                    else:
                        modified_lines.append(line)
                else:
                    modified_lines.append(line)

            # Write modified file
            with open(output_file, "w") as f:
                f.writelines(modified_lines)

            if modifications_made:
                mod_summary = "; ".join(
                    [f"{rid}: {old}->{new}" for rid, old, new in modifications_made]
                )
                message = (
                    f"Successfully modified {len(modifications_made)} residue types: {mod_summary}\n"
                    f"Output file: {output_file}\n"
                    f"Ready for Stage 2 parametrization with --notprotonate"
                )
                logger.info(message)
                return True, message
            else:
                message = (
                    f"No residues were modified (no matches found).\n"
                    f"File copied to: {output_file}\n"
                    f"Proceed with Stage 2 parametrization."
                )
                logger.warning(message)
                return True, message

        except Exception as e:
            error_msg = f"Error modifying residue names: {str(e)}"
            logger.error(error_msg, exc_info=True)
            return False, error_msg


def convert_to_namd(
    top_file: str,
    crd_file: str,
    psf_output: str = "system.psf",
    pdb_output: str = "system.pdb",
) -> None:
    """
    Convert Amber files to NAMD format using ParmEd.

    Args:
        top_file: Amber topology file (.top)
        crd_file: Amber coordinate file (.crd)
        psf_output: Output PSF file name
        pdb_output: Output PDB file name

    Raises:
        BuilderError: If conversion fails
    """
    try:
        import parmed as pmd

        logger.info(f"Converting {top_file} and {crd_file} to NAMD format")

        # Load Amber files
        system = pmd.load_file(top_file, crd_file)

        # Save in NAMD format
        system.save(psf_output, overwrite=True)
        system.save(pdb_output, overwrite=True)

        logger.info(f"Files converted and saved as {psf_output} and {pdb_output}")

    except ImportError:
        raise BuilderError(
            "ParmEd is required for NAMD conversion. "
            "Install with: pip install parmed"
        )
    except Exception as e:
        raise BuilderError(f"Error converting to NAMD format: {e}")


def convert_to_amber(
    top_file: str,
    crd_file: str,
    prmtop_output: str = "system.prmtop",
    inpcrd_output: str = "system.inpcrd",
) -> None:
    """
    Convert Amber top/crd files to prmtop/inpcrd format using ParmEd.

    Args:
        top_file: Amber topology file (.top)
        crd_file: Amber coordinate file (.crd)
        prmtop_output: Output parameter topology file name (.prmtop)
        inpcrd_output: Output restart file name (.inpcrd)

    Raises:
        BuilderError: If conversion fails
    """
    try:
        import parmed as pmd

        logger.info(
            f"Converting {top_file} and {crd_file} to Amber prmtop/inpcrd format"
        )

        # Load Amber files
        system = pmd.load_file(top_file, crd_file)

        # Save in Amber prmtop/inpcrd format
        system.save(prmtop_output, overwrite=True)
        system.save(inpcrd_output, overwrite=True)

        logger.info(f"Files converted and saved as {prmtop_output} and {inpcrd_output}")

    except ImportError:
        raise BuilderError(
            "ParmEd is required for Amber conversion. "
            "Install with: pip install parmed"
        )
    except Exception as e:
        raise BuilderError(f"Error converting to Amber format: {e}")
