# gatewizard/tools/equilibration.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Equilibration tools for molecular dynamics simulations.

This module provides tools for generating equilibration protocols and
input files for various molecular dynamics engines.
"""

import os
import re
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple, Set
from dataclasses import dataclass, field, replace as _dc_replace
from collections import defaultdict
import json
import tempfile

from gatewizard.utils.logger import get_logger
from gatewizard.utils.equilibration_resume import (
    GROMACS_RESUME_SHELL,
    NAMD_RESUME_SHELL,
    OPENMM_RESUME_SHELL,
)
from gatewizard.utils.equilibration_resources import resolve_compute_resources_from_stages
from gatewizard.tools.namd_water import namd_water_model_config_block

logger = get_logger(__name__)


def _gromacs_mdrun_resource_flags(
    cpu_cores: Optional[int] = None,
    use_gpu: bool = False,
    gpu_id: int = 0,
    num_gpus: int = 1,
) -> str:
    """Build ``mdrun`` CPU/GPU flags from equilibration resource settings.

    GROMACS 2026 requires ``-ntmpi`` whenever ``-ntomp`` is set on a
    GPU-capable build, even if GPU offload flags are omitted, because GPUs
    may still be auto-detected. Mapping for GPU runs: one thread-MPI rank per
    GPU; OpenMP threads are UI processors divided across those ranks.

    CPU-only runs (including energy minimisation) force ``-nb cpu -pme cpu``
    so a CUDA ``gmx`` does not attempt PME-on-GPU.
    """
    parts: List[str] = []
    cores = int(cpu_cores or 0)
    if use_gpu:
        ngpu = max(1, int(num_gpus or 1))
        ntmpi = ngpu
        parts.append(f"-ntmpi {ntmpi}")
        if cores > 0:
            ntomp = max(1, cores // ntmpi)
            parts.append(f"-ntomp {ntomp}")
        parts.extend(["-nb gpu", "-pme gpu"])
        gid = int(gpu_id or 0)
        # GROMACS -gpu_id is a compact digit string, e.g. "0" or "01".
        parts.append(f"-gpu_id {''.join(str(gid + i) for i in range(ngpu))}")
    else:
        # Always pair -ntomp with -ntmpi on GPU-capable builds (EM / CPU target).
        parts.append("-ntmpi 1")
        if cores > 0:
            parts.append(f"-ntomp {cores}")
        parts.extend(["-nb cpu", "-pme cpu"])
    return " ".join(parts)


def _gromacs_gpu_info(use_gpu: bool, gpu_id: int, num_gpus: int) -> str:
    if not use_gpu:
        return "No"
    ngpu = max(1, int(num_gpus or 1))
    gid = int(gpu_id or 0)
    if ngpu == 1:
        return f"Yes (GPU {gid})"
    gpu_list = ",".join(str(gid + i) for i in range(ngpu))
    return f"Yes ({ngpu} GPUs: {gpu_list})"


# ---------------------------------------------------------------------------
# Module-level helpers shared by NAMD / GROMACS COM colvars generation
# ---------------------------------------------------------------------------


def _build_com_colvars_config(
    atom_numbers: str,
    x0: float,
    y0: float,
    z0: float,
    com_k: float,
    add_rotation: bool,
    rot_k: float,
    ag: Any,
    engine: str = "namd",
    ref_positions_file: Optional[str] = None,
    rotation_ref_positions_mode: str = "auto",
    ref_positions_col: Optional[str] = None,
    ref_positions_col_value: Optional[float] = None,
) -> str:
    """Return Colvars configuration text for translation (+ optional rotation).

    Generates a ``center`` CV (``distance`` to a dummy atom at the initial
    centroid) to restrain translation, and optionally a ``rotation`` CV
    (``orientation``) to restrain rigid-body rotation.

    Args:
        atom_numbers: Space-separated 1-based atom serial numbers.
        x0, y0, z0: Initial centroid in Ångströms.
        com_k: Translation force constant in kcal/mol/Å².
        add_rotation: Also add a ``rotation``/``orientation`` CV.
        rot_k: Rotation force constant in kcal/mol/Å².
        ag: MDAnalysis AtomGroup (used to write ``refPositions`` inline when
            the inline mode is active).
        engine: ``"namd"`` or ``"gromacs"`` (only affects header comment style).
        ref_positions_file: Optional ``refPositionsFile`` path for the
            orientation CV.  For NAMD, a PDB file is accepted.  For GROMACS,
            **only XYZ format is supported** (GROMACS Colvars manual §3.7.3);
            ``GROMACSEquilibrationManager.generate_com_colvars_config``
            auto-generates a ``.xyz`` file when this argument is not given.
        rotation_ref_positions_mode: How to encode orientation reference
            coordinates: ``"auto"`` (default, uses ``refPositionsFile``),
            ``"refPositions"`` (inline coordinates), or
            ``"refPositionsFile"`` (external file).
        ref_positions_col: PDB column flag for selecting atoms from
            ``refPositionsFile`` (``O``, ``B``, ``X``, ``Y``, or ``Z``).
            GROMACS XYZ files do not support this option.
        ref_positions_col_value: Optional numeric value paired with
            ``ref_positions_col``.

    Returns:
        Multi-line Colvars config string.
    """
    header = "# Colvars - COM restraint\n# Generated by GateWizard\n"
    comment = "#"

    lines = [
        header,
        f"{comment} Translation restraint on the geometric center of {len(ag.atoms)} atoms.",
        f"{comment} Force constant: {com_k} kcal/mol/A^2 (same units as positional restraints).",
        f"{comment} To disable, comment out the 'harmonic' block.",
        "",
    ]

    atom_numbers_list = atom_numbers.split()

    def _atom_number_lines(indent: str = "                ") -> List[str]:
        rows: List[str] = []
        for i in range(0, len(atom_numbers_list), 20):
            rows.append(indent + " ".join(atom_numbers_list[i : i + 20]))
        return rows

    for label, axis_vec in [("x", "(1, 0, 0)"), ("y", "(0, 1, 0)"), ("z", "(0, 0, 1)")]:
        lines += [
            "colvar {",
            f"   name center_{label}",
            "   outputValue",
            "",
            "   distanceZ {",
            f"      axis {axis_vec}",
            "      main {",
            "         atomNumbers {",
        ]
        lines += _atom_number_lines(indent="            ")
        lines += [
            "         }",
            "      }",
            "",
            "      ref {",
            "         dummyAtom {",
            f"            ({x0:.6f}, {y0:.6f}, {z0:.6f})",
            "         }",
            "      }",
            "   }",
            "}",
            "",
            "harmonic {",
            f"   name restraint_{label}",
            f"   colvars center_{label}",
            "   centers 0.0",
            f"   forceConstant {com_k:.4f}",
            "}",
            "",
        ]

    if add_rotation:
        mode_norm = str(rotation_ref_positions_mode).strip().lower()
        if mode_norm == "auto":
            use_ref_positions_file = True
        elif mode_norm in {"refpositions", "inline"}:
            use_ref_positions_file = False
        elif mode_norm in {"refpositionsfile", "file"}:
            use_ref_positions_file = True
        else:
            raise ValueError(
                "Invalid rotation_ref_positions_mode. "
                "Use 'auto', 'refPositions', or 'refPositionsFile'."
            )

        if use_ref_positions_file and not ref_positions_file:
            raise ValueError(
                "Rotation restraint with refPositionsFile mode requires "
                "ref_positions_file."
            )
        if ref_positions_col_value is not None and not ref_positions_col:
            raise ValueError("ref_positions_col_value requires ref_positions_col.")

        ref_block = ""
        if not use_ref_positions_file:
            ref_lines = []
            for atom in ag.atoms:
                pos = atom.position
                ref_lines.append(
                    f"            ({pos[0]:.4f}, {pos[1]:.4f}, {pos[2]:.4f})"
                )
            ref_block = "\n".join(ref_lines)
        lines += [
            f"{comment} Rotation restraint (identity quaternion = no rotation).",
            "colvar {",
            "   name rotation",
            "   outputValue",
            "",
            "   orientation {",
        ]
        lines += ["      atoms {", "         atomNumbers {"]
        lines += _atom_number_lines()
        lines += ["         }", "      }"]
        if use_ref_positions_file:
            lines += [f"      refPositionsFile {ref_positions_file}"]
            if ref_positions_col:
                lines += [f"      refPositionsCol {ref_positions_col}"]
            if ref_positions_col_value is not None:
                lines += [f"      refPositionsColValue {ref_positions_col_value:g}"]
        else:
            lines += [
                "      refPositions {",
                ref_block,
                "      }",
            ]
        lines += [
            "   }",
            "}",
            "",
            "harmonic {",
            "   name rotation",
            "   colvars rotation",
            "   centers (1.0, 0.0, 0.0, 0.0)",
            f"   forceConstant {rot_k:.4f}",
            "}",
            "",
        ]

    return "\n".join(lines)


def _build_com_colvars_activation_block(engine: str, config_filename: str) -> str:
    """Return the input-block snippet that activates COM colvars."""
    if engine == "gromacs":
        return (
            "\n; Colvars activation\n"
            "colvars-active         = yes\n"
            f"colvars-configfile     = {config_filename}\n"
        )
    return (
        "\n# Colvars activation\n" "colvars on\n" f"colvarsConfig {config_filename}\n"
    )


@dataclass
class EquilibrationStage:
    """
    Parameters for a single equilibration stage.

    Accepted by both :class:`OpenMMEquilibrationManager` and
    :class:`NAMDEquilibrationManager` as an alternative to plain dicts.
    Fields prefixed with the engine name are silently ignored by the other engine.

    Args:
        name: Human-readable stage label.
        ensemble: NVT | NPT | NPAT | NPgT.
        time_ns: Simulation length in nanoseconds.
        timestep: Integration timestep in femtoseconds.
        temperature: Target temperature in Kelvin.
        constraints: Force constants (kcal/mol/Å²) keyed by atom class:
            ``protein_backbone``, ``protein_sidechain``, ``lipid_head``,
            ``lipid_tail``, ``water``, ``ions``, ``other``.
        minimize_steps: Energy minimisation steps (first stage only).
        dcd_freq: Trajectory write frequency in steps (OpenMM).
        steps: Explicit step count override (NAMD); computed from
            ``time_ns / timestep`` if None.
        pressure: Target pressure in bar (NAMD, NPT/NPAT/NPgT ensembles).
        surface_tension: Surface tension in dyn/cm (NAMD, NPAT/NPgT ensembles).

    Example::

        from dataclasses import replace
        from gatewizard.tools.equilibration import (
            OpenMMEquilibrationManager, EquilibrationStage
        )

        stages = OpenMMEquilibrationManager.get_default_stage_params("NPT",
                                                                      include_production=True)

        # Attribute access
        stages[-1].time_ns = 100.0

        # Immutable-style copy with dataclasses.replace
        stages[-1] = replace(stages[-1], time_ns=100.0, temperature=303.15)

        manager = OpenMMEquilibrationManager(prepared_folder)
        result = manager.setup_openmm_equilibration(
            stage_params_list=stages, scheme_type="NPT"
        )
    """

    name: str
    ensemble: str
    time_ns: float
    timestep: float
    temperature: float
    constraints: Dict[str, float] = field(default_factory=dict)
    minimize_steps: int = 0
    dcd_freq: int = 5000
    # NAMD-specific optional fields
    steps: Optional[int] = None
    pressure: Optional[float] = None
    surface_tension: Optional[float] = None

    def replace(self, **kwargs) -> "EquilibrationStage":
        """Return a copy of this stage with the given fields overridden.

        Example::

            fast_prod = stages[-1].replace(time_ns=10.0, temperature=303.15)
        """
        return _dc_replace(self, **kwargs)

    def to_dict(self) -> Dict[str, Any]:
        """Serialise to a plain dict compatible with the setup methods."""
        d: Dict[str, Any] = {
            "name": self.name,
            "ensemble": self.ensemble,
            "time_ns": self.time_ns,
            "timestep": self.timestep,
            "temperature": self.temperature,
            "constraints": dict(self.constraints),
            "minimize_steps": self.minimize_steps,
            "dcd_freq": self.dcd_freq,
        }
        if self.steps is not None:
            d["steps"] = self.steps
        if self.pressure is not None:
            d["pressure"] = self.pressure
        if self.surface_tension is not None:
            d["surface_tension"] = self.surface_tension
        return d


class EquilibrationProtocol:
    """Base class for equilibration protocols."""

    def __init__(self, name: str, description: str = ""):
        self.name = name
        self.description = description
        self.stages = []

    def add_stage(self, stage: Dict[str, Any]):
        """Add an equilibration stage to the protocol."""
        self.stages.append(stage)

    def to_dict(self) -> Dict[str, Any]:
        """Convert protocol to dictionary."""
        return {
            "name": self.name,
            "description": self.description,
            "stages": self.stages,
        }

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EquilibrationProtocol":
        """Create protocol from dictionary."""
        protocol = cls(data["name"], data.get("description", ""))
        protocol.stages = data.get("stages", [])
        return protocol


class NAMDEquilibrationManager:
    """Manager for NAMD equilibration simulations."""

    def __init__(self, working_dir: Path, namd_executable: str = "namd3"):
        self.working_dir = Path(working_dir)
        self.namd_executable = namd_executable
        self.logger = get_logger(self.__class__.__name__)

        # Path to NAMD templates (homogeneous layout with other engines)
        self.namd_templates_dir = (
            Path(__file__).parent.parent.parent / "equilibration" / "namd"
        )

    def find_system_files(self) -> Optional[Dict[str, str]]:
        """
        Automatically detect system files in working directory.

        Looks for standard AMBER system files and bilayer PDB with CRYST1 record.

        Returns:
            Dictionary with system file paths, or None if required files not found:
            {
                'prmtop': Path to system.prmtop,
                'inpcrd': Path to system.inpcrd,
                'pdb': Path to system.pdb,
                'bilayer_pdb': Path to bilayer PDB with CRYST1
            }

        Example:
            >>> manager = NAMDEquilibrationManager(Path("work_dir"))
            >>> system_files = manager.find_system_files()
            >>> if system_files:
            ...     result = manager.setup_namd_equilibration(
            ...         system_files=system_files,
            ...         stage_params_list=stages
            ...     )
        """
        system_files = {}

        # Find AMBER topology file
        prmtop_files = list(self.working_dir.glob("*.prmtop"))
        if not prmtop_files:
            self.logger.error("No .prmtop file found in working directory")
            return None
        system_files["prmtop"] = str(prmtop_files[0])
        self.logger.info(f"Found topology: {prmtop_files[0].name}")

        # Find AMBER coordinate file
        inpcrd_files = list(self.working_dir.glob("*.inpcrd"))
        if not inpcrd_files:
            # Try alternative extensions
            inpcrd_files = list(self.working_dir.glob("*.crd"))
            if not inpcrd_files:
                inpcrd_files = list(self.working_dir.glob("*.rst"))

        if not inpcrd_files:
            self.logger.error("No .inpcrd/.crd/.rst file found in working directory")
            return None
        system_files["inpcrd"] = str(inpcrd_files[0])
        self.logger.info(f"Found coordinates: {inpcrd_files[0].name}")

        # Find system PDB file
        system_pdb = self.working_dir / "system.pdb"
        if not system_pdb.exists():
            # Try to find any PDB that's not a bilayer file
            pdb_files = [
                f
                for f in self.working_dir.glob("*.pdb")
                if "bilayer" not in f.name.lower()
            ]
            if pdb_files:
                system_pdb = pdb_files[0]
            else:
                self.logger.error("No system.pdb file found in working directory")
                return None
        system_files["pdb"] = str(system_pdb)
        self.logger.info(f"Found system PDB: {system_pdb.name}")

        # Find bilayer PDB with CRYST1 record
        bilayer_pdb = self._find_bilayer_pdb_with_cryst1()
        if not bilayer_pdb:
            self.logger.error(
                "No bilayer PDB with CRYST1 record found in working directory"
            )
            self.logger.error(
                "Required: bilayer*_lipid.pdb file from packmol-memgen --parametrize"
            )
            return None
        system_files["bilayer_pdb"] = str(bilayer_pdb)
        self.logger.info(f"Found bilayer PDB with CRYST1: {bilayer_pdb.name}")

        return system_files

    def _get_config_name(self, stage_name: str, stage_index: int) -> str:
        """
        Convert GUI display names to valid config file names.

        Maps names like 'Equilibration 1', 'Equilibration 2', etc. to 'step1', 'step2', etc.
        For custom names, uses stage_index to generate sequential step names.

        Args:
            stage_name: Display name from GUI (e.g., "Equilibration 1" or "Initial Equilibration")
            stage_index: Zero-based index of the stage in the protocol

        Returns:
            Valid config file name (e.g., "step1", "step2", etc.)
        """
        # Handle the standard naming convention with spaces "Equilibration N"
        if stage_name.startswith("Equilibration "):
            try:
                stage_num = stage_name.split()[1]
                return f"step{stage_num}"
            except (IndexError, ValueError):
                pass

        # Handle Production stage
        if stage_name == "Production":
            return "step7_production"

        # Handle legacy names (in case they exist) - convert to new convention
        legacy_mapping = {
            "equilibration_1": "step1",
            "equilibration_2": "step2",
            "equilibration_3": "step3",
            "equilibration_4": "step4",
            "equilibration_5": "step5",
            "equilibration_6": "step6",
            "eq1": "step1",  # Convert old naming
            "eq2": "step2",
            "eq3": "step3",
            "eq4": "step4",
            "eq5": "step5",
            "eq6": "step6",
            "production": "step7_production",
        }

        # Check legacy mapping first
        if stage_name.lower() in legacy_mapping:
            return legacy_mapping[stage_name.lower()]

        # For custom stage names, use stage_index to generate sequential step names
        # This ensures proper ordering: step1, step2, step3, etc.
        return f"step{stage_index + 1}"

    def _read_box_dimensions(self, pdb_file: Path) -> Tuple[float, float, float]:
        """
        Read box dimensions from PDB file.

        First tries to read CRYST1 record, then estimates from coordinates.
        Example CRYST1 line: CRYST1   70.335   70.833   85.067  90.00  90.00  90.00 P 1           1

        Args:
            pdb_file: Path to PDB file

        Returns:
            Tuple of (a, b, c) box dimensions in Angstroms
        """
        try:
            with open(pdb_file, "r") as f:
                for line in f:
                    if line.startswith("CRYST1"):
                        # CRYST1 record contains unit cell parameters
                        # Format: CRYST1    a       b       c     alpha  beta  gamma sgroup    z
                        #         CRYST1 70.335  70.833  85.067  90.00  90.00  90.00 P 1        1
                        a = float(line[6:15].strip())
                        b = float(line[15:24].strip())
                        c = float(line[24:33].strip())
                        self.logger.info(
                            f"Read box dimensions from CRYST1 record: {a:.2f} x {b:.2f} x {c:.2f} Å"
                        )
                        return a, b, c
        except (FileNotFoundError, ValueError, IndexError) as e:
            self.logger.warning(f"Could not read CRYST1 from {pdb_file}: {e}")

        # If no CRYST1 record found, estimate from coordinates
        try:
            self.logger.info(
                f"No CRYST1 record found in {pdb_file.name}, estimating box dimensions from coordinates"
            )
            return self._estimate_box_from_coordinates(pdb_file)
        except Exception as e:
            self.logger.warning(
                f"Could not estimate box dimensions from {pdb_file}: {e}"
            )

        # Return default dimensions for membrane systems (approximate)
        self.logger.info("Using default box dimensions (100x100x100 Å)")
        return 100.0, 100.0, 100.0

    def _estimate_box_from_coordinates(
        self, pdb_file: Path
    ) -> Tuple[float, float, float]:
        """
        Estimate box dimensions from coordinate extremes with padding.

        Args:
            pdb_file: Path to PDB file

        Returns:
            Tuple of (a, b, c) box dimensions in Angstroms
        """
        x_coords, y_coords, z_coords = [], [], []

        with open(pdb_file, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        x_coords.append(x)
                        y_coords.append(y)
                        z_coords.append(z)
                    except (ValueError, IndexError):
                        continue

        if not x_coords:
            raise ValueError("No valid coordinates found in PDB file")

        # Calculate ranges and add padding (10 Å on each side)
        padding = 10.0
        x_range = max(x_coords) - min(x_coords) + 2 * padding
        y_range = max(y_coords) - min(y_coords) + 2 * padding
        z_range = max(z_coords) - min(z_coords) + 2 * padding

        self.logger.info(
            f"Estimated box dimensions: {x_range:.2f} x {y_range:.2f} x {z_range:.2f} Å"
        )
        return x_range, y_range, z_range

    def _calculate_pme_grid_size(
        self, box_dimension: float, cutoff: float = 9.0
    ) -> int:
        """
        Calculate optimal PME grid size for a given box dimension.

        The PME grid size should be:
        1. At least 2x the cutoff distance in grid points
        2. Efficiently factorable (preferably powers of 2, 3, 5)
        3. About 1 Å per grid point for good accuracy

        Args:
            box_dimension: Box dimension in Angstroms
            cutoff: Electrostatic cutoff in Angstroms

        Returns:
            Optimal PME grid size
        """
        # Rule of thumb: ~1 Å per grid point, but at least 2x cutoff
        min_grid_size = max(int(box_dimension), int(2 * cutoff))

        # Find the next efficiently factorable number
        # NAMD/FFTW work best with numbers that factor into small primes (2, 3, 5)
        grid_size = self._find_efficient_grid_size(min_grid_size)

        self.logger.debug(
            f"Box dimension: {box_dimension:.2f} Å, "
            f"Min grid size: {min_grid_size}, "
            f"Optimal grid size: {grid_size}"
        )

        return grid_size

    def _find_efficient_grid_size(self, min_size: int) -> int:
        """
        Find the smallest efficiently factorable number >= min_size.

        Efficient numbers for FFT are those that factor into small primes (2, 3, 5).

        Args:
            min_size: Minimum required grid size

        Returns:
            Efficient grid size >= min_size
        """
        if min_size <= 1:
            return 1

        # Generate efficient numbers by multiplying powers of 2, 3, 5
        efficient_numbers = []

        # Generate numbers up to a reasonable limit (min_size * 2)
        limit = min_size * 2

        # Powers of 2
        power_2 = 1
        while power_2 <= limit:
            # Powers of 3
            power_3 = power_2
            while power_3 <= limit:
                # Powers of 5
                power_5 = power_3
                while power_5 <= limit:
                    efficient_numbers.append(power_5)
                    power_5 *= 5
                power_3 *= 3
            power_2 *= 2

        # Sort and find the first number >= min_size
        efficient_numbers.sort()

        for num in efficient_numbers:
            if num >= min_size:
                return num

        # Fallback: if no efficient number found, use next power of 2
        power_of_2 = 1
        while power_of_2 < min_size:
            power_of_2 *= 2

        return power_of_2

    def _find_bilayer_pdb_with_cryst1(self) -> Optional[Path]:
        """
        Find bilayer PDB file that contains CRYST1 record for box dimensions.

        Prioritizes bilayer*_lipid.pdb files generated by packmol-memgen --parametrize,
        which contain essential CRYST1 box information for MD simulations.

        Returns:
            Path to bilayer PDB with CRYST1 record, or None if not found
        """
        # First priority: Look for final prepared files with pattern bilayer*_lipid.pdb
        final_pattern = "bilayer*_lipid.pdb"
        final_files = list(self.working_dir.glob(final_pattern))

        for pdb_file in final_files:
            if self._is_final_prepared_pdb(pdb_file):
                self.logger.info(
                    f"Found final prepared bilayer PDB with CRYST1: {pdb_file}"
                )
                return pdb_file

        # Second priority: Look for other bilayer files with CRYST1
        other_patterns = ["bilayer_*.pdb", "*_bilayer.pdb", "*membrane*.pdb"]

        for pattern in other_patterns:
            bilayer_files = list(self.working_dir.glob(pattern))
            for pdb_file in bilayer_files:
                # Skip lipid-only files unless it's the final prepared pattern
                if "lipid" in pdb_file.name.lower() and not pdb_file.name.endswith(
                    "_lipid.pdb"
                ):
                    continue

                # Check if this file has CRYST1 record and is not intermediate
                if self._is_final_prepared_pdb(pdb_file):
                    self.logger.info(f"Found bilayer PDB with CRYST1: {pdb_file}")
                    return pdb_file

        return None

    def _find_original_bilayer_pdb(self) -> Optional[Path]:
        """
        Find the final prepared bilayer PDB file (from packmol-memgen).

        The correct file should:
        1. Have pattern bilayer*_lipid.pdb (final prepared file)
        2. Contain CRYST1 header (properly formatted)
        3. NOT start with "REMARK   Packmol generated pdb file" (intermediate file)

        Returns:
            Path to the final prepared bilayer PDB file, or None if not found
        """
        # First priority: Look for final prepared files with pattern bilayer*_lipid.pdb
        final_pattern = "bilayer*_lipid.pdb"
        final_files = list(self.working_dir.glob(final_pattern))

        for pdb_file in final_files:
            # Verify this is the correct final file (has CRYST1, not intermediate)
            if self._is_final_prepared_pdb(pdb_file):
                self.logger.info(f"Found final prepared bilayer PDB: {pdb_file}")
                return pdb_file

        # Second priority: Look for other bilayer files with CRYST1 header
        other_patterns = ["bilayer_*.pdb", "*_bilayer.pdb"]

        for pattern in other_patterns:
            bilayer_files = list(self.working_dir.glob(pattern))
            for pdb_file in bilayer_files:
                # Skip if it's an intermediate file or doesn't have CRYST1
                if not self._is_final_prepared_pdb(pdb_file):
                    continue
                self.logger.info(f"Found bilayer PDB with CRYST1: {pdb_file}")
                return pdb_file

        # Fallback: Any bilayer file (warn user about potential issues)
        for pattern in ["bilayer_*.pdb", "*_bilayer.pdb"]:
            bilayer_files = list(self.working_dir.glob(pattern))
            if bilayer_files:
                self.logger.warning(
                    f"Using fallback bilayer PDB (may not have CRYST1): {bilayer_files[0]}"
                )
                return bilayer_files[0]

        return None

    def _is_final_prepared_pdb(self, pdb_file: Path) -> bool:
        """
        Check if a PDB file is the final prepared file (not an intermediate).

        The final prepared file should:
        1. Have CRYST1 header line
        2. NOT start with intermediate file markers like "REMARK   Packmol generated pdb file"

        Args:
            pdb_file: Path to PDB file to check

        Returns:
            True if this is the final prepared file, False otherwise
        """
        try:
            with open(pdb_file, "r") as f:
                lines = f.readlines()

            if not lines:
                return False

            has_cryst1 = False
            is_intermediate = False

            # Check first few lines for indicators
            for i, line in enumerate(lines[:10]):  # Check first 10 lines
                line = line.strip()

                # Check for CRYST1 header (good indicator of final file)
                if line.startswith("CRYST1"):
                    has_cryst1 = True

                # Check for intermediate file markers (bad indicators)
                if (
                    "Packmol generated pdb file" in line
                    and "Packmol Memgen estimated parameters" in line
                ):
                    is_intermediate = True
                    break

                if "charmmlipid2amber.py transformed file" in line:
                    is_intermediate = True
                    break

            # Final file should have CRYST1 and NOT be intermediate
            result = has_cryst1 and not is_intermediate

            if result:
                self.logger.debug(
                    f"✅ {pdb_file.name} is final prepared file (has CRYST1, not intermediate)"
                )
            else:
                if is_intermediate:
                    self.logger.debug(
                        f"❌ {pdb_file.name} is intermediate file (packmol-memgen generated)"
                    )
                elif not has_cryst1:
                    self.logger.debug(f"❌ {pdb_file.name} missing CRYST1 header")

            return result

        except Exception as e:
            self.logger.warning(f"Could not check PDB file {pdb_file}: {e}")
            return False

    def _read_amber_box_dimensions(
        self, coord_file: Path
    ) -> Tuple[float, float, float]:
        """
        Read box dimensions from AMBER coordinate (.inpcrd, .crd or .rst) file.

        Args:
            coord_file: Path to AMBER coordinate file

        Returns:
            Tuple of (a, b, c) box dimensions in Angstroms
        """
        try:
            with open(coord_file, "r") as f:
                lines = f.readlines()

            # Box information is in the last line of AMBER coordinate files
            if len(lines) >= 2:
                last_line = lines[-1].strip()
                box_values = last_line.split()

                if len(box_values) >= 3:
                    a = float(box_values[0])
                    b = float(box_values[1])
                    c = float(box_values[2])
                    self.logger.info(
                        f"Read box dimensions from AMBER file: {a:.2f} x {b:.2f} x {c:.2f} Å"
                    )
                    return a, b, c

        except (FileNotFoundError, ValueError, IndexError) as e:
            self.logger.warning(
                f"Could not read box dimensions from AMBER file {coord_file}: {e}"
            )

        # Fallback to default
        self.logger.info(
            "Using default box dimensions for AMBER system (100x100x100 Å)"
        )
        return 100.0, 100.0, 100.0

    def generate_config_file(
        self,
        stage_name: str,
        stage_params: Dict[str, Any],
        stage_index: int,
        system_files: Dict[str, str],
        previous_stage_name: Optional[str] = None,
    ) -> str:
        """
        Generate NAMD configuration file for a specific equilibration stage using AMBER force field.

        Args:
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths (should include 'prmtop' and 'inpcrd')
            previous_stage_name: Name of the previous stage for restart files (optional)

        Returns:
            Content of the NAMD configuration file
        """

        config_lines = []

        # Header
        config_lines.extend(
            [
                "#############################################################",
                f"## NAMD Configuration File for {stage_params.get('name', stage_name)}",
                f"## Generated by Gatewizard",
                f"## Stage {stage_index + 1}: {stage_params.get('description', '')}",
                f"## Force Field: AMBER",
                "#############################################################",
                "",
            ]
        )

        # Input files - AMBER format (now using local copies)
        config_lines.extend(
            [
                "# Input files - AMBER format",
                "amber              on",
                f"parmfile           system.prmtop",
                f"ambercoor          system.inpcrd",
            ]
        )

        # Restart files for stages after the first
        if stage_index > 0:
            # Use passed previous stage name or build from previous stage index
            if previous_stage_name:
                prev_stage = self._generate_output_name(
                    previous_stage_name, stage_index - 1
                )
            else:
                # Fallback: use the input name generation method
                prev_stage = self._generate_input_name(stage_index, previous_stage_name)

            config_lines.extend(
                [
                    f"bincoordinates     {prev_stage}.coor",
                    f"binvelocities      {prev_stage}.vel",
                    f"extendedsystem     {prev_stage}.xsc",
                ]
            )

        config_lines.append("")

        # AMBER Force field settings
        config_lines.extend(
            [
                "# AMBER Force field settings",
                "exclude            scaled1-4",
                "oneFourScaling         0.833333333",  # = 1/1.2 (SCEE=1.2 in AMBER)
                "scnb               2.0",  # VDW 1-4 scaling factor
                "readexclusions     yes",  # Read exclusions from PARM file
                "switching          off",  # Turn off switching (AMBER default)
                "LJcorrection       on",  # Apply analytical tail correction
                "zeromomentum       on",  # Remove COM drift (netfrc=1 in AMBER)
                "",
            ]
        )

        # Output files
        output_name = self._generate_output_name(stage_name, stage_index)
        dcd_freq = stage_params.get("dcd_freq", 5000)

        config_lines.extend(
            [
                "# Output files",
                f"outputName         {output_name}",
                f"dcdfile            {output_name}.dcd",
                f"dcdfreq            {dcd_freq}",
                f"outputEnergies     {dcd_freq}",
                f"outputPressure     {dcd_freq}",
                f"outputTiming       {dcd_freq}",
                f"xstFreq            {dcd_freq}",
                "",
            ]
        )

        # Basic simulation parameters
        timestep = stage_params.get(
            "timestep", 2.0
        )  # Default 2 fs (NAMD uses femtoseconds)
        steps = stage_params.get("steps", 125000)

        config_lines.extend(
            [
                "# Simulation parameters",
                f"timestep           {timestep}",
                f"nonbondedFreq      1",
                f"fullElectFrequency 1",  # AMBER default
                f"stepspercycle      10",
                "# Note: numsteps parameter is intentionally omitted",
                "# Steps are controlled by the 'run' command at the end",
                "",
                "# AMBER-compatible force field settings",
                "rigidBonds         all",  # SHAKE all bonds (ntc=2, ntf=2)
                "rigidIterations    100",
                f"cutoff             {stage_params.get('cutoff', 9.0)}",  # Default AMBER cutoff
                f"pairlistdist       {stage_params.get('cutoff', 9.0) + 2.0}",  # cutoff + 2.0
            ]
        )
        from gatewizard.tools.namd_water import namd_water_model_config_block

        wm = getattr(self, "water_model", "tip3p")
        for line in namd_water_model_config_block(wm).splitlines():
            stripped = line.strip()
            if stripped and not stripped.startswith("#"):
                config_lines.append(stripped)
        if wm == "tip3p" or not any(
            l.startswith("waterModel") for l in config_lines[-5:]
        ):
            config_lines.extend(
                [
                    "rigidTolerance     1.0e-8",
                    "useSettle          on",
                    "watermodel         tip3",
                ]
            )
        config_lines.append("")

        # Temperature control - Langevin thermostat (corresponds to ntt=3 in AMBER)
        temperature = stage_params.get("temperature", 310.15)
        ensemble = stage_params.get("ensemble", "NPT")

        if ensemble in ["NVT", "NPT", "NPAT", "NPgT"]:
            damping = stage_params.get("langevin_damping", 5.0)  # gamma_ln in AMBER

            # For first stage, set initial temperature
            # For subsequent stages, only set thermostat parameters (velocities come from restart)
            if stage_index == 0:
                config_lines.extend(
                    [
                        "# Temperature control - Langevin thermostat (first stage)",
                        f"temperature        {temperature}",  # tempi in AMBER
                        "langevin           on",  # ntt=3 in AMBER
                        f"langevinTemp       {temperature}",  # temp0 in AMBER
                        f"langevinDamping    {damping}",  # gamma_ln in AMBER (ps^-1)
                        "langevinHydrogen   off",  # AMBER default
                        "",
                    ]
                )
            else:
                config_lines.extend(
                    [
                        "# Temperature control - Langevin thermostat (restart stage)",
                        "langevin           on",  # ntt=3 in AMBER
                        f"langevinTemp       {temperature}",  # temp0 in AMBER
                        f"langevinDamping    {damping}",  # gamma_ln in AMBER (ps^-1)
                        "langevinHydrogen   off",  # AMBER default
                        "",
                    ]
                )

        # Pressure control - Berendsen barostat (corresponds to ntp=1 in AMBER)
        if ensemble in ["NPT", "NPAT", "NPgT"]:
            pressure = stage_params.get("pressure", 1.01325)  # AMBER default pressure
            surface_tension = stage_params.get("surface_tension", 0.0)

            # Use Berendsen pressure control for AMBER compatibility
            config_lines.extend(
                [
                    "# Pressure control - Berendsen barostat",
                    "BerendsenPressure     on",  # ntp=1 in AMBER
                    f"BerendsenPressureTarget {pressure}",  # pres0 in AMBER
                    "BerendsenPressureCompressibility  4.57e-5",  # compressibility in AMBER (1/bar)
                    "BerendsenPressureRelaxationTime 100.0",  # taup in AMBER (fs in NAMD, ps in AMBER)
                    "useGroupPressure      yes",  # needed for rigidBonds
                ]
            )

            # Configure pressure scaling based on ensemble
            if ensemble == "NPAT":
                # Semi-isotropic pressure control for membrane systems
                config_lines.extend(
                    [
                        "useFlexibleCell       yes",  # allow cell shape changes
                        "useConstantRatio      yes",  # keep XY ratio constant
                    ]
                )

                # Add surface tension if specified
                if surface_tension > 0.0:
                    config_lines.extend(
                        [
                            f"# Surface tension control for NPAT ensemble",
                            f"surfaceTensionTarget  {surface_tension}",  # dyn/cm
                        ]
                    )

            elif ensemble == "NPgT":
                # Constant surface tension ensemble
                config_lines.extend(
                    [
                        "useFlexibleCell       yes",  # allow cell shape changes
                        "useConstantRatio      yes",  # keep XY ratio constant
                    ]
                )

                # Surface tension is required for NPgT
                if surface_tension == 0.0:
                    surface_tension = 0.0  # Default value in dyn/cm
                    self.logger.warning(
                        f"NPgT ensemble requires surface tension. Using default: {surface_tension} dyn/cm"
                    )

                config_lines.extend(
                    [
                        f"# Surface tension control for NPgT ensemble",
                        f"surfaceTensionTarget  {surface_tension}",  # dyn/cm
                    ]
                )

            else:  # NPT
                # Isotropic pressure control
                config_lines.extend(
                    [
                        "useFlexibleCell       no",  # isotropic scaling
                        "useConstantArea       no",
                    ]
                )

            config_lines.append("")

        # Restraints/Constraints (if needed)
        constraints = stage_params.get("constraints", {})
        has_restraints = any(float(v) > 0 for v in constraints.values())

        if has_restraints:
            # Calculate constraint scaling based on the maximum restraint force
            max_force = max(float(v) for v in constraints.values() if float(v) > 0)
            constraint_scaling = min(
                10.0, max(1.0, max_force)
            )  # Scale between 1.0 and 10.0

            # Use stage-specific restraints file if available, otherwise general file
            config_name = self._get_config_name(stage_name, stage_index)
            if config_name == "step7_production":
                stage_restraints_file = f"restraints/{config_name}_restraints.pdb"
            else:
                stage_restraints_file = (
                    f"restraints/{config_name}_equilibration_restraints.pdb"
                )
            general_restraints_file = "restraints.pdb"

            config_lines.extend(
                [
                    "# Harmonic restraints",
                    "constraints        on",
                    "consexp            2",  # harmonic restraints
                    f"consref            {stage_restraints_file}",  # reference coordinates
                    f"conskfile          {stage_restraints_file}",  # force constant file
                    "conskcol           B",  # use B-factor column
                    f"constraintScaling  {constraint_scaling}",  # overall scaling factor
                    f"# Note: If {stage_restraints_file} not found, use {general_restraints_file}",
                    "",
                ]
            )

        # Non-bonded interactions - AMBER settings are already included above
        # The cutoff, pairlistdist, switching, exclude, etc. are set in the simulation parameters

        # Get box dimensions for PME grid size calculation
        # PRIORITY: Use bilayer_pdb from system_files for CRYST1 record (most accurate)
        bilayer_pdb_path = self.working_dir / system_files.get("bilayer_pdb", "")
        inpcrd_file = self.working_dir / system_files.get("inpcrd", "system.inpcrd")

        # Try to read from bilayer PDB with CRYST1 first (highest priority)
        if bilayer_pdb_path.exists():
            box_a, box_b, box_c = self._read_box_dimensions(bilayer_pdb_path)
            self.logger.info(
                f"Using box dimensions from bilayer PDB for cell basis: {bilayer_pdb_path.name}"
            )
        elif inpcrd_file.exists():
            box_a, box_b, box_c = self._read_amber_box_dimensions(inpcrd_file)
            self.logger.info(
                f"Using box dimensions from AMBER inpcrd file: {inpcrd_file.name}"
            )
        else:
            self.logger.warning(
                "No coordinate file found, using default box dimensions"
            )
            box_a, box_b, box_c = 100.0, 100.0, 100.0

        # PME electrostatics - AMBER compatible settings
        config_lines.extend(
            [
                "# PME electrostatics - AMBER compatible",
                "PME                yes",
                "PMETolerance       1.0e-6",  # dsum_tol in AMBER
                "PMEInterpOrder     4",  # order=4 in AMBER (cubic spline)
                "PMEGridSpacing     1.0",  # Let NAMD automatically calculate grid sizes
                "",
            ]
        )

        # Periodic boundary conditions
        # Box dimensions already calculated above for PME grid sizing
        config_lines.extend(
            [
                "# Periodic boundary conditions",
                f"# Box dimensions: {box_a:.2f} x {box_b:.2f} x {box_c:.2f} Å",
                f"cellBasisVector1   {box_a:.6f}   0.000000   0.000000",
                f"cellBasisVector2   0.000000   {box_b:.6f}   0.000000",
                f"cellBasisVector3   0.000000   0.000000   {box_c:.6f}",
                "cellOrigin         0.0   0.0   0.0",
                "",
            ]
        )

        # Wrap output
        config_lines.extend(
            ["# Wrap output", "wrapWater          on", "wrapAll            on", ""]
        )

        # Minimization specific settings
        integrator = stage_params.get("integrator", "")
        if "minimization" in stage_name.lower() or integrator == "conjugate_gradient":
            minimize_steps = min(5000, steps)
            config_lines.extend(
                ["# Energy minimization", f"minimize           {minimize_steps}", ""]
            )

        # Run command
        if "minimization" not in stage_name.lower():
            config_lines.extend(["# Run the simulation", f"run               {steps}"])

        return "\n".join(config_lines)

    # ------------------------------------------------------------------ #
    #  MDAnalysis-based restraint selection helpers                        #
    # ------------------------------------------------------------------ #

    # Default MDAnalysis selection strings for each restraint category
    DEFAULT_SELECTIONS = {
        "protein_backbone": "protein and backbone",
        "protein_sidechain": "protein and not backbone",
        "lipid_head": (
            "(resname POPC POPE POPS DPPC DMPC DOPC DSPC PC PE PS PA PG PI SM "
            "OL LA MY ST AR OLE PAL STE LIN CHOL CHL CHOLEST PALM OLEO STEROL) "
            "and (name P O11 O12 O13 O14 O21 O22 O31 O32 O33 O34 "
            "O1P O2P O3P O4P OP1 OP2 OP3 OP4 "
            "N C11 C12 C13 C14 N31 C32 C33 C34 C35 C1 C2 C3 "
            "HN1 HN2 HN3 HO2 HO3 HS)"
        ),
        "lipid_tail": (
            "(resname POPC POPE POPS DPPC DMPC DOPC DSPC PC PE PS PA PG PI SM "
            "OL LA MY ST AR OLE PAL STE LIN CHOL CHL CHOLEST PALM OLEO STEROL) "
            "and not (name P O11 O12 O13 O14 O21 O22 O31 O32 O33 O34 "
            "O1P O2P O3P O4P OP1 OP2 OP3 OP4 "
            "N C11 C12 C13 C14 N31 C32 C33 C34 C35 C1 C2 C3 "
            "HN1 HN2 HN3 HO2 HO3 HS)"
        ),
        "water": "resname TIP3 HOH WAT SOL TIP4 SPC T3P T4P",
        "ions": (
            "resname NA CL K CA MG ZN FE CU SOD CLA POT CAL MAG ZIN IRN COP "
            "Na+ Cl- K+ Ca2+ Mg2+ Zn2+ Fe2+ Fe3+ Cu2+ "
            "NA+ CL- LIT RUB CES BAR"
        ),
        "other": (
            "not (protein or "
            "(resname POPC POPE POPS DPPC DMPC DOPC DSPC PC PE PS PA PG PI SM "
            "OL LA MY ST AR OLE PAL STE LIN CHOL CHL CHOLEST PALM OLEO STEROL) or "
            "(resname TIP3 HOH WAT SOL TIP4 SPC T3P T4P) or "
            "(resname NA CL K CA MG ZN FE CU SOD CLA POT CAL MAG ZIN IRN COP "
            "Na+ Cl- K+ Ca2+ Mg2+ Zn2+ Fe2+ Fe3+ Cu2+ "
            "NA+ CL- LIT RUB CES BAR))"
        ),
    }

    @staticmethod
    def count_selection_atoms(pdb_path: str, selection: str) -> int:
        """
        Count atoms matching an MDAnalysis selection expression.

        Args:
            pdb_path: Path to a PDB file.
            selection: MDAnalysis selection string.

        Returns:
            Number of atoms matching the selection (0 if selection is invalid
            or MDAnalysis is not available).
        """
        try:
            import MDAnalysis as mda
            import warnings

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(str(pdb_path))
                ag = u.select_atoms(selection)
                return len(ag)
        except Exception as exc:
            logger.debug(f"count_selection_atoms failed for '{selection}': {exc}")
            return 0

    @staticmethod
    def get_default_selections(pdb_path: str) -> Dict[str, str]:
        """
        Build the default selection dict, auto-detecting extra ligands / residues.

        The seven standard categories (protein_backbone, protein_sidechain,
        lipid_head, lipid_tail, water, ions, other) are always present.
        Any residue that falls into the *other* category is additionally
        split into individual ``ligand_<RESNAME>`` entries so users can
        assign per-ligand restraint forces.

        Args:
            pdb_path: Path to a PDB file.

        Returns:
            ``{category_name: mda_selection_string, ...}`` including any
            auto-detected ligand entries.
        """
        sels = dict(NAMDEquilibrationManager.DEFAULT_SELECTIONS)

        try:
            import MDAnalysis as mda
            import warnings

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(str(pdb_path))

                # Identify residues in "other" (not protein/lipid/water/ion)
                other_ag = u.select_atoms(sels["other"])
                if len(other_ag) > 0:
                    unique_resnames = sorted(set(other_ag.resnames))
                    for resname in unique_resnames:
                        key = f"ligand_{resname}"
                        sels[key] = f"resname {resname}"
        except Exception as exc:
            logger.debug(f"get_default_selections ligand detection failed: {exc}")

        return sels

    @staticmethod
    def count_all_selections(
        pdb_path: str, selections: Optional[Dict[str, str]] = None
    ) -> Dict[str, int]:
        """
        Count atoms for every selection in the dict.

        Args:
            pdb_path: Path to a PDB file.
            selections: ``{name: mda_selection_string}``.  If *None* the
                default selections (with auto-detected ligands) are used.

        Returns:
            ``{name: atom_count, ...}``
        """
        if selections is None:
            selections = NAMDEquilibrationManager.get_default_selections(pdb_path)

        counts: Dict[str, int] = {}
        try:
            import MDAnalysis as mda
            import warnings

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(str(pdb_path))
                for name, sel_str in selections.items():
                    try:
                        counts[name] = len(u.select_atoms(sel_str))
                    except Exception:
                        counts[name] = 0
        except Exception as exc:
            logger.debug(f"count_all_selections failed: {exc}")
            for name in selections:
                counts[name] = 0
        return counts

    def generate_restraints_file_mda(
        self,
        system_pdb: Path,
        selections_with_forces: Dict[str, Tuple[str, float]],
        output_file: Path,
        stage_name: str = "",
    ) -> None:
        """
        Generate a restraints PDB using MDAnalysis selections.

        Each entry in *selections_with_forces* maps a category name to a
        ``(mda_selection_string, force)`` tuple.  For every ATOM/HETATM line
        the **first matching** selection determines the B-factor written.
        Atoms that match no selection get B-factor 0.0.

        Args:
            system_pdb: Path to the system PDB file.
            selections_with_forces: ``{name: (selection_string, force), ...}``
            output_file: Destination path for the restraints PDB.
            stage_name: Label used in log messages.
        """
        import MDAnalysis as mda
        import warnings

        if not system_pdb.exists():
            self.logger.error(f"System PDB not found: {system_pdb}")
            return

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(str(system_pdb))

        # Pre-compute index sets for each selection
        ordered_selections: List[Tuple[str, set, float]] = []
        stats: Dict[str, int] = {}
        for name, (sel_str, force) in selections_with_forces.items():
            try:
                indices = set(u.select_atoms(sel_str).indices)
            except Exception as exc:
                self.logger.warning(f"Selection '{name}' failed: {exc}")
                indices = set()
            ordered_selections.append((name, indices, force))
            stats[name] = 0

        # Build per-atom force array (first match wins)
        n_atoms = len(u.atoms)
        forces = [0.0] * n_atoms
        for atom_idx in range(n_atoms):
            for name, idx_set, force in ordered_selections:
                if atom_idx in idx_set:
                    forces[atom_idx] = force
                    stats[name] += 1
                    break

        # Read original PDB and replace B-factors
        with open(system_pdb, "r") as fh:
            lines = fh.readlines()

        out_lines: List[str] = []
        atom_serial = 0
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                new_line = line[:60] + f"{forces[atom_serial]:6.2f}" + line[66:]
                out_lines.append(new_line)
                atom_serial += 1
            else:
                out_lines.append(line)

        with open(output_file, "w") as fh:
            fh.writelines(out_lines)

        self.logger.info(f"Generated MDAnalysis restraints: {output_file}")
        self.logger.info(f"Stage: {stage_name} | Total atoms: {n_atoms}")
        for name, count in stats.items():
            if count > 0:
                force_val = selections_with_forces[name][1]
                self.logger.info(
                    f"  {name}: {count} atoms, force = {force_val} kcal/mol/Å²"
                )

    def generate_restraints_file(
        self,
        system_pdb: Path,
        constraints: Dict[str, float],
        output_file: Path,
        stage_name: str = "",
        selections: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        Generate restraints PDB file with B-factors for constraint forces.

        Uses the final system.pdb file for generating restraints to ensure
        consistency with the parametrized system used in simulations.

        When *selections* is provided, each key in *constraints* is resolved
        to the corresponding MDAnalysis selection string and
        ``generate_restraints_file_mda`` is used.  Otherwise the legacy
        residue-name heuristic is applied.

        Args:
            system_pdb: Path to the system.pdb file
            constraints: Dictionary of constraint types and forces
            output_file: Path for output restraints file
            stage_name: Name of the equilibration stage (for documentation)
            selections: Optional ``{constraint_name: mda_selection_string}``
                mapping.  If provided, MDAnalysis is used instead of the
                built-in classification heuristic.
        """

        # If MDAnalysis selections are provided, delegate to the MDA method
        if selections:
            sel_with_forces: Dict[str, Tuple[str, float]] = {}
            for name, force in constraints.items():
                sel_str = selections.get(name)
                if sel_str:
                    sel_with_forces[name] = (sel_str, force)
            if sel_with_forces:
                self.generate_restraints_file_mda(
                    system_pdb, sel_with_forces, output_file, stage_name
                )
                return

        # Use only the final system.pdb file for restraints
        if not system_pdb.exists():
            self.logger.error(f"System PDB file not found: {system_pdb}")
            return

        self.logger.info(f"Using system.pdb for restraints: {system_pdb.name}")

        # Read system PDB file
        with open(system_pdb, "r") as f:
            lines = f.readlines()

        # Process each line and assign restraint forces
        restraint_lines = []
        atom_count = 0
        restraint_stats = {
            "protein_backbone": 0,
            "protein_sidechain": 0,
            "lipid_head": 0,
            "lipid_tail": 0,
            "water": 0,
            "ions": 0,
            "other": 0,
        }

        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                # Parse atom information
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21].strip()

                # Determine restraint force based on atom type
                restraint_force, atom_type = self._get_restraint_force_detailed(
                    atom_name, residue_name, chain_id, constraints
                )

                # Update statistics
                if atom_type in restraint_stats:
                    restraint_stats[atom_type] += 1

                # Replace B-factor with restraint force (columns 61-66)
                new_line = line[:60] + f"{restraint_force:6.2f}" + line[66:]
                restraint_lines.append(new_line)
                atom_count += 1
            else:
                # Keep non-atom lines as is
                restraint_lines.append(line)

        # Write restraints file
        with open(output_file, "w") as f:
            f.writelines(restraint_lines)

        # Log statistics
        self.logger.info(f"Generated restraints file: {output_file}")
        self.logger.info(f"Source PDB: {system_pdb.name}")
        self.logger.info(f"Stage: {stage_name}")
        self.logger.info(f"Total atoms processed: {atom_count}")
        for atom_type, count in restraint_stats.items():
            if count > 0:
                force = constraints.get(atom_type, 0.0)
                self.logger.info(
                    f"  {atom_type}: {count} atoms, force = {force} kcal/mol/Å²"
                )

    def _get_restraint_force_detailed(
        self,
        atom_name: str,
        residue_name: str,
        chain_id: str,
        constraints: Dict[str, float],
    ) -> tuple:
        """
        Determine the appropriate restraint force for an atom with detailed classification.

        Args:
            atom_name: Name of the atom
            residue_name: Name of the residue
            chain_id: Chain identifier
            constraints: Dictionary of constraint forces

        Returns:
            Tuple of (restraint_force, atom_type)
        """

        # Standard protein residues (including protonation states from AMBER)
        protein_residues = {
            # Standard amino acids
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
            "HSE",
            "HSD",
            "HSP",
            "CYX",
            # Protonation states from AMBER/propka (based on PROTONATION_STATES dict)
            "ASH",  # Protonated aspartic acid
            "GLH",  # Protonated glutamic acid
            "HIE",  # Histidine with proton on epsilon nitrogen
            "HID",  # Histidine with proton on delta nitrogen
            "HIP",  # Histidine with both nitrogens protonated
            "LYN",  # Neutral lysine (deprotonated)
            "TYM",  # Deprotonated tyrosine
            "CYM",  # Deprotonated cysteine
            "SEP",  # Phosphorylated serine
            "T2P",  # Phosphorylated threonine
            # Terminal caps
            "ACE",
            "NHE",
            "NME",
            "COO",
        }

        # Lipid residues (include common AMBER/CHARMM lipid names)
        lipid_residues = {
            # Standard lipid names
            "POPC",
            "POPE",
            "POPS",
            "DPPC",
            "DMPC",
            "DOPC",
            "DSPC",
            "CHOL",
            "CHOLEST",
            "PALM",
            "OLEO",
            "STEROL",
            # AMBER-style lipid names (common in packmol-memgen)
            "PC",
            "PE",
            "PS",
            "PA",
            "PG",
            "PI",
            "SM",
            "CHL",
            "CHOL",
            "OL",
            "LA",
            "MY",
            "PA",
            "ST",
            "AR",  # AMBER lipid residue codes
            # Additional lipid variants
            "OLE",
            "PAL",
            "STE",
            "LIN",
            # Note: LYN removed - it's neutral lysine (protein), not a lipid
        }

        # Water residues (include various naming conventions)
        water_residues = {"TIP3", "HOH", "WAT", "SOL", "TIP4", "SPC", "T3P", "T4P"}

        # Ion residues (include various naming conventions)
        ion_residues = {
            "NA",
            "CL",
            "K",
            "CA",
            "MG",
            "ZN",
            "FE",
            "CU",
            "Na+",
            "Cl-",
            "K+",
            "Ca2+",
            "Mg2+",
            "Zn2+",
            "Fe2+",
            "Fe3+",
            "SOD",
            "CLA",
            "POT",
            "CAL",
            "MAG",
            "ZIN",
            "IRN",
            "COP",
        }

        # Protein atoms
        if residue_name in protein_residues:
            # Backbone atoms (including hydrogens)
            backbone_atoms = {
                "N",
                "CA",
                "C",
                "O",
                "OXT",
                "H",
                "HN",
                "HA",
                "HT1",
                "HT2",
                "HT3",
            }

            if atom_name in backbone_atoms:
                return constraints.get("protein_backbone", 0.0), "protein_backbone"
            else:
                return constraints.get("protein_sidechain", 0.0), "protein_sidechain"

        # Lipid molecules
        elif residue_name in lipid_residues:
            # Head group atoms (phosphate, choline, ethanolamine, etc.)
            head_atoms = {
                # Phosphate groups
                "P",
                "O11",
                "O12",
                "O13",
                "O14",
                "O21",
                "O22",
                "O31",
                "O32",
                "O33",
                "O34",
                "O1P",
                "O2P",
                "O3P",
                "O4P",
                "OP1",
                "OP2",
                "OP3",
                "OP4",
                # Choline and ethanolamine heads (specific to head group only)
                "N",
                "C11",
                "C12",
                "C13",
                "C14",  # Removed C15, C16 as they can be tail carbons
                "N31",
                "C32",
                "C33",
                "C34",
                "C35",
                # Glycerol backbone (connects head to tails)
                "C1",
                "C2",
                "C3",
                "O21",
                "O31",
                # Common head group patterns
                "HN1",
                "HN2",
                "HN3",
                "HO2",
                "HO3",
                "HS",
            }

            # Check atom name patterns for head vs tail classification
            if (
                atom_name in head_atoms
                or atom_name.startswith(("P", "O1", "O2", "N"))
                or "P" in atom_name
                or "N3" in atom_name
                or (
                    atom_name.startswith("C")
                    and len(atom_name) <= 2
                    and atom_name in ["C1", "C2", "C3"]
                )
            ):
                return constraints.get("lipid_head", 0.0), "lipid_head"
            else:
                # Tail carbons and other atoms (including C15, C16, etc.)
                return constraints.get("lipid_tail", 0.0), "lipid_tail"

        # Water molecules
        elif residue_name in water_residues:
            return constraints.get("water", 0.0), "water"

        # Ions
        elif residue_name in ion_residues:
            return constraints.get("ions", 0.0), "ions"

        # Other molecules (ligands, cofactors, etc.)
        else:
            return constraints.get("other", 0.0), "other"

    @staticmethod
    def get_default_stage_params(
        scheme_type: str = "NPT",
        temperature: float = 310.15,
        include_production: bool = False,
    ) -> List[Dict[str, Any]]:
        """
        Return default CHARMM-GUI-style equilibration stages for a membrane protein system.

        Six stages with gradually decreasing positional restraints, following the
        standard CHARMM-GUI membrane equilibration schedule. Suitable as a starting
        point that can be further customised before passing to setup_namd_equilibration.

        Args:
            scheme_type: Ensemble for all stages (NVT | NPT | NPAT | NPgT).
            temperature: Simulation temperature in Kelvin (default 310.15).
            include_production: When True, append a 50 ns unrestrained production
                stage (default False).

        Returns:
            List of :class:`EquilibrationStage` objects ready to pass to
            setup_namd_equilibration.  Fields can be edited via attribute
            assignment or the :meth:`~EquilibrationStage.replace` method.

        Example::

            >>> from dataclasses import replace
            >>> stages = NAMDEquilibrationManager.get_default_stage_params("NPT",
            ...                                                              include_production=True)
            >>> stages[-1].time_ns = 100.0          # mutable attribute set
            >>> stages[0] = stages[0].replace(temperature=303.15)  # immutable copy
            >>> manager = NAMDEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_namd_equilibration(stage_params_list=stages)
        """
        valid = {"NVT", "NPT", "NPAT", "NPgT"}
        if scheme_type not in valid:
            raise ValueError(
                f"scheme_type must be one of {sorted(valid)}, got '{scheme_type}'"
            )

        def _steps(time_ns: float, timestep_fs: float) -> int:
            return int(round(time_ns * 1_000_000 / timestep_fs))

        pressure = 1.0 if scheme_type in {"NPT", "NPAT", "NPgT"} else None
        surface_tension = 0.0 if scheme_type in {"NPAT", "NPgT"} else None

        def _stage(name, time_ns, timestep, minimize_steps=0, **constraints_overrides):
            base_constraints = {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "water": 0.0,
                "ions": 0.0,
                "other": 0.0,
            }
            base_constraints.update(constraints_overrides)
            return EquilibrationStage(
                name=name,
                ensemble=scheme_type,
                time_ns=time_ns,
                steps=_steps(time_ns, timestep),
                timestep=timestep,
                temperature=temperature,
                minimize_steps=minimize_steps,
                pressure=pressure,
                surface_tension=surface_tension,
                constraints=base_constraints,
            )

        stages: List[EquilibrationStage] = [
            _stage(
                "Equilibration 1",
                0.125,
                1.0,
                minimize_steps=10000,
                protein_backbone=10.0,
                protein_sidechain=5.0,
                lipid_head=2.5,
            ),
            _stage(
                "Equilibration 2",
                0.125,
                1.0,
                protein_backbone=5.0,
                protein_sidechain=2.5,
                lipid_head=1.0,
            ),
            _stage(
                "Equilibration 3",
                0.125,
                1.0,
                protein_backbone=2.5,
                protein_sidechain=1.0,
                lipid_head=0.5,
            ),
            _stage(
                "Equilibration 4",
                0.25,
                1.0,
                protein_backbone=1.0,
                protein_sidechain=0.5,
            ),
            _stage("Equilibration 5", 0.25, 2.0, protein_backbone=0.5),
            _stage("Equilibration 6", 0.5, 2.0, protein_backbone=0.1),
        ]

        if include_production:
            stages.append(_stage("Production", 50.0, 2.0))

        return stages

    def setup_namd_equilibration(
        self,
        system_files: Optional[Dict[str, str]] = None,
        stage_params_list: Optional[List[Dict[str, Any]]] = None,
        output_name: str = "equilibration",
        scheme_type: Optional[str] = None,
        namd_executable: str = "namd3",
        selections: Optional[Dict[str, str]] = None,
        add_com_restraint: bool = False,
        com_restraint_k: float = 10.0,
        add_rotation_restraint: bool = True,
        rotation_restraint_k: float = 2000.0,
        com_selection: str = "name CA",
        rotation_ref_positions_mode: str = "auto",
        ref_positions_file: Optional[str] = None,
        ref_positions_col: Optional[str] = None,
        ref_positions_col_value: Optional[float] = None,
        water_model: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Complete NAMD equilibration setup - replicates GUI workflow.

        This simplified method performs all steps needed for equilibration:
        1. Auto-detects system files (if not provided)
        2. Creates output directory structure (equilibration/namd/)
        3. Copies system files to NAMD directory
        4. Generates NAMD configuration files for all stages
        5. Generates restraint files for each stage
        6. Creates run script
        7. Creates protocol summary

        Args:
            system_files: Dictionary with source file paths (optional - will auto-detect if None):
                - 'prmtop': AMBER topology file (.prmtop)
                - 'inpcrd': AMBER coordinate file (.inpcrd)
                - 'pdb': System PDB file
                - 'bilayer_pdb': (REQUIRED) Bilayer PDB file with CRYST1 record for box dimensions
                If None, will search working_dir for standard file names.
            stage_params_list: List of stage dictionaries with equilibration parameters.
                Each stage must include 'ensemble' key (NVT, NPT, NPAT, or NPgT).
            output_name: Output directory name (default: "equilibration")
            scheme_type: Equilibration scheme (optional - auto-detected from stages).
                If None, will be extracted from the 'ensemble' field of the first stage.
                Can be explicitly set if needed (NVT, NPT, NPAT, or NPgT).
            namd_executable: NAMD executable path (default: "namd3")

        Returns:
            Dictionary with paths to generated files:
            {
                'output_dir': Path to main output directory,
                'namd_dir': Path to NAMD directory,
                'config_files': List of config file paths,
                'restraints_dir': Path to restraints directory,
                'run_script': Path to run script,
                'summary_file': Path to protocol summary
            }

        Example 1 (Auto-detect files):
            >>> manager = NAMDEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_namd_equilibration(
            ...     stage_params_list=[
            ...         {'name': 'Equilibration 1', 'time_ns': 0.125, ...}
            ...     ],
            ...     scheme_type="NPT"
            ... )

        Example 2 (Explicit file paths):
            >>> manager = NAMDEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_namd_equilibration(
            ...     system_files={
            ...         'prmtop': '/path/to/system.prmtop',
            ...         'inpcrd': '/path/to/system.inpcrd',
            ...         'pdb': '/path/to/system.pdb',
            ...         'bilayer_pdb': '/path/to/bilayer_lipid.pdb'
            ...     },
            ...     stage_params_list=[
            ...         {'name': 'Equilibration 1', 'time_ns': 0.125, ...}
            ...     ],
            ...     scheme_type="NPT"
            ... )
            >>> # Run: cd result['namd_dir'] && ./run_equilibration.sh
        """
        import shutil
        import json

        from gatewizard.tools.namd_water import (
            namd_water_model_config_block,
            normalize_water_model,
            read_water_model_from_builder_status,
        )

        if water_model is None:
            water_model = read_water_model_from_builder_status(self.working_dir)
        self.water_model = normalize_water_model(water_model or "tip3p")
        self.logger.info(f"NAMD water model for config: {self.water_model}")

        self.logger.info("=== Setting up NAMD equilibration ===")

        # Normalise EquilibrationStage objects to plain dicts early so that
        # all subsequent code can safely use dict-API (.get, subscript, etc.)
        if stage_params_list:
            stage_params_list = [
                s.to_dict() if isinstance(s, EquilibrationStage) else s
                for s in stage_params_list
            ]

        # Auto-detect scheme_type from stages if not provided
        if scheme_type is None:
            if stage_params_list and len(stage_params_list) > 0:
                # Extract ensemble from first stage
                first_ensemble = stage_params_list[0].get("ensemble", "NPT")
                scheme_type = first_ensemble
                self.logger.info(
                    f"Auto-detected scheme_type from stages: {scheme_type}"
                )
            else:
                # Default fallback
                scheme_type = "NPT"
                self.logger.info(
                    f"No stages provided, using default scheme_type: {scheme_type}"
                )

        # Validate scheme_type
        valid_schemes = ["NVT", "NPT", "NPAT", "NPgT"]
        if scheme_type not in valid_schemes:
            raise ValueError(
                f"Invalid scheme_type '{scheme_type}'. Must be one of {valid_schemes}"
            )

        # Auto-detect system files if not provided
        if system_files is None:
            self.logger.info("Auto-detecting system files in working directory...")
            system_files = self.find_system_files()
            if system_files is None:
                raise FileNotFoundError(
                    "Could not auto-detect required system files. "
                    "Please provide system_files dictionary explicitly."
                )

        # Validate required files exist
        for file_type in ["prmtop", "inpcrd", "pdb", "bilayer_pdb"]:
            if file_type not in system_files:
                raise ValueError(f"Missing required key '{file_type}' in system_files")
            file_path = Path(system_files[file_type])
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")

        # Create output directory structure
        # output_name may be an absolute Path (mirrors OpenMM behaviour)
        namd_dir = (
            Path(output_name)
            if Path(str(output_name)).is_absolute()
            else Path(self.working_dir) / output_name
        )
        output_dir = namd_dir
        restraints_dir = namd_dir / "restraints"

        for directory in [namd_dir, restraints_dir]:
            directory.mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Created directory: {directory}")

        # Step 1: Copy system files to NAMD directory
        self.logger.info("Copying system files to NAMD directory...")
        copied_files = {}

        for file_type, source_path in system_files.items():
            source_path_obj = Path(source_path)

            # Special handling for bilayer_pdb - keep original name, only used for CRYST1
            if file_type == "bilayer_pdb":
                target_filename = source_path_obj.name  # Keep original filename
                target_path = namd_dir / target_filename
                shutil.copy2(source_path, target_path)
                copied_files[file_type] = target_filename
                self.logger.info(
                    f"  Copied {file_type}: {source_path_obj.name} (for CRYST1 box info only)"
                )
            else:
                # Standard naming: system.prmtop, system.inpcrd, system.pdb
                target_filename = f"system{source_path_obj.suffix}"
                target_path = namd_dir / target_filename
                shutil.copy2(source_path, target_path)
                copied_files[file_type] = target_filename
                self.logger.info(
                    f"  Copied {file_type}: {source_path_obj.name} -> {target_filename}"
                )

        # Step 2: Generate NAMD configuration files for each stage
        self.logger.info("Generating NAMD configuration files...")
        config_files = []
        protocols_dict = {}

        # Use default stages when none are provided
        if not stage_params_list:
            stage_params_list = self.get_default_stage_params(scheme_type)
            self.logger.info(
                f"No stages provided — using default {scheme_type} protocol "
                f"({len(stage_params_list)} stages)"
            )

        # Normalise EquilibrationStage objects to plain dicts
        stage_params_list = [
            s.to_dict() if isinstance(s, EquilibrationStage) else s
            for s in stage_params_list
        ]

        # Convert stage list to protocols dictionary
        for i, stage_params in enumerate(stage_params_list):
            stage_name = stage_params.get("name", f"Equilibration {i+1}")
            protocols_dict[stage_name] = stage_params

        previous_stage_name = None
        for i, (stage_name, stage_params) in enumerate(protocols_dict.items()):
            # Generate config using CHARMM-GUI template system
            config_content = self.generate_charmm_gui_config_file(
                stage_name=stage_name,
                stage_params=stage_params,
                stage_index=i,
                system_files=copied_files,  # Use relative names
                scheme_type=scheme_type,
                previous_stage_name=previous_stage_name,
                all_stage_settings=protocols_dict,
            )

            # Write configuration file
            config_name = self._get_config_name(stage_name, i)
            if config_name == "step7_production":
                config_file = namd_dir / f"{config_name}.conf"
            else:
                config_file = namd_dir / f"{config_name}_equilibration.conf"

            config_file.write_text(config_content)
            config_files.append(config_file)
            self.logger.info(f"  Generated: {config_file.name}")

            previous_stage_name = stage_name

        # Step 3: Generate restraint files for each stage
        self.logger.info("Generating restraint files...")
        system_pdb = namd_dir / copied_files.get("pdb", "system.pdb")

        if system_pdb.exists():
            for i, (stage_name, stage_params) in enumerate(protocols_dict.items()):
                constraints = stage_params.get("constraints", {})
                has_restraints = any(float(v) > 0 for v in constraints.values())

                if has_restraints:
                    config_name = self._get_config_name(stage_name, i)
                    if config_name == "step7_production":
                        restraint_file = (
                            restraints_dir / f"{config_name}_restraints.pdb"
                        )
                    else:
                        restraint_file = (
                            restraints_dir
                            / f"{config_name}_equilibration_restraints.pdb"
                        )

                    self.generate_restraints_file(
                        system_pdb=system_pdb,
                        constraints=constraints,
                        output_file=restraint_file,
                        stage_name=stage_params.get("name", stage_name),
                        selections=selections,
                    )
                    self.logger.info(f"  Generated: {restraint_file.name}")
        else:
            self.logger.warning(
                f"System PDB not found: {system_pdb}, skipping restraints"
            )

        # Step 4: Generate run script
        self.logger.info("Generating run script...")
        run_script_content = self.generate_run_script(protocols_dict, namd_executable)
        run_script = namd_dir / "run_equilibration.sh"
        run_script.write_text(run_script_content)
        run_script.chmod(0o755)
        self.logger.info(f"  Generated: {run_script.name}")

        # Step 5: Create protocol summary
        protocol_summary = {
            "protocol_name": f"{scheme_type} Equilibration Protocol",
            "total_stages": len(protocols_dict),
            "scheme_type": scheme_type,
            "stages": protocols_dict,
            "namd_executable": namd_executable,
            "force_field": "AMBER",
        }

        summary_file = namd_dir / "protocol_summary.json"
        with open(summary_file, "w") as f:
            json.dump(protocol_summary, f, indent=2)
        self.logger.info(f"  Generated: {summary_file.name}")

        # --- COM restraint ---
        com_colvars_path: Optional[Path] = None
        if add_com_restraint:
            pdb_src = system_files.get("pdb") if system_files else None
            if pdb_src and Path(pdb_src).exists():
                dest_pdb = namd_dir / Path(pdb_src).name
                com_colvars_relpath = Path("restraints") / "com_restraint.col"
                com_colvars_path = self.generate_com_colvars_config(
                    pdb_path=dest_pdb if dest_pdb.exists() else Path(pdb_src),
                    output_file=namd_dir / com_colvars_relpath,
                    com_restraint_k=com_restraint_k,
                    add_rotation_restraint=add_rotation_restraint,
                    rotation_restraint_k=rotation_restraint_k,
                    selection=com_selection,
                    rotation_ref_positions_mode=rotation_ref_positions_mode,
                    ref_positions_file=ref_positions_file,
                    ref_positions_col=ref_positions_col,
                    ref_positions_col_value=ref_positions_col_value,
                )
                if com_colvars_path:
                    activation_block = _build_com_colvars_activation_block(
                        "namd", str(com_colvars_relpath)
                    )
                    for config_file in config_files:
                        config_text = config_file.read_text()
                        if "colvars on" not in config_text:
                            config_file.write_text(
                                config_text.rstrip() + activation_block
                            )
                    self.logger.info(
                        "  COM colvars file generated and activated in NAMD configs."
                    )
            else:
                self.logger.warning(
                    "No PDB in system_files; COM colvars file not generated."
                )

        self.logger.info("=== Setup complete ===")
        self.logger.info(f"Output directory: {namd_dir}")
        self.logger.info(f"To run: cd {namd_dir} && ./run_equilibration.sh")

        return {
            "output_dir": output_dir,
            "namd_dir": namd_dir,
            "config_files": config_files,
            "restraints_dir": restraints_dir,
            "run_script": run_script,
            "summary_file": summary_file,
            "com_colvars": com_colvars_path,
        }

    def generate_com_colvars_config(
        self,
        pdb_path: Path,
        output_file: Path,
        com_restraint_k: float = 10.0,
        add_rotation_restraint: bool = True,
        rotation_restraint_k: float = 2000.0,
        selection: str = "name CA",
        rotation_ref_positions_mode: str = "auto",
        ref_positions_file: Optional[str] = None,
        ref_positions_col: Optional[str] = None,
        ref_positions_col_value: Optional[float] = None,
    ) -> Optional[Path]:
        """Generate NAMD Colvars configuration to restrain protein centre of mass.

        Prevents protein translation (and optionally rotation) by applying
        harmonic restraints on the **geometric centre** of selected atoms —
        a centre-of-mass / centre-of-geometry restraint that introduces no
        per-atom positional bias.

        Three 1-D CVs (``distanceX``, ``distanceY``, ``distanceZ``) are used
        for translation; an ``orientation`` CV is added when
        *add_rotation_restraint* is True.  The initial centroid is computed
        from *pdb_path* so that the rest position matches the starting
        structure.

        Include the generated file in NAMD via::

            colvarsConfig com_restraint.col

        Args:
            pdb_path: System PDB file (used to compute the initial centroid).
            output_file: Destination for the Colvars ``.col`` config.
            com_restraint_k: Translation force constant in kcal/mol/Å².
            add_rotation_restraint: When True, also add an orientation CV.
            rotation_restraint_k: Rotation force constant in kcal/mol/Å².
            selection: MDAnalysis selection string for the reference atoms
                (default: ``"name CA"`` for Cα atoms).
            rotation_ref_positions_mode: Orientation reference mode:
                ``"auto"``, ``"refPositions"``, or ``"refPositionsFile"``.
                ``"auto"`` defaults to ``"refPositionsFile"`` for NAMD.
            ref_positions_col: Optional PDB column (``O/B/X/Y/Z``) to select
                reference atoms from ``refPositionsFile``.
            ref_positions_col_value: Optional numeric value paired with
                ``ref_positions_col``.

        Returns:
            Path to the written config file, or ``None`` if MDAnalysis is
            unavailable or the selection matches no atoms.
        """
        try:
            import MDAnalysis as mda  # type: ignore

            u = mda.Universe(str(pdb_path))
            ag = u.select_atoms(selection)
            atoms = getattr(ag, "atoms", ag)
            if len(atoms) == 0:
                self.logger.warning(
                    f"No atoms matched '{selection}' in {pdb_path.name}; "
                    "COM colvars file not generated."
                )
                return None

            com = ag.center_of_geometry()
            x0, y0, z0 = float(com[0]), float(com[1]), float(com[2])
            atom_nums = " ".join(str(int(a.index) + 1) for a in atoms)

            content = _build_com_colvars_config(
                atom_numbers=atom_nums,
                x0=x0,
                y0=y0,
                z0=z0,
                com_k=com_restraint_k,
                add_rotation=add_rotation_restraint,
                rot_k=rotation_restraint_k,
                ag=ag,
                engine="namd",
                ref_positions_file=ref_positions_file or Path(pdb_path).name,
                rotation_ref_positions_mode=rotation_ref_positions_mode,
                ref_positions_col=ref_positions_col,
                ref_positions_col_value=ref_positions_col_value,
            )

            output_file.write_text(content)
            self.logger.info(
                f"  COM colvars config ({len(atoms)} atoms, centroid "
                f"[{x0:.2f},{y0:.2f},{z0:.2f}] Å): {output_file.name}"
            )
            return output_file

        except ImportError:
            self.logger.warning(
                "MDAnalysis not available; COM colvars file not generated. "
            )
            return None
        except Exception as exc:
            self.logger.error(f"COM colvars generation failed: {exc}")
            return None

    def generate_colvars_file(
        self,
        system_pdb: Path,
        output_file: Path,
        stage_params: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Generate NAMD colvars configuration file for bilayer thickness restraint using bilayer_utils.

        Args:
            system_pdb: Path to the system PDB file to analyze for phosphate atoms
            output_file: Path for output colvars configuration file (should be bilayer_thickness.col)
            stage_params: Stage parameters including bilayer thickness and force constant
        """
        try:
            if stage_params and stage_params.get("bilayer_thickness") is not None:
                self._generate_bilayer_thickness_restraint(
                    system_pdb, output_file, stage_params
                )
            else:
                self.logger.info(
                    "Bilayer thickness restraint not enabled for this stage"
                )

        except Exception as e:
            self.logger.error(f"Error generating colvars file: {e}")
            self.logger.info("Falling back to basic approach")
            # Ensure stage_params is not None for fallback
            if stage_params is None:
                stage_params = {}
            self._generate_bilayer_thickness_restraint_fallback(
                system_pdb, output_file, stage_params
            )

    def _generate_bilayer_thickness_restraint(
        self,
        system_pdb: Path,
        output_file: Path,
        stage_params: Optional[Dict[str, Any]],
    ) -> None:
        """
        Generate colvars configuration for bilayer thickness harmonic restraint using bilayer_utils.
        This creates a simple harmonic restraint to maintain bilayer thickness.
        """
        try:
            # Import bilayer utilities
            import sys
            import os

            # Add the gatewizard utils path to sys.path
            gatewizard_root = os.path.dirname(
                os.path.dirname(os.path.dirname(__file__))
            )
            utils_path = os.path.join(gatewizard_root, "utils")
            if utils_path not in sys.path:
                sys.path.insert(0, utils_path)

            from bilayer_utils import BilayerAnalyzer  # type: ignore

            # Initialize bilayer analyzer with all phosphate patterns (including P31)
            analyzer = BilayerAnalyzer()

            # Analyze bilayer to find ALL phosphate atoms
            upper_bilayer, lower_bilayer = analyzer.analyze_bilayer_from_pdb(
                str(system_pdb)
            )

            if not upper_bilayer or not lower_bilayer:
                self.logger.warning(
                    "Could not find phosphate atoms in both bilayers for thickness restraint"
                )
                return

            # Get statistics for logging
            stats = analyzer.get_bilayer_statistics(upper_bilayer, lower_bilayer)
            self.logger.info(
                f"Found {stats['total_phosphorus_atoms']} total phosphorus atoms for bilayer thickness restraint"
            )
            self.logger.info(f"Upper bilayer: {stats['upper_bilayer_count']} atoms")
            self.logger.info(f"Lower bilayer: {stats['lower_bilayer_count']} atoms")
            self.logger.info(f"Atom types found: {', '.join(stats['atom_types'])}")
            self.logger.info(f"Residue types: {', '.join(stats['residue_types'])}")

            # Get parameters from stage_params or use defaults
            if stage_params is None:
                stage_params = {}
            target_thickness = float(stage_params.get("bilayer_thickness", 39.1))  # Å
            force_constant = float(
                stage_params.get("force_constant", 10.0)
            )  # kcal/mol/Å²

            # Generate simple harmonic restraint configuration
            colvars_content = self._create_harmonic_thickness_config(
                upper_bilayer, lower_bilayer, target_thickness, force_constant
            )

            # Write colvars file
            with open(output_file, "w") as f:
                f.write(colvars_content)

            self.logger.info(
                f"Generated bilayer thickness restraint file with {stats['total_phosphorus_atoms']} phosphate atoms: {output_file}"
            )
            self.logger.info(
                f"Target thickness: {target_thickness} Å, Force constant: {force_constant} kcal/mol/Å²"
            )

        except ImportError as e:
            self.logger.error(f"Could not import bilayer_utils: {e}")
            self.logger.info("Falling back to basic approach")
            if stage_params is None:
                stage_params = {}
            self._generate_bilayer_thickness_restraint_fallback(
                system_pdb, output_file, stage_params
            )
        except Exception as e:
            self.logger.error(f"Error in bilayer analysis: {e}")
            self.logger.info("Falling back to basic approach")
            if stage_params is None:
                stage_params = {}
            self._generate_bilayer_thickness_restraint_fallback(
                system_pdb, output_file, stage_params
            )

    def _create_harmonic_thickness_config(
        self,
        upper_bilayer,
        lower_bilayer,
        target_thickness: float,
        force_constant: float,
    ) -> str:
        """
        Create simple harmonic restraint colvars configuration for bilayer thickness.

        Args:
            upper_bilayer: List of PhosphorusAtom objects for upper bilayer
            lower_bilayer: List of PhosphorusAtom objects for lower bilayer
            target_thickness: Target bilayer thickness in Å
            force_constant: Force constant in kcal/mol/Å²

        Returns:
            Complete colvars configuration string with harmonic restraint
        """
        # Extract ALL NAMD indices (0-based) - no limit, include all atoms
        upper_indices = [str(atom.namd_index) for atom in upper_bilayer]
        lower_indices = [str(atom.namd_index) for atom in lower_bilayer]

        # Format indices in readable chunks (10 per line for readability)
        def format_indices_multiline(
            indices: List[str], indent: str = "            "
        ) -> str:
            if not indices:
                return ""

            lines = []
            for i in range(0, len(indices), 10):
                chunk = indices[i : i + 10]
                lines.append(indent + " ".join(chunk))
            return "\n".join(lines)

        upper_formatted = format_indices_multiline(upper_indices)
        lower_formatted = format_indices_multiline(lower_indices)

        # Create simple harmonic restraint configuration
        config = f"""# NAMD Colvars Configuration for Bilayer Thickness Restraint
# Generated by Gatewizard using bilayer analysis
# Upper bilayer atoms: {len(upper_bilayer)} phosphates
# Lower bilayer atoms: {len(lower_bilayer)} phosphates
# Total atoms used: {len(upper_bilayer) + len(lower_bilayer)} phosphates

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNumbers {' '.join(upper_indices)}
        }}
        group2 {{
            atomNumbers {' '.join(lower_indices)}
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant {force_constant}  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""

        return config

    def _generate_bilayer_thickness_restraint_fallback(
        self, system_pdb: Path, output_file: Path, stage_params: Dict[str, Any]
    ) -> None:
        """
        Fallback method using basic phosphate detection for bilayer thickness restraint.
        Used when bilayer_utils is not available.
        """
        self.logger.warning(
            "Using fallback phosphate detection method for bilayer thickness restraint"
        )

        # Find phosphate atoms in the system using basic method
        phosphate_atoms = self._find_phosphate_atoms_basic(system_pdb)

        if len(phosphate_atoms) < 2:
            self.logger.warning(
                "Insufficient phosphate atoms found for bilayer thickness restraint"
            )
            return

        # Separate into upper and lower leaflets (simple z-coordinate based)
        z_coords = [atom["z"] for atom in phosphate_atoms]
        z_center = sum(z_coords) / len(z_coords)

        upper_leaflet = [atom for atom in phosphate_atoms if atom["z"] > z_center]
        lower_leaflet = [atom for atom in phosphate_atoms if atom["z"] < z_center]

        if len(upper_leaflet) == 0 or len(lower_leaflet) == 0:
            self.logger.warning("Could not separate phosphates into leaflets")
            return

        # Get parameters from stage_params or use defaults
        target_thickness = (
            float(stage_params.get("bilayer_thickness", 39.1)) if stage_params else 39.1
        )  # Å
        force_constant = (
            float(stage_params.get("force_constant", 10.0)) if stage_params else 10.0
        )  # kcal/mol/Å²

        # Generate basic harmonic restraint configuration
        colvars_content = self._create_basic_harmonic_thickness_config(
            upper_leaflet, lower_leaflet, target_thickness, force_constant
        )

        # Write colvars file
        with open(output_file, "w") as f:
            f.write(colvars_content)

        self.logger.info(
            f"Generated fallback bilayer thickness restraint file: {output_file}"
        )
        self.logger.info(
            f"Found {len(upper_leaflet)} upper and {len(lower_leaflet)} lower leaflet phosphates"
        )
        self.logger.info(
            f"Target thickness: {target_thickness} Å, Force constant: {force_constant} kcal/mol/Å²"
        )

    def _generate_phosphate_distance_colvars_fallback(
        self, system_pdb: Path, output_file: Path
    ) -> None:
        """
        Fallback method using the original simple phosphate detection.
        Used when bilayer_utils is not available.
        """
        self.logger.warning("Using fallback phosphate detection method")

        # Find phosphate atoms in the system using basic method
        phosphate_atoms = self._find_phosphate_atoms_basic(system_pdb)

        if len(phosphate_atoms) < 2:
            self.logger.warning(
                "Insufficient phosphate atoms found for distance measurement"
            )
            return

        # Separate into upper and lower leaflets (simple z-coordinate based)
        z_coords = [atom["z"] for atom in phosphate_atoms]
        z_center = sum(z_coords) / len(z_coords)

        upper_leaflet = [atom for atom in phosphate_atoms if atom["z"] > z_center]
        lower_leaflet = [atom for atom in phosphate_atoms if atom["z"] < z_center]

        if len(upper_leaflet) == 0 or len(lower_leaflet) == 0:
            self.logger.warning("Could not separate phosphates into leaflets")
            return

        # Generate basic colvars configuration
        colvars_content = self._create_basic_phosphate_distance_config(
            upper_leaflet, lower_leaflet
        )

        # Write colvars file
        with open(output_file, "w") as f:
            f.write(colvars_content)

        self.logger.info(f"Generated fallback colvars file: {output_file}")
        self.logger.info(
            f"Found {len(upper_leaflet)} upper and {len(lower_leaflet)} lower leaflet phosphates"
        )

    def _find_phosphate_atoms_basic(self, system_pdb: Path) -> List[Dict]:
        """Find phosphate atoms (P atoms) in lipid molecules."""
        phosphate_atoms = []

        with open(system_pdb, "r") as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith(("ATOM", "HETATM")):
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()

                    # Look for phosphorus atoms in lipid residues (3-character names)
                    # Support both 'P' and 'P31' atom names
                    if atom_name in ["P", "P31"] and residue_name in [
                        "PC",
                        "PA",
                        "PE",
                        "PS",
                        "PG",
                        "PI",
                    ]:
                        try:
                            atom_id = int(line[6:11].strip())
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())

                            phosphate_atoms.append(
                                {
                                    "id": atom_id,
                                    "name": atom_name,
                                    "residue": residue_name,
                                    "x": x,
                                    "y": y,
                                    "z": z,
                                    "line_num": line_num,
                                }
                            )
                        except (ValueError, IndexError) as e:
                            self.logger.warning(f"Error parsing line {line_num}: {e}")

        return phosphate_atoms

    def _create_basic_phosphate_distance_config(
        self, upper_leaflet: List[Dict], lower_leaflet: List[Dict]
    ) -> str:
        """Create basic colvars configuration for phosphate distance measurement (fallback method)."""

        # Use ALL atoms from each leaflet (no artificial limit)
        upper_atoms = " ".join([str(atom["id"]) for atom in upper_leaflet])
        lower_atoms = " ".join([str(atom["id"]) for atom in lower_leaflet])

        colvars_config = f"""# NAMD Colvars Configuration for Phosphate Distance Measurement (Fallback)
# Generated for ABF simulation of bilayer thickness
# Upper leaflet atoms: {len(upper_leaflet)} phosphates
# Lower leaflet atoms: {len(lower_leaflet)} phosphates

colvarsTrajFrequency    1000
colvarsRestartFrequency 1000

# Define collective variable for phosphate-phosphate distance
colvar {{
    name phosphate_distance
    
    # Distance between center of mass of upper and lower leaflet phosphates
    distance {{
        group1 {{
            atomNumbers {upper_atoms}
        }}
        group2 {{
            atomNumbers {lower_atoms}
        }}
    }}
}}

# ABF bias for phosphate distance
abf {{
    colvars          phosphate_distance
    fullSamples      200
    historyFreq      1000
    inputPrefix      ""
    outputPrefix     "phosphate_distances"
}}

# Metadynamics hills for enhanced sampling (optional)
metadynamics {{
    colvars         phosphate_distance
    hillWeight      0.1
    hillWidth       0.5
    newHillFreq     1000
    writeHillsFreq  1000
    outputPrefix    "phosphate_hills"
}}
"""
        return colvars_config

    def _create_basic_harmonic_thickness_config(
        self,
        upper_leaflet: List[Dict],
        lower_leaflet: List[Dict],
        target_thickness: float,
        force_constant: float,
    ) -> str:
        """Create basic harmonic restraint configuration for bilayer thickness (fallback method)."""

        # Use ALL atoms from each leaflet (no artificial limit)
        upper_atoms = " ".join([str(atom["id"]) for atom in upper_leaflet])
        lower_atoms = " ".join([str(atom["id"]) for atom in lower_leaflet])

        colvars_config = f"""# NAMD Colvars Configuration for Bilayer Thickness Restraint (Fallback)
# Generated for harmonic restraint of bilayer thickness
# Upper leaflet atoms: {len(upper_leaflet)} phosphates
# Lower leaflet atoms: {len(lower_leaflet)} phosphates

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNumbers {upper_atoms}
        }}
        group2 {{
            atomNumbers {lower_atoms}
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant {force_constant}  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
        return colvars_config

    def generate_run_script(
        self, protocols: Dict[str, Dict], namd_executable: Optional[str] = None
    ) -> str:
        """
        Generate bash script to run all equilibration stages.
        Each stage uses its own CPU/GPU settings from the protocol configuration.

        Args:
            protocols: Dictionary of equilibration protocols with per-stage resource settings
            namd_executable: Path to NAMD executable

        Returns:
            Content of the run script
        """

        namd_exe = namd_executable or self.namd_executable

        script_lines = [
            "#!/bin/bash",
            "#############################################################",
            "## NAMD Equilibration Run Script",
            "## Generated by Gatewizard",
            "#############################################################",
            "",
            "# Set NAMD executable",
            f'NAMD="{namd_exe}"',
            "",
            "# Check if NAMD is available",
            "if ! command -v $NAMD &> /dev/null; then",
            '    echo "Error: NAMD executable not found: $NAMD"',
            "    exit 1",
            "fi",
            "",
            'echo "Starting NAMD equilibration protocol..."',
            'echo "Each stage uses individual CPU/GPU settings"',
            'echo ""',
            "",
            NAMD_RESUME_SHELL,
            "",
        ]

        # Add commands for each stage
        for i, (stage_key, stage_data) in enumerate(protocols.items()):
            stage_num = i + 1
            stage_name = stage_data.get("name", stage_key)
            steps = stage_data.get("steps", "N/A")
            timestep = stage_data.get("timestep", "N/A")
            use_gpu = stage_data.get("use_gpu", False)
            cpu_cores = stage_data.get("cpu_cores", 1)
            gpu_id = stage_data.get("gpu_id", 0)
            num_gpus = stage_data.get("num_gpus", 1)

            # Build NAMD command with appropriate flags
            namd_cmd = f"$NAMD"

            # Add processor specification
            namd_cmd += f" +p{cpu_cores}"

            # Add GPU flags if enabled
            if use_gpu:
                if num_gpus == 1:
                    namd_cmd += f" +devices {gpu_id}"  # Single GPU device
                else:
                    # Multiple GPUs: create device list starting from gpu_id
                    device_list = ",".join(str(gpu_id + i) for i in range(num_gpus))
                    namd_cmd += f" +devices {device_list}"

            # Complete command - use config-safe names for file names with new step naming
            config_name = self._get_config_name(stage_key, i)
            if config_name == "step7_production":
                namd_cmd += f" step7_production.conf > step7_production.log 2>&1"
                namd_stem = "step7_production"
            else:
                namd_cmd += f" {config_name}_equilibration.conf > {config_name}_equilibration.log 2>&1"
                namd_stem = f"{config_name}_equilibration"

            # Create detailed resource information
            gpu_info = "No"
            if use_gpu:
                if num_gpus == 1:
                    gpu_info = f"Yes (GPU {gpu_id})"
                else:
                    gpu_list = ",".join(str(gpu_id + i) for i in range(num_gpus))
                    gpu_info = f"Yes ({num_gpus} GPUs: {gpu_list})"

            script_lines.extend(
                [
                    f"# Stage {stage_num}: {stage_name}",
                    f'if [ "$RESUME" = "1" ] && _gw_namd_stage_done "{namd_stem}"; then',
                    f'  echo "RESUME: skipping stage {stage_num} ({namd_stem})"',
                    "else",
                    f'  echo "Running Stage {stage_num}: {stage_name}"',
                    f'  echo "Steps: {steps}, Timestep: {timestep} ps"',
                    f'  echo "Resources: {cpu_cores} CPU cores, GPU: {gpu_info}"',
                    f"  {namd_cmd}",
                    "",
                    "  if [ $? -ne 0 ]; then",
                    f'    echo "Error in Stage {stage_num}: {stage_name}"',
                    "    exit 1",
                    "  fi",
                    f'  echo "Stage {stage_num} completed successfully"',
                    "fi",
                    'echo ""',
                    "",
                ]
            )

        script_lines.extend(
            [
                'echo "All equilibration stages completed successfully!"',
                'echo "Check the log files for detailed output"',
            ]
        )

        return "\n".join(script_lines)

    def run_equilibration(
        self, config_files: List[Path], num_processors: int = 4
    ) -> subprocess.Popen:
        """
        Run NAMD equilibration simulation.

        Args:
            config_files: List of NAMD configuration files to run
            num_processors: Number of processors to use

        Returns:
            Process object for the running simulation
        """

        # Create run script
        script_content = self._create_run_script(config_files, num_processors)

        # Write script to temporary file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".sh", delete=False) as f:
            f.write(script_content)
            script_path = f.name

        # Make script executable
        os.chmod(script_path, 0o755)

        # Run script
        process = subprocess.Popen(
            ["bash", script_path],
            cwd=self.working_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        self.logger.info(f"Started NAMD equilibration with PID: {process.pid}")
        return process

    def _create_run_script(self, config_files: List[Path], num_processors: int) -> str:
        """Create a run script for the given configuration files."""

        script_lines = [
            "#!/bin/bash",
            f"# NAMD Equilibration Script",
            f"# Generated by Gatewizard",
            "",
            f'NAMD="{self.namd_executable}"',
            f"NPROCS={num_processors}",
            "",
        ]

        for i, config_file in enumerate(config_files):
            stage_num = i + 1
            script_lines.extend(
                [
                    f'echo "Running stage {stage_num}: {config_file.name}"',
                    f"$NAMD +p$NPROCS {config_file.name}",
                    f'if [ $? -ne 0 ]; then echo "Error in stage {stage_num}"; exit 1; fi',
                    "",
                ]
            )

        script_lines.append('echo "Equilibration completed successfully!"')

        return "\n".join(script_lines)

    def load_charmm_gui_template(
        self,
        scheme_type: str,
        stage_number: int,
        system_files: Dict[str, str],
        target_thickness: Optional[float] = None,
    ) -> str:
        """
        Load and customize CHARMM-GUI template file for a specific scheme and stage.

        Args:
            scheme_type: Type of scheme (NVT, NPT, NPAT, NPgT)
            stage_number: Stage number (1-12 for equilibration, 13 for production)
            system_files: Dictionary of system file paths
            target_thickness: Target bilayer thickness (deprecated, not used anymore)

        Returns:
            Customized NAMD configuration content
        """
        # Map scheme types to folder names
        scheme_folders = {
            "NVT": "01_NVT",
            "NPT": "02_NPT",
            "NPAT": "03_NPAT",
            "NPgT": "04_NPgT",
        }

        if scheme_type not in scheme_folders:
            raise ValueError(f"Unknown scheme type: {scheme_type}")

        # Build template file path
        scheme_folder = scheme_folders[scheme_type]
        if stage_number <= 12:
            # Use stage 6 template for all equilibration stages 7-12
            if stage_number <= 6:
                template_file = f"step6.{stage_number}_equilibration.inp"
            else:
                template_file = (
                    "step6.6_equilibration.inp"  # Reuse last equilibration template
                )
        elif stage_number == 13:
            template_file = "step7_production.inp"
        else:
            raise ValueError(f"Invalid stage number: {stage_number}")

        template_path = self.namd_templates_dir / scheme_folder / template_file

        if not template_path.exists():
            raise FileNotFoundError(f"CHARMM-GUI template not found: {template_path}")

        # Read template content
        with open(template_path, "r") as f:
            template_content = f.read()

        # Customize template for Gatewizard
        customized_content = self._customize_charmm_gui_template_old(
            template_content, system_files, target_thickness
        )

        return customized_content

    def _customize_charmm_gui_template_old(
        self,
        template_content: str,
        system_files: Dict[str, str],
        target_thickness: Optional[float] = None,
    ) -> str:
        """
        (DEPRECATED - kept for backward compatibility)
        Customize CHARMM-GUI template content for Gatewizard.

        Args:
            template_content: Original template content
            system_files: Dictionary of system file paths
            target_thickness: Target bilayer thickness for restraints

        Returns:
            Customized template content
        """
        lines = template_content.split("\n")
        customized_lines = []

        skip_restraint_section = False

        for line in lines:
            # Replace system file paths (use paths as-is from system_files, no hardcoded ../../)
            # Match both CHARMM-GUI format (step5_input.*) and Gatewizard format (system.*)
            if line.strip().startswith("parmfile") and (
                "step5_input.parm7" in line or "system.prmtop" in line
            ):
                customized_lines.append(
                    f"parmfile                {system_files.get('prmtop', 'system.prmtop')}"
                )
            elif line.strip().startswith("ambercoor") and (
                "step5_input.rst7" in line or "system.inpcrd" in line
            ):
                customized_lines.append(
                    f"ambercoor               {system_files.get('inpcrd', 'system.inpcrd')}"
                )

            # Skip planar and dihedral restraints as requested
            elif (
                "planar restraint" in line.lower()
                or "dihedral restraint" in line.lower()
            ):
                skip_restraint_section = True
                customized_lines.append(f"# {line.strip()} - DISABLED BY GATEWIZARD")
                continue
            elif skip_restraint_section and (
                line.strip().startswith("#") or line.strip() == ""
            ):
                customized_lines.append(line)
                continue
            elif skip_restraint_section:
                if any(
                    keyword in line.lower()
                    for keyword in ["colvars", "extrabonds", "exec sed"]
                ):
                    customized_lines.append(
                        f"# {line.strip()} - DISABLED BY GATEWIZARD"
                    )
                    continue
                else:
                    skip_restraint_section = False

            # Bilayer thickness restraint is now handled in _generate_restraints_block
            else:
                customized_lines.append(line)

        return "\n".join(customized_lines)

    def generate_bilayer_thickness_colvar(
        self,
        target_thickness: float,
        output_path: Path,
        pdb_file: Optional[Path] = None,
    ) -> None:
        """
        Generate collective variable file for bilayer thickness restraint.

        Args:
            target_thickness: Target bilayer thickness in Angstroms
            output_path: Path where to save the colvar file
            pdb_file: Optional PDB file for automatic phosphorus atom detection
        """
        # Try to use the new bilayer analyzer if PDB file is provided
        if pdb_file and pdb_file.exists():
            try:
                # Import the bilayer utilities
                import sys
                import os

                # Add the gatewizard utils path to sys.path
                gatewizard_root = os.path.dirname(
                    os.path.dirname(os.path.dirname(__file__))
                )
                utils_path = os.path.join(gatewizard_root, "utils")
                if utils_path not in sys.path:
                    sys.path.insert(0, utils_path)

                from bilayer_utils import generate_bilayer_thickness_colvar_from_pdb  # type: ignore

                # Generate colvar configuration using automatic analysis
                colvar_config = generate_bilayer_thickness_colvar_from_pdb(
                    str(pdb_file), colvar_name="bilayer_thickness"
                )

                # Add header comment and restraint configuration
                colvar_content = f"""# Collective variable for bilayer thickness control
# Generated by Gatewizard using automatic bilayer analysis
# PDB file: {pdb_file}

{colvar_config}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant 10.0  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""

                self.logger.info(
                    f"Generated bilayer thickness colvar using automatic analysis from {pdb_file}"
                )

            except (ImportError, Exception) as e:
                self.logger.warning(f"Could not use automatic bilayer analysis: {e}")
                self.logger.info("Falling back to generic colvar configuration")

                # Fall back to the old generic method
                colvar_content = f"""# Collective variable for bilayer thickness control
# Generated by Gatewizard (generic configuration)

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNameResidueRange P 1-999999
        }}
        group2 {{
            atomNameResidueRange P 1-999999  
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant 10.0  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
        else:
            # No PDB file provided or doesn't exist, use generic method
            self.logger.info(
                "No PDB file provided for automatic analysis, using generic colvar configuration"
            )
            colvar_content = f"""# Collective variable for bilayer thickness control
# Generated by Gatewizard (generic configuration)

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNameResidueRange P 1-999999
            atomNameResidueRange P 1-999999
        }}
        group2 {{
            atomNameResidueRange P 1-999999  
            atomNameResidueRange P 1-999999
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant 10.0  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""

        # Ensure output directory exists
        from gatewizard.utils.helpers import create_directory_robust

        create_directory_robust(output_path.parent)

        # Write colvar file
        with open(output_path, "w") as f:
            f.write(colvar_content)

        self.logger.info(f"Generated bilayer thickness colvar file: {output_path}")

    def generate_charmm_gui_config_file(
        self,
        stage_name: str,
        stage_params: Dict[str, Any],
        stage_index: int,
        system_files: Dict[str, str],
        scheme_type: str,
        previous_stage_name: Optional[str] = None,
        all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None,
        force_scheme_type: bool = False,
    ) -> str:
        """
        Generate NAMD configuration file using CHARMM-GUI templates with GateWizard customizations.

        Args:
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths
            scheme_type: CHARMM-GUI scheme type (NVT, NPT, NPAT, NPgT) - default ensemble for protocol
            previous_stage_name: Name of the previous stage for restart files
            force_scheme_type: If True, always use scheme_type for template selection,
                             ignoring stage-specific ensemble values (GUI mode)

        Returns:
            Content of the NAMD configuration file
        """
        # Skip minimization stage - it's now incorporated into the first equilibration
        if stage_name == "minimization":
            self.logger.info(
                "Skipping separate minimization stage - now included in eq1_equilibration"
            )
            return ""

        # Map stage names to template files using config name mapping with stage_index
        config_name = self._get_config_name(stage_name, stage_index)

        # Determine which ensemble/scheme to use for template selection
        # Priority: 1) custom_template key, 2) stage's ensemble key (unless forced), 3) global scheme_type
        if force_scheme_type:
            # GUI mode: always use the selected scheme_type for all stages
            stage_ensemble = scheme_type
        else:
            # API mode: allow per-stage ensemble customization
            stage_ensemble = stage_params.get("ensemble") or scheme_type

        if isinstance(stage_ensemble, str):
            stage_ensemble = stage_ensemble.upper()

        # Check if user explicitly specified a custom template
        custom_template = stage_params.get("custom_template", None)

        if custom_template:
            # User explicitly specified which template to use (e.g., 'step6.3_equilibration.inp')
            template_filename = custom_template
            template_scheme = (
                stage_ensemble  # Use stage's ensemble for folder selection
            )
            self.logger.info(
                f"Using custom template for {stage_name}: {custom_template} from {stage_ensemble} ensemble"
            )
        else:
            # Auto-select template based on stage index and ensemble
            # Define template mapping based on config names (6 equilibration stages + production)
            template_mapping = {
                "step1": "step6.1_equilibration.inp",
                "step2": "step6.2_equilibration.inp",
                "step3": "step6.3_equilibration.inp",
                "step4": "step6.4_equilibration.inp",
                "step5": "step6.5_equilibration.inp",
                "step6": "step6.6_equilibration.inp",
                "step7_production": "step7_production.inp",
            }

            # Get template filename based on stage position
            template_filename = template_mapping.get(
                config_name, "step6.1_equilibration.inp"
            )

            # Determine which ensemble scheme to use for template folder selection
            explicit_stage_ensemble = stage_params.get("ensemble")
            if explicit_stage_ensemble and stage_ensemble != scheme_type:
                # Stage uses different ensemble than the protocol default
                self.logger.warning(
                    f"Stage {stage_index + 1} ({stage_name}) uses ensemble '{stage_ensemble}' "
                    f"but protocol default is '{scheme_type}'. Using '{stage_ensemble}' template."
                )
                template_scheme = stage_ensemble
            else:
                template_scheme = scheme_type

        # Load and customize template
        return self._load_and_customize_charmm_gui_template(
            template_scheme,
            template_filename,
            stage_name,
            stage_params,
            stage_index,
            system_files,
            previous_stage_name,
            all_stage_settings,
        )

    def _load_and_customize_charmm_gui_template(
        self,
        scheme_type: str,
        template_filename: str,
        stage_name: str,
        stage_params: Dict[str, Any],
        stage_index: int,
        system_files: Dict[str, str],
        previous_stage_name: Optional[str] = None,
        all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> str:
        """
        Load CHARMM-GUI template and customize with GateWizard parameters.

        Args:
            scheme_type: CHARMM-GUI scheme type (NVT, NPT, NPAT, NPgT)
            template_filename: Template file to load
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths
            previous_stage_name: Name of the previous stage for restart files

        Returns:
            Customized NAMD configuration content
        """
        # Map scheme types to template directories
        scheme_mapping = {
            "NVT": "01_NVT",
            "NPT": "02_NPT",
            "NPAT": "03_NPAT",
            "NPgT": "04_NPgT",
        }

        scheme_folder = scheme_mapping.get(scheme_type, "01_NVT")
        template_path = self.namd_templates_dir / scheme_folder / template_filename

        if not template_path.exists():
            raise FileNotFoundError(f"CHARMM-GUI template not found: {template_path}")

        # Read template content
        with open(template_path, "r") as f:
            template_content = f.read()

        # Customize template with GateWizard parameters
        return self._customize_charmm_gui_template(
            template_content,
            stage_name,
            stage_params,
            stage_index,
            system_files,
            previous_stage_name,
            all_stage_settings,
        )

    def _customize_charmm_gui_template(
        self,
        template_content: str,
        stage_name: str,
        stage_params: Dict[str, Any],
        stage_index: int,
        system_files: Dict[str, str],
        previous_stage_name: Optional[str] = None,
        all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> str:
        """
        Customize CHARMM-GUI template content with GateWizard parameters.

        Args:
            template_content: Raw template content
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths
            previous_stage_name: Name of the previous stage for restart files

        Returns:
            Customized configuration content
        """
        # Temperature
        temperature = stage_params.get("temperature", 303.15)

        # DCD frequency for trajectory output
        dcd_freq = stage_params.get("dcd_freq", 5000)

        # Margin parameter for NPAT simulations
        margin = stage_params.get("margin", 5.0)

        # Time in nanoseconds and timestep
        time_ns = stage_params.get("time_ns", 0.125)  # Default 125 ps = 0.125 ns
        timestep = stage_params.get(
            "timestep", 1.0
        )  # Default 1 fs (NAMD uses femtoseconds)

        # Calculate steps for display (but template uses the equation)
        calculated_steps = int(time_ns * 1e6 / timestep)

        # Minimization steps (only for first stage)
        minimize_steps = (
            stage_params.get("minimize_steps", 10000) if stage_index == 0 else 0
        )

        # Calculate firsttimestep based on previous stages
        first_timestep = self._calculate_first_timestep(
            stage_index, stage_params, all_stage_settings
        )

        # Handle cell basis vectors
        cell_basis_block = self._generate_cell_basis_block(stage_index)

        # Handle PME settings
        pme_block = self._generate_pme_block()

        # Handle restraints
        restraints_block = self._generate_restraints_block(
            stage_name, stage_params, stage_index
        )

        # Handle production steps for step13
        production_steps = stage_params.get("steps", 50000000)  # Default 50M steps

        # Generate input/output names for NAMD TCL variables
        output_name = self._generate_output_name(stage_name, stage_index)
        input_name = self._generate_input_name(stage_index, previous_stage_name)

        # Generate initial temperature directive (only for first stage)
        if stage_index == 0:  # First stage gets initial temperature assignment
            initial_temp_directive = f"temperature            $temp               # Initial temperature assignment for first stage"
        else:  # Subsequent stages don't set initial temperature (read from restart files)
            initial_temp_directive = "# No initial temperature assignment - reading velocities from restart file"

        # Perform replacements
        customized_content = template_content.replace("{TEMPERATURE}", str(temperature))
        customized_content = customized_content.replace("{DCD_FREQ}", str(dcd_freq))
        customized_content = customized_content.replace(
            "{OUTPUT_ENERGIES}", str(dcd_freq)
        )
        customized_content = customized_content.replace("{XST_FREQ}", str(dcd_freq))
        customized_content = customized_content.replace(
            "{OUTPUT_TIMING}", str(dcd_freq)
        )
        customized_content = customized_content.replace("{MARGIN}", str(margin))
        customized_content = customized_content.replace("{TIME_NS}", str(time_ns))
        customized_content = customized_content.replace("{TIMESTEP}", str(timestep))
        customized_content = customized_content.replace(
            "{CELL_BASIS_VECTORS}", cell_basis_block
        )
        customized_content = customized_content.replace("{PME_SETTINGS}", pme_block)
        customized_content = customized_content.replace(
            "{RESTRAINTS_BLOCK}", restraints_block
        )
        customized_content = customized_content.replace(
            "{PRODUCTION_STEPS}", str(production_steps)
        )
        customized_content = customized_content.replace(
            "{RUN_STEPS}", str(calculated_steps)
        )
        customized_content = customized_content.replace(
            "{MINIMIZE_STEPS}", str(minimize_steps)
        )
        customized_content = customized_content.replace(
            "{FIRST_TIMESTEP}", str(first_timestep)
        )
        customized_content = customized_content.replace(
            "{INITIAL_TEMPERATURE_DIRECTIVE}", initial_temp_directive
        )
        customized_content = customized_content.replace(
            "{WATER_MODEL_BLOCK}",
            namd_water_model_config_block(getattr(self, "water_model", "tip3p")),
        )

        # Replace system file paths (parmfile and ambercoor)
        # This allows using either relative paths (when files are copied) or absolute paths
        import re

        if "prmtop" in system_files:
            prmtop_path = system_files["prmtop"]
            customized_content = re.sub(
                r"parmfile\s+[\w/.]+\.prmtop",
                f"parmfile                {prmtop_path}",
                customized_content,
            )
            customized_content = re.sub(
                r"parmfile\s+[\w/.]+\.parm7",
                f"parmfile                {prmtop_path}",
                customized_content,
            )
            customized_content = re.sub(
                r"parmfile\s+[\w/.]+\.top",
                f"parmfile                {prmtop_path}",
                customized_content,
            )

        if "inpcrd" in system_files:
            inpcrd_path = system_files["inpcrd"]
            customized_content = re.sub(
                r"ambercoor\s+[\w/.]+\.inpcrd",
                f"ambercoor               {inpcrd_path}",
                customized_content,
            )
            customized_content = re.sub(
                r"ambercoor\s+[\w/.]+\.rst7",
                f"ambercoor               {inpcrd_path}",
                customized_content,
            )
            customized_content = re.sub(
                r"ambercoor\s+[\w/.]+\.rst",
                f"ambercoor               {inpcrd_path}",
                customized_content,
            )

        # Replace NAMD TCL variable names if they exist as placeholders
        customized_content = customized_content.replace("{OUTPUT_NAME}", output_name)
        customized_content = customized_content.replace("{INPUT_NAME}", input_name)

        # Replace hardcoded outputname and inputname in templates using regex

        # Replace set outputname lines (e.g., "set outputname eq6_equilibration;" -> "set outputname eq7_equilibration;")
        customized_content = re.sub(
            r"set\s+outputname\s+\w+;",
            f"set outputname          {output_name};",
            customized_content,
        )

        # Replace set inputname lines (e.g., "set inputname eq5_equilibration;" -> "set inputname eq6_equilibration;")
        # Only replace if we have an input name (not first stage)
        if input_name:
            customized_content = re.sub(
                r"set\s+inputname\s+\w+;",
                f"set inputname           {input_name};",
                customized_content,
            )
        else:
            # For first stage, remove or comment out the inputname line AND the restart file directives
            customized_content = re.sub(
                r"set\s+inputname\s+\w+;",
                "# set inputname           (not needed for first stage);",
                customized_content,
            )
            # Also comment out the restart file directives that depend on inputname
            customized_content = re.sub(
                r"binCoordinates\s+\$inputname\.coor;",
                "# binCoordinates          $inputname.coor;    # (not needed for first stage)",
                customized_content,
            )
            customized_content = re.sub(
                r"binVelocities\s+\$inputname\.vel;",
                "# binVelocities           $inputname.vel;     # (not needed for first stage)",
                customized_content,
            )
            customized_content = re.sub(
                r"extendedSystem\s+\$inputname\.xsc;",
                "# extendedSystem          $inputname.xsc;     # (not needed for first stage)",
                customized_content,
            )

        return customized_content

    def _calculate_first_timestep(
        self,
        stage_index: int,
        stage_params: Dict[str, Any],
        all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> int:
        """Calculate the first timestep for a stage based on previous stages."""
        if stage_index == 0:
            return 0  # First stage always starts from 0

        if all_stage_settings is None:
            # Fallback to old behavior if all_stage_settings not provided
            current_stage_steps = stage_params.get("steps", 125000)
            cumulative_steps = stage_index * current_stage_steps
            return cumulative_steps

        # Build list of actual stage keys from all_stage_settings in order
        # The keys might be in format "Equilibration 1", "Equilibration 2", ..., "Production"
        # or "equilibration_1", "equilibration_2", ..., "production"
        stage_keys = list(all_stage_settings.keys())

        cumulative_steps = 0
        for i in range(stage_index):
            if i < len(stage_keys):
                stage_key = stage_keys[i]
                stage_config = all_stage_settings[stage_key]

                # Get run steps based on time_ns and timestep
                time_ns = stage_config.get("time_ns", 0.125)
                timestep = stage_config.get("timestep", 1.0)
                run_steps = int(time_ns * 1e6 / timestep)

                # Add minimize steps for first stage only
                if i == 0:
                    minimize_steps = stage_config.get("minimize_steps", 10000)
                    cumulative_steps += minimize_steps

                cumulative_steps += run_steps
            else:
                # Default steps if stage not found
                cumulative_steps += 125000

        return cumulative_steps

    def _generate_output_name(self, stage_name: str, stage_index: int) -> str:
        """Generate output name for NAMD configuration."""
        # Use the new config name mapping function with stage_index
        config_name = self._get_config_name(stage_name, stage_index)
        # Production uses step7_production, equilibration stages use step{N}_equilibration
        if config_name == "step7_production":
            return config_name
        else:
            return f"{config_name}_equilibration"

    def _generate_input_name(
        self, stage_index: int, previous_stage_name: Optional[str] = None
    ) -> str:
        """Generate input name for NAMD restart files."""
        if stage_index == 0:
            # First stage doesn't need input name
            return ""

        if previous_stage_name:
            # Use previous stage output name - convert previous stage name to config name
            # Use stage_index - 1 for the previous stage
            prev_config_name = self._get_config_name(
                previous_stage_name, stage_index - 1
            )
            # Production uses step7_production, equilibration stages use step{N}_equilibration
            if prev_config_name == "step7_production":
                return prev_config_name
            else:
                return f"{prev_config_name}_equilibration"

        # Fallback: generate based on stage index using direct step naming
        # For stage_index 1 (2nd stage) -> need step1_equilibration as input
        # For stage_index 2 (3rd stage) -> need step2_equilibration as input, etc.
        if stage_index > 0:
            prev_step_num = stage_index  # This gives us the previous step number
            return f"step{prev_step_num}_equilibration"

        # Final fallback (should never reach here)
        return ""

    def _generate_cell_basis_block(self, stage_index: int) -> str:
        """Generate cell basis vectors block - only for first stage, others use .xsc files."""
        if stage_index > 0:
            # After first stage, box dimensions come from .xsc files
            return "# Cell dimensions read from .xsc file"

        # For first stage, get box dimensions with proper priority order
        try:
            # PRIORITY 1: Read from bilayer*_lipid.pdb files with CRYST1 records (most accurate)
            bilayer_pdb = self._find_bilayer_pdb_with_cryst1()
            if bilayer_pdb and bilayer_pdb.exists():
                a, b, c = self._read_box_dimensions(bilayer_pdb)
                self.logger.info(
                    f"Using box dimensions from bilayer*_lipid.pdb for cell basis: {bilayer_pdb.name}"
                )
            else:
                # PRIORITY 2: Try system.pdb as fallback
                system_pdb = self.working_dir / "system.pdb"
                if system_pdb.exists():
                    a, b, c = self._read_box_dimensions(system_pdb)
                    self.logger.info(
                        f"Using box dimensions from system.pdb for cell basis: {system_pdb.name}"
                    )
                else:
                    # PRIORITY 3: Default dimensions
                    a, b, c = 100.0, 100.0, 100.0
                    self.logger.warning(
                        "No PDB files found, using default box dimensions for cell basis"
                    )

            zcen = 0.0
        except Exception as e:
            self.logger.warning(
                f"Error reading box dimensions for cell basis: {e}, using defaults"
            )
            a, b, c = 100.0, 100.0, 100.0
            zcen = 0.0

        return f"""cellBasisVector1     {a:.3f}   0.0   0.0
cellBasisVector2     0.0   {b:.3f}   0.0
cellBasisVector3     0.0   0.0   {c:.3f}
cellOrigin           0.0   0.0   {zcen:.3f}"""

    def _generate_pme_block(self) -> str:
        """Generate PME settings block without hardcoded grid sizes."""
        # PME grid sizes should be calculated automatically by NAMD
        # based on the system size and PMEGridSpacing
        return f"""PME                     yes
PMEInterpOrder          6
PMEGridSpacing          1.0"""

    def _generate_restraints_block(
        self, stage_name: str, stage_params: Dict[str, Any], stage_index: int
    ) -> str:
        """Generate restraints block using GateWizard's restraint system."""
        constraints = stage_params.get("constraints", {})

        # Check if any constraints are defined
        has_restraints = any(float(v) > 0 for v in constraints.values())

        restraints_lines = []

        # Add position restraints if defined
        if has_restraints:
            # Map stage names to restraint file names using the correct naming scheme
            # Based on the display name to config name conversion
            config_name = self._get_config_name(stage_name, stage_index)

            if config_name == "step7_production":
                restraint_file = (
                    f"{config_name}_restraints.pdb"  # step7_production_restraints.pdb
                )
            else:
                # For equilibration stages: step1_equilibration_restraints.pdb, step2_equilibration_restraints.pdb, etc.
                restraint_file = f"{config_name}_equilibration_restraints.pdb"

            # Use constraintScaling = 1.0 and keep GUI force values
            restraints_lines.extend(
                [
                    "# Position restraints",
                    "constraints             on",
                    "consexp                 2",
                    f"consref                 restraints/{restraint_file}",
                    f"conskfile               restraints/{restraint_file}",
                    "conskcol                B",
                    "constraintScaling       1.0",
                ]
            )

        # If no restraints at all, add a comment
        if not restraints_lines:
            return "# No restraints defined for this stage"

        return "\n".join(restraints_lines)


class OpenMMEquilibrationManager:
    """Manager for OpenMM equilibration simulations using CHARMM-GUI templates.

    Mirrors the NAMDEquilibrationManager API. The ``constraints`` dict accepts
    the same keys and kcal/mol/Å² units as NAMD; non-zero values control which
    atom types appear in the OpenMM restraint index files (prot_pos.txt,
    lipid_pos.txt). Per-stage force constant magnitudes follow the CHARMM-GUI
    equilibration schedule embedded in the template files. Lipid dihedral
    restraints (fc_ldih) are always disabled; only positional restraints are used.
    """

    SCHEME_MAPPING: Dict[str, str] = {
        "NVT": "01_NVT",
        "NPT": "02_NPT",
        "NPAT": "03_NPAT",
        "NPgT": "04_NPgT",
    }

    TEMPLATE_MAPPING: Dict[str, str] = {
        "step1": "step6.1_equilibration.inp",
        "step2": "step6.2_equilibration.inp",
        "step3": "step6.3_equilibration.inp",
        "step4": "step6.4_equilibration.inp",
        "step5": "step6.5_equilibration.inp",
        "step6": "step6.6_equilibration.inp",
        "step7_production": "step7_production.inp",
    }

    STAGE_INDEX_TO_KEY: Dict[int, str] = {
        1: "step1",
        2: "step2",
        3: "step3",
        4: "step4",
        5: "step5",
        6: "step6",
        7: "step7_production",
    }

    def __init__(self, working_dir: Path):
        self.working_dir = Path(working_dir)
        self.templates_dir = (
            Path(__file__).parent.parent.parent / "equilibration" / "openmm"
        )
        self.scripts_dir = self.templates_dir / "scripts"
        self.logger = get_logger(self.__class__.__name__)

    def find_system_files(self) -> Optional[Dict[str, str]]:
        """
        Automatically detect AMBER system files in the working directory.

        Returns:
            Dict with keys ``prmtop``, ``inpcrd``, ``pdb``, ``bilayer_pdb``
            or None on failure.

        Example:
            >>> manager = OpenMMEquilibrationManager(Path("work_dir"))
            >>> files = manager.find_system_files()
            >>> if files:
            ...     result = manager.setup_openmm_equilibration(system_files=files)
        """
        system_files: Dict[str, Any] = {}

        prmtop_files = list(self.working_dir.glob("*.prmtop"))
        if not prmtop_files:
            self.logger.error("No .prmtop file found in working directory")
            return None
        system_files["prmtop"] = str(prmtop_files[0])
        self.logger.info(f"Found topology: {prmtop_files[0].name}")

        inpcrd_files = list(self.working_dir.glob("*.inpcrd"))
        if not inpcrd_files:
            inpcrd_files = list(self.working_dir.glob("*.crd"))
        if not inpcrd_files:
            inpcrd_files = list(self.working_dir.glob("*.rst"))
        if not inpcrd_files:
            self.logger.error("No .inpcrd/.crd/.rst file found in working directory")
            return None
        system_files["inpcrd"] = str(inpcrd_files[0])
        self.logger.info(f"Found coordinates: {inpcrd_files[0].name}")

        system_pdb = self.working_dir / "system.pdb"
        if system_pdb.exists():
            system_files["pdb"] = str(system_pdb)
        else:
            pdb_files = [
                p
                for p in self.working_dir.glob("*.pdb")
                if "bilayer" not in p.name.lower()
            ]
            if pdb_files:
                system_files["pdb"] = str(pdb_files[0])
            else:
                self.logger.warning(
                    "No .pdb file found; restraint generation will be skipped"
                )
                system_files["pdb"] = None

        if system_files.get("pdb"):
            self.logger.info(f"Found PDB: {Path(system_files['pdb']).name}")

        # Bilayer PDB with CRYST1 record — needed to supply box dimensions when
        # prmtop has IFBOX=0 (common in membrane-system preparations).
        bilayer_pdb: Optional[Path] = None
        for pattern in ("bilayer*_lipid.pdb", "bilayer_*.pdb", "*_bilayer.pdb"):
            candidates = list(self.working_dir.glob(pattern))
            if candidates:
                bilayer_pdb = candidates[0]
                break
        if bilayer_pdb:
            system_files["bilayer_pdb"] = str(bilayer_pdb)
            self.logger.info(f"Found bilayer PDB: {bilayer_pdb.name}")
        else:
            system_files["bilayer_pdb"] = None
            self.logger.warning(
                "No bilayer PDB (bilayer_*_lipid.pdb) found; box dimensions may "
                "not be set correctly if prmtop has IFBOX=0"
            )

        return system_files

    @staticmethod
    def get_default_selections(pdb_path: str) -> Dict[str, str]:
        """
        Auto-detect MDAnalysis selection strings for the system in *pdb_path*.

        Delegates to :meth:`NAMDEquilibrationManager.get_default_selections`,
        which inspects the PDB for standard protein/lipid residues and any
        non-standard residues that should be treated as ligands.

        Returns:
            Dict mapping category keys (``"protein_backbone"``,
            ``"protein_sidechain"``, ``"lipid_head"``, ``"lipid_tail"``,
            ``"ligand_<RESNAME>"``, …) to MDAnalysis selection strings.

        Example:
            >>> sels = OpenMMEquilibrationManager.get_default_selections("system.pdb")
            >>> print(sels["protein_backbone"])
            protein and backbone
        """
        return NAMDEquilibrationManager.get_default_selections(pdb_path)

    def setup_openmm_equilibration(
        self,
        system_files: Optional[Dict[str, str]] = None,
        stage_params_list: Optional[List[Dict[str, Any]]] = None,
        output_name: str = "equilibration",
        scheme_type: Optional[str] = None,
        selections: Optional[Dict[str, str]] = None,
        add_com_restraint: bool = False,
        com_restraint_k: float = 10.0,
        add_rotation_restraint: bool = False,
        rotation_restraint_k: float = 2000.0,
        com_selection: str = "name CA",
    ) -> Dict[str, Any]:
        """
        Complete OpenMM equilibration setup.

        Generates .inp configuration files, restraint index files (prot_pos.txt,
        lipid_pos.txt, dihe.txt), copies the CHARMM-GUI Python runner scripts,
        and produces a bash run script.

        Args:
            system_files: Dict with ``prmtop``, ``inpcrd``, ``pdb``. Auto-detected
                from working_dir if None.
            stage_params_list: List of stage parameter dicts. Supported keys:

                - ``name`` (str): human-readable label (optional)
                - ``ensemble`` (str): NVT | NPT | NPAT | NPgT
                - ``time_ns`` (float): simulation time in nanoseconds
                - ``timestep`` (float): integration timestep in femtoseconds (default 2.0)
                - ``temperature`` (float): temperature in Kelvin (default 310.15)
                - ``dcd_freq`` (int): trajectory write frequency in steps (default 5000)
                - ``minimize_steps`` (int): minimization steps for first stage (default 5000)
                - ``constraints`` (dict): kcal/mol/Å² force constants; keys:
                  ``protein_backbone``, ``protein_sidechain``, ``lipid_head``,
                  ``lipid_tail``, ``water``, ``ions``, ``other``, or any
                  custom key such as ``ligand_ABC``.  Non-zero values generate
                  the corresponding restraint index files and set ``fc_bb``,
                  ``fc_sc``, ``fc_lpos`` in the .inp file (converted to
                  kJ/mol/nm²).  Custom keys write ``restraints/custom_pos.txt``.

            output_name: Subdirectory name under working_dir.
            scheme_type: Override ensemble for all stages. Auto-detected from the
                ``ensemble`` field of the first stage if None.
            selections: Optional ``{key: mda_selection_string}`` dict that
                overrides the auto-detected MDAnalysis selections used when
                generating restraint index files.  Keys should match the
                ``constraints`` dict keys (e.g. ``"protein_backbone"``,
                ``"ligand_ABC"``).  Requires MDAnalysis.

        Returns:
            Dict with keys ``openmm_dir``, ``config_files``, ``run_script``,
            and ``system_files``.

        Example:
            >>> from pathlib import Path
            >>> from gatewizard.tools.equilibration import OpenMMEquilibrationManager
            >>> stages = [
            ...     {"name": "Eq1", "ensemble": "NPT", "time_ns": 0.125,
            ...      "temperature": 310.15, "timestep": 1.0,
            ...      "constraints": {"protein_backbone": 10.0, "lipid_head": 2.5}},
            ... ]
            >>> manager = OpenMMEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_openmm_equilibration(stage_params_list=stages)
            >>> print(result["openmm_dir"])
        """
        self.logger.info("=== Setting up OpenMM equilibration ===")

        if system_files is None:
            system_files = self.find_system_files()
            if system_files is None:
                raise ValueError(
                    "Could not auto-detect system files in working directory"
                )

        if not stage_params_list:
            _scheme_for_defaults = scheme_type or "NPT"
            stage_params_list = self.get_default_stage_params(_scheme_for_defaults)
            self.logger.info(
                f"No stages provided — using default {_scheme_for_defaults} protocol "
                f"({len(stage_params_list)} stages)"
            )

        # Normalise EquilibrationStage objects to plain dicts
        stage_params_list = [
            s.to_dict() if isinstance(s, EquilibrationStage) else s
            for s in stage_params_list
        ]

        if scheme_type is None:
            scheme_type = stage_params_list[0].get("ensemble", "NPT")
            self.logger.info(f"Auto-detected scheme_type: {scheme_type}")

        if scheme_type not in self.SCHEME_MAPPING:
            raise ValueError(
                f"Unknown scheme_type '{scheme_type}'. "
                f"Must be one of {list(self.SCHEME_MAPPING.keys())}"
            )

        openmm_dir = self.working_dir / output_name
        openmm_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Output directory: {openmm_dir}")

        # Copy system files
        self.logger.info("Copying system files...")
        for key in ("prmtop", "inpcrd", "pdb", "bilayer_pdb"):
            src = system_files.get(key)
            if src and Path(src).exists():
                dest = openmm_dir / Path(src).name
                shutil.copy2(src, dest)
                self.logger.info(f"  Copied {Path(src).name}")

        # Copy OpenMM Python runner scripts
        self.logger.info("Copying OpenMM Python scripts...")
        if self.scripts_dir.exists():
            for script in sorted(self.scripts_dir.glob("*.py")):
                shutil.copy2(script, openmm_dir / script.name)
                self.logger.info(f"  Copied {script.name}")
        else:
            self.logger.warning(f"Scripts directory not found: {self.scripts_dir}")

        # Generate restraint index files before configs (configs need has_custom flag)
        system_pdb_path = openmm_dir / Path(system_files.get("pdb", "system.pdb")).name
        restraint_files: Dict[str, Any] = {}
        has_custom_restraints = False
        if system_pdb_path.exists():
            restraint_files = self.generate_openmm_restraint_files(
                system_pdb=system_pdb_path,
                stage_params_list=stage_params_list,
                output_dir=openmm_dir,
                selections=selections,
            )
            has_custom_restraints = bool(restraint_files.get("custom_pos_per_stage"))
        else:
            self.logger.warning(
                f"System PDB not found at {system_pdb_path}; skipping restraint index generation."
            )

        config_files: List[Path] = []
        stage_config_names: List[str] = []
        custom_pos_per_stage: Dict[int, Optional[Path]] = restraint_files.get(
            "custom_pos_per_stage", {}
        )

        for stage_index, stage_params in enumerate(stage_params_list, start=1):
            stage_name = stage_params.get("name", f"Stage {stage_index}")
            self.logger.info(f"Processing stage {stage_index}: {stage_name}")

            config_content = self.generate_openmm_config(
                stage_name=stage_name,
                stage_params=stage_params,
                stage_index=stage_index,
                scheme_type=scheme_type,
                custom_pos_file=custom_pos_per_stage.get(stage_index),
            )

            config_filename = self._get_config_filename(stage_index)
            config_path = openmm_dir / config_filename
            config_path.write_text(config_content)
            config_files.append(config_path)
            stage_config_names.append(config_filename.replace(".inp", ""))
            self.logger.info(f"  Written: {config_filename}")

        prmtop_name = Path(system_files.get("prmtop", "system.prmtop")).name
        inpcrd_name = Path(system_files.get("inpcrd", "system.inpcrd")).name
        bilayer_pdb_src = system_files.get("bilayer_pdb")
        bilayer_pdb_name = Path(bilayer_pdb_src).name if bilayer_pdb_src else None
        compute = resolve_compute_resources_from_stages(stage_params_list)
        run_script_path = self.generate_run_script(
            stage_config_names=stage_config_names,
            openmm_dir=openmm_dir,
            prmtop_name=prmtop_name,
            inpcrd_name=inpcrd_name,
            bilayer_pdb_name=bilayer_pdb_name,
            cpu_cores=compute["cpu_cores"],
            use_gpu=compute["use_gpu"],
            gpu_id=compute["gpu_id"],
            num_gpus=compute["num_gpus"] or 1,
        )
        self.logger.info(f"Run script: {run_script_path.name}")
        self.logger.info("=== OpenMM equilibration setup complete ===")

        # --- COM restraint (writes com_restraint_params.json for omm_restraints.py) ---
        com_restraint_path: Optional[Path] = None
        if add_com_restraint:
            pdb_src = system_files.get("pdb") if system_files else None
            if pdb_src:
                pdb_for_com = openmm_dir / Path(pdb_src).name
                if not pdb_for_com.exists():
                    pdb_for_com = Path(pdb_src)
                com_restraint_path = self._write_openmm_com_params(
                    pdb_path=pdb_for_com,
                    output_dir=openmm_dir,
                    com_restraint_k=com_restraint_k,
                    add_rotation_restraint=add_rotation_restraint,
                    rotation_restraint_k=rotation_restraint_k,
                    selection=com_selection,
                )
            else:
                self.logger.warning(
                    "No PDB in system_files; COM restraint params not written."
                )

        return {
            "openmm_dir": openmm_dir,
            "config_files": config_files,
            "run_script": run_script_path,
            "system_files": system_files,
            "restraint_files": restraint_files,
            "com_restraint_params": com_restraint_path,
        }

    def generate_openmm_config(
        self,
        stage_name: str,
        stage_params: Dict[str, Any],
        stage_index: int,
        scheme_type: str,
        custom_pos_file: Optional[Path] = None,
        has_custom_restraints: bool = False,
    ) -> str:
        """
        Generate an OpenMM .inp configuration file for a single equilibration stage.

        Loads the CHARMM-GUI template for the given ensemble and stage, then
        substitutes runtime parameters including force constants derived from
        the ``constraints`` dict (converted from kcal/mol/Å² to kJ/mol/nm²).

        Args:
            stage_name: Label used for logging only.
            stage_params: Stage parameter dict (see setup_openmm_equilibration).
            stage_index: 1-based position of this stage in the protocol.
            scheme_type: Ensemble type (NVT, NPT, NPAT, or NPgT).
            custom_pos_file: When set, a ``custom_pos_file =`` line is injected
                into the .inp so ``omm_restraints.py`` loads the per-stage file.
                When None, no custom positional restraints are applied.
            has_custom_restraints: Deprecated — ignored; kept for back-compat.

        Returns:
            String content of the generated .inp file.
        """
        template_key = self.STAGE_INDEX_TO_KEY.get(stage_index, "step7_production")
        scheme_folder = self.SCHEME_MAPPING[scheme_type]
        template_filename = self.TEMPLATE_MAPPING[template_key]
        template_path = self.templates_dir / scheme_folder / template_filename

        if not template_path.exists():
            raise FileNotFoundError(
                f"OpenMM template not found: {template_path}. "
                f"Expected in equilibration/openmm/{scheme_folder}/"
            )

        content = template_path.read_text()

        temperature = float(stage_params.get("temperature", 310.15))
        timestep_fs = float(stage_params.get("timestep", 2.0))
        dt_ps = timestep_fs / 1000.0
        time_ns = float(stage_params.get("time_ns", 0.5))
        nstep = max(1, int(round(time_ns * 1000.0 / dt_ps)))

        is_production = template_key == "step7_production"
        dcd_freq_default = 50000 if is_production else 5000
        dcd_freq = int(stage_params.get("dcd_freq", dcd_freq_default))

        # Compute OpenMM force constants from constraints dict
        # Conversion: 1 kcal/mol/Å² = 418.4 kJ/mol/nm²
        _KCAL_TO_KJ = 418.4
        _STD_KEYS = frozenset(
            {
                "protein_backbone",
                "protein_sidechain",
                "lipid_head",
                "lipid_tail",
            }
        )
        constraints = stage_params.get("constraints", {})
        fc_bb_kj = float(constraints.get("protein_backbone", 0.0)) * _KCAL_TO_KJ
        fc_sc_kj = float(constraints.get("protein_sidechain", 0.0)) * _KCAL_TO_KJ
        fc_lpos_kj = (
            max(
                float(constraints.get("lipid_head", 0.0)),
                float(constraints.get("lipid_tail", 0.0)),
            )
            * _KCAL_TO_KJ
        )
        # Check whether this stage uses any custom (non-standard) restraint
        stage_has_custom = custom_pos_file is not None
        rest = (
            "yes"
            if (fc_bb_kj > 0 or fc_sc_kj > 0 or fc_lpos_kj > 0 or stage_has_custom)
            else "no"
        )

        is_first_stage = stage_index == 1
        minimize_steps = (
            int(stage_params.get("minimize_steps", 5000)) if is_first_stage else 0
        )

        content = content.replace("{TEMPERATURE}", f"{temperature:.2f}")
        content = content.replace("{NSTEP}", str(nstep))
        content = content.replace("{DT}", f"{dt_ps:.3f}")
        content = content.replace("{NSTDCD}", str(dcd_freq))
        # Some production templates ship with a literal "rest = no" instead of
        # "{REST}". Make restraint toggling robust for both formats.
        if "{REST}" in content:
            content = content.replace("{REST}", rest)
        else:
            content = re.sub(
                r"(?m)^(rest\s*=\s*)(yes|no)(\s*#.*)?$",
                lambda m: f"{m.group(1)}{rest}{m.group(3) or ''}",
                content,
            )

        # Older production templates may omit FC lines entirely.
        # Ensure OpenMM input always carries user force constants.
        if not re.search(r"(?m)^\s*fc_bb\s*=", content):
            content += (
                "\n"
                "fc_bb       = 0.0                                   # Positional restraint force constant for protein backbone (kJ/mol/nm^2)\n"
            )
        if not re.search(r"(?m)^\s*fc_sc\s*=", content):
            content += "fc_sc       = 0.0                                   # Positional restraint force constant for protein side-chain (kJ/mol/nm^2)\n"
        if not re.search(r"(?m)^\s*fc_lpos\s*=", content):
            content += "fc_lpos     = 0.0                                   # Positional restraint force constant for lipids (kJ/mol/nm^2)\n"
        # Override CHARMM-GUI template force constants with user-specified values
        content = re.sub(
            r"fc_bb\s*=\s*[\d.]+", f"fc_bb       = {fc_bb_kj:.4f}", content
        )
        content = re.sub(
            r"fc_sc\s*=\s*[\d.]+", f"fc_sc       = {fc_sc_kj:.4f}", content
        )
        content = re.sub(
            r"fc_lpos\s*=\s*[\d.]+", f"fc_lpos     = {fc_lpos_kj:.4f}", content
        )
        if custom_pos_file is not None:
            # Inject the per-stage custom restraint file path as a parseable .inp parameter
            # so omm_restraints.py loads the correct file for this stage.
            content += (
                "\n"
                f"custom_pos_file = restraints/{custom_pos_file.name}"
                "                   # Per-stage custom positional restraints (GateWizard)\n"
            )
        if is_first_stage:
            content = content.replace("{MINI_NSTEP}", str(minimize_steps))
            content = content.replace("{GEN_VEL}", "yes")

        self.logger.debug(
            f"Stage {stage_index} ({stage_name}): template={template_filename}, "
            f"T={temperature:.2f}K, nstep={nstep}, dt={dt_ps:.3f}ps, rest={rest}, "
            f"fc_bb={fc_bb_kj:.1f}, fc_sc={fc_sc_kj:.1f}, fc_lpos={fc_lpos_kj:.1f} kJ/mol/nm²"
        )
        return content

    def _write_openmm_com_params(
        self,
        pdb_path: Path,
        output_dir: Path,
        com_restraint_k: float = 10.0,
        add_rotation_restraint: bool = False,
        rotation_restraint_k: float = 2000.0,
        selection: str = "name CA",
    ) -> Optional[Path]:
        """Write ``com_restraint_params.json`` for the OpenMM runtime script.

        The JSON contains the 0-based Cα atom indices and their initial centroid
        so ``omm_restraints.py`` can apply a ``CustomCentroidBondForce`` to
        restrain the protein centre of mass. When rotation restraint is enabled,
        three anchor Cα atoms and their reference coordinates are also stored so
        ``omm_restraints.py`` can apply a rotational restraint.

        Args:
            pdb_path: System PDB file.
            output_dir: Directory to write ``com_restraint_params.json``.
            com_restraint_k: Force constant in kcal/mol/Å².
            add_rotation_restraint: Whether to enable rotational restraint.
            rotation_restraint_k: Rotation force constant in kcal/mol/Å².
            selection: MDAnalysis selection string for atoms used to define
                COM translation and optional rotation anchors.

        Returns:
            Path to the written JSON, or ``None`` on failure.
        """
        try:
            import MDAnalysis as mda  # type: ignore

            u = mda.Universe(str(pdb_path))
            ag = u.select_atoms(selection)
            if len(ag) == 0:
                self.logger.warning(
                    f"No atoms matched '{selection}' for OpenMM COM restraint."
                )
                return None

            com = ag.center_of_geometry()
            rotation_anchor_indices: list[int] = []
            rotation_ref_positions_angstrom: list[list[float]] = []

            if add_rotation_restraint and len(ag) >= 3:
                import numpy as np

                coords = np.asarray(ag.positions, dtype=float)
                centered = coords - coords.mean(axis=0)
                _, _, vh = np.linalg.svd(centered, full_matrices=False)
                pc1 = vh[0]
                pc2 = vh[1] if vh.shape[0] > 1 else np.array([0.0, 1.0, 0.0])

                proj1 = centered @ pc1
                idx_a = int(np.argmin(proj1))
                idx_b = int(np.argmax(proj1))

                proj2 = np.abs(centered @ pc2)
                proj2[[idx_a, idx_b]] = -1.0
                idx_c = int(np.argmax(proj2))

                # Fallback when all points collapse onto the first principal axis.
                if idx_c in (idx_a, idx_b):
                    for i in range(len(ag)):
                        if i not in (idx_a, idx_b):
                            idx_c = i
                            break

                anchor_atoms = [ag[idx_a], ag[idx_b], ag[idx_c]]
                rotation_anchor_indices = [int(a.index) for a in anchor_atoms]
                rotation_ref_positions_angstrom = [
                    [float(a.position[0]), float(a.position[1]), float(a.position[2])]
                    for a in anchor_atoms
                ]
            elif add_rotation_restraint:
                self.logger.warning(
                    "Rotation restraint requested for OpenMM but fewer than 3 Cα atoms were found."
                )

            params = {
                "ca_indices": [int(a.index) for a in ag],
                "centroid_angstrom": [float(com[0]), float(com[1]), float(com[2])],
                "force_constant_kcal_mol_A2": com_restraint_k,
                "add_rotation_restraint": add_rotation_restraint,
                "rotation_force_constant_kcal_mol_A2": rotation_restraint_k,
                "rotation_anchor_indices": rotation_anchor_indices,
                "rotation_ref_positions_angstrom": rotation_ref_positions_angstrom,
            }
            out = output_dir / "com_restraint_params.json"
            out.write_text(json.dumps(params, indent=2))
            self.logger.info(
                f"  COM restraint params ({len(ag)} Cα, "
                f"centroid [{com[0]:.2f},{com[1]:.2f},{com[2]:.2f}] Å): {out.name}"
            )
            self.logger.info(
                "  To activate OpenMM COM restraint, pass --com-restraint flag "
                "to omm_restraints.py (see equilibration/openmm/scripts/omm_restraints.py)."
            )
            return out

        except ImportError:
            self.logger.warning(
                "MDAnalysis not available; OpenMM COM restraint params not written."
            )
            return None
        except Exception as exc:
            self.logger.error(f"OpenMM COM restraint param writing failed: {exc}")
            return None

    def generate_run_script(
        self,
        stage_config_names: List[str],
        openmm_dir: Path,
        prmtop_name: str,
        inpcrd_name: str,
        bilayer_pdb_name: Optional[str] = None,
        cpu_cores: Optional[int] = None,
        use_gpu: Optional[bool] = None,
        gpu_id: int = 0,
        num_gpus: int = 1,
        platform: Optional[str] = None,
    ) -> Path:
        """
        Generate a bash script that runs all equilibration stages sequentially.

        Each stage restarts from the previous stage's .rst file. The script
        exits with an error message if any stage fails.

        When ``bilayer_pdb_name`` is provided it is passed as ``-b`` to every
        ``openmm_run.py`` invocation.  This allows openmm_run.py to recover
        the periodic box dimensions from the CRYST1 record when the prmtop has
        IFBOX=0 (common in membrane-system preparations).

        Args:
            stage_config_names: List of config base names (without .inp extension).
            openmm_dir: Directory where the script is written.
            prmtop_name: Filename of the AMBER topology.
            inpcrd_name: Filename of the AMBER coordinates.
            bilayer_pdb_name: Filename of bilayer PDB with CRYST1 box record (optional).
            cpu_cores: CPU thread count (``Threads`` property / ``--threads``).
            use_gpu: Prefer GPU platforms when True; force CPU when False.
            gpu_id: First GPU device index for ``--device``.
            num_gpus: Number of consecutive GPU devices starting at ``gpu_id``.
            platform: Explicit OpenMM platform name (CUDA / OpenCL / CPU / Metal).

        Returns:
            Path to the generated ``run_equilibration.sh`` script.
        """
        default_platform = (platform or "").strip()
        if not default_platform and use_gpu is False:
            # Explicit CPU request from the UI (compute target = CPU).
            default_platform = "CPU"
        # use_gpu True with no platform → leave empty so openmm_run auto-picks
        # CUDA > OpenCL > CPU (same as previous behaviour).

        device_index = ""
        wants_gpu = use_gpu is True or (
            default_platform
            and default_platform.upper() not in {"", "AUTO", "CPU", "REFERENCE"}
        )
        if wants_gpu:
            ngpu = max(1, int(num_gpus or 1))
            gid = int(gpu_id or 0)
            device_index = ",".join(str(gid + i) for i in range(ngpu))

        threads = int(cpu_cores) if cpu_cores and int(cpu_cores) > 0 else None

        lines = [
            "#!/bin/bash",
            "## OpenMM Equilibration Run Script",
            "## Generated by GateWizard - run from the directory containing this file",
            "",
            "# --- Platform selection ---",
            "# Auto-detects CUDA > OpenCL > CPU by default when PLATFORM is empty.",
            "# Override with: PLATFORM=CPU bash run_equilibration.sh",
            "# Override with: PLATFORM=OpenCL bash run_equilibration.sh",
            f'PLATFORM="${{PLATFORM:-{default_platform}}}"',
            "",
        ]
        if device_index:
            lines += [
                "# GPU device index (OpenMM DeviceIndex / OpenCLDeviceIndex)",
                f'DEVICE_INDEX="${{DEVICE_INDEX:-{device_index}}}"',
                "",
            ]
        if threads is not None:
            lines += [
                "# CPU thread count (used when PLATFORM=CPU)",
                f'THREADS="${{THREADS:-{threads}}}"',
                "",
            ]
        lines += [
            "# Override Python interpreter with: PYTHON=python3 bash run_equilibration.sh",
            'PYTHON="${PYTHON:-python}"',
            f'PRMTOP="{prmtop_name}"',
            f'INPCRD="{inpcrd_name}"',
        ]
        if bilayer_pdb_name:
            lines.append(
                f'BILAYER_PDB="{bilayer_pdb_name}"  # CRYST1 box source when prmtop IFBOX=0'
            )
        resource_bits = []
        if threads is not None:
            resource_bits.append(f"{threads} CPU threads")
        if device_index:
            resource_bits.append(f"GPU device(s) {device_index}")
        resource_echo = ", ".join(resource_bits) if resource_bits else "auto"
        lines += [
            "",
            'echo "Starting OpenMM equilibration protocol..."',
            f'echo "Resources: {resource_echo}"',
            "",
            OPENMM_RESUME_SHELL,
            "",
        ]

        for i, config_name in enumerate(stage_config_names):
            stage_num = i + 1
            inp_file = f"{config_name}.inp"
            rst_out = f"{config_name}.rst"
            dcd_out = f"{config_name}.dcd"
            log_out = f"{config_name}.log"

            lines.append(f"# Stage {stage_num}: {config_name}")
            lines.append(f'if [ "$RESUME" = "1" ] && _gw_openmm_stage_done "{config_name}"; then')
            lines.append(f'  echo "RESUME: skipping stage {stage_num} ({config_name})"')
            lines.append("else")
            cmd = f"$PYTHON openmm_run.py -i {inp_file} -ff amber -p $PRMTOP -c $INPCRD"
            if bilayer_pdb_name:
                cmd += " -b $BILAYER_PDB"
            if i > 0:
                prev_rst = f"{stage_config_names[i - 1]}.rst"
                cmd += f" -irst {prev_rst}"
            cmd += (
                f" -orst {rst_out} -odcd {dcd_out} ${{PLATFORM:+--platform $PLATFORM}}"
            )
            if device_index:
                cmd += " ${DEVICE_INDEX:+--device $DEVICE_INDEX}"
            if threads is not None:
                cmd += " ${THREADS:+--threads $THREADS}"
            lines.append(f"  ({cmd}) 2>&1 | tee {log_out}")
            lines.append(
                f'  if [ ${{PIPESTATUS[0]}} -ne 0 ]; then echo "Stage {stage_num} ({config_name}) failed"; exit 1; fi'
            )
            lines.append("fi")
            lines.append("")

        lines.append('echo "Equilibration complete."')
        lines.append("")

        script_path = openmm_dir / "run_equilibration.sh"
        script_path.write_text("\n".join(lines))
        script_path.chmod(0o755)
        return script_path

    def generate_openmm_restraint_files(
        self,
        system_pdb: Path,
        stage_params_list: List[Dict[str, Any]],
        output_dir: Path,
        selections: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Optional[Path]]:
        """
        Generate OpenMM restraint index files for all active equilibration stages.

        Creates a ``restraints/`` subdirectory under *output_dir* and writes:

        - ``prot_pos.txt`` — 0-based atom indices labelled ``BB`` (backbone) or
          ``SC`` (sidechain); consumed by ``omm_restraints.py`` with ``fc_bb``
          and ``fc_sc``.
        - ``lipid_pos.txt`` — 0-based atom indices for lipid atoms; consumed with
          ``fc_lpos``.
        - ``custom_pos.txt`` — 0-based atom indices with per-atom force constants
          (kJ/mol/nm²) for ligands and any non-standard constraint categories;
          consumed by the GateWizard extension in ``omm_restraints.py``.

        Only files that are needed (i.e. the corresponding constraint force is > 0
        in at least one stage) are created.  MDAnalysis is used when available
        for accurate atom selection; a name-based fallback handles standard
        categories when it is absent.  Custom categories always require MDAnalysis.

        Args:
            system_pdb: Path to the system PDB file.  Atom order must match the
                topology file (prmtop/psf) used by OpenMM.
            stage_params_list: List of stage parameter dicts; the ``constraints``
                sub-dict determines which categories need index files.
            output_dir: Top-level equilibration output directory.  A
                ``restraints/`` sub-directory is created here.
            selections: Optional ``{key: mda_selection_string}`` dict that
                overrides auto-detected MDAnalysis selections.  Keys must match
                the ``constraints`` dict keys (e.g. ``\"protein_backbone\"``,
                ``\"ligand_ABC\"``).

        Returns:
            ``{"prot_pos": Path|None, "lipid_pos": Path|None,
            "custom_pos": Path|None}``
        """
        _KCAL_TO_KJ = 418.4  # 1 kcal/mol/Å² = 418.4 kJ/mol/nm²
        _STD_KEYS = frozenset(
            {
                "protein_backbone",
                "protein_sidechain",
                "lipid_head",
                "lipid_tail",
            }
        )

        # Collect the maximum force per constraint key across all stages
        max_forces: Dict[str, float] = {}
        for sp in stage_params_list:
            for key, force in sp.get("constraints", {}).items():
                if float(force) > max_forces.get(key, 0.0):
                    max_forces[key] = float(force)

        needs_prot = (
            max_forces.get("protein_backbone", 0.0) > 0
            or max_forces.get("protein_sidechain", 0.0) > 0
        )
        needs_lipid = (
            max_forces.get("lipid_head", 0.0) > 0
            or max_forces.get("lipid_tail", 0.0) > 0
        )
        # Set of custom keys that have > 0 force in at least one stage
        all_custom_keys: set = {
            k for k, v in max_forces.items() if k not in _STD_KEYS and v > 0
        }

        result: Dict[str, Any] = {
            "prot_pos": None,
            "lipid_pos": None,
            "custom_pos_per_stage": {},  # {stage_index: Path | None}
        }

        if not (needs_prot or needs_lipid or all_custom_keys):
            self.logger.info(
                "No active restraints found; skipping index file generation."
            )
            return result

        restraints_dir = output_dir / "restraints"
        restraints_dir.mkdir(exist_ok=True)
        self.logger.info(f"Generating OpenMM restraint index files → {restraints_dir}")

        try:
            import MDAnalysis as mda
            import warnings

            if selections is None:
                selections = NAMDEquilibrationManager.get_default_selections(
                    str(system_pdb)
                )

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                u = mda.Universe(str(system_pdb))

            # --- prot_pos.txt ---
            if needs_prot:
                bb_sel = selections.get(
                    "protein_backbone",
                    NAMDEquilibrationManager.DEFAULT_SELECTIONS["protein_backbone"],
                )
                sc_sel = selections.get(
                    "protein_sidechain",
                    NAMDEquilibrationManager.DEFAULT_SELECTIONS["protein_sidechain"],
                )
                bb_atoms = u.select_atoms(bb_sel)
                sc_atoms = u.select_atoms(sc_sel)
                if len(bb_atoms) + len(sc_atoms) > 0:
                    prot_path = restraints_dir / "prot_pos.txt"
                    with open(prot_path, "w") as fh:
                        for atom in bb_atoms:
                            fh.write(f"{atom.index}  BB\n")
                        for atom in sc_atoms:
                            fh.write(f"{atom.index}  SC\n")
                    result["prot_pos"] = prot_path
                    self.logger.info(
                        f"  prot_pos.txt: {len(bb_atoms)} BB + {len(sc_atoms)} SC atoms"
                    )
                else:
                    self.logger.warning(
                        "  prot_pos.txt skipped: no protein atoms matched"
                    )

            # --- lipid_pos.txt ---
            if needs_lipid:
                lh_sel = selections.get(
                    "lipid_head",
                    NAMDEquilibrationManager.DEFAULT_SELECTIONS["lipid_head"],
                )
                head_atoms = u.select_atoms(lh_sel)
                if max_forces.get("lipid_tail", 0.0) > 0:
                    lt_sel = selections.get(
                        "lipid_tail",
                        NAMDEquilibrationManager.DEFAULT_SELECTIONS["lipid_tail"],
                    )
                    all_lipid = head_atoms | u.select_atoms(lt_sel)
                else:
                    all_lipid = head_atoms
                if len(all_lipid) > 0:
                    lipid_path = restraints_dir / "lipid_pos.txt"
                    with open(lipid_path, "w") as fh:
                        for atom in all_lipid:
                            fh.write(f"{atom.index}\n")
                    result["lipid_pos"] = lipid_path
                    self.logger.info(f"  lipid_pos.txt: {len(all_lipid)} atoms")
                else:
                    self.logger.warning(
                        "  lipid_pos.txt skipped: no lipid atoms matched"
                    )

            # --- custom_pos_stage{N}.txt (per-stage files) ---
            # Build atom groups once (expensive MDAnalysis query), then write
            # stage-specific files with the correct per-stage force constant.
            if all_custom_keys:
                # Resolve MDAnalysis selection strings for every custom key
                custom_atom_groups: Dict[str, Any] = {}
                for key in sorted(all_custom_keys):
                    sel_str = selections.get(key)
                    if sel_str is not None:
                        sel_alias = sel_str.strip().lower().replace(" ", "_")
                        if sel_alias in NAMDEquilibrationManager.DEFAULT_SELECTIONS:
                            sel_str = NAMDEquilibrationManager.DEFAULT_SELECTIONS[
                                sel_alias
                            ]
                        elif sel_alias == "ion":
                            sel_str = NAMDEquilibrationManager.DEFAULT_SELECTIONS[
                                "ions"
                            ]
                    elif key in NAMDEquilibrationManager.DEFAULT_SELECTIONS:
                        sel_str = NAMDEquilibrationManager.DEFAULT_SELECTIONS[key]
                    elif key == "ion":
                        sel_str = NAMDEquilibrationManager.DEFAULT_SELECTIONS["ions"]
                    if sel_str is None:
                        self.logger.warning(
                            f"  custom_pos: no selection for '{key}'. "
                            f"Add selections['{key}'] = 'resname ...' to include it."
                        )
                        continue
                    try:
                        ag = u.select_atoms(sel_str)
                    except Exception as exc:
                        self.logger.warning(
                            f"  custom_pos: selection for '{key}' failed: {exc}"
                        )
                        continue
                    if len(ag) == 0:
                        self.logger.warning(
                            f"  custom_pos: '{key}' matched 0 atoms (sel='{sel_str}')"
                        )
                        continue
                    custom_atom_groups[key] = ag

                # Write per-stage files with the correct force constant for each stage
                for stage_i, sp in enumerate(stage_params_list, start=1):
                    stage_constraints = sp.get("constraints", {})
                    stage_custom: Dict[str, float] = {
                        k: float(v)
                        for k, v in stage_constraints.items()
                        if k not in _STD_KEYS
                        and float(v) > 0
                        and k in custom_atom_groups
                    }
                    if not stage_custom:
                        result["custom_pos_per_stage"][stage_i] = None
                        continue
                    custom_lines: List[str] = [
                        "# GateWizard custom positional restraints\n",
                        "# Format: atom_index(0-based)  force_kJ_mol_nm2\n",
                    ]
                    for key in sorted(stage_custom):
                        force_kcal = stage_custom[key]
                        force_kj = force_kcal * _KCAL_TO_KJ
                        ag = custom_atom_groups[key]
                        custom_lines.append(
                            f"# {key}: {force_kcal} kcal/mol/\u00c5\u00b2 = {force_kj:.2f} kJ/mol/nm\u00b2\n"
                        )
                        for atom in ag:
                            custom_lines.append(f"{atom.index}  {force_kj:.4f}\n")
                        self.logger.info(
                            f"  custom_pos_stage{stage_i}.txt: {len(ag)} atoms for '{key}' "
                            f"@ {force_kj:.2f} kJ/mol/nm\u00b2"
                        )
                    custom_path = restraints_dir / f"custom_pos_stage{stage_i}.txt"
                    with open(custom_path, "w") as fh:
                        fh.writelines(custom_lines)
                    result["custom_pos_per_stage"][stage_i] = custom_path

        except ImportError:
            self.logger.warning(
                "MDAnalysis not available — using PDB name-based heuristic for "
                "prot_pos.txt / lipid_pos.txt.  Install MDAnalysis for accurate "
            )
            if needs_prot or needs_lipid:
                self._generate_openmm_restraints_fallback(
                    system_pdb,
                    restraints_dir,
                    max_forces,
                    result,
                    needs_prot,
                    needs_lipid,
                )
            if all_custom_keys:
                self.logger.error(
                    "MDAnalysis is required to generate restraints for custom "
                    f"categories ({', '.join(sorted(all_custom_keys))}). "
                )

        return result

    def _generate_openmm_restraints_fallback(
        self,
        system_pdb: Path,
        restraints_dir: Path,
        max_forces: Dict[str, float],
        result: Dict[str, Optional[Path]],
        needs_prot: bool,
        needs_lipid: bool,
    ) -> None:
        """PDB name-based fallback for prot_pos.txt / lipid_pos.txt (no MDAnalysis)."""
        _BB_NAMES = frozenset({"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA"})
        _PROT_RESNAMES = frozenset(
            {
                "ALA",
                "ARG",
                "ASN",
                "ASP",
                "CYS",
                "GLN",
                "GLU",
                "GLY",
                "HIS",
                "ILE",
                "LEU",
                "LYS",
                "MET",
                "PHE",
                "PRO",
                "SER",
                "THR",
                "TRP",
                "TYR",
                "VAL",
                "HSE",
                "HSD",
                "HSP",
                "CYX",
                "HIE",
                "HID",
                "HIP",
                "ASH",
                "GLH",
                "LYN",
                "TYM",
                "CYM",
                "SEP",
                "T2P",
                "ACE",
                "NHE",
                "NME",
                "COO",
            }
        )
        _LIPID_RESNAMES = frozenset(
            {
                "POPC",
                "POPE",
                "POPS",
                "DPPC",
                "DMPC",
                "DOPC",
                "DSPC",
                "CHOL",
                "CHOLEST",
                "PC",
                "PE",
                "PS",
                "PA",
                "PG",
                "PI",
                "SM",
                "CHL",
                "OL",
                "LA",
                "MY",
                "ST",
                "AR",
                "OLE",
                "PAL",
                "STE",
                "LIN",
                "PALM",
                "OLEO",
                "STEROL",
            }
        )
        _LIPID_HEAD_NAMES = frozenset(
            {
                "P",
                "O11",
                "O12",
                "O13",
                "O14",
                "O21",
                "O22",
                "O31",
                "O32",
                "O33",
                "O34",
                "O1P",
                "O2P",
                "O3P",
                "O4P",
                "OP1",
                "OP2",
                "OP3",
                "OP4",
                "N",
                "C11",
                "C12",
                "C13",
                "C14",
                "N31",
                "C32",
                "C33",
                "C34",
                "C35",
                "C1",
                "C2",
                "C3",
                "HN1",
                "HN2",
                "HN3",
                "HO2",
                "HO3",
                "HS",
            }
        )

        prot_lines: List[str] = []
        lipid_lines: List[str] = []
        atom_idx = 0
        try:
            with open(system_pdb, "r") as fh:
                for line in fh:
                    if not line.startswith(("ATOM", "HETATM")):
                        continue
                    atom_name = line[12:16].strip()
                    resname = line[17:20].strip()
                    if needs_prot and resname in _PROT_RESNAMES:
                        label = "BB" if atom_name in _BB_NAMES else "SC"
                        prot_lines.append(f"{atom_idx}  {label}\n")
                    elif needs_lipid and resname in _LIPID_RESNAMES:
                        is_head = atom_name in _LIPID_HEAD_NAMES
                        include = (is_head and max_forces.get("lipid_head", 0) > 0) or (
                            not is_head and max_forces.get("lipid_tail", 0) > 0
                        )
                        if include:
                            lipid_lines.append(f"{atom_idx}\n")
                    atom_idx += 1
        except Exception as exc:
            self.logger.error(
                f"Fallback restraint generation failed reading PDB: {exc}"
            )
            return

        if needs_prot and prot_lines:
            prot_path = restraints_dir / "prot_pos.txt"
            with open(prot_path, "w") as fh:
                fh.writelines(prot_lines)
            result["prot_pos"] = prot_path
            self.logger.info(
                f"  prot_pos.txt (fallback): {len(prot_lines)} protein atoms"
            )

        if needs_lipid and lipid_lines:
            lipid_path = restraints_dir / "lipid_pos.txt"
            with open(lipid_path, "w") as fh:
                fh.writelines(lipid_lines)
            result["lipid_pos"] = lipid_path
            self.logger.info(
                f"  lipid_pos.txt (fallback): {len(lipid_lines)} lipid atoms"
            )

    @staticmethod
    def get_default_stage_params(
        scheme_type: str = "NPT",
        temperature: float = 310.15,
        include_production: bool = False,
    ) -> List[Dict[str, Any]]:
        """
        Return default CHARMM-GUI-style equilibration stages for a membrane protein system.

        Six stages with gradually decreasing positional restraints, following the
        standard CHARMM-GUI membrane equilibration schedule. Suitable as a starting
        point that can be further customised before passing to setup_openmm_equilibration.

        Args:
            scheme_type: Ensemble for all stages (NVT | NPT | NPAT | NPgT).
            temperature: Simulation temperature in Kelvin (default 310.15).
            include_production: When True, append a 50 ns unrestrained production
                stage (default False).

        Returns:
            List of :class:`EquilibrationStage` objects ready to pass to
            setup_openmm_equilibration.  Fields can be edited via attribute
            assignment or the :meth:`~EquilibrationStage.replace` method.

        Example::

            >>> from dataclasses import replace
            >>> stages = OpenMMEquilibrationManager.get_default_stage_params("NPT",
            ...                                                               include_production=True)
            >>> stages[-1].time_ns = 100.0          # mutable attribute set
            >>> stages[0] = stages[0].replace(temperature=303.15)  # immutable copy
            >>> manager = OpenMMEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_openmm_equilibration(stage_params_list=stages)
        """
        valid = {"NVT", "NPT", "NPAT", "NPgT"}
        if scheme_type not in valid:
            raise ValueError(
                f"scheme_type must be one of {sorted(valid)}, got '{scheme_type}'"
            )

        def _stage(
            name,
            time_ns,
            timestep,
            minimize_steps=0,
            dcd_freq=5000,
            **constraints_overrides,
        ):
            base_constraints = {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "water": 0.0,
                "ions": 0.0,
                "other": 0.0,
            }
            base_constraints.update(constraints_overrides)
            return EquilibrationStage(
                name=name,
                ensemble=scheme_type,
                time_ns=time_ns,
                timestep=timestep,
                temperature=temperature,
                minimize_steps=minimize_steps,
                dcd_freq=dcd_freq,
                constraints=base_constraints,
            )

        stages: List[EquilibrationStage] = [
            _stage(
                "Equilibration 1",
                0.125,
                1.0,
                minimize_steps=5000,
                protein_backbone=10.0,
                protein_sidechain=5.0,
                lipid_head=2.5,
            ),
            _stage(
                "Equilibration 2",
                0.125,
                1.0,
                protein_backbone=5.0,
                protein_sidechain=2.5,
                lipid_head=1.0,
            ),
            _stage(
                "Equilibration 3",
                0.125,
                1.0,
                protein_backbone=2.5,
                protein_sidechain=1.0,
                lipid_head=0.5,
            ),
            _stage(
                "Equilibration 4",
                0.25,
                1.0,
                protein_backbone=1.0,
                protein_sidechain=0.5,
            ),
            _stage("Equilibration 5", 0.25, 2.0, protein_backbone=0.5),
            _stage("Equilibration 6", 0.5, 2.0, protein_backbone=0.1),
        ]

        if include_production:
            stages.append(_stage("Production", 50.0, 2.0, dcd_freq=50000))

        return stages

    def _get_config_filename(self, stage_index: int) -> str:
        """Return the .inp filename for a given 1-based stage index."""
        if stage_index <= 6:
            return f"step{stage_index}_equilibration.inp"
        return "step7_production.inp"


class GROMACSEquilibrationManager:
    """Manager for GROMACS equilibration simulations using CHARMM-GUI MDP templates.

    Mirrors the NAMDEquilibrationManager / OpenMMEquilibrationManager API.
    The ``constraints`` dict uses the same keys and **kcal/mol/Å²** units;
    values are converted to kJ/mol/nm² when writing MDP ``define`` lines.

    GROMACS input files are generated from the MDP templates shipped in
    ``equilibration/gromacs/{ensemble}/``.  If only AMBER files are found in
    the working directory, they are automatically converted to GROMACS format
    using **ParmEd** (must be installed).

    Workflow generated:
        1. ``step0_minimization.mdp`` — energy minimisation
        2. ``step6.1_equilibration.mdp`` … ``step6.6_equilibration.mdp``
        3. ``step7_production.mdp``
        4. ``run_equilibration.sh`` — orchestrates ``gmx grompp`` + ``gmx mdrun``

    Position restraints (POSRES) are handled via MDP ``define`` macros, identical
    to the CHARMM-GUI scheme: ``POSRES_FC_BB``, ``POSRES_FC_SC``,
    ``POSRES_FC_LIPID``, ``POSRES_FC_WATER``, ``POSRES_FC_ION``,
    ``POSRES_FC_OTHER``, plus per-key macros for custom MDAnalysis selections.
    Separate ``.itp`` files are generated and included into the matching
    molecule types in ``topol.top``.

    Example::

        from pathlib import Path
        from gatewizard.tools.equilibration import GROMACSEquilibrationManager

        manager = GROMACSEquilibrationManager(Path("/work/dir"))
        stages = GROMACSEquilibrationManager.get_default_stage_params("NPT",
                                                                       include_production=True)
        stages[-1].time_ns = 100.0
        result = manager.setup_gromacs_equilibration(stage_params_list=stages)
        print(result["gromacs_dir"])
    """

    # kJ/mol/nm² per kcal/mol/Å²
    _KCAL_TO_KJ: float = 418.4

    SCHEME_MAPPING: Dict[str, str] = {
        "NVT": "01_NVT",
        "NPT": "02_NPT",
        "NPAT": "03_NPAT",
        "NPgT": "04_NPgT",
    }

    TEMPLATE_MAPPING: Dict[str, str] = {
        "step0_minimization": "step6.0_minimization.mdp",
        "step1": "step6.1_equilibration.mdp",
        "step2": "step6.2_equilibration.mdp",
        "step3": "step6.3_equilibration.mdp",
        "step4": "step6.4_equilibration.mdp",
        "step5": "step6.5_equilibration.mdp",
        "step6": "step6.6_equilibration.mdp",
        "step7_production": "step7_production.mdp",
    }

    # 1-based stage index → template key
    STAGE_INDEX_TO_KEY: Dict[int, str] = {
        0: "step0_minimization",
        1: "step1",
        2: "step2",
        3: "step3",
        4: "step4",
        5: "step5",
        6: "step6",
        7: "step7_production",
    }

    def __init__(self, working_dir: Path, gmx_executable: str = "gmx"):
        self.working_dir = Path(working_dir)
        self.gmx_executable = gmx_executable
        self.templates_dir = (
            Path(__file__).parent.parent.parent / "equilibration" / "gromacs"
        )
        self.logger = get_logger(self.__class__.__name__)

    # ------------------------------------------------------------------
    # File discovery & conversion
    # ------------------------------------------------------------------

    def find_system_files(self) -> Optional[Dict[str, str]]:
        """Auto-detect GROMACS or AMBER system files in *working_dir*.

        Searches for::

            GROMACS-native  → system.gro (or *.gro), topol.top, index.ndx
            AMBER fallback  → system.prmtop + system.inpcrd + system.pdb

        Returns:
            Dict with keys ``gro``, ``top``, ``ndx`` (GROMACS) **or**
            ``prmtop``, ``inpcrd``, ``pdb``, ``bilayer_pdb`` (AMBER).
            Returns ``None`` when neither set is found.
        """
        d = self.working_dir

        # GROMACS native
        gro_files = sorted(d.glob("*.gro"))
        top_file = d / "topol.top"
        ndx_file = d / "index.ndx"
        if gro_files and top_file.exists():
            result: Dict[str, str] = {
                "gro": str(gro_files[0]),
                "top": str(top_file),
            }
            if ndx_file.exists():
                result["ndx"] = str(ndx_file)
            self.logger.info(
                f"Found GROMACS files: gro={gro_files[0].name}, top={top_file.name}"
            )
            return result

        # AMBER fallback
        prmtop = None
        inpcrd = None
        pdb = None
        bilayer_pdb = None
        for f in sorted(d.glob("*.prmtop")):
            prmtop = str(f)
            break
        for f in sorted(d.glob("*.inpcrd")):
            inpcrd = str(f)
            break
        for f in sorted(d.glob("*.pdb")):
            if "lipid" in f.name.lower():
                bilayer_pdb = str(f)
            elif pdb is None:
                pdb = str(f)
        if prmtop and inpcrd:
            result = {"prmtop": prmtop, "inpcrd": inpcrd}
            if pdb:
                result["pdb"] = pdb
            if bilayer_pdb:
                result["bilayer_pdb"] = bilayer_pdb
            self.logger.info(
                f"Found AMBER files in {d}, will convert to GROMACS automatically"
            )
            return result

        self.logger.warning(f"No GROMACS or AMBER system files found in {d}")
        return None

    def convert_from_amber(
        self,
        prmtop: Path,
        inpcrd: Path,
        output_dir: Path,
        bilayer_pdb: Optional[Path] = None,
    ) -> Dict[str, Path]:
        """Convert AMBER topology/coordinates to GROMACS format using **ParmEd**.

        Args:
            prmtop: AMBER topology file (.prmtop).
            inpcrd: AMBER coordinate file (.inpcrd).
            output_dir: Directory where ``system.gro`` and ``topol.top`` are written.
            bilayer_pdb: PDB with CRYST1 record; provides box dimensions when
                the prmtop has no box (IFBOX=0).

        Returns:
            Dict with keys ``gro`` and ``top`` pointing to the written files.

        Raises:
            ImportError: If ParmEd is not installed.
            RuntimeError: If conversion fails.
        """
        try:
            import parmed as pmd  # type: ignore
        except ImportError as exc:
            raise ImportError(
                "ParmEd is required for AMBER→GROMACS conversion. "
            ) from exc

        self.logger.info(f"Converting AMBER files to GROMACS using ParmEd…")
        struct = pmd.load_file(str(prmtop), str(inpcrd))

        # Apply box from bilayer PDB when prmtop lacks one
        if struct.box is None and bilayer_pdb and Path(bilayer_pdb).exists():
            try:
                ref = pmd.load_file(str(bilayer_pdb))
                if ref.box is not None:
                    struct.box = ref.box
                    self.logger.info(f"  Box from bilayer PDB: {struct.box[:3]} Å")
            except Exception as exc:
                self.logger.warning(f"  Could not read box from bilayer PDB: {exc}")

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        gro_path = output_dir / "system.gro"
        top_path = output_dir / "topol.top"
        struct.save(str(gro_path), overwrite=True)
        struct.save(str(top_path), overwrite=True)
        self.logger.info(f"  Written: {gro_path.name}, {top_path.name}")
        return {"gro": gro_path, "top": top_path}

    def generate_index_ndx(
        self,
        gro_path: Path,
        output_path: Path,
        extra_groups: Optional[str] = None,
    ) -> Path:
        """Generate GROMACS index file with standard groups.

        Calls ``gmx make_ndx`` on the GRO file.  An optional string of extra
        ``make_ndx`` commands (e.g. ``"1|2\\nq"`` for a merged group) can be
        supplied via *extra_groups*.

        Args:
            gro_path: GRO structure file.
            output_path: Output ``.ndx`` file path.
            extra_groups: Additional ``make_ndx`` commands (before ``q``).

        Returns:
            Path to the generated index file.
        """
        cmds = extra_groups + "\nq\n" if extra_groups else "q\n"
        try:
            result = subprocess.run(
                [
                    self.gmx_executable,
                    "make_ndx",
                    "-f",
                    str(gro_path),
                    "-o",
                    str(output_path),
                ],
                input=cmds,
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                self.logger.warning(
                    f"gmx make_ndx exited {result.returncode}:\n{result.stderr[-400:]}"
                )
            else:
                self.logger.info(f"  Generated index: {output_path.name}")
        except FileNotFoundError:
            self.logger.warning(
                f"gmx executable '{self.gmx_executable}' not found. "
                "Index file not generated — supply one manually or set gmx_executable."
            )
        return output_path

    @staticmethod
    def _append_charmm_gui_groups(ndx_path: Path) -> None:
        """Append SOLU/MEMB/SOLV/SOLU_MEMB groups to an existing GROMACS index file.

        The MDP templates use CHARMM-GUI-style group names for ``tc_grps`` and
        ``comm_grps``.  ParmEd-converted AMBER systems provide groups named
        ``Protein`` (for the solute) and ``Other`` (for lipids).  This method
        creates the expected aliases so that ``gmx grompp`` can resolve them.

        Groups appended:

        - ``SOLU``     — all ``Protein`` atoms
        - ``MEMB``     — all ``Other`` atoms (lipids)
        - ``SOLV``     — System minus Protein minus Other (water + ions)
        - ``SOLU_MEMB`` — Protein union Other

        Groups that cannot be constructed (e.g. no Protein or Other found) are
        silently skipped so that pure-protein or pure-membrane setups degrade
        gracefully.

        Args:
            ndx_path: Path to an existing GROMACS index (.ndx) file.
        """
        if not ndx_path.exists():
            return

        content = ndx_path.read_text()

        def _parse_group(name: str) -> List[int]:
            m = re.search(
                rf"\[\s*{re.escape(name)}\s*\]\s*\n([\s\d]+?)(?=\[|\Z)",
                content,
                re.DOTALL,
            )
            return list(map(int, m.group(1).split())) if m else []

        protein = _parse_group("Protein")
        other = _parse_group("Other")
        system = _parse_group("System")

        if not system:
            return  # Malformed or empty index — nothing to do

        protein_set = set(protein)
        other_set = set(other)
        system_set = set(system)
        solv_atoms = sorted(system_set - protein_set - other_set)
        solu_memb_atoms = sorted(protein_set | other_set)

        def _format_group(name: str, atoms: List[int]) -> str:
            lines = [f"[ {name} ]"]
            for i in range(0, len(atoms), 15):
                lines.append("  ".join(str(a) for a in atoms[i : i + 15]))
            return "\n".join(lines) + "\n"

        parts: List[str] = []
        if protein:
            parts.append(_format_group("SOLU", protein))
        if other:
            parts.append(_format_group("MEMB", other))
        if solv_atoms:
            parts.append(_format_group("SOLV", solv_atoms))
        if protein and other:
            parts.append(_format_group("SOLU_MEMB", solu_memb_atoms))

        if parts:
            with ndx_path.open("a") as fh:
                for part in parts:
                    fh.write("\n" + part)

    # Residue / molecule-type name sets (aligned with NAMD DEFAULT_SELECTIONS)
    _WATER_MOLNAMES: frozenset = frozenset(
        {
            "SOL",
            "TIP3",
            "TIP3P",
            "HOH",
            "WAT",
            "TIP4",
            "TIP4P",
            "SPC",
            "SPCE",
            "T3P",
            "T4P",
        }
    )
    _ION_MOLNAMES: frozenset = frozenset(
        {
            "NA",
            "CL",
            "K",
            "CA",
            "MG",
            "ZN",
            "FE",
            "CU",
            "SOD",
            "CLA",
            "POT",
            "CAL",
            "MAG",
            "ZIN",
            "IRN",
            "COP",
            "NA+",
            "CL-",
            "K+",
            "CA2+",
            "MG2+",
            "ZN2+",
            "FE2+",
            "FE3+",
            "CU2+",
            "LIT",
            "RUB",
            "CES",
            "BAR",
        }
    )
    _LIPID_MOLNAMES: frozenset = frozenset(
        {
            "POPC",
            "POPE",
            "POPS",
            "DPPC",
            "DMPC",
            "DOPC",
            "DSPC",
            "PC",
            "PE",
            "PS",
            "PA",
            "PG",
            "PI",
            "SM",
            "OL",
            "LA",
            "MY",
            "ST",
            "AR",
            "OLE",
            "PAL",
            "STE",
            "LIN",
            "CHOL",
            "CHL",
            "CHOLEST",
            "PALM",
            "OLEO",
            "STEROL",
        }
    )
    _STD_POSRES_KEYS: frozenset = frozenset(
        {
            "protein_backbone",
            "protein_sidechain",
            "lipid_head",
            "lipid_tail",
            "water",
            "ions",
            "ion",
            "other",
        }
    )
    _PROTEIN_LIPID_KEYS: frozenset = frozenset(
        {
            "protein_backbone",
            "protein_sidechain",
            "lipid_head",
            "lipid_tail",
        }
    )

    @staticmethod
    def _parse_mol_atom_counts_from_top(top_path: Path) -> List[int]:
        """Return the number of atoms in each ``[ moleculetype ]`` block.

        Parses ``[ atoms ]`` sections and counts non-comment lines that start
        with an integer (the local atom index).  Returns one entry per
        molecule type in order of appearance.

        Args:
            top_path: Path to a GROMACS topology file.

        Returns:
            List of atom counts, one per molecule type, e.g. ``[496, 134, 1, 1, 3]``
            for a system with a 496-atom protein, 134-atom lipid, two ion types,
            and water.
        """
        type_defs = GROMACSEquilibrationManager._parse_moleculetype_defs(top_path)
        return [n for _, n in type_defs]

    @staticmethod
    def _parse_moleculetype_defs(top_path: Path) -> List[Tuple[str, int]]:
        """Return ``[(name, n_atoms), ...]`` for each ``[ moleculetype ]`` block."""
        defs: List[Tuple[str, int]] = []
        current_name: Optional[str] = None
        current_count = 0
        in_atoms = False
        in_mol = False
        awaiting_name = False
        for raw in top_path.read_text().splitlines():
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                if section == "moleculetype":
                    if in_mol and current_name is not None and current_count > 0:
                        defs.append((current_name, current_count))
                    current_name = None
                    current_count = 0
                    in_mol = True
                    in_atoms = False
                    awaiting_name = True
                elif section == "atoms" and in_mol:
                    in_atoms = True
                    awaiting_name = False
                else:
                    in_atoms = False
                    if section != "moleculetype":
                        awaiting_name = False
            elif awaiting_name and in_mol and current_name is None:
                current_name = line.split()[0]
                awaiting_name = False
            elif in_atoms:
                parts = line.split()
                if parts and parts[0].isdigit():
                    current_count += 1
        if in_mol and current_name is not None and current_count > 0:
            defs.append((current_name, current_count))
        return defs

    @staticmethod
    def _parse_molecules_section(top_path: Path) -> List[Tuple[str, int]]:
        """Return ``[(moltype_name, count), ...]`` from the ``[ molecules ]`` block."""
        molecules: List[Tuple[str, int]] = []
        in_molecules = False
        for raw in top_path.read_text().splitlines():
            line = raw.split(";")[0].strip()
            if not line:
                continue
            if line.startswith("["):
                section = line.strip("[] ").lower()
                in_molecules = section == "molecules"
                continue
            if in_molecules:
                parts = line.split()
                if len(parts) >= 2 and parts[1].lstrip("+-").isdigit():
                    molecules.append((parts[0], int(parts[1])))
        return molecules

    @classmethod
    def _build_mol_copies(
        cls, top_path: Path
    ) -> Tuple[List[Tuple[str, int]], List[Dict[str, Any]]]:
        """Build per-copy global atom ranges from a topology.

        Returns:
            ``(type_defs, copies)`` where *type_defs* is
            ``[(name, n_atoms), ...]`` and *copies* is a list of dicts with
            keys ``name``, ``n_atoms``, ``start``, ``end``, ``copy_idx``
            (``start``/``end`` are 1-based global atom indices).
            ``copy_idx`` is cumulative per molecule-type name even when the
            same name appears in multiple ``[ molecules ]`` rows.
        """
        type_defs = cls._parse_moleculetype_defs(top_path)
        type_atoms = {name: n for name, n in type_defs}
        molecules = cls._parse_molecules_section(top_path)
        if not molecules and type_defs:
            molecules = [(name, 1) for name, _ in type_defs]

        copies: List[Dict[str, Any]] = []
        g = 1
        name_copy_counters: Dict[str, int] = defaultdict(int)
        for name, count in molecules:
            n_atoms = type_atoms.get(name)
            if n_atoms is None:
                continue
            for _ in range(count):
                copy_idx = name_copy_counters[name]
                name_copy_counters[name] += 1
                copies.append(
                    {
                        "name": name,
                        "n_atoms": n_atoms,
                        "start": g,
                        "end": g + n_atoms - 1,
                        "copy_idx": copy_idx,
                    }
                )
                g += n_atoms
        return type_defs, copies

    @classmethod
    def _classify_moltype(cls, name: str) -> str:
        """Classify a molecule type name as protein|lipid|water|ions|other."""
        key = name.strip().upper()
        bare = key.rstrip("+-0123456789")
        if key in cls._WATER_MOLNAMES or bare in cls._WATER_MOLNAMES:
            return "water"
        if key in cls._ION_MOLNAMES or bare in cls._ION_MOLNAMES:
            return "ions"
        if key in cls._LIPID_MOLNAMES or bare in cls._LIPID_MOLNAMES:
            return "lipid"
        if key in {"PROTEIN", "PROT", "SYSTEM", "MOL", "PROA", "PROBE"}:
            return "protein"
        return "other"

    @classmethod
    def _pick_protein_lipid_molnames(
        cls,
        type_defs: List[Tuple[str, int]],
        molecules: List[Tuple[str, int]],
    ) -> Tuple[Optional[str], Optional[str]]:
        """Choose protein and lipid molecule-type names for POSRES includes.

        ParmEd membrane topologies often name types ``system1`` / ``system2``
        (not Protein/POPC).  When names are unclassified:

        * protein → largest ``n_atoms`` among non-water/non-ion types
        * lipid → most copies among remaining non-water/non-ion types

        Never falls back to an ion type for lipids.
        """
        copy_counts: Dict[str, int] = defaultdict(int)
        for name, count in molecules:
            copy_counts[name] += int(count)

        by_class = {name: cls._classify_moltype(name) for name, _ in type_defs}
        protein = next(
            (n for n, _ in type_defs if by_class[n] == "protein"),
            None,
        )
        lipid = next(
            (n for n, _ in type_defs if by_class[n] == "lipid"),
            None,
        )

        others = [
            (name, natoms, copy_counts.get(name, 1))
            for name, natoms in type_defs
            if by_class[name] == "other"
        ]
        if protein is None and others:
            protein = max(others, key=lambda x: (x[1], x[2]))[0]
            others = [c for c in others if c[0] != protein]
        if lipid is None and others:
            lipid = max(others, key=lambda x: (x[2], x[1]))[0]
        return protein, lipid

    @staticmethod
    def _parse_ndx_group(index_path: Path, group_idx: int) -> List[int]:
        """Return 1-based atom indices from a GROMACS ``.ndx`` group by index."""
        groups: List[List[int]] = []
        current: Optional[List[int]] = None
        for line in Path(index_path).read_text().splitlines():
            if line.startswith("["):
                current = []
                groups.append(current)
            elif current is not None:
                current.extend(int(x) for x in line.split() if x.lstrip("+-").isdigit())
        if group_idx < 0 or group_idx >= len(groups):
            return []
        return groups[group_idx]

    @staticmethod
    def _local_indices_in_copy(
        global_ids: List[int], copy: Dict[str, Any]
    ) -> Set[int]:
        """Map global 1-based indices into local 1-based indices for *copy*."""
        start = int(copy["start"])
        end = int(copy["end"])
        return {g - start + 1 for g in global_ids if start <= int(g) <= end}

    @staticmethod
    def _posres_macro_name(key: str) -> str:
        """Return a C-preprocessor-safe POSRES force-constant macro name."""
        special = {
            "water": "POSRES_FC_WATER",
            "ions": "POSRES_FC_ION",
            "ion": "POSRES_FC_ION",
            "other": "POSRES_FC_OTHER",
            "protein_backbone": "POSRES_FC_BB",
            "protein_sidechain": "POSRES_FC_SC",
            "lipid_head": "POSRES_FC_LIPID",
            "lipid_tail": "POSRES_FC_LIPID",
        }
        if key in special:
            return special[key]
        safe = re.sub(r"[^A-Za-z0-9_]", "_", key).upper()
        if not safe or safe[0].isdigit():
            safe = "X_" + safe
        return f"POSRES_FC_{safe}"

    @staticmethod
    def _safe_posres_filename(text: str) -> str:
        """Sanitize a constraint/moltype name for use in an .itp filename."""
        safe = re.sub(r"[^A-Za-z0-9_]+", "_", text.strip()).strip("_").lower()
        return safe or "custom"

    @staticmethod
    def _write_posre_itp(
        outfile: Path, local_indices: Set[int], macro: str
    ) -> Path:
        """Write a position-restraints .itp using *macro* as the FC placeholder."""
        lines = [
            f"; Position restraints generated by GateWizard ({macro})",
            "[ position_restraints ]",
            ";  i funct       fcx        fcy        fcz",
        ]
        for idx in sorted(local_indices):
            lines.append(f"{idx:6d}     1  {macro} {macro} {macro}")
        outfile.write_text("\n".join(lines) + "\n")
        return outfile

    def _map_selection_to_local_indices(
        self,
        structure_path: Path,
        selection: str,
        copies: List[Dict[str, Any]],
    ) -> Tuple[Dict[str, Set[int]], bool]:
        """Map an MDAnalysis selection to local indices per molecule type.

        Returns:
            ``(local_by_moltype, partial_copies)`` where *local_by_moltype*
            maps molecule-type name → set of 1-based local atom indices, and
            *partial_copies* is True when only some copies of a multi-copy
            type were selected (GROMACS POSRES will still restrain all copies).
        """
        import MDAnalysis as mda  # type: ignore
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            u = mda.Universe(str(structure_path))
            ag = u.select_atoms(selection)

        if len(ag) == 0:
            return {}, False

        n_atoms_total = int(sum(c["n_atoms"] for c in copies))
        lookup: List[Optional[Dict[str, Any]]] = [None] * n_atoms_total
        for copy in copies:
            for local in range(1, copy["n_atoms"] + 1):
                g0 = copy["start"] + local - 2  # 0-based
                if 0 <= g0 < len(lookup):
                    lookup[g0] = {
                        "name": copy["name"],
                        "local": local,
                        "copy_idx": copy["copy_idx"],
                    }

        n_copies: Dict[str, int] = defaultdict(int)
        for copy in copies:
            n_copies[copy["name"]] = max(n_copies[copy["name"]], copy["copy_idx"] + 1)

        selected_by_type_copy: Dict[Tuple[str, int], Set[int]] = defaultdict(set)
        for atom in ag:
            g0 = int(atom.index)
            if g0 < 0 or g0 >= len(lookup) or lookup[g0] is None:
                continue
            info = lookup[g0]
            assert info is not None
            selected_by_type_copy[(info["name"], info["copy_idx"])].add(info["local"])

        local_by_moltype: Dict[str, Set[int]] = defaultdict(set)
        for (molname, _copy_idx), locals_set in selected_by_type_copy.items():
            local_by_moltype[molname].update(locals_set)

        partial = False
        for molname in local_by_moltype:
            total = n_copies.get(molname, 1)
            if total <= 1:
                continue
            hit_copies = {
                cidx for (name, cidx) in selected_by_type_copy if name == molname
            }
            if len(hit_copies) < total:
                partial = True
                continue
            ref = None
            for cidx in range(total):
                s = selected_by_type_copy.get((molname, cidx), set())
                if ref is None:
                    ref = s
                elif s != ref:
                    partial = True
                    break

        return dict(local_by_moltype), partial

    def _write_first_lipid_gro_from_range(
        self,
        gro_path: Path,
        output_path: Path,
        start_global: int,
        end_global: int,
    ) -> int:
        """Write a single-molecule mini-GRO with atoms *start_global*–*end_global*.

        The atoms are renumbered 1…N so that ``gmx genrestr`` outputs local
        position-restraint indices valid within the molecule type's ``[ atoms ]``
        block (indices must be 1..N_mol_atoms in GROMACS).

        Args:
            gro_path: Full-system GRO file.
            output_path: Destination for the mini-GRO.
            start_global: 1-based first atom index to extract (inclusive).
            end_global: 1-based last atom index to extract (inclusive).

        Returns:
            Number of atoms written (0 on failure).
        """
        lines = gro_path.read_text().splitlines()
        if len(lines) < 3:
            return 0
        try:
            n_total = int(lines[1].strip())
        except ValueError:
            return 0
        atom_lines = lines[2 : 2 + n_total]
        box_line = lines[2 + n_total] if len(lines) > 2 + n_total else ""
        lip_lines = atom_lines[start_global - 1 : end_global]
        n_lip = len(lip_lines)
        if n_lip == 0:
            return 0
        out = [lines[0], f"{n_lip:5d}"]
        for i, line in enumerate(lip_lines):
            out.append(f"    1{line[5:15]}{i + 1:5d}{line[20:]}")
        out.append(box_line)
        output_path.write_text("\n".join(out) + "\n")
        return n_lip

    def generate_posres_itp(
        self,
        gro_path: Path,
        index_path: Path,
        output_dir: Path,
        top_path: Optional[Path] = None,
        pdb_path: Optional[Path] = None,
        selections: Optional[Dict[str, str]] = None,
        constraints_max: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        """Generate GROMACS position-restraint .itp files with macro force constants.

        Creates protein/lipid files via index-local (or genrestr) POSRES and
        water / ions / other / custom constraint files via MDAnalysis selections
        mapped to per-moleculetype local indices.

        Protein posres (groups 4=Backbone, 8=SideChain) are remapped to
        **local** indices of the first protein molecule copy.  Full-system
        ``gmx genrestr`` indices are wrong when the protein has multiple
        copies (ParmEd ``system1`` × N) — GROMACS requires 1…n_atoms_local.

        Lipid posres extract the first lipid molecule into a mini-GRO so
        ``gmx genrestr`` emits local indices 1…N.

        Water, ions, other, and custom keys use *selections* (or
        :attr:`NAMDEquilibrationManager.DEFAULT_SELECTIONS`) and write
        ``posre_<key>_<moltype>.itp`` files included into each affected
        molecule type.  GROMACS POSRES apply to **all copies** of a molecule
        type; a warning is logged when the selection hits only a subset.

        Args:
            gro_path: GRO structure file (full system).
            index_path: GROMACS index (.ndx) file.
            output_dir: Directory where .itp files are written.
            top_path: Topology file; used for molecule ranges and includes.
            pdb_path: Optional PDB for MDAnalysis (preferred over GRO).
            selections: ``{constraint_key: mda_selection_string}`` overrides.
            constraints_max: Max force (kcal/mol/Å²) per key across stages;
                keys with force ≤ 0 are skipped for MDA-based posres.

        Returns:
            Dict with ``backbone``, ``sidechain``, ``lipid`` (Path|None),
            ``includes`` (``{moltype: [Path, ...]}``), and ``macros``
            (list of POSRES_FC_* names used for MDA/custom keys).
        """
        result: Dict[str, Any] = {
            "backbone": None,
            "sidechain": None,
            "lipid": None,
            "includes": {},
            "macros": [],
        }
        includes: Dict[str, List[Path]] = defaultdict(list)
        macros_used: List[str] = []

        def _replace_fc(content: str, macro: str) -> str:
            return re.sub(
                r"(\s+1)(\s+1)(\s+1)\s*$",
                f"  {macro}   {macro}   {macro}",
                content,
                flags=re.MULTILINE,
            )

        def _genrestr(
            group_idx: int,
            outfile: Path,
            macro: str,
            gro_override: Optional[Path] = None,
            ndx_override: Optional[Path] = None,
        ) -> Optional[Path]:
            _gro = str(gro_override or gro_path)
            _ndx = str(ndx_override or index_path)
            try:
                proc = subprocess.run(
                    [
                        self.gmx_executable,
                        "genrestr",
                        "-f",
                        _gro,
                        "-n",
                        _ndx,
                        "-o",
                        str(outfile),
                        "-fc",
                        "1",
                        "1",
                        "1",
                    ],
                    input=f"{group_idx}\n",
                    capture_output=True,
                    text=True,
                )
                if proc.returncode != 0 or not outfile.exists():
                    self.logger.warning(
                        f"gmx genrestr failed for group {group_idx} "
                        f"(gro={Path(_gro).name}): {proc.stderr[-300:]}"
                    )
                    return None
                outfile.write_text(_replace_fc(outfile.read_text(), macro))
                self.logger.info(f"  Generated {outfile.name}")
                return outfile
            except FileNotFoundError:
                self.logger.warning(
                    f"gmx executable not found; skipping {outfile.name}"
                )
                return None

        result["backbone"] = None
        result["sidechain"] = None

        if top_path is None:
            for candidate in (
                gro_path.parent / "topol_posres.top",
                gro_path.parent / "topol.top",
            ):
                if candidate.exists():
                    top_path = candidate
                    break

        type_defs: List[Tuple[str, int]] = []
        copies: List[Dict[str, Any]] = []
        molecules_sec: List[Tuple[str, int]] = []
        protein_name: Optional[str] = None
        lipid_name: Optional[str] = None
        if top_path is not None and top_path.exists():
            type_defs, copies = self._build_mol_copies(top_path)
            molecules_sec = self._parse_molecules_section(top_path)
            protein_name, lipid_name = self._pick_protein_lipid_molnames(
                type_defs, molecules_sec
            )
            self.logger.debug(
                f"Topology moltypes: {type_defs}; {len(copies)} molecule copies; "
                f"protein={protein_name!r} lipid={lipid_name!r}"
            )

        prot_copy = next(
            (
                c
                for c in copies
                if protein_name
                and c["name"] == protein_name
                and c["copy_idx"] == 0
            ),
            None,
        )
        if prot_copy is not None and index_path.exists():
            for group_idx, key, outfile, macro in (
                (
                    4,
                    "backbone",
                    output_dir / "posre_protein_backbone.itp",
                    "POSRES_FC_BB",
                ),
                (
                    8,
                    "sidechain",
                    output_dir / "posre_protein_sidechain.itp",
                    "POSRES_FC_SC",
                ),
            ):
                gids = self._parse_ndx_group(index_path, group_idx)
                local = self._local_indices_in_copy(gids, prot_copy)
                if not local:
                    self.logger.warning(
                        f"No {key} atoms inside first {protein_name} copy "
                        f"(atoms {prot_copy['start']}–{prot_copy['end']}); "
                        "falling back to full-system genrestr."
                    )
                    result[key] = _genrestr(group_idx, outfile, macro)
                else:
                    bad = [i for i in local if i < 1 or i > int(prot_copy["n_atoms"])]
                    if bad:
                        self.logger.warning(
                            f"Invalid local {key} indices {bad[:5]}…; "
                            "falling back to full-system genrestr."
                        )
                        result[key] = _genrestr(group_idx, outfile, macro)
                    else:
                        self._write_posre_itp(outfile, local, macro)
                        result[key] = outfile
                        self.logger.info(
                            f"  Generated {outfile.name} "
                            f"({len(local)} local atoms in {protein_name}, {macro})"
                        )
        else:
            result["backbone"] = _genrestr(
                4, output_dir / "posre_protein_backbone.itp", "POSRES_FC_BB"
            )
            result["sidechain"] = _genrestr(
                8, output_dir / "posre_protein_sidechain.itp", "POSRES_FC_SC"
            )
            if result["backbone"] or result["sidechain"]:
                self.logger.warning(
                    "Protein molecule type not resolved; protein POSRES may use "
                    "global indices and fail grompp if the protein has multiple copies."
                )

        if top_path is None or not top_path.exists():
            self.logger.warning(
                "Topology not found; cannot determine lipid atom range. "
                "Falling back to group-12 genrestr (may produce wrong indices)."
            )
            result["lipid"] = _genrestr(
                12, output_dir / "posre_lipid.itp", "POSRES_FC_LIPID"
            )
        else:
            lipid_copy = next(
                (
                    c
                    for c in copies
                    if lipid_name
                    and c["name"] == lipid_name
                    and c["copy_idx"] == 0
                ),
                None,
            )

            if lipid_copy is None:
                self.logger.warning(
                    "No lipid molecule type found; skipping lipid posres."
                )
            else:
                lip_start = int(lipid_copy["start"])
                lip_end = int(lipid_copy["end"])
                n_lip = int(lipid_copy["n_atoms"])
                self.logger.info(
                    f"  Lipid posres: extracting first {lipid_copy['name']} "
                    f"(atoms {lip_start}–{lip_end}, {n_lip} atoms)"
                )
                with tempfile.TemporaryDirectory() as _tmpdir:
                    _tmpdir_p = Path(_tmpdir)
                    mini_gro = _tmpdir_p / "first_lipid.gro"
                    mini_ndx = _tmpdir_p / "first_lipid.ndx"
                    written = self._write_first_lipid_gro_from_range(
                        gro_path, mini_gro, lip_start, lip_end
                    )
                    if written == 0:
                        self.logger.warning(
                            f"Could not extract lipid atoms {lip_start}–{lip_end} "
                            f"from {gro_path.name}; skipping lipid posres."
                        )
                    else:
                        try:
                            subprocess.run(
                                [
                                    self.gmx_executable,
                                    "make_ndx",
                                    "-f",
                                    str(mini_gro),
                                    "-o",
                                    str(mini_ndx),
                                ],
                                input="q\n",
                                capture_output=True,
                                text=True,
                            )
                            result["lipid"] = _genrestr(
                                0,
                                output_dir / "posre_lipid.itp",
                                "POSRES_FC_LIPID",
                                gro_override=mini_gro,
                                ndx_override=mini_ndx,
                            )
                            if result["lipid"] is not None:
                                includes[str(lipid_copy["name"])].append(
                                    result["lipid"]
                                )
                        except FileNotFoundError:
                            self.logger.warning(
                                "gmx executable not found; skipping lipid posres."
                            )

        max_forces = {
            k: float(v) for k, v in (constraints_max or {}).items() if float(v) > 0
        }
        mda_keys = sorted(k for k in max_forces if k not in self._PROTEIN_LIPID_KEYS)
        if "ion" in mda_keys and "ions" in mda_keys:
            mda_keys = [k for k in mda_keys if k != "ion"]

        if mda_keys and copies:
            structure = (
                Path(pdb_path) if pdb_path and Path(pdb_path).exists() else gro_path
            )
            resolved = dict(NAMDEquilibrationManager.DEFAULT_SELECTIONS)
            if selections:
                resolved.update(selections)
            if "ion" in resolved and "ions" not in resolved:
                resolved["ions"] = resolved["ion"]

            try:
                for key in mda_keys:
                    sel_key = "ions" if key == "ion" else key
                    sel = resolved.get(sel_key) or resolved.get(key)
                    if not sel:
                        self.logger.warning(
                            f"No MDAnalysis selection for constraint '{key}'; "
                            "skipping GROMACS posres for this key."
                        )
                        continue
                    try:
                        local_by_mt, partial = self._map_selection_to_local_indices(
                            structure, sel, copies
                        )
                    except Exception as exc:
                        self.logger.warning(
                            f"MDAnalysis selection failed for '{key}' "
                            f"({sel!r}): {exc}"
                        )
                        continue
                    if not local_by_mt:
                        self.logger.warning(
                            f"Selection for '{key}' matched no atoms; skipping."
                        )
                        continue
                    if partial:
                        self.logger.warning(
                            f"Constraint '{key}' selected only a subset of copies "
                            "of one or more molecule types. GROMACS POSRES will "
                            "restrain ALL copies of those types."
                        )
                    macro = self._posres_macro_name(key)
                    if macro not in macros_used:
                        macros_used.append(macro)
                    key_safe = self._safe_posres_filename(key)
                    for molname, local_indices in local_by_mt.items():
                        mol_safe = self._safe_posres_filename(molname)
                        outfile = output_dir / f"posre_{key_safe}_{mol_safe}.itp"
                        self._write_posre_itp(outfile, local_indices, macro)
                        includes[molname].append(outfile)
                        self.logger.info(
                            f"  Generated {outfile.name} "
                            f"({len(local_indices)} local atoms, {macro})"
                        )
            except ImportError:
                self.logger.warning(
                    "MDAnalysis not available; skipping water/ions/other/"
                    "custom GROMACS position restraints."
                )
        elif mda_keys and not copies:
            self.logger.warning(
                "Cannot map water/ions/other/custom restraints without a "
                "parseable topology [molecules] section."
            )

        result["includes"] = {k: list(v) for k, v in includes.items()}
        result["macros"] = macros_used
        return result

    def _add_posres_to_topology(
        self,
        topol_path: Path,
        posres_files: Dict[str, Any],
        output_path: Optional[Path] = None,
    ) -> Path:
        """Insert ``#ifdef POSRES`` blocks into a GROMACS topology file.

        Inserts protein backbone/sidechain includes at the end of the protein
        molecule type (detected by name or ParmEd ``systemN`` heuristics),
        lipid includes into the lipid molecule type, and any MDA/custom
        includes into the named molecule types listed in
        ``posres_files["includes"]``.

        Args:
            topol_path: Source ``topol.top`` file.
            posres_files: Dict from :meth:`generate_posres_itp`.
            output_path: Destination for the modified topology.

        Returns:
            Path to the written topology.
        """
        if output_path is None:
            output_path = topol_path.parent / "topol_posres.top"

        content = topol_path.read_text()

        moltype_starts = [
            m.start()
            for m in re.finditer(r"^\[ moleculetype \]", content, re.MULTILINE)
        ]
        system_match = re.search(r"^\[ system \]", content, re.MULTILINE)
        end_cap = system_match.start() if system_match else len(content)

        name_to_insert_pos: Dict[str, int] = {}
        ordered_names: List[str] = []
        for i, start in enumerate(moltype_starts):
            block_end = (
                moltype_starts[i + 1] if i + 1 < len(moltype_starts) else end_cap
            )
            block = content[start:block_end]
            name = None
            for raw in block.splitlines()[1:]:
                line = raw.split(";")[0].strip()
                if not line or line.startswith("["):
                    if line.startswith("["):
                        break
                    continue
                name = line.split()[0]
                break
            if name:
                ordered_names.append(name)
                name_to_insert_pos[name] = block_end

        if not moltype_starts:
            shutil.copy2(topol_path, output_path)
            self.logger.warning(
                "Could not locate molecule boundaries in topology; "
                "posres includes NOT added. Add them manually."
            )
            return output_path

        type_defs = self._parse_moleculetype_defs(topol_path)
        molecules_sec = self._parse_molecules_section(topol_path)
        protein_name, lipid_name = self._pick_protein_lipid_molnames(
            type_defs, molecules_sec
        )

        def _block_for(paths: List[Path]) -> str:
            if not paths:
                return ""
            seen = set()
            lines = ["\n#ifdef POSRES"]
            for p in paths:
                name = Path(p).name
                if name in seen:
                    continue
                seen.add(name)
                lines.append(f'#include "{name}"')
            lines.append("#endif\n")
            return "\n".join(lines) if len(lines) > 2 else ""

        inserts: Dict[int, List[str]] = defaultdict(list)

        prot_paths: List[Path] = []
        if posres_files.get("backbone"):
            prot_paths.append(Path(posres_files["backbone"]))
        if posres_files.get("sidechain"):
            prot_paths.append(Path(posres_files["sidechain"]))
        if prot_paths and ordered_names:
            prot_target = protein_name if protein_name in name_to_insert_pos else None
            if prot_target is None:
                prot_target = next(
                    (
                        n
                        for n in ordered_names
                        if self._classify_moltype(n) == "protein"
                    ),
                    ordered_names[0],
                )
            blk = _block_for(prot_paths)
            if blk:
                inserts[name_to_insert_pos[prot_target]].append(blk)

        if posres_files.get("lipid"):
            lip_path = Path(posres_files["lipid"])
            lipid_names: List[str] = []
            if lipid_name and lipid_name in name_to_insert_pos:
                lipid_names = [lipid_name]
            else:
                lipid_names = [
                    n for n in ordered_names if self._classify_moltype(n) == "lipid"
                ]
            for n in lipid_names:
                already = {
                    Path(p).name
                    for p in (posres_files.get("includes") or {}).get(n, [])
                }
                if lip_path.name in already:
                    continue
                blk = _block_for([lip_path])
                if blk:
                    inserts[name_to_insert_pos[n]].append(blk)

        for molname, paths in (posres_files.get("includes") or {}).items():
            if molname not in name_to_insert_pos:
                self.logger.warning(
                    f"Molecule type '{molname}' not found in topology; "
                    f"cannot include {[Path(p).name for p in paths]}"
                )
                continue
            skip_names = set()
            if posres_files.get("backbone"):
                skip_names.add(Path(posres_files["backbone"]).name)
            if posres_files.get("sidechain"):
                skip_names.add(Path(posres_files["sidechain"]).name)
            filtered = [Path(p) for p in paths if Path(p).name not in skip_names]
            blk = _block_for(filtered)
            if blk:
                inserts[name_to_insert_pos[molname]].append(blk)

        new_content = content
        for pos in sorted(inserts.keys(), reverse=True):
            block = "".join(inserts[pos])
            new_content = new_content[:pos] + block + new_content[pos:]

        output_path.write_text(new_content)
        self.logger.info(f"  Topology with posres: {output_path.name}")
        return output_path

    # ------------------------------------------------------------------
    # Default stage parameters
    # ------------------------------------------------------------------

    @staticmethod
    def get_default_stage_params(
        scheme_type: str = "NPT",
        temperature: float = 310.15,
        include_production: bool = False,
    ) -> List["EquilibrationStage"]:
        """Return default CHARMM-GUI-style GROMACS equilibration stages.

        Six stages with gradually decreasing positional restraints, matching
        the standard CHARMM-GUI membrane equilibration schedule (kJ/mol/nm²
        equivalents are computed from kcal/mol/Å² at run time).

        Args:
            scheme_type: Ensemble for all stages (NVT | NPT | NPAT | NPgT).
            temperature: Simulation temperature in Kelvin (default 310.15 K).
            include_production: When True, append a 50 ns unrestrained
                production stage.

        Returns:
            List of :class:`EquilibrationStage` objects ready to pass to
            :meth:`setup_gromacs_equilibration`.

        Example::

            >>> stages = GROMACSEquilibrationManager.get_default_stage_params("NPT")
            >>> stages[-1].time_ns = 5.0
            >>> manager = GROMACSEquilibrationManager(Path("/work"))
            >>> result = manager.setup_gromacs_equilibration(stage_params_list=stages)
        """
        valid = {"NVT", "NPT", "NPAT", "NPgT"}
        if scheme_type not in valid:
            raise ValueError(
                f"scheme_type must be one of {sorted(valid)}, got '{scheme_type}'"
            )

        def _stage(name, time_ns, timestep, minimize_steps=0, **constraints_overrides):
            base = {
                "protein_backbone": 0.0,
                "protein_sidechain": 0.0,
                "lipid_head": 0.0,
                "lipid_tail": 0.0,
                "water": 0.0,
                "ions": 0.0,
                "other": 0.0,
            }
            base.update(constraints_overrides)
            return EquilibrationStage(
                name=name,
                ensemble=scheme_type,
                time_ns=time_ns,
                timestep=timestep,
                temperature=temperature,
                minimize_steps=minimize_steps,
                constraints=base,
            )

        stages: List[EquilibrationStage] = [
            _stage(
                "Minimization",
                0.0,
                1.0,
                minimize_steps=5000,
                protein_backbone=10.0,
                protein_sidechain=5.0,
                lipid_head=2.5,
            ),
            _stage(
                "Equilibration 1",
                0.125,
                1.0,
                protein_backbone=10.0,
                protein_sidechain=5.0,
                lipid_head=2.5,
            ),
            _stage(
                "Equilibration 2",
                0.125,
                1.0,
                protein_backbone=5.0,
                protein_sidechain=2.5,
                lipid_head=1.0,
            ),
            _stage(
                "Equilibration 3",
                0.125,
                1.0,
                protein_backbone=2.5,
                protein_sidechain=1.0,
                lipid_head=1.0,
            ),
            _stage(
                "Equilibration 4",
                0.25,
                1.0,
                protein_backbone=1.0,
                protein_sidechain=0.5,
            ),
            _stage("Equilibration 5", 0.25, 2.0, protein_backbone=0.5),
            _stage("Equilibration 6", 0.5, 2.0, protein_backbone=0.1),
        ]
        if include_production:
            stages.append(_stage("Production", 50.0, 2.0))
        return stages

    # ------------------------------------------------------------------
    # MDP file generation
    # ------------------------------------------------------------------

    @staticmethod
    def _topology_posres_macros(posres_files: Optional[Dict[str, Any]]) -> Set[str]:
        """Return POSRES_FC_* macros actually referenced by generated ``.itp`` files."""
        if not posres_files:
            return set()
        macros: Set[str] = set(posres_files.get("macros") or [])
        if posres_files.get("backbone"):
            macros.add("POSRES_FC_BB")
        if posres_files.get("sidechain"):
            macros.add("POSRES_FC_SC")
        if posres_files.get("lipid"):
            macros.add("POSRES_FC_LIPID")
        return macros

    def generate_mdp_file(
        self,
        stage_name: str,
        stage_params: Dict[str, Any],
        stage_index: int,
        scheme_type: str,
        used_posres_macros: Optional[Set[str]] = None,
    ) -> str:
        """Generate GROMACS MDP file content for a single equilibration stage.

        Loads the appropriate template from ``equilibration/gromacs/{ensemble}/``
        and substitutes runtime parameters including force constants converted
        from kcal/mol/Å² to kJ/mol/nm².

        Args:
            stage_name: Human-readable stage label (used for logging).
            stage_params: Stage parameter dict (see :meth:`setup_gromacs_equilibration`).
            stage_index: 0-based for minimization, 1-based for equilibration steps,
                7 for production.
            scheme_type: Ensemble type (NVT, NPT, NPAT, or NPgT).
            used_posres_macros: POSRES_FC_* macros present in the topology. Unused
                ``-D`` defines are stripped to avoid grompp warnings.

        Returns:
            String content of the generated ``.mdp`` file.
        """
        template_key = self.STAGE_INDEX_TO_KEY.get(stage_index, "step7_production")
        scheme_folder = self.SCHEME_MAPPING[scheme_type]
        template_filename = self.TEMPLATE_MAPPING[template_key]
        template_path = self.templates_dir / scheme_folder / template_filename

        if not template_path.exists():
            raise FileNotFoundError(f"GROMACS MDP template not found: {template_path}")

        content = template_path.read_text()

        temperature = float(stage_params.get("temperature", 310.15))
        timestep_fs = float(stage_params.get("timestep", 2.0))
        dt_ps = timestep_fs / 1000.0
        time_ns = float(stage_params.get("time_ns", 0.5))
        nsteps = max(1, int(round(time_ns * 1_000_000 / timestep_fs)))

        is_minimization = stage_index == 0
        is_production = stage_index == 7
        mini_nsteps = (
            int(stage_params.get("minimize_steps", 5000)) if is_minimization else 5000
        )

        # Force constants (kcal/mol/Å² → kJ/mol/nm²)
        C = self._KCAL_TO_KJ
        constraints = stage_params.get("constraints", {})
        fc_bb = float(constraints.get("protein_backbone", 0.0)) * C
        fc_sc = float(constraints.get("protein_sidechain", 0.0)) * C
        fc_lip = (
            max(
                float(constraints.get("lipid_head", 0.0)),
                float(constraints.get("lipid_tail", 0.0)),
            )
            * C
        )
        fc_water = float(constraints.get("water", 0.0)) * C
        fc_ion = (
            max(
                float(constraints.get("ions", 0.0)),
                float(constraints.get("ion", 0.0)),
            )
            * C
        )
        fc_other = float(constraints.get("other", 0.0)) * C
        custom_fcs: Dict[str, float] = {}
        for key, val in constraints.items():
            if key in self._STD_POSRES_KEYS:
                continue
            fv = float(val) * C
            if fv > 0:
                custom_fcs[self._posres_macro_name(key)] = fv
        # Dihedral restraint: ~SC*0.5 when SC>0, else BB*0.1 (matching CHARMM-GUI pattern)
        fc_dih = fc_sc * 0.5 if fc_sc > 0 else (fc_bb * 0.1 if fc_bb > 0 else 0.0)
        any_posres = (
            fc_bb > 0
            or fc_sc > 0
            or fc_lip > 0
            or fc_water > 0
            or fc_ion > 0
            or fc_other > 0
            or bool(custom_fcs)
        )

        # ------ substitute nsteps / dt ------
        content = re.sub(
            r"nsteps\s*=\s*\d+",
            f"nsteps                  = {mini_nsteps if is_minimization else nsteps}",
            content,
        )
        content = re.sub(r"(?m)^(dt\s*=\s*)[\d.]+", rf"\g<1>{dt_ps:.3f}", content)

        # ------ substitute temperature ------
        # ref_t may have multiple values e.g. "303.15 303.15 303.15"
        content = re.sub(
            r"(ref_t\s*=\s*)[\d.\s]+",
            lambda m: m.group(1)
            + "  ".join(
                [f"{temperature:.2f}"] * m.group(0).count("303")
                if "303" in m.group(0)
                else [f"{temperature:.2f}"]
            ),
            content,
        )
        # Simpler ref_t replacement: count the number of groups from tc_grps
        tc_grps_match = re.search(r"tc_grps\s*=\s*(.+)", content)
        if tc_grps_match:
            n_groups = len(tc_grps_match.group(1).split())
            ref_t_str = "  ".join([f"{temperature:.2f}"] * n_groups)
            content = re.sub(
                r"ref_t\s*=\s*[\d.\s]+",
                f"ref_t                   = {ref_t_str}\n",
                content,
            )
            tau_t_match = re.search(r"tau_t\s*=\s*([\d.\s]+)", content)
            if tau_t_match:
                tau_vals = tau_t_match.group(1).split()
                tau_str = "  ".join([tau_vals[0]] * n_groups)
                content = re.sub(
                    r"tau_t\s*=\s*[\d.\s]+",
                    f"tau_t                   = {tau_str}\n",
                    content,
                )
        content = re.sub(
            r"gen-temp\s*=\s*[\d.]+",
            f"gen-temp                = {temperature:.2f}",
            content,
        )

        # ------ nstxout-compressed ------
        dcd_freq_default = 50000 if is_production else 5000
        dcd_freq = int(stage_params.get("dcd_freq", dcd_freq_default))
        content = re.sub(
            r"nstxout-compressed\s*=\s*\d+",
            f"nstxout-compressed      = {dcd_freq}",
            content,
        )

        # ------ force constants in define line ------
        def _sub_std_macros(text: str) -> str:
            text = re.sub(r"POSRES_FC_BB=[\d.]+", f"POSRES_FC_BB={fc_bb:.1f}", text)
            text = re.sub(r"POSRES_FC_SC=[\d.]+", f"POSRES_FC_SC={fc_sc:.1f}", text)
            text = re.sub(
                r"POSRES_FC_LIPID=[\d.]+", f"POSRES_FC_LIPID={fc_lip:.1f}", text
            )
            text = re.sub(
                r"POSRES_FC_WATER=[\d.]+", f"POSRES_FC_WATER={fc_water:.1f}", text
            )
            text = re.sub(r"POSRES_FC_ION=[\d.]+", f"POSRES_FC_ION={fc_ion:.1f}", text)
            text = re.sub(
                r"POSRES_FC_OTHER=[\d.]+", f"POSRES_FC_OTHER={fc_other:.1f}", text
            )
            text = re.sub(r"DDIHRES_FC=[\d.]+", f"DDIHRES_FC={fc_dih:.1f}", text)
            return text

        content = _sub_std_macros(content)

        def _ensure_macro_on_define(text: str, macro: str, value: float) -> str:
            """Ensure ``-DMACRO=value`` appears on the define line."""
            token = f"-D{macro}="
            if re.search(rf"{re.escape(token)}[\d.]+", text):
                return re.sub(
                    rf"{re.escape(token)}[\d.]+", f"{token}{value:.1f}", text
                )
            if re.search(r"(?m)^define\s*=", text):
                return re.sub(
                    r"(?m)^(define\s*=.*)$",
                    rf"\1 -D{macro}={value:.1f}",
                    text,
                    count=1,
                )
            return text

        # Only ensure macros that the topology actually references (avoids
        # grompp "defined but were not used" warnings for WATER/OTHER/etc.).
        std_extra_macros = (
            ("POSRES_FC_WATER", fc_water),
            ("POSRES_FC_ION", fc_ion),
            ("POSRES_FC_OTHER", fc_other),
        )
        for macro, value in std_extra_macros:
            if used_posres_macros is None or macro in used_posres_macros:
                content = _ensure_macro_on_define(content, macro, value)

        for macro, value in custom_fcs.items():
            if used_posres_macros is None or macro in used_posres_macros:
                content = _ensure_macro_on_define(content, macro, value)

        if not is_minimization:
            # Keep POSRES only when at least one force constant is non-zero.
            if not any_posres:
                content = re.sub(r"^define\s*=.*\n", "", content, flags=re.MULTILINE)
            elif not re.search(r"(?m)^define\s*=", content):
                # CHARMM-GUI production templates may omit the define line;
                # inject one so user-requested production restraints are active.
                std_tokens = []
                for macro, value in (
                    ("POSRES_FC_BB", fc_bb),
                    ("POSRES_FC_SC", fc_sc),
                    ("POSRES_FC_LIPID", fc_lip),
                    ("POSRES_FC_WATER", fc_water),
                    ("POSRES_FC_ION", fc_ion),
                    ("POSRES_FC_OTHER", fc_other),
                ):
                    if used_posres_macros is None or macro in used_posres_macros:
                        std_tokens.append(f"-D{macro}={value:.1f}")
                for macro, value in sorted(custom_fcs.items()):
                    if used_posres_macros is None or macro in used_posres_macros:
                        std_tokens.append(f"-D{macro}={value:.1f}")
                extra = (" " + " ".join(std_tokens)) if std_tokens else ""
                content = f"define                  = -DPOSRES{extra}\n" + content

        # GateWizard does not add dihedral restraints to the topology, so the
        # DIHRES / DIHRES_FC macros would always be unused and cause a warning.
        # Strip them from the define line unconditionally.
        content = re.sub(r"\s*-DDIHRES\b", "", content)
        content = re.sub(r"\s*-DDIHRES_FC=[\d.]+", "", content)

        # Drop POSRES_FC_* defines that are not referenced by any generated ITP.
        if used_posres_macros is not None:
            for match in re.findall(r"-D(POSRES_FC_[A-Za-z0-9_]+)=[\d.]+", content):
                if match not in used_posres_macros:
                    content = re.sub(rf"\s*-D{re.escape(match)}=[\d.]+", "", content)

        self.logger.debug(
            f"Stage {stage_index} ({stage_name}): T={temperature:.2f}K, "
            f"nsteps={nsteps}, dt={dt_ps:.3f}ps, "
            f"fc_bb={fc_bb:.1f} fc_sc={fc_sc:.1f} fc_lip={fc_lip:.1f} "
            f"fc_water={fc_water:.1f} fc_ion={fc_ion:.1f} "
            f"fc_other={fc_other:.1f} kJ/mol/nm²"
        )
        return content

    def _get_mdp_filename(self, stage_index: int) -> str:
        """Return the output ``.mdp`` filename for a given stage index.

        The naming matches the NAMD / OpenMM convention used by GateWizard:
        ``step0_minimization.mdp``, ``step1_equilibration.mdp`` …
        ``step6_equilibration.mdp``, ``step7_production.mdp``.
        (The template files on disk keep their original CHARMM-GUI names such
        as ``step6.1_equilibration.mdp``.)
        """
        if stage_index == 0:
            return "step0_minimization.mdp"
        if 1 <= stage_index <= 6:
            return f"step{stage_index}_equilibration.mdp"
        return "step7_production.mdp"

    # ------------------------------------------------------------------
    # Run script
    # ------------------------------------------------------------------

    def generate_run_script(
        self,
        gromacs_dir: Path,
        gro_name: str,
        top_name: str,
        ndx_name: Optional[str],
        n_stages: int,
        gmx_executable: str = "gmx",
        gmxrc_path: Optional[str] = None,
        cpu_cores: Optional[int] = None,
        use_gpu: bool = False,
        gpu_id: int = 0,
        num_gpus: int = 1,
    ) -> Path:
        """Generate a bash run script for the full GROMACS equilibration protocol.

        Optionally sources a GROMACS ``GMXRC`` (typical for ``/usr/local/gromacs``
        installs). Conda / absolute ``gmx`` binaries usually do not need this.

        Args:
            gromacs_dir: Output directory containing all input files.
            gro_name: Name of the input GRO file.
            top_name: Name of the topology file (``topol_posres.top`` when posres).
            ndx_name: Name of the index file (``None`` → not passed to grompp).
            n_stages: Number of equilibration stages (6 by default).
            gmx_executable: GROMACS executable name/path.
            gmxrc_path: Optional path to ``GMXRC`` to ``source`` before running.
            cpu_cores: OpenMP thread count for ``mdrun -ntomp``.
            use_gpu: When True, enable GPU offload (``-nb gpu -pme gpu``).
            gpu_id: First GPU device index for ``-gpu_id``.
            num_gpus: Number of consecutive GPU devices starting at ``gpu_id``.

        Returns:
            Path to the written ``run_equilibration.sh``.
        """
        # Dynamics stages may use GPU offload; energy minimisation cannot
        # (PME GPU rejects non-dynamical integrators such as steep/cg).
        md_flags = _gromacs_mdrun_resource_flags(
            cpu_cores=cpu_cores,
            use_gpu=use_gpu,
            gpu_id=gpu_id,
            num_gpus=num_gpus,
        )
        em_flags = _gromacs_mdrun_resource_flags(
            cpu_cores=cpu_cores,
            use_gpu=False,
            gpu_id=gpu_id,
            num_gpus=num_gpus,
        )
        mdrun_md = f" {md_flags}" if md_flags else ""
        mdrun_em = f" {em_flags}" if em_flags else ""
        gpu_info = _gromacs_gpu_info(use_gpu, gpu_id, num_gpus)
        cores_label = int(cpu_cores or 0) if cpu_cores else "auto"
        lines = [
            "#!/bin/bash",
            "## GROMACS Equilibration Run Script",
            "## Generated by GateWizard — run from the directory containing this file",
            "",
        ]
        if gmxrc_path:
            lines += [
                "# Source GROMACS environment (system install)",
                f'source "{gmxrc_path}"',
                "",
            ]
        lines += [
            f'GMX="{gmx_executable}"',
            f'GRO="{gro_name}"',
            f'TOP="{top_name}"',
        ]
        if ndx_name:
            lines.append(f'NDX="{ndx_name}"')
        lines += [
            "",
            'echo "Starting GROMACS equilibration protocol…"',
            f'echo "Resources: {cores_label} CPU threads (OpenMP), GPU: {gpu_info}"',
            'echo "Note: step0 energy minimisation runs on CPU (GROMACS PME GPU does not support steep/cg)."',
            "",
            GROMACS_RESUME_SHELL,
            "",
            "# --- Step 0: Energy minimisation (CPU only) ---",
        ]
        grompp_ndx = "-n ${NDX}" if ndx_name else ""
        lines += [
            'if [ "$RESUME" = "1" ] && _gw_gromacs_stage_done "step0_minimization"; then',
            '  echo "RESUME: skipping step0_minimization"',
            "else",
            f"  ${{GMX}} grompp -f step0_minimization.mdp -o step0_minimization.tpr \\",
            f"      -c ${{GRO}} -r ${{GRO}} -p ${{TOP}} {grompp_ndx} -maxwarn 2",
            f"  ${{GMX}} mdrun -v{mdrun_em} -s step0_minimization.tpr -deffnm step0_minimization || {{ echo 'Minimisation failed'; exit 1; }}",
            "fi",
            "",
            "# --- Steps 1–6: Equilibration ---",
        ]
        for i in range(1, n_stages + 1):
            prev = "step0_minimization" if i == 1 else f"step{i - 1}_equilibration"
            curr = f"step{i}_equilibration"
            lines += [
                f'if [ "$RESUME" = "1" ] && _gw_gromacs_stage_done "{curr}"; then',
                f'  echo "RESUME: skipping {curr}"',
                "else",
                f"  ${{GMX}} grompp -f {curr}.mdp -o {curr}.tpr \\",
                f"      -c {prev}.gro -r ${{GRO}} -p ${{TOP}} {grompp_ndx} -maxwarn 2",
                f"  ${{GMX}} mdrun -v{mdrun_md} -s {curr}.tpr -deffnm {curr} || {{ echo 'Stage {i} failed'; exit 1; }}",
                "fi",
                "",
            ]
        lines += [
            "# --- Step 7: Production ---",
            'if [ "$RESUME" = "1" ] && _gw_gromacs_stage_done "step7_production"; then',
            '  echo "RESUME: skipping step7_production"',
            "else",
            f"  ${{GMX}} grompp -f step7_production.mdp -o step7_production.tpr \\",
            f"      -c step{n_stages}_equilibration.gro -p ${{TOP}} {grompp_ndx} -maxwarn 2",
            f"  ${{GMX}} mdrun -v{mdrun_md} -s step7_production.tpr -deffnm step7_production || {{ echo 'Production failed'; exit 1; }}",
            "fi",
            "",
            'echo "GROMACS equilibration complete."',
        ]
        script_path = gromacs_dir / "run_equilibration.sh"
        script_path.write_text("\n".join(lines) + "\n")
        script_path.chmod(0o755)
        self.logger.info(f"Run script: {script_path.name}")
        return script_path

    # ------------------------------------------------------------------
    # COM restraint generation (colvars)
    # ------------------------------------------------------------------

    def generate_com_colvars_config(
        self,
        pdb_path: Path,
        output_file: Path,
        com_restraint_k: float = 10.0,
        add_rotation_restraint: bool = True,
        rotation_restraint_k: float = 2000.0,
        selection: str = "name CA",
        rotation_ref_positions_mode: str = "auto",
        ref_positions_file: Optional[str] = None,
        ref_positions_col: Optional[str] = None,
        ref_positions_col_value: Optional[float] = None,
        ref_base_dir: Optional[Path] = None,
    ) -> Optional[Path]:
        """Generate a GROMACS Colvars configuration to restrain protein COM.

        Prevents lateral drift **and** (optionally) rotation by applying
        harmonic restraints on the geometric center of selected atoms rather
        than per-atom positional restraints.

        GROMACS 2021+ ships the Colvars library.  To activate, add
        ``colvars-active = yes`` and ``colvars-configfile = com_restraint.dat``
        to the MDP file.

        The generated file uses ``distanceZ`` collective variables (translation)
        plus an ``orientation`` CV (rotation).  The initial center and reference
        positions are computed from *pdb_path*.

        .. note::

            **GROMACS Colvars only accepts XYZ format for** ``refPositionsFile``
            (see section 3.7.3 of the GROMACS Colvars manual).  PDB files are
            **not** supported.  When *rotation_ref_positions_mode* is ``"auto"``
            or ``"refPositionsFile"`` and no *ref_positions_file* is given, a
            ``.xyz`` reference file is automatically generated alongside
            *output_file* and referenced from the Colvars config.

        Args:
            pdb_path: System PDB file used to compute initial COM.
            output_file: Path for the written Colvars ``.dat`` config.
            com_restraint_k: Translation force constant in kcal/mol/Å².
            add_rotation_restraint: When True also generate an orientation CV.
            rotation_restraint_k: Rotation force constant in kcal/mol/Å².
            selection: MDAnalysis selection string for the reference atoms
                (default: Cα atoms).
            rotation_ref_positions_mode: Orientation reference mode:
                ``"auto"``, ``"refPositions"``, or ``"refPositionsFile"``.
                ``"auto"`` uses ``refPositionsFile`` with an auto-generated
                XYZ file (required by GROMACS Colvars).  Use
                ``"refPositions"`` to embed coordinates inline instead.
            ref_positions_file: Path to a pre-existing **XYZ** file used when
                mode is ``"refPositionsFile"``.  If not provided, a
                ``.xyz`` file is generated automatically from *pdb_path*.
            ref_positions_col: Not applicable for XYZ files; reserved for
                future PDB-based support in GROMACS Colvars.
            ref_positions_col_value: Optional numeric value paired with
                ``ref_positions_col``.
            ref_base_dir: The directory from which GROMACS runs grompp/mdrun
                (i.e. the directory that contains the MDP and GRO files).
                When provided, the auto-generated XYZ file path written into
                the Colvars config is expressed as a path **relative to
                ref_base_dir**, so GROMACS can resolve it correctly.  If
                ``None``, only the filename is used (works only when the XYZ
                file is in the same directory as grompp's working directory).

        Returns:
            Path to the written config file, or ``None`` on failure.
        """
        try:
            import MDAnalysis as mda  # type: ignore

            u = mda.Universe(str(pdb_path))
            ag = u.select_atoms(selection)
            if len(ag) == 0:
                self.logger.warning(
                    f"No atoms selected by '{selection}' in {pdb_path.name}; "
                    "COM restraint not generated."
                )
                return None

            # Initial centroid (Ångströms)
            com = ag.center_of_geometry()
            x0, y0, z0 = float(com[0]), float(com[1]), float(com[2])
            atom_nums = " ".join(str(int(a.index) + 1) for a in ag)

            # Resolve reference file for the rotation CV.
            # GROMACS Colvars only supports XYZ format for refPositionsFile
            # (section 3.7.3 of the GROMACS Colvars manual);
            # PDB files are NOT supported by GROMACS, unlike NAMD.
            # When the mode is 'auto' or 'refPositionsFile' and no explicit
            # file is provided, we auto-generate a .xyz file alongside the
            # colvars config file.
            mode_norm = str(rotation_ref_positions_mode).strip().lower()
            will_use_file = mode_norm in {"auto", "refpositionsfile", "file"}
            resolved_ref_file = ref_positions_file
            if add_rotation_restraint and will_use_file and not ref_positions_file:
                xyz_path = output_file.with_name(output_file.stem + "_rot_ref.xyz")
                xyz_path.parent.mkdir(parents=True, exist_ok=True)
                xyz_lines = [
                    str(len(ag.atoms)),
                    "Orientation CV reference positions - generated by GateWizard",
                ]
                for atom in ag.atoms:
                    pos = atom.position  # MDAnalysis positions are in Angstrom
                    xyz_lines.append(f"CA  {pos[0]:.6f}  {pos[1]:.6f}  {pos[2]:.6f}")
                xyz_path.write_text("\n".join(xyz_lines) + "\n")
                if ref_base_dir is not None:
                    resolved_ref_file = str(xyz_path.relative_to(ref_base_dir))
                else:
                    resolved_ref_file = xyz_path.name
                self.logger.info(f"  Rotation reference XYZ file: {xyz_path.name}")

            content = _build_com_colvars_config(
                atom_numbers=atom_nums,
                x0=x0,
                y0=y0,
                z0=z0,
                com_k=com_restraint_k,
                add_rotation=add_rotation_restraint,
                rot_k=rotation_restraint_k,
                ag=ag,
                engine="gromacs",
                ref_positions_file=resolved_ref_file,
                rotation_ref_positions_mode=rotation_ref_positions_mode,
                ref_positions_col=ref_positions_col,
                ref_positions_col_value=ref_positions_col_value,
            )

            output_file.write_text(content)
            self.logger.info(
                f"  COM colvars config ({len(ag)} atoms, centroid "
                f"[{x0:.2f},{y0:.2f},{z0:.2f}] Å): {output_file.name}"
            )
            return output_file

        except ImportError:
            self.logger.warning(
                "MDAnalysis not available; COM colvars file not generated. "
            )
            return None
        except Exception as exc:
            self.logger.error(f"COM colvars generation failed: {exc}")
            return None

    # ------------------------------------------------------------------
    # Full setup orchestration
    # ------------------------------------------------------------------

    def setup_gromacs_equilibration(
        self,
        system_files: Optional[Dict[str, str]] = None,
        stage_params_list: Optional[List[Dict[str, Any]]] = None,
        output_name: str = "equilibration",
        scheme_type: Optional[str] = None,
        selections: Optional[Dict[str, str]] = None,
        gmx_executable: str = "gmx",
        gmxrc_path: Optional[str] = None,
        add_com_restraint: bool = False,
        com_restraint_k: float = 10.0,
        add_rotation_restraint: bool = True,
        rotation_restraint_k: float = 2000.0,
        com_selection: str = "name CA",
        rotation_ref_positions_mode: str = "auto",
        ref_positions_file: Optional[str] = None,
        ref_positions_col: Optional[str] = None,
        ref_positions_col_value: Optional[float] = None,
    ) -> Dict[str, Any]:
        """Complete GROMACS equilibration setup.

        Converts AMBER files to GROMACS format (if needed), generates index
        and position-restraint files, writes MDP files for all stages, and
        produces a bash run script.

        Args:
            system_files: Dict with GROMACS keys (``gro``, ``top``, ``ndx``)
                *or* AMBER keys (``prmtop``, ``inpcrd``, ``pdb``,
                ``bilayer_pdb``).  Auto-detected from *working_dir* if None.
            stage_params_list: List of stage dicts or
                :class:`EquilibrationStage` objects.  Keys:

                - ``name`` (str)
                - ``ensemble`` (str): NVT | NPT | NPAT | NPgT
                - ``time_ns`` (float)
                - ``timestep`` (float, fs)
                - ``temperature`` (float, K)
                - ``constraints`` (dict, kcal/mol/Å²):
                  ``protein_backbone``, ``protein_sidechain``,
                  ``lipid_head``, ``lipid_tail``, ``water``, ``ions``,
                  ``other``, plus any custom keys with MDAnalysis selections

                Defaults to the 7-stage CHARMM-GUI protocol when *None*.
            output_name: Subdirectory name under *working_dir*.
            scheme_type: Override ensemble; auto-detected from first stage.
            selections: MDAnalysis selection overrides for water/ions/other
                and custom constraint keys (same dict used by NAMD/OpenMM).
            gmx_executable: GROMACS executable (default: ``"gmx"``).  Can
                include full path, e.g. ``"/usr/local/gromacs/bin/gmx"``.
            gmxrc_path: Optional ``GMXRC`` to source in the run script (system
                installs). Leave ``None`` for conda / self-contained binaries.
            add_com_restraint: When True, generate a Colvars COM restraint
                file (requires GROMACS 2021+ with Colvars support).
            com_restraint_k: COM translation force constant in kcal/mol/Å².
            add_rotation_restraint: Add orientation restraint in colvars.
            rotation_restraint_k: Rotation force constant in kcal/mol/Å².
            rotation_ref_positions_mode: How to encode orientation reference
                coordinates: ``"auto"`` (default, generates a ``.xyz`` file),
                ``"refPositions"`` (inline), or ``"refPositionsFile"``
                (supply *ref_positions_file*).
            ref_positions_file: Path to an existing **XYZ** file used as the
                orientation reference when mode is ``"refPositionsFile"``.
                If not given and mode is ``"auto"``, a ``.xyz`` file is
                generated automatically.  **Note:** GROMACS Colvars only
                supports XYZ format — PDB files will cause a runtime error.

        Returns:
            Dict with keys:
            ``gromacs_dir``, ``mdp_files``, ``run_script``, ``system_files``,
            ``posres_files``, ``gro``, ``top``, ``ndx``.

        Example::

            >>> from pathlib import Path
            >>> from gatewizard.tools.equilibration import GROMACSEquilibrationManager
            >>> manager = GROMACSEquilibrationManager(Path("/work/membrane"))
            >>> stages = GROMACSEquilibrationManager.get_default_stage_params("NPT",
            ...                                                                include_production=True)
            >>> stages[-1].time_ns = 5.0   # quick test run
            >>> result = manager.setup_gromacs_equilibration(stage_params_list=stages)
            >>> # cd result["gromacs_dir"] && bash run_equilibration.sh
        """
        self.logger.info("=== Setting up GROMACS equilibration ===")

        if system_files is None:
            system_files = self.find_system_files()
            if system_files is None:
                raise ValueError(
                    "Could not auto-detect system files in working directory"
                )

        if not stage_params_list:
            _scheme = scheme_type or "NPT"
            stage_params_list = self.get_default_stage_params(_scheme)
            self.logger.info(
                f"No stages provided — using default {_scheme} protocol "
                f"({len(stage_params_list)} stages)"
            )

        # Normalise EquilibrationStage → dict
        stage_params_list = [
            s.to_dict() if isinstance(s, EquilibrationStage) else s
            for s in stage_params_list
        ]

        if scheme_type is None:
            scheme_type = stage_params_list[0].get("ensemble", "NPT")
            self.logger.info(f"Auto-detected scheme_type: {scheme_type}")

        if scheme_type not in self.SCHEME_MAPPING:
            raise ValueError(
                f"Unknown scheme_type '{scheme_type}'. "
                f"Must be one of {list(self.SCHEME_MAPPING.keys())}"
            )

        # GROMACS uses a dedicated step0 minimization stage before equilibration
        # (unlike NAMD/OpenMM where minimization is embedded in Equilibration 1).
        # If the first stage is not explicitly named "Minimization" but carries
        # minimize_steps > 0, auto-split it so GROMACS always emits the full
        # step0 + step1…step6 + step7 sequence.
        first_stage = stage_params_list[0]
        if (
            first_stage.get("name", "").lower()
            not in ("minimization", "energy minimization", "energy_minimization")
            and int(first_stage.get("minimize_steps", 0)) > 0
        ):
            mini_stage = dict(first_stage)
            mini_stage["name"] = "Minimization"
            mini_stage["time_ns"] = 0.0

            eq1_stage = dict(first_stage)
            eq1_stage["minimize_steps"] = 0

            stage_params_list = [mini_stage, eq1_stage] + stage_params_list[1:]
            self.logger.info(
                "  Auto-inserted GROMACS Minimization stage (split from "
                f"'{first_stage.get('name', 'first stage')}')"
            )

        # output_name may be an absolute Path (mirrors OpenMM behaviour)
        gromacs_dir = (
            Path(output_name)
            if Path(str(output_name)).is_absolute()
            else self.working_dir / output_name
        )
        gromacs_dir.mkdir(parents=True, exist_ok=True)
        restraints_dir = gromacs_dir / "restraints"
        restraints_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"Output directory: {gromacs_dir}")

        # --- Obtain GRO + TOP ---
        gro_path: Optional[Path] = None
        top_path: Optional[Path] = None
        ndx_path: Optional[Path] = None

        if "gro" in system_files:
            # Native GROMACS files — copy them
            gro_src = Path(system_files["gro"])
            top_src = Path(system_files["top"])
            gro_path = gromacs_dir / gro_src.name
            top_path = gromacs_dir / top_src.name
            shutil.copy2(gro_src, gro_path)
            shutil.copy2(top_src, top_path)
            if "ndx" in system_files:
                ndx_src = Path(system_files["ndx"])
                ndx_path = gromacs_dir / ndx_src.name
                shutil.copy2(ndx_src, ndx_path)
            self.logger.info(f"  Copied GROMACS files")
        else:
            # AMBER → convert with ParmEd
            prmtop = system_files.get("prmtop")
            inpcrd = system_files.get("inpcrd")
            bilayer_pdb_src = system_files.get("bilayer_pdb")
            if not prmtop or not inpcrd:
                raise ValueError(
                    "system_files must contain 'gro'+'top' or 'prmtop'+'inpcrd'"
                )
            conv = self.convert_from_amber(
                prmtop=Path(prmtop),
                inpcrd=Path(inpcrd),
                output_dir=gromacs_dir,
                bilayer_pdb=Path(bilayer_pdb_src) if bilayer_pdb_src else None,
            )
            gro_path = conv["gro"]
            top_path = conv["top"]
            # Copy PDB for reference / COM restraints
            pdb_src = system_files.get("pdb")
            if pdb_src and Path(pdb_src).exists():
                shutil.copy2(pdb_src, gromacs_dir / Path(pdb_src).name)

        # --- Generate index.ndx ---
        ndx_path = gromacs_dir / "index.ndx"
        self.generate_index_ndx(gro_path, ndx_path)
        # Add CHARMM-GUI-style aliases (SOLU/MEMB/SOLV/SOLU_MEMB) expected
        # by the MDP templates' tc_grps and comm_grps fields.
        if ndx_path.exists():
            self._append_charmm_gui_groups(ndx_path)
        ndx_name = ndx_path.name if ndx_path.exists() else None

        # --- Generate posres .itp files ---
        posres_files: Dict[str, Any] = {}
        any_constraints = any(
            any(float(v) > 0 for v in s.get("constraints", {}).values())
            for s in stage_params_list
        )
        if any_constraints and gro_path.exists():
            constraints_max: Dict[str, float] = {}
            for s in stage_params_list:
                for k, v in s.get("constraints", {}).items():
                    fv = float(v)
                    if fv > constraints_max.get(k, 0.0):
                        constraints_max[k] = fv
            pdb_for_posres = None
            for cand in (
                gromacs_dir / "system.pdb",
                Path(system_files.get("pdb", "")) if system_files.get("pdb") else None,
            ):
                if cand is not None and cand.exists():
                    pdb_for_posres = cand
                    break
            posres_files = self.generate_posres_itp(
                gro_path=gro_path,
                index_path=ndx_path,
                output_dir=gromacs_dir,
                top_path=top_path,
                pdb_path=pdb_for_posres,
                selections=selections,
                constraints_max=constraints_max,
            )
            # Modify topology to include posres
            top_path = self._add_posres_to_topology(
                topol_path=top_path,
                posres_files=posres_files,
            )
            self.logger.info(f"  Using topology with posres: {top_path.name}")

        # --- Generate MDP files ---
        mdp_files: List[Path] = []
        n_eq_stages = 0
        used_posres_macros = self._topology_posres_macros(posres_files)

        for stage_index, stage_params in enumerate(stage_params_list):
            stage_name = stage_params.get("name", f"Stage {stage_index}")

            # Determine stage_index_key:
            # First stage with minimize_steps > 0 → minimization (idx 0)
            # Then equilibration stages (1..6), then production (7)
            if stage_index == 0 and int(stage_params.get("minimize_steps", 0)) > 0:
                key_idx = 0
            elif stage_params.get("name", "").lower() == "production":
                key_idx = 7
            else:
                key_idx = min(stage_index, 6) if stage_index > 0 else 1
                if key_idx == 0:
                    key_idx = 1

            self.logger.debug(
                f"  Stage {stage_index} → key_idx={key_idx} ({stage_name})"
            )
            content = self.generate_mdp_file(
                stage_name=stage_name,
                stage_params=stage_params,
                stage_index=key_idx,
                scheme_type=scheme_type,
                used_posres_macros=used_posres_macros if posres_files else None,
            )
            mdp_filename = self._get_mdp_filename(key_idx)
            mdp_path = gromacs_dir / mdp_filename
            mdp_path.write_text(content)
            mdp_files.append(mdp_path)
            if 1 <= key_idx <= 6:
                n_eq_stages = max(n_eq_stages, key_idx)
            self.logger.info(f"  Written: {mdp_filename}")

        if n_eq_stages == 0:
            n_eq_stages = 6  # default

        # --- COM colvars restraint ---
        com_colvars_path: Optional[Path] = None
        if add_com_restraint:
            pdb_for_com = gromacs_dir / "system.pdb"
            if not pdb_for_com.exists():
                # try to find any .pdb
                pdbs = list(gromacs_dir.glob("*.pdb"))
                pdb_for_com = pdbs[0] if pdbs else None
            if pdb_for_com and pdb_for_com.exists():
                com_colvars_relpath = Path("restraints") / "com_restraint.dat"
                com_colvars_path = self.generate_com_colvars_config(
                    pdb_path=pdb_for_com,
                    output_file=gromacs_dir / com_colvars_relpath,
                    com_restraint_k=com_restraint_k,
                    add_rotation_restraint=add_rotation_restraint,
                    rotation_restraint_k=rotation_restraint_k,
                    selection=com_selection,
                    rotation_ref_positions_mode=rotation_ref_positions_mode,
                    ref_positions_file=ref_positions_file,
                    ref_positions_col=ref_positions_col,
                    ref_positions_col_value=ref_positions_col_value,
                    ref_base_dir=gromacs_dir,
                )
                if com_colvars_path:
                    activation_block = _build_com_colvars_activation_block(
                        "gromacs", str(com_colvars_relpath)
                    )
                    for mdp_path in mdp_files:
                        mdp_text = mdp_path.read_text()
                        if "colvars-active" not in mdp_text:
                            mdp_path.write_text(mdp_text.rstrip() + activation_block)
                    self.logger.info(
                        "  COM colvars file generated and activated in GROMACS MDP files."
                    )
            else:
                self.logger.warning(
                    "No PDB found in output dir; COM colvars file not generated."
                )

        # --- Run script ---
        compute = resolve_compute_resources_from_stages(stage_params_list)
        run_script = self.generate_run_script(
            gromacs_dir=gromacs_dir,
            gro_name=gro_path.name,
            top_name=top_path.name,
            ndx_name=ndx_name,
            n_stages=n_eq_stages,
            gmx_executable=gmx_executable,
            gmxrc_path=gmxrc_path,
            cpu_cores=compute["cpu_cores"],
            use_gpu=compute["use_gpu"],
            gpu_id=compute["gpu_id"],
            num_gpus=compute["num_gpus"] or 1,
        )

        self.logger.info("=== GROMACS equilibration setup complete ===")
        return {
            "gromacs_dir": gromacs_dir,
            "mdp_files": mdp_files,
            "run_script": run_script,
            "system_files": system_files,
            "posres_files": posres_files,
            "gro": gro_path,
            "top": top_path,
            "ndx": ndx_path,
            "com_colvars": com_colvars_path,
        }

    @staticmethod
    def get_default_selections(pdb_path: str) -> Dict[str, str]:
        """Delegate to :meth:`NAMDEquilibrationManager.get_default_selections`."""
        return NAMDEquilibrationManager.get_default_selections(pdb_path)


class AmberEquilibrationManager:
    """Manager for AMBER equilibration simulations (placeholder)."""

    def __init__(self, working_dir: Path):
        self.working_dir = Path(working_dir)
        self.logger = get_logger(self.__class__.__name__)

    def generate_input_file(self, stage_name: str, stage_params: Dict[str, Any]) -> str:
        """Generate AMBER input file (placeholder)."""
        self.logger.info("AMBER equilibration not yet implemented")
        return ""


class EquilibrationAnalyzer:
    """Analyzer for equilibration simulation results."""

    def __init__(self, working_dir: Path):
        self.working_dir = Path(working_dir)
        self.logger = get_logger(self.__class__.__name__)

    def analyze_energy_convergence(self, log_files: List[Path]) -> Dict[str, Any]:
        """
        Analyze energy convergence from simulation log files.

        Args:
            log_files: List of simulation log files

        Returns:
            Dictionary with convergence analysis results
        """

        results = {
            "converged": False,
            "total_energy": [],
            "temperature": [],
            "pressure": [],
            "volume": [],
        }

        # Placeholder implementation
        self.logger.info("Energy convergence analysis not yet implemented")

        return results

    def generate_plots(
        self, analysis_results: Dict[str, Any], output_dir: Path
    ) -> List[Path]:
        """
        Generate plots for equilibration analysis.

        Args:
            analysis_results: Results from analysis
            output_dir: Directory to save plots

        Returns:
            List of generated plot files
        """

        plots = []

        # Placeholder for plot generation
        self.logger.info("Plot generation not yet implemented")

        return plots
