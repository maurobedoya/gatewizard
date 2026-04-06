# gatewizard/core/mempro.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Membrane protein orientation module using MemPrO.

This module provides a programmatic API for orienting membrane proteins
using the MemPrO tool (https://github.com/pstansfeld/MemPrO).
MemPrO positions proteins correctly in the membrane prior to system
building with packmol-memgen.
"""

import os
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)


class MemProError(Exception):
    """Error raised by MemPrO operations."""


class OrientationResult:
    """
    A single ranked orientation from a MemPrO run.

    Attributes:
        rank: Integer rank (1 = best).
        relative_potential: Relative potential score.
        hits_pct: Percentage of hits.
        rerank_potential: Re-ranked potential.
        rerank_depth: Re-rank minima depth.
        rerank_value: Re-rank value.
        pdb_path: Path to the oriented PDB file for this rank.
    """

    def __init__(
        self,
        rank: int,
        relative_potential: float,
        hits_pct: float,
        rerank_potential: float,
        rerank_depth: float,
        rerank_value: float,
        pdb_path: str,
    ):
        self.rank = rank
        self.relative_potential = relative_potential
        self.hits_pct = hits_pct
        self.rerank_potential = rerank_potential
        self.rerank_depth = rerank_depth
        self.rerank_value = rerank_value
        self.pdb_path = pdb_path

    def __repr__(self) -> str:
        return (
            f"OrientationResult(rank={self.rank}, "
            f"potential={self.relative_potential:.3f}, "
            f"hits={self.hits_pct:.1f}%)"
        )


class MemPrO:
    """
    Main class for orienting membrane proteins using MemPrO.

    This class wraps the ``mempro`` CLI to orient a protein structure
    in a membrane, parse the ranked results, and provide easy access
    to the oriented PDB files.

    Example::

        from gatewizard.core.mempro import MemPrO

        mp = MemPrO()
        results = mp.run("protein.pdb")
        print(results[0])           # best orientation
        print(results[0].pdb_path)  # path to oriented PDB
    """

    def __init__(self):
        """Initialize MemPrO wrapper."""
        self._mempro_cmd = shutil.which("mempro")

    @staticmethod
    def is_available() -> bool:
        """Return True if the ``mempro`` executable is on PATH."""
        return shutil.which("mempro") is not None

    def run(
        self,
        pdb_file: str,
        output_dir: Optional[str] = None,
        n_cpus: Optional[int] = None,
        n_iters: int = 150,
        grid_size: int = 36,
        dual_membrane: bool = False,
        peripheral: bool = False,
        use_weights: bool = False,
        flip: bool = False,
        membrane_thickness: Optional[float] = None,
        extra_args: Optional[List[str]] = None,
    ) -> List[OrientationResult]:
        """
        Run MemPrO orientation on a PDB file.

        Args:
            pdb_file: Path to the input PDB file.
            output_dir: Name of the output directory (default: ``Orient``).
            n_cpus: Number of CPU cores to use (default: all).
            n_iters: Number of minimisation iterations (default: 150).
            grid_size: Number of starting configurations (default: 36).
            dual_membrane: Enable dual membrane orientation.
            peripheral: Enable peripheral protein orientation.
            use_weights: Use B-factors to weight orientation.
            flip: Flip protein in the Z-axis after orientation.
            membrane_thickness: Initial membrane thickness in Å (default: 28).
            extra_args: Additional CLI arguments passed verbatim.

        Returns:
            List of :class:`OrientationResult` sorted by rank.

        Raises:
            MemProError: If ``mempro`` is not installed or execution fails.
            FileNotFoundError: If the input PDB file does not exist.
        """
        pdb_path = Path(pdb_file).resolve()
        if not pdb_path.is_file():
            raise FileNotFoundError(f"PDB file not found: {pdb_path}")

        if not self.is_available():
            raise MemProError(
                "mempro is not installed. "
                "Install with: pip install git+https://github.com/pstansfeld/MemPrO.git"
            )

        cmd = [self._mempro_cmd, "-f", str(pdb_path)]

        if output_dir:
            cmd.extend(["-o", str(output_dir)])
        if n_cpus is not None:
            cmd.extend(["-nc", str(n_cpus)])
        if n_iters != 150:
            cmd.extend(["-ni", str(n_iters)])
        if grid_size != 36:
            cmd.extend(["-ng", str(grid_size)])
        if dual_membrane:
            cmd.append("-dm")
        if peripheral:
            cmd.append("-pr")
        if use_weights:
            cmd.append("-w")
        if flip:
            cmd.append("-flip")
        if membrane_thickness is not None:
            cmd.extend(["-mt", str(membrane_thickness)])
        if extra_args:
            cmd.extend(extra_args)

        # Run from the directory containing the PDB so output lands nearby
        work_dir = str(pdb_path.parent)
        logger.info(f"Running MemPrO: {' '.join(cmd)}")

        try:
            result = subprocess.run(
                cmd,
                cwd=work_dir,
                capture_output=True,
                text=True,
                timeout=3600,
            )
        except subprocess.TimeoutExpired:
            raise MemProError("MemPrO timed out after 1 hour.")
        except Exception as e:
            raise MemProError(f"Failed to run MemPrO: {e}")

        if result.returncode != 0:
            logger.error(f"MemPrO stderr: {result.stderr}")
            raise MemProError(
                f"MemPrO exited with code {result.returncode}.\n{result.stderr}"
            )

        # Determine the output folder
        orient_name = output_dir if output_dir else "Orient"
        orient_path = Path(work_dir) / orient_name
        if not orient_path.is_dir():
            raise MemProError(f"Expected output directory not found: {orient_path}")

        return self.parse_results(str(orient_path))

    @staticmethod
    def parse_results(orient_dir: str) -> List[OrientationResult]:
        """
        Parse MemPrO results from an existing Orient directory.

        Args:
            orient_dir: Path to the Orient output directory.

        Returns:
            List of :class:`OrientationResult` sorted by rank.

        Raises:
            MemProError: If the orientation file cannot be parsed.
        """
        orient_path = Path(orient_dir)
        txt_file = orient_path / "orientation.txt"

        if not txt_file.is_file():
            raise MemProError(f"orientation.txt not found in {orient_path}")

        results: List[OrientationResult] = []

        with open(txt_file, "r") as fh:
            for line in fh:
                line = line.strip()
                # Skip header, comments, and blank lines
                if not line or line.startswith("Generated") or line.startswith("Rank"):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    rank = int(parts[0])
                    rel_pot = float(parts[1])
                    hits = float(parts[2])
                    rr_pot = float(parts[3])
                    rr_depth = float(parts[4])
                    rr_val = float(parts[5])
                except (ValueError, IndexError):
                    continue

                rank_dir = orient_path / f"Rank_{rank}"
                pdb_file = rank_dir / f"oriented_rank_{rank}.pdb"
                pdb_str = str(pdb_file) if pdb_file.is_file() else ""

                results.append(
                    OrientationResult(
                        rank=rank,
                        relative_potential=rel_pot,
                        hits_pct=hits,
                        rerank_potential=rr_pot,
                        rerank_depth=rr_depth,
                        rerank_value=rr_val,
                        pdb_path=pdb_str,
                    )
                )

        results.sort(key=lambda r: r.rank)
        logger.info(f"Parsed {len(results)} orientation results from {orient_path}")
        return results

    @staticmethod
    def get_oriented_pdb(orient_dir: str, rank: int = 1) -> str:
        """
        Get the path to an oriented PDB file for a specific rank.

        Args:
            orient_dir: Path to the Orient output directory.
            rank: Desired rank number (default: 1).

        Returns:
            Absolute path to the oriented PDB file.

        Raises:
            FileNotFoundError: If the PDB file for the given rank does not exist.
        """
        pdb_path = Path(orient_dir) / f"Rank_{rank}" / f"oriented_rank_{rank}.pdb"
        if not pdb_path.is_file():
            raise FileNotFoundError(
                f"Oriented PDB not found for rank {rank}: {pdb_path}"
            )
        return str(pdb_path.resolve())

    def build_command(
        self,
        pdb_file: str,
        output_dir: Optional[str] = None,
        n_cpus: Optional[int] = None,
        n_iters: int = 150,
        grid_size: int = 36,
        dual_membrane: bool = False,
        peripheral: bool = False,
        use_weights: bool = False,
        flip: bool = False,
        membrane_thickness: Optional[float] = None,
        extra_args: Optional[List[str]] = None,
    ) -> List[str]:
        """
        Build the MemPrO command without executing it.

        Accepts the same parameters as :meth:`run`.

        Returns:
            List of command-line tokens.
        """
        if not self._mempro_cmd:
            raise MemProError("mempro executable not found on PATH.")

        cmd = [self._mempro_cmd, "-f", str(Path(pdb_file).resolve())]
        if output_dir:
            cmd.extend(["-o", str(output_dir)])
        if n_cpus is not None:
            cmd.extend(["-nc", str(n_cpus)])
        if n_iters != 150:
            cmd.extend(["-ni", str(n_iters)])
        if grid_size != 36:
            cmd.extend(["-ng", str(grid_size)])
        if dual_membrane:
            cmd.append("-dm")
        if peripheral:
            cmd.append("-pr")
        if use_weights:
            cmd.append("-w")
        if flip:
            cmd.append("-flip")
        if membrane_thickness is not None:
            cmd.extend(["-mt", str(membrane_thickness)])
        if extra_args:
            cmd.extend(extra_args)
        return cmd
