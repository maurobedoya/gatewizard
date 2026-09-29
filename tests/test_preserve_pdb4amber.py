"""Optional pdb4amber integration for preserve_residue_numbers."""

import shutil
from pathlib import Path

import pytest

from gatewizard.core.preparation import PreparationManager
from gatewizard.utils.residue_restore import snapshot_pdb_residues


def _atom(serial: int, name: str, resname: str, chain: str, resid: int, x: float) -> str:
    return (
        f"ATOM  {serial:5d}  {name:<3s} {resname:>3s} {chain}{resid:4d}    "
        f"{x:8.3f}{x:8.3f}{x:8.3f}  1.00  0.00           {name[0]}\n"
    )


def _pdb4amber_available() -> bool:
    return shutil.which("pdb4amber") is not None


@pytest.mark.skipif(not _pdb4amber_available(), reason="pdb4amber not installed")
def test_pdb4amber_preserve_keeps_gap_and_hetero(tmp_path: Path):
    src = tmp_path / "gapped.pdb"
    src.write_text(
        _atom(1, "N", "ALA", "A", 10, 1.0)
        + _atom(2, "CA", "ALA", "A", 10, 2.0)
        + _atom(3, "C", "ALA", "A", 10, 3.0)
        + _atom(4, "O", "ALA", "A", 10, 4.0)
        + _atom(5, "N", "GLY", "A", 20, 5.0)
        + _atom(6, "CA", "GLY", "A", 20, 6.0)
        + _atom(7, "C", "GLY", "A", 20, 7.0)
        + _atom(8, "O", "GLY", "A", 20, 8.0)
        + "HETATM    9  C1  LIG A 900       9.000   9.000   9.000  1.00  0.00           C\n"
        + "END\n",
        encoding="utf-8",
    )
    out = tmp_path / "prepared.pdb"
    manager = PreparationManager()
    result = manager.run_pdb4amber_with_cap_fix(
        input_pdb=str(src),
        output_pdb=str(out),
        fix_caps=False,
        preserve_residue_numbers=True,
        original_pdb=str(src),
    )
    assert result["success"]
    assert result["residue_numbers_preserved"] is True
    snaps = snapshot_pdb_residues(out)
    resids = {(chain, resid) for chain, resid, _icode, resname in snaps}
    assert ("A", 10) in resids
    assert ("A", 20) in resids
    assert ("A", 900) in resids
    hetero = [
        line
        for line in out.read_text(encoding="utf-8").splitlines()
        if line.startswith("HETATM") and "LIG" in line
    ]
    assert hetero
