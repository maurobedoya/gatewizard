"""Write one trajectory timestep as a single-model PDB."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def write_current_frame_pdb(universe: Any, dest: str | Path) -> str:
    """Write the universe's current timestep as one PDB.

    A multi-model PDB or a DCD/XTC stays on disk unchanged. Only the frame
    the universe is positioned on is written, including its periodic box
    when the trajectory has one.
    """
    path = Path(dest).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".writing")
    try:
        universe.atoms.write(str(partial), file_format="PDB")
        partial.replace(path)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    return str(path)
