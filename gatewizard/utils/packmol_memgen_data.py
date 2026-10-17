"""Locate the PACKMOL-Memgen data files tleap needs for Lipid_ext."""

from pathlib import Path

import packmol_memgen


def packmol_memgen_data_dir() -> Path:
    """Directory that holds ``leaprc.lipid_ext``, ``frcmod.lipid_ext``, and ``lipid_ext.lib``.

    ``leaprc.lipid_ext`` loads the frcmod and lib by bare filename, so tleap
    only finds them after ``addPath`` to this directory.
    """
    data = Path(packmol_memgen.__file__).resolve().parent / "data"
    leap = data / "leaprc.lipid_ext"
    if not leap.is_file() or not (data / "frcmod.lipid_ext").is_file() or not (data / "lipid_ext.lib").is_file():
        raise FileNotFoundError(
            "PACKMOL-Memgen Lipid_ext parameters were not found next to the installed "
            f"packmol_memgen package ({data})."
        )
    return data
