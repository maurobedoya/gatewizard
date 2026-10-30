"""Read Maestro ``.mae`` / ``.maegz`` connection-table properties.

The desktop Visualize page uses this module for Glide, MM/GBSA, IFD / Prime,
and other CT-level scores. Geometry materialization (cache PDBs) stays in the
GUI backend; this API returns titles, atom counts, and typed properties.
"""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path
from typing import Any, Optional, Sequence, Union

# Maestro sometimes embeds a full PDB (or similar) in a string property.
MAESTRO_PROP_STRING_MAX = 500

_MAESTRO_PROP_TYPES = {"r": "real", "i": "integer", "b": "boolean", "s": "string"}

# Same preference order as the Visualize entry picker and structure-row headline.
HEADLINE_PROPERTY_NAMES = (
    "r_psp_MMGBSA_dG_Bind",
    "r_i_docking_score",
    "r_i_glide_gscore",
    "r_psp_Prime_Energy",
)

PathLike = Union[str, Path]


class MaestroError(RuntimeError):
    """Raised when a Maestro file cannot be read for properties."""


def is_maestro_path(path: PathLike) -> bool:
    """True for ``.mae``, ``.maegz``, or ``.mae.gz``."""
    p = Path(path)
    suf = p.suffix.lower()
    if suf in {".mae", ".maegz"}:
        return True
    return p.name.lower().endswith(".mae.gz")


def _open_mae_text_path(path: Path) -> tuple[Path, bool]:
    """Return ``(path_to_uncompressed_mae, needs_cleanup)``."""
    name = path.name.lower()
    if name.endswith(".maegz") or name.endswith(".mae.gz"):
        tmp = tempfile.NamedTemporaryFile(suffix=".mae", delete=False)
        tmp_path = Path(tmp.name)
        tmp.close()
        with gzip.open(path, "rb") as src, tmp_path.open("wb") as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
        return tmp_path, True
    return path, False


def maestro_ct_properties(mol: Any) -> list[dict[str, Any]]:
    """CT-level properties from an RDKit molecule.

    Keeps ``r_`` / ``i_`` / ``b_`` / ``s_`` values. Drops RDKit private names
    (leading ``_``) and strings longer than :data:`MAESTRO_PROP_STRING_MAX`.

    Each item is ``{ "name", "type", "value" }`` with
    ``type`` in ``real`` / ``integer`` / ``boolean`` / ``string``.
    """
    if mol is None:
        return []
    try:
        raw = mol.GetPropsAsDict(includePrivate=False, includeComputed=False)
    except Exception:
        raw = {}
        try:
            for name in mol.GetPropNames():
                raw[name] = mol.GetProp(name)
        except Exception:
            return []

    out: list[dict[str, Any]] = []
    for name, value in raw.items():
        if not isinstance(name, str) or len(name) < 3:
            continue
        if name.startswith("_"):
            continue
        kind = _MAESTRO_PROP_TYPES.get(name[0])
        if kind is None or name[1] != "_":
            continue
        if kind == "string":
            text = str(value)
            if len(text) > MAESTRO_PROP_STRING_MAX:
                continue
            value = text
        elif kind == "real":
            try:
                value = float(value)
            except (TypeError, ValueError):
                continue
            if not (value == value):  # NaN
                continue
        elif kind == "integer":
            try:
                value = int(value)
            except (TypeError, ValueError):
                continue
        else:
            value = bool(value)
        out.append({"name": name, "type": kind, "value": value})
    out.sort(key=lambda item: item["name"])
    return out


def pick_headline_property(
    properties: Sequence[dict[str, Any]] | None,
) -> Optional[dict[str, Any]]:
    """First preferred real score, else the first real property, else ``None``."""
    list_props = list(properties or [])
    by_name = {p.get("name"): p for p in list_props if isinstance(p, dict)}
    for name in HEADLINE_PROPERTY_NAMES:
        hit = by_name.get(name)
        if hit and hit.get("type") == "real":
            return hit
    for prop in list_props:
        if isinstance(prop, dict) and prop.get("type") == "real":
            return prop
    return None


def _entry_label(source: Path, index: int, title: str | None) -> str:
    base = source.name
    if title and str(title).strip():
        return f"{str(title).strip()} · {base}"
    return f"CT {index + 1} · {base}"


def list_maestro_entries(path: PathLike) -> dict[str, Any]:
    """List connection tables and their CT-level properties.

    Requires RDKit ``MaeMolSupplier``. Raises :class:`MaestroError` when the
    file cannot be read.

    Returns::

        {
          "sourcePath": str,
          "entries": [
            {
              "index": int,          # supplier index (0-based)
              "title": str | None,
              "label": str,
              "atomCount": int,
              "kind": "maestro_ct",
              "properties": [ {"name", "type", "value"}, ... ],
            },
            ...
          ],
        }
    """
    p = Path(path).expanduser().resolve()
    if not p.is_file():
        raise FileNotFoundError(f"File not found: {p}")
    if not is_maestro_path(p):
        raise MaestroError(f"Not a Maestro file: {p.name}")

    try:
        from rdkit import Chem
    except ImportError as exc:
        raise MaestroError(
            "RDKit is required to read Maestro properties. "
            "Install gatewizard with its scientific dependencies."
        ) from exc
    if not hasattr(Chem, "MaeMolSupplier"):
        raise MaestroError("This RDKit build has no MaeMolSupplier")

    mae_path, cleanup = _open_mae_text_path(p)
    entries: list[dict[str, Any]] = []
    try:
        supplier = Chem.MaeMolSupplier(str(mae_path), sanitize=False, removeHs=False)
        for i, mol in enumerate(supplier):
            if mol is None:
                continue
            title = None
            try:
                if mol.HasProp("_Name"):
                    title = str(mol.GetProp("_Name"))
            except Exception:
                title = None
            n_atoms = int(mol.GetNumAtoms())
            entries.append(
                {
                    "index": i,
                    "title": title,
                    "label": _entry_label(p, i, title),
                    "atomCount": n_atoms,
                    "kind": "maestro_ct",
                    "properties": maestro_ct_properties(mol),
                }
            )
    except MaestroError:
        raise
    except Exception as exc:
        raise MaestroError(f"Could not read Maestro file {p.name}: {exc}") from exc
    finally:
        if cleanup:
            try:
                mae_path.unlink(missing_ok=True)
            except OSError:
                pass

    if not entries:
        raise MaestroError(f"No connection tables found in {p.name}")
    return {"sourcePath": str(p), "entries": entries}


def read_maestro_properties(
    path: PathLike,
    indices: Sequence[int] | None = None,
) -> list[dict[str, Any]]:
    """Return properties for selected CTs (all CTs when ``indices`` is omitted).

    Each item is ``{ "index", "title", "label", "atomCount", "properties" }``.
    Schrödinger energy properties are in kcal/mol.
    """
    info = list_maestro_entries(path)
    entries = info["entries"]
    if indices is None:
        chosen = entries
    else:
        want = {int(i) for i in indices}
        chosen = [e for e in entries if e["index"] in want]
        missing = want - {e["index"] for e in chosen}
        if missing:
            raise ValueError(f"CT index out of range: {sorted(missing)}")
    return [
        {
            "index": e["index"],
            "title": e["title"],
            "label": e["label"],
            "atomCount": e["atomCount"],
            "properties": e["properties"],
        }
        for e in chosen
    ]
