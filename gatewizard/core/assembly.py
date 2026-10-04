"""Build a biological assembly from operators already stored in a PDB or mmCIF.

The coordinates in ``ATOM`` records are the asymmetric unit. The oligomer or
capsid is the set of ``REMARK 350`` / ``BIOMT`` operators (PDB) or
``_pdbx_struct_assembly`` operators (mmCIF). Crystal neighbors from ``CRYST1``
are a different operation and are not built here.
"""

from __future__ import annotations

import math
import os
import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import gemmi
import MDAnalysis as mda
import numpy as np

# Same comfort line the viewer uses for heavy representations. A capsid past
# this asks for confirmation before it is written.
LARGE_ASSEMBLY_ATOM_COUNT = 80_000

ASSEMBLY_MMCIF_TAG = "_gatewizard.biological_assembly"
_NCS_ID = "ncs"
_WATER_NAMES = frozenset({"HOH", "WAT", "DOD"})
# Atoms closer than this are the same site copied by a symmetry operator.
# A real contact is several angstroms further out.
_STACK_CUTOFF = 0.5


class AssemblyError(ValueError):
    """The file has no biological assembly that can be built."""


def _might_be_mmcif(path: Path) -> bool:
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return True
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for _ in range(20):
                line = handle.readline()
                if not line:
                    return False
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                return stripped.lower().startswith("data_")
    except OSError:
        return False
    return False


def is_biological_assembly_mmcif(path: str | Path) -> bool:
    """True when this file was written by :func:`build_assembly`.

    A Save that copied the mmCIF and named it ``.pdb`` still counts: the
    marker is in the text, and the suffix is not trusted on its own.
    """
    file_path = Path(path)
    if not _might_be_mmcif(file_path):
        return False
    try:
        with file_path.open("r", encoding="utf-8", errors="replace") as handle:
            overlap = ""
            while True:
                chunk = handle.read(1 << 20)
                if not chunk:
                    return ASSEMBLY_MMCIF_TAG in overlap
                blob = overlap + chunk
                if ASSEMBLY_MMCIF_TAG in blob:
                    return True
                overlap = blob[-len(ASSEMBLY_MMCIF_TAG) :]
    except OSError:
        return False


def list_assemblies(path: str | Path) -> dict[str, Any]:
    """List biological assemblies stored in a PDB or mmCIF file.

    Each entry has an id, the author or software description, how many chains
    the expansion produces, and an atom estimate with and without crystal
    waters. An assembly whose operators are all the identity is reported as
    already the biological unit and is not expanded.
    """
    structure = _read(path)
    rows: list[dict[str, Any]] = []
    for assembly in structure.assemblies:
        if not _operators(assembly):
            continue
        rows.append(_assembly_row(structure, assembly))
    # Strict NCS is a different expansion. Offer it when the file has no
    # biological assembly that actually moves a chain.
    if not any(not row["already_biological_unit"] for row in rows) and _ncs_worth_expanding(structure):
        rows.append(_ncs_row(structure))
    return {"assemblies": rows, "message": _list_message(rows)}


def build_assembly(
    path: str | Path,
    assembly_id: str,
    include_waters: bool = False,
    drop_overlaps: bool = True,
    dest: str | Path | None = None,
) -> dict[str, Any]:
    """Write one biological assembly as mmCIF.

    Chain copies are named ``A1``, ``A2``, … so a later loader can keep them
    apart. Crystal waters (``HOH``, ``WAT``, ``DOD``) are dropped unless
    ``include_waters`` is set. Ligands stay. Identity operators do not
    duplicate atoms. The result is not a classic PDB: those chain names do
    not fit in one column.

    Ions and waters that sit on a symmetry axis are copied once per operator
    and land on the same point. ``drop_overlaps`` keeps one of each piled-up
    hetatm and leaves protein atoms in place.
    """
    structure = _read(path)
    assembly_key = str(assembly_id).strip()
    if assembly_key == _NCS_ID:
        return _build_ncs(
            structure,
            include_waters=include_waters,
            drop_overlaps=drop_overlaps,
            dest=dest,
            source=path,
        )
    assembly = _find_assembly(structure, assembly_key)
    if _all_identity(assembly):
        return _not_built(
            assembly_key,
            already=True,
            message="Every assembly operator is the identity, so there is nothing to build.",
        )
    model = _expanded_model(structure, assembly, include_waters=include_waters)
    message = ""
    if drop_overlaps:
        model, message = _without_stacked_copies(model)
    written = _write_mmcif(model, source=path, dest=dest)
    return _built_result(written, assembly_key, message=message)


def _read_mmcif_structure(path: Path) -> gemmi.Structure:
    """Read an mmCIF structure even when the file was saved with a ``.pdb`` name."""
    if path.suffix.lower() in {".cif", ".mmcif"}:
        return gemmi.read_structure(str(path))
    document = gemmi.cif.read(str(path))
    return gemmi.make_structure_from_block(document.sole_block())


def universe_from_assembly_mmcif(path: str | Path) -> mda.Universe:
    """Load an assembly mmCIF into MDAnalysis without shortening chain names.

    ``A1`` and ``B3`` stay on ``chainID``. Segment ids are left blank so a
    later chain lookup does not treat them as CHARMM segments and collapse
    them to one letter.
    """
    structure = _read_mmcif_structure(Path(path))
    if len(structure) == 0:
        raise AssemblyError("Assembly file has no model")
    model = structure[0]
    names: list[str] = []
    elements: list[str] = []
    positions: list[tuple[float, float, float]] = []
    resindex: list[int] = []
    chain_ids: list[str] = []
    residue_names: list[str] = []
    residue_ids: list[int] = []
    residue_icodes: list[str] = []
    current = None
    for chain in model:
        chain_name = str(chain.name)
        for residue in chain:
            icode = str(residue.seqid.icode or "").strip()
            key = (chain_name, str(residue.name).strip(), int(residue.seqid.num), icode)
            if key != current:
                current = key
                residue_names.append(key[1])
                residue_ids.append(key[2])
                residue_icodes.append(icode or "")
            residue_index = len(residue_names) - 1
            for atom in residue:
                names.append(str(atom.name).strip())
                element = str(atom.element.name or "").strip()
                elements.append(element or "X")
                point = atom.pos
                positions.append((float(point.x), float(point.y), float(point.z)))
                resindex.append(residue_index)
                chain_ids.append(chain_name)
    if not names:
        raise AssemblyError("Assembly file has no atoms")
    universe = mda.Universe.empty(
        len(names),
        n_residues=len(residue_names),
        atom_resindex=resindex,
        trajectory=True,
    )
    universe.add_TopologyAttr("name", names)
    universe.add_TopologyAttr("type", elements)
    universe.add_TopologyAttr("element", elements)
    universe.add_TopologyAttr("resname", residue_names)
    universe.add_TopologyAttr("resid", residue_ids)
    universe.add_TopologyAttr("icode", residue_icodes)
    universe.add_TopologyAttr("chainID", chain_ids)
    universe.add_TopologyAttr("segid", [""])
    universe.atoms.positions = np.asarray(positions, dtype=np.float64)
    return universe


def _read(path: str | Path) -> gemmi.Structure:
    file_path = Path(path)
    if not file_path.is_file():
        raise AssemblyError(f"Structure not found: {file_path}")
    try:
        structure = gemmi.read_structure(str(file_path))
    except (OSError, ValueError, RuntimeError) as exc:
        raise AssemblyError(f"Could not read structure: {exc}") from exc
    if len(structure) == 0:
        raise AssemblyError("Structure has no model")
    return structure


def _operators(assembly: gemmi.Assembly) -> list[Any]:
    operators = []
    for generator in assembly.generators:
        operators.extend(generator.operators)
    return operators


def _is_identity_transform(transform: gemmi.Transform) -> bool:
    matrix = transform.mat.tolist()
    for row in range(3):
        for col in range(3):
            expected = 1.0 if row == col else 0.0
            if abs(float(matrix[row][col]) - expected) > 1e-4:
                return False
    vector = transform.vec
    return abs(float(vector.x)) < 1e-3 and abs(float(vector.y)) < 1e-3 and abs(float(vector.z)) < 1e-3


def _all_identity(assembly: gemmi.Assembly) -> bool:
    operators = _operators(assembly)
    return bool(operators) and all(_is_identity_transform(op.transform) for op in operators)


def _ncs_worth_expanding(structure: gemmi.Structure) -> bool:
    if len(structure.ncs) == 0:
        return False
    for operator in structure.ncs:
        if getattr(operator, "given", False):
            continue
        if not _is_identity_transform(operator.tr):
            return True
    return False


def _description(assembly: gemmi.Assembly) -> str:
    detail = str(assembly.oligomeric_details or "").strip()
    who: list[str] = []
    if assembly.author_determined:
        who.append("author")
    if assembly.software_determined:
        who.append("software")
    if detail and who:
        return f"{detail} ({', '.join(who)})"
    if detail:
        return detail
    if who:
        return ", ".join(who)
    return ""


def _prepared_model(structure: gemmi.Structure, include_waters: bool) -> gemmi.Model:
    model = structure[0].clone()
    model.remove_alternative_conformations()
    if not include_waters:
        model.remove_waters()
    return model


def _expanded_model(
    structure: gemmi.Structure,
    assembly: gemmi.Assembly,
    include_waters: bool,
) -> gemmi.Model:
    model = _prepared_model(structure, include_waters)
    return gemmi.make_assembly(assembly, model, gemmi.HowToNameCopiedChain.AddNumber)


def _is_water(residue: gemmi.Residue) -> bool:
    return str(residue.name).strip().upper() in _WATER_NAMES


def _is_het(residue: gemmi.Residue) -> bool:
    return str(residue.het_flag).strip().upper() == "H" or _is_water(residue)


def _format_extra(counts: Counter[str]) -> str:
    return ", ".join(f"{name} × {count}" for name, count in sorted(counts.items()))


def _overlap_fields(
    other: Counter[str],
    water: Counter[str],
    protein: set[str],
) -> dict[str, Any]:
    combined = other + water
    return {
        "stacked_extra_atoms": int(sum(other.values())),
        "stacked_extra_atoms_with_waters": int(sum(combined.values())),
        "stacked_detail": _format_extra(other),
        "stacked_detail_with_waters": _format_extra(combined),
        "protein_overlap_chains": sorted(protein),
    }


def _stack_message(other: Counter[str], water: Counter[str], protein: set[str]) -> str:
    combined = other + water
    total = int(sum(combined.values()))
    parts: list[str] = []
    if total:
        noun = "atom" if total == 1 else "atoms"
        parts.append(
            f"Removed {total} extra {noun} that landed on the same point ({_format_extra(combined)})."
        )
    if protein:
        names = ", ".join(sorted(protein))
        parts.append(
            f"Protein atoms from {names} also land on the same point; those copies were kept."
        )
    return " ".join(parts)


def _scan_stacks(
    model: gemmi.Model,
    cutoff: float = _STACK_CUTOFF,
) -> tuple[set[tuple[int, int, int]], Counter[str], Counter[str], set[str]]:
    """Find hetatm atoms that repeat an earlier atom of the same element.

    Protein atoms are never marked for removal. A polymer atom within
    ``cutoff`` of a polymer atom on another chain is reported instead.
    """
    limit = cutoff * cutoff
    kept: dict[tuple[int, int, int], list[tuple[str, str, bool, float, float, float]]] = {}
    drop: set[tuple[int, int, int]] = set()
    extra_other: Counter[str] = Counter()
    extra_water: Counter[str] = Counter()
    protein: set[str] = set()

    def near(x: float, y: float, z: float):
        ix = math.floor(x / cutoff)
        iy = math.floor(y / cutoff)
        iz = math.floor(z / cutoff)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    yield from kept.get((ix + dx, iy + dy, iz + dz), ())

    def store(element: str, chain_name: str, polymer: bool, x: float, y: float, z: float) -> None:
        key = (math.floor(x / cutoff), math.floor(y / cutoff), math.floor(z / cutoff))
        kept.setdefault(key, []).append((element, chain_name, polymer, x, y, z))

    for chain_index, chain in enumerate(model):
        chain_name = str(chain.name)
        for residue_index, residue in enumerate(chain):
            het = _is_het(residue)
            water = _is_water(residue)
            resname = str(residue.name).strip()
            for atom_index, atom in enumerate(residue):
                element = str(atom.element.name or "").strip().upper()
                if not element:
                    element = str(atom.name).strip().upper()[:1] or "X"
                point = atom.pos
                x, y, z = float(point.x), float(point.y), float(point.z)
                same_element = False
                protein_hit = ""
                for kept_element, kept_chain, kept_polymer, kx, ky, kz in near(x, y, z):
                    if kept_chain == chain_name:
                        continue
                    dx = x - kx
                    dy = y - ky
                    dz = z - kz
                    if dx * dx + dy * dy + dz * dz >= limit:
                        continue
                    if element == kept_element:
                        same_element = True
                    if not het and kept_polymer:
                        protein_hit = kept_chain
                if het and same_element:
                    drop.add((chain_index, residue_index, atom_index))
                    if water:
                        extra_water[resname] += 1
                    else:
                        extra_other[resname] += 1
                    continue
                if protein_hit:
                    protein.add(chain_name)
                    protein.add(protein_hit)
                store(element, chain_name, not het, x, y, z)
    return drop, extra_other, extra_water, protein


def _residue_copy(residue: gemmi.Residue, atoms: list[gemmi.Atom]) -> gemmi.Residue:
    copied = gemmi.Residue()
    copied.name = residue.name
    copied.seqid = residue.seqid
    copied.het_flag = residue.het_flag
    copied.subchain = residue.subchain
    copied.entity_id = residue.entity_id
    copied.label_seq = residue.label_seq
    for atom in atoms:
        copied.add_atom(atom)
    return copied


def _model_without(model: gemmi.Model, drop: set[tuple[int, int, int]]) -> gemmi.Model:
    if not drop:
        return model
    fresh = gemmi.Model(1)
    for chain_index, chain in enumerate(model):
        new_chain = gemmi.Chain(chain.name)
        for residue_index, residue in enumerate(chain):
            kept = [
                atom
                for atom_index, atom in enumerate(residue)
                if (chain_index, residue_index, atom_index) not in drop
            ]
            if not kept:
                continue
            if len(kept) == len(residue):
                new_chain.add_residue(residue)
            else:
                new_chain.add_residue(_residue_copy(residue, kept))
        if len(new_chain):
            fresh.add_chain(new_chain)
    return fresh


def _without_stacked_copies(model: gemmi.Model) -> tuple[gemmi.Model, str]:
    drop, other, water, protein = _scan_stacks(model)
    return _model_without(model, drop), _stack_message(other, water, protein)


def _summarize(model: gemmi.Model) -> tuple[int, int, list[str]]:
    """Non-water atoms, all atoms, and chains that still hold a non-water residue."""
    other = 0
    water = 0
    chains: list[str] = []
    for chain in model:
        chain_other = 0
        for residue in chain:
            count = len(residue)
            if _is_water(residue):
                water += count
            else:
                other += count
                chain_other += count
        if chain_other:
            chains.append(str(chain.name))
    return other, other + water, chains


def _assembly_row(structure: gemmi.Structure, assembly: gemmi.Assembly) -> dict[str, Any]:
    identity = _all_identity(assembly)
    overlaps: dict[str, Any] | None = None
    if identity:
        dry_atoms, _ignored, chains = _summarize(_prepared_model(structure, include_waters=False))
        _dry, wet_atoms, _wet_chains = _summarize(_prepared_model(structure, include_waters=True))
    else:
        # Waters stay in this model so one expansion supplies both atom counts.
        expanded = _expanded_model(structure, assembly, include_waters=True)
        dry_atoms, wet_atoms, chains = _summarize(expanded)
        _drop, other, water, protein = _scan_stacks(expanded)
        overlaps = _overlap_fields(other, water, protein)
    return _row(
        assembly_id=str(assembly.name),
        description=_description(assembly),
        chains=chains,
        atom_count=dry_atoms,
        atom_count_with_waters=wet_atoms,
        already=identity,
        overlaps=overlaps,
    )


def _ncs_row(structure: gemmi.Structure) -> dict[str, Any]:
    wet = structure.clone()
    wet[0].remove_alternative_conformations()
    wet.expand_ncs(gemmi.HowToNameCopiedChain.AddNumber)
    dry_atoms, wet_atoms, chains = _summarize(wet[0])
    _drop, other, water, protein = _scan_stacks(wet[0])
    return _row(
        assembly_id=_NCS_ID,
        description="Strict NCS (MTRIX)",
        chains=chains,
        atom_count=dry_atoms,
        atom_count_with_waters=wet_atoms,
        already=False,
        overlaps=_overlap_fields(other, water, protein),
    )


def _row(
    assembly_id: str,
    description: str,
    chains: list[str],
    atom_count: int,
    atom_count_with_waters: int,
    already: bool,
    overlaps: dict[str, Any] | None = None,
) -> dict[str, Any]:
    stacked = overlaps or _overlap_fields(Counter(), Counter(), set())
    return {
        "id": assembly_id,
        "description": description,
        "oligomeric_count": len(chains),
        "atom_count": int(atom_count),
        "atom_count_with_waters": int(atom_count_with_waters),
        "already_biological_unit": already,
        "needs_confirmation": int(atom_count) >= LARGE_ASSEMBLY_ATOM_COUNT,
        "needs_confirmation_with_waters": int(atom_count_with_waters) >= LARGE_ASSEMBLY_ATOM_COUNT,
        **stacked,
    }


def _list_message(rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "This file has no biological assembly records, so there is nothing to build."
    if all(row["already_biological_unit"] for row in rows):
        return "Every assembly operator is the identity, so there is nothing to build."
    return ""


def _find_assembly(structure: gemmi.Structure, assembly_id: str) -> gemmi.Assembly:
    for assembly in structure.assemblies:
        if str(assembly.name) == assembly_id:
            return assembly
    known = ", ".join(str(assembly.name) for assembly in structure.assemblies) or "none"
    raise AssemblyError(f"No biological assembly {assembly_id!r} in this file (found {known}).")


def _not_built(assembly_id: str, already: bool, message: str) -> dict[str, Any]:
    return {
        "built": False,
        "path": "",
        "n_atoms": 0,
        "chains": [],
        "assembly_id": assembly_id,
        "already_biological_unit": already,
        "message": message,
    }


def _built_result(path: Path, assembly_id: str, message: str = "") -> dict[str, Any]:
    structure = gemmi.read_structure(str(path))
    _dry, atoms, chains = _summarize(structure[0])
    return {
        "built": True,
        "path": str(path),
        "n_atoms": int(atoms),
        "chains": chains,
        "assembly_id": assembly_id,
        "already_biological_unit": False,
        "message": message,
    }


def _build_ncs(
    structure: gemmi.Structure,
    include_waters: bool,
    drop_overlaps: bool,
    dest: str | Path | None,
    source: str | Path,
) -> dict[str, Any]:
    if not _ncs_worth_expanding(structure):
        return _not_built(
            _NCS_ID,
            already=True,
            message="Every NCS operator is the identity, so there is nothing to build.",
        )
    work = structure.clone()
    work[0].remove_alternative_conformations()
    if not include_waters:
        work[0].remove_waters()
    work.expand_ncs(gemmi.HowToNameCopiedChain.AddNumber)
    message = ""
    model = work[0]
    if drop_overlaps:
        model, message = _without_stacked_copies(model)
    written = _write_mmcif(model, source=source, dest=dest)
    return _built_result(written, _NCS_ID, message=message)


def _write_mmcif(model: gemmi.Model, source: str | Path, dest: str | Path | None) -> Path:
    out = gemmi.Structure()
    out.add_model(model)
    out.remove_empty_chains()
    if dest is None:
        stem = Path(source).stem or "assembly"
        handle, name = tempfile.mkstemp(prefix=f"{stem}_assembly_", suffix=".cif")
        os.close(handle)
        dest_path = Path(name)
    else:
        dest_path = Path(dest)
    document = out.make_mmcif_document()
    document.sole_block().set_pair(ASSEMBLY_MMCIF_TAG, "yes")
    document.write_file(str(dest_path))
    return dest_path
