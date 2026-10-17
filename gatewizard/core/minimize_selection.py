"""Minimize a chosen set of atoms.

Protein and Amber ions use OpenMM. Every other atom Amber can build stays
fixed. Selected residues Amber cannot build get a short repulsive bump, with
the rest of the structure held fixed. The rest of the file is copied unchanged.
"""

from __future__ import annotations

import os
import tempfile
from typing import Any, Dict, List, Mapping, Optional, Sequence, Set, Tuple

import MDAnalysis as mda
import numpy as np
from openmm import LangevinMiddleIntegrator, unit
from openmm.app import CutoffNonPeriodic, ForceField, Modeller, PDBFile, Simulation

FAST_ITERATIONS = 50
MAX_ITERATIONS = 2000
DEFAULT_TOLERANCE = 10.0
NEIGHBOR_A = 8.0
ResidueKey = Tuple[str, int, str]
_LJ_EPS = 2.0
_SPRING_K = 80.0
_MAX_STEP_A = 0.15
_MIN_STEP_A = 1e-4
# Bondi van der Waals radii, angstroms. Sigma for a pair is the sum.
_VDW = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "FE": 2.00,
    "ZN": 1.39,
    "MG": 1.73,
    "CA": 2.31,
    "NA": 2.27,
    "K": 2.75,
}
_COVALENT = {
    "H": 0.31,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "P": 1.07,
    "S": 1.05,
    "CL": 1.02,
    "BR": 1.20,
    "I": 1.39,
    "FE": 1.32,
    "ZN": 1.22,
    "MG": 1.41,
    "CA": 1.76,
    "NA": 1.66,
    "K": 2.03,
}

_ALIASES = {"HIS": "HID", "HSD": "HID", "HSE": "HIE", "HSP": "HIP"}
# Monoatomic ions in amber14/tip3p.xml. OpenMM's protein file does not include them.
_ION_ALIASES = {
    "K": "K",
    "K+": "K",
    "POT": "K",
    "NA": "NA",
    "NA+": "NA",
    "SOD": "NA",
    "CL": "CL",
    "CL-": "CL",
    "CLA": "CL",
    "MG": "MG",
    "MG2": "MG",
    "CA": "CA",
    "CA2": "CA",
    "ZN": "ZN",
    "ZN2": "ZN",
    "LI": "LI",
    "RB": "RB",
    "CS": "CS",
    "F": "F",
    "BR": "BR",
    "IOD": "IOD",
    "I": "IOD",
}
_AMBER = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "ASH",
    "CYS",
    "CYM",
    "CYX",
    "GLN",
    "GLU",
    "GLH",
    "GLY",
    "HID",
    "HIE",
    "HIP",
    "ILE",
    "LEU",
    "LYS",
    "LYN",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
    "ACE",
    "NME",
}


def minimize_selection(
    input_pdb: str,
    output_pdb: str,
    selection: Optional[str] = None,
    atom_indices: Optional[Sequence[int]] = None,
    residues: Optional[Sequence[Mapping[str, Any]]] = None,
    max_iterations: int = FAST_ITERATIONS,
    tolerance: Optional[float] = None,
    max_rmsd: Optional[float] = None,
) -> Dict[str, Any]:
    """Move the selected atoms. Everything else stays fixed.

    Protein and Amber ions use OpenMM. ``max_iterations`` is the L-BFGS limit
    and the steric step count. ``tolerance`` is the OpenMM force tolerance in
    kJ/mol/nm; the steric bump uses the same number, converted per angstrom.
    ``max_rmsd`` is an optional cap, in angstroms, on how far the selected
    atoms of each pass may move. The default is a short pass with no cap.
    """
    if not os.path.isfile(input_pdb):
        raise FileNotFoundError(f"Structure not found: {input_pdb}")
    steps = max(1, min(int(max_iterations or FAST_ITERATIONS), MAX_ITERATIONS))
    tol = None if tolerance is None else float(tolerance)
    cap = None if max_rmsd is None else float(max_rmsd)
    universe = mda.Universe(input_pdb)
    chosen, bump, relaxed = _chosen_atoms(universe, selection, atom_indices, residues)
    if len(chosen) == 0 and not bump:
        raise RuntimeError("Nothing to minimize. The selection has no atoms.")
    moved = 0
    if len(chosen):
        context = _amber_context(universe)
        movable = {int(index) for index in chosen.indices}
        positions = _cap_rmsd(universe, _minimize_context(context, movable, steps, tol), cap)
        moved += _apply_positions(universe, positions)
    if bump:
        positions = _steric_bump(universe, bump, steps, tol, cap)
        moved += _apply_positions(universe, positions)
    _write_pdb(input_pdb, output_pdb, universe)
    return {
        "output_file": os.path.abspath(output_pdb),
        "atoms_moved": moved,
        "movable": (len(chosen) if len(chosen) else 0) + len(bump),
        "skipped": [],
        "relaxed": relaxed,
        "max_iterations": steps,
        "tolerance": DEFAULT_TOLERANCE if tol is None else tol,
        "max_rmsd": cap,
    }


def _amber_resname(resname: str) -> Optional[str]:
    raw = str(resname or "").strip().upper()
    ion = _ION_ALIASES.get(raw)
    if ion:
        return ion
    mapped = _ALIASES.get(raw, raw)
    if mapped in _AMBER:
        return mapped
    return None


def _fragment_chain(atom: Any, resname: str) -> str:
    """Keep ions off the protein chain so the amino acid terminus still matches Amber."""
    if resname in set(_ION_ALIASES.values()):
        return "X"
    chain = (_residue_key(atom)[0] or "A")[:1]
    return chain or "A"


def _fragment_atom_name(resname: str, atom_name: str) -> str:
    """Amber ion templates use one atom whose name matches the residue."""
    if resname in set(_ION_ALIASES.values()):
        return resname
    return str(atom_name)


def _chosen_atoms(
    universe: mda.Universe,
    selection: Optional[str],
    atom_indices: Optional[Sequence[int]],
    residues: Optional[Sequence[Mapping[str, Any]]],
) -> Tuple[mda.AtomGroup, List[int], List[Dict[str, Any]]]:
    if atom_indices:
        valid = [int(index) for index in atom_indices if 0 <= int(index) < len(universe.atoms)]
        group = universe.atoms[valid] if valid else universe.atoms[[]]
    elif residues:
        parts = []
        for item in residues:
            chain = str(item.get("chain") or "").strip()
            resid = int(item.get("resid"))
            if chain:
                parts.append(f"(chainID {chain} and resid {resid})")
            else:
                parts.append(f"(resid {resid})")
        group = universe.select_atoms(" or ".join(parts)) if parts else universe.atoms[[]]
    elif selection and str(selection).strip() and str(selection).strip().lower() != "all":
        group = universe.select_atoms(str(selection).strip())
    else:
        raise RuntimeError("Choose residues, atoms, or a molecule before minimizing.")
    kept: List[int] = []
    bump: List[int] = []
    relaxed: List[Dict[str, Any]] = []
    seen: Set[ResidueKey] = set()
    for atom in group.atoms:
        key = _residue_key(atom)
        if _amber_resname(atom.resname):
            kept.append(int(atom.index))
            continue
        bump.append(int(atom.index))
        if key not in seen:
            seen.add(key)
            relaxed.append(
                {
                    "chain": key[0],
                    "resid": key[1],
                    "resname": key[2],
                    "reason": "Relaxed sterically",
                }
            )
    return universe.atoms[kept] if kept else universe.atoms[[]], bump, relaxed


def _amber_context(universe: mda.Universe) -> mda.AtomGroup:
    """Protein and ion atoms Amber can build. Unknown residues stay out of the model."""
    indices: List[int] = []
    for residue in universe.residues:
        if _amber_resname(str(residue.resname)) is None:
            continue
        indices.extend(int(index) for index in residue.atoms.indices)
    return universe.atoms[indices] if indices else universe.atoms[[]]


def _minimize_context(
    context: mda.AtomGroup,
    movable: Set[int],
    steps: int,
    tolerance: Optional[float],
) -> Dict[int, np.ndarray]:
    fragment = _fragment_pdb(context)
    _append_missing_template_hydrogens(fragment)
    try:
        pdb = PDBFile(fragment)
        _restore_amber_residue_names(pdb.topology, fragment)
        forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
        modeller = Modeller(pdb.topology, pdb.positions)
        # A structure with no hydrogens gets a full set here. A prepared
        # structure keeps the hydrogens it already has; any that are missing
        # were added to the temporary file above and are not written back.
        if not _fragment_has_hydrogen(pdb.topology):
            modeller.addHydrogens(forcefield)
        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=CutoffNonPeriodic,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=None,
            ignoreExternalBonds=True,
        )
    except Exception as exc:
        raise RuntimeError(f"OpenMM could not build this selection: {exc}") from exc
    finally:
        try:
            os.remove(fragment)
        except OSError:
            pass

    by_name = _original_by_name(context)
    for atom in modeller.topology.atoms():
        original = _take_original(by_name, atom)
        if original is None or original not in movable:
            system.setParticleMass(atom.index, 0.0)
    integrator = LangevinMiddleIntegrator(
        300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds
    )
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    if tolerance is None or tolerance <= 0:
        simulation.minimizeEnergy(maxIterations=steps)
    else:
        simulation.minimizeEnergy(
            tolerance=float(tolerance) * unit.kilojoule_per_mole / unit.nanometer,
            maxIterations=steps,
        )
    state = simulation.context.getState(getPositions=True)
    coords = state.getPositions(asNumpy=True).value_in_unit(unit.angstrom)
    by_name = _original_by_name(context)
    updated: Dict[int, np.ndarray] = {}
    for atom in modeller.topology.atoms():
        original = _take_original(by_name, atom)
        if original is not None and original in movable:
            updated[original] = np.asarray(coords[atom.index], dtype=np.float64)
    return updated


def _cap_rmsd(
    universe: mda.Universe,
    positions: Mapping[int, np.ndarray],
    max_rmsd: Optional[float],
) -> Dict[int, np.ndarray]:
    """Scale a displacement so the selected atoms sit on the RMSD cap."""
    if max_rmsd is None or not positions:
        return dict(positions)
    before = []
    delta = []
    keys = []
    for index, xyz in positions.items():
        start = np.asarray(universe.atoms[int(index)].position, dtype=np.float64)
        shift = np.asarray(xyz, dtype=np.float64) - start
        keys.append(int(index))
        before.append(start)
        delta.append(shift)
    shifts = np.vstack(delta)
    rmsd = float(np.sqrt(np.mean(np.sum(shifts * shifts, axis=1))))
    limit = max(0.0, float(max_rmsd))
    if rmsd <= limit or rmsd == 0.0:
        return dict(positions)
    scale = limit / rmsd
    return {index: start + scale * shift for index, start, shift in zip(keys, before, delta)}


def _steric_bump(
    universe: mda.Universe,
    movable: Sequence[int],
    steps: int,
    tolerance: Optional[float],
    max_rmsd: Optional[float] = None,
) -> Dict[int, np.ndarray]:
    """Repulsive Lennard-Jones plus bond springs. No attractive well."""
    coords = np.asarray(universe.atoms.positions, dtype=np.float64).copy()
    mov = np.asarray([int(index) for index in movable], dtype=int)
    mov_set = set(int(index) for index in mov)
    radii = np.asarray([_vdw_radius(atom) for atom in universe.atoms], dtype=np.float64)
    bonds = _movable_bonds(universe, mov_set, coords)
    excluded = _excluded_pairs(universe, mov_set, bonds)
    angles = _pairs_two_bonds_apart(bonds)
    bond_rest = {
        pair: float(np.linalg.norm(coords[pair[0]] - coords[pair[1]]))
        for pair in bonds
    }
    angle_rest = {
        pair: float(np.linalg.norm(coords[pair[0]] - coords[pair[1]]))
        for pair in angles
    }
    start = coords.copy()
    force_limit = (DEFAULT_TOLERANCE if tolerance is None or tolerance <= 0 else float(tolerance)) / 10.0
    step = 0.05
    for _ in range(steps):
        energy, forces = _steric_forces(coords, mov, mov_set, radii, bonds, bond_rest, excluded)
        peak = float(np.max(np.linalg.norm(forces[mov], axis=1))) if len(mov) else 0.0
        if peak < force_limit:
            break
        scale = step / max(peak, 1e-8)
        delta = forces[mov] * scale
        lengths = np.linalg.norm(delta, axis=1)
        too_far = lengths > _MAX_STEP_A
        delta[too_far] *= (_MAX_STEP_A / lengths[too_far])[:, None]
        trial = coords.copy()
        trial[mov] += delta
        _restore_geometry(trial, bonds, bond_rest, angles, angle_rest)
        trial_energy, _trial_forces = _steric_forces(trial, mov, mov_set, radii, bonds, bond_rest, excluded)
        if trial_energy <= energy + 1e-6:
            coords = trial
            step = min(step * 1.2, _MAX_STEP_A)
        else:
            step *= 0.5
            if step < _MIN_STEP_A:
                break
    _restore_geometry(coords, bonds, bond_rest, angles, angle_rest)
    if max_rmsd is not None:
        coords = _cap_movable_coords(start, coords, mov, max_rmsd)
        _restore_geometry(coords, bonds, bond_rest, angles, angle_rest)
    return {int(index): coords[int(index)].copy() for index in mov}


def _pairs_two_bonds_apart(bonds: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Atom pairs separated by two bonds. For water this is the H-H pair that sets the angle."""
    graph: Dict[int, Set[int]] = {}
    for left, right in bonds:
        graph.setdefault(left, set()).add(right)
        graph.setdefault(right, set()).add(left)
    found: Set[Tuple[int, int]] = set()
    for neighbors in graph.values():
        ordered = sorted(neighbors)
        for pos, left in enumerate(ordered):
            for right in ordered[pos + 1 :]:
                found.add((left, right))
    return sorted(found)


def _restore_geometry(
    coords: np.ndarray,
    bonds: Sequence[Tuple[int, int]],
    bond_rest: Mapping[Tuple[int, int], float],
    angles: Sequence[Tuple[int, int]],
    angle_rest: Mapping[Tuple[int, int], float],
) -> None:
    """Hold bond lengths and the distances that set bond angles."""
    for _ in range(12):
        bond_error = _project_distances(coords, bonds, bond_rest)
        angle_error = _project_distances(coords, angles, angle_rest)
        if bond_error < 1e-4 and angle_error < 1e-4:
            break


def _project_distances(
    coords: np.ndarray,
    pairs: Sequence[Tuple[int, int]],
    rest: Mapping[Tuple[int, int], float],
) -> float:
    worst = 0.0
    for left, right in pairs:
        vector = coords[left] - coords[right]
        length = float(np.linalg.norm(vector))
        if length < 1e-8:
            continue
        target = rest[(left, right)]
        worst = max(worst, abs(length - target))
        shift = (0.5 * (length - target) / length) * vector
        coords[left] -= shift
        coords[right] += shift
    return worst


def _cap_movable_coords(
    start: np.ndarray,
    coords: np.ndarray,
    movable: np.ndarray,
    max_rmsd: float,
) -> np.ndarray:
    shifts = coords[movable] - start[movable]
    rmsd = float(np.sqrt(np.mean(np.sum(shifts * shifts, axis=1)))) if len(movable) else 0.0
    limit = max(0.0, float(max_rmsd))
    if rmsd <= limit or rmsd == 0.0:
        return coords
    capped = coords.copy()
    capped[movable] = start[movable] + (limit / rmsd) * shifts
    return capped


def _steric_forces(
    coords: np.ndarray,
    mov: np.ndarray,
    mov_set: Set[int],
    radii: np.ndarray,
    bonds: Sequence[Tuple[int, int]],
    rest: Mapping[Tuple[int, int], float],
    excluded: Set[Tuple[int, int]],
) -> Tuple[float, np.ndarray]:
    forces = np.zeros_like(coords)
    energy = 0.0
    for ia in mov:
        delta = coords[ia] - coords
        dist = np.linalg.norm(delta, axis=1)
        close = np.where((dist < NEIGHBOR_A) & (dist > 1e-8))[0]
        for other in close:
            pair = (int(ia), int(other)) if int(ia) < int(other) else (int(other), int(ia))
            if int(other) in mov_set and int(other) <= int(ia):
                continue
            if pair in excluded:
                continue
            separation = max(float(dist[other]), 0.5)
            contact = float(radii[ia] + radii[other])
            if separation >= contact:
                continue
            sigma = contact / (2.0 ** (1.0 / 6.0))
            sr = sigma / separation
            sr6 = sr ** 6
            sr12 = sr6 * sr6
            energy += 4.0 * _LJ_EPS * (sr12 - sr6) + _LJ_EPS
            slope = (24.0 * _LJ_EPS / separation) * (-2.0 * sr12 + sr6)
            push = (-slope / separation) * delta[other] * (separation / float(dist[other]))
            forces[ia] += push
            if int(other) in mov_set:
                forces[int(other)] -= push
    for left, right in bonds:
        vector = coords[left] - coords[right]
        length = float(np.linalg.norm(vector))
        if length < 1e-8:
            continue
        stretch = length - rest[(left, right)]
        energy += 0.5 * _SPRING_K * stretch * stretch
        slope = _SPRING_K * stretch
        pull = (-slope / length) * vector
        forces[left] += pull
        forces[right] -= pull
    return energy, forces


def _element_symbol(atom: Any) -> str:
    element = str(getattr(atom, "element", "") or "").strip().upper()
    if element and element not in {"X"}:
        return element
    name = str(getattr(atom, "name", "") or "").strip().upper()
    for symbol in ("CL", "BR", "FE", "ZN", "MG", "NA"):
        if name.startswith(symbol):
            return symbol
    for char in name:
        if char.isalpha():
            return char
    return "C"


def _vdw_radius(atom: Any) -> float:
    return _VDW.get(_element_symbol(atom), 1.70)


def _covalent_radius(atom: Any) -> float:
    return _COVALENT.get(_element_symbol(atom), 0.77)


def _movable_bonds(
    universe: mda.Universe,
    movable: Set[int],
    coords: np.ndarray,
) -> List[Tuple[int, int]]:
    """Bonds inside the movable set, from CONECT and from a distance guess.

    A file that already lists protein bonds must not hide the water O-H pairs.
    Those missing pairs were still repelling, which stretched the hydrogens.
    """
    found: Set[Tuple[int, int]] = set()
    bonds = getattr(universe, "bonds", None)
    if bonds is not None and len(bonds):
        for bond in bonds:
            left = int(bond.atoms[0].index)
            right = int(bond.atoms[1].index)
            if left in movable and right in movable and left != right:
                found.add((min(left, right), max(left, right)))
    ordered = sorted(movable)
    atoms = list(universe.atoms)
    for pos, left in enumerate(ordered):
        for right in ordered[pos + 1 :]:
            length = float(np.linalg.norm(coords[left] - coords[right]))
            if 0.4 < length <= _bond_limit(atoms[left], atoms[right]):
                found.add((left, right))
    return sorted(found)


def _bond_limit(left: Any, right: Any) -> float:
    """Longest distance still treated as a bond, in angstroms.

    HOH hydrogens in a prepared PDB sit near 1.18 Å. A strict covalent sum
    stops at about 1.16 Å, so those O-H pairs were missed and then repelled.
    """
    total = _covalent_radius(left) + _covalent_radius(right)
    hydrogen = _element_symbol(left) == "H" or _element_symbol(right) == "H"
    if hydrogen and _same_residue(left, right):
        return 1.45 * total
    return 1.2 * total


def _same_residue(left: Any, right: Any) -> bool:
    return (
        int(left.resid) == int(right.resid)
        and str(getattr(left, "chainID", "") or "") == str(getattr(right, "chainID", "") or "")
        and str(left.resname).strip().upper() == str(right.resname).strip().upper()
    )


def _excluded_pairs(
    universe: mda.Universe,
    movable: Set[int],
    guessed: Sequence[Tuple[int, int]],
) -> Set[Tuple[int, int]]:
    graph: Dict[int, Set[int]] = {}

    def _add(left: int, right: int) -> None:
        graph.setdefault(left, set()).add(right)
        graph.setdefault(right, set()).add(left)

    for left, right in guessed:
        _add(left, right)
    bonds = getattr(universe, "bonds", None)
    if bonds is not None and len(bonds):
        for bond in bonds:
            left = int(bond.atoms[0].index)
            right = int(bond.atoms[1].index)
            if left in movable or right in movable:
                _add(left, right)
    excluded: Set[Tuple[int, int]] = set()
    for atom in movable:
        for neighbor in graph.get(atom, ()):
            excluded.add((min(atom, neighbor), max(atom, neighbor)))
            for further in graph.get(neighbor, ()):
                if further != atom:
                    excluded.add((min(atom, further), max(atom, further)))
    return excluded


def _apply_positions(universe: mda.Universe, positions: Mapping[int, np.ndarray]) -> int:
    moved = 0
    for index, xyz in positions.items():
        before = np.asarray(universe.atoms[index].position, dtype=np.float64)
        if float(np.linalg.norm(xyz - before)) > 0.01:
            moved += 1
        universe.atoms[index].position = xyz
    return moved


def _write_pdb(source: str, dest: str, universe: mda.Universe) -> None:
    by_index = {int(atom.index): np.asarray(atom.position, dtype=np.float64) for atom in universe.atoms}
    order = [int(atom.index) for atom in universe.atoms]
    cursor = 0
    lines: List[str] = []
    with open(source, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and cursor < len(order):
                xyz = by_index[order[cursor]]
                cursor += 1
                text = line.rstrip("\n")
                if len(text) < 54:
                    text = text.ljust(54)
                coords = f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                lines.append(text[:30] + coords + text[54:] + "\n")
            else:
                lines.append(line if line.endswith("\n") else line + "\n")
    with open(dest, "w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(lines)


def _chain_continues_after(atom: Any) -> bool:
    """True when a later protein residue on this chain is part of the model."""
    chain, resid, _resname = _residue_key(atom)
    ions = set(_ION_ALIASES.values())
    for residue in atom.universe.residues:
        other = _residue_key(residue.atoms[0])
        if other[0] != chain or other[1] <= resid:
            continue
        mapped = _amber_resname(other[2])
        if mapped and mapped not in ions:
            return True
    return False


def _restore_amber_residue_names(topology: Any, path: str) -> None:
    """Put GLH, ASH, HIE, and CYX back after OpenMM reads them as GLU, ASP, HIS, and CYS."""
    wanted: Dict[Tuple[str, str], str] = {}
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21].strip() or "A"
        wanted[(chain, line[22:26].strip())] = line[17:20].strip()
    for residue in topology.residues():
        chain = str(residue.chain.id).strip() or "A"
        name = wanted.get((chain, str(residue.id).strip()))
        if name and residue.name != name:
            residue.name = name


def _fragment_has_hydrogen(topology: Any) -> bool:
    for atom in topology.atoms():
        element = getattr(atom, "element", None)
        symbol = str(getattr(element, "symbol", "") or "").upper()
        if symbol == "H":
            return True
    return False


def _append_missing_template_hydrogens(path: str) -> None:
    """Add hydrogens a prepared residue is short, only in the temporary model.

    OpenMM refuses a threonine that has every atom except the hydroxyl
    hydrogen. Those atoms are not copied into the written PDB.
    """
    forcefield = ForceField("amber14-all.xml", "amber14/tip3p.xml")
    residues: List[Dict[str, Any]] = []
    other: List[str] = []
    current: Optional[Dict[str, Any]] = None
    for line in open(path, encoding="utf-8", errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            if current is not None:
                residues.append(current)
                current = None
            other.append(line)
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip() or "A"
        resid = int(line[22:26])
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        key = (chain, resid, resname)
        if current is None or current["key"] != key:
            if current is not None:
                residues.append(current)
            current = {"key": key, "names": {}, "serial": int(line[6:11])}
        current["names"][name] = xyz
        current["serial"] = max(int(current["serial"]), int(line[6:11]))
    if current is not None:
        residues.append(current)
    added = False
    for residue in residues:
        names = {name for name in residue["names"] if name != "OXT"}
        template = _template_missing_hydrogens(forcefield, residue["key"][2], names)
        if template is None:
            continue
        present = set(residue["names"])
        if not any(_is_hydrogen_name(name) for name in present):
            continue
        template_names = {atom.name for atom in template.atoms}
        for name in list(residue["names"]):
            if name not in template_names and _is_hydrogen_name(name):
                del residue["names"][name]
                added = True
        present = set(residue["names"])
        missing = [atom for atom in template.atoms if atom.name not in present and atom.element.symbol == "H"]
        for slot, atom in enumerate(missing):
            xyz = _placed_hydrogen(template, atom.name, residue["names"], slot)
            if xyz is None:
                continue
            residue["names"][atom.name] = xyz
            added = True
    if not added:
        return
    end = [line for line in other if line.startswith("END")]
    kept = [line for line in other if not line.startswith("END")]
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(kept)
        handle.writelines(_residue_lines(residues))
        handle.writelines(end or ["END\n"])


def _residue_lines(residues: Sequence[Mapping[str, Any]]) -> List[str]:
    lines: List[str] = []
    serial = 1
    for residue in residues:
        chain, resid, resname = residue["key"]
        for name, xyz in residue["names"].items():
            lines.append(
                f"ATOM  {serial:5d} {_pdb_name(name)} {resname:>3} {chain}{resid:4d}    "
                f"{float(xyz[0]):8.3f}{float(xyz[1]):8.3f}{float(xyz[2]):8.3f}  1.00  0.00\n"
            )
            serial += 1
    return lines


def _is_hydrogen_name(name: str) -> bool:
    text = name.strip().upper()
    return text.startswith("H") or (len(text) > 1 and text[0].isdigit() and "H" in text)


def _template_missing_hydrogens(forcefield: ForceField, resname: str, names: Set[str]) -> Any:
    heavy = {name for name in names if name != "OXT" and not _is_hydrogen_name(name)}
    hydrogens = {name for name in names if _is_hydrogen_name(name)}
    found = []
    for key, template in forcefield._templates.items():
        if template.name != resname and not str(key).endswith(resname):
            continue
        template_names = {atom.name for atom in template.atoms}
        template_heavy = {atom.name for atom in template.atoms if atom.element.symbol != "H"}
        if not heavy <= template_heavy:
            continue
        known = heavy | {name for name in hydrogens if name in template_names}
        missing = [atom for atom in template.atoms if atom.name not in known and atom.name != "OXT"]
        if any(atom.element.symbol != "H" for atom in missing):
            continue
        extra = len([name for name in hydrogens if name not in template_names])
        found.append((len(missing), extra, template.name != resname, template))
    if not found:
        return None
    found.sort(key=lambda item: (item[0], item[1], item[2]))
    if found[0][0] == 0 and found[0][1] == 0:
        return None
    return found[0][3]


def _placed_hydrogen(
    template: Any,
    atom_name: str,
    positions: Mapping[str, np.ndarray],
    slot: int,
) -> Optional[np.ndarray]:
    index = template.getAtomIndexByName(atom_name)
    parent = None
    for left, right in template.bonds:
        other_index = right if left == index else left if right == index else None
        if other_index is None:
            continue
        other = template.atoms[other_index]
        if other.element.symbol != "H" and other.name in positions:
            parent = other
            break
    if parent is None:
        return None
    origin = np.asarray(positions[parent.name], dtype=np.float64)
    neighbors = []
    parent_index = template.getAtomIndexByName(parent.name)
    for left, right in template.bonds:
        other_index = right if left == parent_index else left if right == parent_index else None
        if other_index is None or other_index == index:
            continue
        other = template.atoms[other_index]
        if other.element.symbol != "H" and other.name in positions:
            neighbors.append(np.asarray(positions[other.name], dtype=np.float64))
    if neighbors:
        direction = origin - neighbors[0]
    else:
        direction = np.array([1.0, 0.0, 0.0])
    length = float(np.linalg.norm(direction))
    if length < 1e-6:
        direction = np.array([1.0, 0.0, 0.0])
        length = 1.0
    direction = direction / length
    if slot:
        side = np.array([-direction[1], direction[0], 0.0])
        side_length = float(np.linalg.norm(side))
        if side_length < 1e-6:
            side = np.array([0.0, 1.0, 0.0])
        else:
            side = side / side_length
        direction = direction + side * (0.6 * slot)
        direction = direction / float(np.linalg.norm(direction))
    return origin + 0.96 * direction


def _fragment_pdb(context: mda.AtomGroup) -> str:
    fd, path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    lines = ["REMARK GateWizard minimize selection\n"]
    serial = 1
    protein_by_chain: Dict[str, List[str]] = {}
    ion_lines: List[str] = []
    last_by_chain: Dict[str, Tuple[int, np.ndarray, str]] = {}
    has_oxt: Set[Tuple[str, int]] = set()
    ions = set(_ION_ALIASES.values())
    has_hydrogen = False
    for atom in context.atoms:
        element = str(getattr(atom, "element", "") or "").strip().upper()
        if element == "H" or str(atom.name).strip().upper().startswith("H"):
            has_hydrogen = True
            break
    for atom in context.atoms:
        resname = _amber_resname(str(atom.resname)) or "UNK"
        chain = _fragment_chain(atom, resname)
        name = _pdb_name(_fragment_atom_name(resname, str(atom.name)))
        x, y, z = (float(value) for value in atom.position)
        record = (
            f"ATOM  {serial:5d} {name} {resname:>3} {chain}{int(atom.resid):4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00\n"
        )
        serial += 1
        if resname in ions:
            ion_lines.append(record)
            continue
        protein_by_chain.setdefault(chain, []).append(record)
        if str(atom.name).strip() == "C" and (not has_hydrogen or not _chain_continues_after(atom)):
            last_by_chain[chain] = (int(atom.resid), np.asarray(atom.position, dtype=np.float64), resname)
        if str(atom.name).strip() == "OXT":
            has_oxt.add((chain, int(atom.resid)))
    for chain, records in protein_by_chain.items():
        lines.extend(records)
        if chain not in last_by_chain:
            continue
        resid, carbon, resname = last_by_chain[chain]
        if (chain, resid) in has_oxt:
            continue
        oxygen = carbon + np.array([1.2, 0.0, 0.0])
        lines.append(
            f"ATOM  {serial:5d}  OXT {resname:>3} {chain}{resid:4d}    "
            f"{oxygen[0]:8.3f}{oxygen[1]:8.3f}{oxygen[2]:8.3f}  1.00  0.00\n"
        )
        serial += 1
    lines.extend(ion_lines)
    lines.append("END\n")
    with open(path, "w", encoding="utf-8", newline="\n") as handle:
        handle.writelines(lines)
    return path


def _original_by_name(shell: mda.AtomGroup) -> Dict[Tuple[str, int, str], List[int]]:
    found: Dict[Tuple[str, int, str], List[int]] = {}
    for atom in shell.atoms:
        mapped = _amber_resname(str(atom.resname)) or ""
        chain = _fragment_chain(atom, mapped)
        atom_name = _fragment_atom_name(mapped, str(atom.name)) if mapped else str(atom.name).strip()
        key = (chain, int(atom.resid), atom_name.strip())
        found.setdefault(key, []).append(int(atom.index))
    return found


def _take_original(by_name: Dict[Tuple[str, int, str], List[int]], atom: Any) -> Optional[int]:
    chain = str(atom.residue.chain.id).strip()[:1] or "A"
    key = (chain, _residue_number(atom.residue.id), str(atom.name).strip())
    bucket = by_name.get(key)
    if not bucket:
        return None
    return bucket.pop(0)


def _residue_key(atom: Any) -> ResidueKey:
    chain = str(getattr(atom, "chainID", "") or "").strip()
    if not chain or chain.upper() in {"SYSTEM", "SEGID"}:
        chain = str(getattr(atom, "segid", "") or "").strip()
    if chain.upper() in {"SYSTEM", "SEGID"}:
        chain = ""
    return chain, int(atom.resid), str(atom.resname).strip().upper()


def _residue_number(value: Any) -> int:
    text = str(value).strip()
    return int(text or "0")


def _pdb_name(name: str) -> str:
    text = name.strip()
    if len(text) >= 4:
        return text[:4]
    if len(text) == 1 or (len(text) <= 3 and text[:1].isalpha() and len(text) > 1 and text[1].isdigit()):
        return f" {text:<3}"[:4]
    return f"{text:<4}"[:4]
