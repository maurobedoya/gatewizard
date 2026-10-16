# gatewizard/utils/ions.py
# SPDX-License-Identifier: MIT
"""Ion residue names (Amber / CHARMM / PDB) and tleap neutralization helpers."""

from __future__ import annotations

from typing import Iterable, Optional

# Builder / packmol-memgen / tleap defaults
DEFAULT_CATION = "K+"
DEFAULT_ANION = "Cl-"

# Allowed Builder neutralization / salt ions (AmberTools naming)
AMBER_CATIONS = ("K+", "Na+", "Li+", "Rb+", "Cs+")
AMBER_ANIONS = ("Cl-", "Br-", "F-", "I-")

# Residue names seen in PDB / Amber / CHARMM / OpenMM exports.
# Matching is case-insensitive via ``is_ion_resname``.
ION_RESNAMES: tuple[str, ...] = (
    # Sodium
    "Na+",
    "NA+",
    "NA",
    "SOD",
    # Potassium
    "K+",
    "K",
    "POT",
    # Chloride
    "Cl-",
    "CL-",
    "CL",
    "Cl",
    "CLA",
    # Other common monoatomics
    "LI",
    "LI+",
    "LIT",
    "RB",
    "RB+",
    "RUB",
    "CS",
    "CS+",
    "CES",
    "MG",
    "MG2",
    "MG2+",
    "CA",
    "CA2",
    "CA2+",
    "CAL",
    "ZN",
    "ZN2",
    "ZN2+",
    "MN",
    "FE",
    "FE2",
    "FE3",
    "CU",
    "CU1",
    "CU2",
    "NI",
    "CD",
    "HG",
    "SR",
    "BA",
    "BAR",
    "F",
    "F-",
    "BR",
    "BR-",
    "I",
    "I-",
    "NH4",
)

_ION_RESNAMES_UPPER = frozenset(n.strip().upper() for n in ION_RESNAMES)


def normalize_ion_resname(name: Optional[str]) -> str:
    return str(name or "").strip().upper()


def is_ion_resname(name: Optional[str]) -> bool:
    """True if *name* is a known monoatomic / simple ion residue."""
    return normalize_ion_resname(name) in _ION_RESNAMES_UPPER


def ion_mda_selection(extra: Optional[Iterable[str]] = None) -> str:
    """MDAnalysis ``resname …`` selection covering common ion spellings."""
    names = list(ION_RESNAMES)
    if extra:
        names.extend(str(x) for x in extra if x)
    # Preserve order, drop empties / dupes (case-sensitive for MDA)
    seen: set[str] = set()
    ordered: list[str] = []
    for n in names:
        n = str(n).strip()
        if not n or n in seen:
            continue
        seen.add(n)
        ordered.append(n)
    return "resname " + " ".join(ordered)


def resolve_amber_cation(cation: Optional[str]) -> str:
    text = str(cation or "").strip()
    if not text:
        return DEFAULT_CATION
    # Accept case variants of known Amber names
    for known in AMBER_CATIONS:
        if text.upper() == known.upper():
            return known
    # Pass through if already looks like Amber (e.g. Na+)
    if text[-1:] == "+" or text in AMBER_CATIONS:
        # Normalize common typos
        upper = text.upper()
        if upper in ("NA+", "NA"):
            return "Na+"
        if upper in ("K+", "K"):
            return "K+"
        return text
    return DEFAULT_CATION


def resolve_amber_anion(anion: Optional[str]) -> str:
    text = str(anion or "").strip()
    if not text:
        return DEFAULT_ANION
    for known in AMBER_ANIONS:
        if text.upper() == known.upper():
            return known
    if text[-1:] == "-" or text in AMBER_ANIONS:
        upper = text.upper()
        if upper in ("CL-", "CL", "CLA"):
            return "Cl-"
        return text
    return DEFAULT_ANION


# Amber ``atomic_ions.lib`` atom name for a residue whose deposited atom
# name differs. Residue FE2 is deposited with atom FE; tleap then adds a
# second atom named FE2 and the original FE has no type. CU1, GD3, IOD,
# NH4, U4+, HE+, and HZ+ have the same mismatch when the atom is the
# residue name instead of the library atom.
AMBER_ION_ATOM_BY_RESNAME = {
    "FE2": "FE2",
    "CU1": "CU",
    "GD3": "GD",
    "IOD": "I",
    "NH4": "N",
    "U4+": "U",
    "HE+": "H",
    "HZ+": "H",
}

# Residue FE3 is Fe3+. This Amber library has no FE3 unit. The Fe3+ unit
# is residue FE, atom FE (charge +3). Residue FE2 stays FE2.
AMBER_ION_RESIDUE_ALIAS = {
    "FE3": "FE",
}


def align_amber_ion_atom_names(pdb_path: str) -> int:
    """Rename deposited ion atoms to the Amber residue template name.

    Residue FE3 is written as Amber residue FE, which is the Fe3+ unit.
    Other residues keep the deposited name. Returns the number of atom
    records changed.
    """
    path = pdb_path
    changed = 0
    lines_out: list[str] = []
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
                resname = line[17:20].strip().upper()
                target = AMBER_ION_RESIDUE_ALIAS.get(resname, resname)
                expected = AMBER_ION_ATOM_BY_RESNAME.get(target)
                if target != resname and expected is None:
                    expected = target
                current = line[12:16].strip().upper()
                if expected and (current != expected or target != resname):
                    atom = f"{expected:<4}"[:4]
                    residue = f"{target:<3}"[:3]
                    line = f"{line[:12]}{atom}{line[16:17]}{residue}{line[20:]}"
                    changed += 1
            lines_out.append(line if line.endswith("\n") else line + "\n")
    if changed:
        with open(path, "w", encoding="utf-8", newline="\n") as handle:
            handle.writelines(lines_out)
    return changed


def tleap_neutralization_lines(
    cation: Optional[str] = None,
    anion: Optional[str] = None,
) -> str:
    """tleap commands to neutralize with the chosen cation then anion.

    ``addIonsRand unit ion 0`` adds enough of *ion* to zero the charge when
    the ion and system have opposite signs; the other ion is a no-op when
    already neutralized.
    """
    cat = resolve_amber_cation(cation)
    an = resolve_amber_anion(anion)
    return (
        f"# Neutralize total charge (cation={cat}, anion={an})\n"
        f"addIonsRand system {cat} 0\n"
        f"addIonsRand system {an} 0"
    )
