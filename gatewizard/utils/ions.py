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
