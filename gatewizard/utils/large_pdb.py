"""Make a packed membrane PDB safe for pdb4amber and tleap.

pdb4amber 22 / ParmEd 4.3 read coordinates from columns 31–38, 39–46 and
47–54. A residue number past 9999 is written in more than four columns, so
those slices become text such as ``1 115.50``. A serial past 99999 is written
against the atom name (``272291Cl-``). The digit in ``1Cl-`` is the spilled
serial, not an atom name. Hydrogen names such as ``1HG`` stay.

Amber's residue field holds four digits. ParmEd treats residue 9999 followed
by any other number as overflow and then numbers every later residue 10000,
10001, …, which shifts the columns again. This writer stays at 9998. Deposited
numbers are kept while they are still unique. A repeated or oversized number
takes the next free value. After 9998 residues, later residues continue on
chains B, C, … and the numbers start again.

tleap's ``loadPDB`` does not copy CRYST1 into the unit, and a CONECT record
uses five-digit serials. On a large system those serials wrap, so the same
CONECT bonds a water hydrogen. CONECT is dropped. Disulfides are added with
``bond system.<resid>.SG``, which selects the first residue with that number.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

_DIGITS = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
_NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)")
_CHAINS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
MAX_RESID = 9998
_DISULFIDE_CUTOFF = 2.5
_KEPT_RECORDS = ("TER", "CRYST1", "REMARK", "TITLE", "HEADER", "END", "SSBOND")


def hy36encode(width: int, number: int) -> str:
    """Encode ``number`` in ``width`` columns, using hybrid-36 past 99999."""
    if number < 10**width:
        return f"{number:{width}d}"
    offset = number - 10**width + 10 * 36 ** (width - 1)
    chars: list[str] = []
    remaining = offset
    for _ in range(width):
        remaining, rem = divmod(remaining, 36)
        chars.append(_DIGITS[rem])
    if remaining or not chars[-1].isalpha():
        raise ValueError(f"{number} does not fit in hybrid-36 width {width}")
    return "".join(reversed(chars))


def _looks_like_atom_name(name: str) -> bool:
    """True for a PDB atom name, including ``1HG``. ``1Cl-`` is a spilled serial."""
    if not 1 <= len(name) <= 4:
        return False
    if name[0].isalpha():
        return True
    return len(name) >= 2 and name[0] in "123456789" and name[1] == "H"


def _split_serial_and_name(token: str) -> tuple[int, str | None]:
    if token.isdigit():
        return int(token), None
    for size in range(min(4, len(token) - 1), 0, -1):
        serial_text, name = token[:-size], token[-size:]
        if serial_text.isdigit() and len(serial_text) >= 5 and _looks_like_atom_name(name):
            return int(serial_text), name
    raise ValueError(f"could not read atom serial from {token!r}")


def _split_glued_numbers(token: str) -> list[str]:
    pieces = _NUMBER.findall(token)
    if len(pieces) > 1 and "".join(pieces) == token:
        return pieces
    return [token]


def _pdb_atom_name(name: str) -> str:
    text = name.strip()
    if len(text) >= 4:
        return text[:4]
    if len(text) == 1 or (len(text) <= 3 and text[:1].isalpha() and text[1:2].isdigit()):
        return f" {text:<3}"[:4]
    return f"{text:<4}"[:4]


def _columns_are_unsafe(line: str) -> bool:
    if len(line) < 54:
        return True
    try:
        float(line[30:38])
        float(line[38:46])
        float(line[46:54])
        resid = int(line[22:26])
    except ValueError:
        return True
    if resid > MAX_RESID or resid < -999:
        return True
    name = line[12:16].strip()
    if name and name[0].isdigit() and not (len(name) >= 2 and name[1] == "H"):
        return True
    return False


def needs_amber_pdb_rewrite(path: Path) -> bool:
    """True when pdb4amber cannot read this PDB's columns or residue numbers."""
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")) and _columns_are_unsafe(line):
                return True
    return False


_CHAIN_RESID = re.compile(r"^[A-Za-z](-?\d+)$")


def _split_record(parts: list[str]) -> tuple[str, list[str]]:
    """Separate the record name from a serial glued onto it, as in ``HETATM10000``."""
    head = parts[0]
    record = "ATOM  " if head.startswith("ATOM") else "HETATM"
    prefix = "ATOM" if head.startswith("ATOM") else "HETATM"
    rest = head[len(prefix) :]
    tail = parts[1:]
    if rest.isdigit():
        return record, [rest, *tail]
    if rest and rest[0].isdigit():
        _serial, glued = _split_serial_and_name(rest)
        if glued:
            return record, [str(_serial), glued, *tail]
        return record, [str(_serial), *tail]
    return record, tail


def _residue_number(token: str) -> int | None:
    """Read a residue number, including a chain letter stuck to it (``A1000``)."""
    if token.lstrip("+-").isdigit():
        return int(token)
    matched = _CHAIN_RESID.match(token)
    if matched:
        return int(matched.group(1))
    return None


def parse_atom(line: str) -> tuple[str, str, str, int, float, float, float, float, float]:
    """Read one atom from a column-shifted or whitespace PDB line."""
    parts = line.split()
    record, fields_in = _split_record(parts)
    _serial, glued_name = _split_serial_and_name(fields_in[0])
    if glued_name is None:
        raw_name = fields_in[1]
        fields = fields_in[2:]
    else:
        raw_name = glued_name
        fields = fields_in[1:]
    resname = fields[0][:3]
    index = 1
    resid = _residue_number(fields[index])
    if resid is None:
        index += 1
        resid = _residue_number(fields[index])
    if resid is None:
        raise ValueError(f"could not read residue number from {line[:90]!r}")
    numbers: list[str] = []
    for token in fields[index + 1 :]:
        numbers.extend(_split_glued_numbers(token))
    if len(numbers) < 3:
        raise ValueError(f"could not read coordinates from {line[:90]!r}")
    x, y, z = (float(numbers[0]), float(numbers[1]), float(numbers[2]))
    occ = float(numbers[3]) if len(numbers) > 3 else 1.0
    bfactor = float(numbers[4]) if len(numbers) > 4 else 0.0
    return record, raw_name, resname, resid, x, y, z, occ, bfactor


def _atom_line(
    record: str,
    serial: int,
    name: str,
    resname: str,
    chain: str,
    resid: int,
    x: float,
    y: float,
    z: float,
    occ: float,
    bfactor: float,
) -> str:
    # Four spaces after the residue number put x in column 31.
    return (
        f"{record}{hy36encode(5, serial)} {_pdb_atom_name(name)} "
        f"{resname:>3} {chain}{resid:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{bfactor:6.2f}\n"
    )


class _ResidueAllocator:
    """Keep deposited residue numbers until they would collide or exceed 9998."""

    def __init__(self) -> None:
        self.used: set[int] = set()
        self.next_free = 1
        self.overflow_slot = 0
        self.overflow_resid = 1

    def assign(self, deposited: int) -> tuple[str, int]:
        if (
            self.overflow_slot == 0
            and 1 <= deposited <= MAX_RESID
            and deposited not in self.used
        ):
            self.used.add(deposited)
            return "A", deposited
        if self.overflow_slot == 0:
            while self.next_free <= MAX_RESID and self.next_free in self.used:
                self.next_free += 1
            if self.next_free <= MAX_RESID:
                chosen = self.next_free
                self.next_free += 1
                self.used.add(chosen)
                return "A", chosen
            self.overflow_slot = 1
            self.overflow_resid = 1
        if self.overflow_resid > MAX_RESID:
            self.overflow_slot += 1
            self.overflow_resid = 1
        chains = _CHAINS[1:]
        if self.overflow_slot > len(chains):
            raise RuntimeError("Ran out of PDB chain ids for residue numbers")
        chain = chains[self.overflow_slot - 1]
        resid = self.overflow_resid
        self.overflow_resid += 1
        return chain, resid


def rewrite_pdb_for_amber(source: Path, dest: Path) -> dict[str, int | str]:
    """Write ``dest`` with Amber column widths. Atom order is unchanged."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    partial = dest.with_name(dest.name + ".writing")
    allocator = _ResidueAllocator()
    serial = 0
    previous: tuple[str, int] | None = None
    chain = "A"
    out_resid = 0
    moved = 0
    try:
        with source.open("r", encoding="utf-8", errors="replace") as inp, partial.open(
            "w", encoding="utf-8", newline="\n"
        ) as out:
            for line in inp:
                if line.startswith(("ATOM", "HETATM")):
                    record, raw_name, resname, resid, x, y, z, occ, bfactor = parse_atom(line)
                    key = (resname, resid)
                    if key != previous:
                        previous = key
                        chain, out_resid = allocator.assign(resid)
                        if chain != "A" or out_resid != resid:
                            moved += 1
                    serial += 1
                    out.write(
                        _atom_line(
                            record, serial, raw_name, resname, chain, out_resid, x, y, z, occ, bfactor
                        )
                    )
                elif line.startswith(_KEPT_RECORDS):
                    # TER separates chains. Dropping it makes tleap peptide-bond
                    # the next chain onto this one. CONECT is not copied.
                    out.write(line if line.endswith("\n") else line + "\n")
        partial.replace(dest)
    except Exception:
        partial.unlink(missing_ok=True)
        raise
    return {
        "atoms": serial,
        "residues_renumbered": moved,
        "last_chain": chain,
        "last_resid": out_resid,
    }


def read_cryst1_box(path: Path) -> tuple[float, float, float, float, float, float] | None:
    """Return ``(a, b, c, alpha, beta, gamma)`` from the first CRYST1 record."""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if not line.startswith("CRYST1") or len(line) < 54:
                    continue
                a = float(line[6:15])
                b = float(line[15:24])
                c = float(line[24:33])
                alpha = float(line[33:40] or 90.0)
                beta = float(line[40:47] or 90.0)
                gamma = float(line[47:54] or 90.0)
                return a, b, c, alpha, beta, gamma
    except (OSError, ValueError):
        return None
    return None


def _strip_conect(path: Path) -> int:
    kept: list[str] = []
    removed = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("CONECT"):
                removed += 1
                continue
            kept.append(line if line.endswith("\n") else line + "\n")
    if removed:
        path.write_text("".join(kept), encoding="utf-8", newline="\n")
    return removed


def _residue_key(resname: str, chain: str, resid: int) -> tuple[str, str, int]:
    return resname.strip(), chain.strip() or "A", resid


def _iter_column_atoms(path: Path):
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_index, line in enumerate(handle):
            if not line.startswith(("ATOM", "HETATM")) or len(line) < 54:
                continue
            name = line[12:16].strip()
            resname = line[17:20]
            chain = line[21:22]
            try:
                resid = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            yield line_index, name, resname, chain, resid, x, y, z


def _disulfide_pairs(path: Path) -> list[tuple[tuple[str, str, int], tuple[str, str, int]]]:
    sulfurs: list[tuple[tuple[str, str, int], float, float, float]] = []
    for _index, name, resname, chain, resid, x, y, z in _iter_column_atoms(path):
        if resname.strip().upper() == "CYX" and name.upper() == "SG":
            sulfurs.append((_residue_key(resname, chain, resid), x, y, z))
    pairs: list[tuple[float, tuple[str, str, int], tuple[str, str, int]]] = []
    for i, (left, x1, y1, z1) in enumerate(sulfurs):
        for right, x2, y2, z2 in sulfurs[i + 1 :]:
            distance = ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5
            if distance <= _DISULFIDE_CUTOFF:
                pairs.append((distance, left, right))
    pairs.sort()
    used: set[tuple[str, str, int]] = set()
    chosen: list[tuple[tuple[str, str, int], tuple[str, str, int]]] = []
    for _distance, left, right in pairs:
        if left in used or right in used or left == right:
            continue
        used.add(left)
        used.add(right)
        chosen.append((left, right))
    return chosen


def _renumber_residue(path: Path, chain: str, old: int, new: int) -> None:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    out: list[str] = []
    for line in lines:
        if (
            line.startswith(("ATOM", "HETATM"))
            and len(line) >= 26
            and line[21:22] == chain
            and line[22:26].strip().lstrip("-").isdigit()
            and int(line[22:26]) == old
        ):
            line = f"{line[:22]}{new:4d}{line[26:]}"
        out.append(line if line.endswith("\n") else line + "\n")
    path.write_text("".join(out), encoding="utf-8", newline="\n")


def _first_occurrence(path: Path) -> dict[int, tuple[str, str, int]]:
    found: dict[int, tuple[str, str, int]] = {}
    seen: set[tuple[str, str, int]] = set()
    for _index, _name, resname, chain, resid, _x, _y, _z in _iter_column_atoms(path):
        key = _residue_key(resname, chain, resid)
        if key in seen:
            continue
        seen.add(key)
        found.setdefault(resid, key)
    return found


def _ensure_bond_residue_numbers(path: Path, pairs: list) -> list[tuple[int, int]]:
    """Return residue numbers tleap can bond. ``system.N`` is the first N."""
    bonds: list[tuple[int, int]] = []
    for left, right in pairs:
        numbers: list[int] = []
        for key in (left, right):
            owner = _first_occurrence(path)
            _resname, chain, resid = key
            if owner.get(resid) == key:
                numbers.append(resid)
                continue
            used_before = set(owner)
            free = next((n for n in range(1, MAX_RESID + 1) if n not in used_before), None)
            if free is None:
                raise RuntimeError(
                    f"Disulfide {key[0]} {chain}:{resid} repeats a residue number "
                    "and no unused number is left for tleap."
                )
            _renumber_residue(path, chain, resid, free)
            numbers.append(free)
        bonds.append((numbers[0], numbers[1]))
    return bonds


def geometry_lines_for_tleap(path: Path) -> list[str]:
    """Bond commands and ``set system box`` for the PDB tleap will load."""
    lines = [
        "# Disulfides and the periodic box. loadPDB does not copy CRYST1,",
        "# and a wrapped CONECT serial bonds the wrong atom.",
    ]
    for left, right in _ensure_bond_residue_numbers(path, _disulfide_pairs(path)):
        lines.append(f"bond system.{left}.SG system.{right}.SG")
    box = read_cryst1_box(path)
    if box is not None:
        a, b, c, alpha, beta, gamma = box
        if abs(alpha - 90.0) < 0.01 and abs(beta - 90.0) < 0.01 and abs(gamma - 90.0) < 0.01:
            lines.append(f"set system box {{ {a:.3f} {b:.3f} {c:.3f} }}")
        else:
            lines.append(
                f"set system box {{ {a:.3f} {b:.3f} {c:.3f} {alpha:.2f} {beta:.2f} {gamma:.2f} }}"
            )
    return lines


def finalize_for_tleap(pdb_path: Path, fragment_path: Path) -> dict[str, int]:
    """Drop CONECT, keep disulfide bonds and the CRYST1 box for tleap."""
    if needs_amber_pdb_rewrite(pdb_path):
        rewritten = pdb_path.with_name(pdb_path.name + ".aligned")
        rewrite_pdb_for_amber(pdb_path, rewritten)
        rewritten.replace(pdb_path)
    removed = _strip_conect(pdb_path)
    fragment_path.write_text(
        "\n".join(geometry_lines_for_tleap(pdb_path)) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    box = read_cryst1_box(pdb_path)
    return {"conect_removed": removed, "has_box": int(box is not None)}


def _prepend_cryst1(pdb_path: Path, source: Path) -> bool:
    if read_cryst1_box(pdb_path) is not None:
        return False
    box_line = None
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("CRYST1"):
                box_line = line if line.endswith("\n") else line + "\n"
                break
    if box_line is None:
        return False
    body = pdb_path.read_text(encoding="utf-8", errors="replace")
    pdb_path.write_text(box_line + body, encoding="utf-8", newline="\n")
    return True


def repair_saved_pdb(pdb_path: Path, cryst1_source: Path) -> dict[str, int | bool]:
    """Make tleap's ``system.pdb`` readable and keep the packed box on it.

    Restraint selection reads this PDB by residue name. Atom order matches
    the topology. tleap's own ``savePDB`` writes residue 100000 again, which
    MDAnalysis cannot read, so the columns are rewritten when they are unsafe.
    """
    rewritten = False
    if needs_amber_pdb_rewrite(pdb_path):
        aligned = pdb_path.with_name(pdb_path.name + ".aligned")
        rewrite_pdb_for_amber(pdb_path, aligned)
        aligned.replace(pdb_path)
        rewritten = True
    copied = _prepend_cryst1(pdb_path, cryst1_source)
    return {"rewritten": rewritten, "cryst1_copied": copied}


def rewrite_if_needed(source: Path, dest: Path) -> bool:
    """Rewrite ``source`` to ``dest`` when Amber cannot read it.

    A readable file leaves no ``dest``, so the packed PDB is used as-is.
    """
    if not needs_amber_pdb_rewrite(source):
        if dest.is_file():
            dest.unlink()
        return False
    rewrite_pdb_for_amber(source, dest)
    return True


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    rewrite = sub.add_parser("rewrite-if-needed")
    rewrite.add_argument("source", type=Path)
    rewrite.add_argument("dest", type=Path)

    finalize = sub.add_parser("finalize-for-tleap")
    finalize.add_argument("pdb", type=Path)
    finalize.add_argument("fragment", type=Path)

    repair = sub.add_parser("repair-saved-pdb")
    repair.add_argument("pdb", type=Path)
    repair.add_argument("cryst1_source", type=Path)

    args = parser.parse_args(argv)
    if args.command == "rewrite-if-needed":
        changed = rewrite_if_needed(args.source, args.dest)
        print("rewritten" if changed else "unchanged")
        return 0
    if args.command == "finalize-for-tleap":
        info = finalize_for_tleap(args.pdb, args.fragment)
        print(
            f"conect_removed {info['conect_removed']} has_box {info['has_box']}"
        )
        return 0
    info = repair_saved_pdb(args.pdb, args.cryst1_source)
    print(f"rewritten {int(info['rewritten'])} cryst1_copied {int(info['cryst1_copied'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
