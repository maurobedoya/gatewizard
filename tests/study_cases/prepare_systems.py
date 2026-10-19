#!/usr/bin/env python3
"""Prepare membrane systems and time every stage.

The builder (packmol-memgen + tleap) runs ``BUILDER_REPLICAS`` times.
Download, chain selection, MemPrO, protonation, and GAFF run once per
PDB and are reused by those replicas. Pass ``--fresh`` to ignore saved
stage times and run them again.

Run this file from the directory where ``work/`` and ``results/`` should
be written. Edit ``SYSTEMS`` before launching. The study-case page
describes each choice. A system with ``assembly`` builds that biological
assembly before chain selection. ``chains`` is the chain list that enters
the build.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import statistics
import subprocess
import sys
import time
import traceback
import urllib.request
from dataclasses import replace
from datetime import datetime
from pathlib import Path

import numpy as np

from gatewizard.core.assembly import build_assembly, universe_from_assembly_mmcif
from gatewizard.core.builder import Builder
from gatewizard.core.mempro import MemPrO, compute_orientation_transform
from gatewizard.core.preparation import (
    PreparationError,
    PreparationManager,
    complete_missing_heavy_atoms,
    resolve_reduce_executable,
    strip_protein_hydrogens,
)
from gatewizard.core.short_gaps import fill_short_gaps
from gatewizard.core.structure_checker import diagnose_structure
from gatewizard.core.structure_manager import StructureManager
from gatewizard.tools.equilibration import (
    AmberEquilibrationManager,
    GROMACSEquilibrationManager,
    NAMDEquilibrationManager,
    OpenMMEquilibrationManager,
    _build_universal_membrane_stages,
)
from gatewizard.tools.ligand_parametrization import (
    forcefield_charge_warning,
    parametrize_ligand_from_system_pdb,
)
from gatewizard.utils.ions import is_ion_resname
from gatewizard.utils.large_pdb import read_cryst1_box
from gatewizard.utils.peptide_residues import PEPTIDE_POLYMER_RESIDUES, remap_pdb_peptide_resnames

# Outputs follow the working directory, so a checkout of this file can be
# launched from any benchmark folder.
WORK = Path.cwd() / "work"
RESULTS = Path.cwd() / "results"

# How many independent packmol-memgen + tleap runs to time per PDB.
BUILDER_REPLICAS = 3

# PropKa target used when residue names are rewritten before pdb4amber.
PH = 7.4

# Membrane MD after a successful build. This stays off, so running the
# script only prepares the system. Turn it on for one PDB by setting
# ``protocol`` on that entry in SYSTEMS. ``engine`` is amber, gromacs,
# openmm, or namd. An empty ``engine_executable`` is looked up on PATH
# (openmm uses this Python). ``stage_ns`` replaces a stage length in
# nanoseconds; None keeps the GateWizard membrane protocol. That protocol
# is minimization (10000 steps), Equilibration 1–6, then production.
# Equilibration 6 is 47.625 ns. A shorter test can set, for example,
# ``"stage_ns": {"Equilibration 6": 0.1}``.
# COM and rotation restraints stay off, matching the Equilibration page.
PROTOCOL = {
    "enabled": False,
    "engine": "openmm",
    "engine_executable": "",
    "openmm_platform": "",
    "minimization": True,
    "equilibration": True,
    "production_ns": 0.0,
    "ensemble": "NPT",
    "temperature": 310.15,
    "mini_steps": 10000,
    "stage_ns": None,
    "replicas": 1,
}

# Residues that stay out of GAFF. A kept name in this set needs
# ``manual_params`` (frcmod + lib) on that system.
NEEDS_MANUAL_PARAMS = frozenset({"OEX", "HEM"})

# Deposited water names. A name in keep_hetero is kept with the protein.
# The builder still adds bulk TIP3P around those waters.
WATER_NAMES = frozenset({"HOH", "WAT", "TIP", "SOL", "DOD"})

# Organic ligands. Amber recommends GAFF with AM1-BCC, or GAFF2 with ABCG2.
# Set atom_type and charge_method on a system to change that pair. A
# parametrize entry can also be a dict with charge, atom_type, and
# charge_method for one residue. GAFF2 with AM1-BCC, and GAFF with ABCG2,
# still run and print a warning.
LIGAND_ATOM_TYPE = "gaff2"
LIGAND_CHARGE_METHOD = "abcg2"

# Cofactors and thylakoid lipids in 3WU2. They belong to the complex.
# This example leaves them out so the protein-and-ions build stays shorter.
# Add a name to keep_hetero when that run should include it.
PSII_ORGANIC = (
    "CLA",
    "PHO",
    "BCR",
    "RRX",
    "PL9",
    "BCT",
    "LHG",
    "DGD",
    "LMG",
    "SQD",
    "HEM",
    "OEX",
)

# Deposited choices, from the smallest packed system to the largest.
# ``parametrize`` maps a kept organic residue name to its net charge.
# atom_type and charge_method above choose GAFF or GAFF2 and ABCG2 or
# AM1-BCC. Ions are kept without a charge entry.
SYSTEMS = [
    {
        "pdb_id": "2MVH",
        "title": "SpoVM",
        # 10 NMR models. RCSB representative conformer is model 1.
        "model": 1,
        "chains": "A",
        "keep_hetero": [],
        "parametrize": {},
        # Amphipathic peptide. Transmembrane MemPrO is the wrong mode.
        "peripheral": True,
        # Smallest system. The membrane protocol is commented out, so this
        # entry only prepares the system. Uncomment the dict to run OpenMM
        # minimization, Equilibration 1–6, and 1 ns of production.
        # "protocol": {
        #     "enabled": True,
        #     "engine": "openmm",
        #     "minimization": True,
        #     "equilibration": True,
        #     "production_ns": 1.0,
        #     "openmm_platform": "CUDA",
        # },
    },
    {
        "pdb_id": "1AFO",
        "title": "Glycophorin A transmembrane dimer",
        # 20 NMR models. Model 1 is the first deposited conformer.
        "model": 1,
        # Biological assembly is the A+B dimer. One chain would be a monomer.
        "chains": "AB",
        "keep_hetero": [],
        "parametrize": {},
        "peripheral": False,
    },
    {
        "pdb_id": "1M0L",
        "title": "Bacteriorhodopsin",
        "model": 1,
        # Deposited asymmetric unit is one chain. The biological assembly
        # is a trimer generated by crystal symmetry; this script does not
        # build that trimer.
        "chains": "A",
        # RET stays. Crystal lipids LI1 and SQU are left out of this example
        # because POPC replaces the bilayer. Keep a lipid that is part of
        # the biology, and give it a charge in parametrize.
        # The deposited formula is neutral retinal (C20H28O), but the oxygen
        # is absent: LINK bonds C15 to Lys216 NZ. After hydrogens are added
        # the residue is C20H27. Charge 0 then has an odd electron count and
        # antechamber refuses a singlet. The protonated Schiff base is +1,
        # which is the charge passed to antechamber. The lysine stays LYS, and the
        # LINK is not rebuilt as a bond.
        "keep_hetero": ["RET"],
        "parametrize": {"RET": 1},
        "peripheral": False,
    },
    {
        "pdb_id": "7AAQ",
        "title": "Sugar/H+ symporter STP10",
        "model": 1,
        "chains": "A",
        # BGC is the bound sugar. Acetate (ACT), PEG detergent (XPE), and
        # monoolein (OLC) are left out of this example. Keep a detergent or
        # a lipid when it is biologically relevant, and parametrize it.
        "keep_hetero": ["BGC"],
        "parametrize": {"BGC": 0},
        "peripheral": False,
    },
    {
        "pdb_id": "3T4D",
        "title": "NaK2K channel Y55F",
        "model": 1,
        # Biological assembly 2 expands deposited chain B, the more complete
        # monomer (96 of 97 SEQRES residues). Assembly 1 expands chain A.
        # Copies are named B1–B4 and then renamed A–D for the one-column PDB.
        "assembly": {"id": "2", "waters": True, "drop_overlaps": True},
        "assembly_from": "B",
        "chains": "ABCD",
        # Crystal waters stay. Symmetry copies of the pore potassium that land
        # on the same point are dropped. The four ions that remain stay.
        "keep_hetero": ["HOH", "K"],
        "parametrize": {},
        "peripheral": False,
    },
    {
        "pdb_id": "6RV2",
        "title": "TASK-1 (K2P3.1)",
        "model": 1,
        # The asymmetric unit has two dimers (A+B and C+D). Assembly 1 is
        # chains A and B. Switch to "CD" for the other dimer.
        "chains": "AB",
        # Selectivity-filter potassium stays. Crystal PC (PC1), detergent
        # (DMU), and cholesterol hemisuccinate (Y01) are left out of this
        # example. Keep a detergent or a crystal lipid when it is
        # biologically relevant.
        "keep_hetero": ["K"],
        "parametrize": {},
        "peripheral": False,
    },
    {
        "pdb_id": "3WU2",
        "title": "Photosystem II",
        "model": 1,
        # Both halves of the dimer (uppercase and lowercase author chains).
        "chains": "AaBbCcDdEeFfHhIiJjKkLlMmOoTtUuVvXxYyZz",
        # Monoatomic ions stay. Chlorophyll, heme, the Mn4Ca cluster, and
        # the thylakoid lipids are part of the complex. They are left out
        # here so this build stays shorter. Add them from PSII_ORGANIC
        # when the question needs them. HEM and OEX also need manual_params.
        "keep_hetero": ["CA", "CL", "MG", "FE2"],
        "parametrize": {},
        "manual_params": {},
        "peripheral": False,
    },
]


def main() -> None:
    args = parse_args()
    replicas = args.replicas if args.replicas is not None else BUILDER_REPLICAS
    if args.report:
        saved = RESULTS / "timings.json"
        rows = json.loads(saved.read_text(encoding="utf-8")) if saved.is_file() else []
        write_report(rows, replicas)
        print(f"Wrote {RESULTS / 'timings.csv'}")
        print(f"Wrote {RESULTS / 'summary.md'}")
        return
    selected = {item.upper() for item in args.systems} if args.systems else None
    rows = []
    for spec in SYSTEMS:
        if selected and spec["pdb_id"].upper() not in selected:
            continue
        rows.extend(run_system(spec, replicas=replicas, fresh=args.fresh))
    merged = merge_saved_rows(rows)
    write_report(merged, replicas)
    print(f"Wrote {RESULTS / 'timings.csv'}")
    print(f"Wrote {RESULTS / 'summary.md'}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "systems",
        nargs="*",
        help="PDB ids to run (default: every entry in SYSTEMS)",
    )
    parser.add_argument(
        "--replicas",
        type=int,
        default=None,
        help=f"Builder repeats per PDB (default: {BUILDER_REPLICAS})",
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Ignore stage JSON files and time every stage again",
    )
    parser.add_argument(
        "--report",
        action="store_true",
        help="Rewrite results/timings.csv from saved stages and build directories. Do not start a stage.",
    )
    return parser.parse_args()


def run_system(spec: dict, replicas: int, fresh: bool) -> list[dict]:
    pdb_id = spec["pdb_id"]
    folder = WORK / pdb_id
    folder.mkdir(parents=True, exist_ok=True)
    print(f"\n=== {pdb_id} {spec['title']} ===")
    rows: list[dict] = []

    def timed(stage: str, fn, replica: int | None = None):
        row = run_stage(folder, pdb_id, stage, fn, replica=replica, fresh=fresh)
        rows.append(row)
        persist_rows(rows, replicas)
        label = stage if replica is None else f"{stage} replica {replica}"
        if row["ok"]:
            print(f"  {label}: {row['seconds']:.1f} s")
        else:
            print(f"  {label} failed: {row['error']}")
        return row

    source = folder / f"{pdb_id}.pdb"
    timed("download", lambda: download_pdb(pdb_id, source))
    if not source.is_file():
        return rows

    gap_spec = dict(spec)
    if spec.get("assembly_from"):
        gap_spec["chains"] = spec["assembly_from"]
    gapped = folder / "short_gaps.pdb"
    row = timed("short_gaps", lambda: complete_short_gaps(source, gapped, gap_spec))
    if not row["ok"]:
        return rows
    subset_source = Path(row["path"])

    if spec.get("assembly"):
        assembled = folder / "assembly.pdb"
        row = timed(
            "assembly",
            lambda: write_biological_assembly(subset_source, assembled, spec),
        )
        if not row["ok"]:
            return rows
        subset_source = Path(row["path"])

    subset = folder / "subset.pdb"
    protein_only = folder / "protein_for_mempro.pdb"
    row = timed(
        "subset",
        lambda: write_subset(subset_source, subset, protein_only, spec),
    )
    if not row["ok"]:
        return rows

    timed(
        "structure_check_before",
        lambda: write_structure_check(
            subset,
            folder / "structure_check_before",
            seqres_from=subset_source if spec.get("assembly") else source,
            chains=set(spec["chains"]),
        ),
    )

    oriented = folder / "oriented.pdb"
    row = timed(
        "mempro",
        lambda: orient_structure(protein_only, subset, oriented, spec, folder / "mempro"),
    )
    if not row["ok"]:
        return rows

    prepared = folder / "prepared.pdb"
    keep_water = bool(WATER_NAMES.intersection(str(name).upper() for name in spec.get("keep_hetero") or []))
    row = timed(
        "protonation",
        lambda: protonate(oriented, prepared, folder / "protonation", keep_water=keep_water),
    )
    if not row["ok"]:
        return rows

    timed(
        "structure_check_after",
        lambda: write_structure_check(prepared, folder / "structure_check_after"),
    )

    row = timed(
        "ligand_parameters",
        lambda: parametrize_kept_ligands(prepared, spec, folder),
    )
    if not row["ok"]:
        return rows
    ligand_params = row.get("ligand_params") or {}

    protein_atoms = count_records(prepared)[0]
    protocol = resolved_protocol(spec)
    for replica in range(1, replicas + 1):
        job_name = f"builder_rep{replica}"
        row = timed(
            "builder",
            lambda replica=replica, job_name=job_name: build_once(
                prepared,
                folder / job_name,
                ligand_params,
                protein_atoms,
            ),
            replica=replica,
        )
        if (
            row.get("ok")
            and protocol["enabled"]
            and replica <= int(protocol.get("replicas") or 1)
        ):
            timed(
                "protocol",
                lambda replica=replica, job_name=job_name, protocol=protocol: run_protocol(
                    folder / job_name,
                    protocol,
                ),
                replica=replica,
            )
    return rows


def run_stage(folder: Path, pdb_id: str, stage: str, fn, replica: int | None, fresh: bool) -> dict:
    marker_name = f"{stage}.json" if replica is None else f"{stage}_rep{replica}.json"
    marker = folder / marker_name
    if marker.is_file() and not fresh:
        saved = json.loads(marker.read_text(encoding="utf-8"))
        if saved.get("ok"):
            return saved
    started = time.perf_counter()
    try:
        payload = fn() or {}
        ok = True
        error = ""
    except Exception as exc:
        payload = {}
        ok = False
        error = f"{exc}\n{traceback.format_exc()}"
    row = {
        "pdb_id": pdb_id,
        "stage": stage,
        "replica": replica,
        "seconds": time.perf_counter() - started,
        "ok": ok,
        "error": error,
    }
    row.update(payload)
    marker.write_text(json.dumps(row, indent=2), encoding="utf-8")
    return row


def download_pdb(pdb_id: str, dest: Path) -> dict:
    dest.parent.mkdir(parents=True, exist_ok=True)
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    urllib.request.urlretrieve(url, dest)
    return {"path": str(dest)}


_CHAIN_LETTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def write_biological_assembly(source: Path, dest: Path, spec: dict) -> dict:
    """Expand one biological assembly and write a one-column PDB.

    This is the Visualize assembly build: crystal waters stay when
    ``waters`` is set, and atoms that land on the same point are dropped
    when ``drop_overlaps`` is set. Copied chains are named ``B1``, ``B2``,
    and so on. Those names do not fit in the PDB chain column, so they are
    renamed to ``A``, ``B``, ``C``, … in that order. SEQRES for
    ``assembly_from`` is copied onto each new chain.
    """
    assembly = spec.get("assembly") or {}
    assembly_id = str(assembly.get("id") or "").strip()
    if not assembly_id:
        raise RuntimeError(f"{spec['pdb_id']} assembly needs an id.")
    cif = dest.with_suffix(".cif")
    built = build_assembly(
        source,
        assembly_id,
        include_waters=bool(assembly.get("waters", True)),
        drop_overlaps=bool(assembly.get("drop_overlaps", True)),
        dest=cif,
    )
    if not built.get("built"):
        raise RuntimeError(built.get("message") or f"Assembly {assembly_id} was not built.")
    universe = universe_from_assembly_mmcif(cif)
    order: list[str] = []
    for atom in universe.atoms:
        chain = str(atom.chainID).strip()
        if chain and chain not in order:
            order.append(chain)
    if len(order) > len(_CHAIN_LETTERS):
        raise RuntimeError(
            f"{spec['pdb_id']} assembly {assembly_id} has {len(order)} chains. "
            "The PDB chain column holds one character."
        )
    letters = {chain: _CHAIN_LETTERS[index] for index, chain in enumerate(order)}
    seqres = _seqres_lines(source, str(spec.get("assembly_from") or ""))
    lines = [f"TITLE     {spec['pdb_id']} biological assembly {assembly_id}\n"]
    for chain in order:
        for line in seqres:
            lines.append(line[:11] + letters[chain] + line[12:])
    serial = 1
    waters = 0
    for atom in universe.atoms:
        chain = letters[str(atom.chainID).strip()]
        resname = str(atom.resname).strip()[:3]
        if resname in WATER_NAMES:
            waters += 1
        raw_name = str(atom.name).strip()
        name = raw_name if len(raw_name) >= 4 else f" {raw_name:<3s}"
        element = str(getattr(atom, "element", "") or "").strip().upper()
        if not element:
            element = "".join(char for char in raw_name if char.isalpha())[:2].upper() or "X"
        element = element[:2]
        record = "HETATM" if resname not in PEPTIDE_POLYMER_RESIDUES else "ATOM  "
        x, y, z = (float(value) for value in atom.position)
        lines.append(
            f"{record}{serial:5d} {name} {resname:>3s} {chain}"
            f"{int(atom.resid):4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
            f"{1.0:6.2f}{0.0:6.2f}          {element:>2s}\n"
        )
        serial += 1
    lines.append("END\n")
    dest.write_text("".join(lines), encoding="utf-8")
    renamed = ",".join(f"{chain}={letters[chain]}" for chain in order)
    message = built.get("message") or ""
    if message:
        print(f"    {message}")
    print(f"    chains {renamed}; {waters} water atoms")
    return {
        "path": str(dest),
        "assembly_id": assembly_id,
        "n_atoms": built["n_atoms"],
        "chain_map": renamed,
        "waters": waters,
        "message": message,
    }


def _seqres_lines(source: Path, chain: str) -> list[str]:
    if not chain:
        return []
    lines = []
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True):
        if line.startswith("SEQRES") and len(line) >= 12 and line[11] == chain:
            lines.append(line if line.endswith("\n") else line + "\n")
    return lines


def complete_short_gaps(source: Path, dest: Path, spec: dict) -> dict:
    """Build numbering gaps of 1 to 6 standard residues before orientation.

    SEQRES is still on the downloaded PDB. The subset written afterward
    drops SEQRES, so the fill has to happen here. Only the new residue
    atoms are inserted into the downloaded PDB. Ligands, waters, and the
    other records stay as deposited. A gap that is longer than 6 residues,
    or whose sequence is not a standard amino acid, is printed and left open.
    """
    report = fill_short_gaps(str(source), str(dest), chains=list(spec["chains"]))
    built = report.get("built") or []
    skipped = report.get("skipped") or []
    for gap in built:
        start = int(gap["after_resid"]) + 1
        end = int(gap["before_resid"]) - 1
        print(f"    chain {gap['chain']} {start}-{end}: {' '.join(gap['sequence'])}")
    for gap in skipped:
        print(
            f"    left open chain {gap['chain']}, {gap['count']} residues "
            f"({gap['reason']})"
        )
    return {
        "path": str(dest),
        "filled": len(built),
        "skipped": len(skipped),
    }


def write_subset(source: Path, subset: Path, protein_only: Path, spec: dict) -> dict:
    keep = {name.upper() for name in spec.get("keep_hetero") or []}
    unknown = sorted(
        name
        for name in keep
        if name in NEEDS_MANUAL_PARAMS and name not in (spec.get("manual_params") or {})
    )
    if unknown:
        raise RuntimeError(
            f"{spec['pdb_id']} keeps {', '.join(unknown)}, which need "
            "manual_params (frcmod and lib). GAFF is not used for those residues."
        )
    missing_charge = sorted(
        name
        for name in keep
        if name not in WATER_NAMES
        and not is_ion_resname(name)
        and name not in (spec.get("parametrize") or {})
        and name not in (spec.get("manual_params") or {})
    )
    if missing_charge:
        raise RuntimeError(
            f"{spec['pdb_id']} keeps {', '.join(missing_charge)} without a "
            "GAFF charge in parametrize or a library in manual_params."
        )

    chains = set(spec["chains"])
    model = int(spec.get("model") or 1)
    kept_lines = []
    protein_lines = []
    current = 1
    saw_model = False
    for line in source.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True):
        if line.startswith("MODEL"):
            saw_model = True
            current = int(line.split()[1])
            continue
        if line.startswith("ENDMDL"):
            continue
        if saw_model and current != model:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21:22]
        if chain not in chains:
            continue
        resname = line[17:20].strip().upper()
        if resname in WATER_NAMES and resname not in keep:
            continue
        if line.startswith("HETATM") and resname not in keep:
            continue
        kept_lines.append(line if line.endswith("\n") else line + "\n")
        if line.startswith("ATOM"):
            protein_lines.append(kept_lines[-1])

    if not protein_lines:
        raise RuntimeError(f"No protein atoms kept for {spec['pdb_id']} chains {spec['chains']}.")
    subset.write_text("".join(kept_lines) + "END\n", encoding="utf-8")
    protein_only.write_text("".join(protein_lines) + "END\n", encoding="utf-8")
    protein_atoms, hetero_atoms = count_records(subset)
    return {
        "protein_atoms": protein_atoms,
        "hetero_atoms": hetero_atoms,
        "path": str(subset),
    }


def orient_structure(
    protein_pdb: Path,
    full_pdb: Path,
    dest: Path,
    spec: dict,
    mempro_dir: Path,
) -> dict:
    mempro_dir.mkdir(parents=True, exist_ok=True)
    results = MemPrO().run(
        str(protein_pdb),
        output_dir=str(mempro_dir),
        peripheral=bool(spec.get("peripheral")),
    )
    if not results:
        raise RuntimeError(f"MemPrO returned no orientations for {spec['pdb_id']}.")
    rotation, translation = compute_orientation_transform(str(protein_pdb), results[0].pdb_path)
    manager = StructureManager()
    manager.load_structure(full_pdb)
    for atom in manager.structure.atoms:
        atom.coord = rotation @ np.asarray(atom.coord, dtype=float) + translation
    manager.structure.write_pdb(dest)
    return {"path": str(dest), "rank": results[0].rank}


def protonate(source: Path, dest: Path, folder: Path, *, keep_water: bool = False) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    manager = PreparationManager()
    named = folder / "protonated.pdb"
    pka_path = manager.run_analysis(str(source), output_dir=str(folder))
    summary_path = manager.extract_summary(pka_path, str(folder / "propka_summary.txt"))
    try:
        manager.parse_summary(summary_path)
    except PreparationError as exc:
        if "No protonable residues" not in str(exc):
            raise
        named.write_text(source.read_text(encoding="utf-8"), encoding="utf-8")
    else:
        manager.apply_protonation_states(str(source), str(named), ph=PH)
    # Same order as Preparation → Prepare. Bridged cysteines become CYX
    # before hydrogens are rebuilt, so reduce does not add a thiol hydrogen
    # that the CYX template cannot accept.
    remap_pdb_peptide_resnames(str(named), str(named))
    n_bridges = manager.apply_disulfide_bonds(str(named), str(named))
    completed = folder / "heavy_atoms_completed.pdb"
    complete_info = complete_missing_heavy_atoms(str(named), str(completed))
    amber_input = completed
    removed_h = 0
    if resolve_reduce_executable():
        stripped = folder / "protein_h_removed.pdb"
        strip_info = strip_protein_hydrogens(str(amber_input), str(stripped))
        removed_h = int(strip_info.get("removed") or 0)
        amber_input = stripped
    result = manager.run_pdb4amber_with_cap_fix(
        str(amber_input),
        str(dest),
        fix_caps=True,
        preserve_residue_numbers=True,
        pdb4amber_options={"reduce": True, "dry": not keep_water},
    )
    if not result.get("success", True) and not dest.is_file():
        raise RuntimeError(result.get("stderr") or "pdb4amber failed")
    protein_atoms, hetero_atoms = count_records(dest)
    return {
        "path": str(dest),
        "protein_atoms": protein_atoms,
        "hetero_atoms": hetero_atoms,
        "ph": PH,
        "disulfide_bonds": n_bridges,
        "heavy_atoms_added": int(complete_info.get("atoms_added") or 0),
        "protein_hydrogens_removed": removed_h,
    }


def write_structure_check(
    pdb_path: Path,
    folder: Path,
    seqres_from: Path | None = None,
    chains: set[str] | None = None,
) -> dict:
    """Write the same report as Visualize → Structure checker.

    The check before preparation keeps SEQRES lines for the selected chains,
    so unmodeled residues are listed. The check after preparation does not:
    protonation renames residues (ASH, HIE, CYX) and those names would look
    like sequence differences. A finding does not stop the later stages.
    """
    folder.mkdir(parents=True, exist_ok=True)
    target = pdb_path
    if seqres_from is not None and chains and seqres_from.is_file():
        header = [
            line if line.endswith("\n") else line + "\n"
            for line in seqres_from.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
            if line.startswith("SEQRES") and len(line) >= 12 and line[11:12] in chains
        ]
        if header:
            target = folder / "with_seqres.pdb"
            target.write_text(
                "".join(header) + pdb_path.read_text(encoding="utf-8", errors="replace"),
                encoding="utf-8",
            )
    report = diagnose_structure(str(target))
    text = render_structure_check(pdb_path, report)
    (folder / "report.txt").write_text(text, encoding="utf-8")
    (folder / "report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    summary = report["summary"]
    print(
        "    "
        f"{summary['missing_atoms']} missing atoms, "
        f"{summary['missing_residues']} missing residues, "
        f"{summary['gaps']} breaks, "
        f"{summary['altlocs']} alternate conformations, "
        f"{summary['chain_differences']} chain differences, "
        f"{summary['other']} other"
    )
    print(f"    {folder / 'report.txt'}")
    return {"path": str(pdb_path), "report": str(folder / "report.txt"), **summary}


def _place(chain: str, resid: object, icode: str = "", resname: str = "") -> str:
    return f"{chain or '?'}{resid}{icode or ''}{f' {resname}' if resname else ''}"


def render_structure_check(pdb_path: Path, report: dict) -> str:
    """Plain-text form of the Structure checker panel."""
    summary = report["summary"]
    lines = [
        f"Structure checker: {pdb_path}",
        (
            f"{summary['missing_atoms']} missing atoms · "
            f"{summary['missing_residues']} missing residues · "
            f"{summary['gaps']} breaks · "
            f"{summary['altlocs']} alternate conformations · "
            f"{summary['chain_differences']} chain differences · "
            f"{summary['other']} other"
        ),
        "",
    ]
    quiet = not any(summary.get(key) for key in (
        "missing_atoms",
        "missing_residues",
        "gaps",
        "altlocs",
        "chain_differences",
        "other",
    ))
    if quiet:
        lines.append(
            "Nothing to report. The loaded model has no missing heavy atoms, "
            "gaps, alternate conformations, or chain differences."
        )
        return "\n".join(lines) + "\n"

    missing_atoms = report.get("missing_atoms") or []
    if missing_atoms:
        lines.append(f"Missing atoms ({len(missing_atoms)})")
        for row in missing_atoms:
            lines.append(
                f"  {_place(row['chain'], row['resid'], row.get('icode') or '', row['resname'])} "
                f"is missing {row['atom']}"
            )
        lines.append("")

    missing_residues = report.get("missing_residues") or []
    if missing_residues:
        lines.append(f"Missing residues ({len(missing_residues)})")
        for row in missing_residues:
            if row.get("source") == "seqres":
                names = " ".join(row.get("resnames") or [])
                extra = f": {names}" if names else ""
                lines.append(
                    f"  Chain {row.get('chain') or '?'} SEQRES has {row['count']} "
                    f"names without coordinates{extra}."
                )
            else:
                plural = "residue" if row["count"] == 1 else "residues"
                listed = ", ".join(str(item) for item in (row.get("missing") or []))
                extra = f" ({listed})" if listed else ""
                lines.append(
                    f"  Chain {row['chain']} is missing {row['count']} {plural} "
                    f"between {row['after_resid']} and {row['before_resid']}{extra}."
                )
        lines.append("")

    gaps = report.get("gaps") or []
    if gaps:
        lines.append(f"Backbone breaks ({len(gaps)})")
        for row in gaps:
            lines.append(
                f"  Chain {row['chain']}: {row['resname_a']} {row['resid_a']} to "
                f"{row['resname_b']} {row['resid_b']} is {row['distance']} Å"
            )
        lines.append("")

    altlocs = report.get("altlocs") or []
    if altlocs:
        lines.append(f"Alternate conformations ({len(altlocs)})")
        for group in altlocs:
            lines.append(
                f"  {_place(group['chain'], group['resid'], group.get('icode') or '', group['resname'])}"
            )
            for option in group.get("options") or []:
                likely = " · most likely" if option.get("likely") else ""
                lines.append(f"    {option['altloc']} {float(option['occupancy']):.2f}{likely}")
        lines.append("")

    pairs = report.get("chain_pairs") or []
    if pairs:
        lines.append(f"Chains ({len(pairs)})")
        for pair in pairs:
            lines.append(
                f"  Chain {pair['chain_a']} ({pair['length_a']}) and chain "
                f"{pair['chain_b']} ({pair['length_b']})"
            )
            note = "These chains match." if pair.get("same") else pair.get("length_note") or ""
            lines.append(f"    {note} Identity {round(float(pair['identity']) * 100)}%.")
            for sub in pair.get("substitutions") or []:
                lines.append(
                    f"    {pair['chain_a']}{sub['resid_a']} {sub['resname_a']} differs from "
                    f"{pair['chain_b']}{sub['resid_b']} {sub['resname_b']}"
                )
            for diff in pair.get("length_differences") or []:
                lines.append(
                    f"    {diff['count_a']} residues on {diff['chain_a']} line up with "
                    f"{diff['count_b']} on {diff['chain_b']}."
                )
        lines.append("")

    issues = report.get("issues") or []
    if issues:
        lines.append(f"Other ({len(issues)})")
        for issue in issues:
            lines.append(f"  {issue.get('message') or ''}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def ligand_charge_choice(spec: dict, name: str, entry) -> tuple[int, str, str]:
    """Net charge, GAFF atom type, and antechamber charge method for one residue.

    ``entry`` is an integer charge, or a dict with ``charge`` and optional
    ``atom_type`` (``gaff2`` or ``gaff``) and ``charge_method`` (``abcg2``
    or ``bcc``). Missing fields use the system entry, then the script defaults.
    """
    atom_type = str(spec.get("atom_type") or LIGAND_ATOM_TYPE)
    charge_method = str(spec.get("charge_method") or LIGAND_CHARGE_METHOD)
    if isinstance(entry, dict):
        if "charge" not in entry:
            raise RuntimeError(f"{spec['pdb_id']} {name} parametrize dict needs charge.")
        atom_type = str(entry.get("atom_type") or atom_type)
        charge_method = str(entry.get("charge_method") or charge_method)
        charge = int(entry["charge"])
    else:
        charge = int(entry)
    if atom_type not in ("gaff2", "gaff"):
        raise RuntimeError(
            f"{spec['pdb_id']} {name} atom_type must be gaff2 or gaff, not {atom_type}."
        )
    if charge_method not in ("abcg2", "bcc"):
        raise RuntimeError(
            f"{spec['pdb_id']} {name} charge_method must be abcg2 or bcc, not {charge_method}."
        )
    return charge, atom_type, charge_method


def parametrize_kept_ligands(pdb_path: Path, spec: dict, folder: Path) -> dict:
    ligand_params = {}
    for name, library in (spec.get("manual_params") or {}).items():
        ligand_params[name] = dict(library)
    charges = spec.get("parametrize") or {}
    for name, entry in charges.items():
        charge, atom_type, charge_method = ligand_charge_choice(spec, name, entry)
        warning = forcefield_charge_warning(atom_type, charge_method)
        if warning:
            print(f"    {warning}")
        produced = parametrize_ligand_from_system_pdb(
            str(pdb_path),
            name,
            str(folder),
            charge=charge,
            charge_method=charge_method,
            atom_type=atom_type,
        )
        ligand_params[name] = {
            "frcmod": produced["frcmod"],
            "lib": produced["lib"],
            "atom_type": atom_type,
            "charge_method": charge_method,
        }
    return {"ligand_params": ligand_params, "n_ligands": len(ligand_params)}


# Same map as gatewizard.utils.ions. Kept here so a copied script still
# fixes the ions when the installed package is older than this checkout.
# Amber atomic_ions.lib: FE2 atom is FE2; FE3 is not a unit (Fe3+ is FE/FE).
_ION_ATOM_BY_RESNAME = {
    "FE2": "FE2",
    "CU1": "CU",
    "GD3": "GD",
    "IOD": "I",
    "NH4": "N",
    "U4+": "U",
    "HE+": "H",
    "HZ+": "H",
}
_ION_RESIDUE_ALIAS = {"FE3": "FE"}


def align_deposited_ion_names(pdb_path: Path) -> None:
    """Match deposited monoatomic ions to Amber's ion-library atom names.

    Residue FE2 keeps that name and its atom becomes FE2. Residue FE3
    (Fe3+) becomes residue FE, atom FE. CA, CL, MG, and K are already
    the names Amber expects and are left as deposited.
    """
    lines = pdb_path.read_text(encoding="utf-8", errors="replace").splitlines(keepends=True)
    changed = False
    out = []
    for line in lines:
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 20:
            resname = line[17:20].strip().upper()
            target = _ION_RESIDUE_ALIAS.get(resname, resname)
            expected = _ION_ATOM_BY_RESNAME.get(target)
            if target != resname and expected is None:
                expected = target
            current = line[12:16].strip().upper()
            if expected and (current != expected or target != resname):
                atom = f"{expected:<4}"[:4]
                residue = f"{target:<3}"[:3]
                line = f"{line[:12]}{atom}{line[16:17]}{residue}{line[20:]}"
                changed = True
        out.append(line if line.endswith("\n") else line + "\n")
    if changed:
        pdb_path.write_text("".join(out), encoding="utf-8")


def resolved_protocol(spec: dict) -> dict:
    """PROTOCOL with this system's ``protocol`` dict laid on top."""
    merged = dict(PROTOCOL)
    merged.update(spec.get("protocol") or {})
    return merged


def protocol_stages(protocol: dict):
    """Minimization, equilibration, and production selected for this system."""
    production_ns = float(protocol.get("production_ns") or 0.0)
    stages = _build_universal_membrane_stages(
        str(protocol.get("ensemble") or "NPT"),
        float(protocol.get("temperature") or 310.15),
        include_minimization=bool(protocol.get("minimization")),
        include_production=production_ns > 0,
        mini_steps=int(protocol.get("mini_steps") or 10000),
        production_ns=production_ns or 1.0,
    )
    if not protocol.get("equilibration"):
        mini = [stage for stage in stages if stage.stage_kind == "minimization"]
        production = [stage for stage in stages if stage.stage_kind == "production"]
        engine = str(protocol.get("engine") or "").lower()
        if mini and engine in {"openmm", "namd"}:
            # Those engines fold minimization into the next MD stage.
            holder = next(stage for stage in stages if stage.name == "Equilibration 1")
            stages = mini + [replace(holder, time_ns=0.0, steps=0)] + production
        else:
            stages = mini + production
    overrides = protocol.get("stage_ns") or {}
    updated = []
    for stage in stages:
        if stage.name in overrides:
            ns = float(overrides[stage.name])
            steps = int(round(ns * 1_000_000 / stage.timestep)) if ns > 0 else 0
            stage = replace(stage, time_ns=ns, steps=steps)
        updated.append(stage)
    if not updated:
        raise RuntimeError("protocol is enabled but no stages were selected")
    return updated


def resolve_engine_executable(protocol: dict) -> str:
    override = str(protocol.get("engine_executable") or "").strip()
    if override:
        return override
    engine = str(protocol.get("engine") or "openmm").lower()
    if engine == "openmm":
        return sys.executable
    name = {"amber": "sander", "gromacs": "gmx", "namd": "namd3"}.get(engine)
    if name is None:
        raise RuntimeError("engine must be amber, gromacs, openmm, or namd")
    found = shutil.which(name)
    if not found:
        raise RuntimeError(
            f"{engine} executable {name!r} was not found on PATH. "
            "Set engine_executable in this system's protocol."
        )
    return found


def collect_job_files(job_dir: Path) -> dict:
    prmtop = next(job_dir.glob("*.prmtop"), None)
    inpcrd = next(job_dir.glob("*.inpcrd"), None)
    if prmtop is None or inpcrd is None:
        raise RuntimeError(f"No system.prmtop / system.inpcrd in {job_dir}")
    files = {"prmtop": str(prmtop), "inpcrd": str(inpcrd)}
    pdb = job_dir / "system.pdb"
    if pdb.is_file():
        files["pdb"] = str(pdb)
    bilayer = None
    candidates = [
        *job_dir.glob("bilayer*_lipid.pdb"),
        *(
            path
            for path in job_dir.glob("bilayer_*.pdb")
            if path.name != "bilayer_for_pdb4amber.pdb"
        ),
    ]
    for path in candidates:
        if read_cryst1_box(path) is not None:
            bilayer = path
            break
    if bilayer is None:
        for path in candidates:
            bilayer = path
            break
    if bilayer is not None:
        files["bilayer_pdb"] = str(bilayer)
    return files


def find_built_job(job_root: Path) -> Path:
    if (job_root / "system.prmtop").is_file():
        return job_root
    matches = sorted(path.parent for path in job_root.glob("02_build*/system.prmtop"))
    if not matches:
        raise RuntimeError(f"No parametrized build under {job_root}")
    return matches[-1]


def run_protocol(job_root: Path, protocol: dict) -> dict:
    """Write and run the GateWizard membrane protocol for one build."""
    job_dir = find_built_job(job_root)
    files = collect_job_files(job_dir)
    stages = [stage.to_dict() for stage in protocol_stages(protocol)]
    engine = str(protocol.get("engine") or "openmm").lower()
    executable = resolve_engine_executable(protocol)
    scheme = str(protocol.get("ensemble") or "NPT")
    output_name = "protocol"
    if engine == "openmm":
        manager = OpenMMEquilibrationManager(job_dir)
        result = manager.setup_openmm_equilibration(
            system_files=files,
            stage_params_list=stages,
            output_name=output_name,
            scheme_type=scheme,
            add_com_restraint=False,
            add_rotation_restraint=False,
        )
    elif engine == "amber":
        manager = AmberEquilibrationManager(job_dir, amber_executable=executable)
        result = manager.setup_amber_equilibration(
            system_files=files,
            stage_params_list=stages,
            output_name=output_name,
            scheme_type=scheme,
            amber_executable=executable,
        )
    elif engine == "gromacs":
        manager = GROMACSEquilibrationManager(job_dir)
        result = manager.setup_gromacs_equilibration(
            system_files=files,
            stage_params_list=stages,
            output_name=output_name,
            scheme_type=scheme,
            gmx_executable=executable,
            add_com_restraint=False,
            add_rotation_restraint=False,
        )
    elif engine == "namd":
        manager = NAMDEquilibrationManager(job_dir, executable)
        result = manager.setup_namd_equilibration(
            system_files=files,
            stage_params_list=stages,
            output_name=output_name,
            scheme_type=scheme,
            namd_executable=executable,
            add_com_restraint=False,
            add_rotation_restraint=False,
        )
    else:
        raise RuntimeError("engine must be amber, gromacs, openmm, or namd")
    script = Path(result.get("run_script") or (job_dir / output_name / "run_equilibration.sh"))
    if not script.is_file():
        raise RuntimeError(f"run_equilibration.sh was not written in {job_dir / output_name}")
    env = os.environ.copy()
    env["PYTHON"] = sys.executable
    platform = str(protocol.get("openmm_platform") or "").strip()
    if engine == "openmm" and platform:
        env["PLATFORM"] = platform
    print(f"  running {engine} protocol: {script}")
    completed = subprocess.run(
        ["bash", str(script)],
        cwd=str(script.parent),
        env=env,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"{script.name} exited {completed.returncode}")
    return {
        "engine": engine,
        "script": str(script),
        "stages": [stage["name"] for stage in stages],
        "production_ns": protocol.get("production_ns"),
    }


def build_once(pdb_path: Path, job_root: Path, ligand_params: dict, protein_atoms: int) -> dict:
    job_root.mkdir(parents=True, exist_ok=True)
    align_deposited_ion_names(pdb_path)
    builder = Builder()
    ok, message, job_dir = builder.prepare_system(
        pdb_file=str(pdb_path),
        working_dir=str(job_root),
        upper_lipids=["POPC"],
        lower_lipids=["POPC"],
        lipid_ratios="1//1",
        protein_ff="ff14SB",
        lipid_ff="lipid21",
        water_model="tip3p",
        salt_concentration=0.15,
        cation="K+",
        anion="Cl-",
        preoriented=True,
        parametrize=True,
        notprotonate=True,
        # reduce names the GLH alpha hydrogen HCA. The Builder checkbox
        # "Remove protein hydrogens" drops those names; tleap puts HA back.
        # Ligand hydrogens, including RET, stay.
        remove_protein_h=True,
        ligand_params=ligand_params,
        wait=True,
        wait_verbose=True,
    )
    if not ok or job_dir is None:
        raise RuntimeError(message)
    system_pdb = find_system_pdb(Path(job_dir))
    _protein, _hetero, system_atoms = count_records(system_pdb, total=True)
    saved = job_root / "system.pdb"
    if system_pdb.resolve() != saved.resolve():
        saved.write_bytes(system_pdb.read_bytes())
    return {
        "protein_atoms": protein_atoms,
        "system_atoms": system_atoms,
        "system_pdb": str(saved),
        "job_dir": str(job_dir),
        "message": message,
    }


def find_system_pdb(job_dir: Path) -> Path:
    direct = job_dir / "system.pdb"
    if direct.is_file():
        return direct
    matches = sorted(job_dir.rglob("system.pdb"))
    if not matches:
        raise RuntimeError(f"system.pdb was not written under {job_dir}")
    return matches[0]


def count_records(path: Path, total: bool = False):
    protein = 0
    hetero = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM"):
                protein += 1
            elif line.startswith("HETATM"):
                hetero += 1
    if total:
        return protein, hetero, protein + hetero
    return protein, hetero


def merge_saved_rows(rows: list[dict]) -> list[dict]:
    """Keep timings for PDBs that this run did not touch."""
    path = RESULTS / "timings.json"
    if not path.is_file():
        return rows
    try:
        previous = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return rows
    if not isinstance(previous, list):
        return rows
    touched = {row["pdb_id"] for row in rows}
    order = []
    for row in previous + rows:
        pdb_id = row.get("pdb_id")
        if pdb_id and pdb_id not in order:
            order.append(pdb_id)
    merged = []
    for pdb_id in order:
        source = rows if pdb_id in touched else previous
        merged.extend(row for row in source if row.get("pdb_id") == pdb_id)
    return merged


def persist_rows(rows: list[dict], replicas: int) -> None:
    """Write the timing table after each finished stage or replica."""
    write_report(merge_saved_rows(rows), replicas)


BUILDER_STEPS = ("packmol", "pdb4amber", "tleap")


def format_seconds(value: float) -> str:
    """Keep a short stage visible. One decimal hides anything under 0.05 s."""
    if value < 1:
        return f"{value:.3f}"
    if value < 10:
        return f"{value:.2f}"
    return f"{value:.1f}"


def builder_job_dir(row: dict) -> Path | None:
    """Directory that holds status.json for one builder replica."""
    stored = row.get("job_dir")
    if stored and Path(stored).is_dir():
        return Path(stored)
    replica = row.get("replica")
    pdb_id = row.get("pdb_id")
    if replica and pdb_id:
        local = WORK / str(pdb_id) / f"builder_rep{replica}" / "02_build_prepared"
        if local.is_dir():
            return local
    return None


def _file_time(folder: Path, name: str) -> datetime | None:
    path = folder / name
    if not path.is_file():
        return None
    return datetime.fromtimestamp(path.stat().st_mtime)


def builder_step_rows(row: dict) -> list[dict]:
    """Split one finished builder replica into packmol, pdb4amber, and tleap.

    The shell script records a start time and an end time. The files it
    writes mark the boundary of each step: packmol.log, the PDB handed to
    tleap, and system.prmtop. A replica that is still packing has no end
    time and is left as the single builder row.
    """
    if row.get("stage") != "builder" or not row.get("ok"):
        return []
    folder = builder_job_dir(row)
    if folder is None:
        return []
    status_path = folder / "status.json"
    if not status_path.is_file():
        return []
    try:
        status = json.loads(status_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return []
    if status.get("status") != "completed" or not status.get("start_time") or not status.get("end_time"):
        return []
    start = datetime.fromisoformat(status["start_time"])
    pack_end = _file_time(folder, "packmol.log") or _file_time(folder, "bilayer_prepared_lipid.pdb")
    amber_end = _file_time(folder, "system_for_tleap.pdb") or _file_time(folder, "bilayer_for_pdb4amber.pdb")
    tleap_end = _file_time(folder, "system.prmtop")
    bounds = [
        ("packmol", start, pack_end),
        ("pdb4amber", pack_end, amber_end),
        ("tleap", amber_end, tleap_end),
    ]
    steps = []
    for stage, began, ended in bounds:
        if began is None or ended is None:
            continue
        seconds = (ended - began).total_seconds()
        if seconds < 0:
            continue
        steps.append({
            "pdb_id": row["pdb_id"],
            "stage": stage,
            "replica": row.get("replica"),
            "seconds": seconds,
            "ok": True,
            "error": "",
        })
    return steps


def attach_builder_steps(rows: list[dict]) -> list[dict]:
    """Add packmol, pdb4amber, and tleap after each builder row."""
    stored = {
        (row.get("pdb_id"), row.get("replica"), row.get("stage")): row
        for row in rows
        if row.get("stage") in BUILDER_STEPS
    }
    attached = []
    for row in rows:
        if row.get("stage") in BUILDER_STEPS:
            continue
        attached.append(row)
        if row.get("stage") != "builder":
            continue
        fresh = builder_step_rows(row)
        if fresh:
            attached.extend(fresh)
            continue
        for stage in BUILDER_STEPS:
            saved = stored.get((row.get("pdb_id"), row.get("replica"), stage))
            if saved:
                attached.append(saved)
    return attached


def write_report(rows: list[dict], replicas: int) -> None:
    rows = attach_builder_steps(rows)
    RESULTS.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pdb_id",
        "stage",
        "replica",
        "seconds",
        "ok",
        "protein_atoms",
        "hetero_atoms",
        "system_atoms",
        "error",
    ]
    with (RESULTS / "timings.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    (RESULTS / "timings.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
    (RESULTS / "summary.md").write_text(render_summary(rows, replicas), encoding="utf-8")


def render_summary(rows: list[dict], replicas: int) -> str:
    lines = [
        "# System preparation timings",
        "",
        f"Builder replicas requested: {replicas}.",
        "Mean and sample standard deviation use the successful replicas only.",
        "",
        "| PDB | Protein atoms | System atoms (each replica) | Builder seconds | Mean | Std |",
        "| --- | ---: | --- | --- | ---: | ---: |",
    ]
    pdb_ids = []
    for row in rows:
        if row["pdb_id"] not in pdb_ids:
            pdb_ids.append(row["pdb_id"])
    for pdb_id in pdb_ids:
        system_rows = [row for row in rows if row["pdb_id"] == pdb_id]
        builds = [row for row in system_rows if row["stage"] == "builder"]
        ok_builds = [row for row in builds if row.get("ok")]
        protein = ""
        for row in system_rows:
            if row.get("protein_atoms"):
                protein = str(row["protein_atoms"])
        atom_list = [str(row.get("system_atoms") or "") for row in ok_builds]
        seconds = [float(row["seconds"]) for row in ok_builds]
        mean = f"{statistics.mean(seconds):.1f}" if seconds else ""
        std = f"{statistics.stdev(seconds):.1f}" if len(seconds) >= 2 else ""
        individual = ", ".join(f"{value:.1f}" for value in seconds)
        lines.append(
            f"| {pdb_id} | {protein} | {', '.join(atom_list)} | {individual} | {mean} | {std} |"
        )
    lines.extend(["", "## Every stage", ""])
    lines.append("| PDB | Stage | Replica | Seconds | OK |")
    lines.append("| --- | --- | ---: | ---: | --- |")
    for row in rows:
        replica = "" if row.get("replica") is None else str(row["replica"])
        lines.append(
            f"| {row['pdb_id']} | {row['stage']} | {replica} | {format_seconds(float(row['seconds']))} | {row['ok']} |"
        )
    checks = [
        row
        for row in rows
        if row.get("stage") in {"structure_check_before", "structure_check_after"} and row.get("ok")
    ]
    if checks:
        lines.extend(["", "## Structure check", ""])
        for row in checks:
            when = "before preparation" if row["stage"].endswith("before") else "after preparation"
            lines.append(
                f"- {row['pdb_id']} {when}: "
                f"{row.get('missing_atoms', 0)} missing atoms, "
                f"{row.get('missing_residues', 0)} missing residues, "
                f"{row.get('gaps', 0)} breaks, "
                f"{row.get('altlocs', 0)} alternate conformations, "
                f"{row.get('chain_differences', 0)} chain differences, "
                f"{row.get('other', 0)} other. "
                f"Full list: `{row.get('report') or ''}`."
            )
    failed = [row for row in rows if not row.get("ok")]
    if failed:
        lines.extend(["", "## Failures", ""])
        for row in failed:
            replica = "" if row.get("replica") is None else f" replica {row['replica']}"
            first = (row.get("error") or "").strip().splitlines()
            message = first[0] if first else "failed"
            lines.append(f"- {row['pdb_id']} {row['stage']}{replica}: {message}")
    lines.append("")
    return "\n".join(lines)


if __name__ == "__main__":
    main()
