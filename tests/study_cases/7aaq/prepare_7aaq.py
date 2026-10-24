#!/usr/bin/env python3
"""Preparation step for 7aaq_protein_oriented.pdb (same files as the GUI Prepare job).

Input (from orient_7aaq.py):
  7aaq_protein_oriented.pdb

Expected outputs in 01_preparation_7aaq_protein_oriented/:
  7aaq_protein_oriented.pdb
  7aaq_protein_oriented.pka
  7aaq_protein_oriented_summary_of_prediction.txt
  7aaq_protein_oriented_protonated.pdb
  7aaq_protein_oriented_protonated_nonprot.pdb
  7aaq_protein_oriented_protonated_renum.txt
  7aaq_protein_oriented_protonated_sslink

Run:
  conda activate gatewizard
  python prepare_7aaq.py
"""

from pathlib import Path
import shutil
import tempfile

from gatewizard.core.preparation import (
    PreparationManager,
    complete_missing_heavy_atoms,
    resolve_reduce_executable,
    strip_protein_hydrogens,
)
from gatewizard.utils.peptide_residues import remap_pdb_peptide_resnames

WORK = Path(__file__).resolve().parent
INPUT_PDB = WORK / "7aaq_protein_oriented.pdb"
OUTPUT_DIR = WORK / "01_preparation_7aaq_protein_oriented"
PH = 6.0

# Sequential Amber numbering, as this study case requests.
# The desktop app keeps original residue numbers unless that box is turned off.
PRESERVE_RESIDUE_NUMBERS = False


def main() -> None:
    if not INPUT_PDB.is_file():
        raise FileNotFoundError(
            f"Missing {INPUT_PDB}. Run orient_7aaq.py first."
        )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copy2(INPUT_PDB, OUTPUT_DIR / INPUT_PDB.name)

    manager = PreparationManager()
    pka_file = manager.run_analysis(str(INPUT_PDB), output_dir=str(OUTPUT_DIR))
    summary_file = manager.extract_summary(pka_file)
    residues = manager.parse_summary(summary_file)

    with tempfile.TemporaryDirectory(prefix="gw_prep_") as tmp_name:
        tmp_dir = Path(tmp_name)
        named = tmp_dir / "protonated_names.pdb"
        manager.apply_protonation_states(
            input_pdb=str(INPUT_PDB),
            output_pdb=str(named),
            ph=PH,
            residues=residues,
        )
        remap_pdb_peptide_resnames(str(named), str(named))
        manager.apply_disulfide_bonds(str(named), str(named))

        completed = tmp_dir / "heavy_atoms_completed.pdb"
        complete_missing_heavy_atoms(str(named), str(completed))
        amber_input = completed
        if resolve_reduce_executable():
            stripped = tmp_dir / "protein_h_removed.pdb"
            strip_protein_hydrogens(str(amber_input), str(stripped))
            amber_input = stripped

        protonated = OUTPUT_DIR / "7aaq_protein_oriented_protonated.pdb"
        result = manager.run_pdb4amber_with_cap_fix(
            input_pdb=str(amber_input),
            output_pdb=str(protonated),
            fix_caps=True,
            preserve_residue_numbers=PRESERVE_RESIDUE_NUMBERS,
            pdb4amber_options={"reduce": True, "dry": True},
        )

    if not result.get("success", True):
        raise RuntimeError(result.get("stderr") or "pdb4amber failed")

    print(f"pKa:      {pka_file}")
    print(f"summary:  {summary_file}")
    print(f"prepared: {protonated}")
    print(f"reduce:   used={result.get('reduce_used')} skipped={result.get('reduce_skipped')}")
    for name in (
        "7aaq_protein_oriented.pdb",
        "7aaq_protein_oriented.pka",
        "7aaq_protein_oriented_summary_of_prediction.txt",
        "7aaq_protein_oriented_protonated.pdb",
        "7aaq_protein_oriented_protonated_nonprot.pdb",
        "7aaq_protein_oriented_protonated_renum.txt",
        "7aaq_protein_oriented_protonated_sslink",
    ):
        path = OUTPUT_DIR / name
        print(f"  {'OK' if path.is_file() else 'MISSING'}  {path.name}")


if __name__ == "__main__":
    main()
