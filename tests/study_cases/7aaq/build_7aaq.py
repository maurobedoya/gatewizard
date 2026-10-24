#!/usr/bin/env python3
"""Builder step for 7aaq_protein_oriented_protonated.pdb (same files as the GUI Build job).

Input (from prepare_7aaq.py):
  01_preparation_7aaq_protein_oriented/7aaq_protein_oriented_protonated.pdb

Output folder (created automatically from the PDB stem):
  02_build_7aaq_protein_oriented_protonated/

Settings match the study case:
  POPC//POPC, ff14SB, lipid21, TIP3P, 0.15 M KCl, dist=12, dist_wat=26,
  preoriented=True, parametrize=True, notprotonate=True.
  MemPrO already placed the protein, so the build skips MEMEMBED.
  A folder that already has a build is refused unless OVERWRITE is True.

Run:
  conda activate gatewizard
  python build_7aaq.py
"""

from pathlib import Path

from gatewizard.core.builder import Builder

WORK = Path(__file__).resolve().parent

# A folder that already has a build is refused. Delete it, or set this True to reuse it.
OVERWRITE = False

PROTONATED_PDB = WORK / "01_preparation_7aaq_protein_oriented" / "7aaq_protein_oriented_protonated.pdb"
JOB_NAME = "02_build_7aaq_protein_oriented_protonated"

EXPECTED_FILES = (
    "7aaq_protein_oriented_protonated.pdb",
    "bilayer_7aaq_protein_oriented_protonated.pdb",
    "bilayer_7aaq_protein_oriented_protonated.pdb_FORCED",
    "bilayer_7aaq_protein_oriented_protonated_lipid.crd",
    "bilayer_7aaq_protein_oriented_protonated_lipid.pdb",
    "bilayer_7aaq_protein_oriented_protonated_lipid.top",
    "bilayer_for_pdb4amber.pdb",
    "leap.log",
    "leap_.log",
    "leap_parametrize.in",
    "packmol-memgen.json",
    "packmol-memgen.log",
    "packmol.inp",
    "packmol.log",
    "process.pid",
    "run_preparation.sh",
    "status.json",
    "status_utils.py",
    "system.inpcrd",
    "system.pdb",
    "system.prmtop",
    "system_for_tleap.pdb",
    "system_for_tleap_nonprot.pdb",
    "system_for_tleap_renum.txt",
    "system_for_tleap_sslink",
    "tleap_geometry.in",
)


def main() -> None:
    if not PROTONATED_PDB.is_file():
        raise FileNotFoundError(
            f"Missing {PROTONATED_PDB}. Run orient_7aaq.py, then prepare_7aaq.py."
        )

    builder = Builder()
    builder.set_configuration(
        water_model="tip3p",
        protein_ff="ff14SB",
        lipid_ff="lipid21",
        preoriented=True,
        parametrize=True,
        notprotonate=True,
        salt_concentration=0.15,
        cation="K+",
        anion="Cl-",
        dist=12,
        dist_wat=26,
        nloop=20,
        nloop_all=100,
        tolerance=2.0,
    )

    ok, msg, job_dir = builder.prepare_system(
        pdb_file=str(PROTONATED_PDB),
        working_dir=str(WORK),
        upper_lipids=["POPC"],
        lower_lipids=["POPC"],
        lipid_ratios="1//1",
        output_folder_name=JOB_NAME,
        overwrite=OVERWRITE,
        wait=True,
        wait_verbose=True,
    )
    print(ok, msg, job_dir)
    if not ok or job_dir is None:
        raise SystemExit(1)

    print("Expected outputs:")
    for name in EXPECTED_FILES:
        path = Path(job_dir) / name
        print(f"  {'OK' if path.exists() else 'MISSING'}  {name}")
    logs_dir = Path(job_dir) / "logs"
    print(f"  {'OK' if logs_dir.is_dir() else 'MISSING'}  logs/")


if __name__ == "__main__":
    main()
