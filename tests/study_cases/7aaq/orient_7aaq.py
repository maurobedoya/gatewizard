#!/usr/bin/env python3
"""Orient 7aaq_protein.pdb with MemPrO before preparation.

The desktop Apply button uses the same two calls. MemPrO ranks poses of the
protein. StructureManager.apply_mempro_orientation copies that rigid move onto
every atom, so BGC moves with the protein.

Output:
  7aaq_protein_oriented.pdb

Run:
  conda activate gatewizard
  python orient_7aaq.py
"""

from pathlib import Path

from gatewizard.core.mempro import MemPrO
from gatewizard.core.structure_manager import StructureManager

WORK = Path(__file__).resolve().parent
INPUT_PDB = WORK / "7aaq_protein.pdb"
OUTPUT_PDB = WORK / "7aaq_protein_oriented.pdb"
OUTPUT_DIR = WORK / "mempro_7aaq_protein"
# Same seed as tests/study_cases/prepare_systems.py. The desktop dialog can
# leave Seed empty. This script sets 7 so the pose repeats.
SEED = 7


def main() -> None:
    if not INPUT_PDB.is_file():
        raise FileNotFoundError(INPUT_PDB)
    if not MemPrO.is_available():
        raise RuntimeError(
            'mempro is not installed. Install with: pip install "gatewizard[full]"'
        )

    results = MemPrO().run(
        str(INPUT_PDB),
        output_dir=str(OUTPUT_DIR),
        n_iters=150,
        grid_size=36,
        device="cpu",
        seed=SEED,
    )
    if not results:
        raise RuntimeError("MemPrO returned no orientations.")

    manager = StructureManager()
    manager.load_structure(str(INPUT_PDB))
    n_atoms = manager.apply_mempro_orientation(results[0].pdb_path)
    manager.structure.write_pdb(str(OUTPUT_PDB))
    print(f"rank:     {results[0].rank}")
    print(f"atoms:    {n_atoms}")
    print(f"oriented: {OUTPUT_PDB}")


if __name__ == "__main__":
    main()
