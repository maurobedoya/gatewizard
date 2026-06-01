import os
from gatewizard.core.mempro import MemPrO, MemProError

orient_dir = os.environ.get("MEMPRO_ORIENT_DIR", "Orient")
try:
    # Parse pre-computed results
    results = MemPrO.parse_results(orient_dir)

    # Load the best orientation into the viewer
    if results and results[0].pdb_path:
        from gatewizard.core.structure_manager import StructureManager

        viewer = StructureManager()
        info = viewer.load_structure(results[0].pdb_path)
        print(f"Loaded rank 1: {info['n_atoms']} atoms, {info['n_chains']} chains")
except (MemProError, FileNotFoundError):
    print(f"Orient directory not found at '{orient_dir}', skipping.")
