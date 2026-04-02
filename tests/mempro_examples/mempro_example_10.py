import os
from gatewizard.core.mempro import MemPrO, MemProError

orient_dir = os.environ.get("MEMPRO_ORIENT_DIR", "Orient")
try:
    results = MemPrO.parse_results(orient_dir)
    if results:
        best = results[0]
        print(f"Rank: {best.rank}")
        print(f"Relative potential: {best.relative_potential}")
        print(f"Hits: {best.hits_pct}%")
        print(f"Re-rank potential: {best.rerank_potential}")
        print(f"Re-rank depth: {best.rerank_depth}")
        print(f"Re-rank value: {best.rerank_value}")
        print(f"PDB path: {best.pdb_path}")
except (MemProError, FileNotFoundError):
    print(f"Orient directory not found at '{orient_dir}', skipping.")
