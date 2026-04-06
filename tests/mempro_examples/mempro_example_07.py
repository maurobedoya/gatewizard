import os
from gatewizard.core.mempro import MemPrO, MemProError

orient_dir = os.environ.get("MEMPRO_ORIENT_DIR", "Orient")
try:
    results = MemPrO.parse_results(orient_dir)
    for r in results:
        print(
            f"Rank {r.rank}: potential={r.relative_potential:.2f}, "
            f"hits={r.hits_pct:.1f}%, pdb={r.pdb_path}"
        )
except (MemProError, FileNotFoundError):
    print(f"Orient directory not found at '{orient_dir}', skipping.")
