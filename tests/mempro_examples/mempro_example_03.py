from pathlib import Path
from gatewizard.core.mempro import MemPrO

pdb_file = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
mp = MemPrO()
if MemPrO.is_available():
    results = mp.run(pdb_file)
    for r in results:
        print(
            f"Rank {r.rank}: potential={r.relative_potential:.2f}, hits={r.hits_pct:.1f}%"
        )
        print(f"  PDB: {r.pdb_path}")
