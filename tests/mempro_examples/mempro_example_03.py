from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run("protein.pdb")
    for r in results:
        print(
            f"Rank {r.rank}: potential={r.relative_potential:.2f}, hits={r.hits_pct:.1f}%"
        )
        print(f"  PDB: {r.pdb_path}")
