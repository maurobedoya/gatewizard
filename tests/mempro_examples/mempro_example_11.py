from pathlib import Path
from gatewizard.core.mempro import MemPrO, MemProError

pdb_file = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
mp = MemPrO()

# Check availability
if not MemPrO.is_available():
    print("Install: pip install git+https://github.com/pstansfeld/MemPrO.git")
else:
    # Run orientation
    try:
        results = mp.run(pdb_file, n_cpus=4)
        print(f"Found {len(results)} orientations\n")

        # Display results table
        print(f"{'Rank':>4}  {'Potential':>10}  {'Hits%':>6}  {'Re-rank':>10}")
        print("-" * 36)
        for r in results:
            print(
                f"{r.rank:>4d}  {r.relative_potential:>10.3f}  "
                f"{r.hits_pct:>5.1f}%  {r.rerank_potential:>10.3f}"
            )

        # Load the best orientation
        if results and results[0].pdb_path:
            print(f"\nBest PDB: {results[0].pdb_path}")
    except MemProError as e:
        print(f"Error: {e}")
