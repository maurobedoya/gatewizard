from pathlib import Path
from gatewizard.core.mempro import MemPrO

pdb_file = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
mp = MemPrO()
if MemPrO.is_available():
    results = mp.run(
        pdb_file,
        output_dir="my_orient",
        n_cpus=4,
        n_iters=200,
        grid_size=72,
    )
    print(f"Found {len(results)} orientations")
