from pathlib import Path
from gatewizard.core.mempro import MemPrO

pdb_file = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
mp = MemPrO()
if MemPrO.is_available():
    results = mp.run(pdb_file, dual_membrane=True)
    print(f"Dual membrane: {len(results)} orientations")
