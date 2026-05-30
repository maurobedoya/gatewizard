from pathlib import Path
from gatewizard.core.mempro import MemPrO

pdb_file = str(Path(__file__).parent.parent / "6RV3_AB.pdb")
mp = MemPrO()
if MemPrO.is_available():
    results = mp.run(pdb_file, peripheral=True)
    print(f"Peripheral: {len(results)} orientations")
