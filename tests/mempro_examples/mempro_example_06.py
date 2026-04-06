from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run("protein.pdb", peripheral=True)
    print(f"Peripheral: {len(results)} orientations")
