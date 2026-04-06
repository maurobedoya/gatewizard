from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run("protein.pdb", dual_membrane=True)
    print(f"Dual membrane: {len(results)} orientations")
