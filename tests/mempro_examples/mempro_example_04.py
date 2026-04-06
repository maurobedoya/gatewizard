from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    results = mp.run(
        "protein.pdb",
        output_dir="my_orient",
        n_cpus=4,
        n_iters=200,
        grid_size=72,
    )
    print(f"Found {len(results)} orientations")
