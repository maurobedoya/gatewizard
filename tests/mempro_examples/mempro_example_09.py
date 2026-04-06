from gatewizard.core.mempro import MemPrO

mp = MemPrO()
if MemPrO.is_available():
    cmd = mp.build_command("protein.pdb", n_cpus=4, dual_membrane=True)
    print("Command:", " ".join(cmd))
