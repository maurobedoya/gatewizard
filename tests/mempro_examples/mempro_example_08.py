from gatewizard.core.mempro import MemPrO

try:
    pdb = MemPrO.get_oriented_pdb("Orient", rank=1)
    print(f"Best orientation: {pdb}")
except FileNotFoundError as e:
    print(f"Not found: {e}")
