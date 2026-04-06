from gatewizard.core.mempro import MemPrO

if MemPrO.is_available():
    print("MemPrO is installed and ready to use")
else:
    print("Install MemPrO: pip install git+https://github.com/pstansfeld/MemPrO.git")
