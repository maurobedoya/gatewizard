import shutil
import sys

if shutil.which("propka3") is None:
    sys.exit("Skip: PropKa 3 is not on PATH")

from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
