import shutil
import sys

if shutil.which("propka3") is None:
    sys.exit("Skip: PropKa 3 is not on PATH")

from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)

curves = analyzer.get_ph_titration_curve(ph_range=(4, 10), ph_step=1.0)

for residue_id, curve in curves.items():
    print(f"\n{residue_id}:")
    for ph, state in curve:
        print(f"  pH {ph:.1f}: {state}")
