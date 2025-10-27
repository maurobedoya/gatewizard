from gatewizard.core.preparation import PreparationManager

analyzer = PreparationManager()
pka_file = analyzer.run_analysis("protein.pdb")
summary_file = analyzer.extract_summary(pka_file)
residues = analyzer.parse_summary(summary_file)

for res in residues:
    states = [analyzer.get_default_protonation_state(res, ph=p) for p in [2.0, 7.0, 11.0]]
    print(f"{res['residue']}{res['res_id']} (pKa={res['pka']:.2f}): "
          f"pH2={states[0]}, pH7={states[1]}, pH11={states[2]}")
