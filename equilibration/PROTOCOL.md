# Membrane equilibration protocol

Same schedule for NAMD, GROMACS, OpenMM, and Amber. Still testing — not a finished production protocol.

The membrane is packed under NPgT with γ=0. The ensemble you pick (NVT, NPT, NPAT, or NPgT) is used only in production. Eq1–6 do not change with that choice.

The GUI default is `gatewizard-gui/resources/protocols/base.json`. The API/CLI default is `_build_universal_membrane_stages()` in `gatewizard/tools/equilibration.py`. If you edit one, update the other.

## Stages

Minimization and production are extra. The six MD stages below add up to 20 ns. There is no Eq7.

| Stage | Ensemble | dt | Time | Save every | Restraints (kcal/mol/Å²) |
|-------|----------|-----|------|------------|--------------------------|
| Mini | NVT | — | 10k steps | — | BB 10, SC 5, lipid 2.5, ions 5 |
| Eq1 heat | NVT | 1 fs | 0.125 ns | 5000 | BB 10, SC 5, lipid 2.5 |
| Eq2 scaffold | NVT | 1 fs | 0.5 ns | 5000 | BB 10, SC 5, lipid 5 |
| Eq3 pack | NPgT γ=0 | 1 fs | 0.25 ns | 5000 | BB 5, SC 2.5, heads 2.5, tails free |
| Eq4 pack | NPgT | 1 fs | 0.5 ns | 5000 | BB 2.5, SC 1, lipids free |
| Eq5 pack | NPgT | 2 fs | 1.0 ns | 5000 | BB 1.0 |
| Eq6 pack | NPgT γ=0 | 2 fs | 17.625 ns | 50000 | BB 0.1 |
| Production | chosen ensemble | 2 fs | you set | 50000 | none |

Barostat starts at Eq3. Backbone k goes 2.5 → 1.0 → 0.1, then 0 in production.

Trajectory files follow the engine (DCD, XTC, or NetCDF), not the ensemble. Eq6 writes as often as production (every 50 000 steps, 100 ps/frame at 2 fs).

## Files

Each engine has one `eq/` folder (mini + Eq1–6) and four production templates:

```
equilibration/{engine}/
  eq/
  production/NVT|NPT|NPAT|NPgT/
```

Eq stages always load `eq/`. Production loads `production/{ensemble}/`. Filenames stay CHARMM-GUI style (`step6.1_equilibration.inp`, `step7_production.mdp`, …).

## Notes

- NVT and NPAT production keep the packed XY box. There is no extra step that eases into that ensemble after Eq6, so the barostat can jump.
- Judge thickness and area per lipid on late production, not the whole trajectory. Packing frames make APL look too high (and the membrane too thin) because the box is still wide.
- packmol-memgen POPC often starts around 90 Å²/lipid. Long NPT with Lipid21 or CHARMM36 is often in the low-to-mid 60s. Getting to ~65 Å² is not guaranteed just by packing longer.
- Amber Eq3 (first barostat) should run on CPU `pmemd`. GPU PME often dies with “box dimensions changed too much” while the box first collapses. Later stages, including Eq6, can use `pmemd.cuda`.
- For Amber, set `IFBOX=1` if any stage uses a barostat — even if production is NVT. Trajectory write interval is `ntwx` from `dcd_freq`; format is NetCDF (`ioutfm=1`).
- Eq6 is long (~17 ns). Production length is separate.
