# STP10 with the API

This walkthrough prepares *Arabidopsis thaliana* STP10 (PDB 7AAQ) with the Python API: MemPro orientation, protonation, a pre-oriented POPC membrane, and OpenMM equilibration inputs for a 100 ns production run. The same choices in the desktop app are [STP10 in the desktop app](7aaq-gui.md).

The scripts and the cleaned protein are in [`tests/study_cases/7aaq/`](https://github.com/maurobedoya/gatewizard/tree/main/tests/study_cases/7aaq). Copy that folder to the directory where the outputs should be written. Each script uses the folder it lives in, not the shell’s current directory.

| Stage | Script | Main input | Main output |
| --- | --- | --- | --- |
| Orientation | [`orient_7aaq.py`](https://github.com/maurobedoya/gatewizard/blob/main/tests/study_cases/7aaq/orient_7aaq.py) | `7aaq_protein.pdb` | `7aaq_protein_oriented.pdb` |
| Preparation | [`prepare_7aaq.py`](https://github.com/maurobedoya/gatewizard/blob/main/tests/study_cases/7aaq/prepare_7aaq.py) | `7aaq_protein_oriented.pdb` | `7aaq_protein_oriented_protonated.pdb` |
| Builder | [`build_7aaq.py`](https://github.com/maurobedoya/gatewizard/blob/main/tests/study_cases/7aaq/build_7aaq.py) | `7aaq_protein_oriented_protonated.pdb` | `system.prmtop`, `system.inpcrd`, `system.pdb` |
| Equilibration | [`equilibrate_7aaq.py`](https://github.com/maurobedoya/gatewizard/blob/main/tests/study_cases/7aaq/equilibrate_7aaq.py) | the built system | `step*.inp`, `run_equilibration.sh` |

## Environment

Python, GateWizard, AmberTools, and OpenMM need to be installed first. The [installation page](../installation.md) covers the `gatewizard` conda environment. OpenMM is the engine that runs the equilibration.

```bash
conda activate gatewizard
cd /path/to/7aaq
python orient_7aaq.py
python prepare_7aaq.py
```

`7aaq_protein.pdb` is chain A of 7AAQ with the bound sugar BGC kept and the crystallization additives removed. The desktop walkthrough produces the same kind of file and names it `7aaq_ligand.pdb`.

Run `equilibrate_7aaq.py` with the same `python` that can import `openmm`. Compute target **Auto** leaves the platform empty, and OpenMM chooses it when the job starts.

## Orientation

`orient_7aaq.py` runs MemPro on `7aaq_protein.pdb`, then applies rank 1 to the whole file. BGC receives the same rigid move as the protein. The desktop steps are the MemPro figures in the [GUI study case](7aaq-gui.md). The structure checker and the sequence bar on that page are desktop-only.

| Parameter | Value | Meaning |
| --- | --- | --- |
| Iterations | `150` | MemPro minimisation steps |
| Grid size | `36` | Starting orientations |
| Device | `cpu` | Matches the desktop dialog |
| `SEED` | `7` | Same seed as the prepare-systems study case, so the pose repeats |

The desktop dialog can leave Seed empty. This script sets the seed so a second run keeps the same pose.

```bash
python orient_7aaq.py
```

The oriented structure is `7aaq_protein_oriented.pdb`. Ranked MemPro files stay in `mempro_7aaq_protein/`.

## Preparation

`prepare_7aaq.py` uses `PreparationManager` to predict pKa values and assign protonation at pH 6.

| Parameter | Value | Meaning |
| --- | --- | --- |
| `INPUT_PDB` | `7aaq_protein_oriented.pdb` | Oriented structure |
| `PH` | `6.0` | Working pH |
| `PRESERVE_RESIDUE_NUMBERS` | `False` | Sequential Amber residue numbers |
| `OUTPUT_DIR` | `01_preparation_7aaq_protein_oriented` | Output directory |

The desktop app leaves **Keep original residue numbers and chains** on. This script turns that off on purpose. Set `PRESERVE_RESIDUE_NUMBERS = True` when the API run should match a desktop folder that kept the original numbers.

```bash
python prepare_7aaq.py
```

Inside `01_preparation_7aaq_protein_oriented/` the script prints `OK` or `MISSING` for each expected name:

- `7aaq_protein_oriented.pka` — pKa predictions
- `7aaq_protein_oriented_summary_of_prediction.txt` — prediction summary
- `7aaq_protein_oriented_protonated.pdb` — structure for the Builder

Also expect the non-protein PDB, the renumbering map, and the disulfide-link file. Continue only when `7aaq_protein_oriented_protonated.pdb` is present and readable.

PropKa assigns HIP when the pH is below the histidine pKa and HIE otherwise. The desktop walkthrough then sets those HIP rows to HID by hand. This script keeps the automatic assignment.

## Builder

`build_7aaq.py` packs a POPC bilayer, TIP3P water, and 0.15 M KCl around the prepared protein, then parametrizes with tleap.

| Choice | Value |
| --- | --- |
| `water_model` | `tip3p` |
| `protein_ff` | `ff14SB` |
| `lipid_ff` | `lipid21` |
| `upper_lipids` / `lower_lipids` | `POPC` / `POPC` |
| `lipid_ratios` | `1//1` |
| `salt_concentration` | `0.15` |
| `cation` / `anion` | `K+` / `Cl-` |
| `dist` | `12` Å between the solute and the boundary |
| `dist_wat` | `26` Å water layer |
| `parametrize` | `True` |
| `notprotonate` | `True` so PropKa is not run again |
| `preoriented` | `True` because MemPro already oriented the protein |

```bash
python build_7aaq.py
```

The output folder is `02_build_7aaq_protein_oriented_protonated`. A folder that already contains a build is refused. Delete it, or set `OVERWRITE = True` in the script to reuse it. The script prints `OK` or `MISSING` for the files it expects. These four are the ones to check before equilibration:

- `system.prmtop` — Amber parameters and topology
- `system.inpcrd` — Amber coordinates
- `system.pdb` — the packed system
- `bilayer_7aaq_protein_oriented_protonated_lipid.pdb` — bilayer PDB used for the box

`packmol-memgen.log` and `leap.log` are the logs to read when a build fails.

## Equilibration

`equilibrate_7aaq.py` writes OpenMM inputs. It does not run OpenMM. The stage list starts with the same Minimization card as the desktop protocol (10,000 steps). `OpenMMEquilibrationManager.get_default_stage_params` then adds Equilibration 1 and 2 under NVT, Equilibration 3 through 6 under NPgT, and production under `SCHEME`. OpenMM writes the minimization steps into `step1_equilibration.inp`.

| Choice | Value |
| --- | --- |
| Engine | OpenMM |
| `SCHEME` | `NPT` for production |
| `TEMPERATURE` | `303.15` K |
| `PRODUCTION_NS` | `100.0` |
| Minimization | `10000` steps, written into step 1 |
| `cpu_cores` | `1` |
| `gpu_id` | `0` |
| `num_gpus` | `1` |
| Compute target | `auto` |
| `add_com_restraint` | `False` |
| `add_rotation_restraint` | `False` |

Those resource numbers are the desktop defaults after the engine is set to OpenMM.

```bash
python equilibrate_7aaq.py
```

The output folder is `03_equilibration_02_build_7aaq_protein_oriented_protonated`. The script checks `run_equilibration.sh`, `run_equilibration_cluster.sh`, `openmm_run.py`, `step1_equilibration.inp` through `step6_equilibration.inp`, `step7_production.inp`, and the copied system files.

## Run OpenMM

From the equilibration folder:

```bash
cd 03_equilibration_02_build_7aaq_protein_oriented_protonated
bash run_equilibration.sh
```

The script uses the `python` on `PATH` unless `PYTHON` is set. With compute target Auto, OpenMM picks the platform at startup. One CPU and GPU 0 are the counts written into the script.

After the run, read the log and confirm the coordinate, trajectory, and restart files exist. Temperature, energy, pressure, and RMSD are the usual checks that the system stayed stable. Those analyses are the last two sections of the [desktop walkthrough](7aaq-gui.md).
