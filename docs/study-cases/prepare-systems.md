# Prepare systems

This walkthrough times a membrane build for several proteins. One script prepares each system in POPC and writes how long each stage took. The membrane build itself is repeated three times. The six entries below are the systems filled in for the manuscript. Edit `SYSTEMS` to use other PDB ids, chains, or replica counts.

Download or open the script from GitHub: [`tests/study_cases/prepare_systems.py`](https://github.com/maurobedoya/gatewizard/blob/main/tests/study_cases/prepare_systems.py). Run it from the directory where `work/` and `results/` should be written. AmberTools, MemPrO, and a GateWizard environment that can import `gatewizard` need to be available there.

## What the script builds

Every system uses the same membrane settings:

| Choice | Value |
| --- | --- |
| Leaflets | POPC / POPC |
| Protein | ff14SB |
| Lipids | Lipid21 |
| Water | TIP3P |
| Salt | 0.15 M |
| Ions added by the builder | K+ and Cl− |
| Protonation pH | 7.4 |
| Builder replicas | 3 |

Packmol-memgen still uses its usual solute padding (`dist` 12 Å, water layer `dist_wat` 26 Å). Those are the builder defaults. Change them in the `prepare_system` call inside `build_once` if a study needs a different box.

The builder is the stage that is repeated. Download, chain selection, MemPrO, PropKa, pdb4amber, and GAFF run once per PDB. The three replicas share that prepared PDB, so the mean and the standard deviation describe packing plus tleap.

```python
BUILDER_REPLICAS = 3
PH = 7.4
```

A stopped run keeps finished stages. `python prepare_systems.py --fresh` times them again. `python prepare_systems.py 2MVH` runs one PDB and leaves the other rows in `results/timings.csv` and `results/summary.md`. Running that PDB again replaces only its rows.

## Decisions that are already written into the script

Each PDB is one entry in `SYSTEMS`. `chains` is the author chain list that is kept. `keep_hetero` is the residue names that stay with the protein. Anything else that is a HETATM, including crystal water, is left out, because the builder adds TIP3P and POPC.

```python
{
    "pdb_id": "6RV2",
    "model": 1,
    "chains": "AB",
    "keep_hetero": ["K"],
    "parametrize": {},
    "peripheral": False,
},
```

| PDB | What it is | What the script keeps | What you should confirm |
| --- | --- | --- | --- |
| 2MVH | SpoVM, NMR, 26 residues, 10 models, chain A | Model 1, chain A, no hetero groups. MemPrO runs in peripheral mode. | Model 1 is the RCSB representative. Peripheral mode matches an amphipathic peptide. |
| 1AFO | Glycophorin A transmembrane dimer, NMR, 20 models, chains A and B | Model 1 and both chains. | The biological assembly is the dimer. One chain would be a monomer control. |
| 1M0L | Bacteriorhodopsin, 1.47 Å, chain A, retinal | Chain A and RET (GAFF2, charge 0). Crystal lipids LI1 and SQU are omitted. | The deposited file is the monomer. The biological assembly is a symmetry trimer, and this script does not build that trimer. RET is parametrized as its own residue; the lysine–retinal LINK is not rebuilt as a bond. |
| 7AAQ | STP10, 1.81 Å, chain A | Chain A and the bound sugar BGC (GAFF2, charge 0). Acetate, the PEG detergent XPE, and monoolein OLC are omitted. | Confirm that the sugar stays and the crystallization additives go. |
| 6RV2 | TASK-1, 3.0 Å, four copies of the subunit | Chains A and B, plus the potassium ions named K (written as K+). Crystal PC, decyl maltoside, and cholesterol hemisuccinate are omitted. | The asymmetric unit contains two dimers. Assembly 1 is A+B (512 modeled residues). Assembly 2 is C+D. Change `chains` to `CD` for the second dimer. |
| 3WU2 | Photosystem II, 1.9 Å, 38 protein chains | Both halves of the dimer, chains `AaBbCcDdEeFfHhIiJjKkLlMmOoTtUuVvXxYyZz`, plus CA, CL, MG, and FE2. | Chlorophyll, pheophytin, carotenes, plastoquinone, heme, the Mn4Ca cluster, detergents, glycerol, and the thylakoid lipids are omitted. See below. |

2MVH is the only peripheral protein:

```python
results = MemPrO().run(
    str(protein_pdb),
    output_dir=str(mempro_dir),
    peripheral=bool(spec.get("peripheral")),
)
```

MemPrO is run on protein atoms. The same rigid transform is then applied to the kept hetero groups, so a ligand moves with the protein.

## Photosystem II cofactors

A POPC rebuild with ff14SB, Lipid21, and TIP3P can finish when the protein and the monoatomic ions are the solute. The organic cofactors are a separate decision. They are named in `PSII_ORGANIC` (`CLA`, `PHO`, `BCR`, `RRX`, `PL9`, `BCT`, the crystal lipids, `HEM`, `OEX`).

To carry an organic residue into the build, add its name to `keep_hetero` and its net charge to `parametrize`. The script then runs GAFF2 once and reuses the library for all three builder replicas.

`HEM` and `OEX` are refused by that path. They need `manual_params` with an existing `frcmod` and `lib`. The oxygen-evolving cluster is not a GAFF molecule.

## The builder call

Shared settings live in one function. `notprotonate=True` keeps the PropKa residue names written earlier. `preoriented=True` tells packmol-memgen that MemPrO already placed the protein.

```python
builder.prepare_system(
    pdb_file=str(pdb_path),
    working_dir=str(job_root),
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
    protein_ff="ff14SB",
    lipid_ff="lipid21",
    water_model="tip3p",
    salt_concentration=0.15,
    cation="K+",
    anion="Cl-",
    preoriented=True,
    parametrize=True,
    notprotonate=True,
    ligand_params=ligand_params,
    wait=True,
    wait_verbose=True,
)
```

`lipid_ratios="1//1"` is required with `POPC//POPC`. This packmol-memgen treats the double slash as two leaflet compositions and stops unless a ratio is given. `ligand_params` is empty when the PDB has no organic hetero group. BGC and RET fill it for 7AAQ and 1M0L.

## What is timed

Stages, in order:

1. `download` — RCSB PDB file
2. `subset` — chosen NMR model, chosen chains, waters and omitted hetero groups removed
3. `mempro` — orientation, then the transform applied to the kept atoms
4. `protonation` — PropKa at pH 7.4, residue names, then pdb4amber with `reduce`
5. `ligand_parameters` — GAFF2 for each name in `parametrize`
6. `builder` — one row per replica

Protein atoms are counted on `prepared.pdb` (ATOM records, hydrogens included). System atoms are counted on each replica `system.pdb` (ATOM plus HETATM). `results/summary.md` lists the individual builder times, the mean, and the sample standard deviation (n − 1). A failed replica stays in `results/timings.csv` and is left out of the mean.

Successful stage times are stored under `work/<PDB>/*.json`. That is what lets a long run continue after an interruption. Delete that PDB’s folder, or pass `--fresh`, when the table must come from a single uninterrupted pass.
