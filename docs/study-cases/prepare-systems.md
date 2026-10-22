# Prepare systems

This walkthrough times a membrane build for several proteins. One script prepares each system in POPC and writes how long each stage took. The membrane build itself is repeated three times. The entries below are the systems currently listed in `SYSTEMS`. Edit that list to use other PDB ids, chains, or replica counts.

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

Packmol-memgen on GateWizard was set to use a solute padding `dist` 12 Å, and water layer `dist_wat` 26 Å. Those are the builder defaults. Change them in the `prepare_system` call inside `build_once` if a study needs a different box.

The builder is the stage that is repeated. Download, chain selection, MemPrO, PropKa, pdb4amber, and GAFF run once per PDB, and every builder replica uses that same prepared PDB. Prepare one orientation, then pack the replicas from it. A second MemPrO run can shift the protein slightly and change the number of atoms. The same coordinates give packmol-memgen a reproducible system, so the replicas stay the same size and are easier to compare. The mean and the standard deviation are the packing plus tleap time.

```python
BUILDER_REPLICAS = 3
PH = 7.4
```

A stopped run keeps finished stages. `python prepare_systems.py --fresh` times them again. `python prepare_systems.py 2MVH` runs one PDB and leaves the other rows in `results/timings.csv` and `results/summary.md`. Running that PDB again replaces only its rows.

## Decisions that are already written into the script

Each PDB is one entry in `SYSTEMS`. `chains` is the author chain list that is kept. `keep_hetero` is the residue names that stay with the protein. A name left off that list is left out of that example. Waters, bound lipids, and detergents belong on the list when they are part of the question. The builder still adds bulk TIP3P and POPC around what you keep.

Experimental water stays when its residue name is listed in `keep_hetero`. In most PDB files that name is `HOH`. Water is not a ligand, so it does not go in `parametrize`. The water must also sit on a chain in `chains`. A blank chain column is dropped. Those waters then move with the protein during orientation, and the builder adds the rest of the solvent around them. Whether to keep them depends on the system: a water in a channel or a binding site is often part of the biology, and a loosely ordered crystal water may not be.

```python
{
    "pdb_id": "6RV2",
    "model": 1,
    "chains": "AB",
    "keep_hetero": ["K", "HOH"],
    "parametrize": {},
    "peripheral": False,
},
```

The systems on this page other than 3T4D do not list `HOH`, so those timing runs use TIP3P as the only solvent. 3T4D lists `HOH` because the crystal waters are part of that example. Add `"HOH"` when the crystallographic waters should stay. `WAT` is the other name this script recognizes.

| PDB | What it is | What this example keeps | Recommendation for a similar system |
| --- | --- | --- | --- |
| 2MVH | SpoVM, NMR, 26 residues, 10 models, chain A | Model 1, chain A, no hetero groups. MemPrO runs in peripheral mode. | For an amphipathic peptide, keep one NMR model and turn peripheral mode on. Transmembrane mode looks for a spanning axis that this kind of peptide does not have. |
| 1AFO | Glycophorin A transmembrane dimer, NMR, 20 models, chains A and B | Model 1 and both chains. | For an NMR dimer, keep one model and both chains of the biological assembly. One chain is a monomer, which is a different system. Leave peripheral mode off when the helices cross the bilayer. |
| 1M0L | Bacteriorhodopsin, 1.47 Å, chain A, retinal | Chain A and RET (GAFF2 with ABCG2, charge +1). Crystal lipids LI1 and SQU are left out of this example. | A deposited monomer can be built on its own when the biological assembly is a symmetry trimer. A retinal linked to a lysine as a protonated Schiff base needs charge +1. Keep a crystal lipid when it is biologically relevant, and give it a charge. Leave it out when POPC is standing in for the bilayer. |
| 7AAQ | STP10, 1.81 Å, chain A | Chain A and the bound sugar BGC (GAFF2 with ABCG2, charge 0). Acetate, the PEG detergent XPE, and monoolein OLC are left out of this example. | Keep a bound ligand and give antechamber its net charge. Acetate, the PEG, and the monoolein are crystallization additives in this deposit. Keep a detergent or a lipid when it is biologically relevant. |
| 3T4D | NaK2K channel Y55F, 1.70 Å, one monomer in the asymmetric unit | Biological assembly 2, the tetramer of the more complete chain. Crystal waters and the pore potassium stay. | Build the biological assembly when the deposit is only the asymmetric unit. Keep crystal water that belongs in the pore. The assembly step drops symmetry copies that land on the same point. |
| 6RV2 | TASK-1, 3.0 Å, four copies of the subunit | Chains A and B, plus the potassium ions named K. Crystal PC, decyl maltoside, and cholesterol hemisuccinate are left out of this example. | When the crystal holds two copies of a dimer, orient one copy. Chains A and B are one dimer; C and D are the other. Keep ions that belong in the protein, such as selectivity-filter potassium. Leave a detergent or a crystal lipid out only when it is not biologically relevant. |
| 3WU2 | Photosystem II, 1.9 Å, 38 protein chains | Both halves of the dimer, chains `AaBbCcDdEeFfHhIiJjKkLlMmOoTtUuVvXxYyZz`, plus CA, CL, MG, and FE2. Chlorophyll, heme, the metal cluster, and the thylakoid lipids are left out of this example. | Those cofactors and lipids are part of the complex. This example leaves them out so the protein-and-ions build is shorter. Include them when the question needs them and the parameter files exist. Match deposited ion atom names to the Amber templates before tleap. |

2MVH is the only peripheral protein:

```python
results = MemPrO().run(
    str(protein_pdb),
    output_dir=str(mempro_dir),
    peripheral=bool(spec.get("peripheral")),
    seed=MEMPRO_SEED,
)
```

`MEMPRO_SEED` is `7`. MemPrO 0.3.1 writes the same orientation again for the same protein, seed, and CPU count.

MemPrO is run on protein atoms. The same rigid transform is then applied to the kept hetero groups, so a ligand moves with the protein.

## Atom counts and times

Every system is a symmetric POPC bilayer. The rows below, and the `SYSTEMS` list, follow the packed system size from smallest to largest. The protein atom count is `prepared.pdb` after protonation. The total atom count is every atom in `system.pdb` after packmol-memgen and tleap: protein, lipids, water, and ions. Those counts were measured before ACE and NME caps. A new run adds the caps unless the PDB already has them, so both counts will be larger. Builder times are packing plus parametrization. MemPrO, protonation, and GAFF run once and are not inside the builder mean.

| PDB | Prepared protein atoms | Total atoms | Builder mean | Note |
| --- | ---: | --- | ---: | --- |
| 2MVH | 447 | 38,242 in each replica | 20.0 min (std 8 s) | Short peptide. A practical first system. |
| 1AFO | 1,291 | 45,542 in each replica | 20.0 min (std 5 s) | NMR dimer. One model, both chains. |
| 1M0L | 3,568 | 58,954 in each replica | 23.6 min (std 13 s) | Linked retinal uses charge +1. Long termini stay unbuilt. |
| 7AAQ | 7,580 | 84,695, 84,683, 84,689 | 31.6 min (std 9 s) | Bound sugar stays. Crystallization additives do not. |
| 3T4D | 6,143 | 88,914 in each replica | 30.9 min (std 9 s) | Assembly 2. Four potassium ions stay. This run’s crystal waters were removed by pdb4amber `--dry` before packing. |
| 6RV2 | 8,266 | 129,701, 129,698, 129,707 | 43.6 min (std 2.0 min) | One dimer, chains A and B. |
| 3WU2 | 80,898 | 904,059, 904,068, 904,086 | 5.04 h (std 1.6 min) | Protein and monoatomic ions. All three replicas are in the mean. |

Prepare one orientation and pack every replica from it. A second MemPrO run can move the protein enough to change the atom count. The same prepared protein gives packmol-memgen matching systems that are easier to compare.

Orientation was the slow step before packing only for Photosystem II (18.5 min). The other proteins oriented in 83–116 s. GAFF ran for two ligands: retinal took 25 s and the bound sugar took 1.7 s. 3T4D has no organic ligand. Its assembly step took 0.09 s, and MemPrO took 116 s.

### 2MVH — SpoVM

SpoVM is an amphipathic peptide, not a transmembrane helix bundle. The deposited NMR entry has 10 models and one chain. This example keeps model 1 and chain A, the RCSB representative, and drops every hetero group.

MemPrO needs peripheral mode for an amphipathic peptide. Transmembrane mode looks for a membrane-spanning axis that this peptide does not have.

```python
{
    "pdb_id": "2MVH",
    "model": 1,
    "chains": "A",
    "keep_hetero": [],
    "parametrize": {},
    "peripheral": True,
}
```

The structure checker reports nothing before or after preparation: no missing heavy atoms, no gaps, no alternate conformations. Protonation adds one atom (446 coordinates become 447). Each built system has 38,242 atoms. The OpenMM protocol dict in the script is commented out. Uncomment it to run minimization, Equilibration 1–6, and 1 ns of production on a box of this size.

### 1AFO — glycophorin A

This is an NMR transmembrane dimer, 20 models, chains A and B. Sending every model to MemPrO stacks them into one solute. This example keeps model 1 and both chains.

```python
{
    "pdb_id": "1AFO",
    "model": 1,
    "chains": "AB",
    "keep_hetero": [],
    "parametrize": {},
    "peripheral": False,
}
```

Leave `peripheral` off for a helix pair that spans the bilayer. One chain builds a monomer, which is a different system.

Before preparation the checker is empty. After preparation, chains A and B each have 40 residues and match at 100% identity. Hydrogens are present on 76 of 80 residues, which is the Amber protonation, not a missing-atom error. The prepared protein has 1,291 atoms (the subset had 1,322, because the NMR file already contained hydrogens that preparation replaces). Each built system has 45,542 atoms.

### 1M0L — bacteriorhodopsin

The deposited asymmetric unit is one chain at 1.47 Å, with retinal and two crystal lipids. The biological assembly is a symmetry trimer. This example builds the monomer. `LI1` and `SQU` are left out here because POPC stands in for the bilayer. Keep either lipid when it is biologically relevant, add its name to `keep_hetero`, and set its net charge in `parametrize`.

Retinal stays. The deposited formula looks like neutral retinal, but the oxygen is absent and a `LINK` record joins C15 to Lys216. After hydrogens are added the residue is C20H27. Charge 0 then has an odd electron count and antechamber refuses a singlet. Charge +1 is the protonated Schiff base. The lysine stays `LYS`, and the `LINK` is not rebuilt as a covalent bond.

```python
{
    "pdb_id": "1M0L",
    "model": 1,
    "chains": "A",
    "keep_hetero": ["RET"],
    "parametrize": {"RET": 1},
    "peripheral": False,
}
```

The subset has 1,756 protein atoms and 20 retinal atoms. Preparation writes 3,568 protein atoms and 47 retinal atoms. GAFF for that retinal takes 25 s and is reused by all three replicas. Each built system has 58,954 atoms.

Before preparation, SEQRES lists 17 residues at one end and 18 at the other that have no coordinates. Both stretches are longer than 6 residues, so the short-gap fill leaves them open. The report after preparation does not repeat SEQRES, which is why that count becomes zero. Retinal remains listed as a non-standard residue. That is expected.

### 7AAQ — STP10

STP10 is a sugar/proton symporter at 1.81 Å, one protein chain. The bound sugar `BGC` stays at charge 0. Acetate (`ACT`), the PEG detergent `XPE`, and monoolein (`OLC`) are left out of this example. In this deposit they come from crystallization, and POPC is the membrane that gets packed. Keep a detergent or a lipid that is biologically relevant, and parametrize it the same way as the sugar.

```python
{
    "pdb_id": "7AAQ",
    "model": 1,
    "chains": "A",
    "keep_hetero": ["BGC"],
    "parametrize": {"BGC": 0},
    "peripheral": False,
}
```

The subset has 3,771 protein atoms and 12 sugar atoms. Preparation writes 7,580 protein atoms and 24 sugar atoms. GAFF for the sugar takes 1.7 s. The three built systems have 84,695, 84,683, and 84,689 atoms.

SEQRES has 20 names without coordinates at one end and 7 at the other. Seven is already past the short-gap limit, so neither stretch is built. After preparation the checker only notes that `BGC` is non-standard and that hydrogens are present on 481 of 487 residues.

### 3T4D — NaK2K channel

The deposit is the asymmetric unit of the Y55F mutant, one copy of the subunit in chain A and one in chain B. The channel is a tetramer. Assembly 1 expands chain A. Assembly 2 expands chain B. Chain B has 96 of the 97 SEQRES residues (the N-terminal methionine has no coordinates). Chain A has 93. This example builds assembly 2, the more complete tetramer.

The build is the same one as Visualize → Build assembly, with crystal waters included and extra copies that land on the same point removed. The cutoff is 0.5 Å. On this structure that removes 12 potassium atoms. Each operator copies the pore ions onto one site. Four potassium ions remain, one set for the tetramer. The 19 crystal waters on chain B are copied onto the four subunits and do not land on each other, so all 76 water atoms stay.

The expanded chains are named B1, B2, B3, and B4. A PDB chain column holds one character, so the script renames them A, B, C, and D before MemPrO. SEQRES from deposited chain B is written onto each of those chains. `chains` is `ABCD` because the biological unit is the whole tetramer.

```python
{
    "pdb_id": "3T4D",
    "model": 1,
    "assembly": {"id": "2", "waters": True, "drop_overlaps": True},
    "assembly_from": "B",
    "chains": "ABCD",
    "keep_hetero": ["HOH", "K"],
    "parametrize": {},
    "peripheral": False,
}
```

`assembly_from` is the deposited chain the operators apply to. Short-gap filling runs on that chain before the assembly is built. This entry has no internal gap of 1–6 residues. `waters` keeps `HOH`. `drop_overlaps` is the checkbox that removes the piled-up copies. Potassium is an ion, so it stays in `keep_hetero` without a GAFF charge.

The assembly file has 3,020 protein atoms, 76 crystal water atoms, and 4 potassium atoms. MemPrO takes 116 s. Preparation writes 6,143 protein atoms. Before preparation the checker lists one missing N-terminal methionine on each of the four chains. After preparation the four chains still match at 100% identity, and hydrogens are present on 380 of 384 residues.

Each packed system has 88,914 atoms. The builder mean is 30.9 min (sample standard deviation 9 s). Packmol is 30.6, 30.4, and 30.4 min of that. pdb4amber is about 20 s, and tleap is 2.3 s.

That timed run called pdb4amber with `--dry`, so the 76 crystal water atoms were removed before packing. The four potassium ions stayed. The script now leaves a water whose name is in `keep_hetero`. A new run of this PDB keeps those waters and will be larger than 88,914 atoms.

Run this PDB on its own when the other systems should be left as they are:

```bash
python tests/study_cases/prepare_systems.py 3T4D
```

### 6RV2 — TASK-1

TASK-1 (K2P3.1) at 3.0 Å deposits two dimers in the asymmetric unit: A+B and C+D. The asymmetric unit is the crystal packing, not one channel. This example keeps one dimer, chains A and B. Set `chains` to `"CD"` for the other dimer.

```python
{
    "pdb_id": "6RV2",
    "model": 1,
    "chains": "AB",
    "keep_hetero": ["K"],
    "parametrize": {},
    "peripheral": False,
}
```

Selectivity-filter potassium stays as residue `K` and atom `K`. pdb4amber in the builder renames those to `K+`, which is the same step the GUI uses. Crystal phosphatidylcholine (`PC1`), decyl maltoside (`DMU`), and cholesterol hemisuccinate (`Y01`) are left out of this example. Leave a detergent or a crystal lipid out only when it is not biologically relevant. A lipid that belongs in the channel goes in `keep_hetero` with a net charge in `parametrize`.

The subset has 4,066 protein atoms and 4 potassium atoms. Preparation fills 82 missing side-chain atoms and writes 8,266 protein atoms. The three built systems have 129,701, 129,698, and 129,707 atoms. This is the first case where the builder replicas spread by more than a few seconds (mean 43.6 min, sample standard deviation 2.0 min).

Chain B is 4 residues longer than chain A (261 versus 257) and the sequences match at 98% identity. That length difference is in the deposit. It is still there after preparation. SEQRES also lists 7 residues on chain A and 3 on chain B with no coordinates. The after report has no internal numbering gap, and it does not repeat those SEQRES-only residues.

### 3WU2 — photosystem II

Photosystem II at 1.9 Å is a dimer. The author chains use both cases, and both halves are kept: `AaBbCcDdEeFfHhIiJjKkLlMmOoTtUuVvXxYyZz`. Monoatomic ions stay: `CA`, `CL`, `MG`, and `FE2`. Chlorophyll, pheophytin, carotenes, plastoquinone, heme, the Mn4Ca cluster, detergents, glycerol, and the thylakoid lipids are part of the biology of this complex. This example leaves them out so the protein-and-ions system packs faster. Include them when the question needs them. A POPC bilayer with ff14SB, Lipid21, and TIP3P can still be built from the protein and those ions while you decide.

```python
{
    "pdb_id": "3WU2",
    "model": 1,
    "chains": "AaBbCcDdEeFfHhIiJjKkLlMmOoTtUuVvXxYyZz",
    "keep_hetero": ["CA", "CL", "MG", "FE2"],
    "parametrize": {},
    "manual_params": {},
    "peripheral": False,
}
```

The subset has 41,092 protein atoms and 15 ion atoms. Preparation writes 80,898 protein atoms and the same 15 ions. MemPrO takes 18.5 min. Protonation takes 86 s. Short-gap filling takes 15 s. The three systems have 904,059, 904,068, and 904,086 atoms. Their builder times are 5.01 h, 5.06 h, and 5.05 h (mean 5.04 h, sample standard deviation 1.6 min). Packmol is 4.95 h, 5.00 h, and 4.99 h of that. pdb4amber is 169 s, 172 s, and 168 s. tleap is 25 s, 24 s, and 24 s. Expect the large-system PDB rewrite: serials pass 99,999, and tleap needs the realigned file.

Before preparation the checker reports 314 missing heavy atoms, 63 alternate conformations, and 107 residues that are absent from the coordinates. Preparation fills the missing protein atoms, keeps the higher-occupancy alternate conformation, and leaves two residues open: homoserine lactone `HSK` at residue 336 on chains D and d. `HSK` is not a standard amino acid, so the short-gap builder does not invent it. The chain-difference list is long because the checker compares every pair of chains in the dimer. Different subunits are supposed to differ. `FE2` stays listed as a non-standard residue.

The deposited iron is residue `FE2` and atom `FE`. Before the build, that atom is renamed `FE2`, which is the name in Amber's ion template. A deposited `FE3` residue would be written as Amber residue `FE` (atom `FE`), the Fe3+ unit in the same library. `CU1`, `GD3`, `IOD`, `NH4`, and `U4+` get the library atom name when the deposited atom is only the element. `CA`, `CL`, `MG`, and `K` already match and stay as deposited.

## Photosystem II cofactors

Chlorophyll, pheophytin, carotenes, plastoquinone, heme, the oxygen-evolving cluster, and the thylakoid lipids are part of photosystem II. This example leaves them out for a shorter protein-and-ions pack. Add a name when that run should include it. They are listed in `PSII_ORGANIC` (`CLA`, `PHO`, `BCR`, `RRX`, `PL9`, `BCT`, the thylakoid lipids, `HEM`, `OEX`).

To carry an organic residue into the build, add its name to `keep_hetero` and its net charge to `parametrize`. The script then parametrizes it once and reuses the library for all three builder replicas. Expect a larger and slower pack than the protein-and-ions solute.

`HEM` and `OEX` need `manual_params` with an existing `frcmod` and `lib`. The oxygen-evolving cluster is not a GAFF molecule. Add those names when the libraries exist.

## The builder call

Shared settings live in one function. `notprotonate=True` keeps the PropKa residue names written earlier and stops packmol-memgen from running reduce, which renames `HA` to `HCA` and breaks tleap. `remove_protein_h=True` is the Builder option “Remove protein hydrogens”: it drops protein hydrogens left with those names, and tleap fills `HA` from the residue template. Ligand hydrogens stay. `preoriented=True` tells packmol-memgen that MemPrO already placed the protein.

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
    remove_protein_h=True,
    ligand_params=ligand_params,
    wait=True,
    wait_verbose=True,
)
```

`lipid_ratios="1//1"` is required with `POPC//POPC`. This packmol-memgen treats the double slash as two leaflet compositions and stops unless a ratio is given. `ligand_params` is empty when the PDB has no organic hetero group. BGC (charge 0) and RET (charge +1) fill it for 7AAQ and 1M0L. The integer in `parametrize` is the net charge passed to antechamber.

The default ligand model is GAFF2 with ABCG2, which is one of the two pairings Amber recommends. The other is GAFF with AM1-BCC. Set `atom_type` to `gaff2` or `gaff`, and `charge_method` to `abcg2` or `bcc`, on the system entry. One residue can override both by using a dict instead of an integer. GAFF2 with AM1-BCC, and GAFF with ABCG2, still run, and the script prints Amber's warning for those pairs. Delete that PDB's `ligand_parameters.json` after changing the model, or the saved library is reused.

```python
"atom_type": "gaff",
"charge_method": "bcc",
"parametrize": {"RET": 1},
```

```python
"parametrize": {
    "RET": {"charge": 1, "atom_type": "gaff2", "charge_method": "abcg2"},
    "BGC": {"charge": 0, "atom_type": "gaff", "charge_method": "bcc"},
},
```

After packmol-memgen, the generated `run_preparation.sh` checks the packed PDB. Residue numbers past 9998, or a serial written against the atom name, are rewritten before `pdb4amber`. A readable membrane is left as packmol wrote it. tleap then receives disulfide bonds for CYX SG pairs and `set system box` from the CRYST1 record. `loadPDB` does not copy that box. The `system.pdb` written by `savePDB` is rewritten when its columns are shifted, so the equilibration restraints can read it. The original `bilayer_*_lipid.pdb` is kept for the box as well.

## Running minimization and equilibration

`PROTOCOL` at the top of the script is off. A system runs the membrane protocol only when its own `protocol` dict sets `enabled` to `True`. That dict can choose the engine (`amber`, `gromacs`, `openmm`, or `namd`), whether minimization and equilibration run, and `production_ns`. `stage_ns` replaces individual stage lengths. `None` keeps the GateWizard schedule: 10000 minimization steps, Equilibration 1–6 (Equilibration 6 is 47.625 ns), then production. COM and rotation restraints stay off, as on the Equilibration page. The run is `protocol/run_equilibration.sh` inside the builder replica, and only the first replica runs it unless `replicas` is raised.

2MVH is the small system. Its protocol dict is commented out, so a run only prepares that system. Uncomment the dict to run OpenMM minimization, every equilibration stage, and 1 ns of production.

```python
"protocol": {
    "enabled": True,
    "engine": "openmm",
    "minimization": True,
    "equilibration": True,
    "production_ns": 1.0,
    "openmm_platform": "CUDA",
},
```

## What is timed

Stages, in order:

1. `download` — RCSB PDB file
2. `short_gaps` — numbering gaps of 1 to 6 standard residues are built from SEQRES before orientation. Only those new atoms are inserted into the downloaded PDB, and later atom serials move up. Ligands and the other records stay as deposited. Longer gaps stay open. A sequence that is not a standard amino acid, such as HSK, stays open. When nothing is built, the downloaded PDB is passed on unchanged. A system with `assembly_from` fills gaps on that deposited chain.
3. `assembly` — only when the system sets `assembly`. Biological assembly operators expand that chain, crystal waters stay when `waters` is set, and atoms that land on the same point are dropped when `drop_overlaps` is set. Copied chains are renamed to one letter so the later PDB stages can read them.
4. `subset` — chosen NMR model, chosen chains, and hetero groups that are not in `keep_hetero` removed
5. `structure_check_before` — the same report as Visualize → Structure checker, on the selected chains. SEQRES for those chains is included, so residues with no coordinates are listed. The report is `work/<PDB>/structure_check_before/report.txt`. Findings do not stop the run.
6. `mempro` — orientation, then the transform applied to the kept atoms
7. `protonation` — ACE and NME are added first, unless the PDB already has a terminal cap, then PropKa at pH 7.4 and Amber residue names, peptide names remapped, disulfide cysteines renamed to CYX, missing protein heavy atoms filled from the Amber templates, protein hydrogens removed, then `pdb4amber --reduce` adds Amber hydrogens. Ligands and ions keep their hydrogens. Capping matches Preparation with Cap protein termini on. A finished `protonation` stage is reused, so delete that PDB’s `protonation.json` when this step must run again.
8. `structure_check_after` — the same checker on `prepared.pdb`. SEQRES is left out here because protonation renames residues. Numbering gaps and any heavy atoms the templates did not fill are still listed, in `work/<PDB>/structure_check_after/report.txt`.
9. `ligand_parameters` — GAFF2 with ABCG2 for each name in `parametrize`, unless that system sets `atom_type` and `charge_method`
10. `builder` — one row per replica, the whole packmol-memgen plus tleap run
11. `packmol`, `pdb4amber`, and `tleap` — the same replica split by the files that build writes. `packmol.log` ends packing, `system_for_tleap.pdb` ends pdb4amber, and `system.prmtop` ends tleap. A replica that is still packing stays as the single builder row.
12. `protocol` — only when that system's `protocol` is enabled, and only for the first replica unless `replicas` is higher. This is minimization, equilibration, and production.

`results/timings.csv` and `results/summary.md` are rewritten after every finished stage and every builder replica. Times under one second stay in the summary (they are not rounded to 0.0). Stopping the run keeps the rows that already finished, including protein atom counts from `prepared.pdb` and system atom counts from each finished `system.pdb`. A failed replica stays in the table and is left out of the mean. `python prepare_systems.py --report` rewrites that table from the saved stages and the build directories. It does not start a stage.

Successful stage times are stored under `work/<PDB>/*.json`. That is what lets a long run continue after an interruption. Delete that PDB’s folder, or pass `--fresh`, when the table must come from a single uninterrupted pass.

## Recommendations

- Orient the protein once, then build the replicas from that prepared PDB. Repeating MemPrO can shift the protein and change the atom count. The same coordinates give packmol-memgen systems of equal size, which are easier to compare.
- Start with a short peptide such as SpoVM when learning the script. It is the peripheral case, the structure checker is empty, and a bilayer is ready in about 20 minutes.
- For an NMR entry, set `model` to one conformer. Glycophorin A has 20 models and SpoVM has 10. The whole file is not one protein.
- For a crystal that contains two copies of an oligomer, orient one copy. TASK-1 chains A and B are one dimer. Chains C and D are the other.
- When the deposit is the asymmetric unit, set `assembly` to the biological assembly that should be built. NaK2K uses assembly 2, keeps `HOH`, and drops the potassium copies that land on one point. The script renames the expanded chains to one letter.
- Put an organic residue in both `keep_hetero` and `parametrize`. The integer is the net charge antechamber receives. Retinal is `{"RET": 1}`. The sugar is `{"BGC": 0}`. Ions such as `K` and `FE2`, and waters such as `HOH`, stay in `keep_hetero` only. The default model is GAFF2 with ABCG2. Set `atom_type` to `gaff` and `charge_method` to `bcc` for GAFF with AM1-BCC.
- Keep experimental water, a bound lipid, or a detergent when it is biologically relevant. Add the residue name to `keep_hetero`, and add a net charge to `parametrize` when the residue is organic. Leave a detergent or a crystal lipid out only when it is not biologically relevant.
- Photosystem II chlorophyll, heme, the metal cluster, and the thylakoid lipids are part of the complex. The example leaves them out so the protein-and-ions pack is shorter. Add the names from `PSII_ORGANIC` when the run should include them. `HEM` and `OEX` need `manual_params`.
- Treat a SEQRES stretch longer than 6 residues as unmodeled. The script will not build it. Homoserine lactone on photosystem II stays open even though it is a single missing residue, because it is not a standard amino acid.
- Read `structure_check_after/report.txt` for gaps that are still open. Hydrogens on most residues, a named ligand, and chain pairs that are different subunits are listed there and are not failed builds.
- Plan on the builder, not MemPrO, for the wall time. Only photosystem II spends more than a couple of minutes in orientation.
