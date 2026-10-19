# User Guide

This guide matches the **[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)** desktop app (Electron). The `gatewizard` pip package is the Python API only — it does not open this window.

## Pages

The left **activity bar** switches pages:

| Page | What it does |
|------|----------------|
| **Visualize** | 3D structures, trajectories, representations, sequence and short loops, structure checker, minimize, mutate / split / superimpose, figures and movies |
| **Preparation** | PropKa pKa table, protonation, disulfide bonds, ACE/NME caps, short loops, alternate conformations, inline 3D |
| **Builder** | packmol-memgen leaflets, salt, ligand force field and charges, tleap, job cards |
| **Equilibration** | NAMD / GROMACS / OpenMM / Amber schemes, Generate, Watch / Pull, cluster |
| **Tools** | Fix PBC and other trajectory utilities |
| **Analysis** | Structural (RMSD/RMSF/APL) and energetic log plots |

![Visualize workspace](images/user-guide/visualize_workspace.png)

*Visualize workspace*

## Visualization

1. **Open** a PDB/mmCIF (or paste an RCSB id).
2. Use the toolbar: **Open / Select / Tools / Save Image**.
3. **Representations** dock: cartoon, licorice, surface, coloring, selections.
4. Orbit with the left mouse button; **Ctrl+scroll** zooms more slowly.

![Create representations](images/user-guide/visualize_workspace_7aaq.png)

### Trajectories

Load a topology + DCD/XTC/TRR (or PSF+DCD). The play slider and optional **Traj-smooth** cards appear once frames are available.

![Trajectory](images/user-guide/visualize_trajectory.png)

### Animation

Turn **Animate** on for a keyframe timeline. Export GIF/WebM/PNG (transparent background is available for those formats; MP4/MOV stay opaque).

![Animation](images/user-guide/visualize_animation.png)

### Tools on Visualize

**Tools ▾** includes mutate residue, split chain, superimpose, and related edits. Mutate writes a new structure into the session.

![Mutate / tools](images/user-guide/visualize_tools_mutate.png)

**Save Image** opens the figure dialog (size, scale, transparent PNG, yellow crop overlay).

### Sequence bar and short loops

**Sequence** opens the protein sequence by chain. Click a residue to add it as licorice; click it again to remove only that residue. Drag to add a range. Right-click adds another representation.

A numbering hole of **1 to 6** residues is marked on that bar. **Fill short gaps** builds it with PDBFixer from SEQRES, or from a sequence you type, and saves a new structure. A longer hole is listed and left open. A residue that is not a standard amino acid, such as HSK, stays open. The same checkbox is on Preparation.

### Structure checker and alternate conformations

**Structure checker** lists missing atoms, gaps, alternate conformations, and chain differences. An alternate conformation is the PDB altLoc letter. The more likely letter is the higher occupancy. A tie keeps the earlier letter.

A ligand, water, or other molecule is its own group, even when it reuses a protein residue number. Atoms that overlap and have no altLoc letter are listed too. Tick one letter to draw only that pose. Tick both to draw both, with the more likely pose in yellow.

**Build assembly** expands the biological assembly already stored in the file. Crystal waters can be included. Copies that land on the same point can be dropped.

### Minimize a selection

The minimize control on the Representations toolbar relaxes the residues, atoms, or molecule already selected. The icon uses the last settings. A new session starts at 50 steps. The corner mark sets the step count, the convergence, and how far the atoms may move, without starting a run.

This is a short clash fix, not a finished energy minimum. Atoms Amber can build stay fixed except the selection. A water or another molecule Amber cannot build is relaxed without breaking its bonds. Potassium and the other Amber monoatomic ions can be minimized. The result is a new PDB, `name_min_001.pdb` and then `name_min_002.pdb`, in a `minimized` folder. The structure you started from stays.

## Preparation

PropKa table plus an inline 3D view of the working PDB.

![Preparation](images/user-guide/preparation.png)

1. Set the working directory (top bar) and load a PDB.
2. Choose pH and optional **protein capping** (ACE/NME).
3. Run analysis. The table lists residue, Old ID, New ID, chain, pKa, and state at pH.
4. Auto-detect disulfide bonds if needed, then **Prepare** (missing protein heavy atoms from Amber templates, PropKa names, AmberTools hydrogens). Missing loop residues are not built. Stale CONECT/LINK records are dropped. If the Amber `reduce` binary is missing, Prepare keeps tleap hydrogens and still writes a prepared PDB.

**Keep original residue numbers and chains** is on by default. PropKa and AmberTools still add hydrogens; residue and chain ids stay as in the input (a missing loop that jumps from 200 to 205 stays that gap). Uncheck it for sequential Amber numbering (tleap / MD).

**Fill short gaps (1–6 residues)** builds those holes with PDBFixer when the sequence is in SEQRES, or from a sequence you type. A longer gap is listed and left open. A residue that is not a standard amino acid stays open. Visualize shows the same marks on the sequence bar and can run the same fill as a new structure.

**Keep the most likely alternate conformations** is on by default. Each residue keeps the higher-occupancy letter, including ligands, waters, and other molecules. A tie keeps the earlier letter. Turn it off to leave the choice to pdb4amber, which keeps the first copy of each atom. Choose the letters on Visualize when the assignment matters.

Capping is for fragments and termini that should not stay charged in MD. The capper still writes ACE as residue 1 and shifts the chain for its own file; **Prepare** with preserve on puts the protein numbers back. New ACE/NME then get N-terminus − 1 (residue **0** if the chain starts at 1) and C-terminus + 1. The capper mapping file is kept for later Builder/PropKa steps.

## Builder

Leaflet composition, water/ions, and parametrization. Each run is a **job card** under the project folder.

![Builder](images/user-guide/builder.png)

1. Working file: prepared protein (or uncheck protein for bilayer-only and set membrane XY).
2. Upper/lower lipids and ratios. GateWizard supports **4028 lipids** from the packmol-memgen catalog (`--available_lipids_all`, version 2025.1.29). 258 of them are Amber lipid21 or lipid17. The other 3770 are parametrized with PACKMOL-Memgen's Lipid_ext library on top of that Amber leaprc; those parameters are not in Amber's `lipid21.lib`. The picker marks those rows Lipid_ext. The catalog lists 4055 names; the 27 SIRAH lipids are coarse-grained and are not supported. A short familiar subset is in the [builder reference](api/builder.md).

![Lipid ratios](images/user-guide/builder_lipid_ratios.png)

3. Water model, protein FF, lipid FF, salt, water thickness.
4. Validate, then generate inputs / start preparation. Watch the job card; **Cancel** stops a running `run_preparation.sh`.

Keep **Skip protonation** on when the PDB already has PropKa names (ASH, GLH, HIP, CYX).

Each ligand and free molecule can use **GAFF** or **GAFF2**, and **ABCG2** or **AM1-BCC**. The default is GAFF2 with ABCG2. GAFF2 with AM1-BCC, and GAFF with ABCG2, show Amber's warning before parametrization. The Final drawing uses the bonds in the parametrized mol2. A ring bond stored as a triple bond is drawn from its length, so a piperazine and an aromatic ring stay ordinary rings.

## Equilibration

Pick an engine (NAMD, GROMACS, OpenMM, Amber), a scheme (NVT/NPT/…), then **Generate**.

![Equilibration](images/user-guide/equilibration.png)

### Jobs

Progress cards: **Watch** (local logs), **Pull** (sync from a remote job), and per-stage status.

![Equilibration progress](images/user-guide/equilibration_progress.png)

### Cluster

**Run on cluster** uses Settings → Clusters (SSH, workdir, Slurm). The dialog writes a cluster run script next to the local equilibration folder.

![Cluster dialog](images/user-guide/equilibration_cluster.png)

## Tools

**Fix PBC** unwraps/centers a trajectory for analysis or Visualize. Other tools on this page stay job-based (status JSON in the project folder).

![Fix PBC](images/user-guide/tools_fix_pbc.png)

Folder-wide job scans pause while you are on Visualize so playback stays smooth.

## Analysis

Two modes on the same page:

- **Structural** — RMSD, RMSF, distances, Rg, area-per-lipid / thickness (FATSLiM, GridMAT, or lipyphilic).
- **Energetic** — NAMD / OpenMM / GROMACS / Amber log plots. Engine selector chooses the parser.

![Structural analysis](images/user-guide/analysis_structural.png)

![Energetic analysis](images/user-guide/analysis_energetic.png)

Selections use MDAnalysis syntax (`protein and name CA`). Multi-file runs need per-file times so the x-axis is in ns, not concatenated steps. Details: [Analysis features](analysis.md).

## Settings and clusters

Gear → **Clusters**: host, user, key, remote workdir, optional scratch. GPU safe mode is session-only unless you change it here.

![Settings clusters](images/user-guide/settings_clusters.png)

## Status bar

The footer shows backend connectivity and **system RAM** (not polled while you stay on Visualize).

## Python API

For scripts, use the [API reference](api/index.md) and [method dictionary](api/dictionary.md). Worked examples are the `tests/<module>_examples/` scripts in the API repository.
