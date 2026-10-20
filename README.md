# GateWizard (API)
<img src="./docs/images/splash.png" alt="GateWizard" width="140" align="right" />

[![PyPI version](https://img.shields.io/pypi/v/gatewizard.svg)](https://pypi.org/project/gatewizard/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://maurobedoya.github.io/gatewizard/)
[![DOI](https://zenodo.org/badge/1073861334.svg)](https://doi.org/10.5281/zenodo.18264074)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/maurobedoya/gatewizard.svg)](https://github.com/maurobedoya/gatewizard/stargazers)

Python library for membrane protein preparation, system building, equilibration setup, and MD trajectory analysis.

**This repository is the API only.** It does not include a graphical interface.

If you want a desktop app, use **[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)**, a separate project that calls this library. You can use the API without the GUI, and the GUI without cloning this repo (it installs the API via pip).

📖 [Documentation](https://maurobedoya.github.io/gatewizard/) · [Manual (PDF)](docs/manual/GateWizard-manual.pdf) · 💬 [Discussions](https://github.com/maurobedoya/gatewizard/discussions) (API + GUI community forum)

## Features

- **Preparation** — clean structures, pKa (PROPKA), protonation, and termini capping
- **Membrane builder** — orient protein, pack lipids and solvent, Amber parametrization (tleap). GateWizard supports **4028 lipid types** from packmol-memgen `--available_lipids_all` (version 2025.1.29). 3770 of them parametrize through PACKMOL-Memgen's Lipid_ext library on top of lipid21 or lipid17; those parameters are not in Amber's `lipid21.lib`. The catalog lists 4055 names; the 27 SIRAH lipids are coarse-grained and are not supported. See the [builder reference](docs/api/builder.md).
- **Equilibration** — the same CHARMM-GUI-style protocols for four MD engines: NAMD, GROMACS, OpenMM, and Amber (NVT, NPT, NPAT, NPγT)
- **Analysis** — trajectory metrics, structural analysis, and plotting helpers
- **Force fields** — Amber protein models (e.g. ff14SB, ff19SB) and common water/lipid setups

## Install

```bash
conda create -n gatewizard -c conda-forge python sqlite ambertools=26 parmed=4.3.0 openmm cudatoolkit gromacs -y
conda activate gatewizard
conda install -c conda-forge -c bioconda reduce -y
pip install "gatewizard[full]"
```

That conda line is for Linux or WSL with an NVIDIA driver. On macOS, leave out `cudatoolkit`. NAMD is not in this environment; install it separately and put `namd3` on `PATH` if you run NAMD equilibration.

`gatewizard[full]` installs the CPU build of [MemPrO](https://github.com/pstansfeld/MemPrO) (`mempro` ≥ 0.3.0). A checkout’s `environment.yml` does the same.

For the library without the OpenMM extra (NAMD, GROMACS, and Amber still use their own programs):

```bash
pip install gatewizard
```

From a git checkout (`environment.yml` already runs `pip install -e ".[full]"`, including CPU GROMACS, OpenMM, and `cudatoolkit`):

```bash
git clone https://github.com/maurobedoya/gatewizard.git
cd gatewizard
conda env create -f environment.yml
conda activate gatewizard
```

`environment.yml` already installs CPU MemPrO with the rest of the environment.

Check the version:

```bash
python -c "import gatewizard; print(gatewizard.__version__)"
```

## Dependencies

GateWizard is split across two repositories ([API](https://github.com/maurobedoya/gatewizard) + [GUI](https://github.com/franciscoadasme/gatewizard-gui)). Below is the **full dependency stack**; subsections mark what **this repo** installs automatically.

### Python — core (`pip install gatewizard`) · *this repo*

Installed automatically:

| Package | Role |
|---------|------|
| Python ≥ 3.8 | Runtime |
| NumPy | Numerical arrays |
| Matplotlib | Plots and analysis figures |
| MDAnalysis | Trajectories and topologies |
| lipyphilic | Membrane / lipid analysis (`AreaPerLipid`, leaflets, thickness; git extra for exclude until PyPI) |
| PROPKA | pKa and protonation |
| RDKit | Ligand 2D structures |
| Pillow | Image I/O |
| psique | Structure / sequence tools |
| gemmi | mmCIF and biological assembly |
| freud-analysis | Voronoi area-per-lipid (EVAPL, GridMAT, VTMC) |
| requests | HTTP helpers |

### Python — optional `[full]` · *this repo*

| Package | Role |
|---------|------|
| ParmEd | Topology conversion (GROMACS, etc.) |
| OpenMM | Equilibration / MD (Python; use with conda `cudatoolkit` for GPU) |
| PDBFixer | Short loops of 1–6 residues |
| MemPrO | Membrane orientation ([GitHub](https://github.com/pstansfeld/MemPrO)) |

Install: `pip install "gatewizard[full]"`. That extra includes the CPU MemPrO package. A checkout’s `environment.yml` already includes it.

For protein-aware **upstream** `AreaPerLipid` (`exclude_sel` before a PyPI release with PR #164): `pip install -r requirements-lipyphilic-git.txt`. GateWizard’s default APL with occupants remains **EVAPL**.

### Python — GUI backend · *[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)*

| Package | Role |
|---------|------|
| FastAPI | Local HTTP API for the desktop app |
| Uvicorn | ASGI server |
| Paramiko | SSH for Run on cluster |
| Biopython | Residue mutation |
| SciPy | Numerical helpers |
| Biotite | Chain superposition |
| PDBFixer | Short loops (also in `gatewizard[full]`) |
| gatewizard[full] | This API (installed into the embedded runtime) |

### Conda (building + OpenMM GPU) · *this repo via `environment.yml`; GUI embeds on Linux/WSL/macOS*

| Package | Role |
|---------|------|
| AmberTools 26 | `tleap`, `antechamber`, `packmol`, `pdb4amber`, … |
| reduce (bioconda) | Amber hydrogens for `pdb4amber --reduce`. AmberTools 26 does not ship it. |
| ParmEd 4.3.0 | From conda-forge (AmberTools-compatible) |
| OpenMM ≥ 8.0 | MD engine (Python) |
| cudatoolkit | OpenMM CUDA platform (Linux/WSL + NVIDIA driver) |
| GROMACS (CPU) | Equilibration. Prefer this conda package over a CUDA GROMACS build in the same env. |
| sqlite | Conda runtime library |

### Desktop app · *[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)*

| Package | Role |
|---------|------|
| Electron | Desktop shell |
| Node.js 20+ | App and build runtime |
| Svelte 5 | UI framework |
| Vite | Frontend bundler |
| Tailwind CSS | Styling |
| Three.js + Threlte | 3D structure viewer |
| electron-updater | In-app update checks |

### External MD engines (install separately · neither repo)

| Tool | Role |
|------|------|
| [NAMD](https://www.ks.uiuc.edu/Research/namd/) | Run NAMD equilibration — `namd3` / `namd2` on `PATH` |
| [GROMACS](https://www.gromacs.org/) | Run GROMACS equilibration — `gmx` on `PATH` (prefer conda-forge **CPU** `gromacs`; CUDA conda builds often hang next to OpenMM — see [Installation](docs/installation.md#gromacs-cpu-vs-cuda-conda-forge)) |
| [Amber](https://ambermd.org/) | Run Amber equilibration — `pmemd` / `pmemd.cuda`, or `sander` from AmberTools |

OpenMM needs no separate binary — use **`openmm` + `cudatoolkit`** (conda) with `gatewizard[full]`. Auto-selects **CUDA → OpenCL → CPU**; override with `PLATFORM=CUDA bash run_equilibration.sh`.

```bash
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"
```

### Provided by this repository

| Component | Included |
|-----------|----------|
| `gatewizard` Python package (core + `[full]` extras) | yes |
| `environment.yml` (AmberTools, OpenMM, cudatoolkit, …) | yes |
| Desktop app (Electron, Svelte, viewer) | no — [gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui) |
| FastAPI backend | no — GUI repo |
| NAMD, GROMACS, and Amber `pmemd` binaries | no — user install (`sander` comes with AmberTools) |

## Desktop GUI

**[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)** — Electron app for the same workflow with a visual interface. It manages its own Python runtime and can update the API from the app. First-launch bootstrap installs **CPU** GROMACS by default (CUDA GROMACS is opt-in); details in [Installation → Desktop GUI runtime](docs/installation.md#desktop-gui-runtime-gromacs--cuda).

## Community and support

| Need | Where |
|------|--------|
| Questions, workflows, methodology | [Discussions](https://github.com/maurobedoya/gatewizard/discussions) — central forum for **API and GUI** (use labels like `api`, `gui`, `api:installation`, `gui:installation`) |
| API bug or feature | [Issues](https://github.com/maurobedoya/gatewizard/issues/new/choose) in this repo |
| GUI bug or feature | [Issues](https://github.com/franciscoadasme/gatewizard-gui/issues/new/choose) in gatewizard-gui |
| Docs and troubleshooting | [Documentation](https://maurobedoya.github.io/gatewizard/) |
| Manual (PDF) | [Download](docs/manual/GateWizard-manual.pdf) |

New to the project? Read the [welcome post](https://github.com/maurobedoya/gatewizard/discussions/36) and say hello in the comments.

