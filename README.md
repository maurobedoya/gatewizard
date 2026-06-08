# GateWizard (API)
<img src="./resources/gatewizard_logo_black_back-white.png" alt="GateWizard" width="140" align="right" />

[![PyPI version](https://img.shields.io/pypi/v/gatewizard.svg)](https://pypi.org/project/gatewizard/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://maurobedoya.github.io/gatewizard/)
[![DOI](https://zenodo.org/badge/1073861334.svg)](https://doi.org/10.5281/zenodo.18264074)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/maurobedoya/gatewizard.svg)](https://github.com/maurobedoya/gatewizard/stargazers)

Python library for membrane protein preparation, system building, equilibration setup, and MD trajectory analysis.

**This repository is the API only.** It does not include a graphical interface.

If you want a desktop app, use **[gatewizard-gui](https://github.com/maurobedoya/gatewizard-gui)**, a separate project that calls this library. You can use the API without the GUI, and the GUI without cloning this repo (it installs the API via pip).

📖 [Documentation](https://maurobedoya.github.io/gatewizard/)

<br clear="left" />

## Features

- **Preparation** — clean structures, pKa (PROPKA), protonation, and termini capping
- **Membrane builder** — orient protein, pack lipids and solvent, Amber parametrization (tleap)
- **Equilibration** — CHARMM-GUI-style protocols for NAMD, GROMACS, and OpenMM (NVT, NPT, NPAT, NPγT)
- **Analysis** — trajectory metrics, structural analysis, and plotting helpers
- **Force fields** — Amber protein models (e.g. ff14SB, ff19SB) and common water/lipid setups

## Install

```bash
conda create -n gatewizard -c conda-forge python sqlite ambertools=24 parmed=4.3.0 openmm cudatoolkit -y
conda activate gatewizard
pip install "gatewizard[full]"
pip install git+https://github.com/pstansfeld/MemPrO.git
```

For optional extras without OpenMM (NAMD/GROMACS only):

```bash
pip install gatewizard
```

From source (includes `openmm` + `cudatoolkit` via conda):

```bash
git clone https://github.com/maurobedoya/gatewizard.git
cd gatewizard
conda env create -f environment.yml
conda activate gatewizard
pip install -e .
```

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
| lipyphilic | Membrane / lipid analysis |
| PROPKA | pKa and protonation |
| RDKit | Ligand 2D structures |
| Pillow | Image I/O |
| psique | Structure / sequence tools |
| requests | HTTP helpers |

### Python — optional `[full]` · *this repo*

| Package | Role |
|---------|------|
| ParmEd | Topology conversion (GROMACS, etc.) |
| OpenMM | Equilibration / MD (Python; use with conda `cudatoolkit` for GPU) |

Install: `pip install "gatewizard[full]"`

Install MemPrO separately for membrane orientation:

`pip install git+https://github.com/pstansfeld/MemPrO.git`

### Python — GUI backend · *[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)*

| Package | Role |
|---------|------|
| FastAPI | Local HTTP API for the desktop app |
| Uvicorn | ASGI server |
| gatewizard[full] | This API (installed into the embedded runtime) |

### Conda (building + OpenMM GPU) · *this repo via `environment.yml`; GUI embeds on Linux/WSL/macOS*

| Package | Role |
|---------|------|
| AmberTools 24 | `tleap`, `antechamber`, `packmol`, `pdb4amber`, … |
| ParmEd 4.3.0 | From conda-forge (AmberTools-compatible) |
| OpenMM ≥ 8.0 | MD engine (Python) |
| cudatoolkit | OpenMM CUDA platform (Linux/WSL + NVIDIA driver) |

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
| [GROMACS](https://www.gromacs.org/) | Run GROMACS equilibration — `gmx` on `PATH` |

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
| NAMD / GROMACS binaries | no — user install |

## Desktop GUI

**[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)** — Electron app for the same workflow with a visual interface. It manages its own Python runtime and can update the API from the app.
