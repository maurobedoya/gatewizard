# GateWizard (API)
<img src="./docs/images/splash.png" alt="GateWizard" width="140" align="right" />

[![PyPI version](https://img.shields.io/pypi/v/gatewizard.svg)](https://pypi.org/project/gatewizard/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://maurobedoya.github.io/gatewizard/)
[![DOI](https://zenodo.org/badge/1073861334.svg)](https://doi.org/10.5281/zenodo.18264074)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/maurobedoya/gatewizard.svg)](https://github.com/maurobedoya/gatewizard/stargazers)

Python library for membrane protein preparation, system building, equilibration setup, and MD trajectory analysis.

**This repository is the API only.** For a desktop app, use **[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui)** (it installs this library via pip).

📖 [Documentation](https://maurobedoya.github.io/gatewizard/) · [Start with the API](https://maurobedoya.github.io/gatewizard/QUICK_REFERENCE/) · [Roadmap](docs/roadmap.md) · [Manual (PDF)](docs/manual/GateWizard-manual.pdf) · 💬 [Discussions](https://github.com/maurobedoya/gatewizard/discussions)

## Features

- **Preparation** — clean structures, pKa (PROPKA), protonation, and termini capping
- **Membrane builder** — orient protein, pack lipids and solvent, Amber parametrization. Lipid catalog: [builder reference](https://maurobedoya.github.io/gatewizard/api/builder/#available-lipids)
- **Equilibration** — NAMD, GROMACS, OpenMM, and Amber (NVT, NPT, NPAT, NPγT)
- **Analysis** — trajectory metrics, structural analysis, and plotting helpers
- **Force fields** — Amber protein models (e.g. ff14SB, ff19SB) and common water/lipid setups

## Install

```bash
conda create -n gatewizard -c conda-forge python sqlite ambertools=26 parmed=4.3.0 openmm cudatoolkit gromacs -y
conda activate gatewizard
conda install -c conda-forge -c bioconda reduce -y
pip install "gatewizard[full]"
```

That line is for Linux or WSL with an NVIDIA driver. On macOS, leave out `cudatoolkit`. Full steps, `environment.yml`, and the library-only install: [Installation](https://maurobedoya.github.io/gatewizard/installation/).

New to the API? Start with the [quick reference](https://maurobedoya.github.io/gatewizard/QUICK_REFERENCE/). It downloads PDB 6RV2 and runs the first calls.

## Dependencies

`pip install gatewizard` installs the library. `pip install "gatewizard[full]"` adds membrane building and OpenMM. AmberTools, `reduce`, and GROMACS come from conda. Versions: [Installation → Dependencies](https://maurobedoya.github.io/gatewizard/installation/#dependencies).

For protein-aware upstream `AreaPerLipid` (`exclude_sel` before a PyPI release with PR #164), see [Analysis → Bilayer](https://maurobedoya.github.io/gatewizard/analysis/#bilayer-analysis-lipyphilic--fatslim).

## Community and support

| Need | Where |
|------|--------|
| Questions, workflows, methodology | [Discussions](https://github.com/maurobedoya/gatewizard/discussions) — API and GUI (labels `api`, `gui`, `api:installation`, `gui:installation`) |
| API bug or feature | [Issues](https://github.com/maurobedoya/gatewizard/issues/new/choose) in this repo |
| GUI bug or feature | [Issues](https://github.com/franciscoadasme/gatewizard-gui/issues/new/choose) in gatewizard-gui |
| Docs and troubleshooting | [Documentation](https://maurobedoya.github.io/gatewizard/) |
| Contributing | [CONTRIBUTING.md](CONTRIBUTING.md) |
| Manual (PDF) | [Download](docs/manual/GateWizard-manual.pdf) |

New to the project? Read the [welcome post](https://github.com/maurobedoya/gatewizard/discussions/36) and say hello in the comments.

## Citation

If you use GateWizard in your research, please cite:

- **Software:** [https://doi.org/10.5281/zenodo.18264074](https://doi.org/10.5281/zenodo.18264074)
- **Publication:** to be added

## License

GateWizard is licensed under the [MIT License](LICENSE).

Copyright (c) 2025-2026 Constanza González, Francisco Adasme, Mauricio González, Janin Riedelsberger, and Mauricio Bedoya
