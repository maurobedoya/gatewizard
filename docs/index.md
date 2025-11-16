# GateWizard Documentation

Welcome to GateWizard, a comprehensive tool for membrane protein preparation and molecular dynamics analysis.

## Overview

GateWizard is a Python-based GUI application designed to streamline the workflow of preparing protein structures for molecular dynamics simulations, with a special focus on membrane proteins. It integrates various computational tools to provide a unified interface for structure preparation, analysis, and visualization.

## Key Features

- **Protein Structure Preparation**: Clean PDB files, add missing atoms, and optimize protein structures
- **Propka Integration**: Automatic pKa calculations with protonation state assignment
- **Protein Capping**: Add ACE/NME caps to protein termini for proper MD simulations
- **Force Field Support**: Compatible with Amber force fields (ff14SB, ff19SB, etc.)
- **Membrane System Building**: Automated membrane protein insertion and equilibration protocols
- **Trajectory Analysis**: RMSD, RMSF, distances, and radius of gyration calculations
- **Energy Analysis**: NAMD log file parsing and visualization
- **Modern GUI**: Built with CustomTkinter for an intuitive, cross-platform user experience

## Quick Start

### Installation

```bash
# Create conda environment with dependencies
conda create -n gatewizard -c conda-forge python sqlite ambertools=24 parmed=4.3.0 -y

# Activate environment
conda activate gatewizard

# Install GateWizard from PyPI
pip install gatewizard
```

### Launch

```bash
gatewizard
```

### Upgrade

```bash
pip install --upgrade gatewizard
```

For more detailed installation instructions, see the [Installation Guide](installation.md).

## Documentation Sections

- [Installation Guide](installation.md) - Detailed installation instructions for all platforms
- [User Guide](user-guide.md) - Complete guide to using GateWizard features
- [Analysis Features](analysis.md) - Trajectory and energy analysis capabilities
- [API Reference](api/index.md) - Developer documentation
- [Troubleshooting](troubleshooting.md) - Common issues and solutions

## System Requirements

- **Operating System**: Linux, macOS, or Windows (with WSL recommended)
- **Python**: 3.8 or higher
- **Conda**: For managing dependencies (Miniconda or Anaconda)
- **RAM**: Minimum 4GB (8GB+ recommended for large systems)
- **Disk Space**: ~2GB for installation plus space for your projects

## Getting Help

If you encounter issues or have questions:

1. Check the [Troubleshooting Guide](troubleshooting.md)
2. Review the [User Guide](user-guide.md) for usage instructions
3. Submit an issue on GitHub (if repository is public)

## Contributing

GateWizard is open-source software. Contributions are welcome! See our development guidelines for more information.

## License

GateWizard is licensed under the [MIT License](https://opensource.org/licenses/MIT).

Copyright (c) 2025 Constanza González and Mauricio Bedoya

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## Authors

- **Constanza González** - constanza.gonzalez.villagra@gmail.com
- **Mauricio Bedoya** - mbedoya@ucm.cl

## Citation

If you use GateWizard in your research, please cite:

```
[Citation information to be added]
```

---

*Last updated: October 2025*
