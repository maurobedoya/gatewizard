# GateWizard

[![PyPI version](https://img.shields.io/pypi/v/gatewizard.svg)](https://pypi.org/project/gatewizard/)
[![Python](https://img.shields.io/pypi/pyversions/gatewizard.svg)](https://pypi.org/project/gatewizard/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://maurobedoya.github.io/gatewizard/)
[![DOI](https://zenodo.org/badge/1073861334.svg)](https://doi.org/10.5281/zenodo.18264074)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub stars](https://img.shields.io/github/stars/maurobedoya/gatewizard.svg)](https://github.com/maurobedoya/gatewizard/stargazers)

A library and GUI application tool for membrane protein preparation and molecular dynamics analysis.

📖 **[Read the Documentation](https://maurobedoya.github.io/gatewizard/)**

## Features

- **Protein Structure Preparation**: Clean PDB files, add missing atoms, and optimize structures
- **Propka Integration**: pKa calculations with automatic protonation state assignment
- **Protein Capping**: Add ACE/NME caps to protein termini before analysis
- **Force Field Support**: Compatible with Amber force fields (ff14SB, ff19SB, etc.)
- **Membrane System Building**: Automated membrane protein insertion and equilibration
- **Modern GUI**: Built with CustomTkinter for an intuitive user experience

## Installation

### Quick Installation from PyPI (Recommended)

1. **Create a conda environment with scientific dependencies:**
   ```bash
   conda create -n gatewizard -c conda-forge python sqlite ambertools=24 parmed=4.3.0 -y
   ```

2. **Activate the environment:**
   ```bash
   conda activate gatewizard
   ```

3. **Install GateWizard from PyPI:**
   ```bash
   pip install gatewizard
   ```

### Alternative: Development Installation

For developers or to install from source:

1. **Clone the repository:**
   ```bash
   git clone https://github.com/maurobedoya/gatewizard.git
   cd gatewizard
   ```

2. **Create environment from file:**
   ```bash
   conda env create -f environment.yml
   conda activate gatewizard
   ```

3. **Install in development mode:**
   ```bash
   pip install -e .
   ```

## Dependencies

### Core Dependencies (Automatically Installed)
- **Python** ≥ 3.8
- **CustomTkinter** ≥ 5.0.0 - Modern GUI framework
- **NumPy** ≥ 1.21.0 - Numerical computing
- **Matplotlib** ≥ 3.5.0 - Plotting and visualization
- **MDAnalysis** ≥ 2.0.0 - Molecular analysis toolkit
- **Propka** ≥ 3.2.0 - pKa calculations
- **BioPython** - PDB file handling

### Scientific Computing Dependencies (Via Conda)
- **AmberTools 24** - Molecular dynamics preparation and analysis
- **Parmed 4.3.0** - Parameter/topology file manipulation (must be from conda-forge for compatibility)

### External Requirements
- **NAMD 3.0.1** - Required for molecular dynamics simulations and equilibration
  - Download from: [NAMD Official Website](https://www.ks.uiuc.edu/Research/namd/)
  - Must be installed separately and accessible in your system PATH

## Usage

### Launch the GUI
```bash
gatewizard
```

### Command Line Options
```bash
gatewizard --help              # Show all options
gatewizard --screen 1          # Launch on secondary monitor
gatewizard --debug             # Enable debug logging
gatewizard --version           # Show version
```

## Upgrading

To upgrade to the latest version:

```bash
# Activate your environment
conda activate gatewizard

# Upgrade GateWizard
pip install --upgrade gatewizard
```

To check your current version:
```bash
gatewizard --version
```

## Troubleshooting

### Common Issues

**ImportError with numpy.compat:**
This indicates a version conflict between NumPy and Parmed. Make sure to install Parmed via conda-forge as shown in the installation instructions.

**pdb4amber command not found:**
Ensure AmberTools is installed via conda-forge and the gatewizard environment is activated.

**GUI not launching:**
Check that CustomTkinter is properly installed. Try reinstalling with `pip install --force-reinstall customtkinter`.

## Development

### Setting up for Development
```bash
# Clone the repository
git clone https://github.com/maurobedoya/gatewizard.git
cd gatewizard

# Create and activate environment
conda env create -f environment.yml
conda activate gatewizard

# Install in development mode
pip install -e .
```

### Project Structure

```
gatewizard/
├── gatewizard/          # Main source code
│   ├── gui/             # GUI components
│   ├── analysis/        # Analysis modules
│   └── ...
├── docs/                # Documentation (GitHub Pages)
├── tests/               # Test suite
├── environment.yml      # Conda environment
├── pyproject.toml      # Project configuration
└── README.md           # This file
```

### Running Tests

All tests are in the `tests/` directory. See [tests/README.md](tests/README.md) for detailed testing documentation.

```bash
# Run all tests
python -m pytest tests/

# Run with verbose output
python -m pytest tests/ -v

# Run specific test file
python -m pytest tests/test_propka_improvements.py

# Run with coverage
python -m pytest tests/ --cov=gatewizard --cov-report=html
```

### Documentation

Documentation is built with MkDocs and hosted on GitHub Pages.

```bash
# Install MkDocs
pip install mkdocs mkdocs-material

# Serve documentation locally
mkdocs serve

# Build documentation
mkdocs build

# Deploy to GitHub Pages
mkdocs gh-deploy
```

View documentation at: `http://localhost:8000` (when serving locally)

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Add/update tests
5. Update documentation
6. Commit your changes (`git commit -m 'Add amazing feature'`)
7. Push to the branch (`git push origin feature/amazing-feature`)
8. Open a Pull Request

### Code Style

- Follow PEP 8 guidelines
- Add docstrings to functions and classes
- Include type hints where appropriate
- Write tests for new features

## License

GateWizard is licensed under the [MIT License](LICENSE).

Copyright (c) 2025 Constanza González and Mauricio Bedoya

## Authors

- Constanza González
- Mauricio Bedoya

