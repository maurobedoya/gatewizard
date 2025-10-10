# Gatewizard

A library and GUI application tool for membrane protein preparation and molecular dynamics analysis.

## Features

- **Protein Structure Preparation**: Clean PDB files, add missing atoms, and optimize structures
- **Propka Integration**: pKa calculations with automatic protonation state assignment
- **Protein Capping**: Add ACE/NME caps to protein termini before analysis
- **Force Field Support**: Compatible with Amber force fields (ff14SB, ff19SB, etc.)
- **Membrane System Building**: Automated membrane protein insertion and equilibration
- **Modern GUI**: Built with CustomTkinter for an intuitive user experience

## Installation

### Quick Installation (Recommended)


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

## Key Features

### Propka Analysis
- Automatic pKa calculations for protein residues
- Protonation state assignment at specified pH
- Optional protein capping with ACE/NME groups
- Export results for molecular dynamics simulations

### System Building
- Automated protein preparation workflow
- Integration with AmberTools for force field assignment
- Support for different membrane types
- Equilibration protocol generation
- Analysis of the equilibrated and production MDs

## Troubleshooting

### Common Issues


## Development

### Setting up for Development
```bash
git clone <repository>
cd gatewizard
conda env create -f environment.yml
conda activate gatewizard
pip install -e .
```


### Documentation

Documentation is and hosted on GitHub Pages.

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


## License

MIT

## Authors

- Constanza González 
- Mauricio Bedoya

