# API Reference Overview

Welcome to the GateWizard API documentation. This reference provides detailed information about all modules, classes, methods, and functions available in GateWizard.

## Available Modules

### [Preparation Module](preparation.md)
Module for predicting pKa values and managing protonation states in protein structures.

**Key Features:**
- pKa prediction and analysis
- Protonation state assignment based on pH
- Disulfide bond detection and application
- Protein capping with ACE/NME groups
- pH-dependent protein structure preparation

**Main Classes:** `PreparationManager`, `ProteinCapper`

---

### [Builder Module](builder.md)
Module for building and preparing molecular dynamics simulation systems.

**Key Features:**
- System configuration and setup
- Integration with CHARMM-GUI
- Membrane protein system preparation

**Main Classes:** `Builder`

---

### [MemPrO Module](mempro.md)
Module for orienting membrane proteins using MemPrO.

**Key Features:**
- Membrane protein orientation
- Ranked orientation results with scoring
- Oriented PDB file access by rank
- Command-line builder

**Main Classes:** `MemPrO`, `OrientationResult`

---

### [Equilibration Module](equilibration.md)
Module for managing equilibration protocols and workflows for NAMD, GROMACS, and OpenMM.

**Key Features:**
- NAMD equilibration protocol generation (CHARMM-GUI template integration)
- GROMACS equilibration protocol generation
- OpenMM equilibration protocol generation
- Flexible restraint system (protein backbone/sidechain, lipid head/tail, ligand, ions)
- MDAnalysis-based atom selection for restraint files
- Progressive force-constant schedule across stages
- Custom stage parameters via `EquilibrationStage` dataclass

**Main Classes:** `NAMDEquilibrationManager`, `GROMACSEquilibrationManager`, `OpenMMEquilibrationManager`, `EquilibrationStage`

---

### [Analysis Module](analysis.md)
Module for analyzing simulation results and monitoring equilibration progress.

**Key Features:**
- NAMD log file parsing and energy analysis
- OpenMM StateDataReporter log parsing
- MDAnalysis-based trajectory analysis (RMSD, RMSF, distances, radius of gyration)
- Lipid bilayer analysis via lipyphilic (area per lipid, membrane thickness)
- Dark-theme matplotlib plots
- Multi-stage log concatenation with per-file time override
- Comprehensive 2×2 energy summary plots

**Main Classes:** `EnergyAnalyzer`, `TrajectoryAnalyzer`, `BilayerTrajectoryAnalyzer`, `OpenMMLogAnalyzer`

---

## Quick Links

- **[Quick Reference](../QUICK_REFERENCE.md)** - Common patterns and snippets
- **[User Guide](../user-guide.md)** - Step-by-step tutorials
- **[Examples Directory](https://github.com/yourusername/gatewizard/tree/main/examples)** - Complete workflow examples

## Getting Help

If you need additional help:

1. Check the [Troubleshooting Guide](../troubleshooting.md)
2. Review the [User Guide](../user-guide.md) for practical examples
3. See complete workflow examples in the `examples/` directory

## Module Organization

Each module documentation page includes:

- **Import statements** - How to import the module
- **Class/Function signatures** - Detailed method signatures
- **Parameters** - Complete parameter descriptions
- **Returns** - Return types and values
- **Examples** - Practical code examples
- **Workflows** - Complete usage patterns
