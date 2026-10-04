# API Reference Overview

Welcome to the GateWizard API documentation. This reference provides detailed information about all modules, classes, methods, and functions available in GateWizard.

## Available Modules

### [Method dictionary](dictionary.md)
A–Z index of public classes and functions. Click a name to jump to its heading.

### [Structure Manager](structure_manager.md)
Load, inspect, edit, and save PDB structures (selections, chain/residue edits, MemPrO orientation).

**Key Features:**

- Load from disk or PDB ID 
- MDAnalysis selections and criteria helpers  
- Rename / renumber / delete atoms 
- Secondary structure assignment 
- Apply a MemPrO orientation to the full system 

**Main Classes:** `StructureManager`

---

### [Biological assembly](assembly.md)
Expand the oligomer or capsid from operators already stored in a PDB or mmCIF file.

**Key Features:**

- List biological assemblies and their chain and atom counts
- Build an assembly as mmCIF with chain names such as `A1` and `A2`
- Leave out crystal waters unless requested
- Leave a file that is already the biological unit unchanged

---

### [Preparation Module](preparation.md)
Module for predicting pKa values and managing protonation states in protein structures.

**Key Features:**

- pKa prediction and analysis
- Protonation state assignment based on pH
- Disulfide bond detection and application
- Protein capping with ACE/NME groups
- pH-dependent protein structure preparation
- Apply Amber protonation names, then complete missing protein atoms (`complete_missing_heavy_atoms`, tleap templates including ASH/GLH protons); `pdb4amber --reduce` can rebuild hydrogens when `reduce` is installed
- Restore original residue numbers after `pdb4amber` (`preserve_residue_numbers`; ACE/NME get N-terminus − 1 / C-terminus + 1)

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

### [Hydration Module](hydration.md)
Module for cavity hydration with standalone PACKMOL (AmberTools TIP3P waters).

**Key Features:**

- PACKMOL availability check
- Hydrogen status detection (heavy-atom-safe mode)
- Cavity volume estimation inside a 3D box
- PACKMOL input generation and execution
- Custom PACKMOL input support

**Main Functions:** `check_packmol_available`, `estimate_cavity_volume`, `hydrate_cavity`

---

### [Equilibration Module](equilibration.md)
Module for managing equilibration protocols and workflows for NAMD, GROMACS, OpenMM, and Amber.

**Key Features:**

- NAMD equilibration protocol generation (CHARMM-GUI template integration)
- GROMACS equilibration protocol generation
- OpenMM equilibration protocol generation
- Amber equilibration protocol generation (`pmemd` / `sander` mdin + GROUP restraints)
- Flexible restraint system (protein backbone/sidechain, lipid head/tail, ligand, ions)
- MDAnalysis-based atom selection for restraint files
- Progressive force-constant schedule across stages
- Custom stage parameters via `EquilibrationStage` dataclass

**Main Classes:** `NAMDEquilibrationManager`, `GROMACSEquilibrationManager`, `OpenMMEquilibrationManager`, `AmberEquilibrationManager`, `EquilibrationStage`

---

### [Analysis Module](analysis.md)
Module for analyzing simulation results and monitoring equilibration progress.

**Key Features:**

- NAMD log file parsing and energy analysis
- OpenMM StateDataReporter log parsing
- MDAnalysis-based trajectory analysis (RMSD, RMSF, distances, radius of gyration)
- Lipid bilayer analysis: **FATSLiM** CLI default (`apl_method='fatslim'`); EVAPL experimental; official lipyphilic `AreaPerLipid` when `apl_method='lipyphilic'`; lipyphilic for leaflets and membrane thickness
- Dark-theme matplotlib plots
- Multi-stage log concatenation with per-file time override
- Comprehensive 2×2 energy summary plots

**Main Classes:** `EnergyAnalyzer` (`namd_analysis`), `TrajectoryAnalyzer` (`trajectory_analysis`), `BilayerTrajectoryAnalyzer` (`lipid_bilayer_analysis`), `OpenMMLogAnalyzer`

---

## Quick Links

- **[Quick Reference](../QUICK_REFERENCE.md)** - Common patterns and snippets
- **[User Guide](../user-guide.md)** - Step-by-step tutorials
- **[Method dictionary](dictionary.md)** - Clickable A–Z of public methods
- **[Example scripts](https://github.com/maurobedoya/gatewizard/tree/main/tests)** - `tests/<module>_examples/` (included in these pages and run by pytest)

## Getting Help

If you need additional help:

1. Check the [Troubleshooting Guide](../troubleshooting.md)
2. Review the [User Guide](../user-guide.md) for practical examples
3. Ask in [GitHub Discussions](https://github.com/maurobedoya/gatewizard/discussions) or report bugs via [issues](https://github.com/maurobedoya/gatewizard/issues/new/choose)
4. See complete workflow examples in the `tests/` example directories

## Module Organization

Each module documentation page includes:

- **Import statements** - How to import the module
- **Class/Function signatures** - Detailed method signatures
- **Parameters** - Complete parameter descriptions
- **Returns** - Return types and values
- **Examples** - Practical code examples
- **Workflows** - Complete usage patterns
