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

**Main Classes:** `PropkaAnalyzer`, `ProteinCapper`

---

### [Builder Module](builder.md)
Module for building and preparing molecular dynamics simulation systems.

**Key Features:**
- System configuration and setup
- Integration with CHARMM-GUI
- Membrane protein system preparation

**Main Classes:** `SystemBuilder`

---

### [Equilibration Module](equilibration.md)
Module for managing equilibration protocols and workflows.

**Key Features:**
- NAMD equilibration protocol generation
- Configuration file management
- CHARMM-GUI template integration

**Main Classes:** `NAMDEquilibrationManager`

---

### [Analysis Module](analysis.md)
Module for analyzing simulation results and monitoring progress.

**Key Features:**
- Equilibration progress monitoring
- NAMD log file parsing
- Result analysis and reporting

**Main Functions:** `get_equilibration_progress()`, `parse_namd_log()`

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
