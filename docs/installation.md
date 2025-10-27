# Installation Guide

This guide provides detailed instructions for installing GateWizard on different platforms.

## Prerequisites

Before installing GateWizard, ensure you have:

- **Conda** (Miniconda or Anaconda): [Download here](https://docs.conda.io/en/latest/miniconda.html)
- **Git** (for cloning the repository): [Download here](https://git-scm.com/downloads)
- **Python 3.8+** (will be installed via conda)

## Quick Installation (Recommended)

This is the fastest way to get GateWizard up and running:

```bash
# 1. Create conda environment with all dependencies
conda create -n gatewizard -c conda-forge python sqlite ambertools=24 parmed=4.3.0 -y

# 2. Activate the environment
conda activate gatewizard

# 3. Clone the repository (if not already done)
git clone <repository-url>
cd gatewizard

# 4. Install GateWizard in development mode
pip install -e .

# 5. Launch GateWizard
gatewizard
```

## Alternative: Using environment.yml

If you prefer to use the provided environment file:

```bash
# 1. Clone the repository
git clone <repository-url>
cd gatewizard

# 2. Create environment from file
conda env create -f environment.yml

# 3. Activate the environment
conda activate gatewizard

# 4. Install GateWizard
pip install -e .

# 5. Launch GateWizard
gatewizard
```

## Platform-Specific Instructions

### Linux

GateWizard works natively on Linux. Follow the Quick Installation steps above.

**Additional Notes:**
- Ensure you have `tkinter` support: `sudo apt-get install python3-tk` (Ubuntu/Debian)
- For display issues, ensure X11 is properly configured

### macOS

GateWizard works on macOS with some considerations:

```bash
# Follow the Quick Installation steps, then:

# If you encounter display issues:
conda install -c conda-forge python.app

# Run using pythonw instead of python
pythonw -m gatewizard
```

**Known Issues:**
- On Apple Silicon (M1/M2), some dependencies may require Rosetta 2
- CustomTkinter may have scaling issues on Retina displays

### Windows

GateWizard is best used on Windows with WSL (Windows Subsystem for Linux):

**Option 1: WSL (Recommended)**
```bash
# 1. Install WSL2 with Ubuntu
wsl --install

# 2. Inside WSL, follow the Linux installation instructions
# 3. Install X server on Windows (e.g., VcXsrv, X410)
# 4. Set DISPLAY environment variable in WSL
export DISPLAY=:0
```

**Option 2: Native Windows (Experimental)**
```bash
# Use Anaconda Prompt
conda create -n gatewizard -c conda-forge python sqlite ambertools=24 parmed=4.3.0 -y
conda activate gatewizard
pip install -e .
gatewizard
```

## Dependencies

### Core Python Dependencies

Automatically installed via pip:

- **customtkinter** >= 5.0.0 - Modern GUI framework
- **numpy** >= 1.21.0 - Numerical computing
- **matplotlib** >= 3.5.0 - Plotting and visualization
- **MDAnalysis** >= 2.0.0 - Molecular analysis
- **propka** >= 3.2.0 - pKa calculations
- **biopython** - PDB file handling

### Scientific Computing Dependencies

Must be installed via conda:

- **AmberTools 24** - Molecular dynamics preparation and analysis
- **Parmed 4.3.0** - Parameter/topology manipulation (must be from conda-forge)

## Verifying Installation

After installation, verify that everything works:

```bash
# Activate environment
conda activate gatewizard

# Check Python version
python --version  # Should be 3.8+

# Check key dependencies
python -c "import customtkinter; print('CustomTkinter:', customtkinter.__version__)"
python -c "import MDAnalysis; print('MDAnalysis:', MDAnalysis.__version__)"
python -c "import parmed; print('ParmEd:', parmed.__version__)"

# Check AmberTools
which pdb4amber  # Should show path in conda environment

# Launch GateWizard
gatewizard --version
```

## Troubleshooting Installation

### Issue: ImportError with numpy.compat

**Cause:** Version conflict between NumPy and Parmed.

**Solution:**
```bash
# Reinstall Parmed from conda-forge
conda install -c conda-forge parmed=4.3.0 --force-reinstall
```

### Issue: pdb4amber command not found

**Cause:** AmberTools not properly installed or environment not activated.

**Solution:**
```bash
# Ensure environment is activated
conda activate gatewizard

# Reinstall AmberTools
conda install -c conda-forge ambertools=24
```

### Issue: CustomTkinter GUI not appearing

**Cause:** Display configuration or tkinter missing.

**Solution:**
```bash
# On Linux
sudo apt-get install python3-tk

# Reinstall CustomTkinter
pip install --force-reinstall customtkinter

# Check display
echo $DISPLAY  # Should show something like :0 or localhost:10.0
```

### Issue: Module not found errors

**Cause:** Installation incomplete or wrong environment.

**Solution:**
```bash
# Ensure you're in the correct environment
conda activate gatewizard

# Reinstall in development mode
cd /path/to/gatewizard
pip install -e . --force-reinstall
```

## Updating GateWizard

To update to the latest version:

```bash
# Update repository
cd /path/to/gatewizard
git pull origin main

# Update dependencies if needed
conda env update -f environment.yml

# Reinstall
pip install -e . --force-reinstall

# Restart GateWizard
gatewizard
```

## Uninstallation

To completely remove GateWizard:

```bash
# Remove conda environment
conda deactivate
conda env remove -n gatewizard

# Remove repository (if desired)
rm -rf /path/to/gatewizard
```

## Next Steps

After successful installation:

1. Read the [User Guide](user-guide.md) to learn about features
2. Try the example workflows
3. Check the [Troubleshooting Guide](troubleshooting.md) if you encounter issues

---

*For additional help, please contact the developers or check the GitHub repository.*
