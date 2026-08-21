# Installation Guide

This guide provides detailed instructions for installing GateWizard on different platforms.

## Prerequisites

Before installing GateWizard, ensure you have:

- **Conda** (Miniconda or Anaconda): [Download here](https://docs.conda.io/en/latest/miniconda.html)
- **Git** (for cloning the repository): [Download here](https://git-scm.com/downloads)
- **Python 3.8+** (will be installed via conda)

## Quick Installation from PyPI (Recommended)

This is the fastest way to get GateWizard up and running:

```bash
# 1. Create conda environment with scientific dependencies
conda create -n gatewizard -c conda-forge python sqlite ambertools=24 parmed=4.3.0 -y

# 2. Activate the environment
conda activate gatewizard

# 3. Install GateWizard from PyPI
pip install gatewizard

# 4. Verify the API package
python -c "import gatewizard; print(gatewizard.__version__)"
```

For the **desktop app**, install [gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui/releases) separately. On first launch the GUI embeds its own micromamba runtime (Python, AmberTools, OpenMM, GROMACS, and `gatewizard` via pip). See [Desktop GUI runtime (GROMACS / CUDA)](#desktop-gui-runtime-gromacs--cuda) below.

## Alternative: Development Installation

For developers or to install from source:

```bash
# 1. Clone the repository
git clone https://github.com/maurobedoya/gatewizard.git
cd gatewizard

# 2. Create environment from file
conda env create -f environment.yml

# 3. Activate the environment
conda activate gatewizard

# 4. Install in development mode
pip install -e .

# 5. Verify the API package
python -c "import gatewizard; print(gatewizard.__version__)"
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

```

**Known Issues:**
- On Apple Silicon (M1/M2), some dependencies may require Rosetta 2

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
pip install gatewizard
python -c "import gatewizard; print(gatewizard.__version__)"
```

For the desktop app on Windows, use the [gatewizard-gui Windows installer](https://github.com/franciscoadasme/gatewizard-gui/releases) or the Linux/WSL build for the full MD workflow.

If Windows Task Manager shows `VmmemWSL` using several GB after you close the app, cap the VM — see [Troubleshooting: VmmemWSL RAM](troubleshooting.md#issue-vmmemwsl-keeps-using-a-lot-of-ram-after-gatewizard-closes).

## Dependencies

### Core Python Dependencies

Automatically installed via pip:

- **Python** ≥ 3.8
- **CustomTkinter** ≥ 5.0.0 - Modern GUI framework
- **NumPy** ≥ 1.21.0 - Numerical computing
- **Matplotlib** ≥ 3.5.0 - Plotting and visualization
- **MDAnalysis** ≥ 2.0.0 - Molecular analysis toolkit
- **Propka** ≥ 3.2.0 - pKa calculations
- **RDKit** ≥ 2023.3.1 - Ligand 2D structure visualization

### Scientific Computing Dependencies

Must be installed via conda:

- **AmberTools 24** - Molecular dynamics preparation and analysis
- **Parmed 4.3.0** - Parameter/topology manipulation (must be from conda-forge)

### External Requirements
- **NAMD** - Required for NAMD equilibration (separate install; not on conda-forge)
  - Recommended: **NAMD 3.0.1 or later** (`namd3`)
  - Download from: [NAMD Official Website](https://www.ks.uiuc.edu/Research/namd/)
  - Must be on PATH or selectable via the GUI engine picker
- **OpenMM** (optional for API users) - Python MD engine via conda-forge
  - `conda install -c conda-forge openmm cudatoolkit` on Linux/WSL with an NVIDIA driver for the CUDA platform
  - macOS uses Metal/OpenCL (no `cudatoolkit`)
  - Auto-selects **CUDA → OpenCL → CPU**; override with `PLATFORM=CUDA` when running equilibration scripts
- **GROMACS** (optional) - Recommended via conda-forge; the GUI also installs a build into its embedded runtime
  - **CPU (recommended default):** `conda install -c conda-forge gromacs`
  - **CUDA (Linux only, advanced):** `conda install -c conda-forge 'gromacs=*=nompi_cuda*'`
  - System installs under `/usr/local/gromacs` (with `GMXRC`) are auto-detected
  - See the notes below before choosing the CUDA conda build

### GROMACS: CPU vs CUDA (conda-forge)

GateWizard equilibration works with either a **CPU** or **CUDA** `gmx`. Prefer the CPU package unless you specifically need a conda CUDA GROMACS binary.

| Goal | Command / behavior |
|------|--------------------|
| Reliable conda install | `conda install -c conda-forge gromacs` |
| CUDA GROMACS (opt-in) | `conda install -c conda-forge 'gromacs=*=nompi_cuda*'` |
| GPU OpenMM (separate) | `conda install -c conda-forge openmm cudatoolkit` |

**Why CUDA GROMACS often fails or hangs**

- Installing `gromacs=*=nompi_cuda*` into an environment that already has **OpenMM + `cudatoolkit`** can stall for a very long time in the dependency **solver** (little or no log output after `Pinned packages:`).
- Log lines about the CUDA Toolkit EULA or Anaconda Terms of Service are **informational**. Non-interactive installs (`conda`/`micromamba` with `-y`) do **not** wait for you to type “yes”.
- If the CUDA solve is cancelled or times out, installing the **CPU** `gromacs` package normally succeeds quickly.
- For GPU MD with GateWizard, **OpenMM + `cudatoolkit`** is usually enough. Use a **system** CUDA GROMACS (or GMXRC) if you need `gmx` on the GPU, or keep conda GROMACS on CPU.

### Desktop GUI runtime (GROMACS / CUDA)

[gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui) bootstraps an embedded micromamba env on first launch (see that repo’s README and `runtime-install.log`).

- **OpenMM:** on Linux/WSL, if `nvidia-smi` sees a GPU, the GUI installs `openmm` + `cudatoolkit` into the embedded env (unless `GATEWIZARD_SKIP_CONDA_CUDA=1`).
- **GROMACS (default):** the GUI installs the **CPU** conda-forge build. This avoids multi-hour CUDA solver hangs next to OpenMM’s toolkit.
- **GROMACS CUDA (opt-in):** set `GATEWIZARD_CONDA_GROMACS_CUDA=1` before launching the GUI. The bootstrap tries a frozen CUDA install with a short timeout, then falls back to CPU GROMACS if the solve stalls.
- Quitting mid-install kills leftover `micromamba` processes for that env so a relaunch is not blocked by an orphaned solver.
- Linux/WSL log: `~/.config/gatewizard-gui/runtime-install.log` (look for `[gromacs] starting…` / `failed` / `installed`).

Manual API environments should follow the same preference: install **CPU** `gromacs` by default; only add the CUDA matchspec if you accept a long or fragile solve.
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

# Verify GateWizard API
python -c "import gatewizard; print(gatewizard.__version__)"
```

## Upgrading GateWizard

To upgrade to the latest version from PyPI:

```bash
# Activate environment
conda activate gatewizard

# Upgrade GateWizard
pip install --upgrade gatewizard

# Verify new version
python -c "import gatewizard; print(gatewizard.__version__)"
```

**Check for updates:**
You can check if a newer version is available:
```bash
pip list --outdated | grep gatewizard
```

**Force reinstall (if needed):**
If you encounter issues after an update:
```bash
pip install --force-reinstall gatewizard
```

## Troubleshooting Installation

### Issue: conda GROMACS CUDA install hangs after “Pinned packages”

**Cause:** The `gromacs=*=nompi_cuda*` solve often conflicts with an environment that already has OpenMM/`cudatoolkit`. The solver can run for a very long time with almost no new log lines. EULA / Terms-of-Service messages in the log are not interactive prompts.

**Solution:**
```bash
# Cancel the stuck install (Ctrl+C), then install CPU GROMACS:
conda install -c conda-forge gromacs -y

# Optional: use a system CUDA GROMACS / GMXRC instead of the conda CUDA build
```

For **gatewizard-gui**, leave the default (CPU GROMACS). Only set `GATEWIZARD_CONDA_GROMACS_CUDA=1` if you explicitly want the timed CUDA attempt. See [Desktop GUI runtime (GROMACS / CUDA)](installation.md#desktop-gui-runtime-gromacs--cuda).

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

**Tip:** GateWizard follows [semantic versioning](https://semver.org/). Check the [releases page](https://github.com/maurobedoya/gatewizard/releases) for changelog and new features.

### For Development Installation

If you installed from source:

```bash
# Update repository
cd /path/to/gatewizard
git pull origin main

# Update dependencies if needed
conda env update -f environment.yml

# Reinstall
pip install -e . --force-reinstall

# Re-verify import
python -c "import gatewizard; print(gatewizard.__version__)"
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
