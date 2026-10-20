# Installation Guide

This guide provides detailed instructions for installing GateWizard on different platforms.

## Prerequisites

Before installing GateWizard, ensure you have:

- **Conda** (Miniconda or Anaconda): [Download here](https://docs.conda.io/en/latest/miniconda.html)
- **Git** (for cloning the repository): [Download here](https://git-scm.com/downloads)
- **Python 3.8+** (will be installed via conda)

## Full API installation (Linux or WSL)

Use this when the new machine should build membranes, orient with MemPrO, and run OpenMM or GROMACS. NAMD stays a separate download.

```bash
conda create -n gatewizard -c conda-forge python sqlite ambertools=26 parmed=4.3.0 openmm cudatoolkit gromacs -y
conda activate gatewizard
conda install -c conda-forge -c bioconda reduce -y
pip install "gatewizard[full]"
python -c "import gatewizard; print(gatewizard.__version__)"
which pdb4amber
which mempro
```

On macOS, omit `cudatoolkit`. Prefer the CPU `gromacs` package. A CUDA GROMACS build in this same environment often stalls in the solver; see [GROMACS: CPU vs CUDA](#gromacs-cpu-vs-cuda-conda-forge).

AmberTools 26 does not ship `reduce`. `conda install -c bioconda reduce` puts that binary in the same environment so `pdb4amber --reduce` can add Amber hydrogens.

The same stack from a git checkout:

```bash
git clone https://github.com/maurobedoya/gatewizard.git
cd gatewizard
conda env create -f environment.yml
conda activate gatewizard
```

`environment.yml` already installs the package with `pip install -e ".[full]"`, including CPU MemPrO (`mempro` ≥ 0.3.0). `pip install "gatewizard[full]"` installs that same CPU package. A GPU build is optional and goes into this conda env only; see [MemPrO on a GPU](#mempro-on-a-gpu).

For the **desktop app**, install [gatewizard-gui](https://github.com/franciscoadasme/gatewizard-gui/releases) separately. On first launch the GUI embeds its own micromamba runtime. That runtime is not this conda env. See [Desktop GUI runtime (GROMACS / CUDA)](#desktop-gui-runtime-gromacs-cuda) below.

## API package only

This installs the Python library without OpenMM or MemPrO:

```bash
conda create -n gatewizard -c conda-forge python sqlite ambertools=26 parmed=4.3.0 -y
conda activate gatewizard
conda install -c conda-forge -c bioconda reduce -y
pip install gatewizard
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
conda create -n gatewizard -c conda-forge python sqlite ambertools=26 parmed=4.3.0 -y
conda activate gatewizard
conda install -c conda-forge -c bioconda reduce -y
pip install gatewizard
python -c "import gatewizard; print(gatewizard.__version__)"
```

For the desktop app on Windows, use the [gatewizard-gui Windows installer](https://github.com/franciscoadasme/gatewizard-gui/releases) or the Linux/WSL build for the full MD workflow.

If Windows Task Manager shows `VmmemWSL` using several GB after you close the app, cap the VM — see [Troubleshooting: VmmemWSL RAM](troubleshooting.md#issue-vmmemwsl-keeps-using-a-lot-of-ram-after-gatewizard-closes).

## Dependencies

### Core Python Dependencies

Installed by `pip install gatewizard`:

- **Python** ≥ 3.8
- **NumPy** ≥ 1.21.0
- **Matplotlib** ≥ 3.5.0
- **requests** ≥ 2.25.0
- **MDAnalysis** ≥ 2.0.0
- **lipyphilic** ≥ 0.12.0
- **freud-analysis** ≥ 2.5.0
- **Propka** ≥ 3.2.0
- **Pillow** ≥ 9.0.0
- **RDKit** ≥ 2023.3.1
- **psique** ≥ 1.2.1
- **gemmi** ≥ 0.6.5

`pip install "gatewizard[full]"` also installs ParmEd, OpenMM, PDBFixer, and the CPU build of MemPrO. AmberTools, `reduce`, the conda ParmEd build, OpenMM’s CUDA toolkit, and CPU GROMACS come from conda, not from that pip extra.

### Scientific Computing Dependencies

Must be installed via conda:

- **AmberTools 26** - Molecular dynamics preparation and analysis
- **reduce** (bioconda) - Amber hydrogens for `pdb4amber --reduce`. AmberTools 26 does not ship this binary.
- **Parmed 4.3.0** - Parameter/topology manipulation (must be from conda-forge)
- **PDBFixer** ≥ 1.9 - short loops of 1–6 residues. Installed by `gatewizard[full]`, not by `pip install gatewizard` alone.

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
- The embedded env is `~/.config/gatewizard-gui/mamba-env`. `sudo apt remove gatewizard-gui-linux` leaves it. `sudo apt purge gatewizard-gui-linux` deletes that directory for the user who ran sudo. It does not remove a separate conda env such as `~/miniconda3/envs/gatewizard`.

Manual API environments should follow the same preference: install **CPU** `gromacs` by default; only add the CUDA matchspec if you accept a long or fragile solve.

### MemPrO on a GPU

`gatewizard[full]` and the desktop app both install the CPU package (`mempro` with `jax[cpu]`). Neither installs `mempro[gpu]`. The API conda env and the GUI micromamba env are separate Pythons. Installing the GPU extra in one does not change the other. The extra is `jax[cuda12]`, for Linux with CUDA 12 drivers. macOS, Windows, and a machine with no NVIDIA driver stay on CPU. Asking for a GPU still fails when the card is not visible, including WSL without GPU passthrough.

**API.** Use the conda env from this page, not a `pip` that belongs to another env:

```bash
conda activate gatewizard
python -m pip install "mempro[gpu]"
```

`MemPrO.run(..., device="gpu")` then uses that env. `device="cpu"` and `device="auto"` still run on CPU when JAX sees no GPU.

**Desktop app.** Quit GateWizard first. The GPU choice in Visualize reads JAX inside the embedded env, not the conda env above. On Linux or WSL:

```bash
~/.config/gatewizard-gui/mamba-env/bin/python -m pip install "mempro[gpu]"
```

On macOS:

```bash
~/Library/gatewizard-gui/mamba-env/bin/python -m pip install "mempro[gpu]"
```

Start the app again. The GPU choice stays off until that env’s JAX can see a GPU. There is no native Windows installer; on a Windows PC this is the WSL command above, inside the Linux distro, not Windows Python and not `%APPDATA%\gatewizard-gui\mamba-env`.

## Verifying Installation

After installation, verify that everything works:

```bash
# Activate environment
conda activate gatewizard

# Check Python version
python --version  # Should be 3.8+

# Check key dependencies
python -c "import MDAnalysis; print('MDAnalysis:', MDAnalysis.__version__)"
python -c "import parmed; print('ParmEd:', parmed.__version__)"
python -c "import openmm; print([openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])"

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

For **gatewizard-gui**, leave the default (CPU GROMACS). Only set `GATEWIZARD_CONDA_GROMACS_CUDA=1` if you explicitly want the timed CUDA attempt. See [Desktop GUI runtime (GROMACS / CUDA)](installation.md#desktop-gui-runtime-gromacs-cuda).

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
conda install -c conda-forge ambertools=26
conda install -c conda-forge -c bioconda reduce
```

### Issue: Module not found errors

**Cause:** Installation incomplete or wrong environment.

**Solution:**
```bash
# Ensure you're in the correct environment
conda activate gatewizard

# Reinstall the API
cd /path/to/gatewizard
pip install -e ".[full]" --force-reinstall
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
pip install -e ".[full]" --force-reinstall

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
