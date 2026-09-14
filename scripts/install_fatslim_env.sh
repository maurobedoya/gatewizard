#!/usr/bin/env bash
# Install a companion conda env for the archived FATSLiM CLI (GateWizard default APL).
#
# FATSLiM is GPLv3 and needs Python ≤3.8 — do NOT install it into GateWizard's
# Python (≥3.9) or the GUI embedded mamba-env. GateWizard calls `fatslim` via
# subprocess only (MIT-safe).
#
# Usage:
#   bash scripts/install_fatslim_env.sh
#   export GATEWIZARD_FATSLIM="$(conda run -n fatslim-py38 which fatslim)"
#
# Platforms: Linux and WSL (primary); macOS best-effort (needs a C compiler /
# OpenMP for the Cython extensions).

set -euo pipefail

ENV_NAME="${FATSLIM_CONDA_ENV:-fatslim-py38}"

if ! command -v conda >/dev/null 2>&1; then
  echo "conda not found on PATH. Install Miniconda/Mambaforge first." >&2
  exit 1
fi

# shellcheck disable=SC1091
if [[ -f "${CONDA_PREFIX:-}/../etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "${CONDA_PREFIX}/../etc/profile.d/conda.sh"
elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "${HOME}/miniconda3/etc/profile.d/conda.sh"
elif [[ -f "${HOME}/mambaforge/etc/profile.d/conda.sh" ]]; then
  # shellcheck source=/dev/null
  source "${HOME}/mambaforge/etc/profile.d/conda.sh"
fi

if conda env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  echo "Conda env ${ENV_NAME} already exists — updating packages."
else
  echo "Creating conda env ${ENV_NAME} (Python 3.8)..."
  conda create -n "${ENV_NAME}" python=3.8 -y
fi

# shellcheck disable=SC1091
eval "$(conda shell.bash hook)"
conda activate "${ENV_NAME}"

pip install --upgrade pip
pip install 'numpy<1.24' 'cython<3'
pip install fatslim
# Pin so MDAnalysis does not pull Python≥3.9-only wheels into this env
pip install 'gsd<3.0' 'MDAnalysis==2.4.3'

FATSLIM_BIN="$(command -v fatslim)"
echo
echo "FATSLiM installed at: ${FATSLIM_BIN}"
fatslim --version || true
echo
echo "For GateWizard API / GUI (WSL backend), export:"
echo "  export GATEWIZARD_FATSLIM=\"${FATSLIM_BIN}\""
echo
echo "Or leave GATEWIZARD_FATSLIM unset if ${ENV_NAME}/bin is on PATH / discoverable."
echo "This installer lives in the gatewizard API source tree (GitHub clone), not in a GUI /usr/bin install."
echo "Docs: docs/analysis.md (FATSLiM section)."
