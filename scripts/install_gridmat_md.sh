#!/usr/bin/env bash
# Clone original GridMAT-MD.pl for GateWizard apl_method='gridmat_md'.
#
# GridMAT-MD is GPLv3 — GateWizard (MIT) never vendors the .pl into the package.
# This script clones https://github.com/jalemkul/gridmat-md and prints the
# GATEWIZARD_GRIDMAT_MD export to use with the backend / API.
#
# Usage (from anywhere):
#   bash /path/to/gatewizard/scripts/install_gridmat_md.sh
#   # optional install dir:
#   GRIDMAT_MD_DIR=$HOME/src/gridmat-md bash scripts/install_gridmat_md.sh
#
# Requires: git, perl (on PATH when running analysis).

set -euo pipefail

DEST="${GRIDMAT_MD_DIR:-${HOME}/gridmat-md}"
REPO_URL="${GRIDMAT_MD_REPO:-https://github.com/jalemkul/gridmat-md.git}"

if ! command -v git >/dev/null 2>&1; then
  echo "git not found on PATH." >&2
  exit 1
fi

if [[ -d "${DEST}/.git" ]]; then
  echo "Updating existing clone at ${DEST}"
  git -C "${DEST}" pull --ff-only || true
elif [[ -e "${DEST}" ]]; then
  echo "Path exists but is not a git clone: ${DEST}" >&2
  exit 1
else
  echo "Cloning ${REPO_URL} → ${DEST}"
  git clone --depth 1 "${REPO_URL}" "${DEST}"
fi

SCRIPT="${DEST}/GridMAT-MD.pl"
if [[ ! -f "${SCRIPT}" ]]; then
  echo "Clone succeeded but GridMAT-MD.pl not found at ${SCRIPT}" >&2
  exit 1
fi

if ! command -v perl >/dev/null 2>&1; then
  echo "Warning: perl not on PATH — install perl before running gridmat_md." >&2
fi

echo
echo "GridMAT-MD.pl ready."
echo "  export GATEWIZARD_GRIDMAT_MD=\"${SCRIPT}\""
echo "Optional: export GATEWIZARD_PERL=\"\$(command -v perl)\""
echo "This installer lives in the gatewizard API source tree (GitHub clone), not in a GUI install."
echo "See docs/analysis.md (gatewizard API repo)."
