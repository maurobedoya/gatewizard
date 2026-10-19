#!/usr/bin/env bash
# Print the documentation to one PDF.
# Use the environment that has requirements-docs.txt installed.
set -eo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
exec python scripts/build_docs_pdf.py "$@"
