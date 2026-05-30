#!/usr/bin/env python3

"""
Gatewizard package entry point.

The graphical user interface has moved to the gatewizard-gui repository.
This package provides the Python API only.
"""

import sys

from gatewizard import __version__


def main():
    """Entry point — prints a notice that the GUI has moved."""
    if len(sys.argv) > 1 and sys.argv[1] in ("-V", "--version"):
        print(f"Gatewizard {__version__}")
        sys.exit(0)
    print(
        f"Gatewizard {__version__} — API-only package.\n"
        "The graphical user interface has moved to the gatewizard-gui repository."
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
