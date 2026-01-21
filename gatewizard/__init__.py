# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Gatewizard - A tool for membrane protein preparation and analysis.
"""

__version__ = "1.0.14"
__author__ = "Constanza González, Mauricio Bedoya"
__email__ = ""
__license__ = "MIT"

# Import main classes for easier access
from gatewizard.core.preparation import (
    run_propka,
    extract_summary_section,
    parse_summary_section,
    modify_pdb_based_on_summary,
)

from gatewizard.core.builder import Builder
from gatewizard.core.job_monitor import JobMonitor

# GUI imports (optional, only if GUI dependencies are available)
try:
    from gatewizard.gui.app import ProteinViewerApp
    GUI_AVAILABLE = True
except ImportError:
    GUI_AVAILABLE = False
    ProteinViewerApp = None

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    "run_propka",
    "extract_summary_section", 
    "parse_summary_section",
    "modify_pdb_based_on_summary",
    "Builder",
    "JobMonitor",
    "GUI_AVAILABLE",
]

# Add ProteinViewerApp to __all__ if GUI is available
if GUI_AVAILABLE:
    __all__.append("ProteinViewerApp")