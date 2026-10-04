# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Gatewizard - A tool for membrane protein preparation and analysis.
"""

__version__ = "1.0.53"
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
from gatewizard.core.titration import render_titration_png, run_titration, titration_figure
from gatewizard.core.assembly import build_assembly, list_assemblies

from gatewizard.core.builder import Builder
from gatewizard.core.job_monitor import JobMonitor
from gatewizard.core.structure_manager import StructureManager
from gatewizard.core.mempro import MemPrO
from gatewizard.tools.equilibration import (
    AmberEquilibrationManager,
    NAMDEquilibrationManager,
    OpenMMEquilibrationManager,
)
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer

__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "__license__",
    "run_propka",
    "extract_summary_section",
    "parse_summary_section",
    "modify_pdb_based_on_summary",
    "render_titration_png",
    "run_titration",
    "titration_figure",
    "list_assemblies",
    "build_assembly",
    "Builder",
    "JobMonitor",
    "StructureManager",
    "MemPrO",
    "AmberEquilibrationManager",
    "NAMDEquilibrationManager",
    "OpenMMEquilibrationManager",
    "OpenMMLogAnalyzer",
]
