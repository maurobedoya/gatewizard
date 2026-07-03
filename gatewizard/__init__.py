# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Gatewizard - A tool for membrane protein preparation and analysis.
"""

__version__ = "1.0.44"
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
from gatewizard.core.structure_manager import StructureManager
from gatewizard.core.mempro import MemPrO
from gatewizard.tools.equilibration import (
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
    "Builder",
    "JobMonitor",
    "StructureManager",
    "MemPrO",
    "NAMDEquilibrationManager",
    "OpenMMEquilibrationManager",
    "OpenMMLogAnalyzer",
]
