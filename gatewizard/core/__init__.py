"""
Core functionality module for Gatewizard.

This module contains the core business logic for protein preparation,
system building, and job monitoring. It is independent of the GUI
and can be used programmatically.
"""

from gatewizard.core.preparation import (
    run_propka,
    extract_summary_section,
    parse_summary_section,
    modify_pdb_based_on_summary,
)

from gatewizard.core.builder import Builder
from gatewizard.core.job_monitor import JobMonitor, JobStatus
from gatewizard.core.file_manager import FileManager

__all__ = [
    "run_propka",
    "extract_summary_section",
    "parse_summary_section", 
    "modify_pdb_based_on_summary",
    "Builder",
    "JobMonitor",
    "JobStatus",
    "FileManager",
]