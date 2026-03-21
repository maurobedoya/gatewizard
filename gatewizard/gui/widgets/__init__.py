"""
GUI widgets module for Gatewizard.

This module contains reusable GUI widgets and components.
"""

from gatewizard.gui.widgets.stage_tabs import StageTabsContainer
from gatewizard.gui.widgets.leaflet_frame import LeafletFrame
from gatewizard.gui.widgets.searchable_combobox import SearchableComboBox
from gatewizard.gui.widgets.progress_tracker import ProgressTracker
from gatewizard.gui.widgets.collapsible_section import CollapsibleSection

try:
    from gatewizard.gui.widgets.vtk_frame import VTKFrame
    VTK_AVAILABLE = True
except ImportError:
    VTK_AVAILABLE = False
    VTKFrame = None

__all__ = [
    "StageTabsContainer",
    "LeafletFrame",
    "SearchableComboBox", 
    "ProgressTracker",
    "CollapsibleSection",
    "VTKFrame",
    "VTK_AVAILABLE",
]