"""
GUI module for Gatewizard.

This module contains the graphical user interface components,
including the main application window, frames, and widgets.
"""

try:
    from gatewizard.gui.app import ProteinViewerApp
    GUI_AVAILABLE = True
except ImportError:
    ProteinViewerApp = None
    GUI_AVAILABLE = False

__all__ = [
    "ProteinViewerApp",
    "GUI_AVAILABLE",
]
