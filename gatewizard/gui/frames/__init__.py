"""
GUI frames module for Gatewizard.

This module contains the main content frames for different analysis stages.
"""

from gatewizard.gui.frames.visualize import VisualizeFrame
from gatewizard.gui.frames.preparation_frame import PreparationFrame
from gatewizard.gui.frames.builder_frame import BuilderFrame
from gatewizard.gui.frames.equilibration import EquilibrationFrame
from gatewizard.gui.frames.analysis import AnalysisFrame

__all__ = [
    "VisualizeFrame",
    "PreparationFrame", 
    "BuilderFrame",
    "EquilibrationFrame",
    "AnalysisFrame",
]