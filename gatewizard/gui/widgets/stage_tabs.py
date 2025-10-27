# gatewizard/gui/widgets/stage_tabs.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Stage tabs widget for navigation between analysis stages.

This module provides a tab-like interface for switching between
different analysis stages in the Gatewizard application.
"""

from typing import Callable, Optional, List

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, BUTTON_COLORS, PREPARATION_STAGES, WIDGET_SIZES
)
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class StageTabsContainer(ctk.CTkFrame):
    def set_stage_names(self, stage_names):
        """Set the names of the stages dynamically and update buttons."""
        self.stages = stage_names
        # Remove old buttons
        for btn in self.stage_buttons.values():
            btn.destroy()
        self.stage_buttons = {}
        # Create new buttons
        for stage in self.stages:
            button = ctk.CTkButton(
                self.tabs_frame,
                text=stage,
                width=WIDGET_SIZES['button_width'] + 20,
                height=WIDGET_SIZES['button_height'] + 10,
                font=FONTS['body'],
                fg_color=BUTTON_COLORS['inactive'],
                hover_color=BUTTON_COLORS['hover_inactive'],
                command=lambda s=stage: self._on_stage_clicked(s)
            )
            self.stage_buttons[stage] = button
        self._setup_layout()

    """
    Container widget for stage navigation tabs.
    
    This widget provides a tab-like interface for switching between
    different analysis stages with visual feedback for the active stage.
    """
    
    def __init__(
        self,
        parent,
        stage_changed_callback: Optional[Callable[[str], None]] = None,
        stages: Optional[List[str]] = None
    ):
        """
        Initialize the stage tabs container.
        
        Args:
            parent: Parent widget
            stage_changed_callback: Callback function called when stage changes
            stages: List of stage names (uses default if None)
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['canvas'])
        
        self.stage_changed_callback = stage_changed_callback
        self.stages = stages or ["Visualize", "Preparation", "Builder", "Equilibration", "Analysis"]
        
        # State variables
        self.active_stage = None
        self.stage_buttons = {}
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
    
    def _create_widgets(self):
        """Create stage tab buttons."""
        self.tabs_frame = ctk.CTkFrame(self, fg_color="transparent")
        
        # Create title label
        self.title_label = ctk.CTkLabel(
            self.tabs_frame,
            text="Analysis Stages:",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Create stage buttons
        for i, stage in enumerate(self.stages):
            button = ctk.CTkButton(
                self.tabs_frame,
                text=stage,
                width=WIDGET_SIZES['button_width'] + 20,
                height=WIDGET_SIZES['button_height'] + 10,
                font=FONTS['body'],
                fg_color=BUTTON_COLORS['inactive'],
                hover_color=BUTTON_COLORS['hover_inactive'],
                command=lambda s=stage: self._on_stage_clicked(s)
            )
            
            self.stage_buttons[stage] = button
    
    def _setup_layout(self):
        """Setup the layout of tab buttons."""
        self.tabs_frame.pack(fill="x", padx=10, pady=5)
        
        # Pack title
        self.title_label.pack(side="left", padx=(0, 20))
        
        # Pack stage buttons
        for stage in self.stages:
            self.stage_buttons[stage].pack(side="left", padx=2)
    
    def _on_stage_clicked(self, stage: str):
        """Handle stage button click."""
        if stage != self.active_stage:
            self.set_active_stage(stage)
            
            if self.stage_changed_callback:
                self.stage_changed_callback(stage)
            
            logger.debug(f"Stage changed to: {stage}")
    
    def set_active_stage(self, stage: str):
        """
        Set the active stage and update button appearance.
        
        Args:
            stage: Name of the stage to set as active
        """
        if stage not in self.stages:
            logger.warning(f"Unknown stage: {stage}")
            return
        
        # Update button colors
        for stage_name, button in self.stage_buttons.items():
            if stage_name == stage:
                # Active button
                button.configure(
                    fg_color=BUTTON_COLORS['active'],
                    hover_color=BUTTON_COLORS['hover_active'],
                    text_color="white"
                )
            else:
                # Inactive button
                button.configure(
                    fg_color=BUTTON_COLORS['inactive'],
                    hover_color=BUTTON_COLORS['hover_inactive'],
                    text_color=COLOR_SCHEME['text']
                )
        
        self.active_stage = stage
    
    def get_active_stage(self) -> Optional[str]:
        """
        Get the currently active stage.
        
        Returns:
            Name of the active stage or None if no stage is active
        """
        return self.active_stage
    
    def enable_stage(self, stage: str, enabled: bool = True):
        """
        Enable or disable a specific stage.
        
        Args:
            stage: Name of the stage to enable/disable
            enabled: True to enable, False to disable
        """
        if stage in self.stage_buttons:
            button = self.stage_buttons[stage]
            button.configure(state="normal" if enabled else "disabled")
    
    def set_stage_tooltip(self, stage: str, tooltip: str):
        """
        Set a tooltip for a specific stage button.
        
        Args:
            stage: Name of the stage
            tooltip: Tooltip text
        """
        # Note: CustomTkinter doesn't have built-in tooltips
        # This is a placeholder for potential future implementation
        pass
    
    def add_stage_indicator(self, stage: str, indicator_type: str = "dot"):
        """
        Add a visual indicator to a stage button.
        
        Args:
            stage: Name of the stage
            indicator_type: Type of indicator ("dot", "badge", etc.)
        """
        # This could be implemented to show completion status,
        # warnings, or other visual indicators
        pass
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the stage tabs and allow dynamic sizing of tab buttons."""
        try:
            # Update title label
            if hasattr(self, 'title_label'):
                self.title_label.configure(font=scaled_fonts['body'])
            # Update stage buttons
            for button in self.stage_buttons.values():
                button.configure(font=scaled_fonts['body'], width=140, height=36)
        except Exception as e:
            logger.warning(f"Error updating fonts in StageTabsContainer: {e}")

class StageProgressIndicator(ctk.CTkFrame):
    """
    Visual progress indicator for stages.
    
    This widget shows the overall progress through the analysis stages
    with visual indicators for completed, current, and pending stages.
    """
    
    def __init__(self, parent, stages: Optional[List[str]] = None):
        """
        Initialize the stage progress indicator.
        
        Args:
            parent: Parent widget
            stages: List of stage names
        """
        super().__init__(parent, fg_color="transparent")
        
        self.stages = stages or ["Visualize", "Preparation", "Builder", "Collective Vars"]
        self.current_stage_index = 0
        self.completed_stages = set()
        
        self._create_widgets()
        self._setup_layout()
    
    def _create_widgets(self):
        """Create progress indicator widgets."""
        self.progress_frame = ctk.CTkFrame(self, fg_color="transparent")
        
        self.stage_indicators = {}
        self.stage_connectors = {}
        
        for i, stage in enumerate(self.stages):
            # Stage indicator circle
            indicator = ctk.CTkButton(
                self.progress_frame,
                text=str(i + 1),
                width=30,
                height=30,
                corner_radius=15,
                font=("Arial", 12, "bold"),
                fg_color=COLOR_SCHEME['inactive'],
                state="disabled"
            )
            self.stage_indicators[stage] = indicator
            
            # Stage label
            label = ctk.CTkLabel(
                self.progress_frame,
                text=stage,
                font=FONTS['small']
            )
            self.stage_indicators[f"{stage}_label"] = label
            
            # Connector line (except for last stage)
            if i < len(self.stages) - 1:
                connector = ctk.CTkFrame(
                    self.progress_frame,
                    width=40,
                    height=2,
                    fg_color=COLOR_SCHEME['inactive']
                )
                self.stage_connectors[f"{stage}_connector"] = connector
    
    def _setup_layout(self):
        """Setup the layout of progress indicators."""
        self.progress_frame.pack(fill="x", padx=10, pady=5)
        
        # Layout indicators and connectors horizontally
        for i, stage in enumerate(self.stages):
            # Stage indicator
            indicator = self.stage_indicators[stage]
            indicator.grid(row=0, column=i*2, padx=5, pady=2)
            
            # Stage label
            label = self.stage_indicators[f"{stage}_label"]
            label.grid(row=1, column=i*2, padx=5, pady=2)
            
            # Connector (except for last stage)
            if i < len(self.stages) - 1:
                connector = self.stage_connectors[f"{stage}_connector"]
                connector.grid(row=0, column=i*2+1, padx=2, pady=2, sticky="ew")
    
    def set_current_stage(self, stage: str):
        """
        Set the current stage and update visual indicators.
        
        Args:
            stage: Name of the current stage
        """
        if stage not in self.stages:
            return
        
        self.current_stage_index = self.stages.index(stage)
        self._update_indicators()
    
    def mark_stage_completed(self, stage: str):
        """
        Mark a stage as completed.
        
        Args:
            stage: Name of the completed stage
        """
        self.completed_stages.add(stage)
        self._update_indicators()
    
    def _update_indicators(self):
        """Update the visual appearance of stage indicators."""
        for i, stage in enumerate(self.stages):
            indicator = self.stage_indicators[stage]
            
            if stage in self.completed_stages:
                # Completed stage - green
                indicator.configure(
                    fg_color=COLOR_SCHEME['active'],
                    text="✓"
                )
            elif i == self.current_stage_index:
                # Current stage - blue
                indicator.configure(
                    fg_color=COLOR_SCHEME['highlight'],
                    text=str(i + 1)
                )
            else:
                # Future stage - gray
                indicator.configure(
                    fg_color=COLOR_SCHEME['inactive'],
                    text=str(i + 1)
                )
    
    def reset_progress(self):
        """Reset all progress indicators."""
        self.current_stage_index = 0
        self.completed_stages.clear()
        self._update_indicators()