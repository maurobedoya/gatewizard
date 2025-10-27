# gatewizard/gui/frames/analysis.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Analysis frame for molecular dynamics trajectory analysis using MDAnalysis.

This module provides the GUI for setting up and running trajectory analysis
using MDAnalysis for calculations like RMSD, RMSF, distances, angles, etc.
"""

import tkinter as tk
from tkinter import messagebox, filedialog, simpledialog
from typing import Optional, Callable, List, Dict, Any
from pathlib import Path
from datetime import datetime
import json
import subprocess
import threading
import os
import re
import numpy as np

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis import rms, align
    MDANALYSIS_AVAILABLE = True
except ImportError:
    MDANALYSIS_AVAILABLE = False

try:
    import matplotlib.pyplot as plt
    import matplotlib.figure
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.backends._backend_tk import NavigationToolbar2Tk
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from gatewizard.gui.constants import (
    COLOR_SCHEME, WIDGET_SIZES, LAYOUT
)
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class AnalysisFrame(ctk.CTkFrame):
    """
    Frame for trajectory analysis setup and execution using MDAnalysis.
    
    This frame provides an interface for setting up various molecular dynamics
    trajectory analyses using the MDAnalysis library.
    """
    
    def __init__(
        self,
        parent,
        get_current_pdb: Optional[Callable[[], str]] = None,
        status_callback: Optional[Callable[[str], None]] = None,
        initial_directory: str = "",
        **kwargs
    ):
        """
        Initialize the Analysis frame.
        
        Args:
            parent: Parent widget
            get_current_pdb: Function to get current PDB file path
            status_callback: Function to call for status updates
            initial_directory: Initial directory for file dialogs
            **kwargs: Additional arguments for CTkFrame
        """
        super().__init__(parent, **kwargs)
        
        self.get_current_pdb = get_current_pdb
        self.status_callback = status_callback or self._default_status_callback
        self.initial_directory = initial_directory or os.getcwd()
        
        # File paths - separate lists for structural and energetic tabs
        self.topology_file = tk.StringVar()
        self.structural_trajectory_files = []  # Files for structural analysis
        self.energetic_log_files = []  # Files for energetic analysis
        self.trajectory_files = []  # Active file list (points to one of the above)
        self.discovered_log_files = {}  # Store discovered log files with basename as key
        self.output_directory = tk.StringVar(value=os.getcwd())  # Auto-assign to current directory
        
        # Time tracking for files (simulation time in ns)
        self.structural_file_times = {}  # {filepath: time_in_ns}
        self.energetic_file_times = {}  # {filepath: time_in_ns}
        self.file_times = {}  # Active time dict (points to one of the above)
        
        # Analysis parameters
        self.analysis_type = tk.StringVar(value="-")
        self.atom_selection = tk.StringVar(value="protein")
        self.custom_atom_selection = tk.StringVar(value="")
        self.use_custom_selection = tk.BooleanVar(value=False)
        self.reference_frame = tk.StringVar(value="0")  # Changed to StringVar to avoid IntVar empty string errors
        self.align_rmsd = tk.BooleanVar(value=True)  # Align structures before RMSD calculation
        
        # Distance analysis parameters
        self.distance_selection1 = tk.StringVar(value="protein and resid 1")
        self.distance_selection2 = tk.StringVar(value="protein and resid 50")
        
        # Radius of gyration parameters
        self.rg_selection = tk.StringVar(value="protein")
        
        # Current analysis tab (structural or energetic)
        self.current_analysis_tab = tk.StringVar(value="Structural")
        
        # NAMD log analysis parameters
        self.namd_log_columns = tk.StringVar(value="")  # Selected columns for NAMD log analysis
        self.namd_available_columns = []  # Store available NAMD columns
        self.column_checkboxes = {}  # Dictionary to store checkbox variables for multi-column selection
        
        # Plot parameters
        self.plot_enabled = tk.BooleanVar(value=True)
        self.plot_title = tk.StringVar(value="RMSD Analysis")
        self.plot_xlabel = tk.StringVar(value="Time (ns)")  # Updated to match default time unit
        self.plot_ylabel = tk.StringVar(value="RMSD (Å)")   # Updated to match default distance unit
        self.plot_color = tk.StringVar(value="blue")
        self.plot_background_color = tk.StringVar(value="#2b2b2b")  # Dark gray background
        self.plot_line_color = tk.StringVar(value="blue")  # Line color (separate from overall color scheme)
        self.plot_xlim_min = tk.StringVar(value="")
        self.plot_xlim_max = tk.StringVar(value="")
        self.plot_ylim_min = tk.StringVar(value="")
        self.plot_ylim_max = tk.StringVar(value="")
        self.plot_units = tk.StringVar(value="Å")  # Distance unit selection (Y-axis) - Default to Angstroms
        self.plot_time_units = tk.StringVar(value="ns")  # Time unit selection (X-axis) - Default to nanoseconds
        self.plot_energy_units = tk.StringVar(value="kcal/mol")  # Energy unit selection (Y-axis for energetic analysis)
        
        # RMSF-specific plot options
        self.rmsf_xaxis_type = tk.StringVar(value="residue_number")  # Options: "residue_number", "residue_type_number", "atom_index"
        self.rmsf_show_residue_labels = tk.BooleanVar(value=True)  # Show residue labels on X-axis
        self.rmsf_residue_name_format = tk.StringVar(value="single")  # Options: "single" (1-letter), "triple" (3-letter)
        self.rmsf_label_frequency = tk.StringVar(value="auto")  # Options: "all", "auto", "every_2", "every_5", "every_10", "every_20"
        
        # Figure background color (the frame around the plot)
        self.plot_figure_bg_color = tk.StringVar(value="#212121")  # Dark background for figure frame
        
        # Text and axes color control
        self.plot_text_color = tk.StringVar(value="Auto")  # Options: "Auto" (based on luminance), preset colors, "Custom RGB"
        
        # Grid visibility control
        self.plot_show_grid = tk.BooleanVar(value=True)  # Show grid lines by default
        
        # Dynamic Y-axis units based on property type
        self.plot_yaxis_units = tk.StringVar(value="kcal/mol")  # General Y-axis units (auto-updates based on data type)
        self.current_property_type = "energy"  # Track current property type: energy, temperature, pressure, volume
        
        # Define unit systems for different property types
        self.unit_systems = {
            'energy': {
                'units': ["kcal/mol", "kJ/mol"],
                'default': "kcal/mol",
                'conversions': {
                    'kcal/mol': 1.0,
                    'kJ/mol': 4.184
                }
            },
            'temperature': {
                'units': ["K", "°C", "°F"],
                'default': "K",
                'conversions': {}  # Special handling needed
            },
            'pressure': {
                'units': ["atm", "bar", "Pa", "kPa", "MPa"],
                'default': "atm",
                'conversions': {
                    'atm': 1.0,
                    'bar': 1.01325,
                    'Pa': 101325.0,
                    'kPa': 101.325,
                    'MPa': 0.101325
                }
            },
            'volume': {
                'units': ["Å³", "nm³", "mL", "L"],
                'default': "Å³",
                'conversions': {
                    'Å³': 1.0,
                    'nm³': 0.001,
                    'mL': 1.66054e-24,  # Å³ to mL
                    'L': 1.66054e-27    # Å³ to L
                }
            }
        }
        
        # Analysis state
        self.current_trajectory = None
        self.analysis_results: Dict[str, Any] = {}
        self.running_analysis = False
        self.plot_frame = None
        self.canvas = None
        self.toolbar = None
        
        # Separate state storage for each tab
        self.structural_state = {
            'analysis_results': {},
            'figure': None,
            'analysis_type': '-',
            'plot_color': 'Blue (#1f77b4)',
            'plot_background_color': 'Dark (#2b2b2b)',
            'plot_title': 'RMSD Analysis',
            'plot_xlabel': 'Time (ps)',
            'plot_ylabel': 'RMSD (nm)'
        }
        self.energetic_state = {
            'analysis_results': {},
            'figure': None,
            'analysis_type': '-',
            'available_columns': [],
            'plot_color': 'Blue (#1f77b4)',
            'plot_background_color': 'Dark (#2b2b2b)',
            'plot_title': 'Energy Analysis',
            'plot_xlabel': 'Time (ps)',
            'plot_ylabel': 'Energy (kcal/mol)'
        }
        
        # UI setup
        self.setup_ui()
        
        # Add trace callbacks to monitor parameter changes
        self._setup_parameter_tracking()
        
        # Check MDAnalysis availability
        if MDANALYSIS_AVAILABLE:
            self.status_callback("Analysis frame initialized with MDAnalysis support")
        else:
            self.status_callback("Analysis frame initialized - MDAnalysis not available")
    
    def _default_status_callback(self, message: str):
        """Default status callback."""
        logger.info(f"Analysis: {message}")
    
    def _setup_parameter_tracking(self):
        """Set up trace callbacks to monitor parameter changes."""
        # Track atom selection changes
        self.atom_selection.trace_add('write', lambda *args: self._on_parameter_change())
        self.custom_atom_selection.trace_add('write', lambda *args: self._on_parameter_change())
        
        # Track reference frame changes with validation
        self.reference_frame.trace_add('write', lambda *args: self._on_reference_frame_change())
        
        # Track distance analysis selections
        self.distance_selection1.trace_add('write', lambda *args: self._on_parameter_change())
        self.distance_selection2.trace_add('write', lambda *args: self._on_parameter_change())
        
        # Track radius of gyration selection
        self.rg_selection.trace_add('write', lambda *args: self._on_parameter_change())
        
        # Track output directory changes
        self.output_directory.trace_add('write', lambda *args: self._on_parameter_change())
    
    def _on_reference_frame_change(self):
        """Called when reference frame changes. Only marks as changed, doesn't auto-fill."""
        try:
            # Get the current value as string
            value = self.reference_frame.get().strip()
            
            # If empty, just return (don't auto-fill during typing)
            if not value:
                return
            
            # Try to parse as integer to validate
            int(value)
            
            # If valid, mark as changed
            self._on_parameter_change()
        except ValueError:
            # If not a valid integer, just ignore (don't auto-correct during typing)
            pass
    
    def _validate_reference_frame(self, event=None):
        """Validate reference frame when focus is lost. Sets to 0 if empty or invalid."""
        try:
            value = self.reference_frame.get().strip()
            
            # If empty, set to "0"
            if not value:
                self.reference_frame.set("0")
                return
            
            # Try to parse as integer
            int_value = int(value)
            
            # Ensure it's non-negative
            if int_value < 0:
                self.reference_frame.set("0")
        except ValueError:
            # If not a valid integer, set to "0"
            self.reference_frame.set("0")
    
    def _on_parameter_change(self):
        """Called when any analysis parameter changes."""
        # Only mark as changed if we're in the Structural tab (has status button)
        if self.current_analysis_tab.get() == "Structural":
            self._mark_analysis_changed()
    
    def setup_ui(self):
        """Set up the user interface."""
        # Configure grid weights for better layout
        self.grid_columnconfigure(0, weight=1)  # Left panel for settings (smaller)
        self.grid_columnconfigure(1, weight=2)  # Right panel for plot (larger)
        self.grid_rowconfigure(1, weight=1)    # Make main content expandable
        
        ## Title
        #title_label = ctk.CTkLabel(
        #    self,
        #    text="Trajectory Analysis with MDAnalysis",
        #    font=ctk.CTkFont(size=20, weight="bold"),
        #    text_color=COLOR_SCHEME["highlight"]
        #)
        #title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20), sticky="ew")
        
        # MDAnalysis status
        if not MDANALYSIS_AVAILABLE:
            warning_label = ctk.CTkLabel(
                self,
                text="⚠️ MDAnalysis not available. Please install with: conda install -c conda-forge mdanalysis",
                font=ctk.CTkFont(size=12),  # Increased from 10 to 12
                text_color="orange"
            )
            warning_label.grid(row=1, column=0, columnspan=2, pady=(0, 10), sticky="ew")
            current_row = 2
        else:
            current_row = 1
        
        # Create main content area
        self.create_main_content(current_row)
    
    def create_main_content(self, row):
        """Create main content with scrollable left panel and plot on right."""
        # Configure the main content row to expand
        self.grid_rowconfigure(row, weight=1)
        
        # Left panel - scrollable settings area with fixed width
        self.left_panel = ctk.CTkScrollableFrame(
            self, 
            width=380,  # Fixed width for settings
            height=500,  # Minimum height to ensure content visibility
            fg_color=COLOR_SCHEME["content_bg"]
        )
        self.left_panel.grid(row=row, column=0, sticky="nsew", padx=(10, 5), pady=10)
        self.left_panel.grid_columnconfigure(0, weight=1)
        # Configure rows - give more space to content (row 1)
        self.left_panel.grid_rowconfigure(0, weight=0)  # Tab buttons - fixed height
        self.left_panel.grid_rowconfigure(1, weight=1)  # Tab content - expandable
        
        # Enable mouse wheel scrolling for the left panel
        self._enable_scrolling_for_panel(self.left_panel)
        
        # Create subtabs for Structural vs Energetic Analysis
        self.create_analysis_subtabs()
        
        # Right panel - plot display area that expands
        self.right_panel = ctk.CTkFrame(self, fg_color=COLOR_SCHEME["content_bg"])
        self.right_panel.grid(row=row, column=1, sticky="nsew", padx=(5, 10), pady=10)
        self.right_panel.grid_columnconfigure(0, weight=1)
        self.right_panel.grid_rowconfigure(1, weight=1)
        
        # Add plot area to right panel
        self.create_plot_display()
    
    def create_analysis_subtabs(self):
        """Create subtabs for Structural and Energetic Analysis."""
        # Tab buttons frame
        tab_frame = ctk.CTkFrame(self.left_panel, fg_color="transparent")
        tab_frame.grid(row=0, column=0, sticky="ew", padx=10, pady=(5, 10))
        tab_frame.grid_columnconfigure(0, weight=1)
        tab_frame.grid_columnconfigure(1, weight=1)
        
        # Structural Analysis tab button
        self.structural_tab_btn = ctk.CTkButton(
            tab_frame,
            text="Structural Analysis",
            command=lambda: self.switch_analysis_tab("Structural"),
            height=35,
            fg_color=COLOR_SCHEME["highlight"],
            hover_color=COLOR_SCHEME["hover"]
        )
        self.structural_tab_btn.grid(row=0, column=0, padx=(0, 5), sticky="ew")
        
        # Energetic Analysis tab button
        self.energetic_tab_btn = ctk.CTkButton(
            tab_frame,
            text="Energetic Analysis",
            command=lambda: self.switch_analysis_tab("Energetic"),
            height=35,
            fg_color=COLOR_SCHEME["buttons"],
            hover_color=COLOR_SCHEME["hover"]
        )
        self.energetic_tab_btn.grid(row=0, column=1, padx=(5, 0), sticky="ew")
        
        # Container for tab content
        self.tab_content_frame = ctk.CTkFrame(self.left_panel, fg_color="transparent")
        self.tab_content_frame.grid(row=1, column=0, sticky="nsew", padx=0, pady=0)
        self.tab_content_frame.grid_columnconfigure(0, weight=1)
        # Give more weight to file section row (row 0)
        self.tab_content_frame.grid_rowconfigure(0, weight=2)  # File section gets more space
        self.tab_content_frame.grid_rowconfigure(1, weight=1)  # Analysis settings
        self.tab_content_frame.grid_rowconfigure(2, weight=1)  # Plot settings
        self.tab_content_frame.grid_rowconfigure(3, weight=1)  # Results
        self.tab_content_frame.grid_rowconfigure(4, weight=1)  # Control buttons
        
        # Start with Structural Analysis tab
        self.switch_analysis_tab("Structural")
    
    def switch_analysis_tab(self, tab_name):
        """Switch between Structural and Energetic analysis tabs."""
        # Save current state before switching
        current_tab = self.current_analysis_tab.get()
        if current_tab == "Structural":
            self.structural_trajectory_files = self.trajectory_files.copy()
            self.structural_file_times = self.file_times.copy()  # Save file times
            # Save current state
            self.structural_state['analysis_results'] = self.analysis_results.copy()
            self.structural_state['analysis_type'] = self.analysis_type.get()
            # Save plot settings
            self.structural_state['plot_color'] = self.plot_color.get()
            self.structural_state['plot_background_color'] = self.plot_background_color.get()
            self.structural_state['plot_title'] = self.plot_title.get()
            self.structural_state['plot_xlabel'] = self.plot_xlabel.get()
            self.structural_state['plot_ylabel'] = self.plot_ylabel.get()
            # Save RMSD data if it exists
            if hasattr(self, 'rmsd_data') and self.rmsd_data is not None:
                self.structural_state['rmsd_data'] = self.rmsd_data.copy() if hasattr(self.rmsd_data, 'copy') else self.rmsd_data
                if hasattr(self, 'time_data'):
                    self.structural_state['time_data'] = self.time_data.copy() if hasattr(self.time_data, 'copy') else self.time_data
        elif current_tab == "Energetic":
            self.energetic_log_files = self.trajectory_files.copy()
            self.energetic_file_times = self.file_times.copy()  # Save file times
            # Save current state
            self.energetic_state['analysis_results'] = self.analysis_results.copy()
            self.energetic_state['analysis_type'] = self.analysis_type.get()
            # Save plot settings
            self.energetic_state['plot_color'] = self.plot_color.get()
            self.energetic_state['plot_background_color'] = self.plot_background_color.get()
            self.energetic_state['plot_title'] = self.plot_title.get()
            self.energetic_state['plot_xlabel'] = self.plot_xlabel.get()
            self.energetic_state['plot_ylabel'] = self.plot_ylabel.get()
            if hasattr(self, 'namd_available_columns'):
                self.energetic_state['available_columns'] = self.namd_available_columns.copy()
        
        # Update current tab
        self.current_analysis_tab.set(tab_name)
        
        # Update button colors
        if tab_name == "Structural":
            self.structural_tab_btn.configure(fg_color=COLOR_SCHEME["highlight"])
            self.energetic_tab_btn.configure(fg_color=COLOR_SCHEME["buttons"])
            # Restore structural files and state
            self.trajectory_files = self.structural_trajectory_files.copy()
            self.file_times = self.structural_file_times.copy()  # Restore file times
            self.analysis_results = self.structural_state['analysis_results'].copy()
            self.analysis_type.set(self.structural_state['analysis_type'])
            # Restore plot settings
            self.plot_color.set(self.structural_state['plot_color'])
            self.plot_background_color.set(self.structural_state['plot_background_color'])
            self.plot_title.set(self.structural_state['plot_title'])
            self.plot_xlabel.set(self.structural_state['plot_xlabel'])
            self.plot_ylabel.set(self.structural_state['plot_ylabel'])
            # Restore RMSD data if it exists
            if 'rmsd_data' in self.structural_state:
                self.rmsd_data = self.structural_state['rmsd_data']
                if 'time_data' in self.structural_state:
                    self.time_data = self.structural_state['time_data']
        else:
            self.structural_tab_btn.configure(fg_color=COLOR_SCHEME["buttons"])
            self.energetic_tab_btn.configure(fg_color=COLOR_SCHEME["highlight"])
            # Restore energetic files and state
            self.trajectory_files = self.energetic_log_files.copy()
            self.file_times = self.energetic_file_times.copy()  # Restore file times
            self.analysis_results = self.energetic_state['analysis_results'].copy()
            self.analysis_type.set(self.energetic_state['analysis_type'])
            # Restore plot settings
            self.plot_color.set(self.energetic_state['plot_color'])
            self.plot_background_color.set(self.energetic_state['plot_background_color'])
            self.plot_title.set(self.energetic_state['plot_title'])
            self.plot_xlabel.set(self.energetic_state['plot_xlabel'])
            self.plot_ylabel.set(self.energetic_state['plot_ylabel'])
            if 'available_columns' in self.energetic_state:
                self.namd_available_columns = self.energetic_state['available_columns'].copy()
        
        # Clear current content
        for widget in self.tab_content_frame.winfo_children():
            widget.destroy()
        
        # Create appropriate content
        if tab_name == "Structural":
            self.create_structural_analysis_content()
        else:
            self.create_energetic_analysis_content()
        
        # Re-enable scrolling after creating new content
        def safe_enable_scrolling():
            try:
                if self.winfo_exists() and self.left_panel.winfo_exists():
                    self._enable_scrolling_for_panel(self.left_panel)
            except:
                pass
        self.tab_content_frame.after(150, safe_enable_scrolling)
        
        # Trigger analysis type change to restore dynamic options
        if self.analysis_type.get() != "-":
            # Need to wait for widgets to be created
            def safe_analysis_type_change():
                try:
                    if self.winfo_exists():
                        self._on_analysis_type_change(self.analysis_type.get())
                except:
                    pass
            self.tab_content_frame.after(200, safe_analysis_type_change)
        
        # Restore the plot if it exists
        self._restore_plot_for_tab(tab_name)
    
    def _restore_plot_for_tab(self, tab_name):
        """Restore the saved plot for the specified tab."""
        if not self.canvas or not self.ax:
            return
        
        # Clear current plot
        self.ax.clear()
        
        if tab_name == "Structural":
            # Restore structural plot if data exists
            if hasattr(self, 'rmsd_data') and self.rmsd_data is not None and hasattr(self, 'time_data'):
                if self.plot_enabled.get():
                    self._plot_rmsd_data(self.rmsd_data, self.time_data)
        else:
            # Restore energetic plot if data exists
            if self.analysis_results and 'timestep' in self.analysis_results:
                if self.plot_enabled.get():
                    self._plot_namd_data()
        
        # If no data to plot, just refresh the canvas
        if not self.ax.has_data():
            if hasattr(self, 'canvas') and self.canvas is not None:
                self.canvas.draw()
    
    def create_structural_analysis_content(self):
        """Create content for Structural Analysis tab."""
        current_row = 0
        
        # File selection section
        self.create_file_section_in_panel(current_row, file_type="dcd")
        current_row += 1
        
        # Analysis settings section
        self.create_analysis_section_in_panel(current_row, analysis_types=["-", "RMSD", "RMSF", "Distances", "Radius of Gyration"])
        current_row += 1
        
        # Plot settings section
        self.create_plot_section_in_panel(current_row)
        current_row += 1
        
        # Control buttons section - moved before Results for better workflow
        self.create_control_section_in_panel(current_row)
        current_row += 1
        
        # Results section
        self.create_results_section_in_panel(current_row)
        
        # Initialize labels to match default units
        self._sync_labels_with_units()
    
    def create_energetic_analysis_content(self):
        """Create content for Energetic Analysis tab."""
        current_row = 0
        
        # File selection section (for log files)
        self.create_file_section_in_panel(current_row, file_type="log")
        current_row += 1
        
        # Analysis settings section
        self.create_analysis_section_in_panel(current_row, analysis_types=["-", "NAMD Log Analysis"])
        current_row += 1
        
        # Plot settings section
        self.create_plot_section_in_panel(current_row)
        current_row += 1
        
        # Control buttons section - moved before Results for better workflow
        self.create_control_section_in_panel(current_row)
        current_row += 1
        
        # Results section
        self.create_results_section_in_panel(current_row)
    
    def create_settings_sections(self):
        """Create all settings sections in the scrollable left panel."""
        current_row = 0
        
        # File selection section
        self.create_file_section_in_panel(current_row)
        current_row += 1
        
        # Analysis settings section
        self.create_analysis_section_in_panel(current_row)
        current_row += 1
        
        # Plot settings section
        self.create_plot_section_in_panel(current_row)
        current_row += 1
        
        # Control buttons section - moved before Results for better workflow
        self.create_control_section_in_panel(current_row)
        current_row += 1
        
        # Results section
        self.create_results_section_in_panel(current_row)
    
    def create_file_section_in_panel(self, row, file_type="dcd"):
        """Create file selection section in the left panel.
        
        Args:
            row: Row position in grid
            file_type: Type of files to work with ("dcd" for structural, "log" for energetic)
        """
        # Store file type for this instance
        self.current_file_type = file_type
        
        # File selection frame
        file_frame = ctk.CTkFrame(self.tab_content_frame, fg_color="transparent")
        file_frame.grid(row=row, column=0, sticky="ew", padx=10, pady=5)
        file_frame.grid_columnconfigure(1, weight=1)
        
        # Title (increased from 14 to 16) - Store as attribute for font updates
        self.file_title_label = ctk.CTkLabel(
            file_frame,
            text="Input Files",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        self.file_title_label.grid(row=0, column=0, columnspan=3, pady=(10, 5), sticky="w", padx=10)
        
        # Topology file (only for structural analysis)
        if file_type == "dcd":
            ctk.CTkLabel(file_frame, text="Topology:").grid(
                row=1, column=0, pady=5, padx=(10, 5), sticky="w"
            )
            
            topology_entry = ctk.CTkEntry(
                file_frame,
                textvariable=self.topology_file,
                width=200
            )
            topology_entry.grid(row=1, column=1, pady=5, padx=5, sticky="ew")
            
            topology_browse_btn = ctk.CTkButton(
                file_frame,
                text="Browse",
                command=self.select_topology_file,
                width=80,
                height=25
            )
            topology_browse_btn.grid(row=1, column=2, pady=5, padx=(5, 10))
            traj_row = 2
        else:
            traj_row = 1
        
        # Trajectories / Log files label with time info - spans full width
        label_text = "Trajectories (with simulation time):" if file_type == "dcd" else "Log Files (with simulation time):"
        ctk.CTkLabel(file_frame, text=label_text).grid(
            row=traj_row, column=0, columnspan=3, pady=5, padx=(10, 5), sticky="w"
        )
        
        # Create a scrollable frame for file list with time inputs - spans full width
        self.file_list_frame = ctk.CTkScrollableFrame(
            file_frame,
            height=350,  # Increased from 150 to 350 for better visibility
            fg_color=("gray90", "gray20")
        )
        self.file_list_frame.grid(row=traj_row+1, column=0, columnspan=3, pady=5, padx=10, sticky="ew")
        self.file_list_frame.grid_columnconfigure(0, weight=1)  # Filename column expands
        self.file_list_frame.grid_columnconfigure(1, weight=0)  # Time column fixed width
        
        # Populate file list
        self._refresh_file_list_display()
        
        # Trajectory buttons - below the file list, aligned to the right
        traj_buttons_frame = ctk.CTkFrame(file_frame, fg_color="transparent")
        traj_buttons_frame.grid(row=traj_row+2, column=0, columnspan=3, pady=5, padx=10, sticky="e")
        
        add_traj_btn = ctk.CTkButton(
            traj_buttons_frame,
            text="Add Files",
            command=lambda: self.add_trajectory_file(file_type),
            width=80,
            height=25
        )
        add_traj_btn.grid(row=0, column=0, pady=2, padx=2)
        
        # Add folder selection button
        folder_btn = ctk.CTkButton(
            traj_buttons_frame,
            text="Add Folder",
            command=lambda: self.add_trajectory_folder(file_type),
            width=80,
            height=25
        )
        folder_btn.grid(row=0, column=1, pady=2, padx=2)
        
        remove_traj_btn = ctk.CTkButton(
            traj_buttons_frame,
            text="Remove",
            command=self.remove_selected_file,
            width=80,
            height=25
        )
        remove_traj_btn.grid(row=0, column=2, pady=2, padx=2)
        
        # Auto-detect time button
        auto_time_btn = ctk.CTkButton(
            traj_buttons_frame,
            text="Auto-Detect Time",
            command=self._auto_detect_file_times,
            width=110,
            height=25
        )
        auto_time_btn.grid(row=0, column=3, pady=2, padx=2)
        
        # Output directory
        output_row = traj_row + 3  # Updated from traj_row + 1
        ctk.CTkLabel(file_frame, text="Output:").grid(
            row=output_row, column=0, pady=(5, 10), padx=(10, 5), sticky="w"
        )
        
        output_entry = ctk.CTkEntry(
            file_frame,
            textvariable=self.output_directory,
            width=200
        )
        output_entry.grid(row=output_row, column=1, pady=(5, 10), padx=5, sticky="ew")
        
        output_browse_btn = ctk.CTkButton(
            file_frame,
            text="Browse",
            command=self.select_output_directory,
            width=80,
            height=25
        )
        output_browse_btn.grid(row=output_row, column=2, pady=(5, 10), padx=(5, 10))
    
    def create_analysis_section_in_panel(self, row, analysis_types=None):
        """Create analysis settings section in the left panel.
        
        Args:
            row: Row position in grid
            analysis_types: List of analysis types to show in dropdown
        """
        if analysis_types is None:
            analysis_types = ["-", "RMSD", "RMSF", "Distances", "Radius of Gyration", "NAMD Log Analysis"]
        
        # Analysis settings frame
        self.settings_frame = ctk.CTkFrame(self.tab_content_frame, fg_color="transparent")
        self.settings_frame.grid(row=row, column=0, sticky="ew", padx=10, pady=5)
        self.settings_frame.grid_columnconfigure(1, weight=1)
        
        # Title (increased from 14 to 16) - Store as attribute for font updates
        self.settings_title_label = ctk.CTkLabel(
            self.settings_frame,
            text="Analysis Settings",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        self.settings_title_label.grid(row=0, column=0, columnspan=2, pady=(10, 5), sticky="w", padx=10)
        
        # Analysis type (always visible) - Store as attribute for font updates
        self.analysis_type_label = ctk.CTkLabel(self.settings_frame, text="Analysis Type:")
        self.analysis_type_label.grid(row=1, column=0, pady=5, padx=(10, 5), sticky="w")
        self.analysis_type_combo = ctk.CTkComboBox(
            self.settings_frame,
            values=analysis_types,
            variable=self.analysis_type,
            width=200,
            command=self._on_analysis_type_change
        )
        self.analysis_type_combo.grid(row=1, column=1, pady=5, padx=(5, 10), sticky="w")
        
        # Dynamic options container (initially empty)
        self.dynamic_options_frame = ctk.CTkFrame(self.settings_frame, fg_color="transparent")
        self.dynamic_options_frame.grid(row=2, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        self.dynamic_options_frame.grid_columnconfigure(1, weight=1)
        
        # Store references to dynamic widgets for cleanup
        self.dynamic_widgets = []
    
    def create_plot_section_in_panel(self, row):
        """Create plot settings section in the left panel."""
        # Plot settings frame
        plot_frame = ctk.CTkFrame(self.tab_content_frame, fg_color="transparent")
        plot_frame.grid(row=row, column=0, sticky="ew", padx=10, pady=5)
        plot_frame.grid_columnconfigure(1, weight=1)
        
        # Title (increased from 14 to 16) - Store as attribute for font updates
        self.plot_title_label = ctk.CTkLabel(
            plot_frame,
            text="Plot Settings",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        self.plot_title_label.grid(row=0, column=0, columnspan=2, pady=(10, 5), sticky="w", padx=10)
        
        # Enable plotting checkbox
        self.plot_checkbox = ctk.CTkCheckBox(
            plot_frame,
            text="Enable plotting",
            variable=self.plot_enabled,
            command=self._toggle_plot_settings
        )
        self.plot_checkbox.grid(row=1, column=0, columnspan=2, pady=5, padx=10, sticky="w")
        
        # Plot title - Store as attribute for font updates
        self.plot_title_label = ctk.CTkLabel(plot_frame, text="Title:")
        self.plot_title_label.grid(row=2, column=0, pady=5, padx=(10, 5), sticky="w")
        plot_title_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_title, width=200)
        plot_title_entry.grid(row=2, column=1, pady=5, padx=(5, 10), sticky="ew")
        
        # Plot line color with custom RGB support
        ctk.CTkLabel(plot_frame, text="Line Color:").grid(
            row=3, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        # Create a frame for color selection
        line_color_frame = ctk.CTkFrame(plot_frame, fg_color="transparent")
        line_color_frame.grid(row=3, column=1, pady=5, padx=(5, 10), sticky="w")
        
        plot_color_combo = ctk.CTkComboBox(
            line_color_frame,
            values=["blue", "red", "green", "orange", "purple", "brown", "pink", "gray", "cyan", "magenta", "black", "Custom RGB"],
            variable=self.plot_line_color,
            width=120,
            command=self._on_line_color_change
        )
        plot_color_combo.grid(row=0, column=0, padx=(0, 5))
        
        # Custom RGB entry for line color (hidden by default)
        self.custom_line_color_entry = ctk.CTkEntry(
            line_color_frame,
            width=80,
            placeholder_text="#RRGGBB"
        )
        self.custom_line_color_entry.grid(row=0, column=1, padx=0)
        self.custom_line_color_entry.grid_remove()  # Hidden initially
        self.custom_line_color_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Plot background color with custom RGB support
        ctk.CTkLabel(plot_frame, text="Plot Area Color:").grid(
            row=4, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        # Create a frame for background color selection
        bg_color_frame = ctk.CTkFrame(plot_frame, fg_color="transparent")
        bg_color_frame.grid(row=4, column=1, pady=5, padx=(5, 10), sticky="w")
        
        plot_bg_combo = ctk.CTkComboBox(
            bg_color_frame,
            values=["Dark Gray (#2b2b2b)", "Black (#000000)", "White (#ffffff)", "Light Gray (#f0f0f0)", "Transparent", "Custom RGB"],
            variable=self.plot_background_color,
            width=120,
            command=self._on_background_color_change
        )
        plot_bg_combo.grid(row=0, column=0, padx=(0, 5))
        
        # Custom RGB entry for background color (hidden by default)
        self.custom_bg_color_entry = ctk.CTkEntry(
            bg_color_frame,
            width=80,
            placeholder_text="#RRGGBB"
        )
        self.custom_bg_color_entry.grid(row=0, column=1, padx=0)
        self.custom_bg_color_entry.grid_remove()  # Hidden initially
        self.custom_bg_color_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Figure background color (the frame around the plot area)
        ctk.CTkLabel(plot_frame, text="Border Color:").grid(
            row=5, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        # Create a frame for figure background color selection
        fig_bg_color_frame = ctk.CTkFrame(plot_frame, fg_color="transparent")
        fig_bg_color_frame.grid(row=5, column=1, pady=5, padx=(5, 10), sticky="w")
        
        plot_fig_bg_combo = ctk.CTkComboBox(
            fig_bg_color_frame,
            values=["Dark Gray (#212121)", "Black (#000000)", "White (#ffffff)", "Light Gray (#f0f0f0)", "Transparent", "Custom RGB"],
            variable=self.plot_figure_bg_color,
            width=120,
            command=self._on_figure_bg_color_change
        )
        plot_fig_bg_combo.grid(row=0, column=0, padx=(0, 5))
        
        # Custom RGB entry for figure background color (hidden by default)
        self.custom_fig_bg_color_entry = ctk.CTkEntry(
            fig_bg_color_frame,
            width=80,
            placeholder_text="#RRGGBB"
        )
        self.custom_fig_bg_color_entry.grid(row=0, column=1, padx=0)
        self.custom_fig_bg_color_entry.grid_remove()  # Hidden initially
        self.custom_fig_bg_color_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Text and axes color control
        ctk.CTkLabel(plot_frame, text="Text & Axes Color:").grid(
            row=6, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        # Create a frame for text color selection
        text_color_frame = ctk.CTkFrame(plot_frame, fg_color="transparent")
        text_color_frame.grid(row=6, column=1, pady=5, padx=(5, 10), sticky="w")
        
        text_color_combo = ctk.CTkComboBox(
            text_color_frame,
            values=["Auto", "Black (#000000)", "White (#ffffff)", "Gray (#808080)", "Custom RGB"],
            variable=self.plot_text_color,
            width=120,
            command=self._on_text_color_change
        )
        text_color_combo.grid(row=0, column=0, padx=(0, 5))
        
        # Custom RGB entry for text color (hidden by default)
        self.custom_text_color_entry = ctk.CTkEntry(
            text_color_frame,
            width=80,
            placeholder_text="#RRGGBB"
        )
        self.custom_text_color_entry.grid(row=0, column=1, padx=0)
        self.custom_text_color_entry.grid_remove()  # Hidden initially
        self.custom_text_color_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Grid visibility control
        ctk.CTkCheckBox(
            plot_frame,
            text="Show Grid",
            variable=self.plot_show_grid,
            command=self._update_plot_from_settings
        ).grid(row=7, column=0, columnspan=2, pady=5, padx=10, sticky="w")
        
        # Units selection (only for Structural Analysis - for Y-axis distance units)
        current_tab = self.current_analysis_tab.get()
        if current_tab == "Structural":
            ctk.CTkLabel(plot_frame, text="Y-Axis Units:").grid(
                row=8, column=0, pady=5, padx=(10, 5), sticky="w"
            )
            plot_units_combo = ctk.CTkComboBox(
                plot_frame,
                values=["nm", "Å"],  # Removed redundant "angstrom"
                variable=self.plot_units,
                width=150,
                command=self._on_distance_unit_change
            )
            plot_units_combo.grid(row=8, column=1, pady=5, padx=(5, 10), sticky="w")
            
            # Time units (X-axis) for RMSD analysis only - will be shown/hidden based on analysis type
            self.time_units_label = ctk.CTkLabel(plot_frame, text="X-Axis Units:")
            self.time_units_label.grid(row=9, column=0, pady=5, padx=(10, 5), sticky="w")
            
            self.plot_time_units_combo = ctk.CTkComboBox(
                plot_frame,
                values=["ps", "ns", "µs"],
                variable=self.plot_time_units,
                width=150,
                command=self._on_time_unit_change
            )
            self.plot_time_units_combo.grid(row=9, column=1, pady=5, padx=(5, 10), sticky="w")
            
            # RMSF-specific options (only shown for RMSF analysis)
            # These will be shown/hidden dynamically based on analysis type
            self.rmsf_xaxis_label = ctk.CTkLabel(plot_frame, text="RMSF X-Axis:")
            self.rmsf_xaxis_label.grid(row=10, column=0, pady=5, padx=(10, 5), sticky="w")
            
            self.rmsf_xaxis_combo = ctk.CTkComboBox(
                plot_frame,
                values=["Residue Number", "Residue Type+Number", "Atom Index"],
                variable=tk.StringVar(value="Residue Number"),
                width=150,
                command=self._on_rmsf_xaxis_change
            )
            self.rmsf_xaxis_combo.grid(row=10, column=1, pady=5, padx=(5, 10), sticky="w")
            
            # Residue name format (1-letter vs 3-letter code)
            self.rmsf_format_label = ctk.CTkLabel(plot_frame, text="Residue Format:")
            self.rmsf_format_label.grid(row=11, column=0, pady=5, padx=(10, 5), sticky="w")
            
            self.rmsf_format_combo = ctk.CTkComboBox(
                plot_frame,
                values=["1-Letter Code", "3-Letter Code"],
                variable=tk.StringVar(value="1-Letter Code"),
                width=150,
                command=self._on_rmsf_format_change
            )
            self.rmsf_format_combo.grid(row=11, column=1, pady=5, padx=(5, 10), sticky="w")
            
            # Label frequency control
            self.rmsf_frequency_label = ctk.CTkLabel(plot_frame, text="Label Frequency:")
            self.rmsf_frequency_label.grid(row=12, column=0, pady=5, padx=(10, 5), sticky="w")
            
            self.rmsf_frequency_combo = ctk.CTkComboBox(
                plot_frame,
                values=["Auto", "All Labels", "Every 2nd", "Every 5th", "Every 10th", "Every 20th"],
                variable=tk.StringVar(value="Auto"),
                width=150,
                command=self._on_rmsf_frequency_change
            )
            self.rmsf_frequency_combo.grid(row=12, column=1, pady=5, padx=(5, 10), sticky="w")
            
            # Checkbox for showing residue labels
            self.rmsf_labels_checkbox = ctk.CTkCheckBox(
                plot_frame,
                text="Show residue labels",
                variable=self.rmsf_show_residue_labels,
                command=self._update_plot_from_settings
            )
            self.rmsf_labels_checkbox.grid(row=13, column=0, columnspan=2, pady=5, padx=10, sticky="w")
            
            # Store RMSF widgets for show/hide (including labels)
            self.rmsf_plot_widgets = [
                self.rmsf_xaxis_label,
                self.rmsf_xaxis_combo,
                self.rmsf_format_label,
                self.rmsf_format_combo,
                self.rmsf_frequency_label,
                self.rmsf_frequency_combo,
                self.rmsf_labels_checkbox
            ]
            
            # Store time units widgets for show/hide
            self.time_units_widgets = [
                self.time_units_label,
                self.plot_time_units_combo
            ]
            
            # Initially hide RMSF widgets and show time units (default for RMSD)
            for widget in self.rmsf_plot_widgets:
                widget.grid_remove()
            
            next_row = 14
        else:
            # For Energetic analysis, skip units rows but add time units
            ctk.CTkLabel(plot_frame, text="X-Axis Units:").grid(
                row=8, column=0, pady=5, padx=(10, 5), sticky="w"
            )
            plot_time_units_combo = ctk.CTkComboBox(
                plot_frame,
                values=["ps", "ns", "µs"],
                variable=self.plot_time_units,
                width=150,
                command=self._on_time_unit_change
            )
            plot_time_units_combo.grid(row=8, column=1, pady=5, padx=(5, 10), sticky="w")
            
            # Add dynamic Y-Axis units for Energetic analysis
            ctk.CTkLabel(plot_frame, text="Y-Axis Units:").grid(
                row=9, column=0, pady=5, padx=(10, 5), sticky="w"
            )
            self.plot_yaxis_units_combo = ctk.CTkComboBox(
                plot_frame,
                values=["kcal/mol", "kJ/mol"],  # Default to energy units
                variable=self.plot_yaxis_units,
                width=150,
                command=self._on_yaxis_unit_change
            )
            self.plot_yaxis_units_combo.grid(row=9, column=1, pady=5, padx=(5, 10), sticky="w")
            next_row = 10
        
        # X axis settings
        ctk.CTkLabel(plot_frame, text="X Label:").grid(
            row=next_row, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_xlabel_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_xlabel, width=200)
        plot_xlabel_entry.grid(row=next_row, column=1, pady=5, padx=(5, 10), sticky="ew")
        plot_xlabel_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Y axis settings
        ctk.CTkLabel(plot_frame, text="Y Label:").grid(
            row=next_row+1, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_ylabel_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_ylabel, width=200)
        plot_ylabel_entry.grid(row=next_row+1, column=1, pady=5, padx=(5, 10), sticky="ew")
        plot_ylabel_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # X range settings
        ctk.CTkLabel(plot_frame, text="X Min:").grid(
            row=next_row+2, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_xlim_min_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_xlim_min, width=100)
        plot_xlim_min_entry.grid(row=next_row+2, column=1, pady=5, padx=(5, 10), sticky="w")
        plot_xlim_min_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        ctk.CTkLabel(plot_frame, text="X Max:").grid(
            row=next_row+3, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_xlim_max_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_xlim_max, width=100)
        plot_xlim_max_entry.grid(row=next_row+3, column=1, pady=5, padx=(5, 10), sticky="w")
        plot_xlim_max_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Y range settings
        ctk.CTkLabel(plot_frame, text="Y Min:").grid(
            row=next_row+4, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_ylim_min_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_ylim_min, width=100)
        plot_ylim_min_entry.grid(row=next_row+4, column=1, pady=5, padx=(5, 10), sticky="w")
        plot_ylim_min_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        ctk.CTkLabel(plot_frame, text="Y Max:").grid(
            row=next_row+5, column=0, pady=(5, 10), padx=(10, 5), sticky="w"
        )
        plot_ylim_max_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_ylim_max, width=100)
        plot_ylim_max_entry.grid(row=next_row+5, column=1, pady=(5, 10), padx=(5, 10), sticky="w")
        plot_ylim_max_entry.bind('<KeyRelease>', self._safe_update_plot_wrapper)
        
        # Reset button to restore default settings
        reset_button = ctk.CTkButton(
            plot_frame,
            text="Reset to Defaults",
            command=self._reset_plot_settings,
            width=200,
            height=28
        )
        reset_button.grid(row=next_row+6, column=0, columnspan=2, pady=(10, 5), padx=10)
        
        # Store references for enabling/disabling
        self.plot_widgets = [
            plot_title_entry, plot_color_combo, plot_bg_combo, plot_xlabel_entry, plot_ylabel_entry,
            plot_xlim_min_entry, plot_xlim_max_entry, plot_ylim_min_entry, plot_ylim_max_entry, reset_button
        ]
        # Add plot_units_combo only if it was created (Structural tab only)
        if current_tab == "Structural":
            try:
                if 'plot_units_combo' in locals():
                    self.plot_widgets.append(plot_units_combo)  # type: ignore[possibly-unbound]
            except NameError:
                pass  # plot_units_combo not defined, skip
        
        # Initial state
        self._toggle_plot_settings()
    
    def create_results_section_in_panel(self, row):
        """Create results section in the left panel."""
        # Results frame
        results_frame = ctk.CTkFrame(self.tab_content_frame, fg_color="transparent")
        results_frame.grid(row=row, column=0, sticky="ew", padx=10, pady=5)
        results_frame.grid_columnconfigure(1, weight=1)
        
        # Title (increased from 14 to 16) - Store as attribute for font updates
        self.results_title_label = ctk.CTkLabel(
            results_frame,
            text="Results & Export",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        self.results_title_label.grid(row=0, column=0, columnspan=2, pady=(10, 5), sticky="w", padx=10)
        
        # Analysis results info
        self.results_info = ctk.CTkTextbox(
            results_frame,
            height=100,
            width=300,
            wrap="word"
        )
        self.results_info.grid(row=1, column=0, columnspan=2, pady=5, padx=10, sticky="ew")
        
        # Enable mouse wheel scrolling for results textbox
        def _on_results_scroll(event):
            if event.delta:  # Windows/Mac
                # Scroll up (negative delta) or down (positive delta)
                direction = -1 if event.delta > 0 else 1
                self.results_info._textbox.yview_scroll(direction, "units")
            elif event.num == 4:  # Linux scroll up
                self.results_info._textbox.yview_scroll(-1, "units")
            elif event.num == 5:  # Linux scroll down
                self.results_info._textbox.yview_scroll(1, "units")
        
        self.results_info.bind("<MouseWheel>", _on_results_scroll)
        self.results_info.bind("<Button-4>", _on_results_scroll)
        self.results_info.bind("<Button-5>", _on_results_scroll)
        
        # Export buttons
        export_frame = ctk.CTkFrame(results_frame, fg_color="transparent")
        export_frame.grid(row=2, column=0, columnspan=2, pady=5, padx=10, sticky="ew")
        
        self.export_csv_btn = ctk.CTkButton(
            export_frame,
            text="Export CSV",
            command=self.export_csv,
            width=80,
            height=28
        )
        self.export_csv_btn.grid(row=0, column=0, padx=(0, 5), pady=5)
        
        self.export_json_btn = ctk.CTkButton(
            export_frame,
            text="Export JSON", 
            command=self.export_json,
            width=80,
            height=28
        )
        self.export_json_btn.grid(row=0, column=1, padx=5, pady=5)
        
        self.export_numpy_btn = ctk.CTkButton(
            export_frame,
            text="Export NumPy",
            command=self.export_numpy,
            width=80,
            height=28
        )
        self.export_numpy_btn.grid(row=0, column=2, padx=(5, 0), pady=5)
    
    def create_control_section_in_panel(self, row):
        """Create control buttons section in the left panel."""
        # Control frame
        control_frame = ctk.CTkFrame(self.tab_content_frame, fg_color="transparent")
        control_frame.grid(row=row, column=0, sticky="ew", padx=10, pady=5)
        control_frame.grid_columnconfigure(0, weight=1)
        
        # Get current tab
        current_tab = self.current_analysis_tab.get()
        
        # Title (only for Structural Analysis) (increased from 14 to 16) - Store as attribute for font updates
        if current_tab == "Structural":
            self.control_title_label = ctk.CTkLabel(
                control_frame,
                text="Analysis Control",
                font=ctk.CTkFont(size=16, weight="bold")
            )
            self.control_title_label.grid(row=0, column=0, pady=(10, 5), sticky="w", padx=10)
        
        # Run analysis button (increased from 14 to 16)
        self.run_button = ctk.CTkButton(
            control_frame,
            text="Run Analysis",
            command=self.run_analysis,
            height=40,
            font=ctk.CTkFont(size=16, weight="bold")
        )
        button_row = 1 if current_tab == "Structural" else 0
        self.run_button.grid(row=button_row, column=0, pady=10, padx=10, sticky="ew")
        
        # Status button (only for Structural Analysis) - clickable for quick status check
        if current_tab == "Structural":
            self.analysis_changed = False  # Track if analysis parameters have changed
            self.status_button = ctk.CTkButton(
                control_frame,
                text="Ready",
                width=100,
                height=30,
                fg_color="green",
                hover_color="darkgreen",
                font=ctk.CTkFont(size=14),
                command=self._on_status_button_click
            )
            self.status_button.grid(row=2, column=0, pady=(0, 10), padx=10)
    
    def create_plot_display(self):
        """Create the plot display area on the right side."""
        # Plot title
        #plot_title_label = ctk.CTkLabel(
        #    self.right_panel,
        #    text="RMSD Analysis Plot",
        #    font=ctk.CTkFont(size=16, weight="bold")
        #)
        #plot_title_label.grid(row=0, column=0, pady=(10, 5), sticky="n")
        
        # Create matplotlib figure and canvas
        self.figure, self.ax = plt.subplots(figsize=(10, 6))  # type: ignore[possibly-unbound]
        self.figure.patch.set_facecolor('#212121')  # Dark background
        self.ax.set_facecolor('#2b2b2b')
        
        # Configure subplot spacing and margins
        # left, right, top, bottom: Position of subplot edges (0-1 as fraction of figure)
        # wspace: Width space between subplots (only matters for multiple columns)
        # hspace: Height space between subplots (only matters for multiple rows)
        self.figure.subplots_adjust(
            left=0.08,    # Left margin (8% of figure width)
            right=0.95,   # Right margin (95% of figure width - leaves 5% on right)
            top=0.93,     # Top margin (93% of figure height - leaves 7% on top)
            bottom=0.10,  # Bottom margin (10% of figure height)
            wspace=0.2,   # Width space between subplots (not used in single subplot)
            hspace=0.2    # Height space between subplots (not used in single subplot)
        )
        
        # Style the plot for dark theme
        self.ax.tick_params(colors='white')
        self.ax.spines['bottom'].set_color('white')
        self.ax.spines['top'].set_color('white')
        self.ax.spines['left'].set_color('white')
        self.ax.spines['right'].set_color('white')
        self.ax.xaxis.label.set_color('white')
        self.ax.yaxis.label.set_color('white')
        self.ax.title.set_color('white')
        
        # Create canvas
        self.canvas = FigureCanvasTkAgg(self.figure, self.right_panel)  # type: ignore[possibly-unbound]
        self.canvas.get_tk_widget().grid(row=1, column=0, sticky="nsew", padx=10, pady=10)
        
        # Navigation toolbar
        toolbar_frame = ctk.CTkFrame(self.right_panel)
        toolbar_frame.grid(row=2, column=0, sticky="ew", padx=10, pady=(0, 10))
        
        self.toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)  # type: ignore[possibly-unbound]
        self.toolbar.update()
        
        # Configure the right panel grid weights
        self.right_panel.grid_rowconfigure(1, weight=1)
    
    def _update_plot_from_settings(self, *args):
        """Update the plot using existing data and current settings based on active tab."""
        try:
            # Safety checks to prevent crashes
            if not hasattr(self, 'plot_enabled') or not self.plot_enabled.get():
                return
            
            # Check if we have a valid figure and canvas
            if not hasattr(self, 'figure') or self.figure is None:
                return
            
            if not hasattr(self, 'canvas') or self.canvas is None:
                return
            
            # Determine which tab is active and plot accordingly
            if not hasattr(self, 'current_analysis_tab'):
                return
                
            current_tab = self.current_analysis_tab.get()
            
            if current_tab == "Structural":
                # Check what type of structural analysis data we have
                if hasattr(self, 'analysis_results') and self.analysis_results:
                    if 'rmsf' in self.analysis_results and 'residue_info' in self.analysis_results:
                        # We have RMSF data
                        try:
                            self._plot_rmsf_data(self.analysis_results['rmsf'], 
                                                self.analysis_results['residue_info'])
                        except Exception as e:
                            logger.warning(f"Error updating RMSF plot: {e}")
                            return
                    elif hasattr(self, 'rmsd_data') and self.rmsd_data is not None and hasattr(self, 'time_data'):
                        # We have RMSD data
                        try:
                            self._plot_rmsd_data(self.rmsd_data, self.time_data)
                        except Exception as e:
                            logger.warning(f"Error updating RMSD plot: {e}")
                            return
            elif current_tab == "Energetic":
                # Check if we have NAMD log data to plot
                if hasattr(self, 'analysis_results') and self.analysis_results and 'timestep' in self.analysis_results:
                    try:
                        self._plot_namd_data()
                    except Exception as e:
                        logger.warning(f"Error updating NAMD plot: {e}")
                        return
        except Exception as e:
            logger.error(f"Error in _update_plot_from_settings: {e}", exc_info=True)
    
    def _safe_update_plot_wrapper(self, event=None):
        """Safe wrapper for plot updates that handles errors gracefully."""
        try:
            if hasattr(self, '_plot_update_in_progress') and self._plot_update_in_progress:
                return  # Prevent recursive calls
            
            self._plot_update_in_progress = True
            self._update_plot_from_settings()
        except Exception as e:
            logger.warning(f"Error in plot update wrapper: {e}")
        finally:
            self._plot_update_in_progress = False
    
    def _configure_subplot_spacing(self, analysis_type="structural"):
        """Configure subplot spacing based on analysis type.
        
        Args:
            analysis_type: Either "structural" or "energetic"
                - structural: Standard spacing for single trajectory plots (RMSD, RMSF)
                - energetic: More space on right for legends when plotting multiple properties
        """
        if not hasattr(self, 'figure') or self.figure is None:
            return
        
        if analysis_type == "energetic":
            # Energetic analysis (NAMD logs) - often has legends on the right
            # Need more space on the right side for legend when plotting multiple columns
            self.figure.subplots_adjust(
                left=0.08,    # Left margin (8% of figure width)
                right=0.79,   # Right margin (85% - leaves 15% for legend on right)
                top=0.93,     # Top margin (93% of figure height)
                bottom=0.10,  # Bottom margin (10% of figure height)
                wspace=0.2,   # Width space between subplots (not used in single subplot)
                hspace=0.2    # Height space between subplots (not used in single subplot)
            )
        else:
            # Structural analysis (RMSD, RMSF) - typically single line, no legend
            # Standard spacing with more plot area
            self.figure.subplots_adjust(
                left=0.08,    # Left margin (8% of figure width)
                right=0.95,   # Right margin (95% - only 5% on right)
                top=0.93,     # Top margin (93% of figure height)
                bottom=0.10,  # Bottom margin (10% of figure height)
                wspace=0.2,   # Width space between subplots (not used in single subplot)
                hspace=0.2    # Height space between subplots (not used in single subplot)
            )
    
    def _on_distance_unit_change(self, value):
        """Callback when distance unit is changed. Auto-updates Y-axis label."""
        # Update the Y label to reflect the new unit
        current_ylabel = self.plot_ylabel.get()
        
        # Replace the unit in the label
        if "nm" in current_ylabel.lower():
            new_ylabel = current_ylabel.replace("nm", value).replace("Nm", value).replace("NM", value)
        elif "å" in current_ylabel.lower() or "angstrom" in current_ylabel.lower():
            new_ylabel = current_ylabel.replace("Å", value).replace("å", value).replace("angstrom", value).replace("Angstrom", value)
        else:
            # If no unit found, append it
            new_ylabel = f"{current_ylabel} ({value})"
        
        self.plot_ylabel.set(new_ylabel)
        
        # Update the plot
        self._update_plot_from_settings()
    
    def _on_line_color_change(self, value):
        """Handle line color selection change.
        
        Shows custom RGB entry field when 'Custom RGB' is selected.
        """
        if value == "Custom RGB":
            # Show the custom color entry
            if hasattr(self, 'custom_line_color_entry'):
                self.custom_line_color_entry.grid()
        else:
            # Hide the custom color entry
            if hasattr(self, 'custom_line_color_entry'):
                self.custom_line_color_entry.grid_remove()
            # Update plot with preset color
            self._update_plot_from_settings()
    
    def _on_background_color_change(self, value):
        """Handle background color selection change.
        
        Shows custom RGB entry field when 'Custom RGB' is selected.
        """
        if value == "Custom RGB":
            # Show the custom color entry
            if hasattr(self, 'custom_bg_color_entry'):
                self.custom_bg_color_entry.grid()
        else:
            # Hide the custom color entry
            if hasattr(self, 'custom_bg_color_entry'):
                self.custom_bg_color_entry.grid_remove()
            # Update plot with preset color
            self._update_plot_from_settings()
    
    def _on_figure_bg_color_change(self, value):
        """Handle figure background color selection change.
        
        The figure background is the frame/border area around the plot.
        Shows custom RGB entry field when 'Custom RGB' is selected.
        """
        if value == "Custom RGB":
            # Show the custom color entry
            if hasattr(self, 'custom_fig_bg_color_entry'):
                self.custom_fig_bg_color_entry.grid()
        else:
            # Hide the custom color entry
            if hasattr(self, 'custom_fig_bg_color_entry'):
                self.custom_fig_bg_color_entry.grid_remove()
            # Update plot with preset color
            self._update_plot_from_settings()
    
    def _on_text_color_change(self, value):
        """Handle text and axes color selection change.
        
        Shows custom RGB entry field when 'Custom RGB' is selected.
        """
        if value == "Custom RGB":
            # Show the custom color entry
            if hasattr(self, 'custom_text_color_entry'):
                self.custom_text_color_entry.grid()
        else:
            # Hide the custom color entry
            if hasattr(self, 'custom_text_color_entry'):
                self.custom_text_color_entry.grid_remove()
            # Update plot with preset or auto color
            self._safe_update_plot_wrapper()
    
    def _on_time_unit_change(self, value):
        """Callback when time unit is changed. Auto-updates X-axis label."""
        # Update the X label to reflect the new unit
        current_xlabel = self.plot_xlabel.get()
        
        # Replace the unit in the label
        if "ps" in current_xlabel.lower():
            new_xlabel = current_xlabel.replace("ps", value).replace("Ps", value).replace("PS", value)
        elif "ns" in current_xlabel.lower():
            new_xlabel = current_xlabel.replace("ns", value).replace("Ns", value).replace("NS", value)
        elif "µs" in current_xlabel.lower() or "us" in current_xlabel.lower():
            new_xlabel = current_xlabel.replace("µs", value).replace("us", value).replace("Us", value).replace("US", value)
        else:
            # If no unit found, append it
            new_xlabel = f"{current_xlabel} ({value})"
        
        self.plot_xlabel.set(new_xlabel)
        
        # Update the plot
        self._safe_update_plot_wrapper()
    
    def _on_energy_unit_change(self, value):
        """Callback when energy unit is changed. Auto-updates Y-axis label and re-plots."""
        # Update the Y label to reflect the new unit
        current_ylabel = self.plot_ylabel.get()
        
        # Replace the unit in the label
        if "kcal/mol" in current_ylabel:
            new_ylabel = current_ylabel.replace("kcal/mol", value)
        elif "kj/mol" in current_ylabel.lower():
            new_ylabel = current_ylabel.replace("kJ/mol", value).replace("kj/mol", value).replace("KJ/mol", value)
        else:
            # If no unit found, append it
            new_ylabel = f"{current_ylabel} ({value})"
        
        self.plot_ylabel.set(new_ylabel)
        
        # Re-plot with new units if we're in energetic analysis and have data
        if self.current_analysis_tab.get() == "Energetic" and self.analysis_results:
            self._safe_update_plot_wrapper()
    
    def _on_yaxis_unit_change(self, value):
        """Callback when Y-axis unit is changed. Auto-updates Y-axis label and re-plots."""
        # Update the Y label to reflect the new unit
        current_ylabel = self.plot_ylabel.get()
        
        # Extract the property name without units
        if "(" in current_ylabel and ")" in current_ylabel:
            property_name = current_ylabel.split('(')[0].strip()
            new_ylabel = f"{property_name} ({value})"
        else:
            new_ylabel = f"{current_ylabel} ({value})"
        
        self.plot_ylabel.set(new_ylabel)
        
        # Re-plot with new units if we're in energetic analysis and have data
        if self.current_analysis_tab.get() == "Energetic" and self.analysis_results:
            self._safe_update_plot_wrapper()
    
    def _on_rmsf_xaxis_change(self, value):
        """Callback when RMSF X-axis type is changed.
        
        Maps display names to internal values, updates X-label, and triggers plot update.
        """
        # Map display names to internal values
        mapping = {
            "Residue Number": "residue_number",
            "Residue Type+Number": "residue_type_number",
            "Atom Index": "atom_index"
        }
        self.rmsf_xaxis_type.set(mapping.get(value, "residue_number"))
        
        # Update X-label based on selection
        xlabel_mapping = {
            "Residue Number": "Residue Number",
            "Residue Type+Number": "Residue",
            "Atom Index": "Atom Index"
        }
        self.plot_xlabel.set(xlabel_mapping.get(value, "Residue Number"))
        
        # Update the plot if we have RMSF data
        if 'rmsf' in self.analysis_results and 'residue_info' in self.analysis_results:
            self._safe_update_plot_wrapper()
    
    def _on_rmsf_format_change(self, value):
        """Callback when RMSF residue name format is changed.
        
        Allows switching between 1-letter code (e.g., 'A') and 3-letter code (e.g., 'ALA').
        This is especially important for terminal caps (ACE, NME) which should always use 3-letter
        code to avoid confusion with alanine (A) or asparagine (N).
        """
        # Map display names to internal values
        mapping = {
            "1-Letter Code": "single",
            "3-Letter Code": "triple"
        }
        self.rmsf_residue_name_format.set(mapping.get(value, "single"))
        
        # Update the plot if we have RMSF data
        if 'rmsf' in self.analysis_results and 'residue_info' in self.analysis_results:
            self._safe_update_plot_wrapper()
    
    def _on_rmsf_frequency_change(self, value):
        """Callback when RMSF label frequency is changed.
        
        Controls how many labels are shown on the X-axis:
        - Auto: Automatically determines spacing based on number of residues
        - All Labels: Shows every residue label
        - Every Nth: Shows labels at regular intervals
        """
        # Map display names to internal values
        mapping = {
            "Auto": "auto",
            "All Labels": "all",
            "Every 2nd": "every_2",
            "Every 5th": "every_5",
            "Every 10th": "every_10",
            "Every 20th": "every_20"
        }
        self.rmsf_label_frequency.set(mapping.get(value, "auto"))
        
        # Update the plot if we have RMSF data
        if 'rmsf' in self.analysis_results and 'residue_info' in self.analysis_results:
            self._safe_update_plot_wrapper()
    
    def _detect_property_type_from_column(self, column_display_name):
        """Detect the property type from the column display name."""
        column_lower = column_display_name.lower()
        
        if 'temperature' in column_lower or 'temp' in column_lower:
            return 'temperature'
        elif 'pressure' in column_lower:
            return 'pressure'
        elif 'volume' in column_lower:
            return 'volume'
        elif 'energy' in column_lower or 'elect' in column_lower or 'vdw' in column_lower or \
             'bond' in column_lower or 'angle' in column_lower or 'dihed' in column_lower or \
             'poteng' in column_lower or 'kineng' in column_lower or 'toteng' in column_lower:
            return 'energy'
        else:
            return 'energy'  # Default to energy
    
    def _update_yaxis_units_for_property(self, property_type):
        """Update the Y-axis units combobox based on the property type."""
        if not hasattr(self, 'plot_yaxis_units_combo'):
            return
        
        if property_type in self.unit_systems:
            # Update the units available in the combobox
            unit_system = self.unit_systems[property_type]
            self.plot_yaxis_units_combo.configure(values=unit_system['units'])
            
            # Set to default unit for this property type
            default_unit = unit_system['default']
            self.plot_yaxis_units.set(default_unit)
            self.current_property_type = property_type
            
            # Update the Y-axis label
            current_ylabel = self.plot_ylabel.get()
            if "(" in current_ylabel and ")" in current_ylabel:
                property_name = current_ylabel.split('(')[0].strip()
            else:
                property_name = current_ylabel
            
            self.plot_ylabel.set(f"{property_name} ({default_unit})")
    
    def _convert_units(self, data, data_key, property_type, target_unit):
        """Convert data from native units to target units based on property type.
        
        Args:
            data: List of values to convert
            data_key: The key identifying the data column
            property_type: Type of property (energy, temperature, pressure, volume)
            target_unit: Target unit to convert to
            
        Returns:
            List of converted values
        """
        # Energy conversion (kcal/mol is native NAMD unit)
        if property_type == 'energy':
            energy_column_keys = [
                'poteng', 'potential', 'kineng', 'kinetic', 'toteng', 'total', 'total3',
                'elect', 'vdw', 'bond', 'angle', 'dihed', 'impr', 'misc', 'boundary'
            ]
            if data_key.lower() in energy_column_keys:
                conversion = self.unit_systems['energy']['conversions'].get(target_unit, 1.0)
                return [val * conversion for val in data]
        
        # Temperature conversion (K is native NAMD unit)
        elif property_type == 'temperature':
            if target_unit == "K":
                return data  # No conversion needed
            elif target_unit == "°C":
                return [val - 273.15 for val in data]
            elif target_unit == "°F":
                return [(val - 273.15) * 9/5 + 32 for val in data]
        
        # Pressure conversion (atm is native NAMD unit)
        elif property_type == 'pressure':
            if 'pressure' in data_key.lower():
                conversion = self.unit_systems['pressure']['conversions'].get(target_unit, 1.0)
                return [val * conversion for val in data]
        
        # Volume conversion (Å³ is native NAMD unit)
        elif property_type == 'volume':
            if 'volume' in data_key.lower():
                conversion = self.unit_systems['volume']['conversions'].get(target_unit, 1.0)
                return [val * conversion for val in data]
        
        # No conversion needed or unknown type
        return data
    
    def _sync_labels_with_units(self):
        """Synchronize axis labels with current unit selections.
        Called when analysis content is created to ensure labels match default units."""
        # Get current units
        time_unit = self.plot_time_units.get()
        distance_unit = self.plot_units.get()
        
        # Update X-axis label (time)
        current_xlabel = self.plot_xlabel.get()
        if "(" in current_xlabel:
            # Extract label base (e.g., "Time" from "Time (ps)")
            label_base = current_xlabel.split("(")[0].strip()
            self.plot_xlabel.set(f"{label_base} ({time_unit})")
        
        # Update Y-axis label (distance) - only for structural analysis
        current_tab = self.current_analysis_tab.get()
        if current_tab == "Structural":
            current_ylabel = self.plot_ylabel.get()
            if "(" in current_ylabel:
                # Extract label base (e.g., "RMSD" from "RMSD (nm)")
                label_base = current_ylabel.split("(")[0].strip()
                self.plot_ylabel.set(f"{label_base} ({distance_unit})")
    
    def _reset_plot_settings(self):
        """Reset all plot settings to default values."""
        current_tab = self.current_analysis_tab.get()
        
        # Default values based on analysis type
        if current_tab == "Structural":
            self.plot_title.set("RMSD Analysis")
            self.plot_xlabel.set("Time (ns)")
            self.plot_ylabel.set("RMSD (Å)")
            self.plot_units.set("Å")  # Default to Angstroms to match initial defaults
        else:  # Energetic
            self.plot_title.set("Energy Analysis")
            self.plot_xlabel.set("Time (ps)")
            self.plot_ylabel.set("Energy (kcal/mol)")
        
        # Common defaults
        self.plot_time_units.set("ns")  # Default to nanoseconds to match initial defaults
        self.plot_color.set("Blue (#1f77b4)")
        self.plot_background_color.set("Dark (#2b2b2b)")
        self.plot_xlim_min.set("")
        self.plot_xlim_max.set("")
        self.plot_ylim_min.set("")
        self.plot_ylim_max.set("")
        
        # Update the plot with new settings
        self._update_plot_from_settings()
    
    def _clear_plot(self):
        """Clear the current plot display."""
        if hasattr(self, 'ax') and self.ax is not None:
            self.ax.clear()
            
            # Add placeholder text
            self.ax.text(0.5, 0.5, 'Select analysis type and run analysis\nto display plot',
                        ha='center', va='center', fontsize=12, color='gray',
                        transform=self.ax.transAxes)
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            
            # Refresh canvas if it exists
            if hasattr(self, 'canvas') and self.canvas is not None:
                self.canvas.draw()
    
    def _try_load_and_plot_existing_data(self):
        """Try to load existing RMSD data from files and plot."""
        try:
            # Look for existing RMSD files in common locations
            possible_files = [
                "rmsd_data.csv",
                "rmsd_analysis.csv", 
                "analysis_results.csv",
                os.path.join(os.getcwd(), "rmsd_data.csv"),
                os.path.join(os.getcwd(), "rmsd_analysis.csv")
            ]
            
            loaded_data = None
            for file_path in possible_files:
                if os.path.exists(file_path):
                    try:
                        # Try loading as CSV
                        data = np.loadtxt(file_path, delimiter=',', skiprows=1)
                        if data.shape[1] >= 2:  # Should have time and RMSD columns
                            loaded_data = {
                                'time': data[:, 0],
                                'rmsd': data[:, 1]
                            }
                            break
                    except:
                        continue
            
            if loaded_data:
                self.rmsd_data = loaded_data['rmsd'] 
                self.time_data = loaded_data['time']
                self._plot_rmsd_data(self.rmsd_data, self.time_data)
                
        except Exception as e:
            print(f"Could not load existing RMSD data: {e}")
    
    def _plot_rmsd_data(self, rmsd_data, time_data):
        """Plot RMSD data with current settings."""
        if not self.plot_enabled.get():
            return
            
        # Clear previous plot
        self.ax.clear()
        
        # Configure subplot spacing for structural analysis
        self._configure_subplot_spacing(analysis_type="structural")
        
        # Set figure background color (the frame around the plot)
        fig_bg_color = self.plot_figure_bg_color.get()
        if fig_bg_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_fig_bg_color_entry'):
                custom_color = self.custom_fig_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    fig_bg_color = custom_color
                else:
                    fig_bg_color = "#212121"  # Fallback
            else:
                fig_bg_color = "#212121"
        elif fig_bg_color == "Transparent":
            fig_bg_color = "none"  # Matplotlib's transparent value
        elif "(" in fig_bg_color and ")" in fig_bg_color:
            # Extract hex code from "Name (#code)" format
            fig_bg_color = fig_bg_color.split("(")[1].split(")")[0]
        
        # Apply figure background
        if hasattr(self, 'figure') and self.figure is not None:
            self.figure.patch.set_facecolor(fig_bg_color)
        
        # Get background color - support both preset, transparent, and custom RGB
        bg_color = self.plot_background_color.get()
        if bg_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_bg_color_entry'):
                custom_color = self.custom_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    bg_color = custom_color
                else:
                    bg_color = "#2b2b2b"  # Fallback to dark gray
            else:
                bg_color = "#2b2b2b"
        elif bg_color == "Transparent":
            bg_color = "none"  # Matplotlib's transparent value
        elif "(" in bg_color and ")" in bg_color:
            # Extract hex code from "Name (#code)" format
            bg_color = bg_color.split("(")[1].split(")")[0]
        
        # Apply styling with user-selected background
        self.ax.set_facecolor(bg_color)
        
        # Get text color - support Auto (based on luminance), preset colors, and custom RGB
        text_color_setting = self.plot_text_color.get()
        if text_color_setting == "Auto":
            # Set text and spine colors based on background luminance
            # For transparent backgrounds, use a default dark text color
            if bg_color == "none":
                text_color = 'black'  # Default for transparent
            else:
                # Calculate luminance to determine if background is light or dark
                try:
                    # Remove '#' and convert hex to RGB
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    # Calculate relative luminance
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    # Fallback for light/dark detection
                    if bg_color in ["#ffffff", "#f0f0f0"]:
                        text_color = 'black'
                    else:
                        text_color = 'white'
        elif text_color_setting == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_text_color_entry'):
                custom_color = self.custom_text_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    text_color = custom_color
                else:
                    text_color = 'black'  # Fallback
            else:
                text_color = 'black'
        elif "(" in text_color_setting and ")" in text_color_setting:
            # Extract hex code from "Name (#code)" format
            text_color = text_color_setting.split("(")[1].split(")")[0]
        else:
            text_color = text_color_setting  # Use as is
            
        self.ax.tick_params(colors=text_color)
        for spine in self.ax.spines.values():
            spine.set_color(text_color)
        
        # Enable/disable grid based on user setting
        if self.plot_show_grid.get():
            self.ax.grid(True, alpha=0.3, color=text_color)
        
        # Convert units if needed
        plot_rmsd = rmsd_data.copy()
        plot_time = time_data.copy().astype(float)  # Time is stored in nanoseconds
        
        # Convert distance units (Y-axis)
        # NOTE: rmsd_data is stored in Angstroms (MDAnalysis default)
        units = self.plot_units.get()
        if units in ["nm", "nanometer"]:
            plot_rmsd /= 10  # Convert Å to nm
        # If Å, keep as is (already in Å)
        
        # Convert time units (X-axis) - time_data is in nanoseconds
        time_units = self.plot_time_units.get()
        if time_units == "ps":
            plot_time = plot_time * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = plot_time / 1000.0  # Convert ns to µs
        # If ns, keep as is (time_data is already in ns)
            
        # Plot data with user-selected line color - support custom RGB
        line_color = self.plot_line_color.get()
        if line_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_line_color_entry'):
                custom_color = self.custom_line_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    line_color = custom_color
                else:
                    line_color = "blue"  # Fallback
            else:
                line_color = "blue"
        
        self.ax.plot(plot_time, plot_rmsd, color=line_color, linewidth=2)
        
        # Set labels with appropriate color
        xlabel = self.plot_xlabel.get() or "Time (ps)"
        ylabel = self.plot_ylabel.get() or f"RMSD ({units})"
        self.ax.set_xlabel(xlabel, color=text_color)
        self.ax.set_ylabel(ylabel, color=text_color)
        
        # Set title with appropriate color
        title = self.plot_title.get() or "RMSD Analysis"
        self.ax.set_title(title, color=text_color, fontweight='bold')
        
        # Set axis limits if specified, otherwise use data range
        # X-axis limits: Use specified values or default to data range
        try:
            xlim_min = self.plot_xlim_min.get().strip()
            xlim_max = self.plot_xlim_max.get().strip()
            
            # Get current data limits as fallback
            current_xlim = self.ax.get_xlim()
            
            # Set X limits: use specified value or current data limit
            x_min = float(xlim_min) if xlim_min else current_xlim[0]
            x_max = float(xlim_max) if xlim_max else current_xlim[1]
            self.ax.set_xlim(x_min, x_max)
        except (ValueError, AttributeError):
            pass
            
        # Y-axis limits: Use specified values or default to data range
        try:
            ylim_min = self.plot_ylim_min.get().strip()
            ylim_max = self.plot_ylim_max.get().strip()
            
            # Get current data limits as fallback
            current_ylim = self.ax.get_ylim()
            
            # Set Y limits: use specified value or current data limit
            y_min = float(ylim_min) if ylim_min else current_ylim[0]
            y_max = float(ylim_max) if ylim_max else current_ylim[1]
            self.ax.set_ylim(y_min, y_max)
        except (ValueError, AttributeError):
            pass
        
        # Refresh canvas
        if hasattr(self, 'canvas') and self.canvas is not None:
            self.canvas.draw()
    
    def _plot_rmsf_data(self, rmsf_data, residue_info):
        """Plot RMSF data with residue-based X-axis.
        
        This function creates a customizable RMSF plot with the following features:
        - Unit conversion (nm ↔ Å)
        - Multiple X-axis types: atom index, residue number, or residue type+number
        - Configurable residue name format: 1-letter (e.g., 'A') or 3-letter (e.g., 'ALA')
        - Special handling for terminal caps (ACE, NME) - always shown in 3-letter code
          to avoid confusion with standard amino acids (A=Alanine, N=Asparagine)
        - Flexible label frequency control: auto, all, or every Nth label
        - Automatic axis limits when not specified by user
        
        Args:
            rmsf_data: Array of RMSF values (one per atom/residue)
            residue_info: List of dicts with residue metadata (residue_number, residue_name, etc.)
        """
        if not self.plot_enabled.get():
            return
            
        # Clear previous plot
        self.ax.clear()
        
        # Configure subplot spacing for structural analysis
        self._configure_subplot_spacing(analysis_type="structural")
        
        # Set figure background color (the frame around the plot)
        fig_bg_color = self.plot_figure_bg_color.get()
        if fig_bg_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_fig_bg_color_entry'):
                custom_color = self.custom_fig_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    fig_bg_color = custom_color
                else:
                    fig_bg_color = "#212121"  # Fallback
            else:
                fig_bg_color = "#212121"
        elif fig_bg_color == "Transparent":
            fig_bg_color = "none"  # Matplotlib's transparent value
        elif "(" in fig_bg_color and ")" in fig_bg_color:
            # Extract hex code from "Name (#code)" format
            fig_bg_color = fig_bg_color.split("(")[1].split(")")[0]
        
        # Apply figure background
        if hasattr(self, 'figure') and self.figure is not None:
            self.figure.patch.set_facecolor(fig_bg_color)
        
        # Get background color - support preset, transparent, and custom RGB
        bg_color = self.plot_background_color.get()
        if bg_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_bg_color_entry'):
                custom_color = self.custom_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    bg_color = custom_color
                else:
                    bg_color = "#2b2b2b"  # Fallback to dark gray
            else:
                bg_color = "#2b2b2b"
        elif bg_color == "Transparent":
            # Use matplotlib's transparent value for the plot background
            bg_color = "none"
        elif "(" in bg_color and ")" in bg_color:
            # Extract hex code from "Name (#code)" format
            bg_color = bg_color.split("(")[1].split(")")[0]
        
        # Apply styling with user-selected background
        self.ax.set_facecolor(bg_color)
        
        # Get text color - support Auto (based on luminance), preset colors, and custom RGB
        text_color_setting = self.plot_text_color.get()
        if text_color_setting == "Auto":
            # Set text and spine colors based on background luminance
            # For transparent backgrounds, default to black text (suitable for most presentation backgrounds)
            if bg_color == "none":
                text_color = 'black'
            else:
                # Calculate luminance to determine if background is light or dark
                try:
                    # Remove '#' and convert hex to RGB
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    # Calculate relative luminance
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    # Fallback for light/dark detection
                    if bg_color in ["#ffffff", "#f0f0f0"]:
                        text_color = 'black'
                    else:
                        text_color = 'white'
        elif text_color_setting == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_text_color_entry'):
                custom_color = self.custom_text_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    text_color = custom_color
                else:
                    text_color = 'black'  # Fallback
            else:
                text_color = 'black'
        elif "(" in text_color_setting and ")" in text_color_setting:
            # Extract hex code from "Name (#code)" format
            text_color = text_color_setting.split("(")[1].split(")")[0]
        else:
            text_color = text_color_setting  # Use as is
            
        self.ax.tick_params(colors=text_color)
        for spine in self.ax.spines.values():
            spine.set_color(text_color)
        
        # Enable/disable grid based on user setting
        if self.plot_show_grid.get():
            self.ax.grid(True, alpha=0.3, color=text_color)
        
        # Convert RMSF units if needed
        # NOTE: rmsf_data is stored in Angstroms (MDAnalysis default)
        plot_rmsf = rmsf_data.copy()
        units = self.plot_units.get()
        if units in ["nm", "nanometer"]:
            plot_rmsf /= 10  # Convert Å to nm
        # If Å, keep as is (already in Å)
        
        # Prepare X-axis data based on user preference
        xaxis_type = self.rmsf_xaxis_type.get()
        
        if xaxis_type == "atom_index":
            # X-axis: Atom indices (0, 1, 2, ...)
            x_values = [res['atom_index'] for res in residue_info]
            x_labels = None
            xlabel = "Atom Index"
        elif xaxis_type == "residue_type_number":
            # X-axis: Residue type + number (e.g., "A1", "ALA1", "ACE1")
            x_values = list(range(len(residue_info)))
            
            # Get user's preferred residue name format
            name_format = self.rmsf_residue_name_format.get()  # "single" or "triple"
            
            x_labels = []
            for res in residue_info:
                res_name = res['residue_name']
                res_num = res['residue_number']
                
                # Special handling for terminal caps (ACE, NME)
                # These should ALWAYS use 3-letter code to avoid confusion:
                # - ACE (acetyl cap) vs A (alanine)
                # - NME (N-methyl cap) vs N (asparagine)
                is_cap = res_name in ['ACE', 'NME']
                
                if is_cap or name_format == "triple":
                    # Use full 3-letter code (ALA, MET, TYR, ACE, NME, etc.)
                    label = f"{res_name}{res_num}"
                else:
                    # Use 1-letter code for standard amino acids
                    label = f"{res_name[0]}{res_num}"
                
                x_labels.append(label)
            
            xlabel = "Residue"
        else:  # residue_number (default)
            # X-axis: Residue numbers from PDB (may not be sequential)
            x_values = [res['residue_number'] for res in residue_info]
            x_labels = None
            xlabel = "Residue Number"
        
        # Plot RMSF data with user-selected line color - support custom RGB
        line_color = self.plot_line_color.get()
        if line_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_line_color_entry'):
                custom_color = self.custom_line_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    line_color = custom_color
                else:
                    line_color = "blue"  # Fallback
            else:
                line_color = "blue"
        
        self.ax.plot(x_values, plot_rmsf, color=line_color, linewidth=2, marker='o', markersize=3)
        
        # Configure X-axis labels if residue type+number is selected
        if x_labels and self.rmsf_show_residue_labels.get():
            n_labels = len(x_labels)
            frequency_setting = self.rmsf_label_frequency.get()
            
            # Determine label spacing based on user's frequency setting
            if frequency_setting == "all":
                # Show all labels
                tick_positions = list(range(n_labels))
                tick_labels = x_labels
            elif frequency_setting == "auto":
                # Auto mode: intelligently choose spacing based on data size
                if n_labels > 50:
                    # Many residues: show ~20 labels
                    step = max(1, n_labels // 20)
                    tick_positions = list(range(0, n_labels, step))
                    tick_labels = [x_labels[i] for i in tick_positions]
                else:
                    # Few residues: show all labels
                    tick_positions = list(range(n_labels))
                    tick_labels = x_labels
            else:
                # User specified frequency (every_2, every_5, every_10, every_20)
                step = int(frequency_setting.split("_")[1])  # Extract number from "every_N"
                tick_positions = list(range(0, n_labels, step))
                tick_labels = [x_labels[i] for i in tick_positions]
            
            # Apply tick positions and labels
            self.ax.set_xticks(tick_positions)
            self.ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        
        # Set axis labels with user-specified or default values
        xlabel_final = self.plot_xlabel.get() or xlabel
        ylabel = self.plot_ylabel.get() or f"RMSF ({units})"
        self.ax.set_xlabel(xlabel_final, color=text_color)
        self.ax.set_ylabel(ylabel, color=text_color)
        
        # Set title with appropriate styling
        title = self.plot_title.get() or "RMSF Analysis"
        self.ax.set_title(title, color=text_color, fontweight='bold')
        
        # Set axis limits if specified, otherwise use data range as default
        # X-axis limits: Use specified values or default to data range
        try:
            xlim_min = self.plot_xlim_min.get().strip()
            xlim_max = self.plot_xlim_max.get().strip()
            
            # Get current data limits as fallback
            current_xlim = self.ax.get_xlim()
            
            # Set X limits: use specified value or current data limit
            x_min = float(xlim_min) if xlim_min else current_xlim[0]
            x_max = float(xlim_max) if xlim_max else current_xlim[1]
            self.ax.set_xlim(x_min, x_max)
        except (ValueError, AttributeError):
            pass
            
        # Y-axis limits: Use specified values or default to data range
        try:
            ylim_min = self.plot_ylim_min.get().strip()
            ylim_max = self.plot_ylim_max.get().strip()
            
            # Get current data limits as fallback
            current_ylim = self.ax.get_ylim()
            
            # Set Y limits: use specified value or current data limit
            y_min = float(ylim_min) if ylim_min else current_ylim[0]
            y_max = float(ylim_max) if ylim_max else current_ylim[1]
            self.ax.set_ylim(y_min, y_max)
        except (ValueError, AttributeError):
            pass
        
        # Adjust layout to prevent label cutoff
        try:
            self.figure.tight_layout()
        except:
            pass
        
        # Refresh canvas to display changes
        if hasattr(self, 'canvas') and self.canvas is not None:
            self.canvas.draw()
    
    def _plot_distance_data(self):
        """Plot distance data between two atom selections."""
        if not self.plot_enabled.get():
            return
        
        # Get data from analysis results
        if 'distance' not in self.analysis_results or 'time' not in self.analysis_results:
            return
        
        distance_data = self.analysis_results['distance']  # In nm
        time_data = self.analysis_results['time']  # In ns
        
        # Clear previous plot
        self.ax.clear()
        
        # Configure subplot spacing for structural analysis
        self._configure_subplot_spacing(analysis_type="structural")
        
        # Set figure background color (the frame around the plot)
        fig_bg_color = self.plot_figure_bg_color.get()
        if fig_bg_color == "Custom RGB":
            if hasattr(self, 'custom_fig_bg_color_entry'):
                custom_color = self.custom_fig_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    fig_bg_color = custom_color
                else:
                    fig_bg_color = "#212121"
            else:
                fig_bg_color = "#212121"
        elif fig_bg_color == "Transparent":
            fig_bg_color = "none"
        elif "(" in fig_bg_color and ")" in fig_bg_color:
            fig_bg_color = fig_bg_color.split("(")[1].split(")")[0]
        
        if hasattr(self, 'figure') and self.figure is not None:
            self.figure.patch.set_facecolor(fig_bg_color)
        
        # Get background color
        bg_color = self.plot_background_color.get()
        if bg_color == "Custom RGB":
            if hasattr(self, 'custom_bg_color_entry'):
                custom_color = self.custom_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    bg_color = custom_color
                else:
                    bg_color = "#2b2b2b"
            else:
                bg_color = "#2b2b2b"
        elif bg_color == "Transparent":
            bg_color = "none"
        elif "(" in bg_color and ")" in bg_color:
            bg_color = bg_color.split("(")[1].split(")")[0]
        
        self.ax.set_facecolor(bg_color)
        
        # Get text color
        text_color_setting = self.plot_text_color.get()
        if text_color_setting == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white' if bg_color not in ["#ffffff", "#f0f0f0"] else 'black'
        elif text_color_setting == "Custom RGB":
            if hasattr(self, 'custom_text_color_entry'):
                custom_color = self.custom_text_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    text_color = custom_color
                else:
                    text_color = 'black'
            else:
                text_color = 'black'
        elif "(" in text_color_setting and ")" in text_color_setting:
            text_color = text_color_setting.split("(")[1].split(")")[0]
        else:
            text_color = text_color_setting
        
        self.ax.tick_params(colors=text_color)
        for spine in self.ax.spines.values():
            spine.set_color(text_color)
        
        # Enable/disable grid
        if self.plot_show_grid.get():
            self.ax.grid(True, alpha=0.3, color=text_color)
        
        # Convert units if needed
        plot_distance = distance_data.copy()
        plot_time = time_data.copy().astype(float)
        
        # Convert distance units (Y-axis)
        # NOTE: distance_data is stored in Angstroms (MDAnalysis default)
        units = self.plot_units.get()
        if units in ["nm", "nanometer"]:
            plot_distance /= 10  # Convert Å to nm
        # If Å, keep as is (already in Å)
        
        # Convert time units (X-axis)
        time_units = self.plot_time_units.get()
        if time_units == "ps":
            plot_time = plot_time * 1000.0
        elif time_units == "µs":
            plot_time = plot_time / 1000.0
        
        # Plot data with user-selected line color
        line_color = self.plot_line_color.get()
        if line_color == "Custom RGB":
            if hasattr(self, 'custom_line_color_entry'):
                custom_color = self.custom_line_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    line_color = custom_color
                else:
                    line_color = "blue"
            else:
                line_color = "blue"
        
        self.ax.plot(plot_time, plot_distance, color=line_color, linewidth=2)
        
        # Set labels
        xlabel = self.plot_xlabel.get() or f"Time ({time_units})"
        ylabel = self.plot_ylabel.get() or f"Distance ({units})"
        self.ax.set_xlabel(xlabel, color=text_color)
        self.ax.set_ylabel(ylabel, color=text_color)
        
        # Set title
        title = self.plot_title.get() or "Distance Analysis"
        self.ax.set_title(title, color=text_color, fontweight='bold')
        
        # Set axis limits if specified
        try:
            xlim_min = self.plot_xlim_min.get().strip()
            xlim_max = self.plot_xlim_max.get().strip()
            current_xlim = self.ax.get_xlim()
            x_min = float(xlim_min) if xlim_min else current_xlim[0]
            x_max = float(xlim_max) if xlim_max else current_xlim[1]
            self.ax.set_xlim(x_min, x_max)
        except (ValueError, AttributeError):
            pass
        
        try:
            ylim_min = self.plot_ylim_min.get().strip()
            ylim_max = self.plot_ylim_max.get().strip()
            current_ylim = self.ax.get_ylim()
            y_min = float(ylim_min) if ylim_min else current_ylim[0]
            y_max = float(ylim_max) if ylim_max else current_ylim[1]
            self.ax.set_ylim(y_min, y_max)
        except (ValueError, AttributeError):
            pass
        
        # Refresh canvas
        if hasattr(self, 'canvas') and self.canvas is not None:
            self.canvas.draw()
    
    def _plot_rg_data(self):
        """Plot radius of gyration data."""
        if not self.plot_enabled.get():
            return
        
        # Get data from analysis results
        if 'rg' not in self.analysis_results or 'time' not in self.analysis_results:
            return
        
        rg_data = self.analysis_results['rg']  # In nm
        time_data = self.analysis_results['time']  # In ns
        
        # Clear previous plot
        self.ax.clear()
        
        # Configure subplot spacing for structural analysis
        self._configure_subplot_spacing(analysis_type="structural")
        
        # Set figure background color (the frame around the plot)
        fig_bg_color = self.plot_figure_bg_color.get()
        if fig_bg_color == "Custom RGB":
            if hasattr(self, 'custom_fig_bg_color_entry'):
                custom_color = self.custom_fig_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    fig_bg_color = custom_color
                else:
                    fig_bg_color = "#212121"
            else:
                fig_bg_color = "#212121"
        elif fig_bg_color == "Transparent":
            fig_bg_color = "none"
        elif "(" in fig_bg_color and ")" in fig_bg_color:
            fig_bg_color = fig_bg_color.split("(")[1].split(")")[0]
        
        if hasattr(self, 'figure') and self.figure is not None:
            self.figure.patch.set_facecolor(fig_bg_color)
        
        # Get background color
        bg_color = self.plot_background_color.get()
        if bg_color == "Custom RGB":
            if hasattr(self, 'custom_bg_color_entry'):
                custom_color = self.custom_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    bg_color = custom_color
                else:
                    bg_color = "#2b2b2b"
            else:
                bg_color = "#2b2b2b"
        elif bg_color == "Transparent":
            bg_color = "none"
        elif "(" in bg_color and ")" in bg_color:
            bg_color = bg_color.split("(")[1].split(")")[0]
        
        self.ax.set_facecolor(bg_color)
        
        # Get text color
        text_color_setting = self.plot_text_color.get()
        if text_color_setting == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white' if bg_color not in ["#ffffff", "#f0f0f0"] else 'black'
        elif text_color_setting == "Custom RGB":
            if hasattr(self, 'custom_text_color_entry'):
                custom_color = self.custom_text_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    text_color = custom_color
                else:
                    text_color = 'black'
            else:
                text_color = 'black'
        elif "(" in text_color_setting and ")" in text_color_setting:
            text_color = text_color_setting.split("(")[1].split(")")[0]
        else:
            text_color = text_color_setting
        
        self.ax.tick_params(colors=text_color)
        for spine in self.ax.spines.values():
            spine.set_color(text_color)
        
        # Enable/disable grid
        if self.plot_show_grid.get():
            self.ax.grid(True, alpha=0.3, color=text_color)
        
        # Convert units if needed
        plot_rg = rg_data.copy()
        plot_time = time_data.copy().astype(float)
        
        # Convert Rg units (Y-axis)
        # NOTE: rg_data is stored in Angstroms (MDAnalysis default)
        units = self.plot_units.get()
        if units in ["nm", "nanometer"]:
            plot_rg /= 10  # Convert Å to nm
        # If Å, keep as is (already in Å)
        
        # Convert time units (X-axis)
        time_units = self.plot_time_units.get()
        if time_units == "ps":
            plot_time = plot_time * 1000.0
        elif time_units == "µs":
            plot_time = plot_time / 1000.0
        
        # Plot data with user-selected line color
        line_color = self.plot_line_color.get()
        if line_color == "Custom RGB":
            if hasattr(self, 'custom_line_color_entry'):
                custom_color = self.custom_line_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    line_color = custom_color
                else:
                    line_color = "blue"
            else:
                line_color = "blue"
        
        self.ax.plot(plot_time, plot_rg, color=line_color, linewidth=2)
        
        # Set labels
        xlabel = self.plot_xlabel.get() or f"Time ({time_units})"
        ylabel = self.plot_ylabel.get() or f"Radius of Gyration ({units})"
        self.ax.set_xlabel(xlabel, color=text_color)
        self.ax.set_ylabel(ylabel, color=text_color)
        
        # Set title
        title = self.plot_title.get() or "Radius of Gyration Analysis"
        self.ax.set_title(title, color=text_color, fontweight='bold')
        
        # Set axis limits if specified
        try:
            xlim_min = self.plot_xlim_min.get().strip()
            xlim_max = self.plot_xlim_max.get().strip()
            current_xlim = self.ax.get_xlim()
            x_min = float(xlim_min) if xlim_min else current_xlim[0]
            x_max = float(xlim_max) if xlim_max else current_xlim[1]
            self.ax.set_xlim(x_min, x_max)
        except (ValueError, AttributeError):
            pass
        
        try:
            ylim_min = self.plot_ylim_min.get().strip()
            ylim_max = self.plot_ylim_max.get().strip()
            current_ylim = self.ax.get_ylim()
            y_min = float(ylim_min) if ylim_min else current_ylim[0]
            y_max = float(ylim_max) if ylim_max else current_ylim[1]
            self.ax.set_ylim(y_min, y_max)
        except (ValueError, AttributeError):
            pass
        
        # Refresh canvas
        if hasattr(self, 'canvas') and self.canvas is not None:
            self.canvas.draw()
    
    def create_file_section(self, row):
        """Create file selection section."""
        # File selection frame - spans both columns
        file_frame = ctk.CTkFrame(self, fg_color=COLOR_SCHEME["content_bg"])
        file_frame.grid(row=row, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        file_frame.grid_columnconfigure(1, weight=1)
        
        # Title (increased from 14 to 16) - Store as attribute for font updates
        self.namd_file_title_label = ctk.CTkLabel(
            file_frame,
            text="Input Files",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        self.namd_file_title_label.grid(row=0, column=0, columnspan=3, pady=(10, 5), sticky="w", padx=10)
        
        # Topology file
        ctk.CTkLabel(file_frame, text="Topology:").grid(
            row=1, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        topology_entry = ctk.CTkEntry(
            file_frame,
            textvariable=self.topology_file,
            placeholder_text="Select PDB/PSF topology file...",
            width=400
        )
        topology_entry.grid(row=1, column=1, pady=5, padx=5, sticky="ew")
        
        topology_btn = ctk.CTkButton(
            file_frame,
            text="Browse",
            command=self.select_topology_file,
            width=80
        )
        topology_btn.grid(row=1, column=2, pady=5, padx=(5, 10))
        
        # Trajectory files
        ctk.CTkLabel(file_frame, text="Trajectory:").grid(
            row=2, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        self.trajectory_listbox = tk.Listbox(
            file_frame,
            height=3,
            font=("Consolas", 9),
            bg=COLOR_SCHEME["content_bg"],
            fg=COLOR_SCHEME["text"],
            selectbackground=COLOR_SCHEME["highlight"]
        )
        self.trajectory_listbox.grid(row=2, column=1, pady=5, padx=5, sticky="ew")
        
        traj_btn_frame = ctk.CTkFrame(file_frame, fg_color="transparent")
        traj_btn_frame.grid(row=2, column=2, pady=5, padx=(5, 10), sticky="n")
        
        add_traj_btn = ctk.CTkButton(
            traj_btn_frame,
            text="Add",
            command=self.add_trajectory_file,
            width=80,
            height=25
        )
        add_traj_btn.grid(row=0, column=0, pady=2)
        
        remove_traj_btn = ctk.CTkButton(
            traj_btn_frame,
            text="Remove",
            command=self.remove_trajectory_file,
            width=80,
            height=25
        )
        remove_traj_btn.grid(row=1, column=0, pady=2)
        
        # Note: Move Up/Down buttons removed - drag-and-drop is now used for reordering
        
        # Output directory
        ctk.CTkLabel(file_frame, text="Output:").grid(
            row=3, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        output_entry = ctk.CTkEntry(
            file_frame,
            textvariable=self.output_directory,
            placeholder_text="Select output directory...",
            width=400
        )
        output_entry.grid(row=3, column=1, pady=(5, 10), padx=5, sticky="ew")
        
        output_btn = ctk.CTkButton(
            file_frame,
            text="Browse",
            command=self.select_output_directory,
            width=80
        )
        output_btn.grid(row=3, column=2, pady=(5, 10), padx=(5, 10))
    
    def create_analysis_section(self, row):
        """Create analysis settings section."""
        # Analysis settings frame - left column only
        settings_frame = ctk.CTkFrame(self, fg_color=COLOR_SCHEME["content_bg"])
        settings_frame.grid(row=row, column=0, sticky="ew", padx=(10, 5), pady=5)
        settings_frame.grid_columnconfigure(1, weight=1)
        
        # Title
        settings_title = ctk.CTkLabel(
            settings_frame,
            text="Analysis Settings",
            font=ctk.CTkFont(size=14, weight="bold")
        )
        settings_title.grid(row=0, column=0, columnspan=2, pady=(10, 5), sticky="w", padx=10)
        
        # Analysis type
        ctk.CTkLabel(settings_frame, text="Analysis Type:").grid(
            row=1, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        analysis_combo = ctk.CTkComboBox(
            settings_frame,
            values=["rmsd", "rmsf", "distances", "angles", "radius_of_gyration"],
            variable=self.analysis_type,
            state="readonly",
            width=200
        )
        analysis_combo.grid(row=1, column=1, pady=5, padx=(5, 10), sticky="w")
        
        # Atom selection
        ctk.CTkLabel(settings_frame, text="Atom Selection:").grid(
            row=2, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        
        selection_combo = ctk.CTkComboBox(
            settings_frame,
            values=["protein", "backbone", "name CA", "protein and not name H", "all"],
            variable=self.atom_selection,
            state="readonly",
            width=200
        )
        selection_combo.grid(row=2, column=1, pady=5, padx=(5, 10), sticky="w")
        
        # Reference frame
        ctk.CTkLabel(settings_frame, text="Reference Frame:").grid(
            row=3, column=0, pady=(5, 10), padx=(10, 5), sticky="w"
        )
        
        ref_frame_entry = ctk.CTkEntry(
            settings_frame,
            textvariable=self.reference_frame,
            width=100
        )
        ref_frame_entry.grid(row=3, column=1, pady=(5, 10), padx=(5, 10), sticky="w")
        
        # Bind validation on focus out and Enter key
        ref_frame_entry.bind('<FocusOut>', self._validate_reference_frame)
        ref_frame_entry.bind('<Return>', self._validate_reference_frame)

    def create_plot_section(self, row):
        """Create plot settings section."""
        # Plot settings frame - right column only
        plot_frame = ctk.CTkFrame(self, fg_color=COLOR_SCHEME["content_bg"])
        plot_frame.grid(row=row, column=1, sticky="ew", padx=(5, 10), pady=5)
        plot_frame.grid_columnconfigure(1, weight=1)
        plot_frame.grid_columnconfigure(3, weight=1)
        
        # Title
        plot_title = ctk.CTkLabel(
            plot_frame,
            text="Plot Settings",
            font=ctk.CTkFont(size=14, weight="bold")
        )
        plot_title.grid(row=0, column=0, columnspan=4, pady=(10, 5), sticky="w", padx=10)
        
        # Enable plotting checkbox
        self.plot_checkbox = ctk.CTkCheckBox(
            plot_frame,
            text="Enable plotting",
            variable=self.plot_enabled,
            command=self._toggle_plot_settings
        )
        self.plot_checkbox.grid(row=1, column=0, columnspan=4, pady=5, padx=10, sticky="w")
        
        # Plot title
        ctk.CTkLabel(plot_frame, text="Title:").grid(
            row=2, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_title_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_title, width=150)
        plot_title_entry.grid(row=2, column=1, pady=5, padx=5, sticky="w")
        
        # Plot color
        ctk.CTkLabel(plot_frame, text="Color:").grid(
            row=2, column=2, pady=5, padx=(10, 5), sticky="w"
        )
        plot_color_combo = ctk.CTkComboBox(
            plot_frame,
            values=["blue", "red", "green", "orange", "purple", "brown", "pink", "gray"],
            variable=self.plot_color,
            width=100
        )
        plot_color_combo.grid(row=2, column=3, pady=5, padx=(5, 10), sticky="w")
        
        # X axis settings
        ctk.CTkLabel(plot_frame, text="X Label:").grid(
            row=3, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_xlabel_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_xlabel, width=150)
        plot_xlabel_entry.grid(row=3, column=1, pady=5, padx=5, sticky="w")
        
        # Y axis settings
        ctk.CTkLabel(plot_frame, text="Y Label:").grid(
            row=3, column=2, pady=5, padx=(10, 5), sticky="w"
        )
        plot_ylabel_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_ylabel, width=150)
        plot_ylabel_entry.grid(row=3, column=3, pady=5, padx=(5, 10), sticky="w")
        
        # Units selection
        ctk.CTkLabel(plot_frame, text="Units:").grid(
            row=4, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_units_combo = ctk.CTkComboBox(
            plot_frame,
            values=["nm", "Å", "angstrom"],
            variable=self.plot_units,
            width=100,
            command=self._update_ylabel_units
        )
        plot_units_combo.grid(row=4, column=1, pady=5, padx=5, sticky="w")
        
        # X range settings
        ctk.CTkLabel(plot_frame, text="X Min:").grid(
            row=5, column=0, pady=5, padx=(10, 5), sticky="w"
        )
        plot_xlim_min_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_xlim_min, width=80)
        plot_xlim_min_entry.grid(row=5, column=1, pady=5, padx=5, sticky="w")
        
        ctk.CTkLabel(plot_frame, text="X Max:").grid(
            row=5, column=2, pady=5, padx=(10, 5), sticky="w"
        )
        plot_xlim_max_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_xlim_max, width=80)
        plot_xlim_max_entry.grid(row=5, column=3, pady=5, padx=(5, 10), sticky="w")
        
        # Y range settings
        ctk.CTkLabel(plot_frame, text="Y Min:").grid(
            row=6, column=0, pady=(5, 10), padx=(10, 5), sticky="w"
        )
        plot_ylim_min_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_ylim_min, width=80)
        plot_ylim_min_entry.grid(row=6, column=1, pady=(5, 10), padx=5, sticky="w")
        
        ctk.CTkLabel(plot_frame, text="Y Max:").grid(
            row=6, column=2, pady=(5, 10), padx=(10, 5), sticky="w"
        )
        plot_ylim_max_entry = ctk.CTkEntry(plot_frame, textvariable=self.plot_ylim_max, width=80)
        plot_ylim_max_entry.grid(row=6, column=3, pady=(5, 10), padx=(5, 10), sticky="w")
        
        # Store references for enabling/disabling
        self.plot_widgets = [
            plot_title_entry, plot_color_combo, plot_xlabel_entry, plot_ylabel_entry, plot_units_combo,
            plot_xlim_min_entry, plot_xlim_max_entry, plot_ylim_min_entry, plot_ylim_max_entry
        ]
        
        # Initial state
        self._toggle_plot_settings()
    
    def _toggle_plot_settings(self):
        """Enable/disable plot settings based on checkbox."""
        state = "normal" if self.plot_enabled.get() else "disabled"
        for widget in self.plot_widgets:
            widget.configure(state=state)
    
    def _toggle_rmsf_plot_widgets(self, show=True):
        """Show or hide RMSF-specific plot widgets."""
        if not hasattr(self, 'rmsf_plot_widgets'):
            return
        
        if show:
            for widget in self.rmsf_plot_widgets:
                try:
                    widget.winfo_exists()
                    widget.grid()
                except (AttributeError, tk.TclError):
                    # Widget doesn't exist or has been destroyed, skip
                    pass
        else:
            for widget in self.rmsf_plot_widgets:
                try:
                    widget.winfo_exists()
                    widget.grid_remove()
                except (AttributeError, tk.TclError):
                    # Widget doesn't exist or has been destroyed, skip
                    pass
    
    def _toggle_time_units_widgets(self, show=True):
        """Show or hide time units widgets (X-Axis Units for RMSD)."""
        if not hasattr(self, 'time_units_widgets'):
            return
        
        if show:
            for widget in self.time_units_widgets:
                try:
                    widget.winfo_exists()
                    widget.grid()
                except (AttributeError, tk.TclError):
                    # Widget doesn't exist or has been destroyed, skip
                    pass
        else:
            for widget in self.time_units_widgets:
                try:
                    widget.winfo_exists()
                    widget.grid_remove()
                except (AttributeError, tk.TclError):
                    # Widget doesn't exist or has been destroyed, skip
                    pass
    
    def _update_ylabel_units(self, selected_unit=None):
        """Update Y-axis label when units change."""
        unit = self.plot_units.get()
        if unit in ["Å", "angstrom"]:
            unit_label = "Å"
        else:
            unit_label = "nm"
        
        # Update the Y-axis label
        current_ylabel = self.plot_ylabel.get()
        if "RMSD" in current_ylabel:
            self.plot_ylabel.set(f"RMSD ({unit_label})")
    
    def create_results_section(self, row):
        """Create results display section."""
        # Results frame - spans both columns initially
        results_frame = ctk.CTkFrame(self, fg_color=COLOR_SCHEME["content_bg"])
        results_frame.grid(row=row, column=0, sticky="nsew", padx=(10, 5), pady=5)
        results_frame.grid_columnconfigure(0, weight=1)
        results_frame.grid_rowconfigure(1, weight=1)
        
        # Title
        results_title = ctk.CTkLabel(
            results_frame,
            text="Analysis Results",
            font=ctk.CTkFont(size=14, weight="bold")
        )
        results_title.grid(row=0, column=0, pady=(10, 5), sticky="w", padx=10)
        
        # Results text area
        self.results_textbox = ctk.CTkTextbox(
            results_frame,
            font=ctk.CTkFont(family="Consolas", size=10),
            wrap="word"
        )
        self.results_textbox.grid(row=1, column=0, sticky="nsew", padx=10, pady=(5, 10))
        
        # Plot frame - initially hidden
        self.plot_main_frame = ctk.CTkFrame(self, fg_color=COLOR_SCHEME["content_bg"])
        self.plot_main_frame.grid(row=row, column=1, sticky="nsew", padx=(5, 10), pady=5)
        self.plot_main_frame.grid_columnconfigure(0, weight=1)
        self.plot_main_frame.grid_rowconfigure(1, weight=1)
        
        # Plot title
        plot_display_title = ctk.CTkLabel(
            self.plot_main_frame,
            text="Plot Display",
            font=ctk.CTkFont(size=14, weight="bold")
        )
        plot_display_title.grid(row=0, column=0, pady=(10, 5), sticky="w", padx=10)
        
        # Plot placeholder frame - make it larger and more visible
        self.plot_display_frame = ctk.CTkFrame(self.plot_main_frame, height=400, width=600)
        self.plot_display_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=(5, 10))
        self.plot_display_frame.grid_propagate(False)  # Maintain fixed size initially
        
        # Add placeholder label for plot
        self.plot_placeholder = ctk.CTkLabel(
            self.plot_display_frame,
            text="Plot will appear here after running analysis...",
            font=ctk.CTkFont(size=12),
            text_color="gray"
        )
        self.plot_placeholder.place(relx=0.5, rely=0.5, anchor="center")
        
        # Show plot frame by default (not hidden)
        # self.plot_main_frame.grid_remove()  # Comment this out to show plot area
        
        # Add placeholder text
        placeholder_text = "Analysis results will appear here...\n\n"
        if MDANALYSIS_AVAILABLE:
            placeholder_text += "MDAnalysis is available for trajectory analysis."
        else:
            placeholder_text += "Please install MDAnalysis to enable analysis: conda install -c conda-forge mdanalysis"
        
        self.results_textbox.insert("0.0", placeholder_text)
        self.results_textbox.configure(state="disabled")
    
    def create_control_section(self, row):
        """Create control buttons section."""
        # Control buttons frame - spans both columns
        control_frame = ctk.CTkFrame(self, fg_color="transparent")
        control_frame.grid(row=row, column=0, columnspan=2, sticky="ew", padx=10, pady=10)
        control_frame.grid_columnconfigure(0, weight=1)
        
        button_frame = ctk.CTkFrame(control_frame, fg_color="transparent")
        button_frame.pack()
        
        # Run analysis button
        self.run_button = ctk.CTkButton(
            button_frame,
            text="Run Analysis",
            command=self.run_analysis,
            font=ctk.CTkFont(size=14, weight="bold"),
            height=40,
            width=150
        )
        self.run_button.grid(row=0, column=0, padx=5)
        
        # Clear button
        clear_button = ctk.CTkButton(
            button_frame,
            text="Clear",
            command=self.clear_all,
            height=40,
            width=100
        )
        clear_button.grid(row=0, column=1, padx=5)
        
        # Save results button
        save_button = ctk.CTkButton(
            button_frame,
            text="Save Results",
            command=self.save_results,
            height=40,
            width=120
        )
        save_button.grid(row=0, column=2, padx=5)
        
        # Export data button
        export_data_button = ctk.CTkButton(
            button_frame,
            text="Export Data",
            command=self.export_data,
            height=40,
            width=120
        )
        export_data_button.grid(row=0, column=3, padx=5)
    
    def load_current_pdb(self):
        """Load the current PDB file from the main application."""
        if self.get_current_pdb:
            pdb_file = self.get_current_pdb()
            if pdb_file and os.path.exists(pdb_file):
                self.topology_file.set(pdb_file)
                self.status_callback(f"Loaded current PDB: {Path(pdb_file).name}")
            else:
                messagebox.showwarning("No PDB", "No valid PDB file available in the current project.")
        else:
            messagebox.showwarning("Not Available", "Current PDB loading not available.")
    
    def select_topology_file(self):
        """Select topology file."""
        filetypes = [
            ("PDB files", "*.pdb"),
            ("PSF files", "*.psf"),
            ("All files", "*.*")
        ]
        
        filename = filedialog.askopenfilename(
            title="Select Topology File",
            filetypes=filetypes,
            initialdir=self.initial_directory
        )
        
        if filename:
            self.topology_file.set(filename)
            self.status_callback(f"Topology file selected: {Path(filename).name}")
    
    def add_trajectory_file(self, file_type="dcd"):
        """Add trajectory file(s) to the list.
        
        Args:
            file_type: Type of files to add ("dcd" for trajectories, "log" for NAMD logs)
        """
        if file_type == "log":
            filetypes = [
                ("Log files", "*.log"),
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ]
            title = "Select NAMD Log File(s)"
        else:
            filetypes = [
                ("DCD files", "*.dcd"),
                ("XTC files", "*.xtc"),
                ("TRR files", "*.trr"),
                ("NetCDF files", "*.nc"),
                ("Log files", "*.log"),
                ("All files", "*.*")
            ]
            title = "Select Trajectory File(s)"
        
        filenames = filedialog.askopenfilenames(
            title=title,
            filetypes=filetypes,
            initialdir=self.initial_directory
        )
        
        if filenames:
            added_count = 0
            duplicate_count = 0
            first_log_added = False
            
            for filename in filenames:
                if filename not in self.trajectory_files:
                    self.trajectory_files.append(filename)
                    # Initialize time to 0 for new files
                    if filename not in self.file_times:
                        self.file_times[filename] = 0.0
                    added_count += 1
                    
                    # If this is a log file and it's the first one, detect columns
                    if file_type == "log" and added_count == 1:
                        first_log_added = True
                        self._detect_namd_columns_from_file(filename)
                else:
                    duplicate_count += 1
            
            # Refresh the file list display
            if hasattr(self, 'file_list_frame'):
                self._refresh_file_list_display()
            
            if added_count > 0:
                file_type_name = "log" if file_type == "log" else "trajectory"
                files_text = "file" if added_count == 1 else "files"
                self.status_callback(f"{added_count} {file_type_name} {files_text} added")
                
                # Mark analysis as changed for trajectory files (not log files in Energetic tab)
                if file_type != "log":
                    self._mark_analysis_changed()
                
                # If columns were detected, update the UI
                if first_log_added and hasattr(self, 'namd_available_columns') and self.namd_available_columns:
                    self.status_callback(f"Detected {len(self.namd_available_columns)} columns from log file")
                    # Refresh the NAMD options if they're currently shown
                    if self.analysis_type.get() == "NAMD Log Analysis":
                        self._refresh_namd_options()
            
            if duplicate_count > 0:
                files_text = "file" if duplicate_count == 1 else "files"
                messagebox.showwarning("Duplicate Files", f"{duplicate_count} {files_text} already in the list.")
    
    def _refresh_file_list_display(self):
        """Refresh the file list display with draggable horizontal bars for each file."""
        # Clear existing widgets
        for widget in self.file_list_frame.winfo_children():
            widget.destroy()
        
        # Header row
        header_frame = ctk.CTkFrame(self.file_list_frame, fg_color="transparent")
        header_frame.grid(row=0, column=0, padx=5, pady=(5, 2), sticky="ew")
        header_frame.grid_columnconfigure(1, weight=1)
        
        ctk.CTkLabel(
            header_frame,
            text=":::",
            font=ctk.CTkFont(size=14, weight="bold"),
            width=25,
            text_color=("gray10", "gray90")
        ).grid(row=0, column=0, padx=2)
        
        ctk.CTkLabel(
            header_frame,
            text="Filename",
            font=ctk.CTkFont(weight="bold", size=12),
            anchor="w",
            text_color=("gray10", "gray90")
        ).grid(row=0, column=1, padx=5, sticky="w")
        
        ctk.CTkLabel(
            header_frame,
            text="Time (ns)",
            font=ctk.CTkFont(weight="bold", size=12),
            anchor="center",
            text_color=("gray10", "gray90")
        ).grid(row=0, column=2, padx=5)
        
        # Add each file as a draggable horizontal bar
        for idx, filepath in enumerate(self.trajectory_files, start=1):
            self._create_draggable_file_bar(idx, filepath)
    
    def _create_draggable_file_bar(self, row_idx, filepath):
        """Create a draggable horizontal bar for a file entry."""
        # Main container frame for the file (horizontal bar)
        file_bar = ctk.CTkFrame(
            self.file_list_frame,
            fg_color=("gray85", "gray25"),
            corner_radius=6,
            height=35
        )
        file_bar.grid(row=row_idx, column=0, padx=5, pady=2, sticky="ew")
        file_bar.grid_columnconfigure(1, weight=1)
        file_bar.grid_propagate(False)
        
        # Store reference to the file bar for selection highlighting
        if not hasattr(self, 'file_bars'):
            self.file_bars = {}
        self.file_bars[filepath] = file_bar
        
        # Drag handle (≡ icon - more compatible than ☰)
        drag_handle = ctk.CTkLabel(
            file_bar,
            text=":::",
            font=ctk.CTkFont(size=16),
            width=25,
            cursor="hand2",
            text_color=("gray40", "gray60")
        )
        drag_handle.grid(row=0, column=0, padx=5, sticky="ns")
        
        # Filename label
        filename_label = ctk.CTkLabel(
            file_bar,
            text=Path(filepath).name,
            anchor="w",
            text_color=("gray10", "gray90"),
            font=ctk.CTkFont(size=11)
        )
        filename_label.grid(row=0, column=1, padx=5, sticky="ew")
        
        # Time input entry
        time_var = tk.StringVar(value=str(self.file_times.get(filepath, 0.0)))
        time_entry = ctk.CTkEntry(
            file_bar,
            textvariable=time_var,
            width=80,
            height=28
        )
        time_entry.grid(row=0, column=2, padx=5, pady=3)
        
        # Update time dict when entry changes
        time_entry.bind('<FocusOut>', lambda e, fp=filepath, tv=time_var: self._update_file_time(fp, tv))
        time_entry.bind('<Return>', lambda e, fp=filepath, tv=time_var: self._update_file_time(fp, tv))
        
        # Enable selection on click (on filename area only, not drag handle)
        filename_label.bind('<Button-1>', lambda e: self._select_file_bar(filepath))
        
        # Highlight if this is the currently selected file
        if hasattr(self, 'selected_file') and self.selected_file == filepath:
            file_bar.configure(fg_color=("lightblue", "darkblue"))
        
        # Enable drag-and-drop
        self._make_draggable(file_bar, drag_handle, filepath, row_idx)
    
    def _select_file_bar(self, filepath):
        """Select a file bar and highlight it."""
        # Deselect previous selection
        if hasattr(self, 'selected_file') and self.selected_file in self.file_bars:
            old_bar = self.file_bars[self.selected_file]
            old_bar.configure(fg_color=("gray85", "gray25"))
        
        # Select new file
        self.selected_file = filepath
        if filepath in self.file_bars:
            self.file_bars[filepath].configure(fg_color=("lightblue", "darkblue"))
        
        self.status_callback(f"Selected: {Path(filepath).name}")
    
    def _make_draggable(self, file_bar, drag_handle, filepath, row_idx):
        """Make a file bar draggable for reordering."""
        drag_data = {'y': 0, 'item': None, 'original_color': None}
        
        def on_drag_start(event):
            drag_data['y'] = event.y_root
            drag_data['item'] = filepath
            # Change color to indicate dragging
            drag_data['original_color'] = file_bar.cget('fg_color')
            file_bar.configure(fg_color=("gray75", "gray35"))
            drag_handle.configure(text_color=("blue", "lightblue"))
        
        def on_drag_motion(event):
            # Calculate how much we've moved
            delta_y = event.y_root - drag_data['y']
            
            # Visual feedback: change opacity or color based on movement
            if abs(delta_y) > 40:  # Moved enough to reorder
                file_bar.configure(fg_color=("gray70", "gray40"))
        
        def on_drag_end(event):
            # Restore original color
            if drag_data['original_color']:
                file_bar.configure(fg_color=drag_data['original_color'])
            drag_handle.configure(text_color=("gray40", "gray60"))
            
            # Calculate new position
            delta_y = event.y_root - drag_data['y']
            
            # Each row is approximately 39 pixels (35 height + 4 padding)
            rows_moved = round(delta_y / 39)
            
            if rows_moved != 0:
                self._reorder_file(filepath, rows_moved)
        
        # Bind drag events to both the handle and the entire bar
        for widget in [drag_handle, file_bar]:
            widget.bind('<Button-1>', on_drag_start)
            widget.bind('<B1-Motion>', on_drag_motion)
            widget.bind('<ButtonRelease-1>', on_drag_end)
    
    def _reorder_file(self, filepath, rows_moved):
        """Reorder a file in the trajectory list."""
        try:
            current_idx = self.trajectory_files.index(filepath)
            new_idx = current_idx + rows_moved
            
            # Clamp to valid range
            new_idx = max(0, min(new_idx, len(self.trajectory_files) - 1))
            
            if new_idx != current_idx:
                # Remove from current position
                self.trajectory_files.pop(current_idx)
                # Insert at new position
                self.trajectory_files.insert(new_idx, filepath)
                
                # Refresh display
                self._refresh_file_list_display()
                
                # Mark analysis as changed
                self._mark_analysis_changed()
                
                self.status_callback(f"Moved {Path(filepath).name} to position {new_idx + 1}")
        except (ValueError, IndexError) as e:
            print(f"Error reordering file: {e}")
    
    def _update_file_time(self, filepath, time_var):
        """Update the time value for a specific file."""
        try:
            time_value = float(time_var.get())
            # Only mark as changed if the value actually changed
            if self.file_times.get(filepath, 0.0) != time_value:
                self.file_times[filepath] = time_value
                self._mark_analysis_changed()
        except ValueError:
            messagebox.showerror("Invalid Time", "Please enter a valid number for simulation time.")
            time_var.set(str(self.file_times.get(filepath, 0.0)))
    
    def remove_selected_file(self):
        """Remove the currently selected file from the list."""
        if not self.trajectory_files:
            messagebox.showwarning("No Files", "No files to remove.")
            return
        
        # Check if a file is selected
        if hasattr(self, 'selected_file') and self.selected_file and self.selected_file in self.trajectory_files:
            filepath = self.selected_file
            filename = Path(filepath).name
            
            # Remove from list
            self.trajectory_files.remove(filepath)
            
            # Remove from time dict
            if filepath in self.file_times:
                del self.file_times[filepath]
            
            # Clear selection
            self.selected_file = None
            
            # Refresh display
            self._refresh_file_list_display()
            
            # Mark analysis as changed
            self._mark_analysis_changed()
            
            self.status_callback(f"Removed: {filename}")
        else:
            messagebox.showinfo("No Selection", "Please click on a file to select it, then click Remove.")
    
    def remove_trajectory_file(self):
        """Legacy method for compatibility - redirects to remove_selected_file."""
        self.remove_selected_file()
    
    def _auto_detect_file_times(self):
        """Auto-detect simulation times from DCD/LOG files."""
        if not self.trajectory_files:
            messagebox.showinfo("No Files", "Please add files first.")
            return
        
        try:
            current_tab = self.current_analysis_tab.get()
            cumulative_time = 0.0
            
            for filepath in self.trajectory_files:
                if current_tab == "Structural":
                    # For DCD files, read number of frames and estimate time
                    try:
                        if not MDANALYSIS_AVAILABLE:
                            raise ImportError("MDAnalysis not available")
                        
                        topology_file = self.topology_file.get()
                        if not topology_file or not os.path.exists(topology_file):
                            raise ValueError("Valid topology file required for auto-detection")
                        
                        u = mda.Universe(topology_file, filepath)  # type: ignore[possibly-unbound]
                        n_frames = len(u.trajectory)
                        # Estimate: assume 2fs timestep, 5000 step output frequency = 10ps per frame
                        time_per_frame_ns = 0.01  # 10ps = 0.01ns
                        file_duration = n_frames * time_per_frame_ns
                        
                        self.file_times[filepath] = cumulative_time + file_duration
                        cumulative_time += file_duration
                        
                    except Exception as e:
                        logger.warning(f"Could not auto-detect time for {Path(filepath).name}: {e}")
                        # Use default estimate
                        self.file_times[filepath] = cumulative_time + 10.0  # Default 10ns
                        cumulative_time += 10.0
                        
                else:  # Energetic
                    # For LOG files, read actual timesteps
                    try:
                        with open(filepath, 'r') as f:
                            lines = f.readlines()
                        
                        timesteps = []
                        for line in lines:
                            if line.startswith("ENERGY:"):
                                parts = line.split()
                                if len(parts) > 1:
                                    try:
                                        timesteps.append(int(parts[1]))
                                    except (ValueError, IndexError):
                                        pass
                        
                        if timesteps:
                            # Calculate duration (timesteps * dt)
                            # Assume 2fs timestep
                            dt_fs = 2.0
                            max_step = max(timesteps)
                            min_step = min(timesteps)
                            duration_fs = (max_step - min_step) * dt_fs
                            duration_ns = duration_fs / 1e6
                            
                            self.file_times[filepath] = cumulative_time + duration_ns
                            cumulative_time += duration_ns
                        else:
                            # Default
                            self.file_times[filepath] = cumulative_time + 10.0
                            cumulative_time += 10.0
                            
                    except Exception as e:
                        logger.warning(f"Could not auto-detect time for {Path(filepath).name}: {e}")
                        self.file_times[filepath] = cumulative_time + 10.0
                        cumulative_time += 10.0
            
            self._refresh_file_list_display()
            messagebox.showinfo("Auto-Detection Complete", "Simulation times have been estimated from file data.")
            
        except Exception as e:
            messagebox.showerror("Auto-Detection Failed", f"Could not auto-detect times:\n\n{str(e)}")
    
    def _calculate_time_array_for_trajectories(self, universe):
        """Calculate proper time array based on user-specified file times and frame counts.
        
        Each file's time value represents the DURATION (simulation time) of that file in ns.
        Files are stacked sequentially: file1 runs from 0 to its duration,
        file2 runs from file1's end to file1's end + file2's duration, etc.
        
        Args:
            universe: MDAnalysis Universe object (can be None, will use self.trajectory_files)
            
        Returns:
            numpy array of time values in nanoseconds
        """
        time_array = []
        cumulative_time_ns = 0.0  # Track cumulative time across all files
        
        # Load each trajectory separately to get frame counts
        topology_file = self.topology_file.get()
        
        for traj_file in self.trajectory_files:
            # Load trajectory to count frames
            temp_u = mda.Universe(topology_file, traj_file)  # type: ignore[possibly-unbound]
            n_frames = len(temp_u.trajectory)
            
            # Get the simulation time for this file (this is the DURATION in ns)
            duration_ns = self.file_times.get(traj_file, 0.0)
            
            if duration_ns > 0 and n_frames > 1:
                # Create linearly spaced time points for this trajectory
                # Start from cumulative_time_ns, end at cumulative_time_ns + duration_ns
                file_times = np.linspace(cumulative_time_ns, cumulative_time_ns + duration_ns, n_frames)
                cumulative_time_ns += duration_ns
            else:
                # Fallback: use frame indices (assume 0.01 ns = 10 ps per frame)
                file_times = cumulative_time_ns + np.arange(n_frames) * 0.01  # 10ps per frame in ns
                cumulative_time_ns += n_frames * 0.01
            
            time_array.extend(file_times)
        
        return np.array(time_array)
        return np.array(time_array)

    def add_trajectory_folder(self, file_type="dcd"):
        """Add trajectory files from a selected folder.
        
        Args:
            file_type: Type of files to search for ("dcd" for trajectories, "log" for NAMD logs)
        """
        title = "Select Equilibration Folder" if file_type == "dcd" else "Select Folder with Log Files"
        directory = filedialog.askdirectory(
            title=title,
            initialdir=self.initial_directory
        )
        
        if directory:
            found_files = []
            
            if file_type == "log":
                # Find all .log files in the directory
                for file_path in Path(directory).glob("*.log"):
                    found_files.append(file_path)
            else:
                # Find all trajectory files (.dcd, .xtc, .trr, .nc) in the directory
                for ext in ["*.dcd", "*.xtc", "*.trr", "*.nc"]:
                    for file_path in Path(directory).glob(ext):
                        found_files.append(file_path)
                
                # Also store discovered log files for later use
                log_files = []
                for file_path in Path(directory).glob("*.log"):
                    log_files.append(file_path)
                self.discovered_log_files = {
                    file_path.stem: str(file_path) for file_path in log_files
                }
            
            # Sort files using natural/numerical sorting
            # This ensures proper order like step4, step5, step6, step7, step7_production
            def natural_sort_key(path):
                """Extract numbers from filename for natural sorting."""
                # Split filename into text and number parts
                parts = re.split(r'(\d+)', path.name.lower())
                # Convert numeric parts to integers for proper sorting
                return [int(part) if part.isdigit() else part for part in parts]
            
            found_files.sort(key=natural_sort_key)
            
            if found_files:
                added_count = 0
                duplicate_count = 0
                
                for file_path in found_files:
                    file_str = str(file_path)
                    if file_str not in self.trajectory_files:
                        self.trajectory_files.append(file_str)
                        # Initialize time to 0 for new files
                        if file_str not in self.file_times:
                            self.file_times[file_str] = 0.0
                        added_count += 1
                    else:
                        duplicate_count += 1
                
                # Refresh the file list display
                if hasattr(self, 'file_list_frame'):
                    self._refresh_file_list_display()
                
                if added_count > 0:
                    file_type_name = "log" if file_type == "log" else "trajectory"
                    files_text = "file" if added_count == 1 else "files"
                    self.status_callback(f"Added {added_count} {file_type_name} {files_text} from folder")
                    
                    # Mark analysis as changed for trajectory files (not log files)
                    if file_type != "log":
                        self._mark_analysis_changed()
                    
                    # Auto-discover topology file if present (only for trajectory files)
                    if file_type == "dcd":
                        self._auto_discover_topology(directory)
                
                if duplicate_count > 0:
                    files_text = "file" if duplicate_count == 1 else "files"
                    messagebox.showinfo("Info", f"{duplicate_count} {files_text} already in the list.")
            else:
                file_ext = ".log" if file_type == "log" else ".dcd/trajectory"
                messagebox.showwarning("No Files Found", f"No {file_ext} files found in the selected directory.")

    def _auto_discover_topology(self, directory):
        """Auto-discover topology file in the selected directory."""
        topology_extensions = ['.pdb', '.psf', '.prmtop']
        
        for ext in topology_extensions:
            topology_files = list(Path(directory).glob(f"*{ext}"))
            if topology_files:
                # Use the first found topology file
                topology_file = str(topology_files[0])
                if not self.topology_file.get():  # Only set if not already set
                    self.topology_file.set(topology_file)
                    self.status_callback(f"Auto-discovered topology: {topology_files[0].name}")
                break
    
    def select_output_directory(self):
        """Select output directory."""
        directory = filedialog.askdirectory(
            title="Select Output Directory",
            initialdir=self.initial_directory
        )
        if directory:
            self.output_directory.set(directory)
            self.status_callback(f"Output directory selected: {Path(directory).name}")
            # Autodetect topology file
            self.autodetect_topology(directory)
    
    def autodetect_topology(self, directory):
        """Automatically detect topology files in the given directory."""
        topology_files = ["system.pdb", "system.psf", "system.inpcrd"]
        
        for topo_file in topology_files:
            full_path = os.path.join(directory, topo_file)
            if os.path.exists(full_path):
                self.topology_file.set(full_path)
                self.status_callback(f"Auto-detected topology file: {topo_file}")
                return
        
        # If none found, check for any .pdb, .psf, or .inpcrd files
        for ext in [".pdb", ".psf", ".inpcrd"]:
            for file in os.listdir(directory):
                if file.endswith(ext):
                    full_path = os.path.join(directory, file)
                    self.topology_file.set(full_path)
                    self.status_callback(f"Auto-detected topology file: {file}")
                    return
    
    def run_analysis(self):
        """Run the analysis based on selected analysis type."""
        analysis_type = self.analysis_type.get()
        
        # Check if analysis type is selected
        if analysis_type == "-":
            messagebox.showwarning("No Analysis Type", "Please select an analysis type first.")
            return
        
        # Route to appropriate analysis method
        if analysis_type == "NAMD Log Analysis":
            self.run_log_analysis()
        else:
            self.run_trajectory_analysis()
    
    def run_log_analysis(self):
        """Run analysis on NAMD log files."""
        # Check if we have log files loaded
        if not self.trajectory_files and not self.analysis_results:
            messagebox.showwarning("No Log Data", "Please add log files first.")
            return
        
        # Load log files if not already loaded
        if self.trajectory_files and not self.analysis_results:
            try:
                self._auto_load_log_files()
            except Exception as e:
                messagebox.showerror("Loading Failed", f"Failed to load log files:\n\n{str(e)}")
                return
        
        # Check if we have any data
        if not self.analysis_results:
            messagebox.showerror("No Data", "No data available. Please check your log files.")
            return
        
        # Ensure NAMD options UI is created (creates the combo box and checkboxes)
        if not hasattr(self, 'namd_column_combo') or not hasattr(self, 'column_checkboxes'):
            # UI hasn't been created yet - trigger it by ensuring analysis type is set
            if self.analysis_type.get() != "NAMD Log Analysis":
                self.analysis_type.set("NAMD Log Analysis")
                # This will trigger _on_analysis_type_change which creates the UI
                self._on_analysis_type_change("NAMD Log Analysis")
                # Wait a moment for UI to be created
                self.after(100, self.run_log_analysis)
                return
        
        # Get selected columns - either from combo or from checkboxes
        selected_columns = self._get_selected_columns()
        
        # If no columns selected, show warning
        if not selected_columns or (len(selected_columns) == 1 and selected_columns[0] == "No data available"):
            messagebox.showwarning("No Column Selected", "Please select which data to plot from the available columns.")
            return
        
        if not self.output_directory.get():
            messagebox.showwarning("Missing Output", "Please select an output directory.")
            return
        
        # Update plot if enabled
        if self.plot_enabled.get():
            self._plot_namd_data()
        
        # Update status
        num_cols = len(selected_columns)
        col_word = "column" if num_cols == 1 else "columns"
        self.status_callback(f"NAMD log analysis completed ({num_cols} {col_word} plotted)")

    def _load_combined_namd_logs(self):
        """Load and combine data from all discovered NAMD log files."""
        combined_data = {
            'timestep': [],
            'temp': [],
            'poteng': [],
            'kineng': [],
            'toteng': [],
            'pressure': [],
            'volume': []
        }
        
        # Track how many data points came from each file and timestep ranges
        self._file_data_counts = {}
        self._file_timestep_ranges = {}  # Store (start_idx, end_idx, min_ts, max_ts) for each file
        
        # Use self.trajectory_files order to maintain consistency with time assignments
        # This ensures that the order matches what the user sees in the UI
        print(f"Loading combined NAMD logs from {len(self.trajectory_files)} files in UI order")
        
        current_idx = 0
        for log_path in self.trajectory_files:
            try:
                basename = Path(log_path).stem
                print(f"Parsing log file: {log_path}")
                file_data = self.parse_namd_log(log_path)
                if file_data:
                    num_points = len(file_data.get('timestep', []))
                    print(f"  Got {num_points} data points from {basename}")
                    
                    # Track data count for this file
                    self._file_data_counts[log_path] = num_points
                    
                    # Track timestep range for this file
                    if num_points > 0 and 'timestep' in file_data:
                        timesteps = file_data['timestep']
                        min_ts = min(timesteps)
                        max_ts = max(timesteps)
                        self._file_timestep_ranges[log_path] = (current_idx, current_idx + num_points, min_ts, max_ts)
                        print(f"  Timestep range: {min_ts} to {max_ts}")
                        current_idx += num_points
                    
                    # Append data from this log file
                    for key in combined_data.keys():
                        if key in file_data and file_data[key]:
                            combined_data[key].extend(file_data[key])
                else:
                    print(f"  No data returned from {basename}")
            except Exception as e:
                print(f"  Error parsing {Path(log_path).stem}: {str(e)}")
                self.status_callback(f"Warning: Could not parse {Path(log_path).name}: {str(e)}")
        
        print(f"Combined data summary: {[(key, len(values)) for key, values in combined_data.items()]}")
        
        if combined_data['timestep']:
            self.analysis_results = combined_data
            self.status_callback(f"Loaded combined NAMD log data ({len(combined_data['timestep'])} total data points)")
        else:
            messagebox.showerror("Error", "No data could be loaded from the discovered log files.")
    
    def run_trajectory_analysis(self):
        """Run the trajectory analysis using MDAnalysis."""
        if not MDANALYSIS_AVAILABLE:
            messagebox.showerror(
                "MDAnalysis Not Available",
                "MDAnalysis is required for analysis. Please install it with:\n\nconda install -c conda-forge mdanalysis"
            )
            return
        
        # Validate inputs
        if not self.topology_file.get():
            messagebox.showerror("Missing Input", "Please select a topology file.")
            return
        
        if not self.trajectory_files:
            messagebox.showerror("Missing Input", "Please add at least one trajectory file.")
            return
        
        if not self.output_directory.get():
            messagebox.showerror("Missing Input", "Please select an output directory.")
            return
        
        if self.running_analysis:
            messagebox.showwarning("Analysis Running", "An analysis is already running.")
            return
        
        # Disable run button
        self.run_button.configure(state="disabled")
        self.running_analysis = True
        
        # Clear results
        if hasattr(self, 'results_info'):
            self.results_info.delete("0.0", "end")
            self.results_info.insert("0.0", "Starting analysis...\n")
        else:
            # Fallback for old widget
            self.results_textbox.configure(state="normal")
            self.results_textbox.delete("0.0", tk.END)
            self.results_textbox.insert("0.0", "Starting analysis...\n")
            self.results_textbox.configure(state="disabled")
        
        # Run analysis in background thread
        thread = threading.Thread(target=self._run_analysis_thread)
        thread.daemon = True
        thread.start()
    
    def _run_analysis_thread(self):
        """Run analysis in background thread."""
        try:
            # Load all trajectories
            self._update_results("Loading trajectories...\n")
            
            topology_file = self.topology_file.get()
            
            # Load all trajectory files using MDAnalysis
            if len(self.trajectory_files) == 1:
                # Single trajectory
                self._update_results(f"Loading trajectory: {Path(self.trajectory_files[0]).name}...\n")
                universe = mda.Universe(topology_file, self.trajectory_files[0])  # type: ignore[possibly-unbound]
                total_frames = len(universe.trajectory)
                self._update_results(f"  Loaded {total_frames} frames\n")
            else:
                # Multiple trajectories - MDAnalysis automatically concatenates them
                self._update_results(f"Loading {len(self.trajectory_files)} trajectory files...\n")
                universe = mda.Universe(topology_file, self.trajectory_files)  # type: ignore[possibly-unbound]
                total_frames = len(universe.trajectory)
                self._update_results(f"Successfully loaded and concatenated {len(self.trajectory_files)} trajectories\n")
            
            self._update_results(f"Total trajectory: {total_frames} frames\n")
            
            # Get atom selection
            selection_text = self.atom_selection.get()
            
            # Determine the actual selection string for MDAnalysis
            actual_selection = None
            
            # Use custom selection if specified
            if selection_text == "Custom" and self.use_custom_selection.get():
                custom_sel = self.custom_atom_selection.get()
                if custom_sel:
                    try:
                        atoms = universe.select_atoms(custom_sel)
                        actual_selection = custom_sel
                        self._update_results(f"Using custom selection: '{custom_sel}'\n")
                    except Exception as e:
                        error_msg = f"ERROR: Invalid atom selection: {e}\n"
                        self._update_results(error_msg)
                        raise ValueError(f"Custom atom selection failed: {e}")
                else:
                    error_msg = "ERROR: Custom selection is empty. Please provide a valid selection.\n"
                    self._update_results(error_msg)
                    raise ValueError("Custom atom selection is empty")
            elif selection_text == "protein":
                actual_selection = "protein"
                atoms = universe.select_atoms(actual_selection)
            elif selection_text == "backbone":
                actual_selection = "backbone"
                atoms = universe.select_atoms(actual_selection)
            elif selection_text == "protein and name C":
                actual_selection = "protein and name C"
                atoms = universe.select_atoms(actual_selection)
            elif selection_text == "name CA":
                actual_selection = "name CA"
                atoms = universe.select_atoms(actual_selection)
            elif selection_text == "protein and not name H":
                actual_selection = "protein and not name H*"
                atoms = universe.select_atoms(actual_selection)
            elif selection_text == "all":
                actual_selection = "all"
                atoms = universe.select_atoms(actual_selection)
            else:
                actual_selection = "all"
                atoms = universe.select_atoms(actual_selection)
            
            self._update_results(f"Selected {len(atoms)} atoms for analysis\n")
            
            # Check if selection returned any atoms
            if len(atoms) == 0:
                error_msg = f"ERROR: Atom selection '{selection_text}' returned no atoms. Please check your selection.\n"
                self._update_results(error_msg)
                raise ValueError(f"Atom selection '{selection_text}' returned no atoms")
            
            # Run specific analysis
            analysis_type = self.analysis_type.get()
            
            if analysis_type.upper() == "RMSD":
                self._run_rmsd_analysis(universe, atoms, actual_selection)
            elif analysis_type.upper() == "RMSF":
                self._run_rmsf_analysis(universe, atoms, actual_selection)
            elif analysis_type.upper() == "DISTANCES":
                self._run_distance_analysis(universe)
            elif analysis_type.upper() == "RADIUS OF GYRATION":
                self._run_rg_analysis(universe, atoms)
            else:
                self._update_results(f"Analysis type '{analysis_type}' not implemented yet.\n")
            
            self._update_results("\nAnalysis completed successfully!\n")
            self.status_callback("Analysis completed")
            # Mark analysis as ready/up-to-date
            self._mark_analysis_ready()
            
        except Exception as e:
            error_msg = f"\n{'='*50}\nANALYSIS FAILED\n{'='*50}\n{str(e)}\n{'='*50}\n"
            self._update_results(error_msg)
            logger.error(f"Analysis error: {e}")
            
            # Clear analysis results to prevent plotting old data
            self.analysis_results = {}
            
            # Capture error message before lambda
            error_text = str(e)
            
            # Show error message to user
            self.after(0, lambda err=error_text: messagebox.showerror(
                "Analysis Failed", 
                f"The analysis could not be completed:\n\n{err}\n\nPlease check your input parameters and try again."
            ))
        
        finally:
            # Re-enable run button
            self.run_button.configure(state="normal")
            self.running_analysis = False
    
    def _update_results(self, text):
        """Update results text box from thread."""
        if hasattr(self, 'results_info'):
            self.results_info.insert("end", text)
        else:
            # Fallback to old method if new widget not available
            self.results_textbox.configure(state="normal")
            self.results_textbox.insert(tk.END, text)
            self.results_textbox.see(tk.END)
            self.results_textbox.configure(state="disabled")
            self.results_textbox.update()
    
    def _run_rmsd_analysis(self, universe, atoms, selection_string):
        """Run RMSD analysis.
        
        Args:
            universe: MDAnalysis Universe object
            atoms: Selected AtomGroup
            selection_string: The selection string used to create the AtomGroup
        """
        self._update_results("Computing RMSD...\n")
        
        # Convert reference frame from string to int, default to 0 if empty or invalid
        try:
            ref_frame = int(self.reference_frame.get() or "0")
        except ValueError:
            ref_frame = 0
            
        if ref_frame >= len(universe.trajectory):
            ref_frame = 0
        
        # Set the reference frame
        universe.trajectory[ref_frame]
        ref_coords = atoms.positions.copy()
        
        # Apply alignment if requested
        if self.align_rmsd.get():
            self._update_results("Aligning structures to reference frame...\n")
            # Use MDAnalysis alignment
            from MDAnalysis.analysis import align as mda_align
            aligner = mda_align.AlignTraj(universe, universe, select=selection_string,
                                         ref_frame=ref_frame, in_memory=True)
            aligner.run()
            self._update_results("Alignment completed (rotation + translation applied).\n")
            
            # Calculate RMSD after alignment
            rmsd_values = []
            for ts in universe.trajectory:
                rmsd = rms.rmsd(atoms.positions, ref_coords, superposition=False)  # type: ignore[possibly-unbound]
                rmsd_values.append(rmsd)
            rmsd_values = np.array(rmsd_values)
        else:
            self._update_results("Calculating RMSD without alignment (raw coordinates)...\n")
            # Calculate raw RMSD without any alignment
            rmsd_values = []
            for ts in universe.trajectory:
                diff = atoms.positions - ref_coords
                rmsd = np.sqrt(np.mean(np.sum(diff * diff, axis=1)))
                rmsd_values.append(rmsd)
            rmsd_values = np.array(rmsd_values)
            self._update_results("Raw RMSD calculated (no centering, no rotation).\n")
        
        # MDAnalysis returns RMSD in Angstroms by default
        # Handle unit conversion for RMSD
        unit = self.plot_units.get()
        if unit in ["Å", "angstrom"]:
            rmsd_values_display = rmsd_values  # Already in Å
            unit_label = "Å"
        else:
            rmsd_values_display = rmsd_values / 10.0  # Convert Å to nm
            unit_label = "nm"
        
        # Create time array based on user-specified file times (returns time in ns)
        time_values = self._calculate_time_array_for_trajectories(None)
        
        # Store data for later use (time in ns, RMSD in Å for internal storage)
        self.rmsd_data = rmsd_values  # Store in Å (MDAnalysis default)
        self.time_data = time_values  # Store in ns
        
        # Save results in multiple formats for easy reloading
        output_dir = self.output_directory.get()
        
        # Save as text file
        output_file = os.path.join(output_dir, "rmsd_results.txt")
        with open(output_file, 'w') as f:
            f.write("# Frame\tTime(ns)\tRMSD(Å)\tRMSD(nm)\n")
            for i, rmsd in enumerate(rmsd_values):
                time_ns = time_values[i]
                rmsd_nm = rmsd / 10.0
                f.write(f"{i}\t{time_ns:.4f}\t{rmsd:.4f}\t{rmsd_nm:.4f}\n")
        
        # Save as CSV for easy reloading
        csv_file = os.path.join(output_dir, "rmsd_data.csv")
        with open(csv_file, 'w') as f:
            f.write("Time(ns),RMSD(Å),RMSD(nm)\n")
            for time_ns, rmsd in zip(time_values, rmsd_values):
                rmsd_nm = rmsd / 10.0
                f.write(f"{time_ns:.4f},{rmsd:.4f},{rmsd_nm:.4f}\n")
        
        # Display statistics
        mean_rmsd = np.mean(rmsd_values_display)
        std_rmsd = np.std(rmsd_values_display)
        max_rmsd = np.max(rmsd_values_display)
        
        self._update_results(f"RMSD Statistics:\n")
        self._update_results(f"  Mean: {mean_rmsd:.4f} {unit_label}\n")
        self._update_results(f"  Std:  {std_rmsd:.4f} {unit_label}\n")
        self._update_results(f"  Max:  {max_rmsd:.4f} {unit_label}\n")
        self._update_results(f"Results saved to: {output_file}\n")
        self._update_results(f"CSV data saved to: {csv_file}\n")
        
        # Store analysis results
        self.analysis_results['rmsd'] = rmsd_values  # Store in Å
        self.analysis_results['time'] = time_values
        self.analysis_results['unit'] = unit_label
        
        # Create plot if enabled
        if self.plot_enabled.get():
            self._plot_rmsd_data(rmsd_values, time_values)
    
    def _export_plot(self):
        """Export current plot to file."""
        if not hasattr(self, 'canvas') or not self.canvas:
            messagebox.showerror("No Plot", "No plot available to export.")
            return
        
        filetypes = [
            ("PNG files", "*.png"),
            ("SVG files", "*.svg"),
            ("PDF files", "*.pdf"),
            ("JPEG files", "*.jpg"),
            ("All files", "*.*")
        ]
        
        filename = filedialog.asksaveasfilename(
            title="Export Plot",
            filetypes=filetypes,
            defaultextension=".png",
            initialdir=self.output_directory.get()
        )
        
        if filename:
            try:
                self.canvas.figure.savefig(filename, dpi=300, bbox_inches='tight')
                self._update_results(f"Plot exported to: {filename}\n")
                self.status_callback(f"Plot exported to {Path(filename).name}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export plot:\n{str(e)}")
                logger.error(f"Plot export error: {e}")
    
    def export_csv(self):
        """Export analysis results to CSV format."""
        if not hasattr(self, 'analysis_results') or not self.analysis_results:
            messagebox.showerror("No Data", "No analysis results available to export.")
            return
        
        filetypes = [("CSV files", "*.csv"), ("All files", "*.*")]
        filename = filedialog.asksaveasfilename(
            title="Export CSV",
            filetypes=filetypes,
            defaultextension=".csv",
            initialdir=self.output_directory.get() if self.output_directory.get() else os.getcwd()
        )
        
        if filename:
            try:
                import csv
                with open(filename, 'w', newline='') as f:
                    writer = csv.writer(f)
                    
                    # Write header
                    if 'time' in self.analysis_results and 'rmsd' in self.analysis_results:
                        writer.writerow(['Time(ps)', 'RMSD(nm)'])
                        # Write data
                        for time, rmsd in zip(self.analysis_results['time'], self.analysis_results['rmsd']):
                            writer.writerow([f"{time:.2f}", f"{rmsd:.4f}"])
                    else:
                        # Generic export for other analysis types
                        for key, values in self.analysis_results.items():
                            if hasattr(values, '__len__') and len(values) > 1:
                                writer.writerow([key] + [str(v) for v in values])
                            else:
                                writer.writerow([key, str(values)])
                
                self._update_results(f"Data exported to CSV: {filename}\n")
                self.status_callback(f"CSV exported to {Path(filename).name}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export CSV:\n{str(e)}")
                logger.error(f"CSV export error: {e}")
    
    def export_json(self):
        """Export analysis results to JSON format."""
        if not hasattr(self, 'analysis_results') or not self.analysis_results:
            messagebox.showerror("No Data", "No analysis results available to export.")
            return
        
        filetypes = [("JSON files", "*.json"), ("All files", "*.*")]
        filename = filedialog.asksaveasfilename(
            title="Export JSON",
            filetypes=filetypes,
            defaultextension=".json",
            initialdir=self.output_directory.get() if self.output_directory.get() else os.getcwd()
        )
        
        if filename:
            try:
                import json
                
                # Prepare data for JSON serialization
                export_data = {}
                for key, value in self.analysis_results.items():
                    if hasattr(value, 'tolist'):  # numpy array
                        export_data[key] = value.tolist()
                    elif hasattr(value, '__len__') and not isinstance(value, str):  # list/array
                        export_data[key] = list(value)
                    else:
                        export_data[key] = value
                
                # Add metadata
                export_data['metadata'] = {
                    'export_time': str(datetime.now()),
                    'analysis_type': self.analysis_type.get(),
                    'atom_selection': self.atom_selection.get(),
                    'reference_frame': self.reference_frame.get()
                }
                
                with open(filename, 'w') as f:
                    json.dump(export_data, f, indent=2)
                
                self._update_results(f"Data exported to JSON: {filename}\n")
                self.status_callback(f"JSON exported to {Path(filename).name}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export JSON:\n{str(e)}")
                logger.error(f"JSON export error: {e}")
    
    def export_numpy(self):
        """Export analysis results to NumPy format."""
        if not hasattr(self, 'analysis_results') or not self.analysis_results:
            messagebox.showerror("No Data", "No analysis results available to export.")
            return
        
        filetypes = [("NumPy files", "*.npy"), ("NumPy compressed", "*.npz"), ("All files", "*.*")]
        filename = filedialog.asksaveasfilename(
            title="Export NumPy",
            filetypes=filetypes,
            defaultextension=".npz",
            initialdir=self.output_directory.get() if self.output_directory.get() else os.getcwd()
        )
        
        if filename:
            try:
                import numpy as np
                
                # Prepare data for NumPy export
                arrays_to_save = {}
                for key, value in self.analysis_results.items():
                    if hasattr(value, '__len__') and not isinstance(value, str):
                        arrays_to_save[key] = np.array(value)
                    elif isinstance(value, (int, float)):
                        arrays_to_save[key] = np.array([value])
                
                if filename.endswith('.npz'):
                    # Save as compressed archive
                    np.savez_compressed(filename, **arrays_to_save)
                elif filename.endswith('.npy') and len(arrays_to_save) == 1:
                    # Save single array
                    key, array = next(iter(arrays_to_save.items()))
                    np.save(filename, array)
                else:
                    # Save as uncompressed archive
                    np.savez(filename, **arrays_to_save)
                
                self._update_results(f"Data exported to NumPy: {filename}\n")
                self.status_callback(f"NumPy exported to {Path(filename).name}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export NumPy:\n{str(e)}")
                logger.error(f"NumPy export error: {e}")
    
    def _run_rmsf_analysis(self, universe, atoms, selection_string):
        """Run RMSF analysis.
        
        RMSF (Root Mean Square Fluctuation) measures the time-averaged fluctuation
        of each atom around its mean position. Returns per-atom values.
        
        Args:
            universe: MDAnalysis Universe object
            atoms: Selected AtomGroup
            selection_string: The selection string used to create the AtomGroup
        """
        self._update_results("Computing RMSF...\n")
        
        # Convert reference frame from string to int, default to 0 if empty or invalid
        try:
            ref_frame = int(self.reference_frame.get() or "0")
        except ValueError:
            ref_frame = 0
            
        if ref_frame >= len(universe.trajectory):
            ref_frame = 0
            self._update_results(f"Reference frame {self.reference_frame.get()} out of range, using frame 0\n")
        
        # Check if selection is appropriate for residue-based plotting
        selection_text = self.atom_selection.get()
        if selection_text not in ["name CA", "name CB"]:
            self._update_results("WARNING: For residue-based plotting, select 'name CA' or 'name CB'.\n")
            self._update_results("         Current selection may have multiple atoms per residue.\n")
        
        # First, align the trajectory to the reference frame
        self._update_results(f"Aligning trajectory to reference frame {ref_frame}...\n")
        from MDAnalysis.analysis import align as mda_align
        aligner = mda_align.AlignTraj(universe, universe, select=selection_string,
                                     ref_frame=ref_frame, in_memory=True)
        aligner.run()
        
        # Calculate RMSF using MDAnalysis
        from MDAnalysis.analysis import rms
        rmsfer = rms.RMSF(atoms)
        rmsfer.run()
        rmsf_values = rmsfer.results.rmsf  # MDAnalysis returns RMSF in Angstroms
        
        # Handle unit conversion for RMSF
        unit = self.plot_units.get()
        if unit in ["Å", "angstrom"]:
            rmsf_values_display = rmsf_values  # Already in Å
            unit_label = "Å"
        else:
            rmsf_values_display = rmsf_values / 10.0  # Convert Å to nm
            unit_label = "nm"
        
        # Map atoms to residues for better visualization
        residue_info = []
        for atom in atoms:
            residue = atom.residue
            residue_info.append({
                'atom_index': atom.index,
                'residue_number': residue.resid,
                'residue_name': residue.resname,
                'residue_index': residue.resindex,
                'chain_id': residue.segid
            })
        
        # Save results with residue information
        output_file = os.path.join(self.output_directory.get(), "rmsf_results.txt")
        with open(output_file, 'w') as f:
            f.write(f"# Atom_Index\tResidue_Number\tResidue_Name\tChain\tRMSF(Å)\tRMSF({unit_label})\n")
            for i, rmsf in enumerate(rmsf_values):
                rmsf_display = rmsf_values_display[i]
                res_info = residue_info[i]
                f.write(f"{res_info['atom_index']}\t{res_info['residue_number']}\t"
                       f"{res_info['residue_name']}\t{res_info['chain_id']}\t"
                       f"{rmsf:.4f}\t{rmsf_display:.4f}\n")
        
        # Save as CSV
        csv_file = os.path.join(self.output_directory.get(), "rmsf_data.csv")
        with open(csv_file, 'w') as f:
            f.write(f"Atom_Index,Residue_Number,Residue_Name,Chain,RMSF(Å),RMSF({unit_label})\n")
            for i, rmsf in enumerate(rmsf_values):
                rmsf_display = rmsf_values_display[i]
                res_info = residue_info[i]
                f.write(f"{res_info['atom_index']},{res_info['residue_number']},"
                       f"{res_info['residue_name']},{res_info['chain_id']},"
                       f"{rmsf:.4f},{rmsf_display:.4f}\n")
        
        # Display statistics
        mean_rmsf = np.mean(rmsf_values_display)
        std_rmsf = np.std(rmsf_values_display)
        max_rmsf = np.max(rmsf_values_display)
        
        self._update_results(f"RMSF Statistics:\n")
        self._update_results(f"  Mean: {mean_rmsf:.4f} {unit_label}\n")
        self._update_results(f"  Std:  {std_rmsf:.4f} {unit_label}\n")
        self._update_results(f"  Max:  {max_rmsf:.4f} {unit_label}\n")
        self._update_results(f"Results saved to: {output_file}\n")
        self._update_results(f"CSV data saved to: {csv_file}\n")
        
        # Store analysis results with residue mapping
        self.analysis_results['rmsf'] = rmsf_values  # Store in Å
        self.analysis_results['atom_indices'] = [a.index for a in atoms]
        self.analysis_results['residue_info'] = residue_info
        self.analysis_results['unit'] = unit_label
        
        # Create plot if enabled
        if self.plot_enabled.get():
            self._plot_rmsf_data(rmsf_values, residue_info)
    
    def _run_distance_analysis(self, universe):
        """Run distance analysis between two atom selections."""
        self._update_results("Computing distances...\n")
        
        # Get user selections
        selection1_text = self.distance_selection1.get().strip()
        selection2_text = self.distance_selection2.get().strip()
        
        if not selection1_text or not selection2_text:
            error_msg = "ERROR: Both selections must be specified for distance analysis.\n"
            self._update_results(error_msg)
            raise ValueError("Distance selections cannot be empty")
        
        try:
            # Parse selections
            atoms1 = universe.select_atoms(selection1_text)
            atoms2 = universe.select_atoms(selection2_text)
            
            self._update_results(f"Selection 1 ('{selection1_text}'): {len(atoms1)} atoms\n")
            self._update_results(f"Selection 2 ('{selection2_text}'): {len(atoms2)} atoms\n")
            
            if len(atoms1) == 0 or len(atoms2) == 0:
                error_msg = "ERROR: One or both selections returned no atoms.\n"
                self._update_results(error_msg)
                raise ValueError("Atom selections must contain at least one atom")
            
            # Calculate geometric centers for each frame
            distances = []
            
            for ts in universe.trajectory:
                # Get coordinates for this frame (in Angstroms)
                coords1 = atoms1.positions  # shape: (n_atoms1, 3)
                coords2 = atoms2.positions  # shape: (n_atoms2, 3)
                
                # Calculate geometric centers
                center1 = np.mean(coords1, axis=0)  # shape: (3,)
                center2 = np.mean(coords2, axis=0)  # shape: (3,)
                
                # Calculate distance between centers
                diff = center2 - center1
                dist = np.sqrt(np.sum(diff * diff))
                distances.append(dist)
            
            distances = np.array(distances)  # MDAnalysis uses Angstroms
            self._update_results("Distance calculation completed.\n")
            
            # Handle unit conversion
            unit = self.plot_units.get()
            if unit in ["Å", "angstrom"]:
                distances_display = distances  # Already in Å
                unit_label = "Å"
            else:
                distances_display = distances / 10.0  # Convert Å to nm
                unit_label = "nm"
            
            # Create time array (in ns)
            time_values = self._calculate_time_array_for_trajectories(universe)
            
            # Store data for plotting (in Å and ns)
            self.analysis_results['distance'] = distances  # In Å
            self.analysis_results['time'] = time_values  # In ns
            
            # Save results
            output_dir = self.output_directory.get()
            
            # Save as text file
            output_file = os.path.join(output_dir, "distance_results.txt")
            with open(output_file, 'w') as f:
                f.write(f"# Distance between geometric centers\n")
                f.write(f"# Selection 1: {selection1_text}\n")
                f.write(f"# Selection 2: {selection2_text}\n")
                f.write("# Frame\tTime(ns)\tDistance(Å)\tDistance(nm)\n")
                for i, dist in enumerate(distances):
                    time_ns = time_values[i]
                    dist_nm = dist / 10.0
                    f.write(f"{i}\t{time_ns:.4f}\t{dist:.4f}\t{dist_nm:.4f}\n")
            
            # Save as CSV
            csv_file = os.path.join(output_dir, "distance_data.csv")
            with open(csv_file, 'w') as f:
                f.write("Time(ns),Distance(Å),Distance(nm)\n")
                for time_ns, dist in zip(time_values, distances):
                    dist_nm = dist / 10.0
                    f.write(f"{time_ns:.4f},{dist:.4f},{dist_nm:.4f}\n")
            
            # Display statistics
            mean_dist = np.mean(distances_display)
            std_dist = np.std(distances_display)
            min_dist = np.min(distances_display)
            max_dist = np.max(distances_display)
            
            self._update_results(f"\nDistance Statistics:\n")
            self._update_results(f"  Mean: {mean_dist:.4f} {unit_label}\n")
            self._update_results(f"  Std:  {std_dist:.4f} {unit_label}\n")
            self._update_results(f"  Min:  {min_dist:.4f} {unit_label}\n")
            self._update_results(f"  Max:  {max_dist:.4f} {unit_label}\n")
            self._update_results(f"Results saved to:\n  {output_file}\n  {csv_file}\n")
            
            # Plot if enabled
            if self.plot_enabled.get():
                self._plot_distance_data()
                
        except Exception as e:
            error_msg = f"ERROR in distance analysis: {str(e)}\n"
            self._update_results(error_msg)
            raise
    
    def _run_rg_analysis(self, universe, atoms):
        """Run radius of gyration analysis."""
        self._update_results("Computing radius of gyration...\n")
        
        # Get user selection
        selection_text = self.rg_selection.get().strip()
        
        if not selection_text:
            error_msg = "ERROR: Atom selection must be specified for Rg analysis.\n"
            self._update_results(error_msg)
            raise ValueError("Rg selection cannot be empty")
        
        try:
            # Parse selection
            rg_atoms = universe.select_atoms(selection_text)
            
            self._update_results(f"Selection ('{selection_text}'): {len(rg_atoms)} atoms\n")
            
            if len(rg_atoms) == 0:
                error_msg = "ERROR: Selection returned no atoms.\n"
                self._update_results(error_msg)
                raise ValueError("Rg atom selection must contain at least one atom")
            
            # Compute radius of gyration for each frame
            rg_values = []
            for ts in universe.trajectory:
                rg = rg_atoms.radius_of_gyration()  # MDAnalysis returns Rg in Angstroms
                rg_values.append(rg)
            
            rg_values = np.array(rg_values)
            self._update_results("Radius of gyration calculation completed.\n")
            
            # Handle unit conversion
            unit = self.plot_units.get()
            if unit in ["Å", "angstrom"]:
                rg_display = rg_values  # Already in Å
                unit_label = "Å"
            else:
                rg_display = rg_values / 10.0  # Convert Å to nm
                unit_label = "nm"
            
            # Create time array (in ns)
            time_values = self._calculate_time_array_for_trajectories(universe)
            
            # Store data for plotting (in Å and ns)
            self.analysis_results['rg'] = rg_values  # In Å
            self.analysis_results['time'] = time_values  # In ns
            
            # Save results
            output_dir = self.output_directory.get()
            
            # Save as text file
            output_file = os.path.join(output_dir, "rg_results.txt")
            with open(output_file, 'w') as f:
                f.write(f"# Radius of Gyration\n")
                f.write(f"# Selection: {selection_text}\n")
                f.write("# Frame\tTime(ns)\tRg(Å)\tRg(nm)\n")
                for i, rg in enumerate(rg_values):
                    time_ns = time_values[i]
                    rg_nm = rg / 10.0
                    f.write(f"{i}\t{time_ns:.4f}\t{rg:.4f}\t{rg_nm:.4f}\n")
            
            # Save as CSV
            csv_file = os.path.join(output_dir, "rg_data.csv")
            with open(csv_file, 'w') as f:
                f.write("Time(ns),Rg(Å),Rg(nm)\n")
                for time_ns, rg in zip(time_values, rg_values):
                    rg_nm = rg / 10.0
                    f.write(f"{time_ns:.4f},{rg:.4f},{rg_nm:.4f}\n")
            
            # Display statistics
            mean_rg = np.mean(rg_display)
            std_rg = np.std(rg_display)
            min_rg = np.min(rg_display)
            max_rg = np.max(rg_display)
            
            self._update_results(f"\nRadius of Gyration Statistics:\n")
            self._update_results(f"  Mean: {mean_rg:.4f} {unit_label}\n")
            self._update_results(f"  Std:  {std_rg:.4f} {unit_label}\n")
            self._update_results(f"  Min:  {min_rg:.4f} {unit_label}\n")
            self._update_results(f"  Max:  {max_rg:.4f} {unit_label}\n")
            self._update_results(f"Results saved to:\n  {output_file}\n  {csv_file}\n")
            
            # Plot if enabled
            if self.plot_enabled.get():
                self._plot_rg_data()
                
        except Exception as e:
            error_msg = f"ERROR in Rg analysis: {str(e)}\n"
            self._update_results(error_msg)
            raise
    
    def clear_all(self):
        """Clear all inputs and results."""
        self.topology_file.set("")
        self.trajectory_files.clear()
        self.trajectory_listbox.delete(0, tk.END)
        self.output_directory.set("")
        self.analysis_type.set("rmsd")
        self.atom_selection.set("protein")
        self.reference_frame.set("0")  # StringVar requires string value
        self.analysis_results.clear()
        
        self.results_textbox.configure(state="normal")
        self.results_textbox.delete("0.0", tk.END)
        placeholder_text = "Analysis results will appear here...\n\n"
        if MDANALYSIS_AVAILABLE:
            placeholder_text += "MDAnalysis is available for trajectory analysis."
        else:
            placeholder_text += "Please install MDAnalysis to enable analysis: conda install -c conda-forge mdanalysis"
        self.results_textbox.insert("0.0", placeholder_text)
        self.results_textbox.configure(state="disabled")
        
        self.status_callback("All fields cleared")
    
    def save_results(self):
        """Save analysis results to file."""
        if not self.results_textbox.get("0.0", tk.END).strip():
            messagebox.showwarning("No Results", "No results to save.")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Save Analysis Results",
            defaultextension=".txt",
            filetypes=[
                ("Text files", "*.txt"),
                ("All files", "*.*")
            ],
            initialdir=self.initial_directory
        )
        
        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write(self.results_textbox.get("0.0", tk.END))
                self.status_callback(f"Results saved to: {Path(filename).name}")
                messagebox.showinfo("Success", f"Results saved to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to save results:\n{str(e)}")
    
    def export_data(self):
        """Export analysis data to various formats."""
        if not self.analysis_results:
            messagebox.showwarning("No Data", "No analysis data to export.")
            return
        
        filetypes = [
            ("CSV files", "*.csv"),
            ("JSON files", "*.json"),
            ("NumPy files", "*.npy"),
            ("All files", "*.*")
        ]
        
        filename = filedialog.asksaveasfilename(
            title="Export Analysis Data",
            filetypes=filetypes,
            defaultextension=".csv",
            initialdir=self.output_directory.get()
        )
        
        if filename:
            try:
                ext = Path(filename).suffix.lower()
                
                if ext == '.csv':
                    self._export_to_csv(filename)
                elif ext == '.json':
                    self._export_to_json(filename)
                elif ext == '.npy':
                    self._export_to_numpy(filename)
                else:
                    # Default to CSV
                    self._export_to_csv(filename)
                
                self.status_callback(f"Data exported to {Path(filename).name}")
                messagebox.showinfo("Success", f"Data exported to:\n{filename}")
                
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export data:\n{str(e)}")
                logger.error(f"Data export error: {e}")
    
    def _export_to_csv(self, filename):
        """Export data to CSV format."""
        import csv
        
        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # Write header
            if 'rmsd' in self.analysis_results and 'time' in self.analysis_results:
                unit_label = self.analysis_results.get('unit', 'nm')
                writer.writerow([f'Time(ps)', f'RMSD({unit_label})'])
                
                # Write data - use display values if available
                time_data = self.analysis_results['time']
                if 'rmsd_display' in self.analysis_results:
                    rmsd_data = self.analysis_results['rmsd_display']
                else:
                    rmsd_data = self.analysis_results['rmsd']
                
                for t, rmsd in zip(time_data, rmsd_data):
                    writer.writerow([f"{t:.2f}", f"{rmsd:.6f}"])
            
            elif 'rmsf' in self.analysis_results:
                writer.writerow(['Atom_Index', 'RMSF(nm)'])
                rmsf_data = self.analysis_results['rmsf']
                for i, rmsf in enumerate(rmsf_data):
                    writer.writerow([i, f"{rmsf:.6f}"])
            
            else:
                # Generic export for other analysis types
                for key, values in self.analysis_results.items():
                    if isinstance(values, np.ndarray):
                        writer.writerow([key])
                        for val in values:
                            writer.writerow([f"{val:.6f}"])
    
    def _export_to_json(self, filename):
        """Export data to JSON format."""
        import json
        
        # Convert numpy arrays to lists for JSON serialization
        export_data = {}
        for key, values in self.analysis_results.items():
            if isinstance(values, np.ndarray):
                export_data[key] = values.tolist()
            else:
                export_data[key] = values
        
        with open(filename, 'w') as f:
            json.dump(export_data, f, indent=2)
    
    def _export_to_numpy(self, filename):
        """Export data to NumPy format."""
        if 'rmsd' in self.analysis_results:
            data_to_save = {
                'rmsd': self.analysis_results['rmsd'],
                'time': self.analysis_results.get('time', np.arange(len(self.analysis_results['rmsd'])))
            }
        else:
            data_to_save = self.analysis_results
        
        np.savez(filename, **data_to_save)

    def _enable_scrolling_for_panel(self, scrollable_frame):
        """Enable mouse wheel scrolling for a CTkScrollableFrame and all its children."""
        def _on_scroll(event):
            # For Windows and MacOS
            if event.delta:
                scrollable_frame._parent_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
            # For Linux
            elif event.num == 4:
                scrollable_frame._parent_canvas.yview_scroll(-1, "units")
            elif event.num == 5:
                scrollable_frame._parent_canvas.yview_scroll(1, "units")
        
        # Bind mouse wheel to the scrollable frame
        scrollable_frame.bind("<MouseWheel>", _on_scroll)  # Windows/Mac
        scrollable_frame.bind("<Button-4>", _on_scroll)    # Linux scroll up
        scrollable_frame.bind("<Button-5>", _on_scroll)    # Linux scroll down
        
        # Recursively bind to all children
        def bind_children(widget):
            for child in widget.winfo_children():
                child.bind("<MouseWheel>", _on_scroll)
                child.bind("<Button-4>", _on_scroll)
                child.bind("<Button-5>", _on_scroll)
                bind_children(child)
        
        # Bind after widgets are created
        scrollable_frame.after(100, lambda: bind_children(scrollable_frame))

    def _bind_mouse_wheel(self, widget):
        """Bind mouse wheel events to scrollable frame."""
        def _on_mouse_wheel(event):
            # Get the widget under the mouse
            current_widget = event.widget
            
            # Check if we're directly on a CTkTextbox or its internal text widget
            current = current_widget
            textbox_widget = None
            
            # First, check if we're inside a CTkTextbox
            while current:
                if isinstance(current, ctk.CTkTextbox):
                    textbox_widget = current
                    break
                # Also check if the widget itself has a _textbox attribute (internal text widget)
                if hasattr(current, 'master') and isinstance(current.master, ctk.CTkTextbox):
                    textbox_widget = current.master
                    break
                current = getattr(current, 'master', None)
            
            # If we found a textbox, scroll it
            if textbox_widget:
                try:
                    if hasattr(textbox_widget, '_textbox'):
                        textbox_widget._textbox.yview_scroll(int(-1 * (event.delta / 120)), "units")
                except Exception:
                    pass
            # Otherwise, use the CTkScrollableFrame's internal scrolling
            elif hasattr(widget, '_parent_canvas'):
                widget._parent_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
        
        def _on_mouse_wheel_linux(event):
            # Get the widget under the mouse
            current_widget = event.widget
            
            # Check if we're directly on a CTkTextbox or its internal text widget
            current = current_widget
            textbox_widget = None
            
            # First, check if we're inside a CTkTextbox
            while current:
                if isinstance(current, ctk.CTkTextbox):
                    textbox_widget = current
                    break
                # Also check if the widget itself has a _textbox attribute (internal text widget)
                if hasattr(current, 'master') and isinstance(current.master, ctk.CTkTextbox):
                    textbox_widget = current.master
                    break
                current = getattr(current, 'master', None)
            
            # If we found a textbox, scroll it
            if textbox_widget:
                try:
                    if hasattr(textbox_widget, '_textbox'):
                        delta = -1 if event.num == 4 else 1
                        textbox_widget._textbox.yview_scroll(delta, "units")
                except Exception:
                    pass
            # Otherwise, scroll the main frame
            elif hasattr(widget, '_parent_canvas'):
                if event.num == 4:
                    widget._parent_canvas.yview_scroll(-1, "units")
                elif event.num == 5:
                    widget._parent_canvas.yview_scroll(1, "units")
        
        # Bind to the widget and all its children
        widget.bind("<MouseWheel>", _on_mouse_wheel)
        widget.bind("<Button-4>", _on_mouse_wheel_linux)
        widget.bind("<Button-5>", _on_mouse_wheel_linux)
        
        # Also bind to all child widgets
        def bind_to_children(parent):
            for child in parent.winfo_children():
                child.bind("<MouseWheel>", _on_mouse_wheel)
                child.bind("<Button-4>", lambda e: widget._parent_canvas.yview_scroll(-1, "units") if hasattr(widget, '_parent_canvas') else None)
                child.bind("<Button-5>", lambda e: widget._parent_canvas.yview_scroll(1, "units") if hasattr(widget, '_parent_canvas') else None)
                bind_to_children(child)
        
        # Use after_idle to ensure all widgets are created
        widget.after_idle(lambda: bind_to_children(widget))

    def _on_analysis_type_change(self, value):
        """Handle analysis type selection change."""
        try:
            # Mark analysis as changed
            self._mark_analysis_changed()
            
            # Clear existing dynamic options
            self._clear_dynamic_options()
            
            # Show appropriate options based on analysis type
            if value == "-":
                # No additional options for default selection
                pass
            elif value in ["RMSD", "RMSF", "Distances", "Radius of Gyration"]:
                self._show_trajectory_analysis_options()
            elif value == "NAMD Log Analysis":
                self._show_namd_log_options()
        except Exception as e:
            logger.error(f"Error in _on_analysis_type_change: {e}", exc_info=True)
            # Don't raise - allow GUI to continue functioning
            return
        
        # Show/hide RMSF-specific plot widgets and time units based on analysis type
        if value == "RMSF":
            # Show RMSF widgets, hide time units
            self._toggle_rmsf_plot_widgets(show=True)
            self._toggle_time_units_widgets(show=False)
            
            # If RMSF data exists, plot it
            if 'rmsf' in self.analysis_results and 'residue_info' in self.analysis_results:
                self._plot_rmsf_data(self.analysis_results['rmsf'], 
                                    self.analysis_results['residue_info'])
            else:
                # Clear plot if no RMSF data available
                self._clear_plot()
        elif value == "RMSD":
            # Hide RMSF widgets, show time units
            self._toggle_rmsf_plot_widgets(show=False)
            self._toggle_time_units_widgets(show=True)
            
            # If RMSD data exists, plot it
            if 'rmsd' in self.analysis_results and 'time' in self.analysis_results:
                self._plot_rmsd_data(self.analysis_results['rmsd'], 
                                    self.analysis_results['time'])
            else:
                # Clear plot if no RMSD data available
                self._clear_plot()
        elif value == "NAMD Log Analysis":
            # Hide RMSF widgets, show time units for energetic analysis
            self._toggle_rmsf_plot_widgets(show=False)
            self._toggle_time_units_widgets(show=True)
            
            # Note: Plot will be triggered by _show_namd_log_options() if data exists
            # No need to explicitly plot here as it's handled in that method
        else:
            # For other analysis types, hide both RMSF widgets and time units
            self._toggle_rmsf_plot_widgets(show=False)
            # Show time units for general structural analysis
            if value in ["Distances", "Radius of Gyration"]:
                self._toggle_time_units_widgets(show=True)
            else:
                self._toggle_time_units_widgets(show=False)
            # Clear plot for other analysis types
            self._clear_plot()
        
        # Update button text and visibility based on analysis type
        self._update_run_button_for_analysis_type(value)
        
        # Update plot labels based on analysis type
        self._update_plot_labels_for_analysis_type(value)
    
    def _update_run_button_for_analysis_type(self, analysis_type):
        """Update run button text and visibility based on analysis type."""
        if hasattr(self, 'run_button'):
            current_tab = self.current_analysis_tab.get()
            
            if current_tab == "Energetic":
                # For Energetic analysis, hide the button (plots are automatic)
                self.run_button.grid_remove()
            elif analysis_type == "NAMD Log Analysis":
                # For NAMD log analysis (shouldn't happen in Energetic, but keep for safety)
                self.run_button.grid()
                self.run_button.configure(text="Plot", command=self._plot_namd_data)
            elif analysis_type in ["RMSD", "RMSF", "Distances", "Radius of Gyration"]:
                # For trajectory analysis, show the RMSD analysis button
                self.run_button.grid()
                self.run_button.configure(text="Run Trajectory Analysis", command=self.run_analysis)
            else:
                # Default state
                self.run_button.grid()
                self.run_button.configure(text="Select Analysis Type", command=None, state="disabled")
    
    def _clear_dynamic_options(self):
        """Clear all dynamic option widgets."""
        for widget in self.dynamic_widgets:
            widget.destroy()
        self.dynamic_widgets.clear()
    
    def _show_trajectory_analysis_options(self):
        """Show options for trajectory-based analysis (RMSD, RMSF, Distances, Radius of Gyration)."""
        row = 0
        analysis_type = self.analysis_type.get()
        
        # Show different options based on analysis type
        if analysis_type == "Distances":
            # Distance analysis: Two atom selections
            # Selection 1
            sel1_label = ctk.CTkLabel(self.dynamic_options_frame, text="Selection 1:")
            sel1_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
            self.dynamic_widgets.append(sel1_label)
            
            sel1_entry = ctk.CTkEntry(
                self.dynamic_options_frame,
                textvariable=self.distance_selection1,
                width=280,
                placeholder_text="e.g., protein and resid 1 to 10"
            )
            sel1_entry.grid(row=row, column=1, columnspan=2, pady=5, padx=(5, 10), sticky="ew")
            self.dynamic_widgets.append(sel1_entry)
            row += 1
            
            # Selection 2
            sel2_label = ctk.CTkLabel(self.dynamic_options_frame, text="Selection 2:")
            sel2_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
            self.dynamic_widgets.append(sel2_label)
            
            sel2_entry = ctk.CTkEntry(
                self.dynamic_options_frame,
                textvariable=self.distance_selection2,
                width=280,
                placeholder_text="e.g., protein and resid 50 to 60"
            )
            sel2_entry.grid(row=row, column=1, columnspan=2, pady=5, padx=(5, 10), sticky="ew")
            self.dynamic_widgets.append(sel2_entry)
            row += 1
            
            # Topology analysis button
            analyze_btn = ctk.CTkButton(
                self.dynamic_options_frame,
                text="Analyze Topology",
                command=self._analyze_topology,
                width=120,
                height=25
            )
            analyze_btn.grid(row=row, column=0, columnspan=3, pady=5, padx=(10, 10), sticky="w")
            self.dynamic_widgets.append(analyze_btn)
            row += 1
            
            # Help text
            help_text = ctk.CTkLabel(
                self.dynamic_options_frame,
                text="Distance is calculated between geometric centers of each selection.\nFor single atoms, use: 'name CA and resid 10'",
                font=ctk.CTkFont(size=10),
                text_color="gray"
            )
            help_text.grid(row=row, column=0, columnspan=3, pady=(0, 5), padx=(10, 10), sticky="w")
            self.dynamic_widgets.append(help_text)
            
        elif analysis_type == "Radius of Gyration":
            # Radius of gyration: Single selection
            rg_label = ctk.CTkLabel(self.dynamic_options_frame, text="Atom Selection:")
            rg_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
            self.dynamic_widgets.append(rg_label)
            
            rg_entry = ctk.CTkEntry(
                self.dynamic_options_frame,
                textvariable=self.rg_selection,
                width=280,
                placeholder_text="e.g., protein, backbone, resid 1 to 100"
            )
            rg_entry.grid(row=row, column=1, columnspan=2, pady=5, padx=(5, 10), sticky="ew")
            self.dynamic_widgets.append(rg_entry)
            row += 1
            
            # Topology analysis button
            analyze_btn = ctk.CTkButton(
                self.dynamic_options_frame,
                text="Analyze Topology",
                command=self._analyze_topology,
                width=120,
                height=25
            )
            analyze_btn.grid(row=row, column=0, columnspan=3, pady=5, padx=(10, 10), sticky="w")
            self.dynamic_widgets.append(analyze_btn)
            row += 1
            
            # Help text
            help_text = ctk.CTkLabel(
                self.dynamic_options_frame,
                text="Calculates radius of gyration for selected atoms over time",
                font=ctk.CTkFont(size=10),
                text_color="gray"
            )
            help_text.grid(row=row, column=0, columnspan=3, pady=(0, 5), padx=(10, 10), sticky="w")
            self.dynamic_widgets.append(help_text)
            
        else:
            # RMSD and RMSF: Original atom selection options
            # Atom selection
            atom_label = ctk.CTkLabel(self.dynamic_options_frame, text="Atom Selection:")
            atom_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
            self.dynamic_widgets.append(atom_label)
            
            atom_combo = ctk.CTkComboBox(
                self.dynamic_options_frame,
                values=["protein", "backbone", "name CA", "protein and name C", "protein and not name H", "all", "Custom"], # Protein and name C is ideal for RMSF analysis as it also considers the ACE and NME caps.
                variable=self.atom_selection,
                width=150,
                command=self._on_atom_selection_change
            )
            atom_combo.grid(row=row, column=1, pady=5, padx=(5, 5), sticky="w")
            self.dynamic_widgets.append(atom_combo)
            
            # Help button for custom selection
            help_btn = ctk.CTkButton(
                self.dynamic_options_frame,
                text="?",
                command=self._show_selection_help,
                width=30,
                height=25
            )
            help_btn.grid(row=row, column=2, pady=5, padx=(0, 5), sticky="w")
            self.dynamic_widgets.append(help_btn)
            
            # Topology analysis button
            analyze_btn = ctk.CTkButton(
                self.dynamic_options_frame,
                text="Analyze Topology",
                command=self._analyze_topology,
                width=120,
                height=25
            )
            analyze_btn.grid(row=row, column=3, pady=5, padx=(0, 10), sticky="w")
            self.dynamic_widgets.append(analyze_btn)
            row += 1
            
            # Custom atom selection entry (hidden by default)
            self.custom_selection_label = ctk.CTkLabel(self.dynamic_options_frame, text="Custom Selection:")
            self.custom_selection_entry = ctk.CTkEntry(
                self.dynamic_options_frame,
                textvariable=self.custom_atom_selection,
                width=280,
                placeholder_text="e.g., protein and resid 1 to 50"
            )
            # Initially hidden
            if self.atom_selection.get() == "Custom":
                self.custom_selection_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
                self.custom_selection_entry.grid(row=row, column=1, columnspan=2, pady=5, padx=(5, 10), sticky="ew")
                self.dynamic_widgets.append(self.custom_selection_label)
                self.dynamic_widgets.append(self.custom_selection_entry)
                row += 1
            
            # Alignment option (only for RMSD)
            if analysis_type == "RMSD":
                align_checkbox = ctk.CTkCheckBox(
                    self.dynamic_options_frame,
                    text="Align structures before RMSD calculation",
                    variable=self.align_rmsd,
                    onvalue=True,
                    offvalue=False
                )
                align_checkbox.grid(row=row, column=0, columnspan=3, pady=5, padx=(10, 10), sticky="w")
                self.dynamic_widgets.append(align_checkbox)
                row += 1
            
            # Reference frame (for both RMSD and RMSF)
            if analysis_type in ["RMSD", "RMSF"]:
                ref_label = ctk.CTkLabel(self.dynamic_options_frame, text="Reference Frame:")
                ref_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
                self.dynamic_widgets.append(ref_label)
                
                ref_entry = ctk.CTkEntry(
                    self.dynamic_options_frame,
                    textvariable=self.reference_frame,
                    width=150
                )
                ref_entry.grid(row=row, column=1, pady=5, padx=(5, 10), sticky="w")
                self.dynamic_widgets.append(ref_entry)
                
                # Bind validation on focus out and Enter key
                ref_entry.bind('<FocusOut>', self._validate_reference_frame)
                ref_entry.bind('<Return>', self._validate_reference_frame)
    
    def _on_atom_selection_change(self, value):
        """Handle atom selection change to show/hide custom selection entry."""
        # Find and remove custom selection widgets if they exist
        for widget in [self.custom_selection_label, self.custom_selection_entry]:
            if widget in self.dynamic_widgets:
                widget.grid_forget()
                self.dynamic_widgets.remove(widget)
        
        if value == "Custom":
            # Find current row position
            row = len([w for w in self.dynamic_widgets if isinstance(w, (ctk.CTkLabel, ctk.CTkCheckBox, ctk.CTkEntry, ctk.CTkComboBox))]) // 2 + 1
            
            self.custom_selection_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
            self.custom_selection_entry.grid(row=row, column=1, columnspan=2, pady=5, padx=(5, 10), sticky="ew")
            self.dynamic_widgets.extend([self.custom_selection_label, self.custom_selection_entry])
            self.use_custom_selection.set(True)
        else:
            self.use_custom_selection.set(False)
    
    def _show_selection_help(self):
        """Show help window for custom atom selection."""
        help_window = ctk.CTkToplevel(self)
        help_window.title("Custom Atom Selection Help")
        help_window.geometry("600x500")
        
        # Title
        title = ctk.CTkLabel(
            help_window,
            text="MDAnalysis Atom Selection Syntax",
            font=ctk.CTkFont(size=16, weight="bold")
        )
        title.pack(pady=(10, 5), padx=10)
        
        # Help text
        help_text = ctk.CTkTextbox(help_window, wrap="word")
        help_text.pack(fill="both", expand=True, padx=10, pady=10)
        
        help_content = """
MDAnalysis uses a powerful and flexible atom selection language. Here are common examples:

BASIC SELECTIONS:
* protein          - All protein atoms
* backbone         - Protein backbone atoms (CA, C, N, O)
* name CA          - All alpha carbon atoms
* all              - All atoms
* water            - All water molecules
* resname ALA      - All alanine residues

COMBINING SELECTIONS:
* protein and backbone        - Protein backbone only
* protein and not name H*     - Protein without hydrogens
* name CA or name CB          - Alpha and beta carbons
* resid 1:50                  - Residues 1 through 50
* protein and resid 10:100    - Protein residues 10-100

RESIDUE SELECTIONS:
* resname ALA GLY VAL         - Specific amino acids
* resid 1 5 10                - Specific residue numbers
* resid 1:50 and name CA      - CA atoms in residues 1-50

SEGMENT/CHAIN SELECTIONS:
* segid A            - Segment A
* segid A B          - Segments A and B

ATOM PROPERTIES:
* type CA            - Atoms of type CA
* mass > 12          - Atoms with mass > 12
* charge < 0         - Negatively charged atoms

EXAMPLES FOR ANALYSIS:
* "protein and backbone" 
  -> Good for RMSD of protein structure
  
* "name CA"
  -> Fast RMSD using only C-alpha atoms
  
* "protein and resid 50:150"
  -> Analyze specific region of protein
  
* "protein and not resname GLY PRO"
  -> Exclude flexible residues
  
* "protein and not name H*"
  -> Protein without hydrogens

For more information, see MDAnalysis documentation:
https://docs.mdanalysis.org/stable/documentation_pages/selections.html
"""
        
        help_text.insert("0.0", help_content)
        help_text.configure(state="disabled")
        
        # Close button
        close_btn = ctk.CTkButton(
            help_window,
            text="Close",
            command=help_window.destroy,
            width=100
        )
        close_btn.pack(pady=(0, 10))
    
    def _analyze_topology(self):
        """Analyze the loaded topology and display fundamental information."""
        if not MDANALYSIS_AVAILABLE:
            messagebox.showerror("MDAnalysis Not Available", "MDAnalysis is required for topology analysis.")
            return
        
        # Check if topology is loaded
        if not self.topology_file.get():
            messagebox.showwarning("No Topology", "Please load a topology file first.")
            return
        
        # Prevent multiple windows from opening
        if hasattr(self, '_topology_window') and self._topology_window:
            try:
                if self._topology_window.winfo_exists():
                    self._topology_window.lift()  # Bring existing window to front
                    self._topology_window.focus()
                    return
            except:
                # Window reference exists but window was destroyed, clear the reference
                self._topology_window = None
        
        try:
            # Load topology using MDAnalysis
            topology_path = self.topology_file.get()
            universe = mda.Universe(topology_path)  # type: ignore[possibly-unbound]
            
            # Type assertion to help linter understand universe is not None
            assert universe is not None, "Failed to load universe"
            
            # Create analysis window (reduced width for better proportions)
            analysis_window = ctk.CTkToplevel(self)
            analysis_window.title("Topology Analysis")
            analysis_window.geometry("800x750")
            
            # Store reference to prevent multiple windows
            self._topology_window = analysis_window
            
            # Initialize font size variable (increased default for better readability)
            self._topology_font_size = 36  # Base font size (increased to 36pt for much better visibility)
            self._topology_original_fonts = {}  # Store original font sizes for accurate scaling
            
            # Add proper cleanup handler
            def on_window_close():
                try:
                    # Clear the window reference
                    if hasattr(self, '_topology_window'):
                        self._topology_window = None
                    # Clean up font storage
                    self._topology_original_fonts = {}
                    # Destroy the window
                    analysis_window.destroy()
                except:
                    pass
            
            analysis_window.protocol("WM_DELETE_WINDOW", on_window_close)
            
            # Create safe lambda functions that check window existence
            def safe_font_increase(event):
                if analysis_window.winfo_exists():
                    self._change_topology_font_size(analysis_window, 6)
            
            def safe_font_decrease(event):
                if analysis_window.winfo_exists():
                    self._change_topology_font_size(analysis_window, -6)
            
            def safe_font_reset(event):
                if analysis_window.winfo_exists():
                    self._reset_topology_font_size(analysis_window)
            
            # Bind keyboard shortcuts for font size control (multiple variants to ensure they work)
            analysis_window.bind('<Control-plus>', safe_font_increase)
            analysis_window.bind('<Control-KP_Add>', safe_font_increase)  # Numpad +
            analysis_window.bind('<Control-equal>', safe_font_increase)  # = key (often same as +)
            analysis_window.bind('<Control-minus>', safe_font_decrease)
            analysis_window.bind('<Control-KP_Subtract>', safe_font_decrease)  # Numpad -
            analysis_window.bind('<Control-underscore>', safe_font_decrease)  # _ key
            analysis_window.bind('<Control-Key-0>', safe_font_reset)  # Reset to default
            
            # Set focus to window so keyboard shortcuts work
            analysis_window.focus_force()
            
            # Initialize font storage as list of (font_object, original_size) tuples
            # Can't use dict because CTkFont objects are not hashable
            self._topology_font_objects = []
            
            # Create shared font objects for all UI elements
            # These fonts will be reused across all widgets for consistent sizing
            title_font = ctk.CTkFont(size=50, weight="bold")
            shortcuts_font = ctk.CTkFont(size=30, slant="italic")
            summary_label_font = ctk.CTkFont(size=35, weight="bold")
            summary_value_font = ctk.CTkFont(size=68, weight="bold")
            header_font = ctk.CTkFont(size=14, weight="bold")
            text_font = ctk.CTkFont(family="Courier", size=15)
            button_font = ctk.CTkFont(size=13)
            # Additional fonts for tab content
            category_header_font = ctk.CTkFont(size=21, weight="bold")
            category_count_font = ctk.CTkFont(size=17, weight="bold")
            residue_name_font = ctk.CTkFont(size=18, weight="bold")
            residue_count_font = ctk.CTkFont(size=17)
            text_small_font = ctk.CTkFont(family="Courier", size=14)
            text_tiny_font = ctk.CTkFont(family="Courier", size=11)
            
            # Store fonts with their original sizes for proportional scaling
            self._topology_font_objects.append((title_font, 50))
            self._topology_font_objects.append((shortcuts_font, 30))
            self._topology_font_objects.append((summary_label_font, 35))
            self._topology_font_objects.append((summary_value_font, 68))
            self._topology_font_objects.append((header_font, 14))
            self._topology_font_objects.append((text_font, 15))
            self._topology_font_objects.append((button_font, 13))
            self._topology_font_objects.append((category_header_font, 21))
            self._topology_font_objects.append((category_count_font, 17))
            self._topology_font_objects.append((residue_name_font, 18))
            self._topology_font_objects.append((residue_count_font, 17))
            self._topology_font_objects.append((text_small_font, 14))
            self._topology_font_objects.append((text_tiny_font, 11))
            
            # Store for later use in creating widgets
            self._topology_fonts = {
                'title': title_font,
                'shortcuts': shortcuts_font,
                'summary_label': summary_label_font,
                'summary_value': summary_value_font,
                'header': header_font,
                'text': text_font,
                'button': button_font,
                'category_header': category_header_font,
                'category_count': category_count_font,
                'residue_name': residue_name_font,
                'residue_count': residue_count_font,
                'text_small': text_small_font,
                'text_tiny': text_tiny_font
            }
            
            # Title
            title = ctk.CTkLabel(
                analysis_window,
                text=f"Topology Analysis: {os.path.basename(topology_path)}",
                font=self._topology_fonts['title']
            )
            title.pack(pady=(15, 5), padx=10)
            
            # Keyboard shortcuts hint
            shortcuts_hint = ctk.CTkLabel(
                analysis_window,
                text="Zoom: Ctrl+Plus / Ctrl+Minus | Reset: Ctrl+0",
                font=self._topology_fonts['shortcuts'],
                text_color="gray"
            )
            shortcuts_hint.pack(pady=(0, 5))
            
            # Store widget references for font size changes
            self._topology_widgets = {'title': title, 'shortcuts_hint': shortcuts_hint}
            
            # Summary info frame with better styling
            summary_frame = ctk.CTkFrame(analysis_window, fg_color=("gray85", "gray20"))
            summary_frame.pack(fill="x", padx=15, pady=(5, 10))
            
            # Create a grid for better layout - configure column weights to expand evenly
            summary_grid = ctk.CTkFrame(summary_frame, fg_color="transparent")
            summary_grid.pack(fill="x", expand=True, pady=10, padx=10)
            
            # Configure grid columns to have equal weight for even distribution
            summary_grid.grid_columnconfigure(0, weight=1)
            summary_grid.grid_columnconfigure(1, weight=1)
            summary_grid.grid_columnconfigure(2, weight=1)
            
            # Atoms info
            atoms_label = ctk.CTkLabel(
                summary_grid,
                text="ATOMS",
                font=self._topology_fonts['summary_label'],
                text_color="gray"
            )
            atoms_label.grid(row=0, column=0, padx=20, sticky="ew")
            
            atoms_value = ctk.CTkLabel(
                summary_grid,
                text=str(universe.atoms.n_atoms),  # type: ignore[union-attr]
                font=self._topology_fonts['summary_value']
            )
            atoms_value.grid(row=1, column=0, padx=20, sticky="ew")
            
            # Residues info
            residues_label = ctk.CTkLabel(
                summary_grid,
                text="RESIDUES",
                font=self._topology_fonts['summary_label'],
                text_color="gray"
            )
            residues_label.grid(row=0, column=1, padx=20, sticky="ew")
            
            residues_value = ctk.CTkLabel(
                summary_grid,
                text=str(universe.atoms.n_residues),  # type: ignore[union-attr]
                font=self._topology_fonts['summary_value']
            )
            residues_value.grid(row=1, column=1, padx=20, sticky="ew")
            
            # Chains info
            chains_label = ctk.CTkLabel(
                summary_grid,
                text="CHAINS/SEGMENTS",
                font=self._topology_fonts['summary_label'],
                text_color="gray"
            )
            chains_label.grid(row=0, column=2, padx=20, sticky="ew")
            
            chains_value = ctk.CTkLabel(
                summary_grid,
                text=str(universe.atoms.n_segments),  # type: ignore[union-attr]
                font=self._topology_fonts['summary_value']
            )
            chains_value.grid(row=1, column=2, padx=20, sticky="ew")
            
            # Store references for font updates
            self._topology_widgets.update({
                'atoms_label': atoms_label,
                'atoms_value': atoms_value,
                'residues_label': residues_label,
                'residues_value': residues_value,
                'chains_label': chains_label,
                'chains_value': chains_value
            })
            
            # Create tabbed interface for different views (transparent background to remove gray bar)
            tabview = ctk.CTkTabview(analysis_window, width=750, fg_color="transparent")
            tabview.pack(fill="both", expand=True, padx=15, pady=(5, 10))
            
            # IMPORTANT: Hide the segmented button gray bar by matching it to the window background
            # The segmented button (tab buttons) has a default gray background that creates a visible bar.
            # To blend it seamlessly, we set its fg_color to match the CTkToplevel window background.
            # Use ("gray90", "gray13") which is the default CTkToplevel background color (light mode, dark mode)
            try:
                # Access the internal segmented button and make its background match the window
                if hasattr(tabview, '_segmented_button'):
                    tabview._segmented_button.configure(fg_color=("gray90", "gray13"))
            except:
                pass
            
            # Tab 1: Residue Summary
            tab_summary = tabview.add("Residue summary")
            # Tab 2: Atom names
            tab_detailed = tabview.add("Atom names")
            # Tab 3: Topology Selections
            tab_examples = tabview.add("Topology selections")
            
            # Analyze residues by category (before creating tabs content)
            residue_categories = {
                'Protein': [],
                'Nucleic': [],
                'Water': [],
                'Ions': [],
                'Lipids': [],
                'Other': []
            }
            
            # Common ion names
            ion_names = {'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'FE', 'CU', 'MN', 'SOD', 'CLA', 'POT', 'CAL'}
            
            # Common lipid residue names
            lipid_names = {'POPC', 'POPE', 'POPS', 'DPPC', 'DMPC', 'DOPC', 'CHOL', 'CHOLESTEROL',
                          'DLPC', 'DLPE', 'DLPS', 'PALM', 'OLEO', 'STEA', 'PA', 'PC', 'OL', 'PE', 'PS', 'PG',
                          'POPG', 'DOPG', 'DPPG', 'DMPG', 'DLPG'}
            
            # Protein residue names for classification
            protein_residues = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                              'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
                              'HIE', 'HID', 'HIP', 'CYX', 'HSD', 'HSE', 'HSP', 'ACE', 'NME'}
            
            # Nucleic acid residue names
            nucleic_residues = {'A', 'T', 'G', 'C', 'U', 'DA', 'DT', 'DG', 'DC', 'ADE', 'THY', 'GUA',
                              'CYT', 'URA', 'RA', 'RU', 'RG', 'RC'}
            
            # Water residue names
            water_residues = {'WAT', 'HOH', 'TIP3', 'TIP4', 'TIP5', 'SPC', 'SOL', 'H2O'}
            
            for residue in universe.residues:  # type: ignore[union-attr]
                res_info = {
                    'name': residue.resname,
                    'index': residue.resindex,
                    'resSeq': residue.resid,
                    'chain': residue.segid,
                    'n_atoms': len(residue.atoms)
                }
                
                # Classify residue based on residue name since MDAnalysis doesn't have is_protein, etc.
                if residue.resname.upper() in protein_residues:
                    residue_categories['Protein'].append(res_info)
                elif residue.resname.upper() in water_residues:
                    residue_categories['Water'].append(res_info)
                elif residue.resname.upper() in ion_names:
                    residue_categories['Ions'].append(res_info)
                elif residue.resname.upper() in lipid_names or 'LIP' in residue.resname.upper():
                    residue_categories['Lipids'].append(res_info)
                elif residue.resname.upper() in nucleic_residues:
                    residue_categories['Nucleic'].append(res_info)
                else:
                    residue_categories['Other'].append(res_info)
            
            # === TAB 1: Residue Summary (Enhanced) ===
            # Create scrollable frame for better navigation
            summary_scroll = ctk.CTkScrollableFrame(tab_summary, fg_color="transparent")
            summary_scroll.pack(fill="both", expand=True, padx=5, pady=5)
            
            # Build summary content for export
            summary_content = "=" * 80 + "\n"
            summary_content += "TOPOLOGY SUMMARY\n"
            summary_content += "=" * 80 + "\n\n"
            
            row_idx = 0
            
            # Category symbols
            category_icons = {
                'Protein': '[P]',
                'Nucleic': '[N]',
                'Water': '[W]',
                'Ions': '[I]',
                'Lipids': '[L]',
                'Other': '[O]'
            }
            
            # Category colors (for better visual distinction)
            category_colors = {
                'Protein': ("blue", "cyan"),
                'Nucleic': ("purple", "magenta"),
                'Water': ("lightblue", "lightblue"),
                'Ions': ("orange", "orange"),
                'Lipids': ("green", "lightgreen"),
                'Other': ("gray", "gray")
            }
            
            for category, residues in residue_categories.items():
                if residues:
                    summary_content += f"\n{'='*80}\n"
                    summary_content += f"{category.upper()} ({len(residues)} residues)\n"
                    summary_content += f"{'='*80}\n\n"
                    
                    # Category header with icon
                    cat_frame = ctk.CTkFrame(summary_scroll, fg_color=("gray85", "gray20"))
                    cat_frame.grid(row=row_idx, column=0, sticky="ew", padx=5, pady=(10, 5))
                    cat_frame.grid_columnconfigure(1, weight=1)
                    row_idx += 1
                    
                    # Symbol and title
                    icon = category_icons.get(category, '[?]')
                    cat_label = ctk.CTkLabel(
                        cat_frame,
                        text=f"{icon} {category.upper()}",
                        font=self._topology_fonts['category_header'],
                        anchor="w"
                    )
                    cat_label.grid(row=0, column=0, sticky="w", padx=10, pady=8)
                    
                    # Count badge
                    count_label = ctk.CTkLabel(
                        cat_frame,
                        text=f"{len(residues)} residues",
                        font=self._topology_fonts['category_count'],
                        text_color=category_colors.get(category, ("gray", "gray")),
                        anchor="e"
                    )
                    count_label.grid(row=0, column=1, sticky="e", padx=10, pady=8)
                    
                    # Group by residue name
                    from collections import defaultdict
                    by_name = defaultdict(list)
                    for res in residues:
                        by_name[res['name']].append(res)
                    
                    # Content frame for residue details
                    content_frame = ctk.CTkFrame(summary_scroll, fg_color=("gray95", "gray15"))
                    content_frame.grid(row=row_idx, column=0, sticky="ew", padx=5, pady=(0, 5))
                    row_idx += 1
                    
                    inner_row = 0
                    for resname in sorted(by_name.keys()):
                        res_list = by_name[resname]
                        summary_content += f"  {resname} ({len(res_list)} residues):\n"
                        
                        # Residue name and count
                        res_name_frame = ctk.CTkFrame(content_frame, fg_color="transparent")
                        res_name_frame.grid(row=inner_row, column=0, sticky="ew", padx=10, pady=5)
                        res_name_frame.grid_columnconfigure(1, weight=1)
                        
                        res_name_label = ctk.CTkLabel(
                            res_name_frame,
                            text=f"● {resname}",
                            font=self._topology_fonts['residue_name'],
                            anchor="w"
                        )
                        res_name_label.grid(row=0, column=0, sticky="w")
                        
                        res_count_label = ctk.CTkLabel(
                            res_name_frame,
                            text=f"({len(res_list)} residues)",
                            font=self._topology_fonts['residue_count'],
                            text_color="gray",
                            anchor="w"
                        )
                        res_count_label.grid(row=0, column=1, sticky="w", padx=(5, 0))
                        inner_row += 1
                        
                        # Show residue indices and resSeq
                        if len(res_list) <= 20:
                            # Show all if few residues
                            details_text = ""
                            for res in res_list:
                                details_text += f"Index: {res['index']:4d}, ResSeq: {res['resSeq']:4d}, Chain: {res['chain']}, Atoms: {res['n_atoms']}\n"
                                summary_content += f"    Index: {res['index']:4d}, ResSeq: {res['resSeq']:4d}, Chain: {res['chain']}, Atoms: {res['n_atoms']}\n"
                            
                            details_label = ctk.CTkLabel(
                                content_frame,
                                text=details_text.strip(),
                                font=self._topology_fonts['text'],
                                anchor="w",
                                justify="left"
                            )
                            details_label.grid(row=inner_row, column=0, sticky="w", padx=20, pady=(0, 8))
                            inner_row += 1
                        else:
                            # Show range for many residues
                            min_idx = min(r['index'] for r in res_list)
                            max_idx = max(r['index'] for r in res_list)
                            min_seq = min(r['resSeq'] for r in res_list)
                            max_seq = max(r['resSeq'] for r in res_list)
                            
                            range_text = f"Index range: {min_idx}-{max_idx}, ResSeq range: {min_seq}-{max_seq}"
                            range_label = ctk.CTkLabel(
                                content_frame,
                                text=range_text,
                                font=self._topology_fonts['text'],
                                anchor="w"
                            )
                            range_label.grid(row=inner_row, column=0, sticky="w", padx=20, pady=(0, 2))
                            inner_row += 1
                            
                            # First 5
                            first_5 = "First 5: " + ", ".join([f"({r['index']},{r['resSeq']})" for r in res_list[:5]])
                            first_label = ctk.CTkLabel(
                                content_frame,
                                text=first_5,
                                font=self._topology_fonts['text_small'],
                                text_color="gray",
                                anchor="w"
                            )
                            first_label.grid(row=inner_row, column=0, sticky="w", padx=20, pady=1)
                            inner_row += 1
                            
                            # Last 5
                            last_5 = "Last 5:  " + ", ".join([f"({r['index']},{r['resSeq']})" for r in res_list[-5:]])
                            last_label = ctk.CTkLabel(
                                content_frame,
                                text=last_5,
                                font=self._topology_fonts['text_small'],
                                text_color="gray",
                                anchor="w"
                            )
                            last_label.grid(row=inner_row, column=0, sticky="w", padx=20, pady=(1, 8))
                            inner_row += 1
                            
                            summary_content += f"    Index range: {min_idx}-{max_idx}, ResSeq range: {min_seq}-{max_seq}\n"
                            summary_content += f"    First 5: "
                            for res in res_list[:5]:
                                summary_content += f"({res['index']},{res['resSeq']}) "
                            summary_content += "\n    Last 5:  "
                            for res in res_list[-5:]:
                                summary_content += f"({res['index']},{res['resSeq']}) "
                            summary_content += "\n"
                        summary_content += "\n"
            
            # === TAB 2: Atom names View ===
            detailed_textbox = ctk.CTkTextbox(tab_detailed, wrap="word", font=self._topology_fonts['text'])
            detailed_textbox.pack(fill="both", expand=True, padx=5, pady=5)
            
            detailed_content = "=" * 80 + "\n"
            detailed_content += "DETAILED ATOM AND RESIDUE INFORMATION\n"
            detailed_content += "=" * 80 + "\n\n"
            
            # Get unique atom names by residue type
            for category, residues in residue_categories.items():
                if residues:
                    detailed_content += f"\n{'='*80}\n"
                    detailed_content += f"{category.upper()} - ATOM NAMES\n"
                    detailed_content += f"{'='*80}\n\n"
                    
                    # Group by residue name
                    from collections import defaultdict
                    by_name = defaultdict(list)
                    for res in residues:
                        by_name[res['name']].append(res)
                    
                    for resname in sorted(by_name.keys()):
                        # Get first residue of this type to show atom names
                        first_res_idx = by_name[resname][0]['index']
                        residue_obj = universe.residues[first_res_idx]  # type: ignore[index]
                        
                        atom_names = [atom.name for atom in residue_obj.atoms]
                        detailed_content += f"  {resname}:\n"
                        detailed_content += f"    Atoms ({len(atom_names)}): {', '.join(atom_names)}\n\n"
            
            detailed_textbox.insert("0.0", detailed_content)
            detailed_textbox.configure(state="disabled")
            
            # === TAB 3: Topology Selections (Interactive) ===
            # Create scrollable frame for selection examples from the current topology
            examples_scroll = ctk.CTkScrollableFrame(tab_examples, fg_color="transparent")
            examples_scroll.pack(fill="both", expand=True, padx=5, pady=5)
            
            # Helper function to create selection button
            def create_selection_button(parent, description, selection, row):
                # Frame for each selection
                sel_frame = ctk.CTkFrame(parent, fg_color=("gray90", "gray17"))
                sel_frame.grid(row=row, column=0, sticky="ew", padx=5, pady=3)
                sel_frame.grid_columnconfigure(0, weight=1)
                
                # Description
                desc_label = ctk.CTkLabel(
                    sel_frame,
                    text=description,
                    font=self._topology_fonts['button'],
                    anchor="w"
                )
                desc_label.grid(row=0, column=0, sticky="w", padx=10, pady=(8, 2))
                
                # Selection text (monospace)
                sel_label = ctk.CTkLabel(
                    sel_frame,
                    text=selection,
                    font=self._topology_fonts['text_tiny'],
                    anchor="w",
                    text_color=("blue", "cyan")
                )
                sel_label.grid(row=1, column=0, sticky="w", padx=10, pady=(2, 8))
                
                # Copy button
                copy_btn = ctk.CTkButton(
                    sel_frame,
                    text="Copy",
                    width=80,
                    height=28,
                    command=lambda s=selection: self._copy_to_clipboard(s, analysis_window)
                )
                copy_btn.grid(row=0, column=1, rowspan=2, padx=10, pady=5)
                
                return sel_frame
            
            # Build examples content for export
            examples_content = "=" * 80 + "\n"
            examples_content += "SELECTION EXAMPLES FOR YOUR SYSTEM\n"
            examples_content += "=" * 80 + "\n\n"
            
            row_counter = 0
            
            # Generate relevant examples
            if residue_categories['Protein']:
                # Protein section header
                protein_header = ctk.CTkLabel(
                    examples_scroll,
                    text="[P] PROTEIN SELECTIONS",
                    font=self._topology_fonts['header'],
                    anchor="w"
                )
                protein_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(10, 5))
                row_counter += 1
                
                examples_content += "PROTEIN SELECTIONS:\n" + "-" * 40 + "\n"
                
                # Show actual residue index ranges
                protein_indices = [r['index'] for r in residue_categories['Protein']]
                protein_resseqs = [r['resSeq'] for r in residue_categories['Protein']]
                min_idx, max_idx = min(protein_indices), max(protein_indices)
                min_seq, max_seq = min(protein_resseqs), max(protein_resseqs)
                
                # All protein
                create_selection_button(examples_scroll, "All protein atoms", "protein", row_counter)
                row_counter += 1
                examples_content += "• All protein atoms:\n  protein\n\n"
                
                # Backbone
                create_selection_button(examples_scroll, "Protein backbone (CA, C, N, O)", "backbone", row_counter)
                row_counter += 1
                examples_content += "• Protein backbone:\n  backbone\n\n"
                
                # Alpha carbons
                create_selection_button(examples_scroll, "Protein alpha carbons only", "name CA", row_counter)
                row_counter += 1
                examples_content += "• Protein alpha carbons:\n  name CA\n\n"
                
                # By residue index
                sel_resid = f"protein and resid {min_idx} to {min(min_idx+50, max_idx)}"
                create_selection_button(
                    examples_scroll, 
                    f"Protein by residue index (0-based, range {min_idx}-{max_idx})", 
                    sel_resid, 
                    row_counter
                )
                row_counter += 1
                examples_content += f"• Protein by residue index:\n  {sel_resid}\n\n"
                
                # By PDB numbering
                sel_resseq = f"protein and resSeq {min_seq} to {min(min_seq+50, max_seq)}"
                create_selection_button(
                    examples_scroll, 
                    f"Protein by PDB numbering (range {min_seq}-{max_seq})", 
                    sel_resseq, 
                    row_counter
                )
                row_counter += 1
                examples_content += f"• Protein by PDB numbering:\n  {sel_resseq}\n\n"
                
                # Without hydrogens
                create_selection_button(examples_scroll, "Protein without hydrogens", "protein and not name H*", row_counter)
                row_counter += 1
                examples_content += "• Protein without hydrogens:\n  protein and not name H*\n\n"
                
                # Note
                note_label = ctk.CTkLabel(
                    examples_scroll,
                    text="TIP: Use 'resid' for MDAnalysis residue index, 'resnum' for PDB numbering",
                    font=self._topology_fonts['text_tiny'],
                    text_color="orange",
                    anchor="w"
                )
                note_label.grid(row=row_counter, column=0, sticky="w", padx=10, pady=5)
                row_counter += 1
            
            if residue_categories['Lipids']:
                # Lipid section header
                lipid_header = ctk.CTkLabel(
                    examples_scroll,
                    text="[L] LIPID/MEMBRANE SELECTIONS",
                    font=self._topology_fonts['header'],
                    anchor="w"
                )
                lipid_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(15, 5))
                row_counter += 1
                
                examples_content += "\nLIPID/MEMBRANE SELECTIONS:\n" + "-" * 40 + "\n"
                
                lipid_names_found = set(r['name'] for r in residue_categories['Lipids'])
                lipid_indices = [r['index'] for r in residue_categories['Lipids']]
                min_lip_idx, max_lip_idx = min(lipid_indices), max(lipid_indices)
                
                # Info label
                info_label = ctk.CTkLabel(
                    examples_scroll,
                    text=f"Lipid types found: {', '.join(sorted(lipid_names_found))}",
                    font=self._topology_fonts['text_tiny'],
                    text_color="gray",
                    anchor="w"
                )
                info_label.grid(row=row_counter, column=0, sticky="w", padx=10, pady=3)
                row_counter += 1
                
                # All lipids
                sel_all_lipids = f"resname {' '.join(sorted(lipid_names_found))}"
                create_selection_button(examples_scroll, "All lipids", sel_all_lipids, row_counter)
                row_counter += 1
                examples_content += f"• All lipids:\n  {sel_all_lipids}\n\n"
                
                # Lipid headgroups (if POPC exists)
                if 'POPC' in lipid_names_found:
                    create_selection_button(
                        examples_scroll, 
                        "POPC headgroups (phosphate and nitrogen)", 
                        "resname POPC and (name P or name N)", 
                        row_counter
                    )
                    row_counter += 1
                    examples_content += "• POPC headgroups:\n  resname POPC and (name P or name N)\n\n"
                
                # Lipids by index
                sel_lip_idx = f"resid {min_lip_idx} to {max_lip_idx}"
                create_selection_button(
                    examples_scroll, 
                    f"Lipids by index range ({min_lip_idx}-{max_lip_idx})", 
                    sel_lip_idx, 
                    row_counter
                )
                row_counter += 1
                examples_content += f"• Lipids by index:\n  {sel_lip_idx}\n\n"
            
            if residue_categories['Water']:
                # Water section header
                water_header = ctk.CTkLabel(
                    examples_scroll,
                    text="[W] WATER SELECTIONS",
                    font=self._topology_fonts['header'],
                    anchor="w"
                )
                water_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(15, 5))
                row_counter += 1
                
                examples_content += "\nWATER SELECTIONS:\n" + "-" * 40 + "\n"
                
                # All water
                create_selection_button(examples_scroll, "All water molecules", "water", row_counter)
                row_counter += 1
                examples_content += "• All water:\n  water\n\n"
                
                # Water oxygens
                create_selection_button(examples_scroll, "Water oxygen atoms only", "water and name O", row_counter)
                row_counter += 1
                examples_content += "• Water oxygens:\n  water and name O\n\n"
            
            if residue_categories['Ions']:
                # Ion section header
                ion_header = ctk.CTkLabel(
                    examples_scroll,
                    text="[I] ION SELECTIONS",
                    font=self._topology_fonts['header'],
                    anchor="w"
                )
                ion_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(15, 5))
                row_counter += 1
                
                examples_content += "\nION SELECTIONS:\n" + "-" * 40 + "\n"
                
                ion_names_found = set(r['name'] for r in residue_categories['Ions'])
                
                # Info label
                ion_info_label = ctk.CTkLabel(
                    examples_scroll,
                    text=f"Ion types found: {', '.join(sorted(ion_names_found))}",
                    font=self._topology_fonts['text_tiny'],
                    text_color="gray",
                    anchor="w"
                )
                ion_info_label.grid(row=row_counter, column=0, sticky="w", padx=10, pady=3)
                row_counter += 1
                
                # All ions
                sel_all_ions = f"resname {' '.join(sorted(ion_names_found))}"
                create_selection_button(examples_scroll, "All ions", sel_all_ions, row_counter)
                row_counter += 1
                examples_content += f"• All ions:\n  {sel_all_ions}\n\n"
            
            if residue_categories['Other']:
                # Other section header
                other_header = ctk.CTkLabel(
                    examples_scroll,
                    text="[O] OTHER MOLECULES (Ligands, cofactors, etc.)",
                    font=self._topology_fonts['header'],
                    anchor="w"
                )
                other_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(15, 5))
                row_counter += 1
                
                examples_content += "\nOTHER MOLECULES:\n" + "-" * 40 + "\n"
                
                other_names_found = set(r['name'] for r in residue_categories['Other'])
                
                # Info label
                other_info_label = ctk.CTkLabel(
                    examples_scroll,
                    text=f"Molecule types found: {', '.join(sorted(other_names_found))}",
                    font=self._topology_fonts['text_tiny'],
                    text_color="gray",
                    anchor="w"
                )
                other_info_label.grid(row=row_counter, column=0, sticky="w", padx=10, pady=3)
                row_counter += 1
                
                for resname in sorted(other_names_found):
                    create_selection_button(examples_scroll, resname, f"resname {resname}", row_counter)
                    row_counter += 1
                    examples_content += f"• {resname}:\n  resname {resname}\n"
            
            # Combined selections section
            combined_header = ctk.CTkLabel(
                examples_scroll,
                text="[+] COMBINED SELECTIONS",
                font=self._topology_fonts['header'],
                anchor="w"
            )
            combined_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(15, 5))
            row_counter += 1
            
            examples_content += "\nCOMBINED SELECTIONS:\n" + "-" * 40 + "\n"
            
            # Protein and lipids
            if residue_categories['Protein'] and residue_categories['Lipids']:
                lipid_names_found = set(r['name'] for r in residue_categories['Lipids'])
                sel_prot_lip = f"protein or resname {' '.join(sorted(lipid_names_found))}"
                create_selection_button(examples_scroll, "Protein and lipids together", sel_prot_lip, row_counter)
                row_counter += 1
                examples_content += f"• Protein and lipids:\n  {sel_prot_lip}\n\n"
            
            # Everything except water
            create_selection_button(examples_scroll, "Everything except water", "not water", row_counter)
            row_counter += 1
            examples_content += "• Everything except water:\n  not water\n\n"
            
            # Important notes section
            notes_header = ctk.CTkLabel(
                examples_scroll,
                text="[!] IMPORTANT NOTES",
                font=self._topology_fonts['header'],
                anchor="w"
            )
            notes_header.grid(row=row_counter, column=0, sticky="w", padx=5, pady=(15, 5))
            row_counter += 1
            
            notes_frame = ctk.CTkFrame(examples_scroll, fg_color=("gray95", "gray15"))
            notes_frame.grid(row=row_counter, column=0, sticky="ew", padx=5, pady=5)
            
            notes_text = """* 'resid' uses MDAnalysis residue indices (0-based within each segment)
* 'resnum' uses PDB residue sequence numbers (may have gaps)
* Use 'resnum' if you want to match PDB numbering
* Use 'resid' for programmatic access
* Segment IDs are used for chain identification
* Click the Copy button to copy any selection to clipboard"""
            
            notes_label = ctk.CTkLabel(
                notes_frame,
                text=notes_text,
                font=self._topology_fonts['text_tiny'],
                anchor="w",
                justify="left"
            )
            notes_label.pack(padx=10, pady=10, anchor="w")
            
            examples_content += "\nIMPORTANT NOTES:\n" + "-" * 40 + "\n"
            examples_content += notes_text
            
            # Enable mouse wheel scrolling for all scrollable frames
            self._enable_scrolling_for_panel(summary_scroll)
            self._enable_scrolling_for_panel(examples_scroll)
            
            # Export button with better styling
            button_frame = ctk.CTkFrame(analysis_window, fg_color="transparent")
            button_frame.pack(pady=(5, 15))
            
            export_btn = ctk.CTkButton(
                button_frame,
                text="Export Analysis to File",
                command=lambda: self._export_topology_analysis(summary_content, detailed_content, examples_content),
                width=200,
                height=35,
                font=self._topology_fonts['button']
            )
            export_btn.pack()
            
        except Exception as e:
            messagebox.showerror("Analysis Error", f"Failed to analyze topology:\n{str(e)}")
            logger.error(f"Topology analysis error: {e}", exc_info=True)
    
    def _store_original_fonts(self, window):
        """Store original font sizes for accurate scaling."""
        def store_fonts(widget):
            try:
                widget_id = id(widget)
                widget_type = type(widget).__name__
                
                if any(t in widget_type for t in ['CTkLabel', 'CTkButton', 'CTkTextbox']):
                    try:
                        # Try _font attribute first
                        font_attr = getattr(widget, '_font', None)
                        if font_attr and hasattr(font_attr, 'cget'):
                            self._topology_original_fonts[widget_id] = font_attr.cget('size')
                        else:
                            # Fallback to configure method
                            current_font = widget.cget('font')
                            if current_font and isinstance(current_font, ctk.CTkFont):
                                self._topology_original_fonts[widget_id] = current_font.cget('size')
                    except:
                        pass
                
                # Process children
                if hasattr(widget, 'winfo_children'):
                    for child in widget.winfo_children():
                        store_fonts(child)
            except:
                pass
        
        store_fonts(window)
    
    def _change_topology_font_size(self, window, delta):
        """Change font size of all widgets in topology window using direct font object updates."""
        try:
            # Check if window still exists
            if not window.winfo_exists():
                return
            
            # Update the base font size
            self._topology_font_size += delta
            
            # Clamp font size to reasonable range
            self._topology_font_size = max(20, min(60, self._topology_font_size))
            
            # Update all stored font objects directly
            # This is more reliable than trying to recursively find and update widgets
            if hasattr(self, '_topology_font_objects') and self._topology_font_objects:
                for font_obj, original_size in self._topology_font_objects:
                    try:
                        # Calculate new size maintaining the ratio to original
                        ratio = original_size / 36  # 36 is our base size
                        new_size = int(self._topology_font_size * ratio)
                        # Clamp to reasonable bounds, but proportionally
                        min_size = max(6, int(original_size * 0.5))  # At least 50% of original or 6pt
                        max_size = int(original_size * 3)  # At most 3x original
                        new_size = max(min_size, min(max_size, new_size))
                        font_obj.configure(size=new_size)
                    except:
                        pass
                
                self.status_callback(f"Font size: {self._topology_font_size}pt")
            else:
                self.status_callback(f"Font size: {self._topology_font_size}pt (no fonts to update)")
            
            # Force window refresh
            window.update_idletasks()
            
        except Exception as e:
            logger.error(f"Error changing topology font size: {e}", exc_info=True)
    
    def _reset_topology_font_size(self, window):
        """Reset topology window font size to default."""
        try:
            # Check if window still exists
            if not window.winfo_exists():
                return
                
            # Calculate the delta needed to get back to 36 (new default)
            delta = 36 - self._topology_font_size
            if delta != 0:
                self._change_topology_font_size(window, delta)
            self.status_callback("Font size reset to default")
        except Exception as e:
            logger.error(f"Error resetting topology font size: {e}", exc_info=True)
    
    def _copy_to_clipboard(self, text, window):
        """Copy text to clipboard and show confirmation."""
        try:
            window.clipboard_clear()
            window.clipboard_append(text)
            window.update()  # Required for clipboard to work
            
            # Show brief confirmation (you could also use a toast notification)
            self.status_callback(f"Copied to clipboard: {text[:50]}{'...' if len(text) > 50 else ''}")
        except Exception as e:
            logger.error(f"Failed to copy to clipboard: {e}")
            messagebox.showwarning("Copy Failed", "Could not copy to clipboard")
    
    def _export_topology_analysis(self, summary, detailed, examples):
        """Export topology analysis to a text file."""
        filename = filedialog.asksaveasfilename(
            title="Export Topology Analysis",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("All files", "*.*")],
            initialdir=self.output_directory.get() if self.output_directory.get() else os.getcwd()
        )
        
        if filename:
            try:
                with open(filename, 'w') as f:
                    f.write(summary)
                    f.write("\n\n")
                    f.write(detailed)
                    f.write("\n\n")
                    f.write(examples)
                messagebox.showinfo("Success", f"Topology analysis exported to:\n{filename}")
            except Exception as e:
                messagebox.showerror("Export Error", f"Failed to export analysis:\n{str(e)}")
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the Analysis frame to support Ctrl+/Ctrl- shortcuts."""
        try:
            # Update section title labels (headings)
            title_labels = [
                'file_title_label', 'settings_title_label', 'plot_title_label',
                'results_title_label', 'control_title_label', 'namd_file_title_label'
            ]
            for widget_name in title_labels:
                if hasattr(self, widget_name):
                    try:
                        getattr(self, widget_name).configure(font=scaled_fonts['heading'])
                    except Exception:
                        pass
            
            # Update all CTkLabel widgets that are attributes of self (body text)
            label_widgets = [
                'time_units_label', 'rmsf_xaxis_label', 'rmsf_format_label', 
                'rmsf_frequency_label', 'progress_label', 'custom_selection_label',
                'analysis_type_label', 'plot_title_label'
            ]
            for widget_name in label_widgets:
                if hasattr(self, widget_name):
                    try:
                        getattr(self, widget_name).configure(font=scaled_fonts['body'])
                    except Exception:
                        pass
            
            # Recursively update all CTkLabel widgets in the left panel that aren't stored as attributes
            def update_child_labels(widget):
                if widget is None:
                    return
                try:
                    for child in widget.winfo_children():
                        if isinstance(child, ctk.CTkLabel) and child not in [
                            getattr(self, attr, None) for attr in dir(self) 
                            if isinstance(getattr(self, attr, None), ctk.CTkLabel)
                        ]:
                            try:
                                child.configure(font=scaled_fonts['body'])
                            except Exception:
                                pass
                        update_child_labels(child)
                except Exception:
                    pass
            
            # Apply recursive update to the main frames
            if hasattr(self, 'left_panel'):
                update_child_labels(self.left_panel)
            
            # Update combo boxes
            combo_widgets = [
                'analysis_type_combo', 'plot_time_units_combo', 'rmsf_xaxis_combo',
                'rmsf_format_combo', 'rmsf_frequency_combo', 'plot_yaxis_units_combo',
                'namd_column_combo'
            ]
            for widget_name in combo_widgets:
                if hasattr(self, widget_name):
                    try:
                        getattr(self, widget_name).configure(font=scaled_fonts['body'])
                    except Exception:
                        pass
            
            # Update checkboxes
            checkbox_widgets = ['plot_checkbox', 'rmsf_labels_checkbox']
            for widget_name in checkbox_widgets:
                if hasattr(self, widget_name):
                    try:
                        getattr(self, widget_name).configure(font=scaled_fonts['body'])
                    except Exception:
                        pass
            
            # Update buttons
            button_widgets = ['run_button', 'namd_log_button', 'status_button']
            for widget_name in button_widgets:
                if hasattr(self, widget_name):
                    try:
                        getattr(self, widget_name).configure(font=scaled_fonts['body'])
                    except Exception:
                        pass
            
            # Update entry fields and other input widgets
            input_widgets = [
                'custom_line_color_entry', 'custom_bg_color_entry', 
                'custom_fig_bg_color_entry', 'custom_text_color_entry',
                'custom_selection_entry', 'topology_file_entry', 'output_directory_entry'
            ]
            for widget_name in input_widgets:
                if hasattr(self, widget_name):
                    try:
                        getattr(self, widget_name).configure(font=scaled_fonts['body'])
                    except Exception:
                        pass
            
            # Recursively update all input widgets in the left panel
            def update_child_inputs(widget):
                if widget is None:
                    return
                try:
                    for child in widget.winfo_children():
                        if isinstance(child, (ctk.CTkEntry, ctk.CTkComboBox, ctk.CTkCheckBox, ctk.CTkButton)):
                            try:
                                child.configure(font=scaled_fonts['body'])
                            except Exception:
                                pass
                        update_child_inputs(child)
                except Exception:
                    pass
            
            # Apply recursive update to input widgets
            if hasattr(self, 'left_panel'):
                update_child_inputs(self.left_panel)
            
            # Update textboxes (results_info, etc.)
            if hasattr(self, 'results_info'):
                try:
                    self.results_info.configure(font=scaled_fonts['small'])
                except Exception:
                    pass
            
            # Update export buttons if they exist
            export_buttons = ['export_csv_btn', 'export_json_btn', 'export_numpy_btn']
            for btn_name in export_buttons:
                if hasattr(self, btn_name):
                    try:
                        getattr(self, btn_name).configure(font=scaled_fonts['small'])
                    except Exception:
                        pass
            
            logger.debug("Fonts updated in AnalysisFrame")
            
        except Exception as e:
            logger.warning(f"Error updating fonts in AnalysisFrame: {e}")
    
    def _on_status_button_click(self):
        """Handle status button click - show current analysis status."""
        try:
            if hasattr(self, 'analysis_changed') and self.analysis_changed:
                messagebox.showinfo(
                    "Analysis Status", 
                    "Settings have changed since last analysis.\nPlease run the analysis again to see updated results."
                )
            else:
                messagebox.showinfo(
                    "Analysis Status", 
                    "Analysis has already been run with the current settings."
                )
        except Exception as e:
            logger.error(f"Error in status button click: {e}")
    
    def _mark_analysis_changed(self):
        """Mark that analysis parameters have changed."""
        if hasattr(self, 'status_button'):
            try:
                # Check if the widget still exists (hasn't been destroyed)
                self.status_button.winfo_exists()
                self.analysis_changed = True
                self.status_button.configure(
                    text="Parameters Changed",
                    fg_color="#F5A300",
                    hover_color="#B87A00",
                    text_color="white"  # Still white text on yellow works fine
                )
            except (AttributeError, tk.TclError):
                # Widget doesn't exist or has been destroyed, ignore
                pass
    
    def _mark_analysis_ready(self):
        """Mark that analysis is ready/up-to-date."""
        if hasattr(self, 'status_button'):
            try:
                # Check if the widget still exists (hasn't been destroyed)
                self.status_button.winfo_exists()
                self.analysis_changed = False
                self.status_button.configure(
                    text="Ready",
                    fg_color="#27A85D",
                    hover_color="#1B7440"
                )
            except (AttributeError, tk.TclError):
                # Widget doesn't exist or has been destroyed, ignore
                pass
    
    def _detect_namd_columns_from_file(self, log_file_path):
        """Detect available columns from a NAMD log file by reading the ETITLE line.
        
        Args:
            log_file_path: Path to the NAMD log file
        """
        try:
            with open(log_file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('ETITLE:'):
                        headers = line.split()[1:]  # Skip 'ETITLE:'
                        
                        # Map common NAMD column names to display names with units
                        column_mapping = {
                            'TEMP': 'Temperature (K)',
                            'TEMPERATURE': 'Temperature (K)',
                            'POTENTIAL': 'Potential Energy (kcal/mol)',
                            'POTENG': 'Potential Energy (kcal/mol)',
                            'KINETIC': 'Kinetic Energy (kcal/mol)',
                            'KINENG': 'Kinetic Energy (kcal/mol)',
                            'TOTAL': 'Total Energy (kcal/mol)',
                            'TOTAL3': 'Total3 Energy (kcal/mol)',
                            'TOTENG': 'Total Energy (kcal/mol)',
                            'PRESSURE': 'Pressure (atm)',
                            'VOLUME': 'Volume (Å³)',
                            'ELECT': 'Electrostatic Energy (kcal/mol)',
                            'VDW': 'Van der Waals Energy (kcal/mol)',
                            'BOND': 'Bond Energy (kcal/mol)',
                            'ANGLE': 'Angle Energy (kcal/mol)',
                            'DIHED': 'Dihedral Energy (kcal/mol)',
                            'IMPR': 'Improper Energy (kcal/mol)',
                            'MISC': 'Miscellaneous Energy (kcal/mol)',
                            'BOUNDARY': 'Boundary Energy (kcal/mol)',
                            'GPRESSURE': 'Gas Pressure (atm)',
                            'PRESSAVG': 'Average Pressure (atm)',
                            'GPRESSAVG': 'Average Gas Pressure (atm)'
                        }
                        
                        # Store detected columns (excluding TS and timestep)
                        self.namd_available_columns = []
                        for header in headers:
                            if header not in ['TS', 'timestep']:
                                display_name = column_mapping.get(header.upper(), f"{header.title()}")
                                self.namd_available_columns.append(display_name)
                        
                        logger.info(f"Detected {len(self.namd_available_columns)} NAMD columns: {self.namd_available_columns}")
                        return
        except Exception as e:
            logger.error(f"Error detecting NAMD columns from {log_file_path}: {e}")
            self.namd_available_columns = []
    
    def _refresh_namd_options(self):
        """Refresh NAMD options UI with newly detected columns."""
        # Clear existing dynamic options
        self._clear_dynamic_options()
        # Recreate NAMD options with new columns
        self._show_namd_log_options()
    
    def _show_namd_log_options(self):
        """Show options for NAMD log analysis."""
        # Auto-load log files if they haven't been loaded yet
        if self.trajectory_files and not self.analysis_results:
            self._auto_load_log_files()
        
        row = 0
        
        # Column selection for NAMD logs
        column_label = ctk.CTkLabel(self.dynamic_options_frame, text="Data to Plot:")
        column_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
        self.dynamic_widgets.append(column_label)
        
        # Get available columns from parsed data
        namd_columns = self._get_available_namd_columns()
        
        column_combo = ctk.CTkComboBox(
            self.dynamic_options_frame,
            values=namd_columns,
            command=self._on_namd_column_change,
            width=250
        )
        column_combo.grid(row=row, column=1, pady=5, padx=5, sticky="w")
        column_combo.set(namd_columns[0] if namd_columns else "No data available")
        self.namd_column_combo = column_combo
        self.dynamic_widgets.append(column_combo)
        
        row += 1
        
        # Multi-column selection
        multi_label = ctk.CTkLabel(self.dynamic_options_frame, text="Add Multiple Columns:")
        multi_label.grid(row=row, column=0, pady=5, padx=(10, 5), sticky="w")
        self.dynamic_widgets.append(multi_label)
        
        # Checkboxes frame for multiple selection
        self.multi_column_frame = ctk.CTkFrame(self.dynamic_options_frame)
        self.multi_column_frame.grid(row=row, column=1, pady=5, padx=5, sticky="w")
        self.dynamic_widgets.append(self.multi_column_frame)
        
        # Create checkboxes for each available column
        self.column_checkboxes = {}
        self._create_column_checkboxes(namd_columns)
        
        # Trigger initial plot if we have data and a valid column selected
        if self.analysis_results and namd_columns and namd_columns[0] != "No data available":
            # Manually trigger the column change callback to initialize the plot
            self._on_namd_column_change(namd_columns[0])
    
    def _get_available_namd_columns(self):
        """Get list of available NAMD columns from detected columns or parsed data."""
        # First, check if we have detected columns from file
        if hasattr(self, 'namd_available_columns') and self.namd_available_columns:
            return self.namd_available_columns
        
        # Fall back to parsed data if available
        if hasattr(self, 'analysis_results') and self.analysis_results:
            # Map common NAMD column names to display names with units
            column_mapping = {
                'TEMP': 'Temperature (K)',
                'TEMPERATURE': 'Temperature (K)',
                'POTENTIAL': 'Potential Energy (kcal/mol)',
                'POTENG': 'Potential Energy (kcal/mol)',
                'KINETIC': 'Kinetic Energy (kcal/mol)',
                'KINENG': 'Kinetic Energy (kcal/mol)',
                'TOTAL': 'Total Energy (kcal/mol)',
                'TOTENG': 'Total Energy (kcal/mol)',
                'PRESSURE': 'Pressure (atm)',
                'VOLUME': 'Volume (Å³)',
                'ELECT': 'Electrostatic Energy (kcal/mol)',
                'VDW': 'Van der Waals Energy (kcal/mol)',
                'BOND': 'Bond Energy (kcal/mol)',
                'ANGLE': 'Angle Energy (kcal/mol)',
                'DIHED': 'Dihedral Energy (kcal/mol)',
                'IMPR': 'Improper Energy (kcal/mol)',
                'MISC': 'Miscellaneous Energy (kcal/mol)',
                'BOUNDARY': 'Boundary Energy (kcal/mol)',
                'GPRESSURE': 'Gas Pressure (atm)',
                'PRESSAVG': 'Average Pressure (atm)',
                'GPRESSAVG': 'Average Gas Pressure (atm)'
            }
            
            available_columns = []
            for key in self.analysis_results.keys():
                # Fix the bug: check if key exists and is not empty before accessing
                # Exclude both 'timestep' and 'TS' as they are X-axis variables, not Y-axis
                if key not in ['timestep', 'TS'] and key in self.analysis_results and self.analysis_results.get(key) is not None:
                    # Also check if it's a list/array with content
                    value = self.analysis_results[key]
                    if (isinstance(value, (list, np.ndarray)) and len(value) > 0) or (not isinstance(value, (list, np.ndarray))):
                        display_name = column_mapping.get(key.upper(), f"{key.title()} (unit)")
                        available_columns.append(display_name)
            
            return available_columns if available_columns else ["No data available"]
        
        # Default columns if no data is loaded yet
        return [
            "Temperature (K)",
            "Potential Energy (kcal/mol)",
            "Kinetic Energy (kcal/mol)", 
            "Total Energy (kcal/mol)",
            "Pressure (atm)",
            "Volume (Å³)"
        ]
    
    def _create_column_checkboxes(self, columns):
        """Create checkboxes for multiple column selection."""
        # Clear existing checkboxes
        for widget in self.multi_column_frame.winfo_children():
            widget.destroy()
        self.column_checkboxes.clear()
        
        # Create checkboxes in a grid layout (2 columns)
        # Show all available columns, not just the first 6
        for i, column in enumerate(columns):
            if column != "No data available":
                var = tk.BooleanVar()
                checkbox = ctk.CTkCheckBox(
                    self.multi_column_frame,
                    text=column.split('(')[0].strip(),  # Short name without units
                    variable=var,
                    command=self._on_multi_column_change,
                    width=120
                )
                checkbox.grid(row=i//2, column=i%2, pady=2, padx=5, sticky="w")
                self.column_checkboxes[column] = var
    
    def _on_multi_column_change(self):
        """Handle multiple column selection changes."""
        # Detect property type from selected columns
        selected_columns = self._get_selected_columns()
        if selected_columns and selected_columns[0] != "No data available":
            # Use the first selected column to determine property type
            property_type = self._detect_property_type_from_column(selected_columns[0])
            self._update_yaxis_units_for_property(property_type)
        
        if self.plot_enabled.get():
            # Ensure canvas exists
            if not hasattr(self, 'canvas') or not self.canvas or not hasattr(self, 'ax') or not self.ax:
                self.create_plot_display()
            
            try:
                self._plot_namd_data()
            except Exception as e:
                logger.error(f"Error updating NAMD plot from multi-column change: {e}", exc_info=True)
    
    def _get_selected_columns(self):
        """Get list of selected columns for multi-plotting."""
        selected = []
        # Always include the main selected column
        if hasattr(self, 'namd_column_combo'):
            main_column = self.namd_column_combo.get()
            if main_column != "No data available":
                selected.append(main_column)
        
        # Add checked columns (only if checkboxes exist)
        if hasattr(self, 'column_checkboxes'):
            for column, var in self.column_checkboxes.items():
                if var.get() and column not in selected:
                    selected.append(column)
        
        return selected
    
    def _show_log_analysis_ui(self):
        """Show UI elements specific to log file analysis."""
        # Hide trajectory section, show log file selection
        if hasattr(self, 'trajectory_listbox'):
            # Hide trajectory-related widgets
            self.trajectory_listbox.grid_remove()
            for widget in self.trajectory_listbox.master.winfo_children():
                if isinstance(widget, ctk.CTkButton) and (widget.cget("text") in ["Add", "Remove"]):
                    widget.grid_remove()
        
        # Show NAMD log file selection button (create if it doesn't exist)
        if not hasattr(self, 'namd_log_button'):
            # Find the parent frame
            file_frame = None
            for child in self.left_panel.winfo_children():
                if isinstance(child, ctk.CTkFrame):
                    for subchild in child.winfo_children():
                        if isinstance(subchild, ctk.CTkLabel) and subchild.cget("text") == "Input Files":
                            file_frame = child
                            break
                    if file_frame:
                        break
            
            if file_frame:
                # Add NAMD log file selection
                ctk.CTkLabel(file_frame, text="NAMD Log:").grid(
                    row=2, column=0, pady=5, padx=(10, 5), sticky="w"
                )
                
                self.namd_log_button = ctk.CTkButton(
                    file_frame,
                    text="Select NAMD Log File",
                    command=self.select_namd_log_file,
                    width=200,
                    height=25
                )
                self.namd_log_button.grid(row=2, column=1, columnspan=2, pady=5, padx=5, sticky="ew")
        else:
            self.namd_log_button.grid()
    
    def _show_trajectory_analysis_ui(self):
        """Show UI elements specific to trajectory analysis."""
        # Show trajectory section, hide log file selection
        if hasattr(self, 'trajectory_listbox'):
            self.trajectory_listbox.grid()
            for widget in self.trajectory_listbox.master.winfo_children():
                if isinstance(widget, ctk.CTkButton) and (widget.cget("text") in ["Add", "Remove"]):
                    widget.grid()
        
        # Hide NAMD log button
        if hasattr(self, 'namd_log_button'):
            self.namd_log_button.grid_remove()
    
    def _update_plot_labels_for_analysis_type(self, analysis_type):
        """Update default plot labels based on analysis type."""
        # Get current units to use in labels
        time_unit = self.plot_time_units.get()
        distance_unit = self.plot_units.get()
        
        label_mapping = {
            "-": {"title": "Analysis Results", "xlabel": "Time", "ylabel": "Value"},
            "RMSD": {"title": "RMSD Analysis", "xlabel": f"Time ({time_unit})", "ylabel": f"RMSD ({distance_unit})"},
            "RMSF": {"title": "RMSF Analysis", "xlabel": "Residue", "ylabel": f"RMSF ({distance_unit})"},
            "Distances": {"title": "Distance Analysis", "xlabel": f"Time ({time_unit})", "ylabel": f"Distance ({distance_unit})"},
            "Radius of Gyration": {"title": "Radius of Gyration", "xlabel": f"Time ({time_unit})", "ylabel": f"Rg ({distance_unit})"},
            "NAMD Log Analysis": {"title": "NAMD Log Analysis", "xlabel": f"Time ({time_unit})", "ylabel": "Value"}
        }
        
        if analysis_type in label_mapping:
            labels = label_mapping[analysis_type]
            self.plot_title.set(labels["title"])
            self.plot_xlabel.set(labels["xlabel"])
            self.plot_ylabel.set(labels["ylabel"])
            
        # Update ylabel for NAMD log analysis based on selected column
        if analysis_type == "NAMD Log Analysis" and hasattr(self, 'namd_log_columns'):
            column = self.namd_log_columns.get()
            if column:
                # Extract unit from column name (e.g., "Temperature (K)" -> "K")
                if "(" in column and ")" in column:
                    unit = column.split("(")[1].split(")")[0]
                    parameter = column.split("(")[0].strip()
                    self.plot_ylabel.set(f"{parameter} ({unit})")
                    self.plot_title.set(f"NAMD {parameter} Analysis")

    def _on_namd_column_change(self, value):
        """Handle NAMD column selection change."""
        # Detect property type and update units
        if value:
            property_type = self._detect_property_type_from_column(value)
            self._update_yaxis_units_for_property(property_type)
        
        # Update plot labels based on selected column
        if value and "(" in value and ")" in value:
            unit = value.split("(")[1].split(")")[0]
            parameter = value.split("(")[0].strip()
            # Use the current Y-axis unit instead of the default from column name
            current_unit = self.plot_yaxis_units.get()
            self.plot_ylabel.set(f"{parameter} ({current_unit})")
            self.plot_title.set(f"NAMD {parameter} Analysis")
        
        # If we have data and plotting is enabled, update the plot
        if self.analysis_results and self.plot_enabled.get():
            # Ensure canvas exists
            if not hasattr(self, 'canvas') or not self.canvas or not hasattr(self, 'ax') or not self.ax:
                self.create_plot_display()
            
            try:
                self._plot_namd_data()
            except Exception as e:
                logger.error(f"Error updating NAMD plot: {e}", exc_info=True)

    def parse_namd_log(self, log_file_path):
        """Parse NAMD log file for energy, temperature, and other data."""
        try:
            data = {'timestep': []}
            column_mapping = {}  # Map column names to indices
            
            with open(log_file_path, 'r') as f:
                lines = f.readlines()
            
            # First pass: find ETITLE line to get column mapping
            etitle_found = False
            for line in lines:
                line = line.strip()
                if line.startswith('ETITLE:'):
                    headers = line.split()[1:]  # Skip 'ETITLE:'
                    # Create mapping of column names to indices
                    for i, header in enumerate(headers):
                        column_mapping[header] = i
                        # Initialize data arrays for all columns found
                        if header not in data:
                            data[header] = []
                    etitle_found = True
                    break
            
            if not etitle_found:
                print("Warning: No ETITLE line found in NAMD log")
                return None
            
            # Second pass: parse ENERGY lines using the column mapping
            timesteps_seen = set()  # Track duplicate timesteps
            previous_timestep = None
            duplicates_found = 0
            
            for line in lines:
                line = line.strip()
                if line.startswith('ENERGY:'):
                    parts = line.split()[1:]  # Skip 'ENERGY:'
                    
                    if len(parts) >= len(column_mapping):
                        try:
                            # Check for timestep first (required)
                            if 'TS' not in column_mapping:
                                continue
                                
                            timestep = int(float(parts[column_mapping['TS']]))
                            
                            # Check for duplicate timesteps
                            if timestep in timesteps_seen:
                                duplicates_found += 1
                                print(f"Warning: Duplicate timestep {timestep} found, skipping (duplicate #{duplicates_found})")
                                continue
                            
                            # Check for non-monotonic timesteps
                            if previous_timestep is not None and timestep < previous_timestep:
                                print(f"Warning: Non-monotonic timestep {timestep} after {previous_timestep}, keeping but noting")
                            
                            timesteps_seen.add(timestep)
                            data['timestep'].append(timestep)
                            previous_timestep = timestep
                            
                            # Extract all available columns
                            for header, index in column_mapping.items():
                                if header != 'TS' and index < len(parts):  # Skip TS as it's already processed
                                    try:
                                        value = float(parts[index])
                                        # Check for unrealistic values that might indicate parsing errors
                                        if abs(value) > 1e10:  # Very large values might be parsing errors
                                            print(f"Warning: Unusually large value {value} for {header} at timestep {timestep}")
                                        data[header].append(value)
                                    except (ValueError, IndexError):
                                        # If we can't parse a value, append None or skip
                                        data[header].append(None)
                        
                        except (ValueError, IndexError) as e:
                            continue
            
            # Post-processing: Check for and report data quality issues
            self._validate_namd_data_quality(data, duplicates_found)
            
            # Ensure all arrays have the same length
            timestep_count = len(data['timestep'])
            for key in data:
                if key != 'timestep':
                    current_length = len(data[key])
                    if current_length < timestep_count:
                        # Pad with None values if needed
                        data[key].extend([None] * (timestep_count - current_length))
                    elif current_length > timestep_count:
                        # Truncate if too long
                        data[key] = data[key][:timestep_count]
            
            # Debug output
            if data['timestep']:
                print(f"Parsed NAMD log: {len(data['timestep'])} data points")
                print(f"Available columns: {list(column_mapping.keys())}")
                if duplicates_found > 0:
                    print(f"Removed {duplicates_found} duplicate timesteps")
                print(f"Timestep range: {min(data['timestep'])} to {max(data['timestep'])}")
            
            return data
            
        except Exception as e:
            print(f"Error parsing NAMD log: {e}")
            return None
    
    def _validate_namd_data_quality(self, data, duplicates_found):
        """Validate NAMD data quality and report issues."""
        if not data or 'timestep' not in data or not data['timestep']:
            return
        
        timesteps = data['timestep']
        
        # Check for gaps in timesteps
        if len(timesteps) > 1:
            timestep_diffs = [timesteps[i+1] - timesteps[i] for i in range(len(timesteps)-1)]
            most_common_diff = max(set(timestep_diffs), key=timestep_diffs.count)
            
            large_gaps = [i for i, diff in enumerate(timestep_diffs) if diff > most_common_diff * 2]
            if large_gaps:
                print(f"Warning: Found {len(large_gaps)} large timestep gaps at indices: {large_gaps[:5]}")
        
        # Check for data completeness
        for key, values in data.items():
            if key != 'timestep':
                none_count = sum(1 for v in values if v is None)
                if none_count > 0:
                    percentage = (none_count / len(values)) * 100
                    print(f"Warning: {key} has {none_count} missing values ({percentage:.1f}%)")
        
        if duplicates_found > 0:
            print(f"Data quality summary: {duplicates_found} duplicate timesteps removed")
    
    def select_namd_log_file(self):
        """Select NAMD log file for analysis."""
        filetypes = [
            ("Log files", "*.log"),
            ("Output files", "*.out"),
            ("Text files", "*.txt"),
            ("All files", "*.*")
        ]
        
        filename = filedialog.askopenfilename(
            title="Select NAMD Log File",
            filetypes=filetypes,
            initialdir=self.initial_directory
        )
        
        if filename:
            # Parse the log file and update analysis results
            data = self.parse_namd_log(filename)
            if data and data['timestep']:
                self.analysis_results = data
                self.status_callback(f"NAMD log loaded: {Path(filename).name} ({len(data['timestep'])} data points)")
                
                # Refresh the UI options with new columns
                if self.analysis_type.get() == "NAMD Log Analysis":
                    self._clear_dynamic_options()
                    self._show_namd_log_options()
                
                # Update plot if enabled
                if self.plot_enabled.get():
                    self._plot_namd_data()
            else:
                messagebox.showerror("Error", "Could not parse NAMD log file or no data found.")
    
    def _auto_load_log_files(self):
        """Automatically load all log files from the trajectory_files list."""
        if not self.trajectory_files:
            raise ValueError("No log files to load")
        
        # Parse all log files and combine them
        all_data = {'timestep': []}
        first_file = True
        total_points = 0
        failed_files = []
        available_columns_set = set()  # Track all unique columns found
        
        # Track file data counts and timestep ranges for time calculation
        self._file_data_counts = {}
        self._file_timestep_ranges = {}
        current_idx = 0
        
        for log_file_path in self.trajectory_files:
            try:
                data = self.parse_namd_log(log_file_path)
                if data and data['timestep']:
                    num_points = len(data['timestep'])
                    
                    # Track available columns from this file
                    for key in data.keys():
                        if key != 'timestep':
                            available_columns_set.add(key)
                    
                    # Track data count for this file
                    self._file_data_counts[log_file_path] = num_points
                    
                    # Track timestep range for this file
                    if num_points > 0:
                        timesteps = data['timestep']
                        min_ts = min(timesteps)
                        max_ts = max(timesteps)
                        self._file_timestep_ranges[log_file_path] = (current_idx, current_idx + num_points, min_ts, max_ts)
                        print(f"File {Path(log_file_path).name}: indices {current_idx}-{current_idx + num_points}, timesteps {min_ts}-{max_ts}")
                        current_idx += num_points
                    
                    if first_file:
                        # First file: initialize all columns
                        all_data = data
                        first_file = False
                        total_points = len(data['timestep'])
                    else:
                        # Subsequent files: append data
                        all_data['timestep'].extend(data['timestep'])
                        for key in data:
                            if key != 'timestep' and key in all_data:
                                all_data[key].extend(data[key])
                        total_points += len(data['timestep'])
                else:
                    failed_files.append(Path(log_file_path).name)
                        
            except Exception as e:
                logger.error(f"Failed to parse log file {log_file_path}: {e}")
                failed_files.append(Path(log_file_path).name)
                continue
        
        if total_points > 0:
            self.analysis_results = all_data
            
            # Populate namd_available_columns with formatted names
            column_mapping = {
                'TEMP': 'Temperature (K)',
                'TEMPERATURE': 'Temperature (K)',
                'POTENTIAL': 'Potential Energy (kcal/mol)',
                'POTENG': 'Potential Energy (kcal/mol)',
                'KINETIC': 'Kinetic Energy (kcal/mol)',
                'KINENG': 'Kinetic Energy (kcal/mol)',
                'TOTAL': 'Total Energy (kcal/mol)',
                'TOTENG': 'Total Energy (kcal/mol)',
                'PRESSURE': 'Pressure (atm)',
                'VOLUME': 'Volume (Å³)',
                'ELECT': 'Electrostatic Energy (kcal/mol)',
                'VDW': 'Van der Waals Energy (kcal/mol)',
                'BOND': 'Bond Energy (kcal/mol)',
                'ANGLE': 'Angle Energy (kcal/mol)',
                'DIHED': 'Dihedral Energy (kcal/mol)',
                'IMPR': 'Improper Energy (kcal/mol)',
            }
            
            self.namd_available_columns = []
            for col in sorted(available_columns_set):
                display_name = column_mapping.get(col, f"{col} (raw)")
                self.namd_available_columns.append(display_name)
            
            print(f"Populated namd_available_columns with {len(self.namd_available_columns)} columns: {self.namd_available_columns}")
            
            num_files = len(self.trajectory_files)
            file_word = "file" if num_files == 1 else "files"
            self.status_callback(f"Loaded {num_files} log {file_word} ({total_points} data points)")
            
            print(f"Created _file_timestep_ranges with {len(self._file_timestep_ranges)} files")
            
            # Warn about failed files if any
            if failed_files:
                logger.warning(f"Failed to load {len(failed_files)} files: {', '.join(failed_files)}")
        else:
            # Clear results if no data was loaded
            self.analysis_results = {}
            error_msg = f"Could not load any data from the {len(self.trajectory_files)} log file(s)"
            if failed_files:
                error_msg += f": {', '.join(failed_files)}"
            raise ValueError(error_msg)

    def get_namd_column_mapping(self):
        """Get mapping of user-friendly column names to NAMD log data keys."""
        return {
            "Temperature (K)": "temp",
            "Potential Energy (kcal/mol)": "poteng",
            "Kinetic Energy (kcal/mol)": "kineng", 
            "Total Energy (kcal/mol)": "toteng",
            "Pressure (atm)": "pressure",
            "Volume (Å³)": "volume"
        }

    def get_selected_namd_data(self):
        """Get the data array for the currently selected NAMD column."""
        if not self.analysis_results or not hasattr(self, 'namd_column_combo'):
            return None, None
            
        selected_display_name = self.namd_column_combo.get()
        
        # Map display name back to data key
        data_key = self._get_data_key_from_display_name(selected_display_name)
        
        if data_key and data_key in self.analysis_results and 'timestep' in self.analysis_results:
            timesteps = self.analysis_results['timestep']
            y_data = self.analysis_results[data_key]
            
            # Convert timesteps to time in picoseconds
            # Assuming 2 fs timestep and output frequency
            timestep_fs = 2.0  # NAMD timestep in femtoseconds
            output_freq = 1000  # How often ENERGY lines are written (steps)
            
            # Convert to picoseconds
            x_data = [ts * timestep_fs * output_freq / 1000.0 for ts in timesteps]
            
            return x_data, y_data
        
        return None, None
    
    def _get_data_key_from_display_name(self, display_name):
        """Convert display name back to data key."""
        if not hasattr(self, 'analysis_results') or not self.analysis_results:
            return None
            
        # Reverse mapping from display name to data key
        reverse_mapping = {
            'Temperature (K)': ['TEMP', 'TEMPERATURE'],
            'Potential Energy (kcal/mol)': ['POTENTIAL', 'POTENG'],
            'Kinetic Energy (kcal/mol)': ['KINETIC', 'KINENG'],
            'Total Energy (kcal/mol)': ['TOTAL', 'TOTENG'],
            'Pressure (atm)': ['PRESSURE'],
            'Volume (Å³)': ['VOLUME'],
            'Electrostatic Energy (kcal/mol)': ['ELECT'],
            'Van der Waals Energy (kcal/mol)': ['VDW'],
            'Bond Energy (kcal/mol)': ['BOND'],
            'Angle Energy (kcal/mol)': ['ANGLE'],
            'Dihedral Energy (kcal/mol)': ['DIHED'],
            'Improper Energy (kcal/mol)': ['IMPR'],
            'Miscellaneous Energy (kcal/mol)': ['MISC'],
            'Boundary Energy (kcal/mol)': ['BOUNDARY'],
            'Gas Pressure (atm)': ['GPRESSURE'],
            'Average Pressure (atm)': ['PRESSAVG'],
            'Average Gas Pressure (atm)': ['GPRESSAVG']
        }
        
        if display_name in reverse_mapping:
            for possible_key in reverse_mapping[display_name]:
                if possible_key in self.analysis_results:
                    return possible_key
        
        # If not found in mapping, try to extract from display name format
        # For custom columns like "Column_name (unit)"
        if '(' in display_name:
            base_name = display_name.split('(')[0].strip()
            # Try exact match first
            for key in self.analysis_results.keys():
                if key.upper() == base_name.upper():
                    return key
        
        return None
    
    def _calculate_time_for_namd_logs(self):
        """Calculate time array for NAMD log files with user-assigned time ranges.
        
        Distributes all data points from each file evenly across the assigned time.
        Ignores the actual timestep values from the TS column - simply spreads
        all points uniformly over the assigned duration.
        
        Returns:
            list: Time values in picoseconds for all data points
        """
        time_array_ps = []
        
        print("\n=== Calculating time array for NAMD logs ===")
        print(f"Number of log files: {len(self.trajectory_files)}")
        print(f"Assigned times: {[(Path(f).name, t) for f, t in self.file_times.items()]}")
        
        # Check if we have the necessary data
        if not hasattr(self, '_file_timestep_ranges') or not self.analysis_results:
            print("WARNING: Using fallback method - file ranges not found")
            # Fallback: distribute evenly
            total_points = len(self.analysis_results.get('timestep', []))
            total_time_ns = sum(self.file_times.values())
            total_time_ps = total_time_ns * 1000.0
            time_per_point = total_time_ps / max(total_points - 1, 1)
            result = [i * time_per_point for i in range(total_points)]
            print(f"Fallback: {total_points} points over {total_time_ps}ps")
            return result
        
        # Track cumulative time across all files
        cumulative_time_ps = 0.0
        
        # Process each file in order
        for filepath in self.trajectory_files:
            print(f"\nProcessing file: {Path(filepath).name}")
            
            if filepath not in self._file_timestep_ranges:
                print(f"  ERROR: File not in _file_timestep_ranges!")
                continue
            
            start_idx, end_idx, min_ts, max_ts = self._file_timestep_ranges[filepath]
            num_points = end_idx - start_idx
            
            print(f"  Data indices: {start_idx} to {end_idx} ({num_points} points)")
            print(f"  Timestep range: {min_ts} to {max_ts}")
            
            if num_points == 0:
                print(f"  WARNING: No points for this file!")
                continue
            
            # Get assigned time for this file in ps
            assigned_time_ps = self.file_times.get(filepath, 0.0) * 1000.0
            print(f"  Assigned time: {assigned_time_ps}ps ({assigned_time_ps/1000.0}ns)")
            
            if assigned_time_ps <= 0:
                # No time assigned, use default spacing (1 ps per point)
                print(f"  Using default spacing (1 ps/point)")
                for i in range(num_points):
                    time_array_ps.append(cumulative_time_ps + i * 1.0)
                cumulative_time_ps += num_points * 1.0
            else:
                # Distribute all points evenly across the assigned time range
                # Use linspace to spread num_points from cumulative_time to cumulative_time + assigned_time
                if num_points == 1:
                    # Single point at the start
                    time_array_ps.append(cumulative_time_ps)
                    print(f"  Single point at {cumulative_time_ps}ps")
                else:
                    # Multiple points distributed evenly
                    file_times = np.linspace(cumulative_time_ps, cumulative_time_ps + assigned_time_ps, num_points)
                    time_array_ps.extend(file_times.tolist())
                    print(f"  Distributed {num_points} points from {cumulative_time_ps}ps to {cumulative_time_ps + assigned_time_ps}ps")
                    print(f"  First few times: {file_times[:3].tolist()}")
                    print(f"  Last few times: {file_times[-3:].tolist()}")
                
                # Update cumulative time by the assigned duration
                cumulative_time_ps += assigned_time_ps
        
        print(f"\nTotal time array: {len(time_array_ps)} points")
        print(f"Time range: {min(time_array_ps) if time_array_ps else 'N/A'} to {max(time_array_ps) if time_array_ps else 'N/A'} ps")
        print("=== Time array calculation complete ===\n")
        
        return time_array_ps
    
    def _plot_namd_data(self):
        """Plot NAMD log data based on selected column(s)."""
        # Check if we have log files
        if not self.trajectory_files:
            messagebox.showerror("No Files", "No log files added. Please add log files before plotting.")
            return
        
        # Check if we need to reload log files
        # (either no data exists, or file list has changed)
        needs_reload = False
        if not self.analysis_results or 'timestep' not in self.analysis_results:
            needs_reload = True
        elif hasattr(self, '_loaded_log_files'):
            # Check if the file list has changed
            if set(self.trajectory_files) != set(self._loaded_log_files):
                needs_reload = True
        else:
            needs_reload = True
        
        # Load/reload log files if needed
        if needs_reload:
            try:
                self._auto_load_log_files()
                # Track which files were loaded
                self._loaded_log_files = self.trajectory_files.copy()
            except Exception as e:
                messagebox.showerror("Loading Failed", f"Failed to load log files:\n\n{str(e)}\n\nPlease check your log files and try again.")
                return
        
        # Verify data was loaded successfully
        if not self.analysis_results or 'timestep' not in self.analysis_results:
            messagebox.showerror("No Data", "Failed to load data from log files. Please check that the files are valid NAMD log files.")
            return
        
        # Get all selected columns
        selected_columns = self._get_selected_columns()
        if not selected_columns:
            messagebox.showwarning("No Column Selected", "Please select at least one column to plot.\n\nUse the 'Data to Plot' dropdown or check boxes to select columns.")
            return
        if selected_columns == ["No data available"]:
            messagebox.showwarning("No Data", "No data available for plotting. Please check your log files.")
            return
        
        # Ensure we have a plot canvas
        if not hasattr(self, 'ax') or self.ax is None:
            self.create_plot_display()
        
        # Clear previous plot
        if self.ax:
            self.ax.clear()
        
        # Configure subplot spacing for energetic analysis (more space for legend)
        self._configure_subplot_spacing(analysis_type="energetic")
        
        # Get user-selected line color - support custom RGB
        line_color = self.plot_line_color.get() if hasattr(self, 'plot_line_color') else 'blue'
        if line_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_line_color_entry'):
                custom_color = self.custom_line_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    line_color = custom_color
                else:
                    line_color = "blue"  # Fallback
            else:
                line_color = "blue"
        
        # Colors for multiple lines (when plotting multiple columns)
        colors = [line_color, 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'gray']
        
        # Get time data (x-axis)
        timesteps = self.analysis_results['timestep']
        
        # Check if user has assigned times for files (any non-zero time)
        has_assigned_times = bool(self.file_times) and any(time > 0 for time in self.file_times.values())
        
        print(f"\n=== Plotting NAMD data ===")
        print(f"Using custom time assignments: {has_assigned_times}")
        print(f"Total data points: {len(timesteps)}")
        
        if has_assigned_times:
            # Distribute data points across user-assigned time ranges for each file
            x_data = self._calculate_time_for_namd_logs()
        else:
            # Use the timestep column from the log files
            timestep_fs = 2.0  # NAMD timestep in femtoseconds
            output_freq = 1000  # How often ENERGY lines are written (steps)
            x_data = [ts * timestep_fs * output_freq / 1000.0 for ts in timesteps]  # ps
            print(f"Using timestep-based x-axis")
        
        # Verify x_data and y_data will have matching lengths
        if len(x_data) != len(timesteps):
            print(f"ERROR: x_data length ({len(x_data)}) != y_data length ({len(timesteps)})")
            messagebox.showerror("Data Mismatch", f"Time array has {len(x_data)} points but data has {len(timesteps)} points. Cannot plot.")
            return
        
        # Convert time units based on user selection
        time_units = self.plot_time_units.get()
        if time_units == "ns":
            x_data = [x / 1000 for x in x_data]  # Convert ps to ns
        elif time_units == "µs":
            x_data = [x / 1000000 for x in x_data]  # Convert ps to µs
        
        print(f"Plot time range: {min(x_data):.2f} to {max(x_data):.2f} {time_units}")
        print(f"=== Plotting complete ===\n")
        
        # Plot each selected column
        lines = []
        labels = []
        
        # Get current Y-axis units and property type
        current_yaxis_unit = self.plot_yaxis_units.get()
        property_type = self.current_property_type
        
        # List of all NAMD energy column keys (case-insensitive matching)
        energy_column_keys = [
            'poteng', 'potential', 'kineng', 'kinetic', 'toteng', 'total', 'total3',
            'elect', 'vdw', 'bond', 'angle', 'dihed', 'impr', 'misc', 'boundary'
        ]
        
        for i, column_display_name in enumerate(selected_columns):
            data_key = self._get_data_key_from_display_name(column_display_name)
            if data_key and data_key in self.analysis_results:
                y_data = self.analysis_results[data_key]
                if y_data:  # Check if data is not empty
                    # Apply unit conversion based on property type
                    y_data = self._convert_units(y_data, data_key, property_type, current_yaxis_unit)
                    
                    color = colors[i % len(colors)]
                    line = self.ax.plot(x_data, y_data, color=color, linewidth=1.5, 
                                      marker='o', markersize=2, label=column_display_name.split('(')[0].strip())
                    lines.extend(line)
                    labels.append(column_display_name.split('(')[0].strip())
        
        # Set labels and title
        self.ax.set_xlabel(self.plot_xlabel.get() or f"Time ({time_units})")
        
        # Set y-label based on user setting or auto-detect from columns
        ylabel = self.plot_ylabel.get()
        if not ylabel:
            # Auto-generate ylabel if not set by user
            if len(selected_columns) == 1:
                # Extract property name from display name
                display_name = selected_columns[0]
                property_name = display_name.split('(')[0].strip()
                
                # Use the current Y-axis unit from the combobox
                current_yaxis_unit = self.plot_yaxis_units_combo.get()
                ylabel = f"{property_name} ({current_yaxis_unit})"
            else:
                ylabel = "Multiple Properties"
        self.ax.set_ylabel(ylabel)
        
        title = self.plot_title.get() or f"NAMD Analysis - {len(selected_columns)} Properties"
        self.ax.set_title(title)
        
        # Add legend for multiple columns
        if len(selected_columns) > 1:
            self.ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        # Set figure background color (the frame around the plot)
        fig_bg_color = self.plot_figure_bg_color.get()
        if fig_bg_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_fig_bg_color_entry'):
                custom_color = self.custom_fig_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    fig_bg_color = custom_color
                else:
                    fig_bg_color = "#212121"  # Fallback
            else:
                fig_bg_color = "#212121"
        elif fig_bg_color == "Transparent":
            fig_bg_color = "none"  # Matplotlib's transparent value
        elif "(" in fig_bg_color and ")" in fig_bg_color:
            # Extract hex code from "Name (#code)" format
            fig_bg_color = fig_bg_color.split("(")[1].split(")")[0]
        
        # Apply figure background
        if hasattr(self, 'figure') and self.figure is not None:
            self.figure.patch.set_facecolor(fig_bg_color)
        
        # Get background color - support preset, transparent, and custom RGB
        bg_color = self.plot_background_color.get()
        if bg_color == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_bg_color_entry'):
                custom_color = self.custom_bg_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    bg_color = custom_color
                else:
                    bg_color = "#2b2b2b"  # Fallback to dark gray
            else:
                bg_color = "#2b2b2b"
        elif bg_color == "Transparent":
            # Use matplotlib's transparent value for the plot background
            bg_color = "none"
        elif "(" in bg_color and ")" in bg_color:
            # Extract hex code from "Name (#code)" format
            bg_color = bg_color.split("(")[1].split(")")[0]
        
        # Apply styling with user-selected background
        self.ax.set_facecolor(bg_color)
        
        # Get text color - support Auto (based on luminance), preset colors, and custom RGB
        text_color_setting = self.plot_text_color.get()
        if text_color_setting == "Auto":
            # Set text and spine colors based on background luminance
            # For transparent backgrounds, default to black text
            if bg_color == "none":
                text_color = 'black'
            else:
                # Calculate luminance to determine if background is light or dark
                try:
                    # Remove '#' and convert hex to RGB
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    # Calculate relative luminance
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    # Fallback for light/dark detection
                    if bg_color in ["#ffffff", "#f0f0f0"]:
                        text_color = 'black'
                    else:
                        text_color = 'white'
        elif text_color_setting == "Custom RGB":
            # Use custom RGB value from entry
            if hasattr(self, 'custom_text_color_entry'):
                custom_color = self.custom_text_color_entry.get().strip()
                if custom_color and custom_color.startswith('#') and len(custom_color) == 7:
                    text_color = custom_color
                else:
                    text_color = 'black'  # Fallback
            else:
                text_color = 'black'
        elif "(" in text_color_setting and ")" in text_color_setting:
            # Extract hex code from "Name (#code)" format
            text_color = text_color_setting.split("(")[1].split(")")[0]
        else:
            text_color = text_color_setting  # Use as is
        
        self.ax.tick_params(colors=text_color)
        for spine in self.ax.spines.values():
            spine.set_color(text_color)
        self.ax.xaxis.label.set_color(text_color)
        self.ax.yaxis.label.set_color(text_color)
        self.ax.title.set_color(text_color)
        
        # Enable/disable grid based on user setting
        if self.plot_show_grid.get():
            self.ax.grid(True, alpha=0.3, color=text_color)
        
        # Set axis limits if specified, otherwise use data range
        # X-axis limits: Use specified values or default to data range
        try:
            xlim_min = self.plot_xlim_min.get().strip()
            xlim_max = self.plot_xlim_max.get().strip()
            
            # Get current data limits as fallback
            current_xlim = self.ax.get_xlim()
            
            # Set X limits: use specified value or current data limit
            x_min = float(xlim_min) if xlim_min else current_xlim[0]
            x_max = float(xlim_max) if xlim_max else current_xlim[1]
            self.ax.set_xlim(x_min, x_max)
        except (ValueError, AttributeError):
            pass
            
        # Y-axis limits: Use specified values or default to data range
        try:
            ylim_min = self.plot_ylim_min.get().strip()
            ylim_max = self.plot_ylim_max.get().strip()
            
            # Get current data limits as fallback
            current_ylim = self.ax.get_ylim()
            
            # Set Y limits: use specified value or current data limit
            y_min = float(ylim_min) if ylim_min else current_ylim[0]
            y_max = float(ylim_max) if ylim_max else current_ylim[1]
            self.ax.set_ylim(y_min, y_max)
        except (ValueError, AttributeError):
            pass
        
        # Refresh the plot
        if hasattr(self, 'canvas') and self.canvas is not None:
            self.canvas.draw()
            self.canvas.flush_events()
    
    def cleanup(self):
        """Cleanup method called when the frame is being destroyed."""
        try:
            # Unbind all events to prevent lambda command errors on shutdown
            # This recursively unbinds all events from all widgets in this frame
            def unbind_all(widget):
                try:
                    # Unbind all event sequences
                    for sequence in ['<FocusOut>', '<Return>', '<Button-1>', '<Button-4>', '<Button-5>', 
                                   '<KeyRelease>', '<MouseWheel>', '<Configure>']:
                        try:
                            widget.unbind(sequence)
                        except:
                            pass
                    
                    # Recursively unbind children
                    if hasattr(widget, 'winfo_children'):
                        for child in widget.winfo_children():
                            unbind_all(child)
                except:
                    pass
            
            unbind_all(self)
        except Exception as e:
            # Silently ignore cleanup errors
            pass
