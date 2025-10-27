# gatewizard/gui/frames/equilibration.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Equilibration frame for molecular dynamics simulations.

This module provides the GUI for setting up and running equilibration protocols
using various molecular simulation engines, starting with NAMD.
"""

import tkinter as tk
from tkinter import messagebox, filedialog
from typing import Optional, Callable, List, Dict, Any
from pathlib import Path
import json
import threading
import shutil

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, WIDGET_SIZES, LAYOUT, ERROR_MESSAGES, SUCCESS_MESSAGES
)
from gatewizard.gui.widgets.progress_tracker import ProgressTracker
from gatewizard.tools.equilibration import NAMDEquilibrationManager, EquilibrationProtocol
from gatewizard.utils.config import set_working_directory
from gatewizard.utils.logger import get_logger
from gatewizard.utils.namd_analysis import get_equilibration_progress, format_timing_info, format_progress_summary

logger = get_logger(__name__)

class EquilibrationFrame(ctk.CTkFrame):
    """
    Frame for molecular dynamics equilibration setup and execution.
    
    This frame handles the setup and execution of equilibration protocols
    for membrane protein systems using various MD engines.
    """
    
    def __init__(
        self,
        parent,
        get_current_pdb: Optional[Callable[[], Optional[str]]] = None,
        status_callback: Optional[Callable[[str], None]] = None,
        initial_directory: Optional[str] = None
    ):
        """
        Initialize the equilibration frame.
        
        Args:
            parent: Parent widget
            get_current_pdb: Callback to get current PDB file
            status_callback: Callback for status updates
            initial_directory: Initial directory for file dialogs
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['content_bg'])
        
        # Initialize logger
        self.logger = get_logger(self.__class__.__name__)
        
        self.get_current_pdb = get_current_pdb
        self.status_callback = status_callback
        self.initial_directory = initial_directory or str(Path.cwd())
        
        # State variables
        self.current_pdb_file = None
        self.working_directory = Path(self.initial_directory)
        self.equilibration_output_name = "equilibration"  # Default output folder name
        self.equilibration_stages = []
        self.custom_protocols = None  # Storage for custom protocol configurations
        
        # Progress monitoring variables
        self.monitoring_active = False
        self.progress_timer = None
        
        # Default AMBER protocol parameters
        self.default_protocols = self._get_default_protocols()
        
        # Create widgets
        self._create_widgets()
        self._create_outputname_section()  # Add output name section
        self._setup_layout()
        self._load_defaults()
        
        # Enable mouse wheel scrolling
        self._bind_mouse_wheel()
        
        # Update status
        if self.status_callback:
            self.status_callback("Equilibration frame initialized")

    def _bind_mouse_wheel(self):
        """Bind mouse wheel events for scrolling."""
        def _on_mousewheel(event):
            # Get the widget under the mouse
            widget = event.widget
            
            # Check if we're directly on a CTkTextbox or its internal text widget
            current = widget
            textbox_widget = None
            
            # First, check if we're inside a CTkTextbox
            while current and current != self:
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
                        return "break"
                except Exception:
                    pass
            
            # Otherwise, traverse up to find if we're inside a stage scrollable frame
            current = widget
            stage_scrollable = None
            main_scrollable = None
            
            while current and current != self:
                if isinstance(current, ctk.CTkScrollableFrame):
                    if hasattr(self, 'main_scroll') and current == self.main_scroll:
                        main_scrollable = current
                    else:
                        # This is a stage scrollable frame
                        stage_scrollable = current
                        break
                current = getattr(current, 'master', None)
            
            # If we're in a stage scrollable frame, let it handle the scrolling
            if stage_scrollable and hasattr(stage_scrollable, '_parent_canvas'):
                stage_scrollable._parent_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
            # Otherwise, scroll the main frame
            elif main_scrollable and hasattr(main_scrollable, '_parent_canvas'):
                main_scrollable._parent_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
            return "break"
        
        def _on_mousewheel_linux(event):
            # Get the widget under the mouse
            widget = event.widget
            
            # Check if we're directly on a CTkTextbox or its internal text widget
            current = widget
            textbox_widget = None
            
            # First, check if we're inside a CTkTextbox
            while current and current != self:
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
                        return "break"
                except Exception:
                    pass
            
            # Traverse up to find if we're inside a stage scrollable frame
            current = widget
            stage_scrollable = None
            main_scrollable = None
            
            while current and current != self:
                if isinstance(current, ctk.CTkScrollableFrame):
                    if hasattr(self, 'main_scroll') and current == self.main_scroll:
                        main_scrollable = current
                    else:
                        # This is a stage scrollable frame
                        stage_scrollable = current
                        break
                current = getattr(current, 'master', None)
            
            # If we're in a stage scrollable frame, let it handle the scrolling
            if stage_scrollable and hasattr(stage_scrollable, '_parent_canvas'):
                if event.num == 4:
                    stage_scrollable._parent_canvas.yview_scroll(-1, "units")
                elif event.num == 5:
                    stage_scrollable._parent_canvas.yview_scroll(1, "units")
            # Otherwise, scroll the main frame
            elif main_scrollable and hasattr(main_scrollable, '_parent_canvas'):
                if event.num == 4:
                    main_scrollable._parent_canvas.yview_scroll(-1, "units")
                elif event.num == 5:
                    main_scrollable._parent_canvas.yview_scroll(1, "units")
            return "break"
        
        # Bind to all widgets recursively
        def bind_to_widget(widget):
            try:
                # Windows and macOS
                widget.bind("<MouseWheel>", _on_mousewheel)
                # Linux
                widget.bind("<Button-4>", _on_mousewheel_linux)
                widget.bind("<Button-5>", _on_mousewheel_linux)
                
                # Recursively bind to all children
                for child in widget.winfo_children():
                    bind_to_widget(child)
            except Exception:
                pass  # Some widgets may not support binding
        
        # Bind after a short delay to ensure widgets are created
        self.after(100, lambda: bind_to_widget(self))
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            if hasattr(self, '_refresh_timer') and self._refresh_timer:
                self.after_cancel(self._refresh_timer)
                self._refresh_timer = None
            
            if hasattr(self, 'progress_timer') and self.progress_timer:
                self.after_cancel(self.progress_timer)
                self.progress_timer = None
                
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")
    
    def _get_default_protocols(self) -> Dict[str, Any]:
        """Get default equilibration protocols for AMBER force field."""
        return {
            "equilibration_1": {
                "name": "Equilibration 1 - Strong Restraints",
                "description": "Initial equilibration with strong restraints",
                "time_ns": 0.25,  # 0.25 ns (250 ps) equivalent to 125000 steps at 2 fs timestep
                "ensemble": "NPT",
                "temperature": 310.15,
                "pressure": 1.0,
                "minimize_steps": 10000,  # Initial minimization steps
                "constraints": {
                    "protein_backbone": 5.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 5.0,
                    "lipid_tail": 5.0,
                    "water": 0.0,
                    "ions": 5.0,
                    "other": 0.0
                },
                "timestep": 1.0,  # 1 fs (early equilibration stages use smaller timestep)
                "dcd_freq": 5000,  # DCD trajectory output frequency
                "use_gpu": True,  # Equilibration stages with GPU
                "cpu_cores": 1,
                "gpu_id": 0,      # GPU device ID
                "num_gpus": 1     # Number of GPUs
            },
            "equilibration_2": {
                "name": "Equilibration 2 - Medium Restraints",
                "description": "Equilibration with medium restraints",
                "time_ns": 0.25,  # 0.25 ns (250 ps) equivalent to 125000 steps at 2 fs timestep
                "ensemble": "NPT",
                "temperature": 310.15,
                "pressure": 1.0,
                "constraints": {
                    "protein_backbone": 5.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 2.5,
                    "lipid_tail": 2.5,
                    "water": 0.0,
                    "ions": 0.0,
                    "other": 0.0
                },
                "timestep": 1.0,  # 1 fs (early equilibration stages use smaller timestep)
                "dcd_freq": 5000,  # DCD trajectory output frequency
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "equilibration_3": {
                "name": "Equilibration 3 - Light Restraints",
                "description": "Equilibration with light restraints",
                "time_ns": 0.25,  # 0.25 ns (250 ps) equivalent to 125000 steps at 2 fs timestep
                "ensemble": "NPAT",  # Use NPAT for membrane systems
                "temperature": 310.15,
                "pressure": 1.0,
                "surface_tension": 0.0,  # dyn/cm for membrane
                "constraints": {
                    "protein_backbone": 5.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 1.0,
                    "lipid_tail": 1.0,
                    "water": 0.0,
                    "ions": 0.0,
                    "other": 0.0
                },
                "timestep": 1.0,  # 1 fs (early equilibration stages use smaller timestep)
                "dcd_freq": 5000,  # DCD trajectory output frequency
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "equilibration_4": {
                "name": "Equilibration 4 - Protein Only",
                "description": "Equilibration with protein-only restraints",
                "time_ns": 0.5,  # 0.5 ns (500 ps) equivalent to 250000 steps at 2 fs timestep
                "ensemble": "NPAT",  # Use NPAT for membrane systems
                "temperature": 310.15,
                "pressure": 1.0,
                "surface_tension": 0.0,  # dyn/cm for membrane
                "constraints": {
                    "protein_backbone": 5.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 0.5,
                    "lipid_tail": 0.1,
                    "water": 0.0,
                    "ions": 0.0,
                    "other": 0.0
                },
                "timestep": 2.0,  # 2 fs (NAMD uses femtoseconds)
                "dcd_freq": 5000,  # DCD trajectory output frequency
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "production_prep": {
                "name": "Production Preparation - ABF",
                "description": "Short ABF simulation for phosphate positioning",
                "time_ns": 0.25,  # 0.25 ns (250 ps) equivalent to 125000 steps at 2 fs timestep
                "ensemble": "NPAT",  # Use NPAT for membrane systems
                "temperature": 310.15,
                "pressure": 1.0,
                "surface_tension": 0.0,  # dyn/cm for membrane
                "constraints": {
                    "protein_backbone": 5.0,
                    "protein_sidechain": 5.0,
                    "lipid_head": 0.5,
                    "lipid_tail": 0.0,
                    "water": 0.0,
                    "ions": 0.0,
                    "other": 0.0
                },
                "timestep": 2.0,  # 2 fs (NAMD uses femtoseconds)
                "dcd_freq": 5000,  # DCD trajectory output frequency
                "abf_coordinate": "phosphate_distance",
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            }
        }
    
    def _get_charmm_gui_protocols(self, scheme_type: str) -> Dict[str, Any]:
        """Get CHARMM-GUI equilibration protocols based on scheme type."""
        # Base protocol structure with 7 stages (6 equilibration + production)
        # Settings match the provided figure with specific restraint values
        base_protocols = {
            "Equilibration 1": {
                "name": "Equilibration 1",
                "description": "Initial equilibration with strong restraints",
                "time_ns": 0.125,  # 125 ps = 0.125 ns
                "steps": 125000,   # steps = time_ns * 1e6 / timestep = 0.125 * 1e6 / 1.0
                "ensemble": "NVT",  # Start with NVT for all schemes
                "temperature": 303.15,
                "pressure": 1.0,
                "minimize_steps": 10000,  # Minimize steps (only for Equilibration 1)
                "constraints": {
                    "protein_backbone": 10.0,  # Backbone
                    "protein_sidechain": 5.0,   # Side Chain
                    "lipid_head": 2.5,          # Lipid head
                    "lipid_tail": 2.5,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 10.0,               # Ions
                    "other": 0.0
                },
                "timestep": 1.0,  # 1.0 fs
                "dcd_freq": 5000,
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "Equilibration 2": {
                "name": "Equilibration 2",
                "description": "Second equilibration stage",
                "time_ns": 0.125,  # 125 ps = 0.125 ns
                "steps": 125000,   # steps = time_ns * 1e6 / timestep = 0.125 * 1e6 / 1.0
                "ensemble": "NVT",
                "temperature": 303.15,
                "pressure": 1.0,
                "constraints": {
                    "protein_backbone": 5.0,    # Backbone
                    "protein_sidechain": 2.5,   # Side Chain
                    "lipid_head": 2.5,          # Lipid head
                    "lipid_tail": 2.5,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 0.0,                # Ions
                    "other": 0.0
                },
                "timestep": 1.0,  # 1.0 fs
                "dcd_freq": 5000,
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "Equilibration 3": {
                "name": "Equilibration 3",
                "description": "Third equilibration stage",
                "time_ns": 0.125,  # 125 ps = 0.125 ns
                "steps": 125000,   # steps = time_ns * 1e6 / timestep = 0.125 * 1e6 / 1.0
                "ensemble": scheme_type,  # Now use the selected scheme
                "temperature": 303.15,
                "pressure": 1.0,
                "surface_tension": 0.0,
                "constraints": {
                    "protein_backbone": 2.5,    # Backbone
                    "protein_sidechain": 1.0,   # Side Chain
                    "lipid_head": 1.0,          # Lipid head
                    "lipid_tail": 1.0,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 0.0,                # Ions
                    "other": 0.0
                },
                "timestep": 1.0,  # 1.0 fs
                "dcd_freq": 5000,
                "margin": 5.0,
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "Equilibration 4": {
                "name": "Equilibration 4",
                "description": "Fourth equilibration stage",
                "time_ns": 0.5,  # 500 ps = 0.5 ns
                "steps": 250000,   # steps = time_ns * 1e6 / timestep = 0.5 * 1e6 / 2.0
                "ensemble": scheme_type,
                "temperature": 303.15,
                "pressure": 1.0,
                "surface_tension": 0.0,
                "constraints": {
                    "protein_backbone": 1,      # Backbone
                    "protein_sidechain": 0.5,   # Side Chain
                    "lipid_head": 0.5,          # Lipid head
                    "lipid_tail": 0.5,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 0.0,                # Ions
                    "other": 0.0
                },
                "timestep": 2.0,  # 2.0 fs
                "dcd_freq": 5000,
                "margin": 5.0,
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "Equilibration 5": {
                "name": "Equilibration 5",
                "description": "Fifth equilibration stage",
                "time_ns": 0.5,  # 500 ps = 0.5 ns
                "steps": 250000,   # steps = time_ns * 1e6 / timestep = 0.5 * 1e6 / 2.0
                "ensemble": scheme_type,
                "temperature": 303.15,
                "pressure": 1.0,
                "surface_tension": 0.0,
                "constraints": {
                    "protein_backbone": 0.5,    # Backbone
                    "protein_sidechain": 0.1,   # Side Chain
                    "lipid_head": 0.1,          # Lipid head
                    "lipid_tail": 0.1,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 0.0,                # Ions
                    "other": 0.0
                },
                "timestep": 2.0,  # 2.0 fs
                "dcd_freq": 5000,
                "margin": 5.0,
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "Equilibration 6": {
                "name": "Equilibration 6",
                "description": "Sixth equilibration stage",
                "time_ns": 0.5,  # 500 ps = 0.5 ns
                "steps": 250000,   # steps = time_ns * 1e6 / timestep = 0.5 * 1e6 / 2.0
                "ensemble": scheme_type,
                "temperature": 303.15,
                "pressure": 1.0,
                "surface_tension": 0.0,
                "constraints": {
                    "protein_backbone": 0.1,    # Backbone
                    "protein_sidechain": 0.0,   # Side Chain
                    "lipid_head": 0.0,          # Lipid head
                    "lipid_tail": 0.0,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 0.0,                # Ions
                    "other": 0.0
                },
                "timestep": 2.0,  # 2.0 fs
                "dcd_freq": 5000,
                "margin": 5.0,
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            },
            "Production": {
                "name": "Production",
                "description": "Production simulation",
                "time_ns": 100.0,  # 100 ns
                "steps": 50000000, # steps = time_ns * 1e6 / timestep = 100.0 * 1e6 / 2.0
                "ensemble": scheme_type,
                "temperature": 303.15,
                "pressure": 1.0,
                "surface_tension": 0.0,
                "constraints": {
                    "protein_backbone": 0.0,    # Backbone
                    "protein_sidechain": 0.0,   # Side Chain
                    "lipid_head": 0.0,          # Lipid head
                    "lipid_tail": 0.0,          # Lipid tail
                    "water": 0.0,               # Water
                    "ions": 0.0,
                    "other": 0.0
                },
                "timestep": 2.0,  # 2 fs (NAMD uses femtoseconds)
                # 2 fs
                "dcd_freq": 50000,  # DCD trajectory output frequency
                "margin": 5.0,  # NAMD margin parameter for NPAT simulations
                "use_gpu": True,
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 1
            }
        }
        
        return base_protocols
    
    def _get_current_protocol_template(self) -> Dict[str, Any]:
        """Get the current CHARMM-GUI protocol template based on scheme selection."""
        # Check if we have custom protocols (from adding stages)
        if hasattr(self, 'custom_protocols') and self.custom_protocols:
            return self.custom_protocols
        
        # Always use CHARMM-GUI protocols
        scheme_type = getattr(self, 'scheme_type_var', None)
        if scheme_type is None:
            # During initialization, use default NPT scheme
            return self._get_charmm_gui_protocols("NPT")
        else:
            return self._get_charmm_gui_protocols(scheme_type.get())
    
    def _create_widgets(self):
        """Create all widgets for the equilibration frame."""
        # Create scrollable main frame
        self.main_scroll = ctk.CTkScrollableFrame(self, fg_color="transparent")
        
        # Title section
        self._create_title_section()
        
        # Working directory section
        self._create_workdir_section()
        
        # Input folder section (for prepared system files)
        self._create_inputfolder_section()
        
        # Output name section
        self._create_outputname_section()
        
        # MD Engine selection
        self._create_engine_section()
        
        # Protocol configuration
        self._create_protocol_section()
        
        # Action buttons
        self._create_action_section()
        
        # Progress tracking
        self._create_progress_section()
    
    def _create_title_section(self):
        """Create the title section."""
        self.title_frame = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.title_label = ctk.CTkLabel(
            self.title_frame,
            text="Molecular Dynamics Equilibration",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.description_label = ctk.CTkLabel(
            self.title_frame,
            text="Setup and run equilibration protocols for membrane protein systems",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
    
    def _create_workdir_section(self):
        """Create the working directory section."""
        self.workdir_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.workdir_label = ctk.CTkLabel(
            self.workdir_section,
            text="Working Directory",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.workdir_frame = ctk.CTkFrame(self.workdir_section, fg_color="transparent")
        
        self.workdir_entry = ctk.CTkEntry(
            self.workdir_frame,
            width=400,
            font=FONTS['body'],
            placeholder_text="Select working directory..."
        )
        
        self.workdir_browse_btn = ctk.CTkButton(
            self.workdir_frame,
            text="Browse",
            width=80,
            command=self._browse_workdir
        )
    
    def _create_inputfolder_section(self):
        """Create the input folder section for prepared system files."""
        self.inputfolder_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.inputfolder_label = ctk.CTkLabel(
            self.inputfolder_section,
            text="Input Folder (Prepared System)",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.inputfolder_sublabel = ctk.CTkLabel(
            self.inputfolder_section,
            text="Folder containing system.prmtop, system.inpcrd, and system.pdb files",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )
        
        self.inputfolder_frame = ctk.CTkFrame(self.inputfolder_section, fg_color="transparent")
        
        self.inputfolder_entry = ctk.CTkEntry(
            self.inputfolder_frame,
            width=400,
            font=FONTS['body'],
            placeholder_text="Select folder with prepared system files..."
        )
        
        self.inputfolder_browse_btn = ctk.CTkButton(
            self.inputfolder_frame,
            text="Browse",
            width=80,
            command=self._browse_inputfolder
        )
    
    def _create_outputname_section(self):
        """Create the output name section."""
        self.outputname_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.outputname_label = ctk.CTkLabel(
            self.outputname_section,
            text="Output Folder Name",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.outputname_frame = ctk.CTkFrame(self.outputname_section, fg_color="transparent")
        
        self.outputname_entry = ctk.CTkEntry(
            self.outputname_frame,
            width=300,
            font=FONTS['body'],
            placeholder_text="Enter output folder name (default: equilibration)"
        )
        
        # Set default value
        self.outputname_entry.insert(0, self.equilibration_output_name)
        
        # Bind change event
        self.outputname_entry.bind('<KeyRelease>', self._on_outputname_changed)
        self.outputname_entry.bind('<FocusOut>', self._on_outputname_changed)
    
    def _create_engine_section(self):
        """Create the MD engine selection section."""
        self.engine_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.engine_label = ctk.CTkLabel(
            self.engine_section,
            text="Molecular Dynamics Engine",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.engine_var = ctk.StringVar(value="NAMD")
        self.engine_menu = ctk.CTkOptionMenu(
            self.engine_section,
            values=["NAMD", "GROMACS", "AMBER"],
            variable=self.engine_var,
            command=self._on_engine_changed
        )
        
        # Engine-specific settings frame
        self.engine_settings_frame = ctk.CTkFrame(self.engine_section, fg_color="transparent")
    
    def _create_protocol_section(self):
        """Create the equilibration protocol configuration section."""
        self.protocol_section = ctk.CTkFrame(
            self.main_scroll, 
            fg_color=COLOR_SCHEME['content_inside_bg'],
            height=650  # Set minimum height for better vertical visibility
        )
        
        self.protocol_label = ctk.CTkLabel(
            self.protocol_section,
            text="CHARMM-GUI Equilibration Protocol",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # CHARMM-GUI scheme selection (always enabled)
        self.charmm_gui_scheme_frame = ctk.CTkFrame(self.protocol_section, fg_color="transparent")
        
        # Scheme type selection (NVT, NPT, NPAT, NPgT)
        self.scheme_type_frame = ctk.CTkFrame(self.charmm_gui_scheme_frame, fg_color="transparent")
        
        self.scheme_type_label = ctk.CTkLabel(
            self.scheme_type_frame,
            text="Equilibration Scheme:",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.scheme_type_var = ctk.StringVar(value="NPT")
        self.scheme_type_menu = ctk.CTkOptionMenu(
            self.scheme_type_frame,
            values=["NVT", "NPT", "NPAT", "NPgT"],
            variable=self.scheme_type_var,
            command=self._on_scheme_type_changed,
            width=120
        )
        
        ## Bilayer thickness control
        #self.thickness_frame = ctk.CTkFrame(self.charmm_gui_scheme_frame, fg_color="transparent")
        #
        #self.thickness_label = ctk.CTkLabel(
        #    self.thickness_frame,
        #    text="Bilayer Thickness (Å):",
        #    font=FONTS['body'],
        #    text_color=COLOR_SCHEME['text']
        #)
        #
        #self.thickness_var = ctk.StringVar(value="39.1")  # Default POPC bilayer thickness
        #self.thickness_entry = ctk.CTkEntry(
        #    self.thickness_frame,
        #    textvariable=self.thickness_var,
        #    width=80,
        #    font=FONTS['body'],
        #    placeholder_text="37.0"
        #)
        
        # Stage configuration frame
        self.stages_config_frame = ctk.CTkFrame(
            self.protocol_section, 
            fg_color="transparent",
            height=550  # Increase minimum height for better visibility of stage columns
        )
        
        # Create stage configuration widgets
        self._create_stage_widgets()
    
    def _create_stage_widgets(self):
        """Create widgets for individual stage configuration."""
        self.stage_widgets = {}
        
        # Clear existing widgets
        for widget in self.stages_config_frame.winfo_children():
            widget.destroy()

        # Create synchronized scrolling option
        self.sync_scroll_frame = ctk.CTkFrame(self.stages_config_frame, fg_color="transparent")
        self.sync_scroll_frame.pack(fill="x", padx=0, pady=(0, 10))
        
        self.sync_scroll_var = ctk.BooleanVar(value=True)
        self.sync_scroll_checkbox = ctk.CTkCheckBox(
            self.sync_scroll_frame, 
            text="Synchronize vertical scrolling across all stages",
            variable=self.sync_scroll_var,
            command=self._on_sync_scroll_changed
        )
        self.sync_scroll_checkbox.pack(side="left", padx=(0, 20))

        # Container for the stage widgets
        self.stages_container = ctk.CTkFrame(self.stages_config_frame, fg_color="transparent")
        self.stages_container.pack(fill="both", expand=True, padx=0, pady=0)
        
        # Create the columns view (only view now)
        self._create_columns_view()
    
    def _on_sync_scroll_changed(self):
        """Handle synchronized scrolling toggle."""
        if hasattr(self, 'column_frames'):
            if self.sync_scroll_var.get():
                self._setup_synchronized_scrolling()
            else:
                self._remove_synchronized_scrolling()
                self._setup_individual_scrolling()
    
    def _setup_individual_scrolling(self):
        """Setup individual scrolling for each column when sync is disabled."""
        if not hasattr(self, 'column_frames'):
            return
        
        for i, column_frame in enumerate(self.column_frames):
            if hasattr(column_frame, '_parent_canvas'):
                canvas = column_frame._parent_canvas
                
                # Create individual scroll handler for this column only
                def create_individual_mousewheel_handler(canvas_ref):
                    def on_mousewheel(event):
                        # Calculate scroll delta
                        if event.delta:
                            delta = -1 * int(event.delta / 120)
                        else:
                            delta = -1 if event.num == 4 else 1
                        
                        # Apply scroll only to this canvas
                        try:
                            canvas_ref.yview_scroll(delta, "units")
                        except Exception:
                            pass  # Silently handle errors
                        
                        return "break"
                    
                    return on_mousewheel
                
                # Create the individual mousewheel handler for this column
                individual_handler = create_individual_mousewheel_handler(canvas)
                
                # Bind to all widgets within this column
                self._bind_scrolling_to_all_widgets(column_frame, individual_handler)
    
    def _bind_scrolling_to_all_widgets(self, parent_widget, mousewheel_handler):
        """Recursively bind scrolling events to all child widgets."""
        try:
            # Bind to the parent widget
            parent_widget.bind("<MouseWheel>", mousewheel_handler, add=True)
            parent_widget.bind("<Button-4>", mousewheel_handler, add=True)
            parent_widget.bind("<Button-5>", mousewheel_handler, add=True)
            
            # Recursively bind to all children
            for child in parent_widget.winfo_children():
                self._bind_scrolling_to_all_widgets(child, mousewheel_handler)
        except Exception:
            # Some widgets might not support binding, so we silently skip them
            pass
    
    def _unbind_scrolling_from_all_widgets(self, parent_widget):
        """Recursively unbind scrolling events from all child widgets."""
        try:
            # Unbind from the parent widget
            parent_widget.unbind("<MouseWheel>")
            parent_widget.unbind("<Button-4>")
            parent_widget.unbind("<Button-5>")
            
            # Recursively unbind from all children
            for child in parent_widget.winfo_children():
                self._unbind_scrolling_from_all_widgets(child)
        except Exception:
            # Some widgets might not support unbinding, so we silently skip them
            pass

    def _setup_synchronized_scrolling(self):
        """Setup synchronized scrolling for all columns."""
        if not hasattr(self, 'column_frames') or len(self.column_frames) < 2:
            return
        
        # Remove any existing bindings first
        self._remove_synchronized_scrolling()
        
        # Store original scroll commands - collect ALL first before setting up sync
        self.original_scroll_commands = []
        for column_frame in self.column_frames:
            if hasattr(column_frame, '_parent_canvas'):
                canvas = column_frame._parent_canvas
                # Store original yview command
                self.original_scroll_commands.append(canvas.yview)
            else:
                # Add placeholder for frames without canvas
                self.original_scroll_commands.append(None)
        
        # Now setup synchronized scrolling for each column
        for i, column_frame in enumerate(self.column_frames):
            if hasattr(column_frame, '_parent_canvas'):
                canvas = column_frame._parent_canvas
                
                # Get the original command for this canvas
                original_cmd = self.original_scroll_commands[i]
                if original_cmd is None:
                    continue
                
                # Create synchronized scroll function
                def create_sync_yview(canvas_index, orig_cmd):
                    def sync_yview(*args):
                        if not self.sync_scroll_var.get():
                            # If sync is disabled, just use original command
                            return orig_cmd(*args)
                        
                        # Apply to all canvases
                        for j, other_frame in enumerate(self.column_frames):
                            if hasattr(other_frame, '_parent_canvas') and self.original_scroll_commands[j] is not None:
                                try:
                                    # Apply same command to all canvases
                                    self.original_scroll_commands[j](*args)
                                except Exception:
                                    pass  # Silently handle errors
                    
                    return sync_yview
                
                # Replace the yview command
                canvas.yview = create_sync_yview(i, original_cmd)
                
                # Create mouse wheel handler for this column
                def create_mousewheel_handler(canvas_index):
                    def on_mousewheel(event):
                        if not self.sync_scroll_var.get():
                            return
                        
                        # Calculate scroll delta
                        if event.delta:
                            delta = -1 * int(event.delta / 120)
                        else:
                            delta = -1 if event.num == 4 else 1
                        
                        # Apply scroll to all canvases
                        for j, other_frame in enumerate(self.column_frames):
                            if hasattr(other_frame, '_parent_canvas'):
                                try:
                                    other_frame._parent_canvas.yview_scroll(delta, "units")
                                except Exception:
                                    pass  # Silently handle errors
                        
                        return "break"
                    
                    return on_mousewheel
                
                # Create the mousewheel handler for this column
                mousewheel_handler = create_mousewheel_handler(i)
                
                # Bind mousewheel to canvas
                canvas.bind("<MouseWheel>", mousewheel_handler, add=True)
                canvas.bind("<Button-4>", mousewheel_handler, add=True)
                canvas.bind("<Button-5>", mousewheel_handler, add=True)
                
                # Bind to all widgets within this column using recursive binding
                self._bind_scrolling_to_all_widgets(column_frame, mousewheel_handler)
    
    def _remove_synchronized_scrolling(self):
        """Remove synchronized scrolling bindings."""
        if hasattr(self, 'column_frames'):
            for i, column_frame in enumerate(self.column_frames):
                try:
                    # Restore original yview command if we have it
                    if (hasattr(self, 'original_scroll_commands') and 
                        i < len(self.original_scroll_commands) and
                        hasattr(column_frame, '_parent_canvas')):
                        column_frame._parent_canvas.yview = self.original_scroll_commands[i]
                    
                    # Unbind mousewheel events from canvas
                    if hasattr(column_frame, '_parent_canvas'):
                        canvas = column_frame._parent_canvas
                        canvas.unbind("<MouseWheel>")
                        canvas.unbind("<Button-4>")
                        canvas.unbind("<Button-5>")
                    
                    # Recursively unbind from all widgets in the column
                    self._unbind_scrolling_from_all_widgets(column_frame)
                    
                except Exception as e:
                    print(f"Unbinding error for column {i}: {e}")
        
        # Clear stored commands
        if hasattr(self, 'original_scroll_commands'):
            self.original_scroll_commands = []
    
    def _create_columns_view(self):
        """Create side-by-side columns view for stages."""
        # Clear existing view
        for widget in self.stages_container.winfo_children():
            widget.destroy()
            
        # Create a scrollable frame for columns
        self.columns_scroll = ctk.CTkScrollableFrame(
            self.stages_container, 
            orientation="horizontal",
            fg_color="transparent",
            height=500  # Set minimum height for better visibility
        )
        self.columns_scroll.pack(fill="both", expand=True)
        
        # Store column frames for synchronized scrolling
        self.column_frames = []
        
        # Determine which protocol to use based on current selection
        current_protocols = self._get_current_protocol_template()
        
        # Create columns for each stage
        for col, (stage_key, stage_data) in enumerate(current_protocols.items()):
            # Create column frame - make it scrollable vertically within the column
            column_frame = ctk.CTkScrollableFrame(
                self.columns_scroll, 
                fg_color=COLOR_SCHEME['content_inside_bg'],
                width=320,  # Increased width for better content fit
                height=480,  # Set a good height for visibility while keeping scrolling
                orientation="vertical"
            )
            column_frame.pack(side="left", fill="both", expand=True, padx=8, pady=8)
            
            # Store reference for synchronized scrolling
            self.column_frames.append(column_frame)
            
            # Stage header
            header_label = ctk.CTkLabel(
                column_frame,
                text=stage_data["name"],
                font=FONTS['heading'],
                text_color=COLOR_SCHEME['text']
            )
            header_label.pack(anchor="w", padx=10, pady=(8, 4))
            
            # Stage description
            desc_label = ctk.CTkLabel(
                column_frame,
                text=stage_data["description"],
                font=FONTS['small'],
                text_color=COLOR_SCHEME['text'],
                wraplength=280
            )
            desc_label.pack(anchor="w", padx=10, pady=(0, 8))
            
            # Create stage widgets in this column
            self._create_stage_content(column_frame, stage_key, stage_data)
        
        # Setup synchronized scrolling if enabled (with delay to ensure widgets are ready)
        if self.sync_scroll_var.get():
            # Use after_idle to ensure all widgets are fully created before setting up sync
            self.stages_container.after_idle(self._setup_synchronized_scrolling)
        else:
            # Setup individual scrolling for each column
            self.stages_container.after_idle(self._setup_individual_scrolling)
    
    def _create_stage_content(self, parent_frame, stage_key, stage_data):
        """Create the content widgets for a single stage."""
        # Store widgets for this stage
        stage_widgets = {}
        
        # Parameters frame
        params_frame = ctk.CTkFrame(parent_frame, fg_color="transparent")
        params_frame.pack(fill="x", padx=10, pady=5)
        
        # Time/Steps dual input - user can input either time or steps
        time_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        time_frame.pack(fill="x", pady=2)
        
        # Time in nanoseconds (first row)
        time_ns_frame = ctk.CTkFrame(time_frame, fg_color="transparent")
        time_ns_frame.pack(fill="x", pady=1)
        
        time_text_label = ctk.CTkLabel(time_ns_frame, text="Time:", width=120, anchor="w", font=FONTS['body'])
        time_text_label.pack(side="left")
        
        # Convert steps to time for initial display
        steps = stage_data.get("steps", 125000)
        timestep = stage_data.get("timestep", 0.001)  # 1 fs default
        time_ns = steps * timestep / 1e6  # Convert to ns
        
        time_var = ctk.StringVar(value=f"{time_ns:.3f}")
        time_entry = ctk.CTkEntry(time_ns_frame, textvariable=time_var, width=100)
        time_entry.pack(side="left", padx=5)
        
        time_unit_label = ctk.CTkLabel(time_ns_frame, text="ns", width=30, anchor="w", font=FONTS['body'])
        time_unit_label.pack(side="left", padx=5)
        
        # Steps (second row)
        steps_frame = ctk.CTkFrame(time_frame, fg_color="transparent")
        steps_frame.pack(fill="x", pady=1)
        
        steps_text_label = ctk.CTkLabel(steps_frame, text="Steps:", width=120, anchor="w", font=FONTS['body'])
        steps_text_label.pack(side="left")
        
        steps_var = ctk.StringVar(value=str(steps))
        steps_entry = ctk.CTkEntry(steps_frame, textvariable=steps_var, width=100)
        steps_entry.pack(side="left", padx=5)
        
        steps_unit_label = ctk.CTkLabel(steps_frame, text="steps", width=50, anchor="w", font=FONTS['body'])
        steps_unit_label.pack(side="left", padx=5)
        
        # Store widgets
        stage_widgets["time_ns"] = time_var
        stage_widgets["time_entry"] = time_entry
        stage_widgets["time_label"] = time_text_label
        stage_widgets["time_unit_label"] = time_unit_label
        stage_widgets["steps"] = steps_var
        stage_widgets["steps_entry"] = steps_entry
        stage_widgets["steps_label"] = steps_text_label
        stage_widgets["steps_unit_label"] = steps_unit_label
        
        # Flag to prevent recursive updates
        updating_time_steps = False
        
        # Add callback to update steps when time changes
        def update_steps_from_time(*args):
            nonlocal updating_time_steps
            if updating_time_steps:
                return
            try:
                updating_time_steps = True
                time_value = float(time_var.get())
                # Try to get timestep from the widget if it exists, otherwise use stage data
                if "timestep" in stage_widgets:
                    timestep_value = float(stage_widgets["timestep"].get())
                else:
                    timestep_value = stage_data.get("timestep", 1.0)
                calculated_steps = int(time_value * 1e6 / timestep_value)
                steps_var.set(str(calculated_steps))
            except (ValueError, ZeroDivisionError):
                steps_var.set("0")
            finally:
                updating_time_steps = False
        
        # Add callback to update time when steps change
        def update_time_from_steps(*args):
            nonlocal updating_time_steps
            if updating_time_steps:
                return
            try:
                updating_time_steps = True
                steps_value = int(float(steps_var.get()))
                # Try to get timestep from the widget if it exists, otherwise use stage data
                if "timestep" in stage_widgets:
                    timestep_value = float(stage_widgets["timestep"].get())
                else:
                    timestep_value = stage_data.get("timestep", 1.0)
                calculated_time = steps_value * timestep_value / 1e6  # Convert to ns
                time_var.set(f"{calculated_time:.3f}")
            except (ValueError, ZeroDivisionError):
                time_var.set("0.000")
            finally:
                updating_time_steps = False
        
        time_var.trace("w", update_steps_from_time)
        steps_var.trace("w", update_time_from_steps)
        
        # Add trace to refresh progress display when time changes
        def refresh_progress_on_time_change(*args):
            # Small delay to avoid rapid updates during typing
            if hasattr(self, '_refresh_timer') and self._refresh_timer:
                self.after_cancel(self._refresh_timer)
            self._refresh_timer = self.after(1000, self._refresh_equilibration_progress)
        
        time_var.trace("w", refresh_progress_on_time_change)
        
        # Minimize Steps - only for Equilibration 1
        if stage_key == "Equilibration 1":
            minimize_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
            minimize_frame.pack(fill="x", pady=2)
            
            minimize_text_label = ctk.CTkLabel(minimize_frame, text="Minimize:", width=120, anchor="w", font=FONTS['body'])
            minimize_text_label.pack(side="left")
            
            minimize_var = ctk.StringVar(value=str(stage_data.get("minimize_steps", 125000)))
            minimize_entry = ctk.CTkEntry(minimize_frame, textvariable=minimize_var, width=100)
            minimize_entry.pack(side="left", padx=5)
            
            minimize_info_label = ctk.CTkLabel(minimize_frame, text="steps", width=120, anchor="w", font=FONTS['small'])
            minimize_info_label.pack(side="left", padx=5)
            
            # Store widgets
            stage_widgets["minimize_steps"] = minimize_var
            stage_widgets["minimize_entry"] = minimize_entry
            stage_widgets["minimize_label"] = minimize_text_label
            stage_widgets["minimize_info_label"] = minimize_info_label
            stage_widgets["minimize_frame"] = minimize_frame
        
        # Ensemble - HIDDEN from GUI (uses template from selected equilibration scheme)
        # The ensemble is determined by the equilibration scheme (NPT, NPAT, NPgT, etc.)
        # and should not be changed per-stage to maintain consistency with CHARMM-GUI protocols
        ensemble_var = ctk.StringVar(value=stage_data["ensemble"])
        stage_widgets["ensemble"] = ensemble_var
        # Note: No ensemble_menu created - ensemble is fixed per scheme
        
        # Display current ensemble (read-only)
        ensemble_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        ensemble_frame.pack(fill="x", pady=2)
        
        ensemble_text_label = ctk.CTkLabel(
            ensemble_frame, 
            text=f"Ensemble: {stage_data['ensemble']}", 
            width=220, 
            anchor="w", 
            font=FONTS['body'],
            text_color=COLOR_SCHEME['inactive']
        )
        ensemble_text_label.pack(side="left")
        
        # Info label explaining ensemble is fixed
        ensemble_info_label = ctk.CTkLabel(
            ensemble_frame,
            text="(defined by equilibration scheme)",
            width=250,
            anchor="w",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )
        ensemble_info_label.pack(side="left", padx=5)
        
        stage_widgets["ensemble_label"] = ensemble_text_label
        stage_widgets["ensemble_info_label"] = ensemble_info_label
        # Container frames to allow show/hide
        stage_widgets["params_frame"] = params_frame
        stage_widgets["ensemble_frame"] = ensemble_frame
        
        # Temperature
        temp_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        temp_frame.pack(fill="x", pady=2)
        
        temp_text_label = ctk.CTkLabel(temp_frame, text="Temperature:", width=120, anchor="w", font=FONTS['body'])
        temp_text_label.pack(side="left")
        temp_var = ctk.StringVar(value=str(stage_data["temperature"]))
        temp_entry = ctk.CTkEntry(temp_frame, textvariable=temp_var, width=100)
        temp_entry.pack(side="left", padx=5)
        temp_unit_label = ctk.CTkLabel(temp_frame, text="K", width=50, anchor="w", font=FONTS['body'])
        temp_unit_label.pack(side="left", padx=5)
        stage_widgets["temperature"] = temp_var
        stage_widgets["temperature_entry"] = temp_entry
        stage_widgets["temperature_label"] = temp_text_label
        stage_widgets["temperature_unit_label"] = temp_unit_label
        
        # Pressure: always create; initial visibility depends on ensemble
        pressure_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        # Pack only if applicable initially
        if stage_data["ensemble"] in ["NPT", "NPAT", "NPgT"]:
            pressure_frame.pack(fill="x", pady=2)
        
        pressure_text_label = ctk.CTkLabel(pressure_frame, text="Pressure:", width=120, anchor="w", font=FONTS['body'])
        pressure_text_label.pack(side="left")
        pressure_var = ctk.StringVar(value=str(stage_data.get("pressure", 1.0)))
        pressure_entry = ctk.CTkEntry(pressure_frame, textvariable=pressure_var, width=100)
        pressure_entry.pack(side="left", padx=5)
        pressure_unit_label = ctk.CTkLabel(pressure_frame, text="bar", width=50, anchor="w", font=FONTS['body'])
        pressure_unit_label.pack(side="left", padx=5)
        stage_widgets["pressure"] = pressure_var
        stage_widgets["pressure_entry"] = pressure_entry
        stage_widgets["pressure_label"] = pressure_text_label
        stage_widgets["pressure_unit_label"] = pressure_unit_label
        stage_widgets["pressure_frame"] = pressure_frame
        # Hide if not applicable initially
        if stage_data["ensemble"] not in ["NPT", "NPAT", "NPgT"]:
            pressure_frame.pack_forget()

        # Surface tension: always create; visible only for NPAT/NPgT
        tension_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        if stage_data["ensemble"] in ["NPAT", "NPgT"]:
            tension_frame.pack(fill="x", pady=2)
        tension_text_label = ctk.CTkLabel(tension_frame, text="Surface Tension:", width=120, anchor="w", font=FONTS['body'])
        tension_text_label.pack(side="left")
        tension_var = ctk.StringVar(value=str(stage_data.get("surface_tension", 0.0)))
        tension_entry = ctk.CTkEntry(tension_frame, textvariable=tension_var, width=100)
        tension_entry.pack(side="left", padx=5)
        tension_unit_label = ctk.CTkLabel(tension_frame, text="dyn/cm", width=50, anchor="w", font=FONTS['body'])
        tension_unit_label.pack(side="left", padx=5)
        stage_widgets["surface_tension"] = tension_var
        stage_widgets["surface_tension_entry"] = tension_entry
        stage_widgets["surface_tension_label"] = tension_text_label
        stage_widgets["surface_tension_unit_label"] = tension_unit_label
        stage_widgets["tension_frame"] = tension_frame
        if stage_data["ensemble"] not in ["NPAT", "NPgT"]:
            tension_frame.pack_forget()

        # Note: Ensemble callback removed since ensemble is now fixed per scheme
        # Pressure and tension visibility is set based on initial ensemble value from protocol
        
        # Timestep
        timestep_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        timestep_frame.pack(fill="x", pady=2)
        
        timestep_text_label = ctk.CTkLabel(timestep_frame, text="Timestep:", width=120, anchor="w", font=FONTS['body'])
        timestep_text_label.pack(side="left")
        timestep_var = ctk.StringVar(value=str(stage_data["timestep"]))
        timestep_entry = ctk.CTkEntry(timestep_frame, textvariable=timestep_var, width=100)
        timestep_entry.pack(side="left", padx=5)
        timestep_unit_label = ctk.CTkLabel(timestep_frame, text="fs", width=50, anchor="w", font=FONTS['body'])
        timestep_unit_label.pack(side="left", padx=5)
        stage_widgets["timestep"] = timestep_var
        stage_widgets["timestep_entry"] = timestep_entry
        stage_widgets["timestep_label"] = timestep_text_label
        stage_widgets["timestep_unit_label"] = timestep_unit_label
        
        # Add callback to update steps when timestep changes (reuse the update functions)
        def update_both_from_timestep(*args):
            update_steps_from_time()
        
        timestep_var.trace("w", update_both_from_timestep)
        
        # DCD frequency
        dcd_freq_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        dcd_freq_frame.pack(fill="x", pady=2)
        
        dcd_freq_text_label = ctk.CTkLabel(dcd_freq_frame, text="DCD Frequency:", width=120, anchor="w", font=FONTS['body'])
        dcd_freq_text_label.pack(side="left")
        dcd_freq_var = ctk.StringVar(value=str(stage_data.get("dcd_freq", 5000)))
        dcd_freq_entry = ctk.CTkEntry(dcd_freq_frame, textvariable=dcd_freq_var, width=100)
        dcd_freq_entry.pack(side="left", padx=5)
        dcd_freq_unit_label = ctk.CTkLabel(dcd_freq_frame, text="steps", width=50, anchor="w", font=FONTS['body'])
        dcd_freq_unit_label.pack(side="left", padx=5)
        stage_widgets["dcd_freq"] = dcd_freq_var
        stage_widgets["dcd_freq_entry"] = dcd_freq_entry
        stage_widgets["dcd_freq_label"] = dcd_freq_text_label
        stage_widgets["dcd_freq_unit_label"] = dcd_freq_unit_label
        stage_widgets["dcd_freq_frame"] = dcd_freq_frame
        
        # Margin parameter (for NPAT simulations)
        margin_frame = ctk.CTkFrame(params_frame, fg_color="transparent")
        margin_frame.pack(fill="x", pady=2)
        
        margin_text_label = ctk.CTkLabel(margin_frame, text="Margin:", width=120, anchor="w", font=FONTS['body'])
        margin_text_label.pack(side="left")
        margin_var = ctk.StringVar(value=str(stage_data.get("margin", 5.0)))
        margin_entry = ctk.CTkEntry(margin_frame, textvariable=margin_var, width=100)
        margin_entry.pack(side="left", padx=5)
        margin_unit_label = ctk.CTkLabel(margin_frame, text="Å", width=50, anchor="w", font=FONTS['body'])
        margin_unit_label.pack(side="left", padx=5)
        stage_widgets["margin"] = margin_var
        stage_widgets["margin_entry"] = margin_entry
        stage_widgets["margin_label"] = margin_text_label
        stage_widgets["margin_unit_label"] = margin_unit_label
        stage_widgets["margin_frame"] = margin_frame
        
        # Computational Resources section
        resources_label = ctk.CTkLabel(
            parent_frame,
            text="Computational Resources:",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        resources_label.pack(anchor="w", padx=10, pady=(10, 5))
        stage_widgets["resources_label"] = resources_label
        
        resources_frame = ctk.CTkFrame(parent_frame, fg_color="transparent")
        resources_frame.pack(fill="x", padx=10, pady=5)
        
        # GPU usage
        gpu_frame = ctk.CTkFrame(resources_frame, fg_color="transparent")
        gpu_frame.pack(fill="x", pady=2)
        
        gpu_text_label = ctk.CTkLabel(gpu_frame, text="Use GPU:", width=120, anchor="w", font=FONTS['body'])
        gpu_text_label.pack(side="left")
        gpu_var = ctk.BooleanVar(value=stage_data.get("use_gpu", True))
        gpu_checkbox = ctk.CTkCheckBox(gpu_frame, text="Enable GPU acceleration", variable=gpu_var,
                                     command=lambda sk=stage_key: self._on_gpu_toggle(sk))
        gpu_checkbox.pack(side="left", padx=5)
        stage_widgets["use_gpu"] = gpu_var
        stage_widgets["gpu_checkbox"] = gpu_checkbox
        stage_widgets["gpu_label"] = gpu_text_label
        
        # CPU cores
        cpu_frame = ctk.CTkFrame(resources_frame, fg_color="transparent")
        cpu_frame.pack(fill="x", pady=2)
        
        cpu_text_label = ctk.CTkLabel(cpu_frame, text="CPU Cores:", width=120, anchor="w", font=FONTS['body'])
        cpu_text_label.pack(side="left")
        cpu_var = ctk.StringVar(value=str(stage_data.get("cpu_cores", 1)))
        cpu_entry = ctk.CTkEntry(cpu_frame, textvariable=cpu_var, width=100)
        cpu_entry.pack(side="left", padx=5)
        cpu_unit_label = ctk.CTkLabel(cpu_frame, text="cores", width=50, anchor="w", font=FONTS['body'])
        cpu_unit_label.pack(side="left", padx=5)
        stage_widgets["cpu_cores"] = cpu_var
        stage_widgets["cpu_cores_entry"] = cpu_entry
        stage_widgets["cpu_label"] = cpu_text_label
        stage_widgets["cpu_unit_label"] = cpu_unit_label
        
        # GPU ID
        gpu_id_frame = ctk.CTkFrame(resources_frame, fg_color="transparent")
        gpu_id_frame.pack(fill="x", pady=2)
        
        gpu_id_text_label = ctk.CTkLabel(gpu_id_frame, text="GPU ID:", width=120, anchor="w", font=FONTS['body'])
        gpu_id_text_label.pack(side="left")
        gpu_id_var = ctk.StringVar(value=str(stage_data.get("gpu_id", 0)))
        gpu_id_entry = ctk.CTkEntry(gpu_id_frame, textvariable=gpu_id_var, width=100)
        gpu_id_entry.pack(side="left", padx=5)
        gpu_id_hint_label = ctk.CTkLabel(gpu_id_frame, text="(device ID)", width=50, anchor="w", font=FONTS['body'])
        gpu_id_hint_label.pack(side="left", padx=5)
        stage_widgets["gpu_id"] = gpu_id_var
        stage_widgets["gpu_id_entry"] = gpu_id_entry  # Store widget reference for enable/disable
        stage_widgets["gpu_id_label"] = gpu_id_text_label
        stage_widgets["gpu_id_hint_label"] = gpu_id_hint_label
        
        # Number of GPUs
        num_gpus_frame = ctk.CTkFrame(resources_frame, fg_color="transparent")
        num_gpus_frame.pack(fill="x", pady=2)
        
        num_gpus_text_label = ctk.CTkLabel(num_gpus_frame, text="Number of GPUs:", width=120, anchor="w", font=FONTS['body'])
        num_gpus_text_label.pack(side="left")
        num_gpus_var = ctk.StringVar(value=str(stage_data.get("num_gpus", 1)))
        num_gpus_entry = ctk.CTkEntry(num_gpus_frame, textvariable=num_gpus_var, width=100)
        num_gpus_entry.pack(side="left", padx=5)
        num_gpus_unit_label = ctk.CTkLabel(num_gpus_frame, text="GPUs", width=50, anchor="w", font=FONTS['body'])
        num_gpus_unit_label.pack(side="left", padx=5)
        stage_widgets["num_gpus"] = num_gpus_var
        stage_widgets["num_gpus_entry"] = num_gpus_entry  # Store widget reference for enable/disable
        stage_widgets["num_gpus_label"] = num_gpus_text_label
        stage_widgets["num_gpus_unit_label"] = num_gpus_unit_label
        
        # Set initial state of GPU widgets based on GPU usage
        gpu_enabled = stage_data.get("use_gpu", True)
        gpu_id_entry.configure(state="normal" if gpu_enabled else "disabled")
        num_gpus_entry.configure(state="normal" if gpu_enabled else "disabled")
        
        # Constraints section
        constraints_label = ctk.CTkLabel(
            parent_frame,
            text="Restraint Forces (kcal/mol/Å²):",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        constraints_label.pack(anchor="w", padx=10, pady=(10, 5))
        stage_widgets["constraints_label"] = constraints_label
        
        constraints_frame = ctk.CTkFrame(parent_frame, fg_color="transparent")
        constraints_frame.pack(fill="x", padx=10, pady=5)
        
        stage_widgets["constraints"] = {}
        stage_widgets["constraint_labels"] = []
        
        for constraint_name, constraint_value in stage_data["constraints"].items():
            const_frame = ctk.CTkFrame(constraints_frame, fg_color="transparent")
            const_frame.pack(fill="x", pady=1)
            
            display_name = constraint_name.replace("_", " ").title()
            const_name_label = ctk.CTkLabel(const_frame, text=f"{display_name}:", width=150, anchor="w")
            const_name_label.pack(side="left")
            
            const_var = ctk.StringVar(value=str(constraint_value))
            const_entry = ctk.CTkEntry(const_frame, textvariable=const_var, width=100)
            const_entry.pack(side="left", padx=5)
            
            stage_widgets["constraints"][constraint_name] = const_var
            # Store entry and label for fonts
            if "constraint_entries" not in stage_widgets:
                stage_widgets["constraint_entries"] = []
            stage_widgets["constraint_entries"].append(const_entry)
            stage_widgets["constraint_labels"].append(const_name_label)
        
        self.stage_widgets[stage_key] = stage_widgets
    
    def _create_action_section(self):
        """Create the action buttons section."""
        self.action_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.action_label = ctk.CTkLabel(
            self.action_section,
            text="Actions",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.action_buttons_frame = ctk.CTkFrame(self.action_section, fg_color="transparent")
        
        self.generate_files_btn = ctk.CTkButton(
            self.action_buttons_frame,
            text="Generate Input Files",
            command=self._generate_input_files,
            width=150
        )
        
        self.run_equilibration_btn = ctk.CTkButton(
            self.action_buttons_frame,
            text="Run Equilibration",
            command=self._run_equilibration,
            width=150
        )
        
        self.save_protocol_btn = ctk.CTkButton(
            self.action_buttons_frame,
            text="Save Protocol",
            command=self._save_protocol,
            width=150
        )
        
        self.load_protocol_btn = ctk.CTkButton(
            self.action_buttons_frame,
            text="Load Protocol",
            command=self._load_protocol,
            width=150
        )
    
    def _create_progress_section(self):
        """Create the progress tracking section for equilibration."""
        self.progress_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])

        self.progress_label = ctk.CTkLabel(
            self.progress_section,
            text="Equilibration Progress",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )

        # Create monitoring controls frame
        self.monitoring_controls = ctk.CTkFrame(self.progress_section, fg_color="transparent")

        # Auto-monitoring toggle
        self.auto_monitoring_var = ctk.BooleanVar(value=False)
        self.auto_monitoring_checkbox = ctk.CTkCheckBox(
            self.monitoring_controls,
            text="Auto-monitor",
            variable=self.auto_monitoring_var,
            command=self._toggle_auto_monitoring
        )

        # Update interval
        self.interval_label = ctk.CTkLabel(
            self.monitoring_controls,
            text="Update every:",
            font=FONTS['body']
        )

        self.interval_var = ctk.StringVar(value="5")
        self.interval_entry = ctk.CTkEntry(
            self.monitoring_controls,
            textvariable=self.interval_var,
            width=50
        )

        self.interval_unit_label = ctk.CTkLabel(
            self.monitoring_controls,
            text="seconds",
            font=FONTS['body']
        )

        # Manual refresh button
        self.refresh_progress_btn = ctk.CTkButton(
            self.monitoring_controls,
            text="Refresh Now",
            command=self._refresh_equilibration_progress,
            width=120
        )
        
        # Background process info button
        self.background_info_btn = ctk.CTkButton(
            self.monitoring_controls,
            text="Process Info",
            command=self._show_background_process_info,
            width=120
        )

        # Create horizontal progress display frame
        self.horizontal_progress_frame = ctk.CTkFrame(self.progress_section, fg_color="transparent")

        # Create stage progress bars and labels
        self.progress_widgets = {}
        self.stage_names = [
            'Equilibration 1',
            'Equilibration 2',
            'Equilibration 3',
            'Equilibration 4',
            'Equilibration 5',
            'Equilibration 6',
            'Production',
        ]

        display_map = {
            'Equilibration 1': 'Equilibration 1',
            'Equilibration 2': 'Equilibration 2',
            'Equilibration 3': 'Equilibration 3',
            'Equilibration 4': 'Equilibration 4',
            'Equilibration 5': 'Equilibration 5',
            'Equilibration 6': 'Equilibration 6',
            'Production': 'Production',
        }

        for stage in self.stage_names:
            stage_frame = ctk.CTkFrame(self.horizontal_progress_frame, fg_color=COLOR_SCHEME['content_inside_bg'])

            # Stage name label
            display_text = display_map.get(stage, stage.replace('_', ' ').title())
            stage_label = ctk.CTkLabel(
                stage_frame,
                text=display_text,
                font=FONTS['body'],
                text_color=COLOR_SCHEME['text']
            )

            # Progress bar
            progress_bar = ctk.CTkProgressBar(
                stage_frame,
                width=150,
                height=15
            )
            progress_bar.set(0)

            # Status label
            status_label = ctk.CTkLabel(
                stage_frame,
                text="Not Started",
                font=FONTS['small'],
                text_color=COLOR_SCHEME['inactive']
            )

            # Steps info label
            steps_label = ctk.CTkLabel(
                stage_frame,
                text="0/0 steps",
                font=FONTS['small'],
                text_color=COLOR_SCHEME['inactive']
            )

            # Pack widgets in stage frame
            stage_label.pack(pady=(5, 2))
            progress_bar.pack(pady=2)
            status_label.pack(pady=1)
            steps_label.pack(pady=(1, 5))

            # Store widgets for updates
            self.progress_widgets[stage] = {
                'frame': stage_frame,
                'label': stage_label,
                'progress_bar': progress_bar,
                'status_label': status_label,
                'steps_label': steps_label,
            }

        # Create summary information frame
        self.summary_frame = ctk.CTkFrame(self.progress_section, fg_color=COLOR_SCHEME['content_inside_bg'])

        self.summary_label = ctk.CTkLabel(
            self.summary_frame,
            text="System Information",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )

        # Create summary text display
        self.summary_text = ctk.CTkTextbox(
            self.summary_frame,
            height=100,
            font=FONTS['code'],
            text_color=COLOR_SCHEME['text']
        )
        self.summary_text.insert("1.0", "No simulation data available yet...")
        self.summary_text.configure(state="disabled")
        
        # Note: Mouse wheel scrolling is handled by the global _bind_mouse_wheel method

        # Timer for progress monitoring
        self.progress_timer = None
        self.monitoring_active = False
    
    def _setup_layout(self):
        """Setup the layout of all widgets."""
        # Pack main scroll frame
        self.main_scroll.pack(fill="both", expand=True, padx=10, pady=10)

        # Title section
        self.title_frame.pack(fill="x", padx=10, pady=(0, 10))
        self.title_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.description_label.pack(anchor="w", padx=15, pady=(0, 10))

        # Working directory section
        self.workdir_section.pack(fill="x", padx=10, pady=(0, 10))
        self.workdir_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.workdir_frame.pack(fill="x", padx=15, pady=(0, 10))
        self.workdir_entry.pack(side="left", fill="x", expand=True, padx=(0, 5))
        self.workdir_browse_btn.pack(side="right")

        # Input folder section
        self.inputfolder_section.pack(fill="x", padx=10, pady=(0, 10))
        self.inputfolder_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.inputfolder_sublabel.pack(anchor="w", padx=15, pady=(0, 5))
        self.inputfolder_frame.pack(fill="x", padx=15, pady=(0, 10))
        self.inputfolder_entry.pack(side="left", fill="x", expand=True, padx=(0, 5))
        self.inputfolder_browse_btn.pack(side="right")

        # Output name section
        self.outputname_section.pack(fill="x", padx=10, pady=(0, 10))
        self.outputname_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.outputname_frame.pack(fill="x", padx=15, pady=(0, 10))
        self.outputname_entry.pack(side="left", padx=(0, 5))

        # Engine section
        self.engine_section.pack(fill="x", padx=10, pady=(0, 10))
        self.engine_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.engine_menu.pack(anchor="w", padx=15, pady=(0, 5))
        self.engine_settings_frame.pack(fill="x", padx=15, pady=(0, 10))

        # Protocol section
        self.protocol_section.pack(fill="x", padx=10, pady=(0, 10))
        self.protocol_label.pack(anchor="w", padx=15, pady=(10, 5))
        
        # CHARMM-GUI scheme selection (always enabled)
        self.charmm_gui_scheme_frame.pack(fill="x", padx=15, pady=(0, 5))
        
        # Scheme type selection
        self.scheme_type_frame.pack(fill="x", pady=(5, 0))
        self.scheme_type_label.pack(side="left", padx=(0, 10))
        self.scheme_type_menu.pack(side="left")
        
        # Bilayer thickness
        #self.thickness_frame.pack(fill="x", pady=(5, 0))
        #self.thickness_label.pack(side="left", padx=(0, 10))
        #self.thickness_entry.pack(side="left")
        
        self.stages_config_frame.pack(fill="both", expand=True, padx=15, pady=(5, 10))

        # Action section
        self.action_section.pack(fill="x", padx=10, pady=(0, 10))
        self.action_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.action_buttons_frame.pack(fill="x", padx=15, pady=(0, 10))
        self.generate_files_btn.pack(side="left", padx=(0, 5))
        self.run_equilibration_btn.pack(side="left", padx=5)
        self.save_protocol_btn.pack(side="left", padx=5)
        self.load_protocol_btn.pack(side="left", padx=5)

        # Progress section
        self.progress_section.pack(fill="x", padx=10, pady=(0, 10))
        self.progress_label.pack(anchor="w", padx=15, pady=(10, 5))

        # Monitoring controls
        self.monitoring_controls.pack(fill="x", padx=15, pady=(0, 10))
        self.auto_monitoring_checkbox.pack(side="left", padx=(0, 10))
        self.interval_label.pack(side="left", padx=(0, 5))
        self.interval_entry.pack(side="left", padx=(0, 5))
        self.interval_unit_label.pack(side="left", padx=(0, 10))
        self.refresh_progress_btn.pack(side="left", padx=(0, 5))
        self.background_info_btn.pack(side="left")

        # Horizontal progress bars
        self.horizontal_progress_frame.pack(fill="x", padx=15, pady=(0, 10))

        # Pack stage frames in grid layout (3 rows, 3 columns for 7 stages)
        for i, stage in enumerate(self.stage_names):
            row = i // 3
            col = i % 3
            self.progress_widgets[stage]['frame'].grid(row=row, column=col, padx=5, pady=5, sticky="ew")

        # Configure grid weights for equal column spacing  
        for col in range(3):
            self.horizontal_progress_frame.grid_columnconfigure(col, weight=1)

        # Summary information
        self.summary_frame.pack(fill="x", padx=15, pady=(0, 10))
        self.summary_label.pack(anchor="w", padx=15, pady=(10, 5))
        self.summary_text.pack(fill="x", padx=15, pady=(0, 10))
    
    def _load_defaults(self):
        """Load default values."""
        self.workdir_entry.delete(0, tk.END)
        self.workdir_entry.insert(0, str(self.working_directory))
        
        # Initialize with CHARMM-GUI protocol (NPT scheme by default)
        self._load_charmm_gui_protocol()
    
    def _browse_workdir(self):
        """Browse for working directory."""
        directory = filedialog.askdirectory(
            title="Select Working Directory",
            initialdir=str(self.working_directory) if self.working_directory else self.initial_directory
        )
        
        if directory:
            self.working_directory = Path(directory)
            set_working_directory(str(self.working_directory))
            self.workdir_entry.delete(0, tk.END)
            self.workdir_entry.insert(0, str(self.working_directory))
            
            # Refresh equilibration progress for new directory
            self._refresh_equilibration_progress()
            
            if self.status_callback:
                self.status_callback(f"Working directory: {self.working_directory}")
    
    def _browse_inputfolder(self):
        """Browse for input folder containing prepared system files."""
        # Use working directory or initial directory as starting point
        initial_dir = str(self.working_directory) if self.working_directory else self.initial_directory
        
        directory = filedialog.askdirectory(
            title="Select Input Folder (Prepared System Files)",
            initialdir=initial_dir
        )
        
        if directory:
            input_folder = Path(directory)
            self.inputfolder_entry.delete(0, tk.END)
            self.inputfolder_entry.insert(0, str(input_folder))
            
            # Validate that required files exist
            required_files = ['system.prmtop', 'system.inpcrd', 'system.pdb']
            missing_files = []
            
            for file in required_files:
                if not (input_folder / file).exists():
                    missing_files.append(file)
            
            if missing_files:
                messagebox.showwarning(
                    "Missing Files",
                    f"Warning: The selected folder is missing the following files:\n\n" +
                    "\n".join(f"  • {f}" for f in missing_files) +
                    f"\n\nPlease ensure the folder contains all required system files:\n" +
                    "\n".join(f"  • {f}" for f in required_files) +
                    "\n\nEquilibration may fail if these files are not present."
                )
            else:
                if self.status_callback:
                    self.status_callback(f"✓ Input folder validated: {input_folder.name}")
    
    def _on_outputname_changed(self, event=None):
        """Handle output name change."""
        new_name = self.outputname_entry.get().strip()
        if new_name and new_name != self.equilibration_output_name:
            self.equilibration_output_name = new_name
            # Refresh progress monitoring with new folder name
            if self.monitoring_active and hasattr(self, '_refresh_equilibration_progress'):
                self._refresh_equilibration_progress()
            
            if self.status_callback:
                self.status_callback(f"Output folder name set to: {new_name}")
    
    def _on_engine_changed(self, engine: str):
        """Handle MD engine selection change."""
        self.logger.info(f"MD engine changed to: {engine}")
        
        # Clear existing engine settings
        for widget in self.engine_settings_frame.winfo_children():
            widget.destroy()
        
        # Add engine-specific settings
        if engine == "NAMD":
            self._create_namd_settings()
        elif engine == "GROMACS":
            self._create_gromacs_settings()
        elif engine == "AMBER":
            self._create_amber_settings()
        
        if self.status_callback:
            self.status_callback(f"Selected MD engine: {engine}")
    
    def _create_namd_settings(self):
        """Create NAMD-specific settings."""
        namd_label = ctk.CTkLabel(
            self.engine_settings_frame,
            text="NAMD Configuration:",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        namd_label.pack(anchor="w", pady=5)
        
        # NAMD executable path
        namd_path_frame = ctk.CTkFrame(self.engine_settings_frame, fg_color="transparent")
        namd_path_frame.pack(fill="x", pady=2)
        
        ctk.CTkLabel(namd_path_frame, text="NAMD executable:", width=120, anchor="w", font=FONTS['body']).pack(side="left")
        self.namd_path_var = ctk.StringVar(value="namd3")
        namd_path_entry = ctk.CTkEntry(namd_path_frame, textvariable=self.namd_path_var, width=200)
        namd_path_entry.pack(side="left", padx=5)
    
    def _create_gromacs_settings(self):
        """Create GROMACS-specific settings."""
        # Placeholder for GROMACS settings
        gromacs_label = ctk.CTkLabel(
            self.engine_settings_frame,
            text="GROMACS settings (To be implemented)",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        gromacs_label.pack(anchor="w", pady=5)
    
    def _create_amber_settings(self):
        """Create AMBER-specific settings."""
        # Placeholder for AMBER settings
        amber_label = ctk.CTkLabel(
            self.engine_settings_frame,
            text="AMBER settings (To be implemented)",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['text']
        )
        amber_label.pack(anchor="w", pady=5)
    
    def _on_scheme_type_changed(self, value):
        """Handle scheme type change (NVT, NPT, NPAT, NPgT)."""
        # Regenerate widgets for the new scheme
        self._load_charmm_gui_protocol()
            
        if self.status_callback:
            self.status_callback(f"Scheme changed to {value}")
    
    def _load_charmm_gui_protocol(self):
        """Load the CHARMM-GUI protocol based on selected scheme."""
        scheme_type = self.scheme_type_var.get()
        
        # Load the CHARMM-GUI protocol configuration
        charmm_gui_protocols = self._get_charmm_gui_protocols(scheme_type)
        
        # Regenerate stage widgets for the new protocol
        self._create_stage_widgets()
        
        # Update the protocol widgets with new values
        self._update_protocol_widgets(charmm_gui_protocols)
        
        if self.status_callback:
            self.status_callback(f"Loaded CHARMM-GUI {scheme_type} protocol with 13 equilibration stages")
    
    def _on_gpu_toggle(self, stage_key: str):
        """Handle GPU enable/disable toggle for a specific stage."""
        if stage_key not in self.stage_widgets:
            return
        
        stage_widgets = self.stage_widgets[stage_key]
        gpu_enabled = stage_widgets["use_gpu"].get()
        
        # Enable/disable GPU ID entry
        if "gpu_id_entry" in stage_widgets:
            stage_widgets["gpu_id_entry"].configure(state="normal" if gpu_enabled else "disabled")
        
        # Enable/disable Number of GPUs entry
        if "num_gpus_entry" in stage_widgets:
            stage_widgets["num_gpus_entry"].configure(state="normal" if gpu_enabled else "disabled")
        
        # Update status message
        if self.status_callback:
            stage_name = self.default_protocols[stage_key]["name"]
            status = "enabled" if gpu_enabled else "disabled"
            self.status_callback(f"GPU settings {status} for {stage_name}")
    
    def _refresh_equilibration_progress(self):
        """Refresh equilibration progress from log files."""
        try:
            # Update working directory from entry field before checking progress
            self._update_working_directory_from_entry()
            
            # Get progress information - look in equilibration/namd subdirectory
            from pathlib import Path
            equilibration_dir = Path(self.working_directory) / self.equilibration_output_name
            progress_dict = get_equilibration_progress(equilibration_dir)
            
            # Update horizontal progress bars
            self._update_horizontal_progress(progress_dict)
            
            # Update summary information
            self._update_summary_info(progress_dict)
            
            if self.status_callback:
                self.status_callback("Progress updated")
                
        except Exception as e:
            self.logger.debug(f" Error refreshing progress: {e}")
            # Update display to show error
            for stage in self.stage_names:
                if stage in self.progress_widgets:
                    self.progress_widgets[stage]['status_label'].configure(text="Error")
                    self.progress_widgets[stage]['progress_bar'].set(0)
            
            self.summary_text.configure(state="normal")
            self.summary_text.delete("1.0", "end")
            self.summary_text.insert("1.0", f"Error updating progress: {str(e)}")
            self.summary_text.configure(state="disabled")
    
    def _get_progress_stage_name(self, display_name: str) -> str:
        """
        Convert GUI display stage names to internal names expected by progress monitoring.
        
        Args:
            display_name: Display name from GUI (e.g., "Equilibration 1")
            
        Returns:
            Internal name for progress monitoring (e.g., "equilibration_1")
        """
        # Map display names to internal names used by namd_analysis.py (6 stages + production)
        display_to_internal = {
            "Equilibration 1": "equilibration_1",
            "Equilibration 2": "equilibration_2",
            "Equilibration 3": "equilibration_3",
            "Equilibration 4": "equilibration_4",
            "Equilibration 5": "equilibration_5",
            "Equilibration 6": "equilibration_6",
            "Production": "production"
        }
        
        return display_to_internal.get(display_name, display_name.lower().replace(" ", "_"))
    
    def _update_horizontal_progress(self, progress_dict):
        """Update the horizontal progress bars for each stage."""
        for stage_name in self.stage_names:
            if stage_name not in self.progress_widgets:
                continue
                
            widgets = self.progress_widgets[stage_name]
            
            # Convert display name to internal name for progress lookup
            internal_name = self._get_progress_stage_name(stage_name)
            
            if internal_name in progress_dict:
                progress = progress_dict[internal_name]
                
                # Update progress bar
                progress_percent = progress.progress_percent / 100.0 if progress.progress_percent else 0.0
                widgets['progress_bar'].set(progress_percent)
                
                # Update status with color coding
                status_text = progress.status.replace('_', ' ').title()
                if progress.status == 'completed':
                    widgets['status_label'].configure(text=status_text, text_color="#4CAF50")  # Green
                    widgets['progress_bar'].configure(progress_color="#4CAF50")
                elif progress.status == 'running':
                    widgets['status_label'].configure(text=status_text, text_color="#FF9800")  # Orange
                    widgets['progress_bar'].configure(progress_color="#FF9800")
                elif progress.status == 'error':
                    widgets['status_label'].configure(text=status_text, text_color="#F44336")  # Red
                    widgets['progress_bar'].configure(progress_color="#F44336")
                else:  # not_started
                    widgets['status_label'].configure(text=status_text, text_color=COLOR_SCHEME['inactive'])
                    widgets['progress_bar'].configure(progress_color=COLOR_SCHEME['highlight'])
                
                # Update steps info with nanosecond information
                if progress.timing:
                    steps_text = f"{progress.timing.steps_completed}/{progress.timing.total_steps} steps"
                    
                    # Add simulated time information if available
                    if progress.timing.simulated_time_ns > 0:
                        # Get target ns for this stage from protocol
                        target_ns = self._get_stage_target_ns(stage_name)
                        if target_ns > 0:
                            steps_text += f" ({progress.timing.simulated_time_ns:.3f}/{target_ns:.3f} ns)"
                        else:
                            steps_text += f" ({progress.timing.simulated_time_ns:.3f} ns)"
                    
                    if progress.timing.ns_per_day > 0:
                        steps_text += f" - {progress.timing.ns_per_day:.4f} ns/day"
                else:
                    # Show target nanoseconds even without timing data
                    target_ns = self._get_stage_target_ns(stage_name)
                    if target_ns > 0:
                        steps_text = f"0/0 steps (target: {target_ns:.3f} ns)"
                    else:
                        steps_text = "0/0 steps"
                
                widgets['steps_label'].configure(text=steps_text)
            else:
                # No data for this stage
                widgets['progress_bar'].set(0)
                widgets['status_label'].configure(text="No Data", text_color=COLOR_SCHEME['inactive'])
                # Show target nanoseconds even for stages with no data
                target_ns = self._get_stage_target_ns(stage_name)
                if target_ns > 0:
                    widgets['steps_label'].configure(text=f"0/0 steps (target: {target_ns:.3f} ns)")
                else:
                    widgets['steps_label'].configure(text="0/0 steps")
    
    def _get_stage_target_ns(self, stage_name: str) -> float:
        """Get the target nanoseconds for a given stage from the current GUI configuration or config files."""
        try:
            # First try to get from current stage widgets (user input)
            if hasattr(self, 'stage_widgets') and stage_name in self.stage_widgets:
                stage_widgets = self.stage_widgets[stage_name]
                if 'time_ns' in stage_widgets:
                    time_ns_str = stage_widgets['time_ns'].get().strip()
                    if time_ns_str:
                        time_ns = float(time_ns_str)
                        
                        # For Equilibration 1, add minimize time to the target
                        if stage_name == 'Equilibration 1':
                            # Get minimize steps and timestep from widgets
                            minimize_steps = 0
                            timestep_fs = 2.0  # default
                            
                            if 'minimize_steps' in stage_widgets:
                                minimize_steps_str = stage_widgets['minimize_steps'].get().strip()
                                if minimize_steps_str:
                                    minimize_steps = int(minimize_steps_str)
                            
                            if 'timestep' in stage_widgets:
                                timestep_str = stage_widgets['timestep'].get().strip()
                                if timestep_str:
                                    timestep_fs = float(timestep_str)
                            
                            # Add minimize time to total target time
                            minimize_time_ns = (minimize_steps * timestep_fs) / 1e6
                            time_ns += minimize_time_ns
                        
                        return time_ns
            
            # Try to read from generated configuration files
            target_ns_from_file = self._read_target_ns_from_config_file(stage_name)
            if target_ns_from_file > 0:
                return target_ns_from_file
            
            # Fallback to protocol template if no user input found
            current_protocol = self._get_current_protocol_template()
            if stage_name in current_protocol:
                time_ns = current_protocol[stage_name].get('time_ns', 0.0)
                
                # For Equilibration 1, add minimize time to the target
                if stage_name == 'Equilibration 1':
                    minimize_steps = current_protocol[stage_name].get('minimize_steps', 10000)
                    timestep_fs = current_protocol[stage_name].get('timestep', 2.0)
                    minimize_time_ns = (minimize_steps * timestep_fs) / 1e6
                    time_ns += minimize_time_ns
                
                return time_ns
                
        except Exception as e:
            self.logger.debug(f" Error getting target ns for {stage_name}: {e}")
        return 0.0
    
    def _read_target_ns_from_config_file(self, stage_name: str) -> float:
        """Read target nanoseconds from the NAMD configuration file."""
        try:
            # Get the working directory
            working_dir = Path(self.working_directory) if self.working_directory else Path.cwd()
            equilibration_dir = working_dir / self.equilibration_output_name / 'namd'
            
            # Map stage names to config file names
            config_file_mapping = {
                'equilibration_1': 'eq1_equilibration_1.conf',
                'equilibration_2': 'eq2_equilibration_2.conf',
                'equilibration_3': 'eq3_equilibration_3.conf',
                'equilibration_4': 'eq4_equilibration_4.conf',
                'equilibration_5': 'eq5_equilibration_5.conf',
                'equilibration_6': 'eq6_equilibration_6.conf',
                'production': 'step7_production.conf'
            }
            
            config_filename = config_file_mapping.get(stage_name)
            if not config_filename:
                return 0.0
                
            config_file = equilibration_dir / config_filename
            if not config_file.exists():
                return 0.0
            
            # Read the config file and extract timestep and run/minimize commands
            with open(config_file, 'r') as f:
                content = f.read()
            
            timestep_fs = 2.0  # Default timestep in fs
            total_steps = 0
            
            # Extract timestep (in fs)
            import re
            timestep_match = re.search(r'timestep\s+([0-9.]+)', content)
            if timestep_match:
                timestep_fs = float(timestep_match.group(1))
            
            # Extract minimize steps (only in first stage)
            minimize_matches = re.findall(r'minimize\s+([0-9]+)', content)
            for match in minimize_matches:
                total_steps += int(match)
            
            # Extract steps from direct run commands with numeric values
            run_matches = re.findall(r'run\s+([0-9]+)', content)
            for match in run_matches:
                total_steps += int(match)
            
            # Extract steps from TCL expressions like: run [expr int($time * 1e6 / $tstep)]
            # First extract the time_ns value from: set time {TIME_NS} or actual value
            time_ns = 0.0
            time_match = re.search(r'set time\s+([0-9.]+)', content)
            if time_match:
                time_ns = float(time_match.group(1))
            elif re.search(r'run\s+\[expr\s+int\(\$time\s*\*\s*1e6\s*/\s*\$tstep\)\]', content):
                # If we find the TCL expression but no set time, try to extract from comments or fallback
                # Look for comment with actual time value
                comment_time_match = re.search(r'#.*([0-9.]+)\s*ns', content)
                if comment_time_match:
                    time_ns = float(comment_time_match.group(1))
            
            # Calculate steps from time_ns if we found it
            if time_ns > 0:
                # Convert ns to fs, then divide by timestep_fs to get steps
                calculated_steps = int(time_ns * 1e6 / timestep_fs)
                total_steps += calculated_steps
            
            # Calculate total nanoseconds
            # timestep_fs is in fs, convert to ns: timestep_fs * steps / 1e6
            total_ns = (timestep_fs * total_steps) / 1e6
            
            return total_ns
            
        except Exception as e:
            self.logger.debug(f" Error reading target ns from config file for {stage_name}: {e}")
            return 0.0
    
    def _update_summary_info(self, progress_dict):
        """Update the summary information display."""
        summary_lines = []
        
        # Find any stage with timing information
        timing_info = None
        for progress in progress_dict.values():
            if progress.timing and progress.timing.atoms > 0:
                timing_info = progress.timing
                break
        
        if timing_info:
            summary_lines.append("=== SYSTEM INFORMATION ===")
            summary_lines.append(f"Atoms: {timing_info.atoms:,}")
            summary_lines.append(f"Processors: {timing_info.processors}")
            summary_lines.append(f"GPUs: {timing_info.gpus}")
            summary_lines.append(f"Timestep: {timing_info.timestep_fs:.2f} fs")
            summary_lines.append("")
            
            summary_lines.append("=== PERFORMANCE ===")
            if timing_info.ns_per_day > 0:
                summary_lines.append(f"Speed: {timing_info.ns_per_day:.2f} ns/day")
            if timing_info.sec_per_step > 0:
                summary_lines.append(f"Time per step: {timing_info.sec_per_step:.3f} s")
            if timing_info.real_time_hours > 0:
                summary_lines.append(f"Runtime: {timing_info.real_time_hours:.2f} hours")
            summary_lines.append("")
            
            if timing_info.hostname:
                summary_lines.append(f"Host: {timing_info.hostname}")
            summary_lines.append("")
        
        # Add stage summary
        summary_lines.append("=== STAGES SUMMARY ===")
        for stage_name in self.stage_names:
            # Convert display name to internal name for progress lookup
            internal_name = self._get_progress_stage_name(stage_name)
            
            if internal_name in progress_dict:
                progress = progress_dict[internal_name]
                status_icon = {
                    'completed': '[OK]',
                    'running': '[>>]',
                    'error': '[XX]',
                    'not_started': '[--]'
                }.get(progress.status, '[??]')
                
                # Use the display name directly (it's already properly formatted)
                summary_lines.append(f"{status_icon} {stage_name}: {progress.status.replace('_', ' ').title()}")
                
                if progress.timing and progress.timing.total_steps > 0:
                    percent = (progress.timing.steps_completed / progress.timing.total_steps) * 100
                    summary_lines.append(f"   Progress: {progress.timing.steps_completed}/{progress.timing.total_steps} steps ({percent:.1f}%)")
        
        if not summary_lines:
            summary_lines = ["No simulation data available yet...", 
                           "Make sure your working directory contains an 'equilibration/namd/' folder with NAMD log files."]
        
        # Update the summary text display
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("1.0", "\n".join(summary_lines))
        self.summary_text.configure(state="disabled")
    
    def _update_working_directory_from_entry(self):
        """Update working directory from the entry field."""
        workdir_text = self.workdir_entry.get().strip()
        if workdir_text and Path(workdir_text).exists():
            self.working_directory = Path(workdir_text)
            set_working_directory(str(self.working_directory))
    
    def _format_equilibration_summary(self, progress_dict):
        """Format equilibration progress for display."""
        lines = ["🧪 EQUILIBRATION PROGRESS MONITOR", "=" * 50, ""]
        
        # Overall summary
        total_stages = len(progress_dict)
        completed_stages = sum(1 for p in progress_dict.values() if p.status == "completed")
        running_stages = sum(1 for p in progress_dict.values() if p.status == "running")
        
        lines.append(f"📊 Overview: {completed_stages}/{total_stages} stages completed, {running_stages} running")
        lines.append("")
        
        # Detailed stage information
        for stage_name, progress in progress_dict.items():
            stage_display = stage_name.replace('_', ' ').title()
            
            # Status icon
            status_icons = {
                "not_started": "[--]",
                "running": "[>>]",
                "completed": "[OK]",
                "error": "[XX]"
            }
            icon = status_icons.get(progress.status, "[??]")
            
            lines.append(f"{icon} {stage_display}")
            lines.append(f"   Status: {progress.status.replace('_', ' ').title()}")
            
            if progress.timing:
                timing = progress.timing
                
                # Progress information
                if timing.total_steps > 0:
                    progress_pct = (timing.steps_completed / timing.total_steps) * 100
                    lines.append(f"   Progress: {timing.steps_completed:,}/{timing.total_steps:,} steps ({progress_pct:.1f}%)")
                elif timing.steps_completed > 0:
                    lines.append(f"   Steps completed: {timing.steps_completed:,}")
                
                # Performance information
                if timing.ns_per_day > 0:
                    lines.append(f"   Performance: {timing.ns_per_day:.2f} ns/day")
                
                if timing.simulated_time_ns > 0:
                    lines.append(f"   Simulated time: {timing.simulated_time_ns:.3f} ns")
                
                if timing.real_time_hours > 0:
                    lines.append(f"   Real time: {timing.real_time_hours:.2f} hours")
                
                # System information
                if timing.processors > 0:
                    gpu_info = f", {timing.gpus} GPUs" if timing.gpus > 0 else ""
                    lines.append(f"   💻 Resources: {timing.processors} CPUs{gpu_info}")
                
                if timing.atoms > 0:
                    lines.append(f"   🧬 System: {timing.atoms:,} atoms")
                
                if timing.hostname:
                    lines.append(f"   🌐 Host: {timing.hostname}")
            
            # Log file information
            if progress.log_file:
                from datetime import datetime
                mod_time = datetime.fromtimestamp(progress.last_updated)
                lines.append(f"   📄 Log: {progress.log_file.name} (updated: {mod_time.strftime('%H:%M:%S')})")
            else:
                lines.append(f"   📄 Log: No log file found")
            
            lines.append("")
        
        # Add timestamp
        from datetime import datetime
        lines.append(f"Last updated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        
        return "\n".join(lines)
    
    def _toggle_auto_refresh(self):
        """Toggle automatic progress monitoring (kept for compatibility)."""
        # Redirect to new monitoring system
        self._toggle_auto_monitoring()
    
    def _toggle_auto_monitoring(self):
        """Toggle automatic progress monitoring."""
        if self.auto_monitoring_var.get():
            self._start_progress_monitoring()
        else:
            self._stop_progress_monitoring()
    
    def _start_auto_refresh(self):
        """Start automatic progress refreshing (kept for compatibility)."""
        self._start_progress_monitoring()
    
    def _start_progress_monitoring(self):
        """Start automatic progress monitoring."""
        if not self.monitoring_active:
            self.monitoring_active = True
            self._schedule_next_refresh()
            if self.status_callback:
                self.status_callback("Progress monitoring started")
    
    def _stop_auto_refresh(self):
        """Stop automatic progress refreshing (kept for compatibility)."""
        self._stop_progress_monitoring()
    
    def _stop_progress_monitoring(self):
        """Stop automatic progress monitoring."""
        self.monitoring_active = False
        if self.progress_timer:
            self.after_cancel(self.progress_timer)
            self.progress_timer = None
        if self.status_callback:
            self.status_callback("Progress monitoring stopped")
    
    def _schedule_next_refresh(self):
        """Schedule the next progress refresh."""
        if self.monitoring_active:
            try:
                interval = int(self.interval_var.get()) * 1000  # Convert to milliseconds
            except ValueError:
                interval = 5000  # Default 5 seconds
            
            self.progress_timer = self.after(interval, self._auto_refresh_callback)
    
    def _auto_refresh_callback(self):
        """Callback for automatic refresh timer."""
        if self.monitoring_active:
            self._refresh_equilibration_progress()
            self._schedule_next_refresh()
    
    def _generate_input_files(self):
        """Generate MD input files."""
        # Update working directory from entry field
        workdir_text = self.workdir_entry.get().strip()
        if workdir_text:
            self.working_directory = Path(workdir_text)
        
        if not self._validate_inputs():
            return
        
        engine = self.engine_var.get()
        
        try:
            if engine == "NAMD":
                self._generate_namd_files()
            elif engine == "GROMACS":
                self._generate_gromacs_files()
            elif engine == "AMBER":
                self._generate_amber_files()
            
            messagebox.showinfo(
                "Success",
                f"Input files generated successfully for {engine}!"
            )
            
            if self.status_callback:
                self.status_callback(f"Generated {engine} input files")
        
        except Exception as e:
            self.logger.error(f"Error generating input files: {e}")
            messagebox.showerror("Error", f"Failed to generate input files: {str(e)}")
    
    def _generate_namd_files(self):
        """Generate NAMD input files for all equilibration stages."""
        output_dir = Path(self.working_directory) / self.equilibration_output_name / "namd"
        
        # Create output directory with robust error handling
        try:
            self.logger.info(f"Creating output directory: {output_dir}")
            
            # Enhanced directory creation with better error handling
            from gatewizard.utils.helpers import create_directory_robust
            create_directory_robust(output_dir)
            
            self.logger.info(f"Successfully created output directory: {output_dir}")
        except OSError as e:
            self.logger.error(f"Failed to create output directory: {e}")
            raise OSError(f"Cannot create output directory: {output_dir}. Error: {e}")
        
        # Copy system files to NAMD directory for self-contained execution
        self._copy_system_files_to_namd_dir(self.working_directory, output_dir)
        
        # Get current protocol parameters
        protocols = self._get_current_protocols()
        
        # Force field is AMBER only
        force_field = "amber"
        working_dir = Path(self.working_directory)
        
        # Check for AMBER topology and coordinate files
        self._check_amber_files(output_dir)
        
        # Initialize NAMD manager with working_dir (where system files are located)
        namd_exe = getattr(self, 'namd_path_var', ctk.StringVar(value="namd3")).get()
        namd_manager = NAMDEquilibrationManager(working_dir, namd_exe)
        
        # Check for AMBER files from preparation
        amber_files = list(working_dir.glob("**/*.top")) + list(working_dir.glob("**/*.inpcrd"))
        if amber_files:
            self.logger.info("Found AMBER files from preparation, using AMBER force field")
        else:
            self.logger.warning("No AMBER files found. Please ensure AMBER topology and coordinate files are available.")
        
        # System files for AMBER - look for .prmtop and .inpcrd files
        system_files = {
            'prmtop': 'system.prmtop',
            'inpcrd': 'system.inpcrd',
            'pdb': 'system.pdb'  # Optional PDB for visualization
        }
        
        # Try to find actual AMBER files in the working directory
        prmtop_files = list(working_dir.glob("**/*.prmtop")) + list(working_dir.glob("**/*.top"))
        inpcrd_files = list(working_dir.glob("**/*.inpcrd")) + list(working_dir.glob("**/*.rst"))
        
        if prmtop_files:
            system_files['prmtop'] = prmtop_files[0].name
            self.logger.info(f"Found topology file: {prmtop_files[0].name}")
        
        if inpcrd_files:
            system_files['inpcrd'] = inpcrd_files[0].name
            self.logger.info(f"Found coordinate file: {inpcrd_files[0].name}")
        
        # Always use CHARMM-GUI scheme
        scheme_type = self.scheme_type_var.get()
        
        # Generate configuration files for each stage (skip minimization)
        previous_stage_key = None
        for i, (stage_key, stage_data) in enumerate(protocols.items()):
            # Skip minimization - now incorporated into first equilibration
            if stage_key == "minimization":
                self.logger.info("Skipping separate minimization stage - now included in eq1_equilibration")
                continue
            
            # Always use CHARMM-GUI template generation
            # Force scheme_type for all stages (GUI enforces scheme-level template selection)
            namd_config = namd_manager.generate_charmm_gui_config_file(
                stage_key, stage_data, i, system_files, scheme_type, previous_stage_key, protocols,
                force_scheme_type=True
            )
            
            # Skip empty configs
            if not namd_config.strip():
                continue
            
            # Write configuration file using config-safe naming with new step convention
            config_name = namd_manager._get_config_name(stage_key, i)
            if config_name == "step7_production":
                config_file = output_dir / f"{config_name}.conf"
            else:
                config_file = output_dir / f"{config_name}_equilibration.conf"
            with open(config_file, 'w') as f:
                f.write(namd_config)
            
            self.logger.info(f"Generated NAMD config: {config_file}")
            
            # Update previous stage key for next iteration
            previous_stage_key = stage_key
        
        # Generate restraints files for each stage (if needed)
        # Use system.pdb from the NAMD output directory (where it was just copied)
        system_pdb = output_dir / "system.pdb"
        if not system_pdb.exists():
            # Fallback: try to find it in the working directory
            self.logger.warning(f"system.pdb not found in output directory: {output_dir}")
            self.logger.info("Searching for system.pdb in working directory...")
            system_pdb = self._find_system_pdb(working_dir)
        
        if system_pdb and system_pdb.exists():
            self.logger.info(f"Using PDB file for restraints generation: {system_pdb}")
            # Create restraints directory with error handling
            restraints_dir = output_dir / "restraints"
            try:
                self.logger.info(f"Creating restraints directory: {restraints_dir}")
                from gatewizard.utils.helpers import create_directory_robust
                create_directory_robust(restraints_dir)
                self.logger.info(f"Successfully created restraints directory: {restraints_dir}")
            except OSError as e:
                self.logger.error(f"Failed to create restraints directory {restraints_dir}: {e}")
                # If restraints directory creation fails, continue without restraints
                self.logger.warning("Continuing without restraints due to directory creation failure")
                restraints_dir = None
            
            # Generate restraints file for each stage with different constraints
            if restraints_dir is not None:
                for i, (stage_key, stage_data) in enumerate(protocols.items()):
                    stage_constraints = stage_data.get('constraints', {})
                    has_restraints = any(float(v) > 0 for v in stage_constraints.values())
                    
                    if has_restraints:
                        # Use config-safe naming for restraints files with new step convention
                        config_name = namd_manager._get_config_name(stage_key, i)
                        if config_name == "step7_production":
                            restraints_file = restraints_dir / f"{config_name}_restraints.pdb"
                        else:
                            restraints_file = restraints_dir / f"{config_name}_equilibration_restraints.pdb"
                        namd_manager.generate_restraints_file(
                            system_pdb, 
                            stage_constraints,
                            restraints_file,
                            stage_data.get('name', stage_key)
                        )
                        self.logger.info(f"Generated restraints for {stage_key}: {restraints_file}")
                    
                    # Note: ABF/colvars generation removed as requested
                
                # Also create a general restraints.pdb for backward compatibility
                # Use constraints from first stage for this
                first_stage = list(protocols.values())[0]
                general_restraints = output_dir / "restraints.pdb"
                namd_manager.generate_restraints_file(
                    system_pdb, 
                    first_stage.get('constraints', {}),
                    general_restraints,
                    "General"
                )
                
                self.logger.info(f"Generated general restraints file: {general_restraints}")
            else:
                self.logger.warning("Skipping restraints file generation due to directory creation failure")
        
        # Generate master run script
        run_script_content = namd_manager.generate_run_script(
            protocols, namd_exe
        )
        
        script_file = output_dir / "run_equilibration.sh"
        with open(script_file, 'w') as f:
            f.write(run_script_content)
        
        # Make script executable
        script_file.chmod(0o755)
        
        self.logger.info(f"Generated run script: {script_file}")
        
        # Create a protocol summary file
        protocol_summary = {
            "protocol_name": "AMBER Equilibration Protocol",
            "total_stages": len(protocols),
            "stages": protocols,
            "namd_executable": namd_exe,
            "force_field": force_field,
            "note": "Individual CPU/GPU settings per stage configured in protocols"
        }
        
        summary_file = output_dir / "protocol_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(protocol_summary, f, indent=2)
        
        self.logger.info(f"Generated protocol summary: {summary_file}")
    
    def _find_system_pdb(self, working_dir: Path) -> Optional[Path]:
        """Find system PDB file in working directory."""
        # First check for system.pdb in the working directory
        system_pdb = working_dir / "system.pdb"
        if system_pdb.exists():
            return system_pdb
        
        # Check for bilayer PDB files (membrane systems)
        bilayer_files = list(working_dir.glob("bilayer_*.pdb"))
        if bilayer_files:
            # Use the main bilayer file (not lipid-only)
            for pdb_file in bilayer_files:
                if "lipid" not in pdb_file.name:
                    return pdb_file
            # If only lipid files, use the first one
            return bilayer_files[0]
        
        # Look for PDB files in preparation directories
        prep_dirs = list(working_dir.glob("membrane_*"))
        for prep_dir in sorted(prep_dirs, key=lambda x: x.stat().st_mtime, reverse=True):
            pdb_files = list(prep_dir.glob("*.pdb"))
            if pdb_files:
                return pdb_files[0]
        
        # Look for any PDB file in the working directory
        # Prioritize system.pdb over other files
        pdb_files = list(working_dir.glob("*.pdb"))
        if pdb_files:
            # First, try to find system.pdb (case-insensitive)
            for pdb_file in pdb_files:
                if pdb_file.name.lower() == "system.pdb":
                    return pdb_file
            
            # If system.pdb not found, exclude protein.pdb and prefer full system files
            non_protein_files = [f for f in pdb_files if "protein" not in f.name.lower()]
            if non_protein_files:
                return non_protein_files[0]
            
            # Last resort: return any PDB file
            return pdb_files[0]
        
        return None
    
    def _find_original_bilayer_pdb(self, working_dir: Path) -> Optional[Path]:
        """
        Find the original bilayer PDB file (from packmol-memgen or similar).
        Prioritizes full system files over lipid-only files.
        
        Args:
            working_dir: Directory to search in
            
        Returns:
            Path to original bilayer PDB file, or None if not found
        """
        # Look for bilayer PDB files
        bilayer_patterns = [
            "bilayer_*.pdb",
            "*_bilayer.pdb"
        ]
        
        for pattern in bilayer_patterns:
            bilayer_files = list(working_dir.glob(pattern))
            for pdb_file in bilayer_files:
                # Skip lipid-only files, prefer full system  
                if "lipid" not in pdb_file.name.lower():
                    self.logger.info(f"Found original bilayer PDB: {pdb_file}")
                    return pdb_file
        
        # If no full system found, try lipid files
        for pattern in bilayer_patterns:
            bilayer_files = list(working_dir.glob(pattern))
            for pdb_file in bilayer_files:
                if "lipid" in pdb_file.name.lower():
                    self.logger.info(f"Using lipid bilayer PDB: {pdb_file}")
                    return pdb_file
        
        return None
    
    def _copy_fallback_files(self, source_dir: Path, output_dir: Path):
        """Copy available files as fallback when conversion fails."""
        import shutil
        
        # Copy PDB files for visualization (optional)
        pdb_files = list(source_dir.glob("*.pdb"))
        if pdb_files:
            shutil.copy2(pdb_files[0], output_dir / "system.pdb")
            self.logger.info(f"Copied PDB file for visualization: {pdb_files[0]}")
        
        # Copy AMBER topology and coordinate files
        top_files = list(source_dir.glob("*.top")) + list(source_dir.glob("*.prmtop"))
        if top_files:
            shutil.copy2(top_files[0], output_dir / "system.prmtop")
            self.logger.info(f"Copied topology file: {top_files[0]}")
        
        inpcrd_files = list(source_dir.glob("*.inpcrd")) + list(source_dir.glob("*.rst"))
        if inpcrd_files:
            shutil.copy2(inpcrd_files[0], output_dir / "system.inpcrd")
            self.logger.info(f"Copied coordinate file: {inpcrd_files[0]}")
    
    def _copy_system_files_to_namd_dir(self, working_dir: Path, namd_dir: Path):
        """Copy system topology and coordinate files to NAMD directory for self-contained execution."""
        import shutil
        
        # Check if user specified an input folder
        input_folder_str = self.inputfolder_entry.get().strip()
        
        if input_folder_str:
            # User specified an input folder - use it directly
            source_dir = Path(input_folder_str)
            self.logger.info(f"Using specified input folder: {source_dir}")
        else:
            # No input folder specified - search in working directory (old behavior)
            source_dir = working_dir
            self.logger.info(f"No input folder specified, searching in working directory: {source_dir}")
        
        # Validate source directory exists
        if not source_dir.exists():
            error_msg = f"Input folder does not exist: {source_dir}"
            self.logger.error(error_msg)
            messagebox.showerror(
                "Input Folder Not Found",
                f"{error_msg}\n\nPlease select a valid input folder containing:\n" +
                "  • system.prmtop\n  • system.inpcrd\n  • system.pdb"
            )
            raise FileNotFoundError(error_msg)
        
        self.logger.info(f"Copying system files from {source_dir} to {namd_dir}")
        
        # Find and copy topology file (.prmtop or .top)
        if input_folder_str:
            # Look directly in the specified folder
            topology_files = list(source_dir.glob("*.prmtop")) + list(source_dir.glob("*.top"))
        else:
            # Search recursively in working directory (old behavior for backward compatibility)
            topology_files = list(source_dir.glob("**/*.prmtop")) + list(source_dir.glob("**/*.top"))
        
        if topology_files:
            # Use the first found topology file
            source_topology = topology_files[0]
            target_topology = namd_dir / "system.prmtop"
            shutil.copy2(source_topology, target_topology)
            self.logger.info(f"Copied topology file: {source_topology} -> {target_topology}")
        else:
            error_msg = f"No topology file (.prmtop or .top) found in: {source_dir}"
            self.logger.error(error_msg)
            messagebox.showerror(
                "Missing Topology File",
                f"{error_msg}\n\nPlease ensure the input folder contains system.prmtop or system.top"
            )
            raise FileNotFoundError(error_msg)
        
        # Find and copy coordinate file (.inpcrd or .rst)
        if input_folder_str:
            # Look directly in the specified folder
            coordinate_files = list(source_dir.glob("*.inpcrd")) + list(source_dir.glob("*.rst"))
        else:
            # Search recursively in working directory (old behavior)
            coordinate_files = list(source_dir.glob("**/*.inpcrd")) + list(source_dir.glob("**/*.rst"))
        
        if coordinate_files:
            # Use the first found coordinate file
            source_coords = coordinate_files[0]
            target_coords = namd_dir / "system.inpcrd"
            shutil.copy2(source_coords, target_coords)
            self.logger.info(f"Copied coordinate file: {source_coords} -> {target_coords}")
        else:
            error_msg = f"No coordinate file (.inpcrd or .rst) found in: {source_dir}"
            self.logger.error(error_msg)
            messagebox.showerror(
                "Missing Coordinate File",
                f"{error_msg}\n\nPlease ensure the input folder contains system.inpcrd or system.rst"
            )
            raise FileNotFoundError(error_msg)
        
        # Optionally copy PDB file for reference/visualization
        if input_folder_str:
            # Look directly in the specified folder
            pdb_files = list(source_dir.glob("*.pdb"))
        else:
            # Search recursively in working directory (old behavior)
            pdb_files = list(source_dir.glob("**/*.pdb"))
        
        if pdb_files:
            # Try to find system.pdb first, otherwise use any PDB
            system_pdb = None
            for pdb_file in pdb_files:
                if pdb_file.name == "system.pdb":
                    system_pdb = pdb_file
                    break
            if not system_pdb:
                system_pdb = pdb_files[0]
            
            target_pdb = namd_dir / "system.pdb"
            shutil.copy2(system_pdb, target_pdb)
            self.logger.info(f"Copied PDB file: {system_pdb} -> {target_pdb}")
        else:
            self.logger.warning(f"No PDB file found in: {source_dir}")
    
    def _generate_gromacs_files(self):
        """Generate GROMACS input files (placeholder)."""
        messagebox.showinfo("Info", "GROMACS file generation not yet implemented")
    
    def _generate_amber_files(self):
        """Generate AMBER input files (placeholder)."""
        messagebox.showinfo("Info", "AMBER file generation not yet implemented")
    
    def _check_running_equilibration(self) -> bool:
        """Check if there's already a running equilibration process."""
        try:
            output_dir = Path(self.working_directory) / self.equilibration_output_name / "namd"
            pid_file = output_dir / "equilibration.pid"
            
            if not pid_file.exists():
                return False
            
            with open(pid_file, 'r') as f:
                pid = int(f.read().strip())
            
            # Check if process is still running
            try:
                import psutil  # type: ignore
                if psutil.pid_exists(pid):
                    process = psutil.Process(pid)
                    if process.is_running():
                        return True
                    else:
                        # Process is dead, remove PID file
                        pid_file.unlink()
                        return False
                else:
                    # PID doesn't exist, remove PID file
                    pid_file.unlink()
                    return False
            except ImportError:
                # psutil not available, use basic OS check
                import os
                import signal
                try:
                    os.kill(pid, 0)  # Send signal 0 to check if process exists
                    return True
                except OSError:
                    # Process doesn't exist, remove PID file
                    pid_file.unlink()
                    return False
                    
        except Exception as e:
            self.logger.warning(f"Error checking running equilibration: {e}")
            return False
    
    def _get_running_equilibration_info(self) -> Optional[Dict[str, Any]]:
        """Get information about the running equilibration process."""
        try:
            output_dir = Path(self.working_directory) / self.equilibration_output_name / "namd"
            pid_file = output_dir / "equilibration.pid"
            log_file = output_dir / "equilibration_background.log"
            
            if not pid_file.exists():
                return None
            
            with open(pid_file, 'r') as f:
                pid = int(f.read().strip())
            
            info = {
                'pid': pid,
                'pid_file': pid_file,
                'log_file': log_file,
                'running': False,
                'command': '',
                'start_time': '',
                'cpu_percent': 0.0,
                'memory_percent': 0.0
            }
            
            # Check if process is still running and get details
            try:
                import psutil  # type: ignore
                if psutil.pid_exists(pid):
                    process = psutil.Process(pid)
                    if process.is_running():
                        info['running'] = True
                        info['command'] = ' '.join(process.cmdline())
                        info['start_time'] = str(process.create_time())
                        info['cpu_percent'] = process.cpu_percent()
                        info['memory_percent'] = process.memory_percent()
            except ImportError:
                # psutil not available, use basic check
                import os
                try:
                    os.kill(pid, 0)
                    info['running'] = True
                except OSError:
                    pass
            
            return info
            
        except Exception as e:
            self.logger.warning(f"Error getting equilibration info: {e}")
            return None
    
    def _show_background_process_info(self):
        """Show information about running background equilibration process."""
        info = self._get_running_equilibration_info()
        
        if info is None:
            messagebox.showinfo(
                "Process Info",
                "No equilibration process information found.\n\n"
                "This could mean:\n"
                "• No equilibration is currently running\n"
                "• The process was started outside of this application\n"
                "• The PID file was removed or corrupted"
            )
            return
        
        if info['running']:
            try:
                import psutil  # type: ignore
                process_info = (
                    f"Equilibration Process Information\n"
                    f"{'=' * 40}\n\n"
                    f"Status: RUNNING [OK]\n"
                    f"Process ID: {info['pid']}\n"
                    f"CPU Usage: {info['cpu_percent']:.1f}%\n"
                    f"Memory Usage: {info['memory_percent']:.1f}%\n"
                    f"Command: {info['command']}\n\n"
                    f"Files:\n"
                    f"• PID file: {info['pid_file']}\n"
                    f"• Log file: {info['log_file']}\n\n"
                    f"The simulation is running in the background and will\n"
                    f"continue even if you close this application."
                )
            except ImportError:
                process_info = (
                    f"Equilibration Process Information\n"
                    f"{'=' * 40}\n\n"
                    f"Status: RUNNING [OK]\n"
                    f"Process ID: {info['pid']}\n\n"
                    f"Files:\n"
                    f"• PID file: {info['pid_file']}\n"
                    f"• Log file: {info['log_file']}\n\n"
                    f"The simulation is running in the background and will\n"
                    f"continue even if you close this application.\n\n"
                    f"Note: Install 'psutil' for detailed process information."
                )
        else:
            process_info = (
                f"Equilibration Process Information\n"
                f"{'=' * 40}\n\n"
                f"Status: NOT RUNNING [XX]\n"
                f"Process ID: {info['pid']} (process ended)\n\n"
                f"Files:\n"
                f"• PID file: {info['pid_file']}\n"
                f"• Log file: {info['log_file']}\n\n"
                f"The equilibration process has finished or was terminated.\n"
                f"Check the log file for completion status and any errors."
            )
        
        messagebox.showinfo("Background Process Info", process_info)
    
    def _run_equilibration(self):
        """Run the equilibration simulation."""
        # Update working directory from entry field
        workdir_text = self.workdir_entry.get().strip()
        if workdir_text:
            self.working_directory = Path(workdir_text)
            
        if not self._validate_inputs():
            return
        
        # Check if there's already a running equilibration
        if self._check_running_equilibration():
            result = messagebox.askyesno(
                "Equilibration Running",
                "There appears to be an equilibration process already running.\n\n"
                "Do you want to start a new one anyway?\n"
                "(This will not stop the existing process)"
            )
            if not result:
                return
        
        # Check if input files exist
        engine = self.engine_var.get().lower()
        output_dir = Path(self.working_directory) / self.equilibration_output_name / engine
        
        if not output_dir.exists():
            messagebox.showerror(
                "Error",
                "Input files not found. Please generate input files first."
            )
            return
        
        # Run in background thread
        thread = threading.Thread(target=self._run_equilibration_thread)
        thread.daemon = True
        thread.start()
        
        if self.status_callback:
            self.status_callback("Starting equilibration run...")
    
    def _run_equilibration_thread(self):
        """Run equilibration in background thread."""
        try:
            engine = self.engine_var.get()
            
            if engine == "NAMD":
                process = self._run_namd_equilibration()
                # Update UI immediately after starting background process
                self.after(0, lambda: self._on_equilibration_started(process))
            else:
                # For other engines, show not implemented message
                engine_name = engine  # Capture the engine name
                self.after(0, lambda: messagebox.showinfo(
                    "Info", 
                    f"{engine_name} execution not yet implemented"
                ))
                return
        
        except Exception as e:
            self.logger.error(f"Error starting equilibration: {e}")
            error_msg = str(e)  # Capture the error message
            self.after(0, lambda: messagebox.showerror(
                "Error", 
                f"Failed to start equilibration: {error_msg}"
            ))
    
    def _run_namd_equilibration(self):
        """Run NAMD equilibration as a detached background process."""
        import subprocess
        import os
        import sys
        
        output_dir = Path(self.working_directory) / self.equilibration_output_name / "namd"
        script_file = output_dir / "run_equilibration.sh"
        
        if not script_file.exists():
            raise FileNotFoundError(f"Run script not found: {script_file}. Please generate input files first.")
        
        # Make script executable
        os.chmod(script_file, 0o755)
        
        # Create log file for the background process
        log_file = output_dir / "equilibration_background.log"
        
        try:
            # Run as detached background process
            if sys.platform == "win32":
                # Windows: Use DETACHED_PROCESS to run independently
                process = subprocess.Popen(
                    ["bash", str(script_file)],
                    cwd=output_dir,
                    stdout=open(log_file, 'w'),
                    stderr=subprocess.STDOUT,
                    creationflags=subprocess.DETACHED_PROCESS | subprocess.CREATE_NEW_PROCESS_GROUP
                )
            else:
                # Unix/Linux: Use nohup and setsid for true background execution
                process = subprocess.Popen(
                    ["nohup", "bash", str(script_file)],
                    cwd=output_dir,
                    stdout=open(log_file, 'w'),
                    stderr=subprocess.STDOUT,
                    preexec_fn=os.setsid  # Create new session to detach from parent
                )
            
            # Save process ID for potential future monitoring
            pid_file = output_dir / "equilibration.pid"
            with open(pid_file, 'w') as f:
                f.write(str(process.pid))
            
            self.logger.info(f"NAMD equilibration started as background process (PID: {process.pid})")
            self.logger.info(f"Log file: {log_file}")
            self.logger.info(f"PID file: {pid_file}")
            
            # Don't wait for process to complete - let it run in background
            return process
            
        except Exception as e:
            self.logger.error(f"Failed to start background equilibration: {e}")
            raise
    
    def _on_equilibration_started(self, process):
        """Handle equilibration startup (not completion)."""
        output_dir = Path(self.working_directory) / self.equilibration_output_name / "namd"
        log_file = output_dir / "equilibration_background.log"
        pid_file = output_dir / "equilibration.pid"
        
        message = (
            f"Equilibration started successfully in the background!\n\n"
            f"Process ID: {process.pid}\n"
            f"Log file: {log_file}\n"
            f"PID file: {pid_file}\n\n"
            f"The simulation will continue running even if you close the application.\n"
            f"Use 'Refresh Now' to monitor progress."
        )
        
        messagebox.showinfo("Equilibration Started", message)
        
        if self.status_callback:
            self.status_callback(f"Equilibration running in background (PID: {process.pid})")
        
        # Start automatic progress monitoring
        if hasattr(self, 'auto_monitoring_var') and not self.auto_monitoring_var.get():
            self.auto_monitoring_var.set(True)
            self._toggle_auto_monitoring()
    
    def _on_equilibration_complete(self):
        """Handle equilibration completion (legacy method, kept for compatibility)."""
        messagebox.showinfo(
            "Success",
            "Equilibration process monitoring completed.\n"
            "Check the progress tracker for current status."
        )
        
        if self.status_callback:
            self.status_callback("Equilibration completed")
        
        # Final progress update
        self._refresh_equilibration_progress()
    
    def _save_protocol(self):
        """Save current protocol to file."""
        protocols = self._get_current_protocols()
        
        file_path = filedialog.asksaveasfilename(
            title="Save Protocol",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            initialdir=str(self.working_directory)
        )
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    json.dump(protocols, f, indent=2)
                
                messagebox.showinfo("Success", "Protocol saved successfully!")
                
                if self.status_callback:
                    self.status_callback(f"Protocol saved: {Path(file_path).name}")
            
            except Exception as e:
                self.logger.error(f"Error saving protocol: {e}")
                messagebox.showerror("Error", f"Failed to save protocol: {str(e)}")
    
    def _load_protocol(self):
        """Load protocol from file."""
        file_path = filedialog.askopenfilename(
            title="Load Protocol",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            initialdir=str(self.working_directory)
        )
        
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    protocols = json.load(f)
                
                # Update widgets with loaded values
                self._update_protocol_widgets(protocols)
                
                messagebox.showinfo("Success", "Protocol loaded successfully!")
                
                if self.status_callback:
                    self.status_callback(f"Protocol loaded: {Path(file_path).name}")
            
            except Exception as e:
                self.logger.error(f"Error loading protocol: {e}")
                messagebox.showerror("Error", f"Failed to load protocol: {str(e)}")
    
    def _get_current_protocols(self) -> Dict[str, Any]:
        """Get current protocol parameters from widgets."""
        # Always get values from widgets if they exist
        # This ensures GUI values take precedence when user makes changes
        protocols = {}
        
        # Get the current CHARMM-GUI protocol template to use as base
        protocol_template = self._get_current_protocol_template()
        
        for stage_key, stage_widgets in self.stage_widgets.items():
            # Use the CHARMM-GUI protocol template as base
            if stage_key in protocol_template:
                stage_data = protocol_template[stage_key].copy()
            else:
                # This should not happen with CHARMM-GUI only, but provide fallback
                stage_data = {}
            
            # Update with widget values (always use GUI values)
            # Get time in nanoseconds and calculate steps
            time_ns = float(stage_widgets['time_ns'].get())
            timestep = float(stage_widgets['timestep'].get())
            calculated_steps = int(time_ns * 1e6 / timestep)
            
            stage_data['time_ns'] = time_ns
            stage_data['steps'] = calculated_steps
            stage_data['ensemble'] = stage_widgets['ensemble'].get()
            stage_data['temperature'] = float(stage_widgets['temperature'].get())
            stage_data['timestep'] = timestep
            
            # Update minimize steps only for Equilibration 1
            if stage_key == "Equilibration 1" and 'minimize_steps' in stage_widgets:
                stage_data['minimize_steps'] = int(stage_widgets['minimize_steps'].get())
            
            # Update DCD frequency
            if 'dcd_freq' in stage_widgets:
                stage_data['dcd_freq'] = int(stage_widgets['dcd_freq'].get())
            
            # Update margin parameter
            if 'margin' in stage_widgets:
                stage_data['margin'] = float(stage_widgets['margin'].get())
            
            if 'pressure' in stage_widgets:
                stage_data['pressure'] = float(stage_widgets['pressure'].get())
            
            if 'surface_tension' in stage_widgets:
                stage_data['surface_tension'] = float(stage_widgets['surface_tension'].get())
            
            # Update constraints
            for constraint_name, constraint_var in stage_widgets['constraints'].items():
                stage_data['constraints'][constraint_name] = float(constraint_var.get())
            
            # Update computational resources
            if 'use_gpu' in stage_widgets:
                stage_data['use_gpu'] = stage_widgets['use_gpu'].get()
            
            if 'cpu_cores' in stage_widgets:
                stage_data['cpu_cores'] = int(stage_widgets['cpu_cores'].get())
            
            if 'gpu_id' in stage_widgets:
                stage_data['gpu_id'] = int(stage_widgets['gpu_id'].get())
            
            if 'num_gpus' in stage_widgets:
                stage_data['num_gpus'] = int(stage_widgets['num_gpus'].get())
            
            # Update bilayer thickness restraint settings if present
            if 'force_constant' in stage_widgets:
                stage_data['force_constant'] = float(stage_widgets['force_constant'].get())
            
            protocols[stage_key] = stage_data
        
        return protocols
    
    def _update_protocol_widgets(self, protocols: Dict[str, Any]):
        """Update protocol widgets with loaded values."""
        for stage_key, stage_data in protocols.items():
            if stage_key in self.stage_widgets:
                stage_widgets = self.stage_widgets[stage_key]
                
                # Update basic parameters
                # Use time_ns if available, otherwise calculate from steps and timestep
                if 'time_ns' in stage_data:
                    time_ns_value = stage_data['time_ns']
                else:
                    # Calculate time_ns from steps and timestep for backward compatibility
                    steps = stage_data.get('steps', 125000)
                    timestep = stage_data.get('timestep', 0.001)
                    time_ns_value = steps * timestep / 1e6
                
                stage_widgets['time_ns'].set(str(time_ns_value))
                stage_widgets['ensemble'].set(stage_data['ensemble'])
                stage_widgets['temperature'].set(str(stage_data['temperature']))
                stage_widgets['timestep'].set(str(stage_data['timestep']))
                
                # Update minimize steps only for Equilibration 1
                if stage_key == "Equilibration 1" and 'minimize_steps' in stage_widgets and 'minimize_steps' in stage_data:
                    stage_widgets['minimize_steps'].set(str(stage_data['minimize_steps']))
                
                # Update DCD frequency
                if 'dcd_freq' in stage_widgets and 'dcd_freq' in stage_data:
                    stage_widgets['dcd_freq'].set(str(stage_data['dcd_freq']))
                
                # Update margin parameter
                if 'margin' in stage_widgets and 'margin' in stage_data:
                    stage_widgets['margin'].set(str(stage_data['margin']))
                
                if 'pressure' in stage_widgets and 'pressure' in stage_data:
                    stage_widgets['pressure'].set(str(stage_data['pressure']))
                
                if 'surface_tension' in stage_widgets and 'surface_tension' in stage_data:
                    stage_widgets['surface_tension'].set(str(stage_data['surface_tension']))
                
                # Update computational resources
                if 'use_gpu' in stage_widgets and 'use_gpu' in stage_data:
                    stage_widgets['use_gpu'].set(stage_data['use_gpu'])
                
                if 'cpu_cores' in stage_widgets and 'cpu_cores' in stage_data:
                    stage_widgets['cpu_cores'].set(str(stage_data['cpu_cores']))
                
                if 'gpu_id' in stage_widgets and 'gpu_id' in stage_data:
                    stage_widgets['gpu_id'].set(str(stage_data['gpu_id']))
                
                if 'num_gpus' in stage_widgets and 'num_gpus' in stage_data:
                    stage_widgets['num_gpus'].set(str(stage_data['num_gpus']))
                
                # Update constraints
                for constraint_name, constraint_value in stage_data['constraints'].items():
                    if constraint_name in stage_widgets['constraints']:
                        stage_widgets['constraints'][constraint_name].set(str(constraint_value))
                
                # Update ABF settings if present
                # Skip ABF-related settings (removed)
    
    def _validate_inputs(self) -> bool:
        """Validate all input parameters."""
        # Update working directory from entry field if not already a Path
        if isinstance(self.working_directory, str):
            self.working_directory = Path(self.working_directory)
        
        # Check working directory
        if not self.working_directory or not self.working_directory.exists():
            messagebox.showerror("Error", "Please select a valid working directory")
            return False
        
        # Check if AMBER system files exist (basic validation)
        # Look for topology and coordinate files
        amber_files_found = False
        
        # Check for topology files
        for ext in ['.prmtop', '.top']:
            if list(self.working_directory.glob(f"**/*{ext}")):
                amber_files_found = True
                break
        
        if not amber_files_found:
            result = messagebox.askyesno(
                "Warning",
                "No AMBER topology files (.prmtop or .top) found in working directory.\n"
                "NAMD requires AMBER topology and coordinate files.\n"
                "Do you want to continue anyway?"
            )
            if not result:
                return False
        
        # Check for coordinate files
        coord_files_found = False
        for ext in ['.crd', '.rst', '.inpcrd']:
            if list(self.working_directory.glob(f"**/*{ext}")):
                coord_files_found = True
                break
        
        if not coord_files_found:
            result = messagebox.askyesno(
                "Warning", 
                "No AMBER coordinate files (.crd, .rst, or .inpcrd) found in working directory.\n"
                "NAMD requires AMBER coordinate files.\n"
                "Do you want to continue anyway?"
            )
            if not result:
                return False
        
        # Validate numeric parameters
        try:
            protocols = self._get_current_protocols()
            for stage_key, stage_data in protocols.items():
                # Check required numeric fields
                int(stage_data['steps'])
                float(stage_data['temperature'])
                float(stage_data['timestep'])
                
                if stage_data['ensemble'] == 'NPT':
                    float(stage_data['pressure'])
                
                # Check constraints
                for constraint_value in stage_data['constraints'].values():
                    float(constraint_value)
        
        except (ValueError, TypeError) as e:
            messagebox.showerror("Error", f"Invalid numeric parameter: {str(e)}")
            return False
        
        return True
    
    def _check_amber_files(self, output_dir: Path):
        """
        Check for required AMBER topology and coordinate files.
        
        Args:
            output_dir: Directory where NAMD files will be generated
        """
        self.logger.info("Checking for AMBER topology and coordinate files")
        
        # Ensure working_directory is a Path object
        working_dir = Path(self.working_directory)
        
        # Required AMBER files for NAMD
        required_files = [
            ("Topology file", [".prmtop", ".top"]),
            ("Coordinate file", [".crd", ".rst", ".inpcrd"])
        ]
        
        missing_files = []
        found_files = []
        
        for file_type, extensions in required_files:
            found = False
            for ext in extensions:
                files = list(working_dir.glob(f"**/*{ext}"))
                if files:
                    found_files.append(f"{file_type}: {files[0].name}")
                    found = True
                    break
            
            if not found:
                missing_files.append(f"{file_type} ({', '.join(extensions)})")
        
        # Create a status file
        status_content = "# AMBER Files Status for NAMD\n\n"
        
        if found_files:
            status_content += "Found files:\n"
            for file_info in found_files:
                status_content += f"[OK] {file_info}\n"
            status_content += "\n"
        
        if missing_files:
            status_content += "Missing files:\n"
            for file_info in missing_files:
                status_content += f"[XX] {file_info}\n"
            status_content += "\n"
            status_content += "Note: NAMD requires AMBER topology (.prmtop/.top) and coordinate (.inpcrd/.pdb) files.\n"
            status_content += "These files should be generated during the system preparation stage.\n"
        else:
            status_content += "All required AMBER files found!\n"
        
        status_content += "\nFor NAMD with AMBER force field, ensure:\n"
        status_content += "1. Topology file contains all force field parameters\n"
        status_content += "2. Coordinate file contains box dimensions (if periodic)\n"
        status_content += "3. Files are compatible with NAMD's AMBER support\n"
        
        status_file = output_dir / "AMBER_FILES_STATUS.txt"
        with open(status_file, 'w') as f:
            f.write(status_content)
        
        logger.info(f"Created AMBER files status: {status_file}")
        
        if missing_files:
            logger.warning(f"Missing AMBER files: {missing_files}")
            logger.warning("NAMD simulation may fail without proper AMBER files")
        else:
            logger.info("All required AMBER files are available")
    
    def on_stage_shown(self):
        """Called when this stage is shown."""
        # Update current PDB if available
        if self.get_current_pdb:
            self.current_pdb_file = self.get_current_pdb()
        
        # Initialize NAMD settings by default
        if not hasattr(self, 'namd_path_var'):
            self._on_engine_changed("NAMD")
    
    def on_pdb_changed(self, pdb_file: Optional[str]):
        """Handle PDB file change."""
        self.current_pdb_file = pdb_file
        
        if self.status_callback:
            if pdb_file:
                self.status_callback(f"PDB file: {Path(pdb_file).name}")
            else:
                self.status_callback("No PDB file selected")
    
    def cleanup(self):
        """Clean up resources when frame is destroyed."""
        # Stop auto-refresh timer
        self._stop_auto_refresh()
        
        # Any other cleanup can be added here
        pass
    

    def update_fonts(self, scaled_fonts):
        """Update all fonts in the equilibration frame."""
        try:
            def safe_cfg(widget, font=None, width=None, height=None):
                try:
                    if widget and hasattr(widget, 'configure'):
                        kwargs = {}
                        if font is not None:
                            kwargs['font'] = font
                        if width is not None:
                            kwargs['width'] = width
                        if height is not None:
                            kwargs['height'] = height
                        if kwargs:
                            widget.configure(**kwargs)
                except Exception:
                    pass

            # Headings
            for attr in ['title_label', 'workdir_label', 'engine_label', 'protocol_label', 'progress_label', 'action_label', 'summary_label', 'inputfolder_label', 'outputname_label']:
                if hasattr(self, attr):
                    safe_cfg(getattr(self, attr), scaled_fonts['heading'])

            # Body text labels and description
            if hasattr(self, 'description_label'):
                safe_cfg(self.description_label, scaled_fonts['body'])
            if hasattr(self, 'inputfolder_sublabel'):
                safe_cfg(self.inputfolder_sublabel, scaled_fonts['body'])

            # Inputs and menus
            if hasattr(self, 'workdir_entry'):
                safe_cfg(self.workdir_entry, scaled_fonts['body'])
            if hasattr(self, 'engine_menu'):
                safe_cfg(self.engine_menu, scaled_fonts['body'])
            if hasattr(self, 'inputfolder_entry'):
                safe_cfg(self.inputfolder_entry, scaled_fonts['body'])
            if hasattr(self, 'outputname_entry'):
                safe_cfg(self.outputname_entry, scaled_fonts['body'])
            
            # CHARMM-GUI scheme controls
            if hasattr(self, 'scheme_type_label'):
                safe_cfg(self.scheme_type_label, scaled_fonts['body'])
            if hasattr(self, 'scheme_type_menu'):
                safe_cfg(self.scheme_type_menu, scaled_fonts['body'])
            # Note: thickness_label and thickness_entry removed as bilayer thickness input is no longer needed

            # Monitoring controls
            for attr in ['auto_monitoring_checkbox', 'interval_label', 'interval_entry', 'interval_unit_label']:
                if hasattr(self, attr):
                    # Use body font for better readability
                    safe_cfg(getattr(self, attr), scaled_fonts['body'])

            # Action buttons (consistent sizing)
            for attr in [
                'generate_files_btn', 'run_equilibration_btn', 'save_protocol_btn', 'load_protocol_btn',
                'workdir_browse_btn', 'inputfolder_browse_btn', 'refresh_progress_btn', 'start_monitor_btn', 
                'stop_monitor_btn', 'background_info_btn', 'open_log_btn', 'open_output_btn', 'open_dir_btn', 
                'reset_btn', 'close_btn']:
                if hasattr(self, attr):
                    safe_cfg(getattr(self, attr), scaled_fonts['body'], width=WIDGET_SIZES.get('button_width', 120), height=WIDGET_SIZES.get('button_height', 32))

            # Engine settings container (labels/entries created dynamically)
            def update_container_fonts(container):
                try:
                    for w in container.winfo_children():
                        cls = w.__class__.__name__
                        if cls in ('CTkLabel',):
                            safe_cfg(w, scaled_fonts['body'])
                        elif cls in ('CTkButton',):
                            safe_cfg(w, scaled_fonts['body'], width=120, height=32)
                        elif cls in ('CTkCheckBox',):
                            safe_cfg(w, scaled_fonts['body'])
                        elif cls in ('CTkEntry', 'CTkOptionMenu'):
                            safe_cfg(w, scaled_fonts['body'])
                        elif cls in ('CTkTextbox',):
                            safe_cfg(w, scaled_fonts.get('small', scaled_fonts['body']))
                        # Recurse into frames/tabviews
                        if hasattr(w, 'winfo_children') and w.winfo_children():
                            update_container_fonts(w)
                except Exception:
                    pass

            if hasattr(self, 'engine_settings_frame'):
                update_container_fonts(self.engine_settings_frame)

            # Stage widgets: labels, entries, menus, checkboxes
            if hasattr(self, 'stage_widgets'):
                for widgets in self.stage_widgets.values():
                    for key in [
                        'desc_label', 'ensemble_label', 'temperature_label', 'temperature_unit_label',
                        'pressure_label', 'pressure_unit_label', 'surface_tension_label', 'surface_tension_unit_label', 
                        'timestep_label', 'timestep_unit_label', 'time_label', 'time_unit_label', 'steps_label', 'steps_unit_label',
                        'minimize_label', 'minimize_info_label', 'dcd_freq_label', 'dcd_freq_unit_label',
                        'margin_label', 'margin_unit_label', 'resources_label', 'gpu_label', 'cpu_label', 'cpu_unit_label',
                        'gpu_id_label', 'gpu_id_hint_label', 'num_gpus_label', 'num_gpus_unit_label',
                        'constraints_label', 'abf_label']:
                        if key in widgets:
                            safe_cfg(widgets[key], scaled_fonts['body'])
                    for key in [
                        'time_entry', 'steps_entry', 'minimize_entry', 'temperature_entry', 'pressure_entry', 'surface_tension_entry',
                        'timestep_entry', 'dcd_freq_entry', 'margin_entry', 'cpu_cores_entry', 'gpu_id_entry', 'num_gpus_entry']:
                        if key in widgets:
                            safe_cfg(widgets[key], scaled_fonts['body'])
                    # Constraint entries and labels
                    if 'constraint_entries' in widgets:
                        for e in widgets['constraint_entries']:
                            safe_cfg(e, scaled_fonts['body'])
                    if 'constraint_labels' in widgets:
                        for l in widgets['constraint_labels']:
                            safe_cfg(l, scaled_fonts['body'])
                    # Note: ensemble_menu removed - ensemble is now read-only per scheme
                    if 'abf_checkbox' in widgets:
                        safe_cfg(widgets['abf_checkbox'], scaled_fonts['body'])

            # Tabview segmented button (tab labels)
            if hasattr(self, 'stages_tabview') and hasattr(self.stages_tabview, 'segmented_button'):  # type: ignore
                try:
                    self.stages_tabview.segmented_button.configure(font=scaled_fonts['body'])  # type: ignore
                except Exception:
                    pass

            # Progress widgets (per stage display)
            if hasattr(self, 'progress_widgets'):
                for widgets in self.progress_widgets.values():
                    for key, widget in widgets.items():
                        if hasattr(widget, 'configure') and 'label' in key:
                            # Stage title/body labels vs status small
                            is_status = key in ('status_label', 'steps_label')
                            safe_cfg(widget, scaled_fonts['small'] if is_status else scaled_fonts['body'])

            # Summary text area
            if hasattr(self, 'summary_text'):
                safe_cfg(self.summary_text, scaled_fonts.get('code', scaled_fonts['small']))

            # Protocol selector checkbox
        except Exception as e:
            logger.warning(f"Error updating fonts in EquilibrationFrame: {e}")
