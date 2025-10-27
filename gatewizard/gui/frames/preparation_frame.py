# gatewizard/gui/frames/preparation_frame.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Preparation frame for protein structure preparation and pKa prediction.

This module provides the GUI for running Propka analysis, viewing results,
and applying protonation states to protein structures.
"""

import tkinter as tk
from tkinter import messagebox, filedialog
from typing import Optional, Callable, List, Dict, Any
from pathlib import Path
import threading
from datetime import datetime

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, WIDGET_SIZES, LAYOUT, ERROR_MESSAGES, SUCCESS_MESSAGES
)
from gatewizard.core.preparation import PreparationManager, PreparationError, PROTONATION_STATES
from gatewizard.core.builder import Builder
from gatewizard.tools.validators import SystemValidator
from gatewizard.utils.logger import get_logger
from gatewizard.utils.protein_capping import ProteinCapper, ProteinCappingError

logger = get_logger(__name__)

class PreparationFrame(ctk.CTkFrame):
    """
    Frame for protein structure preparation and protonation state management.
    
    This frame handles pKa analysis, result visualization, and protonation
    state application for protein structures.
    """
    
    def __init__(
        self,
        parent,
        get_current_pdb: Optional[Callable[[], Optional[str]]] = None,
        status_callback: Optional[Callable[[str], None]] = None,
        initial_directory: Optional[str] = None
    ):
        """
        Initialize the Propka frame.
        
        Args:
            parent: Parent widget
            get_current_pdb: Callback to get current PDB file
            status_callback: Callback for status updates
            initial_directory: Initial directory for file dialogs
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['content_bg'])
        
        self.get_current_pdb = get_current_pdb
        self.status_callback = status_callback
        self.initial_directory = initial_directory or str(Path.cwd())
        
        # State variables
        self.current_pdb_file = None
        self.analysis_file = None  # Track which file was actually used for analysis (original or capped)
        self.residue_mapping = {}  # Track mapping from original to capped residue IDs
        self.capping_used = False  # Track if capping was used in the current analysis
        self.analysis_timestamp = None  # Track when analysis was performed
        self.preparation_manager = PreparationManager()
        self.builder = Builder()
        self.validator = SystemValidator()
        self.protein_capper = ProteinCapper()  # Initialize directly since MDAnalysis is required
        self.protonable_residues = []
        self.custom_states = {}
        self.detected_disulfide_bonds = []
        
        # ComboBox tracking for visual feedback
        self.custom_state_combos = {}  # Dictionary to store ComboBox widgets by residue_id
        
        # Warning system for when options change after analysis
        self.analysis_run = False  # Track if analysis has been run
        self.options_changed_after_analysis = False  # Track if options changed after analysis
        self.last_analysis_settings = {}  # Store settings used in last analysis
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        
        # Enable mouse wheel scrolling
        self._bind_mouse_wheel()
        
        # Load initial file if callback is available
        self._load_initial_file()
    
    def _load_initial_file(self):
        """Load initial file from callback if available."""
        if self.get_current_pdb:
            pdb_file = self.get_current_pdb()
            if pdb_file:
                self.load_pdb_file(pdb_file)
    
    def load_pdb_file(self, file_path: str):
        """Load a PDB file into the working file field."""
        if file_path and Path(file_path).exists():
            self.working_file_entry.delete(0, "end")
            self.working_file_entry.insert(0, file_path)
            self.current_pdb_file = file_path
            self.analysis_file = None  # Reset analysis file when loading new PDB
            self.residue_mapping = {}  # Reset residue mapping when loading new PDB
            self.capping_used = False
            
            # Auto-generate export filename
            export_path = self._generate_export_filename(file_path)
            self.export_file_entry.delete(0, "end")
            self.export_file_entry.insert(0, export_path)
            
            # Enable analysis controls
            self._set_analysis_enabled(True)
            
            if self.status_callback:
                self.status_callback(f"Loaded file: {Path(file_path).name}")
    
    def _create_widgets(self):
        """Create all widgets for the Propka frame."""
        # Main horizontal container
        self.main_container = ctk.CTkFrame(self, fg_color="transparent")
        
        # LEFT PANEL - Controls and Configuration (fixed width)
        self.left_panel = ctk.CTkFrame(self.main_container, fg_color=COLOR_SCHEME['content_inside_bg'], width=500)
        self.left_panel.pack_propagate(False)  # Maintain fixed width
        
        # Create a scrollable frame inside the left panel
        self.left_scroll = ctk.CTkScrollableFrame(self.left_panel, fg_color="transparent")
        self.left_scroll.pack(fill="both", expand=True, padx=5, pady=5)
        
        # File selection section
        self.file_section = ctk.CTkFrame(self.left_scroll, fg_color=COLOR_SCHEME['canvas'])
        
        self.file_label = ctk.CTkLabel(
            self.file_section,
            text="Working File Selection",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Current working file display
        self.working_file_frame = ctk.CTkFrame(self.file_section, fg_color="transparent")
        
        self.working_file_label = ctk.CTkLabel(
            self.working_file_frame,
            text="Current Working File:",
            font=FONTS['body']
        )
        
        self.working_file_entry = ctk.CTkEntry(
            self.working_file_frame,
            placeholder_text="No file selected...",
            height=WIDGET_SIZES['entry_height']
        )
        
        self.browse_file_button = ctk.CTkButton(
            self.working_file_frame,
            text="Browse",
            width=100,
            height=WIDGET_SIZES['button_height'],
            command=self._browse_working_file
        )
        
        # Export file section
        self.export_file_frame = ctk.CTkFrame(self.file_section, fg_color="transparent")
        
        self.export_file_label = ctk.CTkLabel(
            self.export_file_frame,
            text="Export Protonated File:",
            font=FONTS['body']
        )
        
        self.export_file_entry = ctk.CTkEntry(
            self.export_file_frame,
            placeholder_text="Will be auto-generated...",
            height=WIDGET_SIZES['entry_height']
        )
        
        self.export_browse_button = ctk.CTkButton(
            self.export_file_frame,
            text="Browse",
            width=100,
            height=WIDGET_SIZES['button_height'],
            command=self._browse_export_file
        )
        
        # RIGHT PANEL - Analysis Results (expandable)
        self.right_panel = ctk.CTkFrame(self.main_container, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.analysis_results_label = ctk.CTkLabel(
            self.right_panel,  # Back to right_panel as parent
            text="Analysis Results",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Container for analysis results
        self.analysis_results_container = ctk.CTkFrame(self.right_panel, fg_color=COLOR_SCHEME['viewer_bg'])
        
        # Analysis section (moved to left panel)
        self.analysis_section = ctk.CTkFrame(self.left_scroll, fg_color=COLOR_SCHEME['canvas'])
        
        self.analysis_label = ctk.CTkLabel(
            self.analysis_section,
            text="Propka Analysis",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Analysis parameters
        self.params_frame = ctk.CTkFrame(self.analysis_section, fg_color="transparent")
        
        self.ph_label = ctk.CTkLabel(
            self.params_frame,
            text="Target pH:",
            font=FONTS['body']
        )
        
        self.ph_entry = ctk.CTkEntry(
            self.params_frame,
            width=80,
            height=WIDGET_SIZES['entry_height']
        )
        self.ph_entry.insert(0, "7.0")
        self.ph_entry.bind('<KeyRelease>', self._on_option_change)
        self.ph_entry.bind('<FocusOut>', self._on_option_change)
        
        # Protein capping option
        self.cap_protein_var = ctk.BooleanVar(value=False)
        self.cap_protein_check = ctk.CTkCheckBox(
            self.params_frame,
            text="Cap protein termini (ACE/NME)",
            variable=self.cap_protein_var,
            font=FONTS['body'],
            command=self._on_cap_toggle
        )
        
        self.run_button = ctk.CTkButton(
            self.params_frame,
            text="Run Analysis",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._run_analysis
        )
        
        # Disulfide bonds section (moved to left panel)
        self.disulfide_section = ctk.CTkFrame(self.left_scroll, fg_color=COLOR_SCHEME['canvas'])
        
        self.disulfide_label = ctk.CTkLabel(
            self.disulfide_section,
            text="Disulfide Bond Detection",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Disulfide detection controls
        self.disulfide_controls_frame = ctk.CTkFrame(self.disulfide_section, fg_color="transparent")
        
        self.auto_detect_var = ctk.BooleanVar(value=True)
        self.auto_detect_check = ctk.CTkCheckBox(
            self.disulfide_controls_frame,
            text="Auto-detect disulfide bonds",
            variable=self.auto_detect_var,
            font=FONTS['body'],
            command=self._on_option_change
        )
        
        self.distance_label = ctk.CTkLabel(
            self.disulfide_controls_frame,
            text="Max S-S distance (Å):",
            font=FONTS['body']
        )
        
        self.distance_entry = ctk.CTkEntry(
            self.disulfide_controls_frame,
            width=80,
            height=WIDGET_SIZES['entry_height']
        )
        self.distance_entry.insert(0, "2.5")
        self.distance_entry.bind('<KeyRelease>', self._on_option_change)
        self.distance_entry.bind('<FocusOut>', self._on_option_change)
        
        self.detect_button = ctk.CTkButton(
            self.disulfide_controls_frame,
            text="Detect Bonds",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._detect_disulfide_bonds,
            state="disabled"
        )
        
        # Disulfide bonds display
        self.disulfide_display_frame = ctk.CTkFrame(self.disulfide_section, fg_color="transparent")
        
        self.disulfide_list_label = ctk.CTkLabel(
            self.disulfide_display_frame,
            text="Detected S-S Bonds:",
            font=FONTS['body']
        )
        
        self.disulfide_text = ctk.CTkTextbox(
            self.disulfide_display_frame,
            height=80,
            font=FONTS['small']
        )
        self.disulfide_text.insert("1.0", "No disulfide bonds detected")
        self.disulfide_text.configure(state="disabled")
        
        # Action buttons section (moved to left panel)
        self.actions_section = ctk.CTkFrame(self.left_scroll, fg_color=COLOR_SCHEME['canvas'])
        
        self.actions_label = ctk.CTkLabel(
            self.actions_section,
            text="Actions",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.actions_frame = ctk.CTkFrame(self.actions_section, fg_color="transparent")
        
        self.apply_button = ctk.CTkButton(
            self.actions_frame,
            text="Apply States & Run pdb4amber",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._apply_protonation_states,
            state="disabled"
        )
        
        self.export_button = ctk.CTkButton(
            self.actions_frame,
            text="Export Results",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._export_results,
            state="disabled"
        )
        
        self.reset_button = ctk.CTkButton(
            self.actions_frame,
            text="Reset Custom",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._reset_custom_states,
            state="disabled"
        )
        
        # Summary section (moved to left panel)
        self.summary_section = ctk.CTkFrame(self.left_scroll, fg_color=COLOR_SCHEME['canvas'])
        
        self.summary_label = ctk.CTkLabel(
            self.summary_section,
            text="Analysis Summary",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.summary_text = ctk.CTkTextbox(
            self.summary_section,
            height=150,  # Increased from 100 to 150
            font=FONTS['small']
        )
        self.summary_text.insert("1.0", "No analysis performed yet")
        self.summary_text.configure(state="disabled")
        
        # Create detailed results view in right panel
        self._create_results_section()
        
        # Initially disable analysis controls
        self._set_analysis_enabled(False)
        
        # Initialize all buttons to blue color
        self._set_all_buttons_color("blue")
    
    def _create_results_section(self):
        """Create the detailed results section in the right panel."""
        # Detailed results display
        self.results_detailed_section = ctk.CTkFrame(self.analysis_results_container, fg_color=COLOR_SCHEME['canvas'])
        
        self.results_detailed_label = ctk.CTkLabel(
            self.results_detailed_section,
            text="Residue Details",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Results display
        self.results_frame = ctk.CTkFrame(self.results_detailed_section, fg_color="transparent")
        
        # Column headers with sorting functionality
        self.headers_frame = ctk.CTkFrame(self.results_frame, fg_color=COLOR_SCHEME['canvas'])
        
        # Sort state tracking
        self.sort_column = None
        self.sort_ascending = True
        
        self.header_residue = ctk.CTkButton(
            self.headers_frame,
            text="Residue",
            font=FONTS['body'],
            width=80,
            height=30,
            command=lambda: self._sort_results('residue'),
            fg_color=COLOR_SCHEME['canvas'],
            text_color=COLOR_SCHEME['text'],
            hover_color=COLOR_SCHEME['hover']
        )
        
        self.header_id = ctk.CTkButton(
            self.headers_frame,
            text="ID (Orig-New)" if self.cap_protein_var.get() else "ID",
            font=FONTS['body'],
            width=120 if self.cap_protein_var.get() else 60,  # Wider when capping is enabled
            height=30,
            command=lambda: self._sort_results('res_id'),
            fg_color=COLOR_SCHEME['canvas'],
            text_color=COLOR_SCHEME['text'],
            hover_color=COLOR_SCHEME['hover']
        )
        
        self.header_chain = ctk.CTkButton(
            self.headers_frame,
            text="Chain",
            font=FONTS['body'],
            width=60,
            height=30,
            command=lambda: self._sort_results('chain'),
            fg_color=COLOR_SCHEME['canvas'],
            text_color=COLOR_SCHEME['text'],
            hover_color=COLOR_SCHEME['hover']
        )
        
        self.header_pka = ctk.CTkButton(
            self.headers_frame,
            text="pKa",
            font=FONTS['body'],
            width=80,
            height=30,
            command=lambda: self._sort_results('pka'),
            fg_color=COLOR_SCHEME['canvas'],
            text_color=COLOR_SCHEME['text'],
            hover_color=COLOR_SCHEME['hover']
        )
        
        self.header_state_ph = ctk.CTkButton(
            self.headers_frame,
            text="State at pH",
            font=FONTS['body'],
            width=100,
            height=30,
            command=lambda: self._sort_results('state_at_ph'),
            fg_color=COLOR_SCHEME['canvas'],
            text_color=COLOR_SCHEME['text'],
            hover_color=COLOR_SCHEME['hover']
        )
        
        self.header_custom = ctk.CTkLabel(
            self.headers_frame,
            text="Custom State",
            font=FONTS['body'],
            width=120
        )
        
        # Scrollable results area
        self.results_scroll = ctk.CTkScrollableFrame(
            self.results_frame,
            height=400,
            fg_color=COLOR_SCHEME['background']
        )
        
        # Add direct mouse wheel binding to results_scroll
        self._setup_results_scroll_binding()
    
    def _setup_results_scroll_binding(self):
        """Set up direct mouse wheel binding for the results scroll frame."""
        def scroll_results(event):
            # Direct scrolling for the results area
            canvas = None
            for attr_name in ['_parent_canvas', '_scrollable_frame_canvas', 'canvas', '_canvas']:
                if hasattr(self.results_scroll, attr_name):
                    canvas = getattr(self.results_scroll, attr_name)
                    if canvas and hasattr(canvas, 'yview_scroll'):
                        break
            
            if canvas:
                try:
                    if hasattr(event, 'delta'):
                        # Windows/macOS
                        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
                    else:
                        # Linux
                        if event.num == 4:
                            canvas.yview_scroll(-1, "units")
                        elif event.num == 5:
                            canvas.yview_scroll(1, "units")
                except:
                    pass
            return "break"
        
        # Bind directly to the scrollable frame
        try:
            self.results_scroll.bind("<MouseWheel>", scroll_results, add="+")
            self.results_scroll.bind("<Button-4>", scroll_results, add="+")
            self.results_scroll.bind("<Button-5>", scroll_results, add="+")
            
            # Also bind to the internal canvas if accessible
            def bind_to_canvas():
                try:
                    if not self.winfo_exists():
                        return
                    for attr_name in ['_parent_canvas', '_scrollable_frame_canvas', 'canvas', '_canvas']:
                        if hasattr(self.results_scroll, attr_name):
                            canvas = getattr(self.results_scroll, attr_name)
                            if canvas:
                                canvas.bind("<MouseWheel>", scroll_results, add="+")
                                canvas.bind("<Button-4>", scroll_results, add="+")
                                canvas.bind("<Button-5>", scroll_results, add="+")
                                break
                except:
                    pass
            
            # Schedule canvas binding after a short delay
            self.after(50, bind_to_canvas)
        except:
            pass
    
    def _setup_layout(self):
        """Setup the layout with left panel for controls and right panel for analysis results."""
        # Main container fills the entire frame
        self.main_container.pack(fill="both", expand=True, padx=5, pady=5)
        
        # LEFT PANEL - Fixed width for controls
        self.left_panel.pack(side="left", fill="y", padx=(0, 5), pady=0)
        
        # File selection section - very compact
        self.file_section.pack(fill="x", padx=10, pady=(5, 2))
        
        self.file_label.pack(anchor="w", padx=10, pady=(2, 1))
        
        # Working file frame - single line layout
        self.working_file_frame.pack(fill="x", padx=10, pady=(0, 1))
        
        # Pack working file widgets in horizontal layout
        self.working_file_label.pack(side="left", pady=0)
        self.browse_file_button.pack(side="right", padx=(5, 0))
        self.working_file_entry.pack(fill="x", padx=(5, 5))
        
        # Export file frame - single line layout  
        self.export_file_frame.pack(fill="x", padx=10, pady=(0, 2))
        
        # Pack export file widgets in horizontal layout
        self.export_file_label.pack(side="left", pady=0)
        self.export_browse_button.pack(side="right", padx=(5, 0))
        self.export_file_entry.pack(fill="x", padx=(5, 5))
        
        # Analysis section
        self.analysis_section.pack(fill="x", padx=10, pady=5)
        
        self.analysis_label.pack(anchor="w", padx=10, pady=(10, 5))
        
        self.params_frame.pack(fill="x", padx=10, pady=(0, 10))
        
        # Parameters in simple vertical layout
        self.ph_label.pack(anchor="w", pady=(0, 2))
        self.ph_entry.pack(anchor="w", pady=(0, 5))
        
        # Protein capping option
        self.cap_protein_check.pack(anchor="w", pady=(0, 5))
        
        self.run_button.pack(fill="x", pady=(5, 0))
        
        # Disulfide bonds section
        self.disulfide_section.pack(fill="x", padx=10, pady=5)
        
        self.disulfide_label.pack(anchor="w", padx=10, pady=(10, 5))
        
        # Disulfide controls
        self.disulfide_controls_frame.pack(fill="x", padx=10, pady=(0, 5))
        
        self.auto_detect_check.pack(anchor="w", pady=2)
        
        # Distance controls in simple vertical layout  
        self.distance_label.pack(anchor="w", pady=(5, 2))
        self.distance_entry.pack(anchor="w", pady=(0, 5))
        
        self.detect_button.pack(fill="x", pady=(5, 0))
        
        # Disulfide display
        self.disulfide_display_frame.pack(fill="x", padx=10, pady=(0, 10))
        
        self.disulfide_list_label.pack(anchor="w", pady=(0, 2))
        self.disulfide_text.pack(fill="x")
        
        # Actions section
        self.actions_section.pack(fill="x", padx=10, pady=5)
        
        self.actions_label.pack(anchor="w", padx=10, pady=(10, 5))
        
        self.actions_frame.pack(fill="x", padx=10, pady=(0, 10))
        
        self.apply_button.pack(fill="x", pady=2)
        self.export_button.pack(fill="x", pady=2)
        self.reset_button.pack(fill="x", pady=2)
        
        # Summary section - increased height
        self.summary_section.pack(fill="both", expand=True, padx=10, pady=(5, 10))
        
        self.summary_label.pack(anchor="w", padx=10, pady=(10, 5))
        self.summary_text.pack(fill="both", expand=True, padx=10, pady=(0, 10))
        
        # RIGHT PANEL - Expandable analysis results (no whole-panel scrolling)
        self.right_panel.pack(side="right", fill="both", expand=True, padx=0, pady=0)
        
        # Pack right panel widgets directly (no scrolling container)
        self.analysis_results_label.pack(anchor="w", padx=15, pady=(15, 5))
        self.analysis_results_container.pack(fill="both", expand=True, padx=15, pady=(0, 15))
        
        # Results section layout
        self.results_detailed_section.pack(fill="both", expand=True, padx=10, pady=10)
        
        self.results_detailed_label.pack(anchor="w", padx=10, pady=(10, 5))
        
        self.results_frame.pack(fill="both", expand=True, padx=10, pady=(0, 10))
        
        # Headers
        self.headers_frame.pack(fill="x", pady=(0, 5))
        
        self.header_residue.pack(side="left", padx=2)
        self.header_id.pack(side="left", padx=2)
        self.header_chain.pack(side="left", padx=2)
        self.header_pka.pack(side="left", padx=2)
        self.header_state_ph.pack(side="left", padx=2)
        self.header_custom.pack(side="left", padx=2)
        
        # Scrollable results
        self.results_scroll.pack(fill="both", expand=True)
    
    def _bind_mouse_wheel(self):
        """Bind mouse wheel events for scrolling in left panel and residue details only."""
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
            
            # Otherwise, traverse up to find the relevant scrollable frame
            current = widget
            scrollable_frame = None
            
            while current and current != self:
                if isinstance(current, ctk.CTkScrollableFrame):
                    # Check if this is left_scroll or results_scroll (residue details)
                    if current == getattr(self, 'left_scroll', None) or current == getattr(self, 'results_scroll', None):
                        scrollable_frame = current
                        break
                current = current.master
            
            # If we found an appropriate scrollable frame, scroll it
            if scrollable_frame:
                # Try multiple ways to access the canvas (CustomTkinter versions may vary)
                canvas = None
                for attr_name in ['_parent_canvas', '_scrollable_frame_canvas', 'canvas', '_canvas']:
                    if hasattr(scrollable_frame, attr_name):
                        canvas = getattr(scrollable_frame, attr_name)
                        if canvas and hasattr(canvas, 'yview_scroll'):
                            break
                
                if canvas:
                    try:
                        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
                    except:
                        pass  # Ignore errors
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
            
            # Otherwise, traverse up to find the relevant scrollable frame
            current = widget
            scrollable_frame = None
            
            while current and current != self:
                if isinstance(current, ctk.CTkScrollableFrame):
                    # Check if this is left_scroll or results_scroll (residue details)
                    if current == getattr(self, 'left_scroll', None) or current == getattr(self, 'results_scroll', None):
                        scrollable_frame = current
                        break
                current = current.master
            
            # If we found an appropriate scrollable frame, scroll it
            if scrollable_frame:
                # Try multiple ways to access the canvas (CustomTkinter versions may vary)
                canvas = None
                for attr_name in ['_parent_canvas', '_scrollable_frame_canvas', 'canvas', '_canvas']:
                    if hasattr(scrollable_frame, attr_name):
                        canvas = getattr(scrollable_frame, attr_name)
                        if canvas and hasattr(canvas, 'yview_scroll'):
                            break
                
                if canvas:
                    try:
                        if event.num == 4:
                            canvas.yview_scroll(-1, "units")
                        elif event.num == 5:
                            canvas.yview_scroll(1, "units")
                    except:
                        pass  # Ignore errors
            return "break"
        
        # Bind mouse wheel events to all widgets recursively
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
        
        # Start binding from the main frame
        bind_to_widget(self)
        
        # Also add a delayed binding for dynamically created widgets
        def delayed_bind():
            try:
                if not self.winfo_exists():
                    return
                # Re-bind to results_scroll and its children after they're created
                if hasattr(self, 'results_scroll'):
                    bind_to_widget(self.results_scroll)
            except:
                pass
        
        # Schedule delayed binding after 100ms to catch dynamically created widgets
        self.after(100, delayed_bind)
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            # This frame uses several self.after() calls for threading UI updates
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")
    
    def _browse_working_file(self):
        """Browse for the working PDB file to analyze."""
        file_path = filedialog.askopenfilename(
            title="Select PDB File for Propka Analysis",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
            initialdir=self.initial_directory
        )
        
        if file_path:
            self.working_file_entry.delete(0, "end")
            self.working_file_entry.insert(0, file_path)
            self.current_pdb_file = file_path
            self.analysis_file = None  # Reset analysis file when loading new PDB
            self.residue_mapping = {}  # Reset residue mapping when loading new PDB
            self.capping_used = False
            
            # Auto-generate export filename
            export_path = self._generate_export_filename(file_path)
            self.export_file_entry.delete(0, "end")
            self.export_file_entry.insert(0, export_path)
            
            # Enable analysis controls
            self._set_analysis_enabled(True)
            
            # Reset button colors when new file is loaded
            self._set_all_buttons_color("blue")
            
            if self.status_callback:
                self.status_callback(f"Loaded working file: {Path(file_path).name}")
    
    def _browse_export_file(self):
        """Browse for the export file location."""
        current_path = self.export_file_entry.get().strip()
        if not current_path and self.current_pdb_file:
            current_path = self._generate_export_filename(self.current_pdb_file)
        
        file_path = filedialog.asksaveasfilename(
            title="Save Protonated PDB File",
            defaultextension=".pdb",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
            initialfile=Path(current_path).name if current_path else "protonated.pdb",
            initialdir=str(Path(current_path).parent) if current_path else self.initial_directory
        )
        
        if file_path:
            self.export_file_entry.delete(0, "end")
            self.export_file_entry.insert(0, file_path)
            
            # Reset button colors when export path is changed
            if hasattr(self, 'analysis_run') and self.analysis_run:
                self._set_all_buttons_color("blue")
    
    def _generate_export_filename(self, input_file: str) -> str:
        """Generate export filename in the app's initial directory based on input file name."""
        input_path = Path(input_file)
        # Use the app's initial directory instead of the input file's directory
        initial_dir = Path(self.initial_directory)
        export_filename = f"{input_path.stem}_protonated_prepared.pdb"
        return str(initial_dir / export_filename)
    
    def get_export_file_path(self) -> Optional[str]:
        """Get the current export file path."""
        export_path = self.export_file_entry.get().strip()
        return export_path if export_path else None
    
    def _on_cap_toggle(self):
        """Handle capping checkbox toggle - update headers immediately."""
        if self.cap_protein_var.get():
            # When capping is enabled, prepare for mapping display
            self.header_id.configure(text="ID (Orig-New)", width=120)
        else:
            # When capping is disabled, reset to normal display
            self.header_id.configure(text="ID", width=60)
            # Clear any existing residue mapping
            self.residue_mapping = {}
            self.capping_used = False
        
        # Track option change
        self._on_option_change()
        
        # If there are existing results, refresh the display to show/hide mapping
        if hasattr(self, 'protonable_residues') and self.protonable_residues:
            current_ph = float(self.ph_entry.get().strip()) if self.ph_entry.get().strip() else 7.0
            self._refresh_results_display(current_ph)
    
    def _on_option_change(self, event=None):
        """Handle when analysis options change after analysis has been run."""
        if self.analysis_run:
            self.options_changed_after_analysis = True
            self._update_button_colors()
            self._show_options_changed_warning()
        
        # Reset all buttons to blue when any option changes
        self._set_all_buttons_color("blue")
    
    def _get_current_settings(self):
        """Get current analysis settings for comparison (excluding custom states)."""
        return {
            'ph': self.ph_entry.get().strip(),
            'cap_protein': self.cap_protein_var.get(),
            'auto_detect': self.auto_detect_var.get(),
            'distance': self.distance_entry.get().strip()
            # Note: custom_states excluded since they don't require re-analysis
        }
    
    def _update_button_colors(self):
        """Update button colors based on current state and recent actions."""
        if self.options_changed_after_analysis:
            # Change all buttons to red to indicate re-analysis needed
            self._set_all_buttons_color("red")
        # If options haven't changed after analysis, leave buttons as they are
        # (they should be in their success/failure/initial states)
    
    def _set_all_buttons_color(self, color_state: str):
        """Set all main action buttons to a specific color state."""
        color_map = {
            "blue": "#1f6aa5",      # Theme blue for initial/default state
            "green": "#10B981",     # Success green
            "red": "#EF4444",       # Error/warning red
            "gray": "#6B7280"       # Disabled/no results gray
        }
        
        color = color_map.get(color_state, color_map["blue"])
        
        # Update all main action buttons
        if hasattr(self, 'run_button'):
            self.run_button.configure(fg_color=color)
        if hasattr(self, 'detect_button'):
            self.detect_button.configure(fg_color=color)
        if hasattr(self, 'apply_button'):
            self.apply_button.configure(fg_color=color)
        if hasattr(self, 'export_button'):
            self.export_button.configure(fg_color=color)
        if hasattr(self, 'reset_button'):
            self.reset_button.configure(fg_color=color)
    
    def _set_button_color(self, button_name: str, color_state: str):
        """Set a specific button to a specific color state."""
        color_map = {
            "blue": "#1f6aa5",      # Theme blue for initial/default state
            "green": "#10B981",     # Success green
            "red": "#EF4444",       # Error/warning red
            "gray": "#6B7280"       # Disabled/no results gray
        }
        
        color = color_map.get(color_state, color_map["blue"])
        
        button_attr = f"{button_name}_button"
        if hasattr(self, button_attr):
            button = getattr(self, button_attr)
            button.configure(fg_color=color)
    
    def _show_options_changed_warning(self):
        """Show warning message in residue details when options change."""
        if hasattr(self, 'scrollable_frame') and self.options_changed_after_analysis:
            # Add warning message to the top of the residue details
            self._add_warning_message()
    
    def _add_warning_message(self):
        """Add a warning message to the residue details area."""
        if hasattr(self, 'warning_frame'):
            # Remove existing warning
            self.warning_frame.destroy()
        
        # Create warning frame
        self.warning_frame = ctk.CTkFrame(
            self.results_scroll,
            fg_color="#FFF3CD",  # Light yellow background
            corner_radius=5
        )
        self.warning_frame.pack(fill="x", padx=5, pady=(5, 10))
        
        # Warning icon and text
        warning_text = "⚠️ Analysis options have changed. Please run analysis again to update results."
        self.warning_label = ctk.CTkLabel(
            self.warning_frame,
            text=warning_text,
            font=FONTS['body'],
            text_color="#856404",  # Dark yellow text
            wraplength=400
        )
        self.warning_label.pack(pady=10, padx=10)
    
    def _clear_warning_message(self):
        """Clear the warning message from residue details."""
        if hasattr(self, 'warning_frame'):
            self.warning_frame.destroy()
            delattr(self, 'warning_frame')
    
    def _reset_option_tracking(self):
        """Reset option change tracking after analysis is run."""
        self.analysis_run = True
        self.options_changed_after_analysis = False
        self.last_analysis_settings = self._get_current_settings()
        self._update_button_colors()
        self._clear_warning_message()
    
    def _run_analysis(self):
        """Run Propka analysis on the current PDB file."""
        if not self.current_pdb_file:
            messagebox.showwarning("No File", "Please select a working file first.")
            return
        
        # Validate inputs
        ph_str = self.ph_entry.get().strip()
        cap_protein = self.cap_protein_var.get()
        
        valid, error_msg = self.validator.validate_propka_inputs(
            self.current_pdb_file, ph_str, "3"
        )
        
        if not valid:
            messagebox.showerror("Validation Error", error_msg)
            return
        
        # Reset option change tracking since we're running new analysis
        self._reset_option_tracking()
        
        # Run analysis in background thread
        self._run_analysis_thread(ph_str, cap_protein)
    
    def _run_analysis_thread(self, ph_str: str, cap_protein: bool = False):
        """Run Propka analysis in a background thread."""
        def analysis_worker():
            try:
                if self.status_callback:
                    status_msg = "Running Propka analysis"
                    if cap_protein:
                        status_msg += " with protein capping"
                    self.status_callback(f"{status_msg}...")
                
                # Determine which file to use for analysis
                analysis_file = self.current_pdb_file
                
                # Perform capping if requested
                if cap_protein:
                    try:
                        if self.status_callback:
                            self.status_callback("Removing hydrogens and adding caps...")
                        
                        # Create capped version in the target output directory
                        if not self.current_pdb_file:
                            if self.status_callback:
                                self.status_callback("Error: No PDB file selected")
                            return
                        
                        pdb_path = Path(self.current_pdb_file)
                        export_path = Path(self.export_file_entry.get().strip())
                        target_dir = export_path.parent if export_path else Path(self.initial_directory)
                        
                        capped_file = target_dir / f"{pdb_path.stem}_capped{pdb_path.suffix}"
                        
                        analysis_file, self.residue_mapping = self.protein_capper.remove_hydrogens_and_cap(
                            self.current_pdb_file, capped_file, target_dir
                        )
                        
                        self.capping_used = True
                        
                        logger.info(f"Protein capping completed: {analysis_file}")
                        logger.debug(f"Residue mapping: {len(self.residue_mapping)} entries")
                        
                        # Inform user about mapping file creation
                        export_path = Path(self.export_file_entry.get().strip()) if self.export_file_entry.get().strip() else None
                        target_dir = export_path.parent if export_path else Path(self.initial_directory)
                        mapping_filename = f"{capped_file.stem}_gatewizard_residue_mapping.txt"
                        mapping_file = target_dir / mapping_filename
                        if self.status_callback and mapping_file.exists():
                            self.status_callback(f"✅ GateWizard residue mapping file created: {mapping_file.name}")
                        
                    except ProteinCappingError as e:
                        try:
                            if self.winfo_exists():
                                self.after(0, self._show_analysis_error, f"Protein capping failed: {str(e)}")
                        except:
                            pass
                        return
                    except Exception as e:
                        logger.error(f"Unexpected error during protein capping: {e}", exc_info=True)
                        try:
                            if self.winfo_exists():
                                self.after(0, self._show_analysis_error, f"Protein capping failed: {str(e)}")
                        except:
                            pass
                        return
                else:
                    # Clear residue mapping when not capping
                    self.residue_mapping = {}
                    self.capping_used = False
                
                # Store the analysis file for later use in apply_protonation_states
                self.analysis_file = analysis_file
                
                if self.status_callback:
                    self.status_callback("Running Propka analysis...")
                
                # Set propka version to 3 (the only available version)
                self.preparation_manager.propka_version = "3"
                
                # Determine target directory for output files
                export_path = Path(self.export_file_entry.get().strip()) if self.export_file_entry.get().strip() else None
                target_dir = export_path.parent if export_path else Path(self.initial_directory)
                
                # Run analysis on the (possibly capped) file with target directory
                if not analysis_file:
                    if self.status_callback:
                        self.status_callback("Error: Analysis file is not available")
                    return
                pka_file = self.preparation_manager.run_analysis(analysis_file, output_dir=str(target_dir))
                
                # Extract summary to target directory with custom filename
                summary_file = self.preparation_manager.extract_summary(pka_file, output_dir=str(target_dir))
                
                # Parse results
                residues = self.preparation_manager.parse_summary(summary_file)
                
                # Update GUI in main thread
                try:
                    if self.winfo_exists():
                        self.after(0, self._update_results, residues, float(ph_str))
                except:
                    pass
                
            except PreparationError as e:
                try:
                    if self.winfo_exists():
                        self.after(0, self._show_analysis_error, str(e))
                except:
                    pass
            except Exception as e:
                logger.error(f"Unexpected error in Propka analysis: {e}", exc_info=True)
                try:
                    if self.winfo_exists():
                        self.after(0, self._show_analysis_error, f"Unexpected error: {str(e)}")
                except:
                    pass
        
        # Start analysis thread
        analysis_thread = threading.Thread(target=analysis_worker, daemon=True)
        analysis_thread.start()
    
    def _update_headers_for_capping(self):
        """Update column headers to reflect whether capping is enabled."""
        if self.cap_protein_var.get() or self.capping_used:
            # When capping is enabled or was used in current analysis, show both original and current IDs
            self.header_id.configure(text="ID (Orig-New)", width=120)
        else:
            # When capping is disabled, show just the ID
            self.header_id.configure(text="ID", width=60)
    
    def _get_residue_id_display(self, residue: Dict[str, Any]) -> str:
        """
        Get the display text for residue ID, showing both original and current when capping is enabled.
        
        Args:
            residue: Residue dictionary with res_id, residue, and chain info
            
        Returns:
            Formatted string for ID display
        """
        current_resid = residue['res_id']
        chain = residue.get('chain', '')
        resname = residue['residue']
        
        if not self.residue_mapping:
            # No capping - just show current ID
            return str(current_resid)
        
        # Look for original residue ID in the mapping
        original_resid = None
        
        # First try exact match with chain
        for (orig_resname, orig_chain, orig_resid), (new_resname, new_chain, new_resid) in self.residue_mapping.items():
            if new_resid == current_resid and new_chain == chain and new_resname == resname:
                original_resid = orig_resid
                break
        
        # If no exact match, try without chain (in case of chain ID mismatch)
        if original_resid is None:
            for (orig_resname, orig_chain, orig_resid), (new_resname, new_chain, new_resid) in self.residue_mapping.items():
                if new_resid == current_resid and new_resname == resname:
                    original_resid = orig_resid
                    break
        
        if original_resid is not None:
            # Show original → current format
            return f"{original_resid}→{current_resid}"
        else:
            # Fallback to just current ID (might be ACE/NME caps or mapping issue)
            return str(current_resid)
    
    def _update_results(self, residues: List[Dict[str, Any]], ph: float):
        """Update the results display with analysis results."""
        try:
            self.protonable_residues = residues
            self.custom_states = {}
            
            # Store timestamp of when analysis was completed
            self.analysis_timestamp = datetime.now()
            
            # Update headers based on capping status
            self._update_headers_for_capping()
            
            # Reset sort state
            self.sort_column = None
            self.sort_ascending = True
            
            # Clear existing results and refresh display
            self._refresh_results_display(ph)
            
            # Update summary
            self._update_summary(residues, ph)
            
            # Reset sort indicators
            self._update_sort_indicators()
            
            # Enable action buttons
            self.export_button.configure(state="normal")
            self.apply_button.configure(state="normal")
            self.reset_button.configure(state="normal")
            
            if self.status_callback:
                self.status_callback(f"Analysis complete: {len(residues)} protonable residues")
            
            # Set Run Analysis button to green on success
            self._set_button_color("run", "green")
            
            logger.info(f"Propka analysis completed: {len(residues)} residues analyzed")
            
        except Exception as e:
            logger.error(f"Error updating results: {e}", exc_info=True)
            self._show_analysis_error(f"Error displaying results: {str(e)}")
    
    def _create_result_row(self, residue: Dict[str, Any], ph: float, row_index: int, parent_frame):
        """Create a row in the results display for a residue."""
        row_frame = ctk.CTkFrame(parent_frame, fg_color="transparent")
        row_frame.pack(fill="x", pady=1)
        
        # Residue name
        residue_label = ctk.CTkLabel(
            row_frame,
            text=residue['residue'],
            font=FONTS['body'],
            width=80
        )
        residue_label.pack(side="left", padx=2)
        
        # Residue ID (show original→current when capping is enabled)
        id_display = self._get_residue_id_display(residue)
        id_label = ctk.CTkLabel(
            row_frame,
            text=id_display,
            font=FONTS['body'],
            width=120 if (self.cap_protein_var.get() or self.capping_used) else 60
        )
        id_label.pack(side="left", padx=2)
        
        # Chain column
        chain_text = residue.get('chain', '') or '-'
        chain_label = ctk.CTkLabel(
            row_frame,
            text=chain_text,
            font=FONTS['body'],
            width=60
        )
        chain_label.pack(side="left", padx=2)
        
        # pKa value
        pka_label = ctk.CTkLabel(
            row_frame,
            text=f"{residue['pka']:.2f}",
            font=FONTS['body'],
            width=80
        )
        pka_label.pack(side="left", padx=2)
        
        # Default state at pH
        default_state = self.preparation_manager.get_default_protonation_state(residue, ph)
        state_label = ctk.CTkLabel(
            row_frame,
            text=default_state,
            font=FONTS['body'],
            width=100
        )
        state_label.pack(side="left", padx=2)
        
        # Custom state selection
        available_states = self.preparation_manager.get_available_states(residue['residue'])
        
        if available_states:
            # Create unique identifier for this residue (same format as in _on_custom_state_changed)
            chain = residue.get('chain', '')
            if chain:
                residue_id = f"{residue['residue']}{residue['res_id']}_{chain}"
            else:
                residue_id = f"{residue['residue']}{residue['res_id']}"
            
            state_combo = ctk.CTkComboBox(
                row_frame,
                values=list(available_states.values()),
                width=120,
                height=25,
                command=lambda value, r=residue: self._on_custom_state_changed(r, value)
            )
            
            # Set current value (custom state if exists, otherwise default)
            current_custom_state = self.custom_states.get(residue_id)
            if current_custom_state:
                state_combo.set(current_custom_state)
            else:
                state_combo.set(default_state)
            
            # Store reference to ComboBox for color updates
            self.custom_state_combos[residue_id] = state_combo
            
            # Set initial color based on whether state is custom or default
            self._update_combo_color(residue_id, residue, ph)
            
            state_combo.pack(side="left", padx=2)
        else:
            # No alternative states available
            no_states_label = ctk.CTkLabel(
                row_frame,
                text="No alternatives",
                font=FONTS['small'],
                width=120,
                text_color=COLOR_SCHEME['inactive']
            )
            no_states_label.pack(side="left", padx=2)
        
        # Bind mouse wheel events to all widgets in this row
        self._bind_mousewheel_to_widget(row_frame)
    
    def _update_combo_color(self, residue_id: str, residue: Dict[str, Any], ph: float):
        """Update ComboBox color based on whether it has a custom state."""
        if residue_id not in self.custom_state_combos:
            return
        
        combo = self.custom_state_combos[residue_id]
        current_value = combo.get()
        default_state = self.preparation_manager.get_default_protonation_state(residue, ph)
        
        # Set color based on whether current value differs from default
        if current_value != default_state:
            # Custom state - use orange/amber color to indicate modification
            combo.configure(
                fg_color="#F59E0B",      # Amber color for modified state
                button_color="#D97706",  # Darker amber for button
                text_color="white"       # White text for contrast
            )
        else:
            # Default state - use available color scheme colors that exist
            # Use a neutral color scheme that should work with the dark theme
            combo.configure(
                fg_color=COLOR_SCHEME['content_inside_bg'],  # Use existing color key
                button_color=COLOR_SCHEME['buttons'],        # Use existing color key
                text_color=COLOR_SCHEME['text']              # Use existing color key
            )
    
    def _update_all_combo_colors(self, ph: float):
        """Update colors for all ComboBox widgets based on current custom states."""
        for residue_id, combo in self.custom_state_combos.items():
            # Find the corresponding residue data
            for residue in self.protonable_residues:
                chain = residue.get('chain', '')
                if chain:
                    current_residue_id = f"{residue['residue']}{residue['res_id']}_{chain}"
                else:
                    current_residue_id = f"{residue['residue']}{residue['res_id']}"
                
                if current_residue_id == residue_id:
                    self._update_combo_color(residue_id, residue, ph)
                    break

    def _bind_mousewheel_to_widget(self, widget):
        """Bind mouse wheel events to a widget and its children recursively."""
        def _on_mousewheel(event):
            # Find the results_scroll frame
            if hasattr(self, 'results_scroll'):
                scrollable_frame = self.results_scroll
                # Try multiple ways to access the canvas
                canvas = None
                for attr_name in ['_parent_canvas', '_scrollable_frame_canvas', 'canvas', '_canvas']:
                    if hasattr(scrollable_frame, attr_name):
                        canvas = getattr(scrollable_frame, attr_name)
                        if canvas and hasattr(canvas, 'yview_scroll'):
                            break
                
                if canvas:
                    try:
                        canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
                    except:
                        pass
            return "break"
        
        def _on_mousewheel_linux(event):
            # Find the results_scroll frame
            if hasattr(self, 'results_scroll'):
                scrollable_frame = self.results_scroll
                # Try multiple ways to access the canvas
                canvas = None
                for attr_name in ['_parent_canvas', '_scrollable_frame_canvas', 'canvas', '_canvas']:
                    if hasattr(scrollable_frame, attr_name):
                        canvas = getattr(scrollable_frame, attr_name)
                        if canvas and hasattr(canvas, 'yview_scroll'):
                            break
                
                if canvas:
                    try:
                        if event.num == 4:
                            canvas.yview_scroll(-1, "units")
                        elif event.num == 5:
                            canvas.yview_scroll(1, "units")
                    except:
                        pass
            return "break"
        
        try:
            # Bind to the widget itself
            widget.bind("<MouseWheel>", _on_mousewheel)
            widget.bind("<Button-4>", _on_mousewheel_linux)
            widget.bind("<Button-5>", _on_mousewheel_linux)
            
            # Recursively bind to all children
            for child in widget.winfo_children():
                self._bind_mousewheel_to_widget(child)
        except Exception:
            pass  # Some widgets may not support binding
    
    def _sort_results(self, column: str):
        """Sort the results by the specified column."""
        if not hasattr(self, 'protonable_residues') or not self.protonable_residues:
            return
        
        # Toggle sort direction if same column, otherwise ascending
        if self.sort_column == column:
            self.sort_ascending = not self.sort_ascending
        else:
            self.sort_ascending = True
        
        self.sort_column = column
        
        # Sort the data
        try:
            if column == 'residue':
                sorted_residues = sorted(
                    self.protonable_residues,
                    key=lambda x: x['residue'],
                    reverse=not self.sort_ascending
                )
            elif column == 'res_id':
                sorted_residues = sorted(
                    self.protonable_residues,
                    key=lambda x: x['res_id'],
                    reverse=not self.sort_ascending
                )
            elif column == 'chain':
                sorted_residues = sorted(
                    self.protonable_residues,
                    key=lambda x: x.get('chain', ''),
                    reverse=not self.sort_ascending
                )
            elif column == 'pka':
                sorted_residues = sorted(
                    self.protonable_residues,
                    key=lambda x: x['pka'],
                    reverse=not self.sort_ascending
                )
            elif column == 'state_at_ph':
                ph = float(self.ph_entry.get())
                sorted_residues = sorted(
                    self.protonable_residues,
                    key=lambda x: self.preparation_manager.get_default_protonation_state(x, ph),
                    reverse=not self.sort_ascending
                )
            else:
                sorted_residues = self.protonable_residues
            
            # Update the display with sorted data
            self.protonable_residues = sorted_residues
            ph = float(self.ph_entry.get())
            self._refresh_results_display(ph)
            
            # Update header buttons to show sort direction
            self._update_sort_indicators()
            
        except Exception as e:
            logger.error(f"Error sorting results: {e}")
    
    def _update_sort_indicators(self):
        """Update the sort indicators in column headers."""
        # Reset all headers to default
        headers = {
            'residue': self.header_residue,
            'res_id': self.header_id,
            'chain': self.header_chain,
            'pka': self.header_pka,
            'state_at_ph': self.header_state_ph
        }
        
        # Get current base texts (dynamic for ID column based on capping)
        base_texts = {
            'residue': 'Residue',
            'res_id': 'ID (Orig-New)' if (self.cap_protein_var.get() or self.capping_used) else 'ID',
            'chain': 'Chain',
            'pka': 'pKa',
            'state_at_ph': 'State at pH'
        }
        
        for col, header in headers.items():
            if col == self.sort_column:
                arrow = " ▲" if self.sort_ascending else " ▼"
                header.configure(text=f"{base_texts[col]}{arrow}")
            else:
                header.configure(text=base_texts[col])
    
    def _refresh_results_display(self, ph: float):
        """Refresh the results display without changing the data order."""
        # Clear existing results
        for widget in self.results_scroll.winfo_children():
            widget.destroy()
        
        # Clear ComboBox references since widgets are being destroyed
        self.custom_state_combos = {}
        
        # Recreate result rows with current order
        for i, residue in enumerate(self.protonable_residues):
            self._create_result_row(residue, ph, i, self.results_scroll)
    
    def _on_custom_state_changed(self, residue: Dict[str, Any], new_state: str):
        """Handle custom state selection change."""
        # Create unique identifier using the same format as core module
        chain = residue.get('chain', '')
        if chain:
            residue_id = f"{residue['residue']}{residue['res_id']}_{chain}"  # ASP80_A format
        else:
            residue_id = f"{residue['residue']}{residue['res_id']}"
        
        # Get default state
        ph = float(self.ph_entry.get())
        default_state = self.preparation_manager.get_default_protonation_state(residue, ph)
        
        if new_state == default_state:
            # Remove custom state if it matches default
            self.custom_states.pop(residue_id, None)
        else:
            # Store custom state
            self.custom_states[residue_id] = new_state
        
        # Update the ComboBox color to reflect the change
        self._update_combo_color(residue_id, residue, ph)
        
        # For custom state changes, only reset Apply States button to blue
        # Don't call _on_option_change() which would reset all buttons
        self._set_button_color("apply", "blue")
        
        logger.debug(f"Custom state for {residue_id}: {new_state}")
    
    def _update_summary(self, residues: List[Dict[str, Any]], ph: float):
        """Update the analysis summary."""
        try:
            # Calculate statistics
            stats = {}
            for residue in residues:
                res_type = residue['residue']
                stats[res_type] = stats.get(res_type, 0) + 1
            
            # Count states at target pH
            state_counts = {}
            for residue in residues:
                state = self.preparation_manager.get_default_protonation_state(residue, ph)
                state_counts[state] = state_counts.get(state, 0) + 1
            
            # Format summary
            summary_text = f"""Analysis completed at pH {ph:.1f}

Protonable residues found: {len(residues)}

Residue types:
"""
            for res_type, count in sorted(stats.items()):
                summary_text += f"  {res_type}: {count}\n"
            
            summary_text += f"\nProtonation states at pH {ph:.1f}:\n"
            for state, count in sorted(state_counts.items()):
                summary_text += f"  {state}: {count}\n"
            
            if self.custom_states:
                summary_text += f"\nCustom states defined: {len(self.custom_states)}"
            
            # Update summary display
            self.summary_text.configure(state="normal")
            self.summary_text.delete("1.0", "end")
            self.summary_text.insert("1.0", summary_text)
            self.summary_text.configure(state="disabled")
            
        except Exception as e:
            logger.error(f"Error updating summary: {e}")
    
    def _apply_protonation_states(self):
        """Apply protonation states and disulfide bonds to the PDB file, then run pdb4amber."""
        if not self.protonable_residues:
            messagebox.showwarning("No Results", "Please run analysis first.")
            return
        
        # Determine which file to use as input (capped or original)
        input_file = self.analysis_file if self.analysis_file else self.current_pdb_file
        
        if not input_file or not Path(input_file).exists():
            messagebox.showerror("File Error", "Analysis file not found. Please run analysis again.")
            return
        
        # Get output file path from export entry
        output_file = self.export_file_entry.get().strip()
        
        if not output_file:
            messagebox.showwarning("No Export File", "Please specify an export file path.")
            return
        
        try:
            if self.status_callback:
                status_msg = "Applying protonation states and disulfide bonds"
                if input_file != self.current_pdb_file:
                    status_msg += " to capped structure"
                self.status_callback(f"{status_msg}...")
            
            # Get target pH
            ph = float(self.ph_entry.get())
            
            # The user-specified output file should be the final _prepared.pdb file
            final_output_file = output_file
            output_path = Path(final_output_file)
            
            # Create intermediate file for protonation states (without _prepared suffix)
            if output_path.name.endswith('_prepared.pdb'):
                # Remove _prepared from the name for the intermediate file
                intermediate_name = output_path.name.replace('_prepared.pdb', '.pdb')
                intermediate_file = str(output_path.parent / intermediate_name)
            else:
                # If user didn't specify _prepared, add it for the final file
                intermediate_file = output_file
                final_output_file = str(output_path.with_name(f"{output_path.stem}_prepared.pdb"))
            
            # Apply protonation states first (use the analysis file as input)
            modification_results = self.preparation_manager.apply_protonation_states(
                input_file,  # Use the analysis file (original or capped)
                intermediate_file,  # Use intermediate file for protonation states
                ph,
                self.custom_states,
                self.protonable_residues
            )
            
            # Extract counts from the result dictionary
            residue_changes = modification_results['residue_changes']
            record_changes = modification_results['record_changes']
            
            # Apply disulfide bonds if auto-detect is enabled or bonds were detected
            disulfide_changes = 0
            if self.auto_detect_var.get() or self.detected_disulfide_bonds:
                bonds_to_apply = self.detected_disulfide_bonds if self.detected_disulfide_bonds else None
                disulfide_changes = self.preparation_manager.apply_disulfide_bonds(
                    intermediate_file,  # Use the intermediate file as input
                    intermediate_file,  # Overwrite the same file
                    bonds_to_apply,
                    auto_detect=self.auto_detect_var.get()
                )
            
            total_changes = record_changes + disulfide_changes
            
            if self.status_callback:
                self.status_callback("Running pdb4amber to prepare PDB for Amber compatibility...")
            
            # Run pdb4amber with automatic ACE/NME cap fix
            # The fix is automatically applied when capping was used
            fix_caps = self.capping_used or self.cap_protein_var.get()
            
            try:
                pdb4amber_result = self.preparation_manager.run_pdb4amber_with_cap_fix(
                    input_pdb=intermediate_file,
                    output_pdb=final_output_file,
                    fix_caps=fix_caps
                )
                
                if pdb4amber_result['success']:
                    hetatm_fixed = pdb4amber_result['hetatm_fixed']
                    
                    if self.status_callback:
                        status_msg = "pdb4amber processing completed successfully"
                        if hetatm_fixed > 0:
                            status_msg += f" (fixed {hetatm_fixed} cap records)"
                        self.status_callback(f"[OK] {status_msg}")
                    
                    success_message = (
                        f"Modifications applied successfully!\n\n"
                        f"Residue state changes: {residue_changes}\n"
                        f"PDB record changes: {record_changes}\n"
                        f"Disulfide bonds: {disulfide_changes}\n"
                        f"Total PDB changes: {total_changes}\n"
                    )
                    
                    if hetatm_fixed > 0:
                        success_message += f"ACE/NME HETATM fixes: {hetatm_fixed}\n"
                    
                    success_message += (
                        f"\nOriginal modified file: {intermediate_file}\n"
                        f"Amber-ready file: {final_output_file}\n\n"
                        f"The Amber-ready file can be used directly into the Builder tab or Packmol-memgen."
                    )
                    
                    messagebox.showinfo("Success", success_message)
                    
                    # Set Apply States button to green on success
                    self._set_button_color("apply", "green")
                    
                    # Reset option change tracking since states have been applied
                    self._reset_option_tracking()
                    
                    logger.info(f"Applied {residue_changes} residue state changes ({record_changes} PDB records), {disulfide_changes} disulfide bonds, and pdb4amber processing")
                    
            except Exception as pdb4amber_error:
                # pdb4amber failed, but still show success for the propka modifications
                logger.error(f"pdb4amber failed: {pdb4amber_error}")
                
                if self.status_callback:
                    self.status_callback(f"Applied {residue_changes} residue changes ({record_changes} PDB records) and {disulfide_changes} disulfide bonds (pdb4amber failed)")
                
                messagebox.showwarning(
                    "Partial Success",
                    f"Protonation modifications applied successfully!\n\n"
                    f"Residue state changes: {residue_changes}\n"
                    f"PDB record changes: {record_changes}\n"
                    f"Disulfide bonds: {disulfide_changes}\n"
                    f"Total PDB changes: {total_changes}\n"
                    f"Output file: {intermediate_file}\n\n"
                    f"⚠️ Warning: pdb4amber processing failed.\n"
                    f"The modified PDB file is available but may need manual preparation for MD simulations."
                )
                
                # Set Apply States button to red due to pdb4amber failure
                self._set_button_color("apply", "red")
                
                # Reset option change tracking since states have been applied
                self._reset_option_tracking()
                
                logger.warning(f"Applied {residue_changes} residue state changes ({record_changes} PDB records) and {disulfide_changes} disulfide bonds, but pdb4amber failed")
            
        except Exception as e:
            logger.error(f"Error applying modifications: {e}")
            messagebox.showerror("Error", f"Failed to apply modifications: {str(e)}")
            # Set Apply States button to red on failure
            self._set_button_color("apply", "red")
    
    def _export_results(self):
        """Export analysis results to a file."""
        if not self.protonable_residues:
            messagebox.showwarning("No Results", "Please run analysis first.")
            return
        
        # Get export file path
        if self.current_pdb_file:
            current_path = Path(self.current_pdb_file)
            default_name = f"{current_path.stem}_propka_results.txt"
            initial_dir = str(current_path.parent)
        else:
            default_name = "propka_results.txt"
            initial_dir = self.initial_directory
        
        export_file = filedialog.asksaveasfilename(
            title="Export Propka Results",
            defaultextension=".txt",
            filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All files", "*.*")],
            initialfile=default_name,
            initialdir=initial_dir
        )
        
        if not export_file:
            return
        
        try:
            if self.status_callback:
                self.status_callback("Exporting results...")
            
            ph = float(self.ph_entry.get())
            
            # Format analysis timestamp
            if self.analysis_timestamp:
                timestamp_str = self.analysis_timestamp.strftime("%Y-%m-%d %H:%M:%S")
            else:
                timestamp_str = "Unknown"
            
            # Prepare export data
            export_text = f"Propka Analysis Results\n"
            export_text += f"=" * 80 + "\n\n"
            export_text += f"PDB File: {self.current_pdb_file}\n"
            export_text += f"Analysis Date: {timestamp_str}\n"
            export_text += f"Target pH: {ph:.1f}\n"
            
            # Add information about structural modifications
            export_text += f"\nStructural Modifications:\n"
            export_text += f"  Protein Capping (ACE/NME): {'Yes' if self.capping_used else 'No'}\n"
            
            if self.capping_used and self.residue_mapping:
                export_text += f"  Residues affected by capping: {len(self.residue_mapping)}\n"
            
            if self.detected_disulfide_bonds:
                export_text += f"  Disulfide Bonds Detected: {len(self.detected_disulfide_bonds)}\n"
                for i, bond in enumerate(self.detected_disulfide_bonds, 1):
                    (res1, id1), (res2, id2) = bond
                    export_text += f"    {i}. {res1}{id1} ↔ {res2}{id2}\n"
            else:
                export_text += f"  Disulfide Bonds: None detected\n"
            
            export_text += f"\n" + "=" * 80 + "\n\n"
            export_text += "Residue\tID\tChain\tpKa\tDefault_State\tCustom_State\n"
            
            for residue in self.protonable_residues:
                # Use the same ID format as in _on_custom_state_changed
                chain = residue.get('chain', '')
                if chain:
                    residue_id = f"{residue['residue']}{residue['res_id']}_{chain}"
                else:
                    residue_id = f"{residue['residue']}{residue['res_id']}"
                
                # Get ID display (showing original→current if capped)
                id_display = self._get_residue_id_display(residue)
                chain_display = chain if chain else '-'
                    
                default_state = self.preparation_manager.get_default_protonation_state(residue, ph)
                custom_state = self.custom_states.get(residue_id, "")
                
                export_text += f"{residue['residue']}\t{id_display}\t{chain_display}\t{residue['pka']:.2f}\t{default_state}\t{custom_state}\n"
            
            # Add summary statistics
            export_text += f"\n" + "=" * 80 + "\n"
            export_text += f"Summary:\n"
            export_text += f"  Total Protonable Residues: {len(self.protonable_residues)}\n"
            if self.custom_states:
                export_text += f"  Custom States Defined: {len(self.custom_states)}\n"
            
            # Write to file
            with open(export_file, 'w') as f:
                f.write(export_text)
            
            if self.status_callback:
                self.status_callback(f"Results exported: {Path(export_file).name}")
            
            messagebox.showinfo("Export Complete", f"Results exported to:\n{export_file}")
            
            # Set Export Results button to green on success
            self._set_button_color("export", "green")
            
        except Exception as e:
            logger.error(f"Error exporting results: {e}")
            messagebox.showerror("Export Error", f"Failed to export results: {str(e)}")
            # Set Export Results button to red on failure
            self._set_button_color("export", "red")
    
    def _reset_custom_states(self):
        """Reset all custom protonation states."""
        if not self.custom_states:
            messagebox.showinfo("No Custom States", "No custom states to reset.")
            self._set_button_color("reset", "gray")
            return
        
        result = messagebox.askyesno(
            "Reset Custom States",
            f"This will reset {len(self.custom_states)} custom protonation states to defaults.\n\n"
            "Are you sure?"
        )
        
        if result:
            try:
                self.custom_states = {}
                
                # Update display preserving sort order
                if self.protonable_residues:
                    ph = float(self.ph_entry.get())
                    self._refresh_results_display(ph)
                
                if self.status_callback:
                    self.status_callback("Custom states reset")
                
                # Set Reset Custom button to green on success
                self._set_button_color("reset", "green")
                
            except Exception as e:
                logger.error(f"Error resetting custom states: {e}")
                messagebox.showerror("Error", f"Failed to reset custom states: {str(e)}")
                # Set Reset Custom button to red on failure
                self._set_button_color("reset", "red")
    
    def _show_analysis_error(self, error_message: str):
        """Show analysis error message."""
        logger.error(f"Propka analysis error: {error_message}")
        messagebox.showerror("Analysis Error", f"Propka analysis failed:\n\n{error_message}")
        
        # Set Run Analysis button to red on failure
        self._set_button_color("run", "red")
        
        if self.status_callback:
            self.status_callback("Analysis failed")
    
    def _set_analysis_enabled(self, enabled: bool):
        """Enable or disable analysis controls."""
        state = "normal" if enabled else "disabled"
        
        self.run_button.configure(state=state)
        # Keep pH entry always enabled so users can change it anytime
        # self.ph_entry.configure(state=state)
        self.detect_button.configure(state=state)
    
    def _detect_disulfide_bonds(self):
        """Detect disulfide bonds in the current PDB file."""
        if not self.current_pdb_file:
            messagebox.showerror(
                "Error",
                "No PDB file loaded for disulfide bond detection"
            )
            self._set_button_color("detect", "red")
            return
        
        try:
            # Get distance threshold
            distance_threshold = float(self.distance_entry.get())
            
            # Detect disulfide bonds
            self.detected_disulfide_bonds = self.preparation_manager.detect_disulfide_bonds(
                self.current_pdb_file,
                distance_threshold
            )
            
            # Update display
            self._update_disulfide_display()
            
            # Set button color based on detection results
            if self.detected_disulfide_bonds:
                self._set_button_color("detect", "green")
                if self.status_callback:
                    self.status_callback(f"✅ Detected {len(self.detected_disulfide_bonds)} disulfide bonds")
            else:
                self._set_button_color("detect", "gray")
                if self.status_callback:
                    self.status_callback("No disulfide bonds detected")
                
        except ValueError:
            error_msg = "Invalid distance threshold. Please enter a valid number."
            messagebox.showerror("Error", error_msg)
            self._set_button_color("detect", "red")
        except Exception as e:
            error_msg = f"Error detecting disulfide bonds: {str(e)}"
            messagebox.showerror("Error", error_msg)
            self._set_button_color("detect", "red")
            logger.error(f"Error detecting disulfide bonds: {e}")
    
    def _update_disulfide_display(self):
        """Update the disulfide bonds display text."""
        self.disulfide_text.configure(state="normal")
        self.disulfide_text.delete("1.0", "end")
        
        if self.detected_disulfide_bonds:
            bond_text = []
            for i, bond in enumerate(self.detected_disulfide_bonds, 1):
                (res1, id1), (res2, id2) = bond
                bond_text.append(f"{i}. {res1}{id1} ↔ {res2}{id2}")
            
            self.disulfide_text.insert("1.0", "\n".join(bond_text))
        else:
            self.disulfide_text.insert("1.0", "No disulfide bonds detected")
        
        self.disulfide_text.configure(state="disabled")
    
    def on_stage_shown(self):
        """Called when this stage becomes active."""
        # Update current PDB file
        if self.get_current_pdb:
            self.current_pdb_file = self.get_current_pdb()
            self._set_analysis_enabled(self.current_pdb_file is not None)
    
    def on_pdb_changed(self, pdb_file: Optional[str]):
        """Called when PDB file changes."""
        self.current_pdb_file = pdb_file
        self._set_analysis_enabled(pdb_file is not None)
        
        # Reset button colors when PDB file changes
        if pdb_file:
            self._set_all_buttons_color("blue")
        
        # Clear previous results
        if not pdb_file:
            self.protonable_residues = []
            self.custom_states = {}
            
            # Clear results display
            for widget in self.results_scroll.winfo_children():
                widget.destroy()
            
            # Reset summary
            self.summary_text.configure(state="normal")
            self.summary_text.delete("1.0", "end")
            self.summary_text.insert("1.0", "No analysis performed yet")
            self.summary_text.configure(state="disabled")
            
            # Disable action buttons
            self.export_button.configure(state="disabled")
            self.apply_button.configure(state="disabled")
            self.reset_button.configure(state="disabled")
    
    def cleanup(self):
        """Cleanup resources when frame is destroyed."""
        # No specific cleanup needed for Propka frame
        pass
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the Propka frame."""
        try:
            # File section
            if hasattr(self, 'file_label'):
                self.file_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'working_file_label'):
                self.working_file_label.configure(font=scaled_fonts['body'])
            if hasattr(self, 'working_file_entry'):
                self.working_file_entry.configure(font=scaled_fonts['body'])
            if hasattr(self, 'browse_file_button'):
                self.browse_file_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'export_file_label'):
                self.export_file_label.configure(font=scaled_fonts['body'])
            if hasattr(self, 'export_file_entry'):
                self.export_file_entry.configure(font=scaled_fonts['body'])
            if hasattr(self, 'export_browse_button'):
                self.export_browse_button.configure(font=scaled_fonts['body'])

            # Analysis section
            if hasattr(self, 'analysis_label'):
                self.analysis_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'ph_label'):
                self.ph_label.configure(font=scaled_fonts['body'])
            if hasattr(self, 'ph_entry'):
                self.ph_entry.configure(font=scaled_fonts['body'])
            if hasattr(self, 'run_button'):
                self.run_button.configure(font=scaled_fonts['body'])

            # Disulfide bonds section
            if hasattr(self, 'disulfide_label'):
                self.disulfide_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'auto_detect_check'):
                self.auto_detect_check.configure(font=scaled_fonts['body'])
            if hasattr(self, 'distance_label'):
                self.distance_label.configure(font=scaled_fonts['body'])
            if hasattr(self, 'distance_entry'):
                self.distance_entry.configure(font=scaled_fonts['body'])
            if hasattr(self, 'detect_button'):
                self.detect_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'disulfide_list_label'):
                self.disulfide_list_label.configure(font=scaled_fonts['body'])
            if hasattr(self, 'disulfide_text'):
                try:
                    self.disulfide_text.configure(font=scaled_fonts['small'])
                except:
                    pass

            # Results section
            if hasattr(self, 'results_label'):
                try:
                    self.results_label.configure(font=scaled_fonts['heading'])  # type: ignore
                except Exception:
                    pass
            if hasattr(self, 'header_residue'):
                self.header_residue.configure(font=scaled_fonts['body'])
            if hasattr(self, 'header_id'):
                self.header_id.configure(font=scaled_fonts['body'])
            if hasattr(self, 'header_pka'):
                self.header_pka.configure(font=scaled_fonts['body'])
            if hasattr(self, 'header_state_ph'):
                self.header_state_ph.configure(font=scaled_fonts['body'])
            if hasattr(self, 'header_custom'):
                self.header_custom.configure(font=scaled_fonts['body'])
            # Results scrollable area: children are added dynamically, not handled here

            # Action buttons
            if hasattr(self, 'apply_button'):
                self.apply_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'export_button'):
                self.export_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'reset_button'):
                self.reset_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'save_button'):
                try:
                    self.save_button.configure(font=scaled_fonts['body'])  # type: ignore
                except Exception:
                    pass

            # Status label
            if hasattr(self, 'status_label'):
                try:
                    self.status_label.configure(font=scaled_fonts['small'])  # type: ignore
                except Exception:
                    pass

            # Main result text (black box)
            if hasattr(self, 'result_text'):
                try:
                    self.result_text.configure(font=scaled_fonts['small'])  # type: ignore
                except Exception:
                    pass

        except Exception as e:
            logger.warning(f"Error updating fonts in PreparationFrame: {e}")