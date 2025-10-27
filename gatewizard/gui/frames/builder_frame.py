# gatewizard/gui/frames/builder_frame.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Builder frame for membrane protein system building.

This module provides the GUI for setting up membrane protein systems
using packmol-memgen and related tools.
"""

import tkinter as tk
from tkinter import messagebox, filedialog
from typing import Optional, Callable, List
from pathlib import Path
import os

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, WIDGET_SIZES, LAYOUT, ERROR_MESSAGES, SUCCESS_MESSAGES
)
from gatewizard.gui.widgets.leaflet_frame import LeafletFrame
from gatewizard.gui.widgets.progress_tracker import ProgressTracker
from gatewizard.gui.widgets.searchable_combobox import SearchableComboBox
from gatewizard.core.builder import Builder, BuilderError
from gatewizard.tools.force_fields import ForceFieldManager
from gatewizard.tools.validators import SystemValidator
from gatewizard.utils.config import set_working_directory
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class BuilderFrame(ctk.CTkFrame):
    """
    Frame for membrane protein system building.
    
    This frame handles the setup and execution of membrane protein system
    building using packmol-memgen and related tools.
    """
    
    def __init__(
        self,
        parent,
        get_current_pdb: Optional[Callable[[], Optional[str]]] = None,
        status_callback: Optional[Callable[[str], None]] = None,
        initial_directory: Optional[str] = None
    ):
        """
        Initialize the preparation frame.

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

        # Initialize components
        self.system_builder = Builder()
        self.ff_manager = ForceFieldManager()
        self.validator = SystemValidator()

        # State variables
        self.current_pdb_file = None
        self.working_directory = self.initial_directory
        self.preparation_output_name = "preparation"  # Default output folder name

        # Create widgets
        self._create_widgets()
        self._setup_layout()
        self._load_defaults()
        
        # Enable mouse wheel scrolling
        self._bind_mouse_wheel()

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
            
            # Otherwise, scroll the main frame
            if hasattr(self.main_scroll, '_parent_canvas'):
                self.main_scroll._parent_canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")
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
            
            # Otherwise, scroll the main frame (Linux mouse wheel events)
            if hasattr(self.main_scroll, '_parent_canvas'):
                if event.num == 4:
                    self.main_scroll._parent_canvas.yview_scroll(-1, "units")
                elif event.num == 5:
                    self.main_scroll._parent_canvas.yview_scroll(1, "units")
            return "break"
        
        # Bind to the main frame and its children
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
            except:
                pass  # Some widgets may not support binding
        
        # Bind after a short delay to ensure widgets are created
        self.after(100, lambda: bind_to_widget(self))
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            # This frame uses self.after() calls for UI binding
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")

    def _create_widgets(self):
        """Create all widgets for the preparation frame."""
        # Create scrollable main frame
        self.main_scroll = ctk.CTkScrollableFrame(self, fg_color="transparent")
        
        # Working directory section
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
            placeholder_text="Select working directory...",
            height=WIDGET_SIZES['entry_height']
        )
        self.workdir_entry.insert(0, self.working_directory)
        
        self.workdir_browse_button = ctk.CTkButton(
            self.workdir_frame,
            text="Browse",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._browse_working_directory
        )
        
        # Working file section
        self.workfile_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.workfile_label = ctk.CTkLabel(
            self.workfile_section,
            text="Working File",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.workfile_frame = ctk.CTkFrame(self.workfile_section, fg_color="transparent")
        
        self.workfile_entry = ctk.CTkEntry(
            self.workfile_frame,
            placeholder_text="Select PDB file for preparation...",
            height=WIDGET_SIZES['entry_height']
        )
        
        self.workfile_browse_button = ctk.CTkButton(
            self.workfile_frame,
            text="Browse",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._browse_working_file
        )
        
        # Output folder name section
        self._create_outputname_section()
        
        # Lipid composition section
        self.lipids_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.lipids_label = ctk.CTkLabel(
            self.lipids_section,
            text="Lipid Composition",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Upper and lower leaflet frames
        self.leaflets_frame = ctk.CTkFrame(self.lipids_section, fg_color="transparent")
        
        self.upper_leaflet = LeafletFrame(
            self.leaflets_frame,
            title="Upper Leaflet",
            available_lipids=self.ff_manager.get_available_lipids()
        )
        
        self.lower_leaflet = LeafletFrame(
            self.leaflets_frame,
            title="Lower Leaflet", 
            available_lipids=self.ff_manager.get_available_lipids()
        )
        
        # Force field section
        self.ff_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.ff_label = ctk.CTkLabel(
            self.ff_section,
            text="Force Field Parameters",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.ff_frame = ctk.CTkFrame(self.ff_section, fg_color="transparent")
        
        # Create row frame for force field options
        self.ff_row1 = ctk.CTkFrame(self.ff_frame, fg_color="transparent")
        
        # Water model
        self.water_label = ctk.CTkLabel(
            self.ff_row1,
            text="Water Model:",
            font=FONTS['body']
        )
        
        self.water_combo = SearchableComboBox(
            self.ff_row1,
            values=self.ff_manager.get_water_models(),
            width=WIDGET_SIZES['combobox_width'],
            height=WIDGET_SIZES['combobox_height']
        )
        
        # Protein force field
        self.protein_ff_label = ctk.CTkLabel(
            self.ff_row1,
            text="Protein FF:",
            font=FONTS['body']
        )
        
        self.protein_ff_combo = SearchableComboBox(
            self.ff_row1,
            values=self.ff_manager.get_protein_force_fields(),
            width=WIDGET_SIZES['combobox_width'],
            height=WIDGET_SIZES['combobox_height']
        )
        
        # Lipid force field
        self.lipid_ff_label = ctk.CTkLabel(
            self.ff_row1,
            text="Lipid FF:",
            font=FONTS['body']
        )
        
        self.lipid_ff_combo = SearchableComboBox(
            self.ff_row1,
            values=self.ff_manager.get_lipid_force_fields(),
            width=WIDGET_SIZES['combobox_width'],
            height=WIDGET_SIZES['combobox_height']
        )
        
        # System options section
        self.options_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.options_label = ctk.CTkLabel(
            self.options_section,
            text="System Options",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.options_frame = ctk.CTkFrame(self.options_section, fg_color="transparent")
        
        # Checkboxes for options
        self.preoriented_var = ctk.BooleanVar(value=True)
        self.preoriented_checkbox = ctk.CTkCheckBox(
            self.options_frame,
            text="Protein is pre-oriented",
            variable=self.preoriented_var
        )
        
        self.parametrize_var = ctk.BooleanVar(value=True)
        self.parametrize_checkbox = ctk.CTkCheckBox(
            self.options_frame,
            text="Run parametrization with tleap",
            variable=self.parametrize_var
        )
        
        self.add_salt_var = ctk.BooleanVar(value=True)
        self.add_salt_checkbox = ctk.CTkCheckBox(
            self.options_frame,
            text="Add salt",
            variable=self.add_salt_var,
            command=self._toggle_salt_options
        )
        
        # Salt options
        self.salt_frame = ctk.CTkFrame(self.options_section, fg_color="transparent")
        
        self.salt_conc_label = ctk.CTkLabel(
            self.salt_frame,
            text="Concentration (M):",
            font=FONTS['body']
        )
        
        self.salt_conc_entry = ctk.CTkEntry(
            self.salt_frame,
            width=80,
            height=WIDGET_SIZES['entry_height']
        )
        self.salt_conc_entry.insert(0, "0.15")
        
        self.cation_label = ctk.CTkLabel(
            self.salt_frame,
            text="Cation:",
            font=FONTS['body']
        )
        
        self.cation_combo = SearchableComboBox(
            self.salt_frame,
            values=self.ff_manager.get_available_cations(),
            width=100,
            height=WIDGET_SIZES['combobox_height']
        )
        self.cation_combo.set("K+")
        
        self.anion_label = ctk.CTkLabel(
            self.salt_frame,
            text="Anion:",
            font=FONTS['body']
        )
        
        self.anion_combo = SearchableComboBox(
            self.salt_frame,
            values=self.ff_manager.get_available_anions(),
            width=100,
            height=WIDGET_SIZES['combobox_height']
        )
        self.anion_combo.set("Cl-")
        
        # Water layer distance options
        self.water_distance_frame = ctk.CTkFrame(self.options_section, fg_color="transparent")
        
        self.water_distance_label = ctk.CTkLabel(
            self.water_distance_frame,
            text="Water layer thickness (Å):",
            font=FONTS['body']
        )
        
        self.water_distance_entry = ctk.CTkEntry(
            self.water_distance_frame,
            width=80,
            height=WIDGET_SIZES['entry_height']
        )
        self.water_distance_entry.insert(0, "17.5")  # Default value from packmol-memgen
        
        self.water_distance_help = ctk.CTkLabel(
            self.water_distance_frame,
            text="Width of water layer over membrane/protein in z-axis",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )
        
        # Simplified protonation options
        self.protonation_frame = ctk.CTkFrame(self.options_section, fg_color="transparent")
        
        self.notprotonate_var = ctk.BooleanVar(value=True)
        self.notprotonate_checkbox = ctk.CTkCheckBox(
            self.protonation_frame,
            text="Skip protonation (preserve propka results)",
            variable=self.notprotonate_var
        )
        
        # Help text for skip protonation
        self.protonation_help_text = ctk.CTkTextbox(
            self.protonation_frame,
            height=60,
            font=FONTS['small']
        )
        help_text = """Prepare protocol (current):
1. Packing: packmol-memgen builds the membrane system from your selections
2. Optional: if "Run parametrization with tleap" is checked, the packed system is parametrized with AmberTools (tleap), producing .prmtop, .inpcrd and a finalized PDB
3. Protonation: if "Skip protonation (preserve propka results)" is checked, existing residue names (e.g., GLH/ASH/HIP/HIE/HID) are kept and no re-protonation is performed
4. Salt: ion addition (cation/anion and concentration) is applied during parametrization when enabled"""
        self.protonation_help_text.insert("1.0", help_text)
        self.protonation_help_text.configure(state="disabled")
        
        # Action buttons section
        self.actions_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.actions_label = ctk.CTkLabel(
            self.actions_section,
            text="Actions",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.actions_frame = ctk.CTkFrame(self.actions_section, fg_color="transparent")
        
        self.validate_button = ctk.CTkButton(
            self.actions_frame,
            text="Validate Inputs",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._validate_inputs
        )
        
        self.prepare_button = ctk.CTkButton(
            self.actions_frame,
            text="Start Preparation",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._start_preparation,
            state="disabled"
        )
        
        self.load_defaults_button = ctk.CTkButton(
            self.actions_frame,
            text="Load Defaults",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._load_defaults
        )
        
        # Progress tracking section
        self.progress_section = ctk.CTkFrame(self.main_scroll, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.progress_label = ctk.CTkLabel(
            self.progress_section,
            text="Preparation Progress",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.progress_tracker = ProgressTracker(
            self.progress_section,
            working_directory=self.working_directory
        )

    def _create_outputname_section(self):
        """Create the output folder name section."""
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
            placeholder_text="Enter output folder name (default: preparation)"
        )
        
        # Set default value
        self.outputname_entry.insert(0, self.preparation_output_name)
        
        # Bind change event
        self.outputname_entry.bind('<KeyRelease>', self._on_outputname_changed)
        self.outputname_entry.bind('<FocusOut>', self._on_outputname_changed)

    def _toggle_propka_workflow(self):
        """Toggle Propka workflow options."""
        if not hasattr(self, 'use_propka_workflow_var'):
            return
        enabled = self.use_propka_workflow_var.get()  # type: ignore

        if enabled:
            if hasattr(self, 'propka_stage_frame'):
                self.propka_stage_frame.pack(fill="x", pady=LAYOUT['padding_small'])  # type: ignore
            self._update_propka_help_text()
            # Disable conflicting options
            if hasattr(self, 'notprotonate_checkbox'):
                self.notprotonate_checkbox.configure(state="disabled")
            if hasattr(self, 'use_pdb2pqr_checkbox'):
                self.use_pdb2pqr_checkbox.configure(state="disabled")  # type: ignore
        else:
            if hasattr(self, 'propka_stage_frame'):
                self.propka_stage_frame.pack_forget()  # type: ignore
            # Re-enable options
            if hasattr(self, 'notprotonate_checkbox'):
                self.notprotonate_checkbox.configure(state="normal")
            if hasattr(self, 'use_pdb2pqr_checkbox'):
                self.use_pdb2pqr_checkbox.configure(state="normal")  # type: ignore

    def _update_propka_help_text(self):
        """Update help text based on selected stage."""
        if not hasattr(self, 'propka_stage_var') or not hasattr(self, 'propka_help_text'):
            return
        
        stage = self.propka_stage_var.get()  # type: ignore

        if stage == "stage1":
            help_text = """Stage 1 Workflow:
    1. This will pack your protein into the membrane without parametrization
    2. After completion, run Propka analysis on your original PDB file
    3. Modify residue names in the packed PDB file based on Propka results
    4. Then run Stage 2 for parametrization with the modified packed file"""
        else:  # stage2
            help_text = """Stage 2 Workflow:
    1. Select the packed PDB file with Propka-modified residue names
    2. This will parametrize the system using --notprotonate flag
    3. The system will be ready for MD simulations
    4. No lipid selection needed (system already packed)"""

        self.propka_help_text.delete("1.0", "end")  # type: ignore
        self.propka_help_text.insert("1.0", help_text)  # type: ignore

    def _on_outputname_changed(self, event=None):
        """Handle output folder name change."""
        new_name = self.outputname_entry.get().strip()
        if new_name and new_name != self.preparation_output_name:
            self.preparation_output_name = new_name
            
            if self.status_callback:
                self.status_callback(f"Output folder name set to: {new_name}")

    def _browse_packed_pdb(self):
        """Browse for packed PDB file for Stage 2."""
        file_path = filedialog.askopenfilename(
            title="Select Packed PDB File",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
            initialdir=self.working_directory or self.initial_directory
        )

        if file_path and hasattr(self, 'packed_pdb_entry'):
            self.packed_pdb_entry.delete(0, "end")  # type: ignore
            self.packed_pdb_entry.insert(0, file_path)  # type: ignore
    
    def _setup_layout(self):
        """Setup the layout of widgets."""
        # Main scrollable area
        self.main_scroll.pack(fill="both", expand=True)
        
        # Working directory section
        self.workdir_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.workdir_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.workdir_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        self.workdir_entry.pack(side="left", fill="x", expand=True, padx=(0, LAYOUT['padding_small']))
        self.workdir_browse_button.pack(side="right", padx=LAYOUT['padding_small'])
        
        # Working file section
        self.workfile_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.workfile_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.workfile_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        self.workfile_entry.pack(side="left", fill="x", expand=True, padx=(0, LAYOUT['padding_small']))
        self.workfile_browse_button.pack(side="right", padx=LAYOUT['padding_small'])
        
        # Output folder name section
        self.outputname_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.outputname_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.outputname_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        self.outputname_entry.pack(side="left", fill="x", expand=True, padx=LAYOUT['padding_small'])
        
        # Lipids section
        self.lipids_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.lipids_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.leaflets_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        self.upper_leaflet.pack(side="left", fill="both", expand=True, padx=(0, LAYOUT['padding_small']))
        self.lower_leaflet.pack(side="left", fill="both", expand=True, padx=(LAYOUT['padding_small'], 0))
        
        # Force field section
        self.ff_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.ff_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.ff_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        # First row of force field options
        self.ff_row1.pack(fill="x", pady=LAYOUT['padding_small'])
        
        self.water_label.pack(side="left", padx=(0, LAYOUT['padding_small']))
        self.water_combo.pack(side="left", padx=LAYOUT['padding_small'])
        
        self.protein_ff_label.pack(side="left", padx=(LAYOUT['padding_large'], LAYOUT['padding_small']))
        self.protein_ff_combo.pack(side="left", padx=LAYOUT['padding_small'])
        
        self.lipid_ff_label.pack(side="left", padx=(LAYOUT['padding_large'], LAYOUT['padding_small']))
        self.lipid_ff_combo.pack(side="left", padx=LAYOUT['padding_small'])
        
        # Options section
        self.options_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.options_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.options_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        self.preoriented_checkbox.pack(anchor="w", pady=LAYOUT['padding_small'])
        self.parametrize_checkbox.pack(anchor="w", pady=LAYOUT['padding_small'])
        self.add_salt_checkbox.pack(anchor="w", pady=LAYOUT['padding_small'])
        
        # Salt options
        self.salt_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        self.salt_conc_label.pack(side="left", padx=(LAYOUT['padding_large'], LAYOUT['padding_small']))
        self.salt_conc_entry.pack(side="left", padx=LAYOUT['padding_small'])
        
        self.cation_label.pack(side="left", padx=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        self.cation_combo.pack(side="left", padx=LAYOUT['padding_small'])
        
        self.anion_label.pack(side="left", padx=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        self.anion_combo.pack(side="left", padx=LAYOUT['padding_small'])
        
        # Water distance options
        self.water_distance_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        self.water_distance_label.pack(side="left", padx=(LAYOUT['padding_large'], LAYOUT['padding_small']))
        self.water_distance_entry.pack(side="left", padx=LAYOUT['padding_small'])
        self.water_distance_help.pack(side="left", padx=LAYOUT['padding_medium'])
        
        # Protonation options
        self.protonation_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        self.notprotonate_checkbox.pack(anchor="w", pady=LAYOUT['padding_small'])
        self.protonation_help_text.pack(fill="x", pady=LAYOUT['padding_small'])
        
        # Actions section
        self.actions_section.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.actions_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.actions_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        self.validate_button.pack(side="left", padx=LAYOUT['padding_small'])
        self.prepare_button.pack(side="left", padx=LAYOUT['padding_small'])
        self.load_defaults_button.pack(side="left", padx=LAYOUT['padding_small'])
        
        # Progress section
        self.progress_section.pack(fill="both", expand=True, padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_medium'])
        
        self.progress_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.progress_tracker.pack(fill="both", expand=True, padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
    
    def _load_defaults(self):
        """Load default force field settings and POPC lipids."""
        try:
            recommendations = self.ff_manager.get_recommendations("membrane")
            
            self.water_combo.set(recommendations["water_model"])
            self.protein_ff_combo.set(recommendations["protein_ff"])
            self.lipid_ff_combo.set(recommendations["lipid_ff"])
            
            # Set POPC as default lipid for both leaflets
            available_lipids = self.ff_manager.get_available_lipids()
            if "POPC" in available_lipids:
                # Upper leaflet with POPC
                self.upper_leaflet.set_selected_lipids(["POPC"], [1.0])
                # Lower leaflet with POPC  
                self.lower_leaflet.set_selected_lipids(["POPC"], [1.0])
            
            # Force the frame to update and be visible
            self.ff_section.update()
            self.ff_frame.update()
            
            # Debug: print widget info
            logger.info(f"Defaults loaded - Water: {recommendations['water_model']}, Protein: {recommendations['protein_ff']}, Lipid: {recommendations['lipid_ff']}, Default lipid: POPC")
            
            if self.status_callback:
                self.status_callback("Default settings loaded with POPC lipids")
            
        except Exception as e:
            logger.error(f"Error loading defaults: {e}")
    
    def _browse_working_directory(self):
        """Browse for working directory."""
        directory = filedialog.askdirectory(
            title="Select Working Directory",
            initialdir=self.working_directory or self.initial_directory
        )
        
        if directory:
            self.workdir_entry.delete(0, "end")
            self.workdir_entry.insert(0, directory)
            self.working_directory = directory
            set_working_directory(directory)
            
            # Update progress tracker
            self.progress_tracker.set_working_directory(directory)
    
    def _browse_working_file(self):
        """Browse for working PDB file."""
        file_path = filedialog.askopenfilename(
            title="Select PDB File for Preparation",
            filetypes=[("PDB files", "*.pdb"), ("All files", "*.*")],
            initialdir=self.working_directory or self.initial_directory
        )
        
        if file_path:
            self.workfile_entry.delete(0, "end")
            self.workfile_entry.insert(0, file_path)
            self.current_pdb_file = file_path
            
            # Generate automatic output folder name based on PDB file
            pdb_name = Path(file_path).stem
            auto_output_name = f"preparation_{pdb_name}"
            
            # Update output name only if it's still the default
            if self.preparation_output_name == "preparation":
                self.preparation_output_name = auto_output_name
                self.outputname_entry.delete(0, "end")
                self.outputname_entry.insert(0, auto_output_name)
            
            if self.status_callback:
                self.status_callback(f"Selected working file: {Path(file_path).name}")
    
    def get_working_file(self) -> Optional[str]:
        """Get the current working file path."""
        file_path = self.workfile_entry.get().strip()
        return file_path if file_path else self.current_pdb_file
    
    def load_working_file(self, file_path: str):
        """Load a working file from external source (e.g., from propka frame)."""
        if file_path and Path(file_path).exists():
            self.workfile_entry.delete(0, "end")
            self.workfile_entry.insert(0, file_path)
            self.current_pdb_file = file_path
            
            # Generate automatic output folder name based on PDB file
            pdb_name = Path(file_path).stem
            auto_output_name = f"preparation_{pdb_name}"
            
            # Update output name only if it's still the default
            if self.preparation_output_name == "preparation":
                self.preparation_output_name = auto_output_name
                self.outputname_entry.delete(0, "end")
                self.outputname_entry.insert(0, auto_output_name)
            
            if self.status_callback:
                self.status_callback(f"Loaded working file: {Path(file_path).name}")
    
    def _toggle_salt_options(self):
        """Toggle salt options based on checkbox state."""
        enabled = self.add_salt_var.get()
        state = "normal" if enabled else "disabled"
        
        self.salt_conc_entry.configure(state=state)
        self.cation_combo.configure(state=state)
        self.anion_combo.configure(state=state)

    def _collect_inputs(self) -> dict:
        """Collect all input parameters for simplified workflow."""
        # Get working file - prioritize the working file entry over current_pdb_file
        working_file = self.workfile_entry.get().strip()
        if not working_file:
            working_file = self.current_pdb_file or ""
        
        # Get lipid ratios
        upper_lipids_ratios = self.upper_leaflet.get_lipids_with_ratios()
        lower_lipids_ratios = self.lower_leaflet.get_lipids_with_ratios()

        upper_ratio_str = ":".join([str(ratio) for _, ratio in upper_lipids_ratios]) if upper_lipids_ratios else ""
        lower_ratio_str = ":".join([str(ratio) for _, ratio in lower_lipids_ratios]) if lower_lipids_ratios else ""

        if upper_ratio_str and lower_ratio_str:
            lipid_ratios = f"{upper_ratio_str}//{lower_ratio_str}"
        elif upper_ratio_str:
            lipid_ratios = upper_ratio_str
        elif lower_ratio_str:
            lipid_ratios = f"//{lower_ratio_str}"
        else:
            lipid_ratios = ""

        # Determine command flags based on skip protonation option
        skip_protonation = self.notprotonate_var.get()
        
        return {
            'pdb_file': working_file,
            'upper_lipids': self.upper_leaflet.get_selected_lipids(),
            'lower_lipids': self.lower_leaflet.get_selected_lipids(),
            'lipid_ratios': lipid_ratios,
            'water_model': self.water_combo.get(),
            'protein_ff': self.protein_ff_combo.get(),
            'lipid_ff': self.lipid_ff_combo.get(),
            'preoriented': self.preoriented_var.get(),
            'parametrize': self.parametrize_var.get(),  # Use checkbox value for tleap parametrization
            'add_salt': self.add_salt_var.get(),
            'salt_concentration': self.salt_conc_entry.get().strip(),
            'cation': self.cation_combo.get(),
            'anion': self.anion_combo.get(),
            'dist_wat': self.water_distance_entry.get().strip(),
            'notprotonate': skip_protonation,  # Skip protonation if requested
            'simplified_workflow': True  # Flag to indicate simplified workflow
        }

    def _validate_inputs(self):
        """Validate all preparation inputs for simplified workflow."""
        try:
            if self.status_callback:
                self.status_callback("Validating inputs...")

            # Get current inputs
            inputs = self._collect_inputs()

            # Basic validation
            if not inputs['pdb_file']:
                messagebox.showerror("Validation Failed", "Please select a working PDB file")
                self.prepare_button.configure(state="disabled")
                return

            if not os.path.exists(inputs['pdb_file']):
                messagebox.showerror("Validation Failed", f"PDB file not found: {inputs['pdb_file']}")
                self.prepare_button.configure(state="disabled")
                return

            if not inputs['upper_lipids'] and not inputs['lower_lipids']:
                messagebox.showerror("Validation Failed", "Please select lipids for at least one leaflet")
                self.prepare_button.configure(state="disabled")
                return
                    
            # Use validator for detailed validation
            valid, error_msg = self.validator.validate_system_inputs(**inputs)

            if valid:
                self.prepare_button.configure(state="normal")
                messagebox.showinfo("Validation Successful", "All inputs are valid. Ready to start preparation.")
                if self.status_callback:
                    self.status_callback("Validation successful")
            else:
                self.prepare_button.configure(state="disabled")
                messagebox.showerror("Validation Failed", error_msg)
                if self.status_callback:
                    self.status_callback("Validation failed")

        except Exception as e:
            logger.error(f"Error during validation: {e}")
            messagebox.showerror("Validation Error", f"Error during validation: {str(e)}")
            self.prepare_button.configure(state="disabled")

    def _start_preparation(self):
        """Start the system preparation process with simplified workflow."""
        try:
            if self.status_callback:
                self.status_callback("Starting system preparation...")

            # Get working directory
            working_dir = self.workdir_entry.get().strip()
            if not working_dir:
                messagebox.showerror("No Directory", "Please select a working directory.")
                return

            # Collect inputs
            inputs = self._collect_inputs()

            # Configure system builder with simplified settings
            config = {
                'water_model': inputs['water_model'],
                'protein_ff': inputs['protein_ff'],
                'lipid_ff': inputs['lipid_ff'],
                'preoriented': inputs['preoriented'],
                'parametrize': inputs['parametrize'],
                'add_salt': inputs['add_salt'],
                'salt_concentration': float(inputs['salt_concentration']) if inputs['salt_concentration'] else 0.15,
                'cation': inputs['cation'],
                'anion': inputs['anion'],
                'dist_wat': float(inputs['dist_wat']) if inputs['dist_wat'] else 17.5,
                'notprotonate': inputs['notprotonate']
            }

            self.system_builder.set_configuration(**config)

            # Start preparation
            success, message, job_dir = self.system_builder.prepare_system(
                pdb_file=inputs['pdb_file'],
                working_dir=working_dir,
                upper_lipids=inputs['upper_lipids'],
                lower_lipids=inputs['lower_lipids'],
                lipid_ratios=inputs['lipid_ratios'],
                output_folder_name=self.preparation_output_name
            )

            if success:
                messagebox.showinfo("Preparation Started", message)
                if self.status_callback:
                    self.status_callback("System preparation started")

                # Refresh progress tracker
                self.progress_tracker.refresh_jobs()

            else:
                messagebox.showerror("Preparation Failed", message)
                if self.status_callback:
                    self.status_callback("Preparation failed")

        except BuilderError as e:
            logger.error(f"System builder error: {e}")
            messagebox.showerror("Preparation Error", str(e))
            if self.status_callback:
                self.status_callback("Preparation error")

        except Exception as e:
            logger.error(f"Unexpected error during preparation: {e}", exc_info=True)
            messagebox.showerror("Unexpected Error", f"An unexpected error occurred: {str(e)}")
            if self.status_callback:
                self.status_callback("Unexpected error")

    
    def set_current_pdb_file(self, pdb_file: str):
        """Set the current PDB file (for compatibility)."""
        if pdb_file and Path(pdb_file).exists():
            self.load_working_file(pdb_file)
    
    def _load_current_pdb(self):
        """Load current PDB file from callback (for compatibility)."""
        if self.get_current_pdb:
            self.current_pdb_file = self.get_current_pdb()
            if self.current_pdb_file:
                self.load_working_file(self.current_pdb_file)

    def set_pdb_file(self, pdb_file: str):
        """Set PDB file (for compatibility).""" 
        self.current_pdb_file = pdb_file
        if pdb_file:
            self.load_working_file(pdb_file)

    def _validate_stage2_pdb(self, packed_pdb_path: str) -> tuple[bool, str]:
        """
        Validate that the Stage 2 PDB file appears to have Propka modifications.
        
        Args:
            packed_pdb_path: Path to the packed PDB file
            
        Returns:
            Tuple of (is_valid, message)
        """
        try:
            with open(packed_pdb_path, 'r') as f:
                content = f.read()
            
            # Check for common Propka-modified residue names
            propka_residues = ['GLH', 'ASH', 'HIP', 'HIE', 'HID', 'LYN', 'TYM', 'CYM']
            found_modifications = []
            
            for res_name in propka_residues:
                if res_name in content:
                    count = content.count(res_name)
                    found_modifications.append(f"{res_name} ({count} atoms)")
            
            if found_modifications:
                message = f"Found Propka modifications: {', '.join(found_modifications)}"
                return True, message
            else:
                # No obvious modifications found - warn but don't block
                warning_msg = (
                    "No obvious Propka modifications detected in the PDB file. "
                    "Make sure you have applied the necessary residue name changes "
                    "(e.g., GLU→GLH, ASP→ASH, HIS→HIP/HIE/HID) based on your Propka analysis."
                )
                return True, warning_msg
                
        except Exception as e:
            return False, f"Error reading PDB file: {str(e)}"
    
    def _show_propka_workflow_help(self):
        """Show detailed help for the Propka workflow."""
        help_text = """Propka Two-Stage Workflow Guide:

OVERVIEW:
This workflow is designed for proteins that need specific protonation states
based on Propka analysis, following the packmol-memgen tutorial recommendations.

STAGE 1 - PACKING:
1. Select your original PDB file and lipid composition
2. Run Stage 1 to pack the protein into the membrane
3. This creates a packed PDB file (bilayer_*.pdb) without parametrization

BETWEEN STAGES:
1. Run Propka analysis on your ORIGINAL PDB file (not the packed one)
2. Identify residues that need different protonation states at your target pH
3. Modify residue names in the PACKED PDB file based on Propka results
   Example: Change GLU to GLH for protonated glutamate

STAGE 2 - PARAMETRIZATION:
1. Select the modified packed PDB file
2. Run Stage 2 with --notprotonate flag
3. This generates topology files ready for MD simulation

WHY TWO STAGES?
- Packmol-memgen normally re-protonates proteins during packing
- This can override your desired protonation states
- The two-stage approach preserves your Propka-based modifications

EXAMPLE WORKFLOW:
1. Original: protein.pdb with GLU 71
2. Stage 1: Creates bilayer_protein.pdb (GLU 71 renumbered to GLU 49)
3. Propka: Shows GLU 71 should be protonated at pH 7
4. Modify: Change GLU 49 to GLH 49 in bilayer_protein.pdb
5. Stage 2: Parametrize bilayer_protein_modified.pdb with --notprotonate
"""
        
        # Create help window
        help_window = ctk.CTkToplevel(self)
        help_window.title("Propka Workflow Help")
        help_window.geometry("700x600")
        help_window.resizable(True, True)
        help_window.transient(self.winfo_toplevel())
        help_window.grab_set()
        
        # Center the help window
        help_window.update_idletasks()
        x = self.winfo_rootx() + (self.winfo_width() // 2) - (700 // 2)
        y = self.winfo_rooty() + (self.winfo_height() // 2) - (600 // 2)
        help_window.geometry(f"+{x}+{y}")
        
        # Create help content
        help_frame = ctk.CTkFrame(help_window)
        help_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        help_label = ctk.CTkLabel(
            help_frame,
            text="Propka Two-Stage Workflow",
            font=FONTS['heading']
        )
        help_label.pack(pady=(0, 10))
        
        help_textbox = ctk.CTkTextbox(
            help_frame,
            font=FONTS['small'],
            wrap="word"
        )
        help_textbox.pack(fill="both", expand=True, pady=(0, 10))
        help_textbox.insert("1.0", help_text)
        help_textbox.configure(state="disabled")
        
        close_button = ctk.CTkButton(
            help_frame,
            text="Close",
            command=help_window.destroy
        )
        close_button.pack()
    
    def _open_propka_workflow_helper(self):
        """Open the Propka workflow helper dialog."""
        try:
            from gatewizard.gui.widgets.propka_workflow_helper import PropkaWorkflowDialog
            
            # Try to find a packed PDB file in the working directory
            packed_pdb_path = ""
            working_dir = Path(self.working_directory)
            if working_dir.exists():
                packed_files = list(working_dir.glob("**/bilayer_*.pdb"))
                if packed_files:
                    packed_pdb_path = str(packed_files[-1])  # Use most recent
            
            dialog = PropkaWorkflowDialog(self, packed_pdb_path)
            dialog.lift()
            dialog.focus_force()
            
        except ImportError:
            messagebox.showinfo(
                "Helper Not Available",
                "The Propka workflow helper is not available. "
                "Please follow the manual workflow steps in the help guide."
            )
        except Exception as e:
            logger.error(f"Error opening Propka workflow helper: {e}")
            messagebox.showerror(
                "Error",
                f"Could not open Propka workflow helper: {str(e)}"
            )
    
    def _setup_propka_bindings(self):
        """Setup bindings for Propka workflow components."""
        if hasattr(self, 'propka_stage1_radio'):
            self.propka_stage1_radio.configure(command=self._update_propka_help_text)  # type: ignore
        if hasattr(self, 'propka_stage2_radio'):
            self.propka_stage2_radio.configure(command=self._update_propka_help_text)  # type: ignore

    def on_stage_shown(self):
        """Called when this stage becomes active."""
        # Update current PDB file
        if self.get_current_pdb:
            self.current_pdb_file = self.get_current_pdb()
        
        # Update working directory in progress tracker
        self.progress_tracker.set_working_directory(self.working_directory)
        self.progress_tracker.refresh_jobs()
    
    def on_pdb_changed(self, pdb_file: Optional[str]):
        """Called when PDB file changes."""
        self.current_pdb_file = pdb_file
        
        # Reset validation state when PDB changes
        self.prepare_button.configure(state="disabled")
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the frame."""
        try:
            def safe_update_font(widget, font, width=None, height=None):
                try:
                    if widget and hasattr(widget, 'configure'):
                        kwargs = {'font': font}
                        if width is not None:
                            kwargs['width'] = width
                        if height is not None:
                            kwargs['height'] = height
                        widget.configure(**kwargs)
                except Exception:
                    pass

            # Main labels
            safe_update_font(self.workdir_label, scaled_fonts['heading'])
            safe_update_font(self.workfile_label, scaled_fonts['heading'])
            safe_update_font(self.lipids_label, scaled_fonts['heading'])
            safe_update_font(self.ff_label, scaled_fonts['heading'])
            if hasattr(self, 'options_label'):
                safe_update_font(self.options_label, scaled_fonts['heading'])
            if hasattr(self, 'actions_label'):
                safe_update_font(self.actions_label, scaled_fonts['heading'])
            if hasattr(self, 'progress_label'):
                safe_update_font(self.progress_label, scaled_fonts['heading'])
            if hasattr(self, 'job_progress_label'):
                safe_update_font(self.job_progress_label, scaled_fonts['heading'])  # type: ignore

            # Body text labels
            safe_update_font(self.water_label, scaled_fonts['body'])
            safe_update_font(self.protein_ff_label, scaled_fonts['body'])
            safe_update_font(self.lipid_ff_label, scaled_fonts['body'])

            # Searchable comboboxes fonts
            for combo_attr in ['water_combo', 'protein_ff_combo', 'lipid_ff_combo', 'cation_combo', 'anion_combo']:
                if hasattr(self, combo_attr):
                    combo = getattr(self, combo_attr)
                    if hasattr(combo, 'update_fonts'):
                        try:
                            combo.update_fonts(scaled_fonts)
                        except Exception:
                            pass

            # Buttons (make all consistent)
            for btn in [
                'workdir_browse_button', 'workfile_browse_button', 'prepare_button',
                'validate_button', 'defaults_button', 'refresh_button', 'clear_completed_button',
                'load_defaults_button']:
                if hasattr(self, btn):
                    safe_update_font(getattr(self, btn), scaled_fonts['body'], width=120, height=32)

            # Options section
            if hasattr(self, 'add_salt_checkbox'):
                safe_update_font(self.add_salt_checkbox, scaled_fonts['body'])
            if hasattr(self, 'preoriented_checkbox'):
                safe_update_font(self.preoriented_checkbox, scaled_fonts['body'])
            if hasattr(self, 'parametrize_checkbox'):
                safe_update_font(self.parametrize_checkbox, scaled_fonts['body'])
            if hasattr(self, 'notprotonate_checkbox'):
                safe_update_font(self.notprotonate_checkbox, scaled_fonts['body'])
            if hasattr(self, 'box_info_label'):
                safe_update_font(self.box_info_label, scaled_fonts['small'])  # type: ignore
            if hasattr(self, 'sol_label'):
                safe_update_font(self.sol_label, scaled_fonts['body'])  # type: ignore
            if hasattr(self, 'cation_label'):
                safe_update_font(self.cation_label, scaled_fonts['body'])
            if hasattr(self, 'anion_label'):
                safe_update_font(self.anion_label, scaled_fonts['body'])
            if hasattr(self, 'salt_conc_label'):
                safe_update_font(self.salt_conc_label, scaled_fonts['body'])
            if hasattr(self, 'skip_prep_checkbox'):
                safe_update_font(self.skip_prep_checkbox, scaled_fonts['body'])  # type: ignore
            if hasattr(self, 'skip_build_checkbox'):
                safe_update_font(self.skip_build_checkbox, scaled_fonts['body'])  # type: ignore

            # Protonation help text
            if hasattr(self, 'protonation_help_text'):
                try:
                    self.protonation_help_text.configure(font=scaled_fonts['small'])
                except Exception:
                    pass

            # Update leaflet frames
            if hasattr(self, 'upper_leaflet') and hasattr(self.upper_leaflet, 'update_fonts'):
                self.upper_leaflet.update_fonts(scaled_fonts)
            if hasattr(self, 'lower_leaflet') and hasattr(self.lower_leaflet, 'update_fonts'):
                self.lower_leaflet.update_fonts(scaled_fonts)
            
        except Exception as e:
            logger.warning(f"Error updating fonts in BuilderFrame: {e}")
    
    def cleanup(self):
        """Cleanup resources when frame is destroyed."""
        if hasattr(self.progress_tracker, 'cleanup'):
            self.progress_tracker.cleanup()