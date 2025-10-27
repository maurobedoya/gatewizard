# gatewizard/gui/frames/visualize.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Visualization frame for protein structure viewing.

This module provides the GUI for loading and visualizing protein structures
using matplotlib-based 3D visualization.
"""

import tkinter as tk
from tkinter import filedialog, messagebox
from typing import Optional, Callable
import warnings
from pathlib import Path
import numpy as np

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

# Import matplotlib components
try:
    import matplotlib
    matplotlib.use('TkAgg')  # Set backend before importing pyplot
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import Axes3D
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    # Set to None so type checkers know these are unavailable
    FigureCanvasTkAgg = None  # type: ignore
    Figure = None  # type: ignore

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, REPRESENTATIONS, COLOR_SCHEMES, FILE_FILTERS,
    WIDGET_SIZES, LAYOUT
)
from gatewizard.core.file_manager import FileManager
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class VisualizationError(Exception):
    """Custom exception for visualization errors."""
    pass

class VisualizeFrame(ctk.CTkFrame):
    """
    Frame for protein structure visualization.
    
    This frame handles PDB file loading, molecular visualization controls,
    and matplotlib-based 3D structure viewing with an improved layout.
    """
    
    def __init__(
        self,
        parent,
        pdb_changed_callback: Optional[Callable[[Optional[str]], None]] = None,
        status_callback: Optional[Callable[[str], None]] = None,
        initial_directory: Optional[str] = None
    ):
        """
        Initialize the visualization frame.
        
        Args:
            parent: Parent widget
            pdb_changed_callback: Callback for PDB file changes
            status_callback: Callback for status updates
            initial_directory: Initial directory for file dialogs
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['content_bg'])
        
        self.pdb_changed_callback = pdb_changed_callback
        self.status_callback = status_callback
        self.initial_directory = initial_directory or str(Path.cwd())
        
        # State variables
        self.current_pdb_file = None
        self.molecular_viewer = None
        self.file_manager = FileManager()
        
        # Matplotlib components
        self.figure = None
        self.ax = None
        self.canvas = None
        
        # Element colors for visualization
        self.element_colors = {
            'C': '#808080', 'N': '#0000FF', 'O': '#FF0000',
            'S': '#FFFF00', 'P': '#FFA500', 'H': '#FFFFFF',
            'default': '#00FF00'
        }
        
        # Secondary structure colors
        self.ss_colors = {
            'H': '#FF0000', 'E': '#0000FF', 'C': '#00FF00',
            'T': '#FFFF00', 'S': '#FF00FF', 'default': '#FFFFFF'
        }
        
        # Chain colors
        self.chain_colors = ['#FF0000', '#00FF00', '#0000FF', '#FFFF00', 
                           '#FF00FF', '#00FFFF', '#FFA500', '#800080']
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        
        # Initialize matplotlib viewer
        self._initialize_viewer()
        
        # Set initial global status message
        self._safe_status_callback("Ready to visualize molecular structures")
    
    def _create_widgets(self):
        """Create all widgets for the visualization frame."""
        # Main horizontal container
        self.main_container = ctk.CTkFrame(self, fg_color="transparent")
        
        # LEFT PANEL - Controls and Information (fixed width)
        self.left_panel = ctk.CTkFrame(self.main_container, fg_color=COLOR_SCHEME['content_inside_bg'], width=350)
        self.left_panel.pack_propagate(False)  # Maintain fixed width
        
        # File selection section
        self.file_section = ctk.CTkFrame(self.left_panel, fg_color=COLOR_SCHEME['canvas'])
        
        self.file_label = ctk.CTkLabel(
            self.file_section,
            text="PDB File",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.file_frame = ctk.CTkFrame(self.file_section, fg_color="transparent")
        
        self.file_entry = ctk.CTkEntry(
            self.file_frame,
            placeholder_text="Select PDB file...",
            width=250,
            height=WIDGET_SIZES['entry_height']
        )
        
        self.browse_button = ctk.CTkButton(
            self.file_frame,
            text="Browse",
            width=80,
            height=WIDGET_SIZES['button_height'],
            command=self.open_file_dialog
        )
        
        self.load_button = ctk.CTkButton(
            self.file_frame,
            text="Load",
            width=80,
            height=WIDGET_SIZES['button_height'],
            command=self._load_structure,
            state="disabled"
        )
        
        # Structure information section (compact)
        self.info_section = ctk.CTkFrame(self.left_panel, fg_color=COLOR_SCHEME['canvas'])
        
        self.info_label = ctk.CTkLabel(
            self.info_section,
            text="Structure Info",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.info_text = ctk.CTkTextbox(
            self.info_section,
            height=80,  # Reduced height
            width=320,
            font=FONTS['small']
        )
        self.info_text.insert("1.0", "No structure loaded")
        self.info_text.configure(state="disabled")
        
        # Visualization controls section (compact)
        self.controls_section = ctk.CTkFrame(self.left_panel, fg_color=COLOR_SCHEME['canvas'])
        
        self.controls_label = ctk.CTkLabel(
            self.controls_section,
            text="Controls",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Representation controls
        self.repr_frame = ctk.CTkFrame(self.controls_section, fg_color="transparent")
        
        self.repr_label = ctk.CTkLabel(
            self.repr_frame,
            text="Style:",
            font=FONTS['small'],
            width=60
        )
        
        self.repr_combo = ctk.CTkComboBox(
            self.repr_frame,
            values=["Atoms", "Backbone", "Secondary Structure", "Cartoon", "Both"],
            width=180,
            height=25,
            command=self._change_representation
        )
        self.repr_combo.set("Secondary Structure")
        
        # Color scheme controls
        self.color_frame = ctk.CTkFrame(self.controls_section, fg_color="transparent")
        
        self.color_label = ctk.CTkLabel(
            self.color_frame,
            text="Color:",
            font=FONTS['small'],
            width=60
        )
        
        self.color_combo = ctk.CTkComboBox(
            self.color_frame,
            values=["Element", "Chain", "Secondary Structure", "Uniform"],
            width=180,
            height=25,
            command=self._change_color_scheme
        )
        self.color_combo.set("Secondary Structure")
        
        # Atom size controls
        self.size_frame = ctk.CTkFrame(self.controls_section, fg_color="transparent")
        
        self.size_label = ctk.CTkLabel(
            self.size_frame,
            text="Size:",
            font=FONTS['small'],
            width=60
        )
        
        self.atom_size_slider = ctk.CTkSlider(
            self.size_frame,
            from_=5,
            to=100,
            number_of_steps=19,
            width=180,
            command=self._update_atom_size
        )
        self.atom_size_slider.set(30)
        
        # View controls (compact)
        self.view_frame = ctk.CTkFrame(self.controls_section, fg_color="transparent")
        
        self.reset_view_button = ctk.CTkButton(
            self.view_frame,
            text="Reset",
            width=60,
            height=25,
            command=self._reset_view
        )
        
        self.center_button = ctk.CTkButton(
            self.view_frame,
            text="Center",
            width=60,
            height=25,
            command=self._center_view
        )
        
        # Show axes checkbox
        self.axes_var = ctk.BooleanVar(value=False)
        self.show_axes_checkbox = ctk.CTkCheckBox(
            self.view_frame,
            text="Axes",
            variable=self.axes_var,
            command=self._toggle_axes,
            width=60
        )
        
        # Export controls (compact)
        self.export_frame = ctk.CTkFrame(self.controls_section, fg_color="transparent")
        
        self.export_image_button = ctk.CTkButton(
            self.export_frame,
            text="Export Image",
            width=100,
            height=25,
            command=self._export_image
        )
        
        self.save_session_button = ctk.CTkButton(
            self.export_frame,
            text="Save Session",
            width=100,
            height=25,
            command=self._save_session
        )
        
        # RIGHT PANEL - 3D Visualization (expandable)
        self.right_panel = ctk.CTkFrame(self.main_container, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.viz_label = ctk.CTkLabel(
            self.right_panel,
            text="3D Structure View",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Container for matplotlib canvas
        self.viz_container = ctk.CTkFrame(self.right_panel, fg_color=COLOR_SCHEME['viewer_bg'])
        
        # Initially disable controls
        self._set_controls_enabled(False)
    
    def _setup_layout(self):
        """Setup the layout with left panel for controls and right panel for 3D view."""
        # Main container fills full frame space - no bottom padding to avoid pushing status bar off-screen
        self.main_container.pack(fill="both", expand=True, padx=5, pady=(5, 2))
        
        # LEFT PANEL - Fixed width for controls
        self.left_panel.pack(side="left", fill="y", padx=(0, 5), pady=0)
        
        # File selection section
        self.file_section.pack(fill="x", padx=10, pady=(10, 5))
        
        self.file_label.pack(anchor="w", padx=10, pady=(10, 5))
        
        self.file_frame.pack(fill="x", padx=10, pady=(0, 10))
        
        self.file_entry.pack(side="top", fill="x", pady=(0, 5))
        
        button_frame = ctk.CTkFrame(self.file_frame, fg_color="transparent")
        button_frame.pack(fill="x")
        self.browse_button.pack(side="left", padx=(0, 5))
        self.load_button.pack(side="left")
        
        # Structure info section
        self.info_section.pack(fill="x", padx=10, pady=5)
        
        self.info_label.pack(anchor="w", padx=10, pady=(10, 5))
        self.info_text.pack(fill="x", padx=10, pady=(0, 10))
        
        # Controls section
        self.controls_section.pack(fill="both", expand=True, padx=10, pady=5)
        
        self.controls_label.pack(anchor="w", padx=10, pady=(10, 5))
        
        # Representation controls
        self.repr_frame.pack(fill="x", padx=10, pady=2)
        self.repr_label.pack(side="left", padx=(0, 5))
        self.repr_combo.pack(side="left", fill="x", expand=True)
        
        # Color controls
        self.color_frame.pack(fill="x", padx=10, pady=2)
        self.color_label.pack(side="left", padx=(0, 5))
        self.color_combo.pack(side="left", fill="x", expand=True)
        
        # Size controls
        self.size_frame.pack(fill="x", padx=10, pady=2)
        self.size_label.pack(side="left", padx=(0, 5))
        self.atom_size_slider.pack(side="left", fill="x", expand=True)
        
        # View controls
        self.view_frame.pack(fill="x", padx=10, pady=5)
        self.reset_view_button.pack(side="left", padx=(0, 5))
        self.center_button.pack(side="left", padx=(0, 5))
        self.show_axes_checkbox.pack(side="left", padx=(10, 0))
        
        # Export controls - reduced bottom padding to match other tabs
        self.export_frame.pack(fill="x", padx=10, pady=(5, 5))
        self.export_image_button.pack(side="top", fill="x", pady=(0, 5))
        self.save_session_button.pack(side="top", fill="x")
        
        # RIGHT PANEL - Expandable 3D view
        self.right_panel.pack(side="right", fill="both", expand=True, padx=0, pady=0)
        
        self.viz_label.pack(anchor="w", padx=15, pady=(15, 5))
        
        # Viz container takes remaining space - no bottom padding to match other tabs
        self.viz_container.pack(fill="both", expand=True, padx=15, pady=(0, 0))

    def _initialize_viewer(self):
        """Initialize the matplotlib-based molecular viewer."""
        if not MATPLOTLIB_AVAILABLE:
            self._safe_status_callback("Matplotlib not available")
            logger.error("Matplotlib not available for visualization")
            return

        try:
            from gatewizard.tools.molecular_viewer import MolecularViewer

            self.molecular_viewer = MolecularViewer()

            # Create matplotlib figure with dark theme and constrained layout to reduce clipping
            if Figure is not None:
                self.figure = Figure(figsize=(10, 8), facecolor='#212121', dpi=100, constrained_layout=True)
                self.ax = self.figure.add_subplot(111, projection='3d')
                self.ax.set_facecolor('#212121')
            else:
                raise ImportError("Figure class not available")

            # Hide axes by default
            self._hide_axes()

            # Create canvas in the visualization container
            if FigureCanvasTkAgg is not None and self.figure is not None:
                self.canvas = FigureCanvasTkAgg(self.figure, self.viz_container)
                self.canvas.draw()
                # Remove extra padding so the canvas uses all available space
                self.canvas.get_tk_widget().pack(fill="both", expand=True)
            else:
                raise ImportError("FigureCanvasTkAgg class not available")
            # Reduce margins to avoid content getting cut when zoomed
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="The figure layout has changed to tight")
                    self.figure.tight_layout()
            except Exception:
                pass

            # Bind scroll event for zooming
            self.canvas.get_tk_widget().bind("<MouseWheel>", self._on_scroll)
            self.canvas.get_tk_widget().bind("<Button-4>", self._on_scroll)
            self.canvas.get_tk_widget().bind("<Button-5>", self._on_scroll)

            logger.info("Matplotlib molecular viewer initialized successfully")
            self._safe_status_callback("3D viewer initialized")

        except Exception as e:
            logger.error(f"Failed to initialize molecular viewer: {e}")
            self._safe_status_callback("3D viewer initialization failed")

    def _safe_status_callback(self, message: str):
        """Safely call status callback with error handling and update local status."""
        try:
            # Update the global status bar (main app)
            if self.status_callback:
                self.status_callback(message)
                
        except Exception as e:
            logger.warning(f"Status callback failed: {e}")
            # Continue without crashing
    
    def open_file_dialog(self):
        """Open file dialog to select PDB file."""
        file_path = filedialog.askopenfilename(
            title="Select PDB File",
            filetypes=FILE_FILTERS['pdb'] + FILE_FILTERS['all'],
            initialdir=self.initial_directory
        )
        
        if file_path:
            self.file_entry.delete(0, "end")
            self.file_entry.insert(0, file_path)
            self.load_button.configure(state="normal")
    
    def _load_structure(self):
        """Load the selected PDB structure."""
        file_path = self.file_entry.get().strip()
        
        if not file_path:
            messagebox.showwarning("No File", "Please select a PDB file first.")
            return
        
        # Validate file
        valid, error_msg = self.file_manager.validate_pdb_file(file_path)
        if not valid:
            messagebox.showerror("Invalid File", error_msg)
            return
        
        try:
            if self.status_callback:
                self.status_callback("Loading structure...")
            
            # Load structure in molecular viewer
            if self.molecular_viewer:
                success = self.molecular_viewer.load_structure(file_path, "protein")
                
                if success:
                    self.current_pdb_file = file_path
                    
                    # Update the 3D visualization
                    self._update_3d_visualization()
                    
                    # Extract and display file information
                    self._update_file_info(file_path)
                    
                    # Enable controls
                    self._set_controls_enabled(True)
                    
                    # Notify callback
                    if self.pdb_changed_callback:
                        self.pdb_changed_callback(file_path)
                    
                    if self.status_callback:
                        self.status_callback(f"Loaded: {Path(file_path).name}")
                    
                    logger.info(f"Structure loaded successfully: {file_path}")
                    
                else:
                    messagebox.showerror("Load Error", "Failed to load structure")
            else:
                messagebox.showerror("Viewer Error", "3D viewer is not available")
                
        except Exception as e:
            logger.error(f"Error loading structure: {e}")
            messagebox.showerror("Error", f"Failed to load structure: {str(e)}")
            if self.status_callback:
                self.status_callback("Load failed")
    
    def _update_file_info(self, file_path: str):
        """Update the file information display with compact format."""
        try:
            # Extract file information
            info = self.file_manager.extract_pdb_info(file_path)
            
            # Format information text (more compact)
            info_text = f"""File: {info['filename']}
PDB ID: {info['pdb_id'] or 'Unknown'}
Atoms: {info['num_atoms']:,} | Residues: {info['num_residues']:,}
Chains: {', '.join(info['chains']) if info['chains'] else 'None'}"""

            if info['resolution']:
                info_text += f"\nResolution: {info['resolution']:.2f} Å"
            
            # Update text widget
            self.info_text.configure(state="normal")
            self.info_text.delete("1.0", "end")
            self.info_text.insert("1.0", info_text)
            self.info_text.configure(state="disabled")
            
        except Exception as e:
            logger.error(f"Error extracting file info: {e}")
            self.info_text.configure(state="normal")
            self.info_text.delete("1.0", "end")
            self.info_text.insert("1.0", f"Error reading file information: {str(e)}")
            self.info_text.configure(state="disabled")
    
    def _change_representation(self, representation: str):
        """Change the molecular representation."""
        if self.molecular_viewer and self.current_pdb_file:
            self._update_3d_visualization()
            if self.status_callback:
                self.status_callback(f"Representation changed to {representation}")
    
    def _change_color_scheme(self, color_scheme: str):
        """Change the color scheme."""
        if self.molecular_viewer and self.current_pdb_file:
            self._update_3d_visualization()
            if self.status_callback:
                self.status_callback(f"Color scheme changed to {color_scheme}")
    
    def _update_atom_size(self, value=None):
        """Handle atom size slider change."""
        if self.molecular_viewer and self.current_pdb_file:
            self._update_3d_visualization()
    
    def _reset_view(self):
        """Reset the view to default."""
        if self.ax and self.canvas:
            self._set_equal_aspect_3d()
            self.canvas.draw()
            if self.status_callback:
                self.status_callback("View reset")
    
    def _center_view(self):
        """Center the view on the structure."""
        if self.ax and self.canvas and self.molecular_viewer and hasattr(self.molecular_viewer, 'atom_coords'):
            if self.molecular_viewer.atom_coords is not None:
                self._set_equal_aspect_3d()
                self.canvas.draw()
                if self.status_callback:
                    self.status_callback("View centered")
    
    def _toggle_axes(self):
        """Toggle axes visibility."""
        if self.ax and self.canvas:
            if self.axes_var.get():
                self._show_axes()
            else:
                self._hide_axes()
            self.canvas.draw()
    
    def _hide_axes(self):
        """Hide all axes, ticks, labels, and background elements"""
        if not self.ax:
            return
            
        try:
            # Hide ticks and labels
            self.ax.set_xticks([])
            self.ax.set_yticks([])
            if hasattr(self.ax, 'set_zticks') and callable(getattr(self.ax, 'set_zticks', None)):
                self.ax.set_zticks([])  # type: ignore
            self.ax.set_xlabel('')
            self.ax.set_ylabel('')
            if hasattr(self.ax, 'set_zlabel'):
                self.ax.set_zlabel('')
            
            # Hide axis lines and panes (3D specific)
            if hasattr(self.ax, 'xaxis') and hasattr(self.ax.xaxis, 'pane'):
                self.ax.xaxis.pane.fill = False  # type: ignore
                self.ax.yaxis.pane.fill = False  # type: ignore
                self.ax.zaxis.pane.fill = False  # type: ignore
                self.ax.xaxis.pane.set_edgecolor('none')  # type: ignore
                self.ax.yaxis.pane.set_edgecolor('none')  # type: ignore
                self.ax.zaxis.pane.set_edgecolor('none')  # type: ignore
            
            # Hide grid
            self.ax.grid(False)
            
            # Make axis lines invisible (3D specific)
            if hasattr(self.ax, 'xaxis') and hasattr(self.ax.xaxis, 'line'):
                self.ax.xaxis.line.set_color('none')  # type: ignore
                self.ax.yaxis.line.set_color('none')  # type: ignore
                self.ax.zaxis.line.set_color('none')  # type: ignore
            
            # Hide axis spines (3D specific)
            if hasattr(self.ax, 'xaxis') and hasattr(self.ax.xaxis, 'set_pane_color'):
                self.ax.xaxis.set_pane_color((0, 0, 0, 0))  # type: ignore
                self.ax.yaxis.set_pane_color((0, 0, 0, 0))  # type: ignore
                self.ax.zaxis.set_pane_color((0, 0, 0, 0))  # type: ignore
                
        except Exception as e:
            logger.warning(f"Error hiding axes: {e}")
    
    def _show_axes(self):
        """Show axes with proper styling"""
        if not self.ax:
            return
            
        try:
            self.ax.set_xlabel('X', color='white')
            self.ax.set_ylabel('Y', color='white')
            if hasattr(self.ax, 'set_zlabel'):
                self.ax.set_zlabel('Z', color='white')
            self.ax.tick_params(colors='white')
            self.ax.grid(True, alpha=0.3)
            
            # Restore pane colors (3D specific)
            if hasattr(self.ax, 'xaxis') and hasattr(self.ax.xaxis, 'pane'):
                self.ax.xaxis.pane.fill = True  # type: ignore
                self.ax.yaxis.pane.fill = True  # type: ignore
                self.ax.zaxis.pane.fill = True  # type: ignore
                self.ax.xaxis.pane.set_edgecolor('white')  # type: ignore
                self.ax.yaxis.pane.set_edgecolor('white')  # type: ignore
                self.ax.zaxis.pane.set_edgecolor('white')  # type: ignore
                
        except Exception as e:
            logger.warning(f"Error showing axes: {e}")
    
    def _on_scroll(self, event):
        """Handle mouse scroll for zooming"""
        if not self.ax:
            return
            
        if event.delta:
            # Windows
            scale_factor = 1.1 if event.delta > 0 else 0.9
        else:
            # Linux
            scale_factor = 1.1 if event.num == 4 else 0.9
        
        # Get current axis limits
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        zlim = self.ax.get_zlim()
        
        # Calculate center points
        x_center = (xlim[0] + xlim[1]) / 2
        y_center = (ylim[0] + ylim[1]) / 2
        z_center = (zlim[0] + zlim[1]) / 2
        
        # Calculate new limits with padding
        x_range = (xlim[1] - xlim[0]) * scale_factor / 2
        y_range = (ylim[1] - ylim[0]) * scale_factor / 2
        z_range = (zlim[1] - zlim[0]) * scale_factor / 2
        pad_x = x_range * 0.05
        pad_y = y_range * 0.05
        pad_z = z_range * 0.05
        
        # Set new limits
        self.ax.set_xlim(x_center - x_range - pad_x, x_center + x_range + pad_x)
        self.ax.set_ylim(y_center - y_range - pad_y, y_center + y_range + pad_y)
        self.ax.set_zlim(z_center - z_range - pad_z, z_center + z_range + pad_z)
        
        # Update layout to avoid clipping
        try:
            if self.figure:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="The figure layout has changed to tight")
                    self.figure.tight_layout()
        except Exception:
            pass
        
        if self.canvas:
            self.canvas.draw()
    
    def _update_3d_visualization(self):
        """Update the 3D matplotlib visualization."""
        if not self.molecular_viewer or not hasattr(self.molecular_viewer, 'atom_coords'):
            return
        
        if self.molecular_viewer.atom_coords is None or not self.ax:
            return
        
        self.ax.clear()
        self.ax.set_facecolor('#212121')
        
        viz_mode = self.repr_combo.get()
        color_scheme = self.color_combo.get()
        atom_size = self.atom_size_slider.get()
        
        # Plot based on visualization mode
        if viz_mode == "Atoms":
            self._plot_atoms(color_scheme, atom_size)
        elif viz_mode == "Backbone":
            self._plot_backbone(color_scheme)
        elif viz_mode == "Secondary Structure":
            self._plot_secondary_structure(color_scheme)
        elif viz_mode == "Cartoon":
            self._plot_cartoon(color_scheme)
        elif viz_mode == "Both":
            self._plot_atoms(color_scheme, atom_size * 0.5)
            self._plot_backbone(color_scheme)
        
        # Handle axes visibility
        if self.axes_var.get():
            self._show_axes()
        else:
            self._hide_axes()
        
        # Set equal aspect ratio and update
        self._set_equal_aspect_3d()
        try:
            if self.figure:
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", message="The figure layout has changed to tight")
                    self.figure.tight_layout()
        except Exception:
            pass
        if self.canvas:
            self.canvas.draw()
    
    def _plot_atoms(self, color_scheme, size):
        """Plot atoms as spheres"""
        if not self.molecular_viewer or not self.ax:
            return
        if not hasattr(self.molecular_viewer, 'atom_coords') or self.molecular_viewer.atom_coords is None:
            return
            
        atom_coords = self.molecular_viewer.atom_coords
        if len(atom_coords) > 0 and hasattr(self.molecular_viewer, 'atom_elements'):
            colors = self._get_colors(self.molecular_viewer.atom_elements, color_scheme)
            self.ax.scatter(
                atom_coords[:, 0],
                atom_coords[:, 1], 
                atom_coords[:, 2],  # type: ignore
                c=colors,
                s=size,
                alpha=0.8,
                edgecolors='black',
                linewidth=0.5
            )
    
    def _plot_backbone(self, color_scheme):
        """Plot backbone as connected lines"""
        if not self.molecular_viewer or not self.ax:
            return
        if not hasattr(self.molecular_viewer, 'chains'):
            return
            
        for chain_id, chain_data in self.molecular_viewer.chains.items():
            if 'ca_coords' not in chain_data:
                continue
            ca_coords = np.array(chain_data['ca_coords'])
            if len(ca_coords) > 1:
                color = self._get_chain_color(chain_id) if color_scheme == "Chain" else '#FFFFFF'
                self.ax.plot(
                    ca_coords[:, 0],
                    ca_coords[:, 1],
                    ca_coords[:, 2],
                    color=color,
                    linewidth=3,
                    alpha=0.9
                )
    
    def _plot_secondary_structure(self, color_scheme):
        """Plot secondary structure with different representations"""
        if not self.molecular_viewer or not self.ax:
            return
        if not hasattr(self.molecular_viewer, 'chains') or not hasattr(self.molecular_viewer, 'secondary_structure'):
            return
            
        for chain_id, chain_data in self.molecular_viewer.chains.items():
            if 'ca_coords' not in chain_data:
                continue
            ca_coords = np.array(chain_data['ca_coords'])
            if len(ca_coords) < 2:
                continue
            
            # Get secondary structure sequence with fallback
            if hasattr(self.molecular_viewer.secondary_structure, 'get'):
                ss_seq = self.molecular_viewer.secondary_structure.get(chain_id, ['C'] * len(ca_coords))  # type: ignore
            else:
                ss_seq = ['C'] * len(ca_coords)
            
            # Group consecutive residues with same secondary structure
            i = 0
            while i < len(ca_coords) - 1:
                current_ss = ss_seq[i] if i < len(ss_seq) else 'C'
                start_idx = i
                
                # Find end of current secondary structure segment
                while (i < len(ca_coords) - 1 and 
                       i < len(ss_seq) and 
                       ss_seq[i] == current_ss):
                    i += 1
                
                end_idx = i
                segment_coords = ca_coords[start_idx:end_idx + 1]
                
                if len(segment_coords) > 1:
                    color = self.ss_colors.get(current_ss, self.ss_colors['default'])
                    
                    if current_ss == 'H':  # Helix - thicker line
                        self.ax.plot(
                            segment_coords[:, 0],
                            segment_coords[:, 1],
                            segment_coords[:, 2],
                            color=color,
                            linewidth=6,
                            alpha=0.9,
                            solid_capstyle='round'
                        )
                    elif current_ss == 'E':  # Sheet - ribbon-like
                        self.ax.plot(
                            segment_coords[:, 0],
                            segment_coords[:, 1],
                            segment_coords[:, 2],
                            color=color,
                            linewidth=8,
                            alpha=0.7,
                            solid_capstyle='butt'
                        )
                    else:  # Coil - thin line
                        self.ax.plot(
                            segment_coords[:, 0],
                            segment_coords[:, 1],
                            segment_coords[:, 2],
                            color=color,
                            linewidth=2,
                            alpha=0.8,
                            linestyle='-'
                        )
                
                i = max(i, start_idx + 1)  # Ensure progress
    
    def _plot_cartoon(self, color_scheme):
        """Plot cartoon representation"""
        if not self.molecular_viewer or not self.ax:
            return
        if not hasattr(self.molecular_viewer, 'chains'):
            return
            
        for chain_id, chain_data in self.molecular_viewer.chains.items():
            if 'ca_coords' not in chain_data:
                continue
            ca_coords = np.array(chain_data['ca_coords'])
            if len(ca_coords) < 4:
                return
            
            # Smooth the backbone using simple interpolation
            smoothed_coords = self._smooth_backbone(ca_coords)
            
            color = self._get_chain_color(chain_id) if color_scheme == "Chain" else '#00AAFF'
            
            self.ax.plot(
                smoothed_coords[:, 0],
                smoothed_coords[:, 1],
                smoothed_coords[:, 2],
                color=color,
                linewidth=5,
                alpha=0.8,
                solid_capstyle='round'
            )
    
    def _smooth_backbone(self, coords, window_size=3):
        """Apply simple smoothing to backbone coordinates"""
        if len(coords) < window_size:
            return coords
        
        smoothed = np.zeros_like(coords)
        half_window = window_size // 2
        
        for i in range(len(coords)):
            start = max(0, i - half_window)
            end = min(len(coords), i + half_window + 1)
            smoothed[i] = np.mean(coords[start:end], axis=0)
        
        return smoothed
    
    def _get_colors(self, elements, scheme):
        """Get colors based on the selected color scheme"""
        if scheme == "Element":
            return [self.element_colors.get(elem, self.element_colors['default']) for elem in elements]
        elif scheme == "Chain":
            # For atoms, assign colors based on which chain they belong to
            colors = []
            if self.molecular_viewer and hasattr(self.molecular_viewer, 'residues'):
                for residue in self.molecular_viewer.residues:
                    chain_color = self._get_chain_color(residue['chain'])
                    for atom in residue['atoms']:
                        colors.append(chain_color)
            return colors[:len(elements)] if colors else ['#00AAFF'] * len(elements)
        elif scheme == "Secondary Structure":
            colors = []
            if (self.molecular_viewer and 
                hasattr(self.molecular_viewer, 'residues') and 
                hasattr(self.molecular_viewer, 'secondary_structure')):
                for residue in self.molecular_viewer.residues:
                    chain_id = residue['chain']
                    res_idx = len([r for r in self.molecular_viewer.residues[:self.molecular_viewer.residues.index(residue)] if r['chain'] == chain_id])
                    # Get secondary structure sequence with fallback
                    if hasattr(self.molecular_viewer.secondary_structure, 'get'):
                        ss_seq = self.molecular_viewer.secondary_structure.get(chain_id, [])  # type: ignore
                    else:
                        ss_seq = []
                    ss = ss_seq[res_idx] if res_idx < len(ss_seq) else 'C'
                    ss_color = self.ss_colors.get(ss, self.ss_colors['default'])
                    for atom in residue['atoms']:
                        colors.append(ss_color)
            return colors[:len(elements)] if colors else ['#00AAFF'] * len(elements)
        else:  # Uniform
            return ['#00AAFF'] * len(elements)
    
    def _get_chain_color(self, chain_id):
        """Get color for a specific chain"""
        if self.molecular_viewer and hasattr(self.molecular_viewer, 'chains'):
            chain_list = sorted(self.molecular_viewer.chains.keys())
            if chain_id in chain_list:
                idx = chain_list.index(chain_id)
                return self.chain_colors[idx % len(self.chain_colors)]
        return '#FFFFFF'
    
    def _set_equal_aspect_3d(self):
        """Set equal aspect ratio for 3D plot"""
        if not self.ax or not self.molecular_viewer:
            return
        
        if not hasattr(self.molecular_viewer, 'atom_coords') or self.molecular_viewer.atom_coords is None:
            return
        
        coords = self.molecular_viewer.atom_coords
        if len(coords) > 0:
            max_range = np.array([
                coords[:, 0].max() - coords[:, 0].min(),
                coords[:, 1].max() - coords[:, 1].min(),
                coords[:, 2].max() - coords[:, 2].min()
            ]).max() / 2.0
            
            mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
            mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
            mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5
            
            # Add small padding so content doesn't touch edges
            pad = max_range * 0.05 if max_range > 0 else 1.0
            self.ax.set_xlim(mid_x - max_range - pad, mid_x + max_range + pad)
            self.ax.set_ylim(mid_y - max_range - pad, mid_y + max_range + pad)
            self.ax.set_zlim(mid_z - max_range - pad, mid_z + max_range + pad)
            
            # Minimize white-space and prevent clipping of artists
            try:
                if self.figure:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("ignore", message="The figure layout has changed to tight")
                        self.figure.tight_layout()
            except Exception:
                pass
    
    def _export_image(self):
        """Export current view as image."""
        if not self.canvas or not self.current_pdb_file:
            messagebox.showwarning("No Structure", "Please load a structure first.")
            return
        
        # Get export file path
        file_path = filedialog.asksaveasfilename(
            title="Export Image",
            defaultextension=".png",
            filetypes=FILE_FILTERS['image'],
            initialdir=self.initial_directory
        )
        
        if file_path:
            try:
                if self.status_callback:
                    self.status_callback("Exporting image...")
                
                # Save the current figure
                if self.figure:
                    self.figure.savefig(
                        file_path,
                        dpi=300,
                        bbox_inches='tight',
                        facecolor='#212121',
                        edgecolor='none'
                )
                
                if self.status_callback:
                    self.status_callback(f"Image exported: {Path(file_path).name}")
                
                messagebox.showinfo("Export Complete", f"Image saved to:\n{file_path}")
                
            except Exception as e:
                logger.error(f"Error exporting image: {e}")
                messagebox.showerror("Export Error", f"Failed to export image: {str(e)}")
    
    def _save_session(self):
        """Save session data."""
        if not self.current_pdb_file:
            messagebox.showwarning("No Structure", "Please load a structure first.")
            return
        
        # Get session file path
        file_path = filedialog.asksaveasfilename(
            title="Save Session",
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            initialdir=self.initial_directory
        )
        
        if file_path:
            try:
                if self.status_callback:
                    self.status_callback("Saving session...")
                
                import json
                
                session_data = {
                    'pdb_file': self.current_pdb_file,
                    'representation': self.repr_combo.get(),
                    'color_scheme': self.color_combo.get(),
                    'atom_size': self.atom_size_slider.get(),
                    'show_axes': self.axes_var.get(),
                    'timestamp': str(Path().cwd())
                }
                
                with open(file_path, 'w') as f:
                    json.dump(session_data, f, indent=2)
                
                if self.status_callback:
                    self.status_callback(f"Session saved: {Path(file_path).name}")
                
                messagebox.showinfo("Save Complete", f"Session saved to:\n{file_path}")
                
            except Exception as e:
                logger.error(f"Error saving session: {e}")
                messagebox.showerror("Save Error", f"Failed to save session: {str(e)}")
    
    def _set_controls_enabled(self, enabled: bool):
        """Enable or disable visualization controls."""
        state = "normal" if enabled else "disabled"
        
        self.repr_combo.configure(state=state)
        self.color_combo.configure(state=state)
        self.atom_size_slider.configure(state=state)
        self.reset_view_button.configure(state=state)
        self.center_button.configure(state=state)
        self.show_axes_checkbox.configure(state=state)
        self.export_image_button.configure(state=state)
        self.save_session_button.configure(state=state)
    
    def on_stage_shown(self):
        """Called when this stage becomes active."""
        pass
    
    def on_pdb_changed(self, pdb_file: Optional[str]):
        """Called when PDB file changes in another frame."""
        if pdb_file != self.current_pdb_file:
            self.current_pdb_file = pdb_file
            
            if pdb_file:
                self.file_entry.delete(0, "end")
                self.file_entry.insert(0, pdb_file)
                self.load_button.configure(state="normal")
                # Auto-load the structure
                self._load_structure()
            else:
                self.file_entry.delete(0, "end")
                self.load_button.configure(state="disabled")
                self._set_controls_enabled(False)
                self.info_text.configure(state="normal")
                self.info_text.delete("1.0", "end")
                self.info_text.insert("1.0", "No structure loaded")
                self.info_text.configure(state="disabled")
                
                # Clear the visualization
                if self.ax:
                    self.ax.clear()
                    self.ax.set_facecolor('#212121')
                    self._hide_axes()
                    if self.canvas:
                        self.canvas.draw()
    
    def cleanup(self):
        """Cleanup resources when frame is destroyed."""
        if self.molecular_viewer:
            try:
                self.molecular_viewer.cleanup()
            except Exception as e:
                logger.error(f"Error during molecular viewer cleanup: {e}")
            finally:
                self.molecular_viewer = None
        
        # Cleanup matplotlib resources
        if self.figure:
            try:
                import matplotlib.pyplot as plt
                plt.close(self.figure)
            except Exception as e:
                logger.error(f"Error closing matplotlib figure: {e}")
            finally:
                self.figure = None
                self.ax = None
                self.canvas = None
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the visualization frame."""
        try:
            # File section
            if hasattr(self, 'file_label'):
                self.file_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'file_entry'):
                self.file_entry.configure(font=scaled_fonts['body'])
            if hasattr(self, 'browse_button'):
                self.browse_button.configure(font=scaled_fonts['body'], width=100, height=32)
            if hasattr(self, 'load_button'):
                self.load_button.configure(font=scaled_fonts['body'], width=100, height=32)
            if hasattr(self, 'info_label'):
                self.info_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'info_text'):
                try:
                    self.info_text.configure(font=scaled_fonts['small'])
                except:
                    pass

            # 3D Structure View label (if present)
            if hasattr(self, 'viz_label'):
                self.viz_label.configure(font=scaled_fonts['heading'])

            # Controls section
            if hasattr(self, 'controls_label'):
                self.controls_label.configure(font=scaled_fonts['heading'])

            # Representation controls
            if hasattr(self, 'repr_label'):
                self.repr_label.configure(font=scaled_fonts['small'])
            if hasattr(self, 'repr_combo'):
                self.repr_combo.configure(font=scaled_fonts['small'])

            # Color controls
            if hasattr(self, 'color_label'):
                self.color_label.configure(font=scaled_fonts['small'])
            if hasattr(self, 'color_combo'):
                self.color_combo.configure(font=scaled_fonts['small'])

            # Size controls
            if hasattr(self, 'size_label'):
                self.size_label.configure(font=scaled_fonts['small'])
            if hasattr(self, 'atom_size_slider'):
                try:
                    self.atom_size_slider.configure(font=scaled_fonts['small'])
                except Exception:
                    # CTkSlider doesn't support font configuration
                    pass

            # View controls (make all buttons consistent)
            font_compatible_controls = ['reset_view_button', 'center_button', 'export_image_button', 'save_session_button']
            for btn in font_compatible_controls:
                if hasattr(self, btn):
                    try:
                        getattr(self, btn).configure(font=scaled_fonts['small'], width=100, height=32)
                    except Exception:
                        # Skip widgets that don't support font configuration
                        pass

            # Checkbox controls (handled separately due to different behavior)
            if hasattr(self, 'show_axes_checkbox'):
                try:
                    self.show_axes_checkbox.configure(font=scaled_fonts['small'])
                except Exception:
                    # Some CTk widgets may not support font configuration
                    pass

            # Note: Status is now handled by global status bar (no local status section)

        except Exception as e:
            logger.warning(f"Error updating fonts in VisualizeFrame: {e}")