# gatewizard/gui/widgets/leaflet_frame.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Leaflet frame widget for lipid selection.

This module provides a widget for selecting lipids for membrane leaflets
in system preparation.
"""

from typing import List, Callable, Optional

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import COLOR_SCHEME, FONTS, WIDGET_SIZES, LAYOUT
from gatewizard.gui.widgets.searchable_combobox import SearchableComboBox
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class LeafletFrame(ctk.CTkFrame):
    """
    Widget for selecting lipids for a membrane leaflet.
    
    This widget allows users to select multiple lipids for either
    upper or lower membrane leaflets.
    """
    
    def __init__(
        self,
        parent,
        title: str = "Leaflet",
        available_lipids: Optional[List[str]] = None,
        max_lipids: int = 5
    ):
        """
        Initialize the leaflet frame.
        
        Args:
            parent: Parent widget
            title: Title for the leaflet (e.g., "Upper Leaflet")
            available_lipids: List of available lipid names
            max_lipids: Maximum number of lipids that can be selected
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.title = title
        self.available_lipids = available_lipids or []
        self.max_lipids = max_lipids
        self.selected_lipids = []
        self.lipid_widgets = []
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
    
    def _create_widgets(self):
        """Create widgets for lipid selection."""
        # Title label
        self.title_label = ctk.CTkLabel(
            self,
            text=self.title,
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        # Lipid selection area
        self.lipids_frame = ctk.CTkFrame(self, fg_color="transparent")
        
        # Add lipid button
        self.add_button = ctk.CTkButton(
            self,
            text="+ Add Lipid",
            width=WIDGET_SIZES['button_width'],
            height=WIDGET_SIZES['button_height'],
            command=self._add_lipid_selection
        )
        
        # Initially add one lipid selection
        self._add_lipid_selection()
    
    def _setup_layout(self):
        """Setup the layout of widgets."""
        self.title_label.pack(anchor="w", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.lipids_frame.pack(fill="both", expand=True, padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
        
        self.add_button.pack(padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
    
    def _add_lipid_selection(self):
        """Add a new lipid selection widget."""
        if len(self.lipid_widgets) >= self.max_lipids:
            return
        
        # Create frame for this lipid selection
        lipid_frame = ctk.CTkFrame(self.lipids_frame, fg_color="transparent")
        lipid_frame.pack(fill="x", pady=LAYOUT['padding_small'])
        
        # Lipid combobox
        lipid_combo = SearchableComboBox(
            lipid_frame,
            values=self.available_lipids,
            width=WIDGET_SIZES['combobox_width'],
            height=WIDGET_SIZES['combobox_height'],
            command=self._on_lipid_changed
        )
        lipid_combo.pack(side="left", padx=(0, LAYOUT['padding_small']))
        
        # Set POPC as default for the first lipid
        if len(self.lipid_widgets) == 0 and "POPC" in self.available_lipids:
            lipid_combo.set("POPC")
        
        # Ratio label
        ratio_label = ctk.CTkLabel(
            lipid_frame,
            text="Ratio:",
            font=FONTS['small'],
            width=40
        )
        ratio_label.pack(side="left", padx=(LAYOUT['padding_small'], 2))
        
        # Ratio entry
        ratio_entry = ctk.CTkEntry(
            lipid_frame,
            placeholder_text="1.0",
            width=60,
            height=WIDGET_SIZES['combobox_height']
        )
        ratio_entry.pack(side="left", padx=(0, LAYOUT['padding_small']))
        ratio_entry.insert(0, "1.0")  # Default ratio
        
        # Remove button
        remove_button = ctk.CTkButton(
            lipid_frame,
            text="×",
            width=30,
            height=WIDGET_SIZES['combobox_height'],
            command=lambda: self._remove_lipid_selection(lipid_frame)
        )
        remove_button.pack(side="left")
        
        # Store references
        lipid_widget = {
            'frame': lipid_frame,
            'combo': lipid_combo,
            'ratio_entry': ratio_entry,
            'button': remove_button
        }
        self.lipid_widgets.append(lipid_widget)
        
        # Update add button state
        self._update_add_button_state()
        
        # Update selected lipids immediately for default POPC
        if len(self.lipid_widgets) == 1 and lipid_combo.get() == "POPC":
            self._update_selected_lipids()
    
    def _remove_lipid_selection(self, frame_to_remove):
        """Remove a lipid selection widget."""
        # Find and remove the widget
        for i, widget in enumerate(self.lipid_widgets):
            if widget['frame'] == frame_to_remove:
                widget['frame'].destroy()
                self.lipid_widgets.pop(i)
                break
        
        # Update selected lipids
        self._update_selected_lipids()
        
        # Update add button state
        self._update_add_button_state()
        
        # Ensure at least one selection remains
        if not self.lipid_widgets:
            self._add_lipid_selection()
    
    def _on_lipid_changed(self, value=None):
        """Handle lipid selection change."""
        self._update_selected_lipids()
    
    def _update_selected_lipids(self):
        """Update the list of selected lipids."""
        self.selected_lipids = []
        
        for widget in self.lipid_widgets:
            lipid = widget['combo'].get().strip()
            if lipid and lipid in self.available_lipids:
                self.selected_lipids.append(lipid)
        
        logger.debug(f"{self.title} selected lipids: {self.selected_lipids}")
    
    def _update_add_button_state(self):
        """Update the state of the add button."""
        if len(self.lipid_widgets) >= self.max_lipids:
            self.add_button.configure(state="disabled")
        else:
            self.add_button.configure(state="normal")
    
    def get_selected_lipids(self) -> List[str]:
        """
        Get the list of selected lipids.
        
        Returns:
            List of selected lipid names
        """
        self._update_selected_lipids()
        return self.selected_lipids.copy()
    
    def get_lipid_ratios(self) -> List[float]:
        """
        Get the list of lipid ratios corresponding to selected lipids.
        
        Returns:
            List of ratios as floats
        """
        ratios = []
        for widget in self.lipid_widgets:
            lipid = widget['combo'].get().strip()
            if lipid and lipid in self.available_lipids:
                try:
                    ratio = float(widget['ratio_entry'].get().strip())
                    ratios.append(ratio)
                except (ValueError, AttributeError):
                    ratios.append(1.0)  # Default ratio if invalid
        return ratios
    
    def get_lipids_with_ratios(self) -> List[tuple[str, float]]:
        """
        Get lipids paired with their ratios.
        
        Returns:
            List of (lipid_name, ratio) tuples
        """
        lipids_ratios = []
        for widget in self.lipid_widgets:
            lipid = widget['combo'].get().strip()
            if lipid and lipid in self.available_lipids:
                try:
                    ratio = float(widget['ratio_entry'].get().strip())
                except (ValueError, AttributeError):
                    ratio = 1.0
                lipids_ratios.append((lipid, ratio))
        return lipids_ratios
    
    def set_selected_lipids(self, lipids: List[str], ratios: Optional[List[float]] = None):
        """
        Set the selected lipids and their ratios.
        
        Args:
            lipids: List of lipid names to select
            ratios: Optional list of ratios corresponding to lipids
        """
        # Clear existing selections
        for widget in self.lipid_widgets:
            widget['frame'].destroy()
        self.lipid_widgets.clear()
        
        # Add selections for each lipid
        for i, lipid in enumerate(lipids):
            self._add_lipid_selection()
            if self.lipid_widgets:
                self.lipid_widgets[-1]['combo'].set(lipid)
                # Set ratio if provided
                if ratios and i < len(ratios):
                    self.lipid_widgets[-1]['ratio_entry'].delete(0, 'end')
                    self.lipid_widgets[-1]['ratio_entry'].insert(0, str(ratios[i]))
        
        # Ensure at least one selection exists with POPC default
        if not self.lipid_widgets:
            self._add_lipid_selection()
        
        self._update_selected_lipids()
        self._update_add_button_state()
    
    def clear_selections(self):
        """Clear all lipid selections."""
        self.set_selected_lipids([])
    
    def set_available_lipids(self, lipids: List[str]):
        """
        Update the list of available lipids.
        
        Args:
            lipids: List of available lipid names
        """
        self.available_lipids = lipids
        
        # Update all comboboxes
        for widget in self.lipid_widgets:
            widget['combo'].configure(values=lipids)
    
    def validate_selections(self) -> tuple[bool, str]:
        """
        Validate the current lipid selections and ratios.
        
        Returns:
            Tuple of (is_valid, error_message)
        """
        selected = self.get_selected_lipids()
        ratios = self.get_lipid_ratios()
        
        if not selected:
            return False, f"No lipids selected for {self.title.lower()}"
        
        # Check for duplicates
        if len(selected) != len(set(selected)):
            return False, f"Duplicate lipids selected in {self.title.lower()}"
        
        # Check if all selected lipids are valid
        invalid_lipids = [lipid for lipid in selected if lipid not in self.available_lipids]
        if invalid_lipids:
            return False, f"Invalid lipids in {self.title.lower()}: {', '.join(invalid_lipids)}"
        
        # Validate ratios
        if len(ratios) != len(selected):
            return False, f"Mismatch between lipids and ratios in {self.title.lower()}"
        
        # Check for valid ratio values
        for i, ratio in enumerate(ratios):
            if ratio <= 0:
                return False, f"Invalid ratio for {selected[i]} in {self.title.lower()}: must be positive"
        
        return True, ""
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the leaflet frame."""
        try:
            # Update title label
            if hasattr(self, 'title_label'):
                self.title_label.configure(font=scaled_fonts['heading'])
            
            # Update add button
            if hasattr(self, 'add_button'):
                self.add_button.configure(font=scaled_fonts['small'])
            
            # Update lipid widgets (these are created dynamically)
            if hasattr(self, 'lipid_widgets'):
                for widget_dict in self.lipid_widgets:
                    if 'remove_button' in widget_dict:
                        try:
                            widget_dict['remove_button'].configure(font=scaled_fonts['small'])
                        except:
                            pass
            
        except Exception as e:
            logger.warning(f"Error updating fonts in LeafletFrame: {e}")