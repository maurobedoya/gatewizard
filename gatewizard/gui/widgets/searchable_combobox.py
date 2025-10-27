# gatewizard/gui/widgets/searchable_combobox.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Searchable combobox widget with filtering capabilities.

This module provides an enhanced combobox widget that allows users
to search and filter through large lists of options.
"""

import tkinter as tk
from typing import List, Optional, Callable, Any

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import COLOR_SCHEME, FONTS
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class SearchableComboBox(ctk.CTkFrame):
    """
    A searchable combobox widget with filtering capabilities.
    
    This widget combines an entry field with a dropdown list that
    filters options as the user types.
    """
    
    def __init__(
        self,
        parent,
        values: Optional[List[str]] = None,
        command: Optional[Callable[[str], None]] = None,
        width: int = 150,
        height: int = 30,
        placeholder_text: str = "Search...",
        **kwargs
    ):
        """
        Initialize the searchable combobox.
        
        Args:
            parent: Parent widget
            values: List of selectable values
            command: Callback function called when selection changes
            width: Widget width
            height: Widget height
            placeholder_text: Placeholder text for the entry
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent, fg_color="transparent", **kwargs)
        
        self.values = values or []
        self.filtered_values = self.values.copy()
        self.command = command
        self.placeholder_text = placeholder_text
        self.current_value = ""
        self.dropdown_visible = False
        
        # Track user interaction to prevent unwanted dropdown opening
        self.last_user_interaction = 0
        self.focus_from_click = False
        
        # Create widgets
        self._create_widgets(width, height)
        self._setup_layout()
        self._setup_bindings()
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            # This widget uses self.after() calls for dropdown management
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")
    
    def _create_widgets(self, width: int, height: int):
        """Create the combobox widgets."""
        # Main entry field
        self.entry = ctk.CTkEntry(
            self,
            placeholder_text=self.placeholder_text,
            width=width - 30,
            height=height,
            font=FONTS['body']
        )
        
        # Dropdown button
        self.dropdown_button = ctk.CTkButton(
            self,
            text="▼",
            width=25,
            height=height,
            font=("Arial", 10),
            command=self._toggle_dropdown
        )
        
        # Dropdown frame (initially hidden)
        self.dropdown_frame = ctk.CTkFrame(
            self,
            fg_color=COLOR_SCHEME['background'],
            border_width=1,
            border_color=COLOR_SCHEME['inactive']
        )
        
        # Scrollable listbox for options
        self.listbox_frame = ctk.CTkScrollableFrame(
            self.dropdown_frame,
            height=150,
            fg_color="transparent"
        )
        
        self.option_buttons = []
        self._update_dropdown_options()
    
    def _setup_layout(self):
        """Setup the layout of widgets."""
        # Configure grid weights
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=0)
        
        # Pack main widgets
        self.entry.grid(row=0, column=0, sticky="ew")
        self.dropdown_button.grid(row=0, column=1, padx=(2, 0))
        
        # Initially hide dropdown
        self.dropdown_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(2, 0))
        self.dropdown_frame.grid_remove()
        
        self.listbox_frame.pack(fill="both", expand=True, padx=5, pady=5)
    
    def _setup_bindings(self):
        """Setup event bindings."""
        # Entry field events
        self.entry.bind("<KeyRelease>", self._on_text_changed)
        self.entry.bind("<FocusIn>", self._on_entry_focus_in)
        self.entry.bind("<Return>", self._on_return_pressed)
        self.entry.bind("<Escape>", self._hide_dropdown)
        self.entry.bind("<Down>", self._on_down_arrow)
        
        # Track user clicks to distinguish from focus restoration
        self.entry.bind("<Button-1>", self._on_entry_click)

        # Focus out event to hide dropdown
        self.entry.bind("<FocusOut>", self._on_focus_out)

    def _on_focus_out(self, event=None):
        """Handle focus out event."""
        # Small delay to allow clicking on dropdown items
        self.after(100, self._check_and_hide_dropdown)

    def _check_and_hide_dropdown(self):
        """Check if dropdown should be hidden."""
        try:
            # Get the widget that currently has focus
            focused_widget = self.focus_get()

            # Check if focus is on dropdown or its children
            if focused_widget:
                # Get the top-level parent of the focused widget
                current = focused_widget
                while current:
                    if current == self or current == self.dropdown_frame:
                        return  # Don't hide if focus is within our widget
                    try:
                        current = current.master
                    except:
                        break
                    
            # Hide dropdown if focus is elsewhere
            self._hide_dropdown()
        except:
            # If there's any error, just hide the dropdown
            self._hide_dropdown()
    
    def _on_text_changed(self, event=None):
        """Handle text change in entry field."""
        search_text = self.entry.get().lower()
        
        # Filter values based on search text
        if search_text:
            self.filtered_values = [
                value for value in self.values
                if search_text in value.lower()
            ]
        else:
            self.filtered_values = self.values.copy()
        
        # Update dropdown options
        self._update_dropdown_options()
        
        # Show dropdown if there are filtered values
        if self.filtered_values and not self.dropdown_visible:
            self._show_dropdown()
    
    def _on_down_arrow(self, event=None):
        """Handle down arrow key press to show dropdown."""
        if not self.dropdown_visible:
            self._show_dropdown()
        return "break"  # Prevent default behavior

    def _on_entry_click(self, event=None):
        """Handle entry field click."""
        import time
        self.last_user_interaction = time.time()
        self.focus_from_click = True

    def _on_entry_focus_in(self, event=None):
        """Handle entry field focus in."""
        import time
        current_time = time.time()
        
        # Only show dropdown if:
        # 1. It's not already visible, AND
        # 2. This focus event is from a recent user click (within 500ms)
        if (not self.dropdown_visible and 
            self.focus_from_click and 
            (current_time - self.last_user_interaction) < 0.5):
            self._show_dropdown()
        
        # Reset the click flag after processing
        self.focus_from_click = False
    
    def _on_return_pressed(self, event=None):
        """Handle Return key press."""
        if self.filtered_values:
            # Select first filtered option
            self._select_option(self.filtered_values[0])
        else:
            # Use current entry text
            self._select_option(self.entry.get())
    
    def _toggle_dropdown(self):
        """Toggle dropdown visibility."""
        if self.dropdown_visible:
            self._hide_dropdown()
        else:
            self._show_dropdown()
    
    def _show_dropdown(self):
        """Show the dropdown list."""
        if not self.dropdown_visible:
            self.dropdown_frame.grid()
            self.dropdown_visible = True
            self.dropdown_button.configure(text="▲")
    
    def _hide_dropdown(self, event=None):
        """Hide the dropdown list."""
        if self.dropdown_visible:
            self.dropdown_frame.grid_remove()
            self.dropdown_visible = False
            self.dropdown_button.configure(text="▼")
    
    def _update_dropdown_options(self):
        """Update the options in the dropdown list."""
        # Clear existing options
        for button in self.option_buttons:
            button.destroy()
        self.option_buttons.clear()

        # Create buttons for filtered values
        for value in self.filtered_values:
            button = ctk.CTkButton(
                self.listbox_frame,
                text=value,
                height=25,
                font=FONTS['body'],
                fg_color="transparent",
                text_color=COLOR_SCHEME['text'],
                hover_color=COLOR_SCHEME['highlight'],
                anchor="w",
                command=lambda v=value: self._select_option_with_focus(v)
            )
            button.pack(fill="x", pady=1)
            self.option_buttons.append(button)
    
    def _select_option(self, value: str):
        """Select an option and update the entry field."""
        self.current_value = value
        self.entry.delete(0, "end")
        self.entry.insert(0, value)
        self._hide_dropdown()

        # Give focus back to entry after selection
        self.entry.focus_set()

        # Call command callback if provided
        if self.command:
            try:
                self.command(value)
            except Exception as e:
                logger.error(f"Error in searchable combobox command: {e}")

    def _select_option_with_focus(self, value: str):
        """Select option and manage focus properly."""
        # Set focus to entry first to prevent focus issues
        self.entry.focus_set()
        # Then select the option
        self._select_option(value)

    def get(self) -> str:
        """Get the current value."""
        return self.entry.get()
    
    def set(self, value: str):
        """Set the current value."""
        self.current_value = value
        self.entry.delete(0, "end")
        self.entry.insert(0, value)
    
    def configure(self, require_redraw=True, **kwargs):
        """Configure widget properties."""
        # Handle our custom configuration options first
        custom_handled = False
        
        if 'values' in kwargs:
            self.values = kwargs.pop('values')
            self.filtered_values = self.values.copy()
            self._update_dropdown_options()
            custom_handled = True
        
        if 'command' in kwargs:
            self.command = kwargs.pop('command')
            custom_handled = True
        
        if 'state' in kwargs:
            state = kwargs.pop('state')
            self.entry.configure(state=state)
            self.dropdown_button.configure(state=state)
            custom_handled = True
        
        # If we handled custom attributes and there are no remaining kwargs, return early
        if custom_handled and not kwargs:
            return
        
        # Otherwise, delegate to parent
        return super().configure(require_redraw, **kwargs)
    
    def focus(self):
        """Set focus to the entry field."""
        self.entry.focus()
    
    def focus_set(self):
        """Set focus to the entry field."""
        self.entry.focus_set()

    def update_fonts(self, scaled_fonts):
        """Update fonts for the combobox components."""
        try:
            # Entry and dropdown button
            self.entry.configure(font=scaled_fonts.get('body', FONTS['body']))
            try:
                self.dropdown_button.configure(font=scaled_fonts.get('small', FONTS['small']))
            except Exception:
                pass
            # Option buttons in dropdown
            for btn in getattr(self, 'option_buttons', []) or []:
                try:
                    btn.configure(font=scaled_fonts.get('body', FONTS['body']))
                except Exception:
                    pass
        except Exception as e:
            logger.debug(f"SearchableComboBox.update_fonts error: {e}")