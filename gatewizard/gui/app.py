# gatewizard/gui/app.py

"""
Main application window for Gatewizard GUI.

This module contains the main application class that coordinates all GUI
components and manages the overall application state.
"""

import sys
import os
import tkinter as tk
from typing import Optional
from pathlib import Path

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI. Install with: pip install customtkinter")

from gatewizard import __version__
from gatewizard.gui.constants import (
    COLOR_SCHEME, WINDOW_CONFIG, FONTS, FONT_SCALE_OPTIONS, get_scaled_fonts
)
from gatewizard.gui.widgets.stage_tabs import StageTabsContainer
from gatewizard.gui.frames.visualize import VisualizeFrame
from gatewizard.gui.frames.preparation_frame import PreparationFrame
from gatewizard.gui.frames.builder_frame import BuilderFrame
from gatewizard.gui.frames.analysis import AnalysisFrame
from gatewizard.gui.frames.equilibration import EquilibrationFrame
from gatewizard.utils.config import load_config, save_config, get_config_value, set_config_value
from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class ProteinViewerApp(ctk.CTk):
    """
    Main application window for Gatewizard.
    
    This class manages the overall GUI layout, coordinates between different
    frames, and handles application-level events and state.
    """
    
    def __init__(self, target_screen: int = 0):
        """
        Initialize the main application window.
        
        Args:
            target_screen: Target screen number for multi-monitor setups
        """
        super().__init__()
        
        # Initialize configuration
        self.config = load_config()
        self.target_screen = target_screen
        
        # Capture the initial working directory where the app was launched
        self.initial_working_directory = str(Path.cwd())
        
        # Application state
        self.current_pdb_file = None
        self.current_stage = "Visualize"
        
        # Setup GUI
        self._setup_window()
        self._setup_theme()
        self._create_widgets()
        self._setup_layout()
        self._setup_bindings()
        
        # Initialize stages
        self._initialize_stages()
        
        # Apply saved font scaling
        self._apply_saved_font_scaling()
        
        # Initialize position saving timer
        self._position_save_timer = None
        self._scheduled_callbacks = []  # Track all scheduled callbacks for cleanup
        self._shutting_down = False  # Flag to prevent new callbacks during shutdown
        
        # Add emergency recovery system
        self._setup_emergency_recovery()
        
        # Setup safe callback system
        self._setup_safe_callback_system()
        
        # Check if initial window size is appropriate (after a short delay)
        self._schedule_callback(500, self._check_initial_window_size)
        
        logger.info("Gatewizard GUI initialized successfully")
    
    def _schedule_callback(self, delay, callback, *args):
        """Schedule a callback and track it for cleanup with automatic safety checks."""
        # Don't schedule new callbacks if we're shutting down
        if getattr(self, '_shutting_down', False):
            return None
        
        # Wrap the callback to check if widget still exists
        def safe_callback():
            try:
                # Check if app still exists and not shutting down
                if not getattr(self, '_shutting_down', False) and self.winfo_exists():
                    callback(*args)
            except Exception as e:
                # Silently ignore errors during shutdown
                if not getattr(self, '_shutting_down', False):
                    logger.debug(f"Error in scheduled callback: {e}")
            
        try:
            callback_id = self.after(delay, safe_callback)
            self._scheduled_callbacks.append(callback_id)
            return callback_id
        except Exception:
            # If scheduling fails (e.g., during shutdown), return None
            return None
    
    def report_callback_exception(self, exc, val, tb):  # type: ignore
        """Log uncaught exceptions from Tkinter callbacks with full traceback.

        This helps diagnose issues that otherwise surface as silent
        'Exception in Tkinter callback' messages without context.
        Filters out harmless cleanup messages during app shutdown.
        """
        try:
            # Filter out harmless "invalid command name" errors during app shutdown
            error_str = str(val) if val else str(exc)
            
            # Expanded list of harmless callback cleanup patterns
            harmless_patterns = [
                "invalid command name",
                "_on_closing",
                "_check_window_state", 
                "_validate_window_geometry",
                "_check_initial_window_size",
                "_save_window_position",
                "_restore_safe_geometry",
                "check_dpi_scaling",
                "update",
                "<lambda>",
                "_refresh_cycle",
                "_auto_update_logs",
                "_check_and_hide_dropdown",
                "_set_grab",
                "_update_propka_status",
                "_display_propka_results",
                "_show_propka_error"
            ]
            
            # Check if this is a harmless cleanup message
            is_harmless = any(pattern in error_str for pattern in harmless_patterns)
            
            if is_harmless:
                # These are harmless cleanup messages - just log as debug
                logger.debug(f"Harmless callback cleanup message: {error_str}")
                return
            
            # Log real errors normally
            logger.error("Tkinter callback error", exc_info=(exc, val, tb))
        except Exception:
            # As a fallback, print to stderr to avoid losing the error
            import traceback, sys as _sys
            print("Tkinter callback error:", file=_sys.stderr)
            traceback.print_exception(exc, val, tb)

    def _apply_saved_font_scaling(self):
        """Apply saved font scaling on startup."""
        try:
            # Default to Large (1.2) if not set
            scale_factor = get_config_value('gui.font_scale', 1.2)
            # Always apply on startup to sync all frames
            self._update_fonts(scale_factor)
        except Exception as e:
            logger.warning(f"Could not apply saved font scaling: {e}")
    
    def _setup_window(self):
        """Setup main window properties."""
        self.title(f"GateWizard v{__version__} - Membrane Protein Preparation Tool")

        # Get screen dimensions first to validate window size
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        
        # Calculate safe maximum window size (leave margin for taskbar, window decorations)
        margin = 100  # Space for taskbar and window decorations
        max_safe_width = screen_width - margin
        max_safe_height = screen_height - margin
        
        # Get window dimensions from config
        desired_width = get_config_value('gui.window_width', WINDOW_CONFIG['default_width'])
        desired_height = get_config_value('gui.window_height', WINDOW_CONFIG['default_height'])

        # Validate and correct problematic saved dimensions
        if desired_width > max_safe_width or desired_height > max_safe_height:
            logger.warning(f"Saved window size ({desired_width}x{desired_height}) too large for screen, correcting...")
            # Save corrected values back to config to prevent future issues
            try:
                config = load_config()
                if desired_width > max_safe_width:
                    config.gui.window_width = min(desired_width, max_safe_width)
                    desired_width = config.gui.window_width
                if desired_height > max_safe_height:
                    config.gui.window_height = min(desired_height, max_safe_height)
                    desired_height = config.gui.window_height
                save_config(config)
                logger.info(f"Updated config with corrected window size: {desired_width}x{desired_height}")
            except Exception as e:
                logger.warning(f"Could not save corrected config: {e}")

        # Validate that the desired size fits on the screen
        width = min(desired_width, max_safe_width)
        height = min(desired_height, max_safe_height)
        
        # Also ensure we meet minimum requirements
        width = max(width, WINDOW_CONFIG['min_width'])
        height = max(height, WINDOW_CONFIG['min_height'])
        
        # If minimum size is still too big for screen, use smaller values
        if width > max_safe_width:
            width = max(800, max_safe_width)  # Fallback to 800 or screen size
        if height > max_safe_height:
            height = max(600, max_safe_height)  # Fallback to 600 or screen size
            
        logger.info(f"Setting window size to {width}x{height} (screen: {screen_width}x{screen_height})")

        # Set window size and constraints
        self.geometry(f"{width}x{height}")
        
        # Set minimum size constraints, but don't make them larger than what fits on screen
        min_width = min(WINDOW_CONFIG['min_width'], max_safe_width)
        min_height = min(WINDOW_CONFIG['min_height'], max_safe_height)
        self.minsize(min_width, min_height)

        # Position window on target screen AFTER setting size
        self._position_on_screen()

        # Setup window closing behavior
        self.protocol("WM_DELETE_WINDOW", self._on_closing)
        
        # Add window state protocol handler for maximize/minimize events
        try:
            # Try to bind window state change protocol (may not work on all platforms)
            self.bind("<Control-Key-F10>", lambda e: self._handle_maximize_request())
        except Exception:
            pass  # Not all window managers support this
    
    def _position_on_screen(self):
        """Position window on the target screen."""
        try:
            # Update the window to get accurate geometry
            self.update_idletasks()

            # Get window dimensions
            window_width = self.winfo_reqwidth()
            window_height = self.winfo_reqheight()

            # Get screen dimensions (total virtual desktop)
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()

            # Simple centering approach - avoid complex multi-monitor detection for initial positioning
            if self.target_screen == 0 or self.target_screen is None:
                # Center on primary screen (assume left half if ultra-wide or multi-monitor)
                primary_width = screen_width // 2 if screen_width > 2000 else screen_width
                x = (primary_width - window_width) // 2
                y = (screen_height - window_height) // 2
            else:
                # For secondary screen, use a safe offset approach
                # Assume secondary screen is to the right
                primary_width = screen_width // 2 if screen_width > 2000 else 1920  # Common primary width
                x = primary_width + 100  # Offset into second screen
                y = 100

            # Ensure window stays within safe bounds (leave margin for taskbar/decorations)
            margin = 50
            max_x = screen_width - window_width - margin
            max_y = screen_height - window_height - margin
            
            x = max(margin, min(x, max_x))
            y = max(margin, min(y, max_y))

            # Check for saved position (but validate it's safe and reasonable)
            if get_config_value('gui.remember_window_position', True):
                saved_x = get_config_value('gui.last_window_x')
                saved_y = get_config_value('gui.last_window_y')
                saved_width = get_config_value('gui.window_width')
                saved_height = get_config_value('gui.window_height')

                if (saved_x is not None and saved_y is not None and 
                    saved_width is not None and saved_height is not None):
                    
                    # Validate saved position is reasonable
                    # Allow some tolerance for window decorations and multi-monitor differences
                    tolerance = 100
                    if (saved_x > -tolerance and saved_x < screen_width - tolerance and
                        saved_y > -tolerance and saved_y < screen_height - tolerance and
                        saved_width > 200 and saved_height > 200 and
                        saved_width < screen_width + tolerance and 
                        saved_height < screen_height + tolerance):
                        
                        x = saved_x
                        y = saved_y
                        # Also restore saved size
                        self.geometry(f"{saved_width}x{saved_height}")

            # Apply position
            self.geometry(f"+{x}+{y}")
            
            # Add a small delay and then validate the position actually worked
            self._schedule_callback(200, self._validate_final_position)

        except Exception as e:
            logger.warning(f"Could not position window: {e}")
            # Fallback to safe default positioning
            self._restore_safe_geometry()
    
    def _validate_final_position(self):
        """Validate that the window ended up in a reasonable position after initial positioning."""
        try:
            # Give the window manager time to position the window
            actual_x = self.winfo_x()
            actual_y = self.winfo_y()
            actual_width = self.winfo_width()
            actual_height = self.winfo_height()
            
            # If window ended up in a clearly wrong position, fix it
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            
            if (actual_x < -500 or actual_x > screen_width + 500 or
                actual_y < -500 or actual_y > screen_height + 500):
                logger.warning(f"Window positioned incorrectly at {actual_x},{actual_y}, correcting...")
                self._restore_safe_geometry()
                return
            
            # Also check if window is too large for the screen
            margin = 80
            if (actual_width > screen_width - margin or actual_height > screen_height - margin):
                logger.warning(f"Window size too large ({actual_width}x{actual_height}) for screen ({screen_width}x{screen_height}), adjusting...")
                safe_width = min(actual_width, screen_width - margin)
                safe_height = min(actual_height, screen_height - margin)
                self.geometry(f"{safe_width}x{safe_height}")
                
        except Exception as e:
            logger.warning(f"Error validating final position: {e}")
    
    def _setup_theme(self):
        """Setup application theme and colors."""
        # Set appearance mode and color theme
        theme = get_config_value('gui.theme', 'dark')
        color_scheme = get_config_value('gui.color_scheme', 'blue')
        
        ctk.DrawEngine.preferred_drawing_method = "circle_shapes"
        ctk.set_appearance_mode(theme)
        ctk.set_default_color_theme(color_scheme)
        
        # Configure colors
        self.configure(fg_color=COLOR_SCHEME['background'])
    
    def _create_widgets(self):
        """Create main GUI widgets."""
        # Create main container
        self.main_container = ctk.CTkFrame(self, fg_color=COLOR_SCHEME['background'])

        # Create title frame
        self.title_frame = ctk.CTkFrame(
            self.main_container,
            fg_color=COLOR_SCHEME['canvas'],
            height=60
        )

        # Create title content frame for better layout
        self.title_content_frame = ctk.CTkFrame(
            self.title_frame,
            fg_color="transparent"
        )

        # Title label
        self.title_label = ctk.CTkLabel(
            self.title_content_frame,
            text=f"GateWizard",
            font=FONTS['title'],
            text_color=COLOR_SCHEME['text']
        )

        # App descriptor placed inline to the right of the title
        self.app_descriptor_label = ctk.CTkLabel(
            self.title_content_frame,
            text="",
            font=FONTS['subtitle'],
            text_color=COLOR_SCHEME['text']
        )

        # Settings button (top right)
        self.settings_button = ctk.CTkButton(
            self.title_frame,
            text="Settings",
            font=FONTS['small'],
            width=80,
            height=28,
            command=self._show_settings,
            fg_color=COLOR_SCHEME['buttons'],
            hover_color=COLOR_SCHEME['hover']
        )

        # Create stage tabs container
        self.stage_tabs = StageTabsContainer(
            self.main_container,
            stage_changed_callback=self._on_stage_changed
        )

        # Create content frame
        self.content_frame = ctk.CTkFrame(
            self.main_container,
            fg_color=COLOR_SCHEME['content_bg']
        )

        # CREATE STATUS BAR FIRST - before stage frames
        self.status_frame = ctk.CTkFrame(
            self.main_container,
            fg_color=COLOR_SCHEME['canvas'],
            height=30
        )

        self.status_label = ctk.CTkLabel(
            self.status_frame,
            text="Ready",
            font=FONTS['small'],
            text_color=COLOR_SCHEME['text']
        )

        # NOW create frames for each stage (after status_label exists)
        self._create_stage_frames()
    
    def _create_stage_frames(self):
        """Create frames for each analysis stage."""
        self.stage_frames = {}
        
        # Visualize frame
        self.stage_frames['Visualize'] = VisualizeFrame(
            self.content_frame,
            pdb_changed_callback=self._on_pdb_changed,
            status_callback=self._update_status,
            initial_directory=self.initial_working_directory
        )
        
        # Preparation frame
        self.stage_frames['Preparation'] = PreparationFrame(
            self.content_frame,
            get_current_pdb=self._get_current_pdb,
            status_callback=self._update_status,
            initial_directory=self.initial_working_directory
        )
        
        # Builder frame
        self.stage_frames['Builder'] = BuilderFrame(
            self.content_frame,
            get_current_pdb=self._get_current_pdb,
            status_callback=self._update_status,
            initial_directory=self.initial_working_directory
        )
        
        # Equilibration frame
        self.stage_frames['Equilibration'] = EquilibrationFrame(
            self.content_frame,
            get_current_pdb=self._get_current_pdb,
            status_callback=self._update_status,
            initial_directory=self.initial_working_directory
        )
        
        # Collective Variables frame
        self.stage_frames['Analysis'] = AnalysisFrame(
            self.content_frame,
            get_current_pdb=self._get_current_pdb, # type: ignore
            status_callback=self._update_status,
            initial_directory=self.initial_working_directory
        )
        
        # Initially hide all frames
        for frame in self.stage_frames.values():
            frame.pack_forget()
    
    def _setup_layout(self):
        """Setup the main window layout."""
        # Pack main container
        self.main_container.pack(fill="both", expand=True, padx=5, pady=5)

        # Pack title frame at the top
        self.title_frame.pack(fill="x", padx=5, pady=(5, 0))
        self.title_frame.pack_propagate(False)

        # Title content inline: app name + descriptor on one row
        self.title_content_frame.pack(side="left", fill="both", expand=True)
        self.title_label.pack(side="left", padx=(20, 0), pady=(10, 10))
        self.app_descriptor_label.pack(side="left", padx=(12, 0), pady=(10, 10))

        # Settings button on the right
        self.settings_button.pack(side="right", padx=10, pady=15)

        # Stage tabs under the title bar
        self.stage_tabs.pack(fill="x", padx=5, pady=5)

        # Main content below tabs - no bottom padding to preserve status bar space
        self.content_frame.pack(fill="both", expand=True, padx=5, pady=(5, 5)) # < this mk setting was problematic.

        # Status bar at the bottom
        self.status_frame.pack(fill="x", padx=5, pady=(0, 5))
        self.status_frame.pack_propagate(False)
        self.status_label.pack(side="left", padx=10, pady=5)
    
    def _setup_bindings(self):
        """Setup event bindings."""
        # Window resize/move tracking
        self.bind("<Configure>", self._on_window_configure)
        
        # Window state change tracking (for maximize/minimize/restore)
        self.bind("<FocusIn>", self._on_window_state_change)
        self.bind("<Map>", self._on_window_state_change)
        self.bind("<Unmap>", self._on_window_state_change)
        
        # Keyboard shortcuts 
        self.bind("<Control-o>", self._handle_open_file_shortcut)
        self.bind("<Control-q>", lambda e: self._on_closing())
        self.bind("<F1>", lambda e: self._show_help())
        
        # Stage switching shortcuts - multiple alternatives for different keyboards
        # Main number row (top of keyboard)
        self.bind("<Control-1>", self._handle_stage_shortcut_1)
        self.bind("<Control-2>", self._handle_stage_shortcut_2)
        self.bind("<Control-3>", self._handle_stage_shortcut_3)
        self.bind("<Control-4>", self._handle_stage_shortcut_4)
        self.bind("<Control-5>", self._handle_stage_shortcut_5)
        
        # Numeric keypad
        self.bind("<Control-KP_1>", self._handle_stage_shortcut_1)
        self.bind("<Control-KP_2>", self._handle_stage_shortcut_2)
        self.bind("<Control-KP_3>", self._handle_stage_shortcut_3)
        self.bind("<Control-KP_4>", self._handle_stage_shortcut_4)
        self.bind("<Control-KP_5>", self._handle_stage_shortcut_5)
        
        # Alternative key codes for international keyboards
        self.bind("<Control-Key-1>", self._handle_stage_shortcut_1)
        self.bind("<Control-Key-2>", self._handle_stage_shortcut_2)
        self.bind("<Control-Key-3>", self._handle_stage_shortcut_3)
        self.bind("<Control-Key-4>", self._handle_stage_shortcut_4)
        self.bind("<Control-Key-5>", self._handle_stage_shortcut_5)
        
        # Alternative with keysym names
        self.bind("<Control-exclam>", self._handle_stage_shortcut_1)      # Shift+1 on some layouts
        self.bind("<Control-at>", self._handle_stage_shortcut_2)          # Shift+2 on some layouts
        self.bind("<Control-numbersign>", self._handle_stage_shortcut_3)  # Shift+3 on some layouts
        self.bind("<Control-dollar>", self._handle_stage_shortcut_4)      # Shift+4 on some layouts
        self.bind("<Control-percent>", self._handle_stage_shortcut_5)     # Shift+5 on some layouts
        
        # Function key alternatives (F2-F6) for more reliable stage switching
        self.bind("<F2>", self._handle_stage_shortcut_1)   # F2 = Visualize
        self.bind("<F3>", self._handle_stage_shortcut_2)   # F3 = Preparation
        self.bind("<F4>", self._handle_stage_shortcut_3)   # F4 = Builder
        self.bind("<F5>", self._handle_stage_shortcut_4)   # F5 = Equilibration
        self.bind("<F6>", self._handle_stage_shortcut_5)   # F6 = Analysis
        
        # Alt+number combinations as another alternative
        self.bind("<Alt-1>", self._handle_stage_shortcut_1)
        self.bind("<Alt-2>", self._handle_stage_shortcut_2)
        self.bind("<Alt-3>", self._handle_stage_shortcut_3)
        self.bind("<Alt-4>", self._handle_stage_shortcut_4)
        self.bind("<Alt-5>", self._handle_stage_shortcut_5)
        
        # Font scaling shortcuts
        self.bind("<Control-plus>", self._increase_font_size)
        self.bind("<Control-equal>", self._increase_font_size)  # For keyboards without separate +
        self.bind("<Control-minus>", self._decrease_font_size)
        
        # Numeric keypad font scaling
        self.bind("<Control-KP_Add>", self._increase_font_size)      # Numeric keypad plus
        self.bind("<Control-KP_Subtract>", self._decrease_font_size)  # Numeric keypad minus
        
        # Allow window to receive focus for shortcuts, but don't steal focus from input widgets
        # Only set focus when clicking on non-input areas
        self.bind("<Button-1>", self._handle_main_window_click)
        
        # Debug key presses to help diagnose keyboard issues (temporarily disabled to fix input widgets)
        # self.bind("<KeyPress>", self._debug_keypress)
        
        # Shortcut tester dialog
        self.bind("<Control-Shift-T>", lambda e: self._create_shortcut_test_dialog())
    
    def _initialize_stages(self):
        """Initialize the first stage."""
        # Show the visualize frame by default and trigger proper stage change
        self._on_stage_changed("Visualize")
        self.stage_tabs.set_active_stage("Visualize")
    
    def _on_stage_changed(self, stage_name: str):
        """
        Handle stage change events.
        
        Args:
            stage_name: Name of the new active stage
        """
        logger.debug(f"Stage changed to: {stage_name}")
        self.current_stage = stage_name
        self._show_stage_frame(stage_name)
        self._update_status(f"Switched to {stage_name} stage")
    
    def _show_stage_frame(self, stage_name: str):
        """
        Show the frame for the specified stage.
        
        Args:
            stage_name: Name of the stage to show
        """
        # Hide all frames
        for frame in self.stage_frames.values():
            frame.pack_forget()
        
        # Show the requested frame
        if stage_name in self.stage_frames:
            self.stage_frames[stage_name].pack(fill="both", expand=True, padx=10, pady=(10, 0))
            
            # Update frame if it has an update method
            frame = self.stage_frames[stage_name]
            if hasattr(frame, 'on_stage_shown'):
                frame.on_stage_shown()
    
    def _on_pdb_changed(self, pdb_file: Optional[str]):
        """
        Handle PDB file change events.
        
        Args:
            pdb_file: Path to the new PDB file (None if cleared)
        """
        self.current_pdb_file = pdb_file
        logger.info(f"PDB file changed: {pdb_file}")
        
        # Notify all frames about the PDB change
        for frame in self.stage_frames.values():
            if hasattr(frame, 'on_pdb_changed'):
                frame.on_pdb_changed(pdb_file)
        
        # Update status
        if pdb_file:
            filename = Path(pdb_file).name
            self._update_status(f"Loaded: {filename}")
        else:
            self._update_status("No PDB file loaded")
    
    def _get_current_pdb(self) -> Optional[str]:
        """
        Get the currently loaded PDB file.
        
        Returns:
            Path to current PDB file or None
        """
        return self.current_pdb_file
    
    def _update_status(self, message: str):
        """Update the status bar message and notify current frame."""
        # Update global status bar
        if hasattr(self, 'status_label') and self.status_label:
            self.status_label.configure(text=message)
            logger.debug(f"Status: {message}")
        else:
            # Fallback - just log if status_label doesn't exist yet
            logger.info(f"Status (no GUI): {message}")
        
        # Also update the current frame's local status display if it has one
        if hasattr(self, 'current_stage') and hasattr(self, 'stage_frames'):
            current_frame = self.stage_frames.get(self.current_stage)
            if current_frame and hasattr(current_frame, '_safe_status_callback'):
                try:
                    # Call the frame's local status update (but don't call global again)
                    if hasattr(current_frame, 'status_text_label') and current_frame.status_text_label:
                        current_frame.status_text_label.configure(text=message)
                except Exception as e:
                    logger.warning(f"Failed to update frame status: {e}")
    
    def _open_file_dialog(self):
        """Open file dialog for PDB selection."""
        if hasattr(self.stage_frames['Visualize'], 'open_file_dialog'):
            self.stage_frames['Visualize'].open_file_dialog()
    
    def _handle_open_file_shortcut(self, event):
        """Handle Ctrl+O shortcut - only works in Visualize tab."""
        logger.info(f"Ctrl+O pressed - current stage: {self.current_stage}")
        if self.current_stage == "Visualize":
            logger.info("Opening file dialog")
            self._open_file_dialog()
            return "break"  # Prevent event propagation
        else:
            # Show a helpful message
            logger.info(f"Ctrl+O blocked - not in Visualize tab")
            self._update_status(f"Ctrl+O only works in Visualize tab (currently in {self.current_stage})")
            return "break"
    
    def _handle_stage_shortcut_1(self, event):
        """Handle Ctrl+1 shortcut - switch to Visualize."""
        try:
            logger.info("Ctrl+1 pressed - switching to Visualize")
            # Use the stage change callback to properly switch tabs and frames
            self._on_stage_changed("Visualize")
            # Also update the button appearance
            self.stage_tabs.set_active_stage("Visualize")
            logger.debug("Switched to Visualize via Ctrl+1")
        except Exception as e:
            logger.error(f"Error switching to Visualize: {e}")
        return "break"
    
    def _handle_stage_shortcut_2(self, event):
        """Handle Ctrl+2 shortcut - switch to Preparation."""
        try:
            logger.info("Ctrl+2 pressed - switching to Preparation")
            # Use the stage change callback to properly switch tabs and frames
            self._on_stage_changed("Preparation")
            # Also update the button appearance
            self.stage_tabs.set_active_stage("Preparation")
            logger.debug("Switched to Preparation via Ctrl+2")
        except Exception as e:
            logger.error(f"Error switching to Preparation: {e}")
        return "break"
    
    def _handle_stage_shortcut_3(self, event):
        """Handle Ctrl+3 shortcut - switch to Builder."""
        try:
            logger.info("Ctrl+3 pressed - switching to Builder")
            # Use the stage change callback to properly switch tabs and frames
            self._on_stage_changed("Builder")
            # Also update the button appearance
            self.stage_tabs.set_active_stage("Builder")
            logger.debug("Switched to Builder via Ctrl+3")
        except Exception as e:
            logger.error(f"Error switching to Builder: {e}")
        return "break"
    
    def _handle_stage_shortcut_4(self, event):
        """Handle Ctrl+4 shortcut - switch to Equilibration."""
        try:
            logger.info("Ctrl+4 pressed - switching to Equilibration")
            # Use the stage change callback to properly switch tabs and frames
            self._on_stage_changed("Equilibration")
            # Also update the button appearance
            self.stage_tabs.set_active_stage("Equilibration")
            logger.debug("Switched to Equilibration via Ctrl+4")
        except Exception as e:
            logger.error(f"Error switching to Equilibration: {e}")
        return "break"
    
    def _handle_stage_shortcut_5(self, event):
        """Handle Ctrl+5 shortcut - switch to Analysis."""
        try:
            logger.info("Ctrl+5 pressed - switching to Analysis")
            # Use the stage change callback to properly switch tabs and frames
            self._on_stage_changed("Analysis")
            # Also update the button appearance
            self.stage_tabs.set_active_stage("Analysis")
            logger.debug("Switched to Analysis via Ctrl+5")
        except Exception as e:
            logger.error(f"Error switching to Analysis: {e}")
        return "break"
    
    def _increase_font_size(self, event):
        """Handle Ctrl++ shortcut - increase font size."""
        try:
            current_scale = get_config_value('gui.font_scale', 1.2)
            
            # Find current scale in options and move to next
            scale_options = list(FONT_SCALE_OPTIONS.values())
            scale_options.sort()
            
            # Find current position
            current_index = -1
            for i, scale in enumerate(scale_options):
                if abs(scale - current_scale) < 0.01:  # Float comparison
                    current_index = i
                    break
            
            # Move to next scale if not at maximum
            if current_index < len(scale_options) - 1:
                new_scale = scale_options[current_index + 1]
                set_config_value('gui.font_scale', new_scale)
                save_config(self.config)
                self._update_fonts(new_scale)
                
                # Find the name for display
                scale_name = "Unknown"
                for name, value in FONT_SCALE_OPTIONS.items():
                    if abs(value - new_scale) < 0.01:
                        scale_name = name
                        break
                
                self._update_status(f"Font size increased to {scale_name}")
                logger.info(f"Font size increased to {new_scale} ({scale_name})")
            else:
                self._update_status("Font size is already at maximum")
                
        except Exception as e:
            logger.error(f"Error increasing font size: {e}")
            self._update_status("Error increasing font size")
        
        return "break"
    
    def _decrease_font_size(self, event):
        """Handle Ctrl+- shortcut - decrease font size."""
        try:
            current_scale = get_config_value('gui.font_scale', 1.2)
            
            # Find current scale in options and move to previous
            scale_options = list(FONT_SCALE_OPTIONS.values())
            scale_options.sort()
            
            # Find current position
            current_index = -1
            for i, scale in enumerate(scale_options):
                if abs(scale - current_scale) < 0.01:  # Float comparison
                    current_index = i
                    break
            
            # Move to previous scale if not at minimum
            if current_index > 0:
                new_scale = scale_options[current_index - 1]
                set_config_value('gui.font_scale', new_scale)
                save_config(self.config)
                self._update_fonts(new_scale)
                
                # Find the name for display
                scale_name = "Unknown"
                for name, value in FONT_SCALE_OPTIONS.items():
                    if abs(value - new_scale) < 0.01:
                        scale_name = name
                        break
                
                self._update_status(f"Font size decreased to {scale_name}")
                logger.info(f"Font size decreased to {new_scale} ({scale_name})")
            else:
                self._update_status("Font size is already at minimum")
                
        except Exception as e:
            logger.error(f"Error decreasing font size: {e}")
            self._update_status("Error decreasing font size")
        
        return "break"
    
    def _debug_keypress(self, event):
        """Debug method to log all keypress events for troubleshooting."""
        # Only log Control, Alt, and Function key combinations to avoid spam
        if event.state & 0x4 or event.state & 0x8 or event.keysym.startswith('F'):  # Control or Alt pressed or Function key
            key_info = {
                'keysym': event.keysym,
                'keycode': event.keycode,
                'char': repr(event.char),
                'state': event.state,
                'state_hex': hex(event.state)
            }
            logger.info(f"Key debug: {key_info}")
            
            # Also show in status for immediate feedback
            if event.state & 0x4:  # Control pressed
                self._update_status(f"Ctrl+{event.keysym} pressed (keycode: {event.keycode})")
            elif event.state & 0x8:  # Alt pressed
                self._update_status(f"Alt+{event.keysym} pressed (keycode: {event.keycode})")
            elif event.keysym.startswith('F'):
                self._update_status(f"{event.keysym} pressed")
                
        # Don't return "break" - let other handlers process the event
        return None
    
    def _handle_main_window_click(self, event):
        """Handle main window clicks - only set focus if not clicking on input widgets."""
        # Get the widget that was clicked
        clicked_widget = event.widget
        
        # Check if the clicked widget or its parents are input widgets
        widget_chain = []
        current = clicked_widget
        while current:
            widget_chain.append(current)
            if hasattr(current, 'master'):
                current = current.master
            else:
                break
        
        # Check if any widget in the chain is an input widget
        input_widget_types = ('CTkEntry', 'CTkTextbox', 'CTkComboBox', 'Entry', 'Text', 'Combobox')
        
        for widget in widget_chain:
            widget_class = widget.__class__.__name__
            if any(input_type in widget_class for input_type in input_widget_types):
                # Don't steal focus from input widgets
                return "break"
        
        # Only set focus to main window if we're not clicking on input areas
        self.focus_set()
        return None
    
    def _create_shortcut_test_dialog(self):
        """Create a dialog to test keyboard shortcuts."""
        test_window = ctk.CTkToplevel(self)
        test_window.title("Keyboard Shortcut Tester")
        test_window.geometry("500x400")
        test_window.transient(self)
        
        # Instructions
        instructions = ctk.CTkLabel(
            test_window,
            text="Test your keyboard shortcuts here.\nPress keys and see what's detected:",
            font=FONTS['body']
        )
        instructions.pack(pady=10)
        
        # Status display
        status_var = tk.StringVar(value="Ready - press a key combination...")
        status_label = ctk.CTkLabel(test_window, textvariable=status_var, font=FONTS['small'])
        status_label.pack(pady=5)
        
        # Key info display
        key_info_text = ctk.CTkTextbox(test_window, height=200, width=460)
        key_info_text.pack(pady=10, padx=20)
        key_info_text.insert("0.0", "Key detection log will appear here...\n")
        
        def test_keypress(event):
            key_info = f"""
Key pressed: {event.keysym}
Key code: {event.keycode}
Character: {repr(event.char)}
State: {event.state} (0x{event.state:x})
Control: {'Yes' if event.state & 0x4 else 'No'}
Alt: {'Yes' if event.state & 0x8 else 'No'}
Shift: {'Yes' if event.state & 0x1 else 'No'}
---
"""
            key_info_text.insert("end", key_info)
            key_info_text.see("end")
            
            # Update status
            modifiers = []
            if event.state & 0x4: modifiers.append("Ctrl")
            if event.state & 0x8: modifiers.append("Alt")
            if event.state & 0x1: modifiers.append("Shift")
            
            mod_str = "+".join(modifiers)
            if mod_str:
                status_var.set(f"Detected: {mod_str}+{event.keysym}")
            else:
                status_var.set(f"Detected: {event.keysym}")
        
        # Bind keypress to test window
        test_window.bind("<KeyPress>", test_keypress)
        test_window.focus_set()
        
        # Buttons for manual stage testing
        button_frame = ctk.CTkFrame(test_window)
        button_frame.pack(pady=10)
        
        stages = ["Visualize", "Preparation", "Builder", "Equilibration", "Analysis"]
        for i, stage in enumerate(stages):
            btn = ctk.CTkButton(
                button_frame,
                text=f"{i+1}. {stage}",
                command=lambda s=stage: self.stage_tabs.set_active_stage(s),
                width=80
            )
            btn.grid(row=0, column=i, padx=2)
        
        # Close button
        close_btn = ctk.CTkButton(test_window, text="Close", command=test_window.destroy)
        close_btn.pack(pady=10)
    
    def _show_help(self):
        """Show help dialog."""
        help_text = f"""
GateWizard v{__version__} - Help

Keyboard Shortcuts:
- Ctrl+O: Open PDB file (only in Visualize tab)

Stage Navigation (multiple options):
- Ctrl+1 or Alt+1 or F2: Switch to Visualize tab
- Ctrl+2 or Alt+2 or F3: Switch to Preparation tab
- Ctrl+3 or Alt+3 or F4: Switch to Builder tab
- Ctrl+4 or Alt+4 or F5: Switch to Equilibration tab
- Ctrl+5 or Alt+5 or F6: Switch to Analysis tab

Font Size:
- Ctrl++/= or Ctrl+Numpad+: Increase font size
- Ctrl+- or Ctrl+Numpad-: Decrease font size

Other:
- Ctrl+Q: Quit application
- F1: Show this help
- Ctrl+Shift+T: Open shortcut tester

Note: Stage shortcuts work with both main number row 
and numeric keypad. If Ctrl+numbers don't work on your 
keyboard, try Alt+numbers or F2-F6 instead.

Stages:
1. Visualize: Load and view protein structures
2. Preparation: Prepare protein structures (pKa, mutations, capping)
3. Builder: Build membrane systems
4. Equilibration: Setup MD protocols
5. Analysis: Analyze simulation results

For more information, visit the documentation.
        """
        
        # Create help window
        help_window = ctk.CTkToplevel(self)
        help_window.title("Gatewizard Help")
        help_window.geometry("500x400")
        help_window.resizable(False, False)
        
        # Center the help window
        help_window.transient(self)
        help_window.grab_set()
        
        # Create help text widget
        help_label = ctk.CTkLabel(
            help_window,
            text=help_text,
            font=FONTS['body'],
            justify="left"
        )
        help_label.pack(fill="both", expand=True, padx=20, pady=20)
        
        # Close button
        close_button = ctk.CTkButton(
            help_window,
            text="Close",
            command=help_window.destroy
        )
        close_button.pack(pady=10)
    
    def _show_settings(self):
        """Show settings dialog."""
        settings_window = ctk.CTkToplevel(self)
        settings_window.title("Gatewizard Settings")
        settings_window.geometry("500x550")
        settings_window.transient(self)
        
        # Center the settings window
        settings_window.update_idletasks()
        x = (settings_window.winfo_screenwidth() - 500) // 2
        y = (settings_window.winfo_screenheight() - 550) // 2
        settings_window.geometry(f"500x550+{x}+{y}")
        
        # Main scrollable frame
        scrollable_frame = ctk.CTkScrollableFrame(settings_window)
        scrollable_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        # Title
        title_label = ctk.CTkLabel(
            scrollable_frame,
            text="Settings",
            font=FONTS['heading']
        )
        title_label.pack(pady=(10, 20))
        
        # Font size section
        font_frame = ctk.CTkFrame(scrollable_frame)
        font_frame.pack(fill="x", pady=10)
        
        font_title = ctk.CTkLabel(
            font_frame,
            text="Display Settings",
            font=FONTS['heading']
        )
        font_title.pack(pady=(10, 5))
        
        font_label = ctk.CTkLabel(
            font_frame,
            text="Font Size:",
            font=FONTS['body']
        )
        font_label.pack(side="left", padx=10, pady=10)
        
        # Get current font scale (default to Large)
        current_scale = get_config_value('gui.font_scale', 1.2)
        current_option = 'Normal'
        for option, scale in FONT_SCALE_OPTIONS.items():
            if abs(scale - current_scale) < 0.01:
                current_option = option
                break
        
        font_var = ctk.StringVar(value=current_option)
        font_dropdown = ctk.CTkComboBox(
            font_frame,
            values=list(FONT_SCALE_OPTIONS.keys()),
            variable=font_var,
            width=120,
            command=self._on_font_size_changed
        )
        font_dropdown.pack(side="right", padx=10, pady=10)
        
        # Logging section
        logging_frame = ctk.CTkFrame(scrollable_frame)
        logging_frame.pack(fill="x", pady=10)
        
        logging_title = ctk.CTkLabel(
            logging_frame,
            text="Logging Settings",
            font=FONTS['heading']
        )
        logging_title.pack(pady=(10, 5))
        
        # Get current logging config
        config = load_config()
        
        # Enable logging checkbox
        logging_enabled_var = ctk.BooleanVar(value=config.logging.enable_file_logging)
        enable_logging_checkbox = ctk.CTkCheckBox(
            logging_frame,
            text="Enable file logging",
            variable=logging_enabled_var,
            font=FONTS['body']
        )
        enable_logging_checkbox.pack(anchor="w", padx=10, pady=(5, 0))
        
        # Log file path (auto-set to launch directory)
        path_frame = ctk.CTkFrame(logging_frame, fg_color="transparent")
        path_frame.pack(fill="x", padx=10, pady=5)
        
        path_label = ctk.CTkLabel(
            path_frame,
            text="Log folder (auto-set to app launch directory):",
            font=FONTS['body']
        )
        path_label.pack(anchor="w", pady=(5, 0))
        
        # Show the current working directory (where app was launched)
        current_log_path = self.initial_working_directory
        path_var = ctk.StringVar(value=current_log_path)
        path_entry = ctk.CTkEntry(
            path_frame,
            textvariable=path_var,
            font=FONTS['body'],
            width=300,
            state="readonly"  # Make it read-only to indicate it's auto-managed
        )
        path_entry.pack(side="left", fill="x", expand=True, pady=2)
        
        # Info label instead of browse button
        info_label = ctk.CTkLabel(
            path_frame,
            text="📁 Auto",
            font=FONTS['small'],
            width=80
        )
        info_label.pack(side="right", padx=(5, 0), pady=2)
        
        # Log file name
        name_frame = ctk.CTkFrame(logging_frame, fg_color="transparent")
        name_frame.pack(fill="x", padx=10, pady=5)
        
        name_label = ctk.CTkLabel(
            name_frame,
            text="Log file name:",
            font=FONTS['body']
        )
        name_label.pack(anchor="w", pady=(5, 0))
        
        name_var = ctk.StringVar(value=config.logging.log_file_name)
        name_entry = ctk.CTkEntry(
            name_frame,
            textvariable=name_var,
            font=FONTS['body']
        )
        name_entry.pack(fill="x", pady=2)
        
        # Log level
        level_frame = ctk.CTkFrame(logging_frame, fg_color="transparent")
        level_frame.pack(fill="x", padx=10, pady=5)
        
        level_label = ctk.CTkLabel(
            level_frame,
            text="Log level:",
            font=FONTS['body']
        )
        level_label.pack(side="left", pady=(5, 0))
        
        level_var = ctk.StringVar(value=config.logging.log_level)
        level_dropdown = ctk.CTkComboBox(
            level_frame,
            values=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
            variable=level_var,
            width=120,
            font=FONTS['body']
        )
        level_dropdown.pack(side="right", pady=2)
        
        # Max file size
        size_frame = ctk.CTkFrame(logging_frame, fg_color="transparent")
        size_frame.pack(fill="x", padx=10, pady=5)
        
        size_label = ctk.CTkLabel(
            size_frame,
            text="Max file size (MB):",
            font=FONTS['body']
        )
        size_label.pack(side="left", pady=(5, 0))
        
        size_var = ctk.StringVar(value=str(config.logging.max_file_size_mb))
        size_entry = ctk.CTkEntry(
            size_frame,
            textvariable=size_var,
            font=FONTS['body'],
            width=80
        )
        size_entry.pack(side="right", pady=2)
        
        # Backup count
        backup_frame = ctk.CTkFrame(logging_frame, fg_color="transparent")
        backup_frame.pack(fill="x", padx=10, pady=5)
        
        backup_label = ctk.CTkLabel(
            backup_frame,
            text="Backup files count:",
            font=FONTS['body']
        )
        backup_label.pack(side="left", pady=(5, 0))
        
        backup_var = ctk.StringVar(value=str(config.logging.backup_count))
        backup_entry = ctk.CTkEntry(
            backup_frame,
            textvariable=backup_var,
            font=FONTS['body'],
            width=80
        )
        backup_entry.pack(side="right", pady=(5, 10))
        
        # Buttons frame
        button_frame = ctk.CTkFrame(scrollable_frame, fg_color="transparent")
        button_frame.pack(fill="x", pady=(20, 10))
        
        # Apply button
        apply_button = ctk.CTkButton(
            button_frame,
            text="Apply",
            command=lambda: self._apply_settings(
                font_var.get(), 
                logging_enabled_var.get(),
                path_var.get(),
                name_var.get(),
                level_var.get(),
                size_var.get(),
                backup_var.get(),
                settings_window
            )
        )
        apply_button.pack(side="right", padx=(0, 10))
        
        # Close button
        close_button = ctk.CTkButton(
            button_frame,
            text="Close",
            command=settings_window.destroy
        )
        close_button.pack(side="right")
        
        # Wait for window to be fully rendered before grabbing
        def setup_grab():
            try:
                if settings_window.winfo_exists():
                    settings_window.grab_set()
                    settings_window.focus()
            except:
                pass  # Ignore grab errors
        
        # Use tracked callback for settings window
        settings_window.after(100, setup_grab)
    
    def _on_font_size_changed(self, selection):
        """Handle font size dropdown change (preview)."""
        pass  # Could add live preview here if desired
    
    # def _browse_log_folder(self, path_var):
    #     """Open folder browser for log file path."""
    #     # No longer used - log folder is auto-set to app launch directory
    #     import tkinter.filedialog as fd
    #     
    #     folder_path = fd.askdirectory(
    #         title="Select Log Folder",
    #         initialdir=path_var.get() or os.path.expanduser("~")
    #     )
    #     
    #     if folder_path:
    #         path_var.set(folder_path)
    
    def _apply_settings(self, font_size_option, logging_enabled, log_path, log_name, 
                       log_level, max_size, backup_count, settings_window):
        """Apply settings changes."""
        try:
            # Save font scale
            scale_factor = FONT_SCALE_OPTIONS[font_size_option]
            
            # Load config and update all settings
            config = load_config()
            config.gui.font_scale = scale_factor
            
            # Update logging settings (except log_file_path which is auto-managed)
            config.logging.enable_file_logging = logging_enabled
            # Note: log_file_path is auto-set to current working directory on each launch, don't save it
            config.logging.log_file_name = log_name.strip() or "gatewizard.log"
            config.logging.log_level = log_level
            
            # Validate and set numeric values
            try:
                config.logging.max_file_size_mb = max(1, int(max_size))
            except ValueError:
                config.logging.max_file_size_mb = 10
                
            try:
                config.logging.backup_count = max(0, int(backup_count))
            except ValueError:
                config.logging.backup_count = 5
            
            # Save all changes
            save_config(config)
            
            # Update all fonts in the application
            self._update_fonts(scale_factor)
            
            # Show confirmation
            status_msg = f"Settings applied. Font: {font_size_option}"
            if logging_enabled:
                log_file = os.path.join(self.initial_working_directory, log_name)
                status_msg += f", Logging: enabled ({log_file})"
            else:
                status_msg += ", Logging: disabled"
            
            self._show_status_message(status_msg)
            
            # Close settings window
            settings_window.destroy()
            
        except Exception as e:
            logger.error(f"Error applying settings: {e}")
            self._show_status_message("Error applying settings")
    
    def _update_fonts(self, scale_factor):
        """Update all fonts in the application."""
        try:
            # Get scaled fonts
            scaled_fonts = get_scaled_fonts(scale_factor)
            
            # Update main window fonts
            self._update_widget_font(self.title_label, scaled_fonts['title'])
            # App descriptor inline with title
            self._update_widget_font(self.app_descriptor_label, scaled_fonts['subtitle'])
            self._update_widget_font(self.settings_button, scaled_fonts['small'])
            self._update_widget_font(self.status_label, scaled_fonts['small'])
            
            # Update stage tab fonts
            if hasattr(self.stage_tabs, 'update_fonts'):
                self.stage_tabs.update_fonts(scaled_fonts)  # type: ignore
            
            # Update stage frame fonts
            for frame in self.stage_frames.values():
                if hasattr(frame, 'update_fonts'):
                    frame.update_fonts(scaled_fonts)
            
            logger.info(f"Fonts updated with scale factor: {scale_factor}")
            
        except Exception as e:
            logger.error(f"Error updating fonts: {e}")
    
    def _update_widget_font(self, widget, font):
        """Safely update a widget's font."""
        try:
            widget.configure(font=font)
        except Exception as e:
            logger.debug(f"Could not update font for widget {type(widget).__name__}: {e}")
    
    def _show_status_message(self, message, duration=3000):
        """Show a temporary status message."""
        if hasattr(self, 'status_label'):
            original_text = self.status_label.cget("text")
            self.status_label.configure(text=message)
            # Reset after duration with safety check
            def reset_status():
                try:
                    if self.winfo_exists() and hasattr(self, 'status_label') and self.status_label.winfo_exists():
                        self.status_label.configure(text=original_text)
                except:
                    pass
            self._schedule_callback(duration, reset_status)
    
    def _on_window_configure(self, event):
        """Handle window configuration changes."""
        if event.widget == self:
            try:
                # Check if window is in a valid state
                current_state = self.state()
                
                # Don't save position if window is maximized, iconified, or withdrawn
                if current_state not in ('normal', None):
                    return
                
                # Save window position if enabled (batched to avoid spam)
                if get_config_value('gui.remember_window_position', True):
                    # Cancel previous timer if exists
                    if self._position_save_timer:
                        self.after_cancel(self._position_save_timer)
                    
                    # Schedule save after 1 second of no movement
                    self._position_save_timer = self._schedule_callback(1000, self._save_window_position)
                    
            except Exception as e:
                # Log the error but don't crash the app
                logger.warning(f"Error in window configure handler: {e}")
    
    def _on_window_state_change(self, event):
        """Handle window state changes (maximize, minimize, restore)."""
        if event.widget == self:
            try:
                # Schedule a state check after a short delay to ensure state is stable
                self._schedule_callback(100, self._check_window_state)
            except Exception as e:
                logger.warning(f"Error in window state change handler: {e}")
    
    def _check_window_state(self):
        """Check and handle current window state safely."""
        try:
            current_state = self.state()
            
            # If window was maximized and then moved to another monitor,
            # we might need to handle geometry recovery
            if current_state == 'zoomed':
                # Window is maximized - this is fine, just don't save position
                return
            elif current_state == 'normal':
                # Window is in normal state - validate geometry is reasonable
                self._validate_window_geometry()
            
        except Exception as e:
            logger.warning(f"Error checking window state: {e}")
            # If we can't determine state, try to restore to a safe position
            self._restore_safe_geometry()
    
    def _validate_window_geometry(self):
        """Validate that the current window geometry is reasonable and visible."""
        try:
            # Get current window position and size
            x = self.winfo_x()
            y = self.winfo_y()
            width = self.winfo_width()
            height = self.winfo_height()
            
            # Get screen dimensions (this gets total virtual screen across all monitors)
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            
            # Check for the specific problematic coordinates we see in multi-monitor setups
            # Values like -32730 are clearly invalid window manager artifacts
            problematic_coords = [-32730, -32709, 32000, -32000]
            
            coords_invalid = False
            for coord in problematic_coords:
                if abs(x - coord) < 100 or abs(y - coord) < 100:
                    logger.warning(f"Detected problematic coordinate: x={x}, y={y}")
                    coords_invalid = True
                    break
            
            # Check if window is completely off-screen or in an invalid position
            margin = 50  # Allow some margin for window decorations
            
            # More strict validation for clearly invalid coordinates
            if (coords_invalid or
                x < -width + margin or x > screen_width - margin or
                y < -height + margin or y > screen_height - margin or
                width < 100 or height < 100 or
                x < -5000 or x > screen_width + 5000 or  # Clearly off-screen
                y < -5000 or y > screen_height + 5000):
                
                logger.warning(f"Window geometry appears invalid (x={x}, y={y}, w={width}, h={height}), restoring safe position")
                self._restore_safe_geometry()
                
        except Exception as e:
            logger.warning(f"Error validating window geometry: {e}")
            self._restore_safe_geometry()
    
    def _restore_safe_geometry(self):
        """Restore window to a safe, visible position and size."""
        try:
            # Use default safe dimensions
            safe_width = WINDOW_CONFIG['default_width']
            safe_height = WINDOW_CONFIG['default_height']
            
            # Position in center of primary screen
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            
            # For multi-monitor setups, use primary screen (left half)
            primary_width = screen_width // 2 if screen_width > 2000 else screen_width
            
            safe_x = (primary_width - safe_width) // 2
            safe_y = (screen_height - safe_height) // 2
            
            # Ensure positive coordinates
            safe_x = max(50, safe_x)
            safe_y = max(50, safe_y)
            
            # Apply safe geometry
            self.geometry(f"{safe_width}x{safe_height}+{safe_x}+{safe_y}")
            
            logger.info(f"Restored window to safe geometry: {safe_width}x{safe_height}+{safe_x}+{safe_y}")
            
        except Exception as e:
            logger.error(f"Failed to restore safe geometry: {e}")
            # Last resort - use minimal safe position
            try:
                self.geometry("800x600+100+100")
            except Exception:
                pass
    
    def _save_window_position(self):
        """Save window position (called by timer to batch saves)."""
        try:
            # Only save position if window is in normal state
            current_state = self.state()
            if current_state not in ('normal', None):
                self._position_save_timer = None
                return
            
            x = self.winfo_x()
            y = self.winfo_y()
            width = self.winfo_width()
            height = self.winfo_height()
            
            # Validate that the position is reasonable before saving
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            
            # Check for specific problematic coordinates that cause crashes
            problematic_coords = [-32730, -32709, 32000, -32000]
            coords_invalid = False
            for coord in problematic_coords:
                if abs(x - coord) < 100 or abs(y - coord) < 100:
                    coords_invalid = True
                    break
            
            # Don't save obviously invalid positions
            if (coords_invalid or
                x < -width or x > screen_width or 
                y < -height or y > screen_height or
                width < 100 or height < 100 or
                x < -5000 or x > screen_width + 5000 or  # Clearly off-screen  
                y < -5000 or y > screen_height + 5000):
                logger.warning(f"Skipping save of invalid window position: {x},{y} {width}x{height}")
                self._position_save_timer = None
                return
            
            # Load config, update values, save once
            config = load_config()
            config.gui.last_window_x = x
            config.gui.last_window_y = y
            config.gui.window_width = width
            config.gui.window_height = height
            save_config(config)
            
            logger.debug(f"Saved window position: {width}x{height}+{x}+{y}")
            self._position_save_timer = None
            
        except Exception as e:
            logger.warning(f"Error saving window position: {e}")
            self._position_save_timer = None
    
    def _setup_emergency_recovery(self):
        """Setup emergency recovery system for window state issues."""
        # Override the state method to add safety checks
        original_state = self.state
        original_geometry = self.geometry
        
        def safe_state(newstate=None):
            try:
                if newstate is None:
                    # Getting current state
                    return original_state()
                else:
                    # Setting new state
                    if newstate == 'zoomed':
                        # Before maximizing, ensure window size is reasonable
                        current_width = self.winfo_width()
                        current_height = self.winfo_height()
                        screen_width = self.winfo_screenwidth()
                        screen_height = self.winfo_screenheight()
                        
                        # If window is larger than screen, resize it first
                        if (current_width > screen_width - 50 or current_height > screen_height - 50):
                            logger.info(f"Window too large for maximize ({current_width}x{current_height} on {screen_width}x{screen_height}), resizing first")
                            safe_width = min(current_width, screen_width - 100)
                            safe_height = min(current_height, screen_height - 100)
                            original_geometry(f"{safe_width}x{safe_height}")
                            # Give time for resize to complete
                            self.update_idletasks()
                        
                        # Save current position before maximizing
                        self._save_window_position()
                        logger.info("Maximizing window...")
                    elif newstate == 'normal':
                        logger.info("Restoring window to normal state...")
                    
                    result = original_state(newstate)
                    
                    # After state change, validate position
                    if newstate == 'normal':
                        self._schedule_callback(100, self._validate_window_geometry)
                    
                    return result
                    
            except Exception as e:
                logger.error(f"Error in window state change: {e}")
                # Try to recover to a safe state
                try:
                    original_state('normal')
                    self._schedule_callback(200, self._restore_safe_geometry)
                except Exception:
                    pass
                raise
        
        def safe_geometry(newGeometry=None):
            try:
                if newGeometry is None:
                    # Getting current geometry
                    return original_geometry()
                else:
                    # Setting new geometry - validate it first
                    if isinstance(newGeometry, str):
                        # Parse and validate geometry string
                        if '+' in newGeometry or '-' in newGeometry:
                            # Has position component
                            parts = newGeometry.replace('+', ' +').replace('-', ' -').split()
                            if len(parts) >= 3:
                                try:
                                    x_part = parts[1] if parts[1].startswith(('+', '-')) else '+' + parts[1]
                                    y_part = parts[2] if parts[2].startswith(('+', '-')) else '+' + parts[2]
                                    x = int(x_part)
                                    y = int(y_part)
                                    
                                    # Validate coordinates are reasonable
                                    screen_width = self.winfo_screenwidth()
                                    screen_height = self.winfo_screenheight()
                                    
                                    # Check for clearly invalid coordinates (like -32730)
                                    if (x < -10000 or x > screen_width + 10000 or 
                                        y < -10000 or y > screen_height + 10000):
                                        logger.warning(f"Rejecting invalid geometry: {newGeometry}")
                                        logger.warning(f"Coordinates ({x},{y}) are outside reasonable bounds")
                                        return  # Don't apply invalid geometry
                                        
                                except (ValueError, IndexError):
                                    pass  # If parsing fails, let tkinter handle it
                    
                    return original_geometry(newGeometry)
                    
            except Exception as e:
                logger.warning(f"Error in geometry operation: {e}")
                raise
        
        # Replace the methods
        self.state = safe_state  # type: ignore
        self.geometry = safe_geometry  # type: ignore
        
        # Set up periodic health check
        self._setup_health_check()
    
    def _setup_health_check(self):
        """Setup periodic health check for window state."""
        def health_check():
            try:
                # Check if window is still responsive and in a valid state
                x = self.winfo_x()
                y = self.winfo_y()
                width = self.winfo_width()
                height = self.winfo_height()
                state = self.state()
                
                # Basic sanity checks
                if state == 'normal':
                    screen_width = self.winfo_screenwidth()
                    screen_height = self.winfo_screenheight()
                    
                    # If window is way off screen, bring it back
                    if (x < -width + 50 or x > screen_width - 50 or
                        y < -height + 50 or y > screen_height - 50):
                        logger.warning(f"Window appears off-screen ({x},{y}), correcting...")
                        self._restore_safe_geometry()
                
            except Exception as e:
                logger.warning(f"Health check failed: {e}")
            
            # Schedule next health check in 30 seconds
            self._schedule_callback(30000, health_check)
        
        # Start first health check after 10 seconds
        self._schedule_callback(10000, health_check)
    
    def _setup_safe_callback_system(self):
        """Setup safe callback system to prevent untracked callbacks."""
        # Store original after method
        self._original_after = self.after
        
        def safe_after(delay, callback=None, *args):
            """Safe after method that tracks callbacks or blocks during shutdown."""
            if getattr(self, '_shutting_down', False):
                # Don't schedule new callbacks during shutdown
                return None
            
            # If this is called without arguments, it's probably from a child widget
            # In that case, use the original method but track it
            if callback is None:
                # This is a tk.call('after', delay) type call - just return None for now
                return None
            
            # For normal callback scheduling, use our tracking system
            return self._schedule_callback(delay, callback, *args)
        
        # Replace the after method
        self.after = safe_after  # type: ignore
        
        # Also try to intercept Tk-level after calls if possible
        try:
            if hasattr(self, 'tk') and self.tk:
                original_tk_call = self.tk.call
                
                def safe_tk_call(command, *args):
                    """Intercept tk.call to catch direct 'after' calls."""
                    if command == 'after' and getattr(self, '_shutting_down', False):
                        # Block after calls during shutdown
                        return None
                    return original_tk_call(command, *args)
                
                self.tk.call = safe_tk_call
        except Exception as e:
            # If we can't intercept tk.call, that's okay - the global cleanup will handle it
            logger.debug(f"Could not intercept tk.call: {e}")
    
    def _check_initial_window_size(self):
        """Check if the initial window size is appropriate for the screen."""
        try:
            # Give the window time to be properly positioned
            self.update_idletasks()
            
            current_width = self.winfo_width()
            current_height = self.winfo_height()
            screen_width = self.winfo_screenwidth()
            screen_height = self.winfo_screenheight()
            
            # Check if window is too large for the screen
            margin = 80  # Space for taskbar and window decorations
            max_safe_width = screen_width - margin
            max_safe_height = screen_height - margin
            
            needs_resize = False
            new_width = current_width
            new_height = current_height
            
            if current_width > max_safe_width:
                new_width = max_safe_width
                needs_resize = True
                
            if current_height > max_safe_height:
                new_height = max_safe_height
                needs_resize = True
                
            if needs_resize:
                logger.info(f"Initial window size too large ({current_width}x{current_height}), adjusting to {new_width}x{new_height}")
                self.geometry(f"{new_width}x{new_height}")
                # Also update the saved config to prevent this issue in future
                try:
                    config = load_config()
                    config.gui.window_width = new_width
                    config.gui.window_height = new_height
                    save_config(config)
                except Exception as e:
                    logger.warning(f"Could not save adjusted window size: {e}")
                    
        except Exception as e:
            logger.warning(f"Error checking initial window size: {e}")
    
    def _handle_maximize_request(self):
        """Handle explicit maximize requests safely."""
        try:
            current_state = self.state()
            if current_state == 'zoomed':
                # Currently maximized, restore to normal
                self.state('normal')
                # Validate the restored position
                self._schedule_callback(100, self._validate_window_geometry)
            else:
                # Currently normal, check if window size is reasonable for maximizing
                current_width = self.winfo_width()
                current_height = self.winfo_height()
                screen_width = self.winfo_screenwidth()
                screen_height = self.winfo_screenheight()
                
                # If window is too large, resize it first
                if (current_width > screen_width - 50 or current_height > screen_height - 50):
                    logger.info(f"Resizing oversized window before maximize: {current_width}x{current_height} -> smaller")
                    safe_width = min(current_width, screen_width - 100)
                    safe_height = min(current_height, screen_height - 100)
                    self.geometry(f"{safe_width}x{safe_height}")
                    self.update_idletasks()
                
                # Save current position before maximizing
                self._save_window_position()
                self.state('zoomed')
        except Exception as e:
            logger.warning(f"Error handling maximize request: {e}")
            # If maximize fails, ensure we're in a good state
            try:
                self.state('normal')
                self._restore_safe_geometry()
            except Exception:
                pass
    
    def _on_closing(self):
        """
        Handle application closing.
        
        ⚠️ CRITICAL WARNING FOR FUTURE DEVELOPERS:
        DO NOT call self.destroy() in this method! CustomTkinter creates lambda
        callbacks during widget destruction that cause "invalid command name" errors.
        Only use self.quit() + sys.exit() for clean shutdown. See comments at end
        of this method for detailed explanation.
        """
        logger.info("Application closing...")
        
        # Set shutdown flag to prevent new callbacks
        self._shutting_down = True
        
        try:
            # Cancel any pending position saves
            if self._position_save_timer:
                self.after_cancel(self._position_save_timer)
            
            # Cancel all tracked scheduled callbacks
            if hasattr(self, '_scheduled_callbacks'):
                for callback_id in self._scheduled_callbacks:
                    try:
                        self.after_cancel(callback_id)
                    except Exception:
                        pass  # Callback may have already executed
                self._scheduled_callbacks.clear()
            
            # Save current window size and final position in one go - but only if in normal state
            try:
                current_state = self.state()
                if current_state in ('normal', None):
                    config = load_config()
                    config.gui.window_width = self.winfo_width()
                    config.gui.window_height = self.winfo_height()
                    config.gui.last_window_x = self.winfo_x()
                    config.gui.last_window_y = self.winfo_y()
                    save_config(config)
                else:
                    logger.info(f"Not saving window position - window state is {current_state}")
            except Exception as e:
                logger.warning(f"Could not save window position on close: {e}")
            
            # Cleanup frames
            for frame in self.stage_frames.values():
                if hasattr(frame, 'cleanup'):
                    try:
                        frame.cleanup()
                    except Exception as e:
                        logger.warning(f"Error cleaning up frame: {e}")
                
                # Also call cleanup_callbacks if available
                if hasattr(frame, 'cleanup_callbacks'):
                    try:
                        frame.cleanup_callbacks()
                    except Exception as e:
                        logger.warning(f"Error cleaning up callbacks in frame: {e}")
            
        except Exception as e:
            logger.error(f"Error during application cleanup: {e}")
        
        finally:
            # Cleanup to prevent callback errors
            try:
                # First, hide the window to prevent any visual issues
                self.withdraw()
            except Exception:
                pass
            
            # IMPORTANT: Give CustomTkinter time to finish any pending animations
            # CustomTkinter uses .after() for button animations and other effects
            # We need to let these complete before destroying the window
            try:
                self.update()  # Process any pending events
            except Exception:
                pass
            
            try:
                # Cancel ALL after callbacks globally - this will catch callbacks from child widgets too
                if hasattr(self, 'tk') and self.tk:
                    # Get all scheduled callbacks and cancel them
                    self.tk.call('after', 'info')  # This returns list of pending after callbacks
            except Exception:
                pass
            
            try:
                # Recursive cleanup of all child widgets
                self._cleanup_widget_callbacks(self)
            except Exception:
                pass
            
            try:
                # Force garbage collection to help clean up references
                import gc
                gc.collect()
            except Exception:
                pass
            
            try:
                # ============================================================
                # CRITICAL: Use quit() ONLY - DO NOT call destroy()!
                # ============================================================
                # CustomTkinter creates lambda callbacks during button click animations
                # (typically 50-100ms animations). When destroy() is called, it triggers
                # cleanup of widgets which in turn creates MORE lambda callbacks as part
                # of CustomTkinter's animation system cleanup process.
                # 
                # These late-stage lambda callbacks get scheduled in the event queue
                # but their command names become invalid when widgets are destroyed,
                # resulting in errors like: "invalid command name '138960674249088<lambda>'"
                #
                # SOLUTION: 
                # - Use quit() to stop the mainloop cleanly
                # - Use sys.exit(0) to terminate the process
                # - DO NOT call destroy() - let the OS handle window cleanup
                # - DO NOT add delays/sleeps - they don't help, destroy() is the problem
                #
                # This approach prevents the lambda callbacks from ever being created,
                # since destroy() never runs. The quit() method simply stops the event
                # loop, and sys.exit() terminates cleanly.
                # ============================================================
                self.quit()
            except Exception:
                pass
            
            # Exit cleanly - the OS will handle window cleanup automatically
            sys.exit(0)
    
    def _cleanup_widget_callbacks(self, widget):
        """Recursively cleanup callbacks from all child widgets."""
        try:
            # Check if this widget has any specific cleanup methods
            if hasattr(widget, 'cleanup_callbacks'):
                try:
                    widget.cleanup_callbacks()
                except Exception as e:
                    logger.debug(f"Error calling cleanup_callbacks on {type(widget).__name__}: {e}")
            
            # Check for common callback attributes and cancel them
            callback_attrs = [
                '_refresh_timer', 'refresh_after_id', 'auto_update_id', 
                'progress_timer', '_position_save_timer', '_after_id',
                'callback_id', 'timer_id', '_timer', '_callback_timer'
            ]
            
            for attr_name in callback_attrs:
                if hasattr(widget, attr_name):
                    callback_id = getattr(widget, attr_name)
                    if callback_id:
                        try:
                            widget.after_cancel(callback_id)
                            setattr(widget, attr_name, None)
                            logger.debug(f"Cancelled callback {attr_name} on {type(widget).__name__}")
                        except Exception as e:
                            logger.debug(f"Could not cancel {attr_name} on {type(widget).__name__}: {e}")
            
            # Look for callback ID lists
            callback_list_attrs = ['_callback_ids', '_scheduled_callbacks', '_timers']
            for attr_name in callback_list_attrs:
                if hasattr(widget, attr_name):
                    callback_list = getattr(widget, attr_name)
                    if isinstance(callback_list, list):
                        for callback_id in callback_list[:]:  # Copy to avoid modification during iteration
                            try:
                                widget.after_cancel(callback_id)
                                callback_list.remove(callback_id)
                            except Exception as e:
                                logger.debug(f"Could not cancel callback from {attr_name}: {e}")
            
            # Recursively process children
            if hasattr(widget, 'winfo_children'):
                for child in widget.winfo_children():
                    self._cleanup_widget_callbacks(child)
                    
        except Exception as e:
            logger.debug(f"Error cleaning up widget {type(widget).__name__}: {e}")

def main():
    """Main function to run the application."""
    try:
        app = ProteinViewerApp()
        app.mainloop()
    except KeyboardInterrupt:
        logger.info("Application interrupted by user")
    except Exception as e:
        logger.error(f"Unexpected error in main application: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()
