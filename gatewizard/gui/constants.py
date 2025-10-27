# gatewizard/gui/constants.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
GUI constants and configuration for Gatewizard.

This module contains all GUI-related constants including colors,
sizes, and other visual configuration parameters.
"""

from typing import Dict, List

# Color schemes using HEX codes
COLOR_SCHEME = {
    'background': '#1E1E1E',
    'canvas': '#2D2D2D', 
    'text': '#FFFFFF',
    'buttons': '#404040',
    'button_text': '#FFFFFF',
    'highlight': '#1f6aa5',
    'active': '#28a745',
    'hover': '#2580c9',
    'inactive': '#6c757d',
    'viewer_bg': '#000000',
    'content_bg': '#2D2D2D',  # content areas of the stages
    'content_inside_bg': '#393E41'  # lighter content for containers
}

# Button color constants
BUTTON_COLORS = {
    'active': "#1f6aa5",     # Dark blue for active button
    'inactive': "#404040",   # Gray for inactive button
    'hover_active': "#2580c9",  # Light blue for active hover
    'hover_inactive': "#4d4d4d"  # Light gray for inactive hover
}

# Window dimensions and constraints
WINDOW_CONFIG = {
    'default_width': 1100,
    'default_height': 700,
    'min_width': 1000,
    'min_height': 600,
}

# Font configurations
FONTS = {
    'title': ('Helvetica', 20, 'bold'),
    'subtitle': ('Helvetica', 16, 'bold'),
    'heading': ('Helvetica', 14, 'bold'),
    'body': ('Helvetica', 12),
    'small': ('Helvetica', 10),
    'code': ('Courier', 10),
}

# Font scaling options
FONT_SCALE_OPTIONS = {
    'Very Small': 0.7,
    'Small': 0.85,
    'Normal': 1.0,
    'Large': 1.2,
    'Very Large': 1.4,
    'Extra Large': 1.6
}

def get_scaled_fonts(scale_factor: float = 1.0) -> dict:
    """Get fonts scaled by the given factor."""
    scaled_fonts = {}
    for font_name, font_config in FONTS.items():
        name, size, *style = font_config
        scaled_size = max(8, int(size * scale_factor))  # Minimum size of 8
        if style:
            scaled_fonts[font_name] = (name, scaled_size, *style)
        else:
            scaled_fonts[font_name] = (name, scaled_size)
    return scaled_fonts

# Molecular viewer settings
VIEWER_CONFIG = {
    'default_representation': 'cartoon',
    'default_color_scheme': 'rainbow',
    'chain_colors': [
        'red', 'green', 'blue', 'yellow', 'cyan', 'magenta',
        'orange', 'pink', 'wheat', 'violet'
    ],
    'update_interval': 1/60,  # 60 FPS for smooth interaction
}

# Available representations for molecular viewer
REPRESENTATIONS = [
    "cartoon",
    "lines", 
    "sticks",
    "spheres",
    "surface",
    "mesh",
    "dots",
    "ribbon"
]

# Available color schemes for molecular viewer
COLOR_SCHEMES = [
    "rainbow",
    "spectrum", 
    "chainbows",
    "element",
    "hydrophobicity",
    "secondary_structure"
]

# Stage configuration for preparation workflow
PREPARATION_STAGES = [
    "Preparation",
    "Ligand Prep", 
    "System Build",
    "Equilibration"
]

# Progress tracker configuration
PROGRESS_CONFIG = {
    'refresh_interval': 1000,  # milliseconds
    'max_displayed_jobs': 50,
    'job_cleanup_interval': 300,  # seconds
}

# File dialog filters
FILE_FILTERS = {
    'pdb': [("PDB files", "*.pdb"), ("Entry files", "*.ent")],
    'image': [("PNG images", "*.png"), ("JPEG images", "*.jpg")],
    'all': [("All files", "*.*")]
}

# Default paths and directories
DEFAULT_PATHS = {
    'working_dir': '.',
    'temp_dir': 'temp',
    'output_dir': 'output',
    'logs_dir': 'logs'
}

# Validation limits
VALIDATION_LIMITS = {
    'max_filename_length': 255,
    'max_ph_value': 14.0,
    'min_ph_value': 0.0,
    'max_salt_concentration': 2.0,
    'min_salt_concentration': 0.0,
    'max_temperature': 500.0,  # Kelvin
    'min_temperature': 200.0,  # Kelvin
}

# Error messages
ERROR_MESSAGES = {
    'no_file_selected': "Please select a file",
    'file_not_found': "Selected file does not exist",
    'invalid_file_format': "Invalid file format",
    'permission_denied': "Permission denied accessing file",
    'invalid_ph': "pH must be between 0 and 14",
    'invalid_salt_conc': "Salt concentration must be between 0 and 2.0 M",
    'no_lipids_selected': "At least one lipid must be selected",
    'invalid_ratios': "Invalid lipid ratios format",
    'incompatible_ff': "Incompatible force field combination",
}

# Success messages  
SUCCESS_MESSAGES = {
    'file_loaded': "File loaded successfully",
    'preparation_started': "System preparation started",
    'analysis_complete': "Analysis completed successfully",
    'settings_saved': "Settings saved successfully",
}

# Layout spacing and padding
LAYOUT = {
    'padding_small': 5,
    'padding_medium': 10,
    'padding_large': 20,
    'spacing_small': 2,
    'spacing_medium': 5,
    'spacing_large': 10,
}

# Widget sizes
WIDGET_SIZES = {
    'button_width': 100,
    'button_height': 30,
    'entry_width': 300,
    'entry_height': 30,
    'combobox_width': 150,
    'combobox_height': 30,
    'progress_height': 200,
}

# Animation and transition settings
ANIMATION = {
    'fade_duration': 300,  # milliseconds
    'slide_duration': 250,
    'button_hover_delay': 100,
    'tooltip_delay': 500,
}

# Icon and image settings (if using icons)
ICONS = {
    'size_small': 16,
    'size_medium': 24,
    'size_large': 32,
    'format': 'png',
}

# Tooltips and help text
TOOLTIPS = {
    'pdb_file': "Select a PDB structure file for analysis",
    'working_dir': "Directory where output files will be saved",
    'ph_value': "pH value for protonation state calculations",
    'salt_conc': "Salt concentration in molar (M)",
    'force_fields': "Force field parameters for molecular dynamics",
    'lipid_ratios': "Relative amounts of lipids (format: upper//lower)",
}

# Help text for complex features
HELP_TEXT = {
    'lipid_ratios': """
Lipid ratios specify the relative amounts of different lipids in each leaflet.

Format: upper_leaflet_ratios//lower_leaflet_ratios

Examples:
- Single lipid per leaflet: 1//1
- Multiple lipids in upper: 0.7:0.3//1
- Multiple lipids in both: 0.5:0.3:0.2//0.6:0.4

Ratios will be normalized automatically.
    """,
    
    'force_fields': """
Force fields define the mathematical model for molecular interactions.

Recommended combinations:
- Membrane proteins: ff19SB + lipid21 + tip3p
- Soluble proteins: ff19SB + opc
- General purpose: ff14SB + tip3p

Always ensure compatibility between force field components.
    """,
    
    'propka': """
PROPKA predicts pKa values of ionizable residues in proteins.

The analysis considers:
- Local electrostatic environment
- Hydrogen bonding patterns
- Solvent accessibility
- Protein structure effects

Results help determine proper protonation states at target pH.
    """
}

# Status messages for different operations
STATUS_MESSAGES = {
    'idle': "Ready",
    'loading': "Loading...",
    'analyzing': "Analyzing structure...",
    'preparing': "Preparing system...",
    'validating': "Validating inputs...",
    'saving': "Saving results...",
    'complete': "Operation completed",
    'error': "Operation failed",
}

# Job monitoring constants
JOB_STATUS = {
    'running': "Running",
    'completed': "Completed", 
    'error': "Error",
    'unknown': "Unknown"
}

# System preparation steps
PREPARATION_STEPS = [
    "MEMEMBED",
    "Packmol", 
    "Parameterization",
    "Building Input"
]