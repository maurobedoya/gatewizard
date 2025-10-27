# gatewizard/tools/validators.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Validation utilities for system preparation inputs.

This module provides validation for all types of inputs
used in molecular dynamics system preparation.
"""

import os
import re
from typing import List, Tuple, Optional, Dict, Any, Union
from pathlib import Path

from gatewizard.utils.logger import get_logger
from gatewizard.tools.force_fields import ForceFieldManager

logger = get_logger(__name__)

class ValidationError(Exception):
    """Custom exception for validation errors."""
    pass

class SystemValidator:
    """
    Validator for system preparation inputs.
    
    This class provides validation methods for PDB files, force fields,
    lipid compositions, salt parameters, and other system preparation inputs.
    """
    
    def __init__(self):
        """Initialize the validator."""
        self.ff_manager = ForceFieldManager()
    
    def validate_pdb_file(self, file_path: str) -> Tuple[bool, str]:
        """
        Validate a PDB file.
        
        Args:
            file_path: Path to the PDB file
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not file_path or not file_path.strip():
            return False, "No PDB file path provided"
        
        file_path = file_path.strip()
        
        # Check if file exists
        if not os.path.exists(file_path):
            return False, f"PDB file does not exist: {file_path}"
        
        # Check if it's a file
        if not os.path.isfile(file_path):
            return False, f"Path is not a file: {file_path}"
        
        # Check file extension
        if not file_path.lower().endswith(('.pdb', '.ent')):
            return False, "File must have .pdb or .ent extension"
        
        # Check if file is readable and has basic PDB content
        try:
            with open(file_path, 'r') as f:
                lines = [f.readline().strip() for _ in range(20)]
            
            # Look for ATOM or HETATM records
            has_atoms = any(
                line.startswith(('ATOM', 'HETATM')) for line in lines if line
            )
            
            if not has_atoms:
                return False, "File does not contain ATOM or HETATM records"
            
            return True, ""
            
        except PermissionError:
            return False, f"Permission denied reading file: {file_path}"
        except Exception as e:
            return False, f"Error reading file: {str(e)}"
    
    def validate_working_directory(self, dir_path: str) -> Tuple[bool, str]:
        """
        Validate a working directory.
        
        Args:
            dir_path: Path to the directory
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not dir_path or not dir_path.strip():
            return False, "No directory path provided"
        
        dir_path = dir_path.strip()
        dir_obj = Path(dir_path)
        
        # Check if directory exists or can be created
        if dir_obj.exists():
            if not dir_obj.is_dir():
                return False, f"Path exists but is not a directory: {dir_path}"
            
            # Check if writable
            if not os.access(dir_path, os.W_OK):
                return False, f"Directory is not writable: {dir_path}"
        else:
            # Try to create directory
            try:
                dir_obj.mkdir(parents=True, exist_ok=True)
            except PermissionError:
                return False, f"Permission denied creating directory: {dir_path}"
            except Exception as e:
                return False, f"Error creating directory: {str(e)}"
        
        return True, ""
    
    def validate_lipid_composition(
        self, 
        upper_lipids: List[str], 
        lower_lipids: List[str]
    ) -> Tuple[bool, str]:
        """
        Validate lipid composition for both leaflets.
        
        Args:
            upper_lipids: List of lipids for upper leaflet
            lower_lipids: List of lipids for lower leaflet
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Check that at least one lipid is specified
        if not upper_lipids and not lower_lipids:
            return False, "At least one lipid must be selected"
        
        # Validate each lipid
        all_lipids = upper_lipids + lower_lipids
        available_lipids = set(self.ff_manager.get_available_lipids())
        
        for lipid in all_lipids:
            if not lipid or not lipid.strip():
                return False, "Empty lipid name found"
            
            if lipid not in available_lipids:
                return False, f"Unknown lipid: {lipid}"
        
        return True, ""
    
    def validate_lipid_ratios(
        self, 
        ratios_string: str, 
        upper_count: int, 
        lower_count: int
    ) -> Tuple[bool, str]:
        """
        Validate lipid ratios string.
        
        Args:
            ratios_string: Ratio string (format: upper_ratios//lower_ratios)
            upper_count: Number of upper leaflet lipids
            lower_count: Number of lower leaflet lipids
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # If only one lipid per leaflet, ratios are optional
        if upper_count <= 1 and lower_count <= 1 and not ratios_string.strip():
            return True, ""
        
        # If multiple lipids, ratios are required
        if (upper_count > 1 or lower_count > 1) and not ratios_string.strip():
            return False, (
                "Lipid ratios are required when using multiple lipids.\n"
                "Format: upper_ratios//lower_ratios (e.g., 0.7:0.3//1)"
            )
        
        if not ratios_string.strip():
            return True, ""
        
        try:
            # Split by leaflet separator
            leaflet_ratios = ratios_string.split("//")
            if len(leaflet_ratios) > 2:
                return False, "Invalid ratio format: too many leaflet separators (//)"
            
            # Validate upper leaflet ratios
            if upper_count > 0:
                upper_ratios_str = leaflet_ratios[0] if leaflet_ratios else ""
                valid, message = self._validate_ratio_numbers(
                    upper_ratios_str, upper_count, "upper"
                )
                if not valid:
                    return False, message
            
            # Validate lower leaflet ratios
            if lower_count > 0:
                if len(leaflet_ratios) < 2:
                    return False, "Missing ratios for lower leaflet (use // to separate upper and lower ratios)"
                
                lower_ratios_str = leaflet_ratios[1]
                valid, message = self._validate_ratio_numbers(
                    lower_ratios_str, lower_count, "lower"
                )
                if not valid:
                    return False, message
            
            return True, ""
            
        except Exception as e:
            return False, f"Error validating ratios: {str(e)}"
    
    def _validate_ratio_numbers(
        self, 
        ratios_str: str, 
        expected_count: int, 
        leaflet_name: str
    ) -> Tuple[bool, str]:
        """
        Validate individual ratio numbers.
        
        Args:
            ratios_str: String with ratios separated by ':'
            expected_count: Expected number of ratios
            leaflet_name: Name of the leaflet for error messages
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not ratios_str.strip():
            if expected_count > 0:
                return False, f"Missing ratios for {leaflet_name} leaflet"
            return True, ""
        
        try:
            ratios = [float(r.strip()) for r in ratios_str.split(":") if r.strip()]
            
            if len(ratios) != expected_count:
                return False, (
                    f"Number of {leaflet_name} leaflet ratios ({len(ratios)}) "
                    f"doesn't match number of lipids ({expected_count})"
                )
            
            if any(r <= 0 for r in ratios):
                return False, f"All ratios for {leaflet_name} leaflet must be positive numbers"
            
            return True, ""
            
        except ValueError:
            return False, f"Invalid ratio format in {leaflet_name} leaflet: must be numbers separated by ':'"
    
    def validate_force_fields(
        self, 
        water_model: str, 
        protein_ff: str, 
        lipid_ff: str
    ) -> Tuple[bool, str]:
        """
        Validate force field combination.
        
        Args:
            water_model: Water model name
            protein_ff: Protein force field name
            lipid_ff: Lipid force field name
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        return self.ff_manager.validate_combination(water_model, protein_ff, lipid_ff)
    
    def validate_salt_parameters(
        self, 
        concentration: str, 
        cation: str, 
        anion: str
    ) -> Tuple[bool, str]:
        """
        Validate salt parameters.
        
        Args:
            concentration: Salt concentration as string
            cation: Cation name
            anion: Anion name
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Validate concentration
        try:
            conc = float(concentration)
            if conc < 0:
                return False, "Salt concentration must be positive"
            if conc > 2.0:
                return False, "Salt concentration should not exceed 2.0 M"
        except ValueError:
            return False, "Invalid salt concentration: must be a number"
        
        # Validate cation
        valid_cation, cation_charge = self.ff_manager.validate_ion(cation)
        if not valid_cation:
            return False, f"Unknown cation: {cation}"
        
        if cation_charge <= 0:
            return False, f"Invalid cation charge: {cation} has charge {cation_charge}"
        
        # Validate anion
        valid_anion, anion_charge = self.ff_manager.validate_ion(anion)
        if not valid_anion:
            return False, f"Unknown anion: {anion}"
        
        if anion_charge >= 0:
            return False, f"Invalid anion charge: {anion} has charge {anion_charge}"
        
        return True, ""
    
    def validate_ph_value(self, ph_str: str) -> Tuple[bool, str]:
        """
        Validate pH value.
        
        Args:
            ph_str: pH value as string
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            ph = float(ph_str)
            if ph < 0 or ph > 14:
                return False, "pH must be between 0 and 14"
            return True, ""
        except ValueError:
            return False, "Invalid pH value: must be a number"
    
    def validate_system_inputs(
        self,
        pdb_file: str,
        upper_lipids: List[str],
        lower_lipids: List[str],
        lipid_ratios: str = "",
        **kwargs
    ) -> Tuple[bool, str]:
        """
        Validation of all system preparation inputs.
        
        Args:
            pdb_file: Path to PDB file
            upper_lipids: List of upper leaflet lipids
            lower_lipids: List of lower leaflet lipids
            lipid_ratios: Lipid ratios string
            **kwargs: Additional parameters to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Validate PDB file
        valid, error = self.validate_pdb_file(pdb_file)
        if not valid:
            return False, error
        
        # Validate lipid composition
        valid, error = self.validate_lipid_composition(upper_lipids, lower_lipids)
        if not valid:
            return False, error
        
        # Validate lipid ratios
        upper_count = len(upper_lipids)
        lower_count = len(lower_lipids)
        valid, error = self.validate_lipid_ratios(lipid_ratios, upper_count, lower_count)
        if not valid:
            return False, error
        
        # Validate force fields if provided
        water_model = kwargs.get('water_model')
        protein_ff = kwargs.get('protein_ff')
        lipid_ff = kwargs.get('lipid_ff')
        
        if water_model and protein_ff and lipid_ff:
            valid, error = self.validate_force_fields(water_model, protein_ff, lipid_ff)
            if not valid:
                return False, error
        
        # Validate salt parameters if salt is enabled
        if kwargs.get('add_salt', False):
            salt_conc = kwargs.get('salt_concentration', '0.15')
            cation = kwargs.get('cation', 'K+')
            anion = kwargs.get('anion', 'Cl-')
            
            valid, error = self.validate_salt_parameters(
                str(salt_conc), cation, anion
            )
            if not valid:
                return False, error
        
        # Validate pH if PDB2PQR is enabled
        if kwargs.get('use_pdb2pqr', False):
            ph = kwargs.get('ph', '7.0')
            valid, error = self.validate_ph_value(str(ph))
            if not valid:
                return False, error
        
        return True, ""
    
    def validate_propka_inputs(
        self, 
        pdb_file: str, 
        ph: str, 
        version: str = "3"
    ) -> Tuple[bool, str]:
        """
        Validate inputs for Propka analysis.
        
        Args:
            pdb_file: Path to PDB file
            ph: pH value as string
            version: Propka version
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        # Validate PDB file
        valid, error = self.validate_pdb_file(pdb_file)
        if not valid:
            return False, error
        
        # Validate pH
        valid, error = self.validate_ph_value(ph)
        if not valid:
            return False, error
        
        # Validate version
        if version not in ["3", "31"]:
            return False, f"Unsupported Propka version: {version}"
        
        return True, ""
    
    def validate_email_address(self, email: str) -> Tuple[bool, str]:
        """
        Validate an email address.
        
        Args:
            email: Email address to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not email or not email.strip():
            return False, "Email address is required"
        
        email = email.strip()
        
        # Basic email validation regex
        pattern = r'^[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}$'
        
        if not re.match(pattern, email):
            return False, "Invalid email address format"
        
        if len(email) > 254:  # RFC 5321 limit
            return False, "Email address is too long"
        
        return True, ""
    
    def validate_filename(self, filename: str) -> Tuple[bool, str]:
        """
        Validate a filename for illegal characters.
        
        Args:
            filename: Filename to validate
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not filename or not filename.strip():
            return False, "Filename cannot be empty"
        
        filename = filename.strip()
        
        # Check for illegal characters
        illegal_chars = r'<>:"/\\|?*'
        for char in illegal_chars:
            if char in filename:
                return False, f"Filename contains illegal character: {char}"
        
        # Check for control characters
        if any(ord(c) < 32 for c in filename):
            return False, "Filename contains control characters"
        
        # Check for reserved names (Windows)
        reserved_names = [
            'CON', 'PRN', 'AUX', 'NUL', 'COM1', 'COM2', 'COM3', 'COM4',
            'COM5', 'COM6', 'COM7', 'COM8', 'COM9', 'LPT1', 'LPT2',
            'LPT3', 'LPT4', 'LPT5', 'LPT6', 'LPT7', 'LPT8', 'LPT9'
        ]
        
        name_without_ext = os.path.splitext(filename)[0].upper()
        if name_without_ext in reserved_names:
            return False, f"Filename uses reserved name: {filename}"
        
        # Check length
        if len(filename) > 255:
            return False, "Filename is too long (max 255 characters)"
        
        return True, ""
    
    def validate_positive_number(
        self, 
        value_str: str, 
        name: str = "value", 
        max_value: Optional[float] = None
    ) -> Tuple[bool, str]:
        """
        Validate a positive number.
        
        Args:
            value_str: String representation of the number
            name: Name of the value for error messages
            max_value: Optional maximum value
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            value = float(value_str)
            
            if value <= 0:
                return False, f"{name} must be positive"
            
            if max_value is not None and value > max_value:
                return False, f"{name} must not exceed {max_value}"
            
            return True, ""
            
        except ValueError:
            return False, f"Invalid {name}: must be a number"
    
    def validate_temperature(self, temp_str: str) -> Tuple[bool, str]:
        """
        Validate temperature value.
        
        Args:
            temp_str: Temperature value as string
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            temp = float(temp_str)
            if temp < 200:
                return False, "Temperature must be at least 200 K"
            if temp > 500:
                return False, "Temperature should not exceed 500 K"
            return True, ""
        except ValueError:
            return False, "Invalid temperature: must be a number"
    
    def validate_pressure(self, pressure_str: str) -> Tuple[bool, str]:
        """
        Validate pressure value.
        
        Args:
            pressure_str: Pressure value as string
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            pressure = float(pressure_str)
            if pressure <= 0:
                return False, "Pressure must be positive"
            if pressure > 10:
                return False, "Pressure should not exceed 10 bar"
            return True, ""
        except ValueError:
            return False, "Invalid pressure: must be a number"
    
    def validate_simulation_time(self, time_str: str) -> Tuple[bool, str]:
        """
        Validate simulation time value.
        
        Args:
            time_str: Simulation time as string
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            time_val = float(time_str)
            if time_val <= 0:
                return False, "Simulation time must be positive"
            if time_val > 1000:
                return False, "Simulation time should not exceed 1000 ns"
            return True, ""
        except ValueError:
            return False, "Invalid simulation time: must be a number"
    
    def validate_timestep(self, timestep_str: str) -> Tuple[bool, str]:
        """
        Validate timestep value.
        
        Args:
            timestep_str: Timestep value as string
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            timestep = float(timestep_str)
            if timestep <= 0:
                return False, "Timestep must be positive"
            if timestep > 4:
                return False, "Timestep should not exceed 4 fs"
            return True, ""
        except ValueError:
            return False, "Invalid timestep: must be a number"
    
    def validate_membrane_dimensions(
        self, 
        x_str: str, 
        y_str: str, 
        z_str: Optional[str] = None
    ) -> Tuple[bool, str]:
        """
        Validate membrane/box dimensions.
        
        Args:
            x_str: X dimension as string
            y_str: Y dimension as string
            z_str: Z dimension as string (optional)
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        try:
            x = float(x_str)
            y = float(y_str)
            
            if x <= 0 or y <= 0:
                return False, "Membrane dimensions must be positive"
            
            if x < 50 or y < 50:
                return False, "Membrane dimensions should be at least 50 Å"
            
            if x > 500 or y > 500:
                return False, "Membrane dimensions should not exceed 500 Å"
            
            if z_str is not None:
                z = float(z_str)
                if z <= 0:
                    return False, "Z dimension must be positive"
                if z > 500:
                    return False, "Z dimension should not exceed 500 Å"
            
            return True, ""
            
        except ValueError:
            return False, "Invalid membrane dimensions: must be numbers"