# gatewizard/core/file_manager.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
File management utilities for Gatewizard.

This module provides utilities for handling PDB files, working directories,
and other file operations commonly used in molecular dynamics workflows.
"""

import os
import shutil
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
import re

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class FileManagerError(Exception):
    """Custom exception for file management errors."""
    pass

class FileManager:
    """
    Utility class for file management operations.
    
    Provides methods for validating, copying, and manipulating files
    commonly used in molecular dynamics workflows.
    """
    
    @staticmethod
    def validate_pdb_file(file_path: str) -> Tuple[bool, str]:
        """
        Validate a PDB file.
        
        Args:
            file_path: Path to the PDB file
            
        Returns:
            Tuple of (is_valid, error_message)
        """
        if not file_path or not file_path.strip():
            return False, "No file path provided"
        
        file_path = file_path.strip()
        
        # Check if file exists
        if not os.path.exists(file_path):
            return False, f"File does not exist: {file_path}"
        
        # Check if it's a file (not directory)
        if not os.path.isfile(file_path):
            return False, f"Path is not a file: {file_path}"
        
        # Check file extension
        if not file_path.lower().endswith(('.pdb', '.ent')):
            return False, "File must have .pdb or .ent extension"
        
        # Check if file is readable
        try:
            with open(file_path, 'r') as f:
                # Read first few lines to check format
                lines = [f.readline().strip() for _ in range(10)]
        except PermissionError:
            return False, f"Permission denied reading file: {file_path}"
        except Exception as e:
            return False, f"Error reading file: {str(e)}"
        
        # Basic PDB format validation
        has_atom_records = any(
            line.startswith(('ATOM', 'HETATM')) for line in lines if line
        )
        
        if not has_atom_records:
            return False, "File does not contain ATOM or HETATM records"
        
        return True, ""
    
    @staticmethod
    def validate_working_directory(dir_path: str) -> Tuple[bool, str]:
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
    
    @staticmethod
    def clean_pdb_file(
        input_file: str, 
        output_file: Optional[str] = None,
        remove_waters: bool = True,
        remove_ligands: bool = False,
        remove_hetero: bool = False
    ) -> str:
        """
        Clean a PDB file by removing unwanted records.
        
        Args:
            input_file: Input PDB file path
            output_file: Output file path (if None, overwrites input)
            remove_waters: Remove water molecules
            remove_ligands: Remove ligand molecules
            remove_hetero: Remove all HETATM records
            
        Returns:
            Path to the cleaned file
            
        Raises:
            FileManagerError: If cleaning fails
        """
        if not os.path.exists(input_file):
            raise FileManagerError(f"Input file does not exist: {input_file}")
        
        output_file = output_file or input_file
        
        try:
            with open(input_file, 'r') as f:
                lines = f.readlines()
            
            cleaned_lines = []
            
            for line in lines:
                # Keep all non-coordinate records
                if not line.startswith(('ATOM', 'HETATM')):
                    cleaned_lines.append(line)
                    continue
                
                # Always keep ATOM records unless specifically removing hetero
                if line.startswith('ATOM'):
                    if not (remove_hetero and len(line) > 17 and line[17:20].strip() in ['HOH', 'WAT']):
                        cleaned_lines.append(line)
                    continue
                
                # Handle HETATM records
                if line.startswith('HETATM'):
                    if remove_hetero:
                        continue
                    
                    residue_name = line[17:20].strip() if len(line) > 20 else ""
                    
                    # Remove waters
                    if remove_waters and residue_name in ['HOH', 'WAT', 'TIP', 'SOL']:
                        continue
                    
                    # Remove common ligands (basic list)
                    if remove_ligands and residue_name in [
                        'ATP', 'ADP', 'GDP', 'GTP', 'NAD', 'FAD', 'HEM', 'CLR'
                    ]:
                        continue
                    
                    cleaned_lines.append(line)
            
            # Write cleaned file
            with open(output_file, 'w') as f:
                f.writelines(cleaned_lines)
            
            logger.info(f"PDB file cleaned: {input_file} -> {output_file}")
            return output_file
            
        except Exception as e:
            raise FileManagerError(f"Error cleaning PDB file: {str(e)}")
    
    @staticmethod
    def extract_pdb_info(file_path: str) -> Dict[str, Any]:
        """
        Extract basic information from a PDB file.
        
        Args:
            file_path: Path to the PDB file
            
        Returns:
            Dictionary containing PDB information
            
        Raises:
            FileManagerError: If extraction fails
        """
        if not os.path.exists(file_path):
            raise FileManagerError(f"PDB file does not exist: {file_path}")
        
        info = {
            'filename': os.path.basename(file_path),
            'file_size': os.path.getsize(file_path),
            'pdb_id': '',
            'title': '',
            'num_atoms': 0,
            'num_residues': 0,
            'chains': [],
            'has_waters': False,
            'has_ligands': False,
            'resolution': None
        }
        
        try:
            with open(file_path, 'r') as f:
                chains_seen = set()
                residues_seen = set()
                
                for line in f:
                    # Header information
                    if line.startswith('HEADER'):
                        info['pdb_id'] = line[62:66].strip()
                    elif line.startswith('TITLE'):
                        info['title'] += line[10:].strip() + ' '
                    elif line.startswith('REMARK   2') and 'RESOLUTION' in line:
                        # Extract resolution
                        match = re.search(r'(\d+\.?\d*)\s*ANGSTROM', line)
                        if match:
                            info['resolution'] = float(match.group(1))
                    
                    # Coordinate records
                    elif line.startswith(('ATOM', 'HETATM')):
                        info['num_atoms'] += 1
                        
                        if len(line) > 26:
                            chain_id = line[21:22].strip()
                            res_num = line[22:26].strip()
                            res_name = line[17:20].strip()
                            
                            chains_seen.add(chain_id)
                            residues_seen.add((chain_id, res_num, res_name))
                            
                            # Check for waters and ligands
                            if res_name in ['HOH', 'WAT', 'TIP', 'SOL']:
                                info['has_waters'] = True
                            elif line.startswith('HETATM') and res_name not in ['HOH', 'WAT', 'TIP', 'SOL']:
                                info['has_ligands'] = True
            
            info['title'] = info['title'].strip()
            info['chains'] = sorted(list(chains_seen))
            info['num_residues'] = len(residues_seen)
            
            return info
            
        except Exception as e:
            raise FileManagerError(f"Error extracting PDB info: {str(e)}")
    
    @staticmethod
    def backup_file(file_path: str, backup_suffix: str = '.bak') -> str:
        """
        Create a backup of a file.
        
        Args:
            file_path: Path to the file to backup
            backup_suffix: Suffix to add to backup file
            
        Returns:
            Path to the backup file
            
        Raises:
            FileManagerError: If backup fails
        """
        if not os.path.exists(file_path):
            raise FileManagerError(f"File does not exist: {file_path}")
        
        backup_path = file_path + backup_suffix
        
        try:
            shutil.copy2(file_path, backup_path)
            logger.info(f"File backed up: {file_path} -> {backup_path}")
            return backup_path
            
        except Exception as e:
            raise FileManagerError(f"Error creating backup: {str(e)}")
    
    @staticmethod
    def create_temp_file(suffix: str = '.tmp', prefix: str = 'gatewizard_') -> str:
        """
        Create a temporary file.
        
        Args:
            suffix: File suffix
            prefix: File prefix
            
        Returns:
            Path to the temporary file
        """
        fd, temp_path = tempfile.mkstemp(suffix=suffix, prefix=prefix)
        os.close(fd)  # Close file descriptor
        return temp_path
    
    @staticmethod
    def create_temp_directory(prefix: str = 'gatewizard_') -> str:
        """
        Create a temporary directory.
        
        Args:
            prefix: Directory prefix
            
        Returns:
            Path to the temporary directory
        """
        return tempfile.mkdtemp(prefix=prefix)
    
    @staticmethod
    def safe_remove(file_path: str) -> bool:
        """
        Safely remove a file or directory.
        
        Args:
            file_path: Path to remove
            
        Returns:
            True if successful, False otherwise
        """
        try:
            path_obj = Path(file_path)
            if path_obj.is_file():
                path_obj.unlink()
            elif path_obj.is_dir():
                shutil.rmtree(file_path)
            
            logger.debug(f"Removed: {file_path}")
            return True
            
        except Exception as e:
            logger.warning(f"Failed to remove {file_path}: {e}")
            return False
    
    @staticmethod
    def get_available_space(directory: str) -> int:
        """
        Get available disk space in bytes.
        
        Args:
            directory: Directory to check
            
        Returns:
            Available space in bytes
        """
        try:
            stat = shutil.disk_usage(directory)
            return stat.free
        except Exception:
            return 0
    
    @staticmethod
    def find_files(
        directory: str, 
        pattern: str = "*", 
        recursive: bool = False
    ) -> List[str]:
        """
        Find files matching a pattern.
        
        Args:
            directory: Directory to search
            pattern: File pattern (glob style)
            recursive: Search recursively
            
        Returns:
            List of matching file paths
        """
        try:
            dir_path = Path(directory)
            if recursive:
                files = list(dir_path.rglob(pattern))
            else:
                files = list(dir_path.glob(pattern))
            
            return [str(f) for f in files if f.is_file()]
            
        except Exception as e:
            logger.warning(f"Error finding files: {e}")
            return []