# gatewizard/tools/molecular_viewer.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Molecular visualization tools using matplotlib.

This module provides a high-level interface for molecular visualization
operations using matplotlib for 3D molecular rendering in the GUI.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
import requests
import os
from typing import Optional, List, Tuple, Dict, Any
from pathlib import Path

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class MolecularViewerError(Exception):
    """Custom exception for molecular viewer errors."""
    pass

class MolecularViewer:
    """
    Matplotlib-based molecular visualization.
    
    This class provides methods for loading structures, changing representations,
    applying color schemes, and rendering using matplotlib.
    """
    
    def __init__(self):
        """Initialize the molecular viewer."""
        self.current_structure = None
        self.atom_coords = None
        self.atom_elements = None
        self.backbone_coords = None
        self.residues = []
        self.secondary_structure = None
        self.chains = {}
        
        # Element colors
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
        
        logger.info("Matplotlib molecular viewer initialized successfully")
    
    def load_structure(self, file_path: str, object_name: str = "protein") -> bool:
        """Load a molecular structure from file."""
        try:
            self.parse_pdb_file(file_path)
            self.current_structure = object_name
            logger.info(f"Structure loaded: {file_path}")
            return True
        except Exception as e:
            logger.error(f"Error loading structure: {e}")
            return False
    
    def parse_pdb_file(self, file_path):
        """Parse PDB file and extract coordinates"""
        atoms = []
        elements = []
        backbone_atoms = []
        residues = []
        chains = {}
        current_residue = None
        
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        element = line[76:78].strip() or line[12:14].strip()[0]
                        atom_name = line[12:16].strip()
                        chain_id = line[21].strip()
                        res_id = int(line[22:26].strip())
                        res_name = line[17:20].strip()
                        
                        coord = [x, y, z]
                        atoms.append(coord)
                        elements.append(element)
                        
                        # Initialize chain if not exists
                        if chain_id not in chains:
                            chains[chain_id] = {'residues': [], 'ca_coords': []}
                        
                        # Handle residue grouping
                        if (current_residue is None or 
                            current_residue['id'] != res_id or 
                            current_residue['chain'] != chain_id):
                            
                            if current_residue is not None:
                                residues.append(current_residue)
                                chains[current_residue['chain']]['residues'].append(current_residue)
                            
                            current_residue = {
                                'name': res_name, 'id': res_id, 'chain': chain_id,
                                'atoms': [], 'ca_coord': None
                            }
                        
                        # Add atom to current residue
                        current_residue['atoms'].append({'coord': coord, 'element': element})
                        
                        if atom_name == 'CA':
                            current_residue['ca_coord'] = coord
                            chains[chain_id]['ca_coords'].append(coord)
                        
                        if atom_name in ['CA', 'N', 'C', 'O']:
                            backbone_atoms.append(coord)
                            
                    except (ValueError, IndexError):
                        continue
            
            # Add last residue
            if current_residue is not None:
                residues.append(current_residue)
                chains[current_residue['chain']]['residues'].append(current_residue)
        
        self.atom_coords = np.array(atoms)
        self.atom_elements = elements
        self.backbone_coords = np.array(backbone_atoms) if backbone_atoms else None
        self.residues = residues
        self.chains = chains
        
        # Calculate secondary structure
        self.calculate_secondary_structure()
    
    def calculate_secondary_structure(self):
        """Calculate secondary structure using simple geometric criteria"""
        self.secondary_structure = {}
        
        for chain_id, chain_data in self.chains.items():
            ca_coords = np.array(chain_data['ca_coords'])
            if len(ca_coords) < 4:
                continue
            
            ss_assignment = []
            for i in range(len(ca_coords)):
                if i < 2 or i >= len(ca_coords) - 2:
                    ss_assignment.append('C')
                else:
                    v1 = ca_coords[i-1] - ca_coords[i-2]
                    v2 = ca_coords[i] - ca_coords[i-1]
                    v3 = ca_coords[i+1] - ca_coords[i]
                    
                    angle1 = self.calculate_angle(v1, v2)
                    angle2 = self.calculate_angle(v2, v3)
                    
                    if 80 < angle1 < 120 and 80 < angle2 < 120:
                        ss_assignment.append('H')  # Helix
                    elif angle1 > 150 or angle2 > 150:
                        ss_assignment.append('E')  # Sheet
                    else:
                        ss_assignment.append('C')  # Coil
            
            self.secondary_structure[chain_id] = ss_assignment
    
    def calculate_angle(self, v1, v2):
        """Calculate angle between two vectors in degrees"""
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        return np.degrees(np.arccos(cos_angle))
    
    # Add other methods as needed for compatibility
    def change_representation(self, representation: str, selection: str = "all") -> bool:
        logger.debug(f"Representation change to {representation} (matplotlib viewer)")
        return True
    
    def change_color_scheme(self, scheme: str, selection: str = "all") -> bool:
        logger.debug(f"Color scheme change to {scheme} (matplotlib viewer)")
        return True
    
    def reset_view(self):
        logger.debug("View reset (matplotlib viewer)")
    
    def cleanup(self):
        logger.info("Matplotlib viewer cleaned up")