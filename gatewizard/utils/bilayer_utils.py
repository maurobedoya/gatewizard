#!/usr/bin/env python3

#!/usr/bin/env python3
"""
Bilayer analysis utilities for gatewizard.

This module provides functions to analyze membrane systems and generate
colvar configurations for bilayer thickness measurements in ABF simulations.
"""

import os
import re
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass

@dataclass
class PhosphorusAtom:
    """Data class for phosphorus atom information."""
    atom_number: int        # PDB atom number (1-based)
    namd_index: int        # NAMD atom index (0-based) 
    atom_name: str         # Atom name (e.g., P31)
    residue_name: str      # Residue name (e.g., PC)
    residue_number: int    # Residue number
    chain_id: str          # Chain identifier
    x: float              # X coordinate
    y: float              # Y coordinate  
    z: float              # Z coordinate

class BilayerAnalyzer:
    """
    Utility class for analyzing bilayer systems and generating colvar configurations.
    """
    
    # Default phosphorus atom patterns (extensible)
    DEFAULT_PHOSPHORUS_PATTERNS = {
        'P31': 'Standard phosphatidylcholine phosphorus',
        'P': 'Generic phosphorus',
        'P1': 'Alternative phosphorus naming', 
        'P2': 'Alternative phosphorus naming',
        'P3': 'Alternative phosphorus naming',
        'P4': 'Alternative phosphorus naming',
        'PO4': 'Phosphate group',
    }
    
    # Default lipid residue patterns (extensible)
    DEFAULT_LIPID_PATTERNS = {
        'PC': 'Phosphatidylcholine',
        'PE': 'Phosphatidylethanolamine',
        'PS': 'Phosphatidylserine', 
        'PA': 'Phosphatidic acid',
        'PG': 'Phosphatidylglycerol',
        'PI': 'Phosphatidylinositol',
        'PIP': 'Phosphatidylinositol phosphate',
        'PIP2': 'Phosphatidylinositol bisphosphate',
        'POPC': 'Palmitoyloleoylphosphatidylcholine',
        'POPE': 'Palmitoyloleoylphosphatidylethanolamine',
        'POPS': 'Palmitoyloleoylphosphatidylserine',
        'DPPC': 'Dipalmitoylphosphatidylcholine',
        'DMPC': 'Dimyristoylphosphatidylcholine',
        'CHOL': 'Cholesterol',
        'CARD': 'Cardiolipin',
        'SM': 'Sphingomyelin',
    }
    
    def __init__(self, phosphorus_patterns: Dict[str, str] = None, 
                 lipid_patterns: Dict[str, str] = None):
        """
        Initialize bilayer analyzer with custom patterns.
        
        Args:
            phosphorus_patterns: Custom phosphorus atom name patterns
            lipid_patterns: Custom lipid residue name patterns
        """
        self.phosphorus_patterns = phosphorus_patterns or self.DEFAULT_PHOSPHORUS_PATTERNS.copy()
        self.lipid_patterns = lipid_patterns or self.DEFAULT_LIPID_PATTERNS.copy()
    
    def add_phosphorus_pattern(self, pattern: str, description: str = ""):
        """Add a new phosphorus atom pattern."""
        self.phosphorus_patterns[pattern] = description
    
    def add_lipid_pattern(self, pattern: str, description: str = ""):
        """Add a new lipid residue pattern."""
        self.lipid_patterns[pattern] = description
    
    def is_phosphorus_atom(self, atom_name: str) -> bool:
        """Check if atom name matches phosphorus patterns."""
        # Exact match
        if atom_name in self.phosphorus_patterns:
            return True
        
        # Pattern match: starts with P and reasonable length
        if atom_name.startswith('P') and len(atom_name) <= 4:
            return True
            
        return False
    
    def is_lipid_residue(self, residue_name: str) -> bool:
        """Check if residue name matches lipid patterns."""
        return residue_name in self.lipid_patterns
    
    def parse_pdb_atom_line(self, line: str) -> Optional[PhosphorusAtom]:
        """
        Parse PDB ATOM/HETATM line and return PhosphorusAtom if it's a lipid phosphorus.
        
        Args:
            line: PDB line string
            
        Returns:
            PhosphorusAtom object if line contains a lipid phosphorus atom, None otherwise
        """
        if not (line.startswith('ATOM') or line.startswith('HETATM')):
            return None
        
        try:
            atom_number = int(line[6:11].strip())
            atom_name = line[12:16].strip()
            residue_name = line[17:20].strip()
            chain_id = line[21:22].strip()
            residue_number = int(line[22:26].strip())
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            
            # Check if this is a phosphorus atom in a lipid
            if self.is_phosphorus_atom(atom_name) and self.is_lipid_residue(residue_name):
                return PhosphorusAtom(
                    atom_number=atom_number,
                    namd_index=atom_number - 1,  # Convert to 0-based for NAMD
                    atom_name=atom_name,
                    residue_name=residue_name,
                    residue_number=residue_number,
                    chain_id=chain_id,
                    x=x, y=y, z=z
                )
        except (ValueError, IndexError):
            pass
        
        return None
    
    def analyze_bilayer_from_pdb(self, pdb_file: str, z_threshold: float = 0.0) -> Tuple[List[PhosphorusAtom], List[PhosphorusAtom]]:
        """
        Analyze PDB file to identify phosphorus atoms in upper and lower bilayers.
        
        Args:
            pdb_file: Path to PDB file
            z_threshold: Z coordinate threshold to separate bilayers (default: 0.0)
            
        Returns:
            Tuple of (upper_bilayer_atoms, lower_bilayer_atoms)
        """
        upper_bilayer = []
        lower_bilayer = []
        
        if not os.path.exists(pdb_file):
            raise FileNotFoundError(f"PDB file not found: {pdb_file}")
        
        with open(pdb_file, 'r') as f:
            for line in f:
                phosphorus_atom = self.parse_pdb_atom_line(line)
                if phosphorus_atom:
                    if phosphorus_atom.z > z_threshold:
                        upper_bilayer.append(phosphorus_atom)
                    else:
                        lower_bilayer.append(phosphorus_atom)
        
        # Sort by atom number for consistent ordering
        upper_bilayer.sort(key=lambda x: x.atom_number)
        lower_bilayer.sort(key=lambda x: x.atom_number)
        
        return upper_bilayer, lower_bilayer
    
    def generate_bilayer_thickness_colvar(self, upper_bilayer: List[PhosphorusAtom], 
                                        lower_bilayer: List[PhosphorusAtom],
                                        colvar_name: str = "bilayer_thickness") -> str:
        """
        Generate NAMD colvar configuration for bilayer thickness measurement.
        
        Args:
            upper_bilayer: List of upper bilayer phosphorus atoms
            lower_bilayer: List of lower bilayer phosphorus atoms
            colvar_name: Name for the colvar (default: "bilayer_thickness")
            
        Returns:
            Colvar configuration string
        """
        if not upper_bilayer or not lower_bilayer:
            raise ValueError("Both upper and lower bilayers must contain phosphorus atoms")
        
        # Extract NAMD indices (0-based)
        upper_indices = [str(atom.namd_index) for atom in upper_bilayer]
        lower_indices = [str(atom.namd_index) for atom in lower_bilayer]
        
        # Format atom numbers in a single line
        def format_indices(indices: List[str], indent: str = " ") -> str:
            return indent + " ".join(indices)
        
        upper_formatted = format_indices(upper_indices)
        lower_formatted = format_indices(lower_indices)
        
        config = f"""colvar {{
    name {colvar_name}
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNumbers {upper_formatted}
        }}
        group2 {{
            atomNumbers {lower_formatted}
        }}
    }}
}}"""
        
        return config
    
    def get_bilayer_statistics(self, upper_bilayer: List[PhosphorusAtom], 
                             lower_bilayer: List[PhosphorusAtom]) -> Dict[str, any]:
        """
        Get statistical information about the bilayer system.
        
        Args:
            upper_bilayer: List of upper bilayer phosphorus atoms
            lower_bilayer: List of lower bilayer phosphorus atoms
            
        Returns:
            Dictionary with bilayer statistics
        """
        all_atoms = upper_bilayer + lower_bilayer
        
        stats = {
            'total_phosphorus_atoms': len(all_atoms),
            'upper_bilayer_count': len(upper_bilayer),
            'lower_bilayer_count': len(lower_bilayer),
            'residue_types': sorted(set(atom.residue_name for atom in all_atoms)),
            'atom_types': sorted(set(atom.atom_name for atom in all_atoms)),
        }
        
        if upper_bilayer:
            z_vals = [atom.z for atom in upper_bilayer]
            stats['upper_z_range'] = (min(z_vals), max(z_vals))
        
        if lower_bilayer:
            z_vals = [atom.z for atom in lower_bilayer]
            stats['lower_z_range'] = (min(z_vals), max(z_vals))
        
        return stats

def generate_bilayer_thickness_colvar_from_pdb(pdb_file: str, 
                                             output_file: str = None,
                                             colvar_name: str = "bilayer_thickness",
                                             z_threshold: float = 0.0,
                                             custom_phosphorus_patterns: Dict[str, str] = None,
                                             custom_lipid_patterns: Dict[str, str] = None) -> str:
    """
    Convenience function to generate bilayer thickness colvar from PDB file.
    
    Args:
        pdb_file: Path to input PDB file
        output_file: Optional output file path  
        colvar_name: Name for the colvar
        z_threshold: Z coordinate threshold to separate bilayers
        custom_phosphorus_patterns: Custom phosphorus atom patterns
        custom_lipid_patterns: Custom lipid residue patterns
        
    Returns:
        Colvar configuration string
    """
    analyzer = BilayerAnalyzer(custom_phosphorus_patterns, custom_lipid_patterns)
    upper_bilayer, lower_bilayer = analyzer.analyze_bilayer_from_pdb(pdb_file, z_threshold)
    
    if not upper_bilayer or not lower_bilayer:
        raise ValueError(f"Could not find phosphorus atoms in both bilayers in {pdb_file}")
    
    config = analyzer.generate_bilayer_thickness_colvar(upper_bilayer, lower_bilayer, colvar_name)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(config)
    
    return config

# Example usage
if __name__ == "__main__":
    # Example of how to use the module
    pdb_file = "system.pdb"
    
    # Create analyzer with default patterns
    analyzer = BilayerAnalyzer()
    
    # Add custom patterns if needed
    analyzer.add_phosphorus_pattern("P5", "Custom phosphorus variant")
    analyzer.add_lipid_pattern("CUSTOM", "Custom lipid type")
    
    try:
        # Analyze bilayer
        upper, lower = analyzer.analyze_bilayer_from_pdb(pdb_file)
        
        # Get statistics
        stats = analyzer.get_bilayer_statistics(upper, lower)
        print(f"Found {stats['total_phosphorus_atoms']} phosphorus atoms")
        print(f"Upper bilayer: {stats['upper_bilayer_count']}")
        print(f"Lower bilayer: {stats['lower_bilayer_count']}")
        
        # Generate colvar
        config = analyzer.generate_bilayer_thickness_colvar(upper, lower)
        
        # Save to file
        with open("bilayer_thickness.col", "w") as f:
            f.write(config)
            
    except FileNotFoundError:
        print(f"PDB file {pdb_file} not found")
    except ValueError as e:
        print(f"Error: {e}")
