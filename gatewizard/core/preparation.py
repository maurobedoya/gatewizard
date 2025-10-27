# gatewizard/core/preparation.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Protein preparation module using PROPKA for predicting pKa values and adjusting protonation states.

This module provides functionality to run PROPKA analysis on protein structures
and modify PDB files based on predicted pKa values and desired pH conditions.
"""

import subprocess
import os
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

# Protonation state mappings for different residue types (AMBER compatible)
PROTONATION_STATES = {
    "ASP": {"protonated": "ASH", "deprotonated": "ASP"},  # Aspartic acid
    "GLU": {"protonated": "GLH", "deprotonated": "GLU"},  # Glutamic acid
    "HIS": {
        "neutral_epsilon": "HIE",  # Histidine with proton on epsilon nitrogen (AMBER default)
        "neutral_delta": "HID",    # Histidine with proton on delta nitrogen
        "protonated": "HIP"        # Histidine with both nitrogens protonated
    },
    "LYS": {"neutral": "LYN", "protonated": "LYS"},       # Lysine
    "TYR": {"neutral": "TYR", "deprotonated": "TYM"},     # Tyrosine
    "CYS": {
        "free": "CYS",             # Free cysteine with SH side chain
        "deprotonated": "CYM",     # Deprotonated cysteine
        "disulfide": "CYX"         # Cysteine in disulfide bond (no H on sulfur)
    },
    "ARG": {"protonated": "ARG"},                         # Arginine (no neutral form in AMBER)
    # Terminal caps
    "NTE": {"capped": "ACE", "protonated": "NHE"},        # N-terminus caps
    "CTE": {"capped": "NME", "deprotonated": "COO"}       # C-terminus caps
}

class PreparationError(Exception):
    """Custom exception for preparation-related errors."""
    pass

class PreparationManager:
    """
    Class for managing protein preparation, PROPKA analysis and protonation state predictions.
    
    This class provides methods for running PROPKA, parsing results,
    and applying protonation states to PDB files.
    """
    
    def __init__(self, propka_version: str = "3"):
        """
        Initialize the preparation manager.
        
        Args:
            propka_version: Version of Propka to use (default: "3")
        """
        self.propka_version = propka_version
        self.last_analysis_file = None
        self.last_summary_file = None
        self.protonable_residues = []
    
    def run_analysis(self, pdb_file: str, output_dir: Optional[str] = None) -> str:
        """
        Run Propka analysis on a PDB file.
        
        Args:
            pdb_file: Path to the input PDB file
            output_dir: Optional directory to save the .pka file. If None, saves in same directory as input file.
            
        Returns:
            Path to the generated .pka file
            
        Raises:
            PreparationError: If Propka execution fails
            FileNotFoundError: If input file doesn't exist
        """
        # Verify input file exists and convert to absolute path
        pdb_path = Path(pdb_file).resolve()
        if not pdb_path.is_file():
            raise FileNotFoundError(f"The file {pdb_file} does not exist.")
        
        logger.info(f"Running Propka {self.propka_version} on {pdb_file}")
        
        # Determine output directory and file paths
        if output_dir:
            output_directory = Path(output_dir).resolve()
            output_directory.mkdir(parents=True, exist_ok=True)
            expected_output_file = output_directory / f"{pdb_path.stem}.pka"
        else:
            output_directory = pdb_path.parent
            expected_output_file = pdb_path.with_suffix('.pka')
        
        # Build command to execute Propka with absolute path
        command = [
            f"propka{self.propka_version}",
            str(pdb_path)  # Use absolute path
        ]
        
        try:
            # Execute the command in the target output directory
            result = subprocess.run(
                command, 
                check=True, 
                capture_output=True, 
                text=True,
                cwd=str(output_directory)
            )
            
            logger.info("Propka executed successfully")
            logger.debug(f"Propka output: {result.stdout}")
            
            # Verify output file was created
            if not expected_output_file.exists():
                raise PreparationError(f"Propka output file not found: {expected_output_file}")
            
            self.last_analysis_file = str(expected_output_file)
            logger.info(f"Propka analysis completed. Output: {expected_output_file}")
            
            return str(expected_output_file)
            
        except subprocess.CalledProcessError as e:
            error_msg = f"Error executing Propka: {e.stderr if e.stderr else str(e)}"
            logger.error(error_msg)
            raise PreparationError(error_msg)
            
        except FileNotFoundError:
            error_msg = f"Propka{self.propka_version} command not found. Please ensure Propka is installed and in PATH."
            logger.error(error_msg)
            raise PreparationError(error_msg)
            
        except Exception as e:
            error_msg = f"Unexpected error occurred: {e}"
            logger.error(error_msg)
            raise PreparationError(error_msg)
    
    def extract_summary(self, propka_file: str, output_file: Optional[str] = None, output_dir: Optional[str] = None) -> str:
        """
        Extract the 'SUMMARY OF THIS PREDICTION' section from a Propka output file.
        
        Args:
            propka_file: Path to the Propka output file (.pka)
            output_file: Path where to save the extracted summary (optional)
            output_dir: Directory where to save the summary file (optional)
            
        Returns:
            Path to the output file containing the summary
            
        Raises:
            FileNotFoundError: If input file doesn't exist
            PreparationError: If summary section is not found
        """
        if not os.path.isfile(propka_file):
            raise FileNotFoundError(f"Propka output file {propka_file} does not exist.")
        
        if output_file is None:
            propka_path = Path(propka_file)
            # Use the source PDB basename for the summary filename
            source_basename = propka_path.stem  # Remove .pka extension
            summary_filename = f"{source_basename}_summary_of_prediction.txt"
            
            if output_dir:
                output_file = str(Path(output_dir) / summary_filename)
            else:
                output_file = str(propka_path.with_name(summary_filename))
        
        logger.info(f"Extracting summary section from {propka_file}")
        
        # Open file and process content
        with open(propka_file, "r", encoding='utf-8') as file:
            lines = file.readlines()

        # Variables to control extraction
        start_delimiter = "SUMMARY OF THIS PREDICTION"
        end_delimiter = "--------------------------------------------------------------------------------------------------------"
        is_summary = False
        summary_content = []

        # Process line by line
        for line in lines:
            if start_delimiter in line:
                is_summary = True  # Start capturing data
            elif end_delimiter in line and is_summary:
                is_summary = False  # Stop capturing data
            elif is_summary:
                summary_content.append(line.strip())  # Save relevant lines

        if not summary_content:
            raise PreparationError(f"Summary section not found in {propka_file}")

        # Save extracted content to new file
        with open(output_file, "w", encoding='utf-8') as file:
            file.write("\n".join(summary_content))

        self.last_summary_file = output_file
        logger.info(f"Summary section saved to {output_file}")
        return output_file
    
    def parse_summary(self, summary_file: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Parse the summary file and extract ionizable residues and their pKa values.
        
        Args:
            summary_file: Path to the summary file (uses last extracted if None)
            
        Returns:
            List of dictionaries containing residue information
            
        Raises:
            FileNotFoundError: If summary file doesn't exist
            PreparationError: If no protonable residues found
        """
        if summary_file is None:
            summary_file = self.last_summary_file
        
        if not summary_file or not os.path.isfile(summary_file):
            raise FileNotFoundError(f"Summary file {summary_file} does not exist.")
        
        logger.info(f"Parsing summary section from {summary_file}")
        
        residues = []
        with open(summary_file, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            
            for line in lines:
                # Skip header lines
                if line.startswith("Group") or not line.strip():
                    continue
                    
                try:
                    parts = line.split()
                    if len(parts) >= 5:
                        residue_name = parts[0]
                        residue_id_or_atom = parts[1]  # Can be residue number (protein) or atom name (ligand)
                        chain_id = parts[2]
                        pka_value = float(parts[3])
                        model_pka = float(parts[4]) if len(parts) > 4 else None
                        # parts[5+] is atom-type classification (for ligands)
                        atom_type = " ".join(parts[5:]).strip() if len(parts) > 5 else ""
                        
                        # Try to determine if this is a protein residue or ligand
                        # Protein residues have integer res_id, ligands have atom names
                        try:
                            res_id = int(residue_id_or_atom)
                            # This is a protein residue
                            atom_name = ""  # Protein residues don't have specific atom in summary
                        except ValueError:
                            # This is a ligand with an atom name
                            res_id = 0  # Ligands don't have meaningful residue numbers in this context
                            atom_name = residue_id_or_atom
                        
                        residue_info = {
                            "residue": residue_name,
                            "res_id": res_id,
                            "chain": chain_id,
                            "pka": pka_value,
                            "atom": atom_name,  # Atom name for ligands, empty for protein residues
                            "atom_type": atom_type,  # Atom type classification (N31, OCO, OP, etc.)
                            "model_pka": model_pka
                        }
                        residues.append(residue_info)
                        
                except (ValueError, IndexError) as e:
                    logger.debug(f"Skipping line due to parsing error: {line.strip()} - {e}")
                    continue
        
        if not residues:
            raise PreparationError("No protonable residues found in summary file")
        
        self.protonable_residues = residues
        logger.info(f"Parsed {len(residues)} protonable residues")
        return residues
    
    def _extract_chain_info_from_pdb(self, pdb_file: str) -> Dict[str, List[str]]:
        """
        Extract chain information for each residue from a PDB file.
        
        Args:
            pdb_file: Path to the PDB file
            
        Returns:
            Dictionary mapping "res_id:res_name" to list of chain IDs
        """
        chain_info = {}
        
        try:
            with open(pdb_file, 'r', encoding='utf-8') as file:
                for line in file:
                    if line.startswith(("ATOM", "HETATM")):
                        try:
                            res_id = int(line[22:26].strip())
                            chain_id = line[21:22].strip()
                            res_name = line[17:20].strip()
                            
                            # Only store for protonable residue types
                            if res_name in ["ASP", "GLU", "HIS", "LYS", "TYR", "CYS", "ARG"]:
                                # Create key using res_id:res_name
                                key = f"{res_id}:{res_name}"
                                if key not in chain_info:
                                    chain_info[key] = []
                                if chain_id not in chain_info[key]:
                                    chain_info[key].append(chain_id)
                                
                        except (ValueError, IndexError):
                            continue
                            
        except Exception as e:
            logger.warning(f"Could not extract chain information from {pdb_file}: {e}")
            
        return chain_info
    
    def get_default_protonation_state(self, residue: Dict[str, Any], ph: float) -> str:
        """
        Determine the default protonation state for a residue based on pKa and pH.
        
        Args:
            residue: Residue dictionary with 'residue', 'res_id', and 'pka' keys
            ph: Target pH value
            
        Returns:
            Three-letter code for the appropriate protonation state
        """
        residue_type = residue['residue']
        pka = residue['pka']
        
        if residue_type not in PROTONATION_STATES:
            return residue_type
        
        states = PROTONATION_STATES[residue_type]
        
        if residue_type in ['ASP', 'GLU']:  # Acidic residues
            return states['protonated'] if ph < pka else states['deprotonated']
        
        elif residue_type == 'HIS':  # Special case for histidine
            if ph < pka:
                return states['protonated']  # HIP
            else:
                return states['neutral_epsilon']  # HIE by default
        
        elif residue_type == 'LYS':  # Lysine
            return states['protonated'] if ph < pka else states['neutral']
        
        elif residue_type == 'ARG':  # Arginine (always protonated in AMBER)
            return states['protonated']  # No neutral form available
        
        elif residue_type == 'TYR':  # Tyrosine
            return states['neutral'] if ph < pka else states['deprotonated']
        
        elif residue_type == 'CYS':  # Cysteine (default to free form)
            return states['free'] if ph < pka else states['deprotonated']
        
        # Default case - return first available state
        return list(states.values())[0]
    
    def get_available_states(self, residue_type: str) -> Dict[str, str]:
        """
        Get available protonation states for a residue type.
        
        Args:
            residue_type: Three-letter residue code
            
        Returns:
            Dictionary of available states
        """
        return PROTONATION_STATES.get(residue_type, {})
    
    def apply_protonation_states(
        self,
        input_pdb: str,
        output_pdb: str,
        ph: float,
        custom_states: Optional[Dict[str, str]] = None,
        residues: Optional[List[Dict[str, Any]]] = None
    ) -> Dict[str, int]:
        """
        Apply protonation states to a PDB file based on pH and/or custom states.
        
        Args:
            input_pdb: Path to input PDB file
            output_pdb: Path for output PDB file
            ph: Target pH value
            custom_states: Optional dictionary mapping residue IDs to custom states
            residues: Optional list of residue data (uses last parsed if None)
            
        Returns:
            Dictionary with counts: {'residue_changes': int, 'record_changes': int}
            
        Raises:
            FileNotFoundError: If input PDB file doesn't exist
            PreparationError: If no residue data available
        """
        if not os.path.isfile(input_pdb):
            raise FileNotFoundError(f"PDB file {input_pdb} does not exist.")
        
        if residues is None:
            residues = self.protonable_residues
        
        if not residues:
            raise PreparationError("No residue data available. Run parse_summary first.")
        
        logger.info(f"Applying protonation states to {input_pdb} for pH {ph}")
        
        custom_states = custom_states or {}
        
        with open(input_pdb, 'r', encoding='utf-8') as file:
            pdb_lines = file.readlines()

        new_pdb = []
        record_modifications = 0
        modified_residues = set()
        
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                res_id = int(line[22:26].strip())
                chain = line[21:22].strip()
                residue_id = f"{res_name}{res_id}"
                residue_chain_id = f"{res_name}{res_id}_{chain}"  # Chain-specific ID for custom states

                # Find matching residue in Propka results (considering chain)
                matching_res = next(
                    (res for res in residues 
                     if res["residue"] == res_name and res["res_id"] == res_id and res["chain"] == chain), 
                    None
                )
                
                if matching_res:
                    # Check for custom state first (try chain-specific, then fallback to chain-agnostic)
                    if residue_chain_id in custom_states:
                        new_name = custom_states[residue_chain_id]
                    elif residue_id in custom_states:
                        new_name = custom_states[residue_id]
                    else:
                        # Use pH-based prediction
                        new_name = self.get_default_protonation_state(matching_res, ph)
                    
                    # Only modify if state changed
                    if new_name != res_name:
                        line = line[:17] + f"{new_name:<3}" + line[20:]
                        record_modifications += 1
                        modified_residues.add(residue_chain_id)  # Track unique residues
            
            new_pdb.append(line)

        # Write modified PDB file
        with open(output_pdb, 'w', encoding='utf-8') as file:
            file.writelines(new_pdb)
        
        residue_changes = len(modified_residues)
        
        logger.info(f"PDB file modified: {residue_changes} residue state changes, {record_modifications} PDB record changes")
        logger.info(f"Modified PDB saved as {output_pdb}")
        
        return {
            'residue_changes': residue_changes,
            'record_changes': record_modifications
        }
    
    def get_residue_statistics(self) -> Dict[str, int]:
        """
        Get statistics about the protonable residues.
        
        Returns:
            Dictionary with residue type counts
        """
        if not self.protonable_residues:
            return {}
        
        stats = {}
        for residue in self.protonable_residues:
            res_type = residue['residue']
            stats[res_type] = stats.get(res_type, 0) + 1
        
        return stats
    
    def get_ph_titration_curve(
        self, 
        ph_range: Tuple[float, float] = (0, 14), 
        ph_step: float = 0.5
    ) -> Dict[str, List[Tuple[float, str]]]:
        """
        Generate titration curves for all protonable residues.
        
        Args:
            ph_range: pH range as (min, max) tuple
            ph_step: Step size for pH values
            
        Returns:
            Dictionary mapping residue IDs to list of (pH, state) tuples
        """
        if not self.protonable_residues:
            return {}
        
        curves = {}
        ph_min, ph_max = ph_range
        ph_values = []
        
        current_ph = ph_min
        while current_ph <= ph_max:
            ph_values.append(current_ph)
            current_ph += ph_step
        
        for residue in self.protonable_residues:
            residue_id = f"{residue['residue']}{residue['res_id']}"
            curve = []
            
            for ph in ph_values:
                state = self.get_default_protonation_state(residue, ph)
                curve.append((ph, state))
            
            curves[residue_id] = curve
        
        return curves
    
    def detect_disulfide_bonds(
        self, 
        pdb_file: str, 
        distance_threshold: float = 2.5
    ) -> List[Tuple[Tuple[str, int], Tuple[str, int]]]:
        """
        Detect potential disulfide bonds between cysteine residues based on distance.
        
        Args:
            pdb_file: Path to PDB file
            distance_threshold: Maximum distance (Å) between sulfur atoms for disulfide bond
            
        Returns:
            List of tuples containing pairs of (residue_name, residue_id) for potential S-S bonds
        """
        if not os.path.isfile(pdb_file):
            raise FileNotFoundError(f"PDB file {pdb_file} does not exist.")
        
        logger.info(f"Detecting potential disulfide bonds in {pdb_file}")
        
        # Extract cysteine sulfur atoms
        cys_sulfurs = []
        
        with open(pdb_file, 'r', encoding='utf-8') as file:
            for line in file:
                if line.startswith(("ATOM", "HETATM")):
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    
                    # Look for sulfur atoms in cysteine residues
                    if res_name == "CYS" and atom_name == "SG":
                        res_id = int(line[22:26].strip())
                        chain = line[21:22].strip()
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        
                        cys_sulfurs.append({
                            'residue': res_name,
                            'res_id': res_id,
                            'chain': chain,
                            'coords': (x, y, z)
                        })
        
        # Calculate distances and find potential bonds
        disulfide_bonds = []
        
        for i in range(len(cys_sulfurs)):
            for j in range(i + 1, len(cys_sulfurs)):
                cys1 = cys_sulfurs[i]
                cys2 = cys_sulfurs[j]
                
                # Calculate distance
                x1, y1, z1 = cys1['coords']
                x2, y2, z2 = cys2['coords']
                distance = ((x2-x1)**2 + (y2-y1)**2 + (z2-z1)**2)**0.5
                
                if distance <= distance_threshold:
                    bond = (
                        (cys1['residue'], cys1['res_id']),
                        (cys2['residue'], cys2['res_id'])
                    )
                    disulfide_bonds.append(bond)
                    logger.info(f"Potential disulfide bond: CYS{cys1['res_id']} - CYS{cys2['res_id']} (distance: {distance:.2f} Å)")
        
        return disulfide_bonds
    
    def apply_disulfide_bonds(
        self,
        input_pdb: str,
        output_pdb: str,
        disulfide_bonds: Optional[List[Tuple[Tuple[str, int], Tuple[str, int]]]] = None,
        auto_detect: bool = True
    ) -> int:
        """
        Apply disulfide bond assignments by changing CYS to CYX for bonded cysteines.
        
        Args:
            input_pdb: Path to input PDB file
            output_pdb: Path for output PDB file
            disulfide_bonds: List of disulfide bond pairs, if None will auto-detect
            auto_detect: Whether to automatically detect disulfide bonds if not provided
            
        Returns:
            Number of actual disulfide bonds applied
        """
        if not os.path.isfile(input_pdb):
            raise FileNotFoundError(f"PDB file {input_pdb} does not exist.")
        
        if disulfide_bonds is None and auto_detect:
            disulfide_bonds = self.detect_disulfide_bonds(input_pdb)
        elif disulfide_bonds is None:
            disulfide_bonds = []
        
        logger.info(f"Applying {len(disulfide_bonds)} disulfide bonds to {input_pdb}")
        
        # Create set of residue IDs that should be CYX
        cyx_residues = set()
        for bond in disulfide_bonds:
            (_, res_id1), (_, res_id2) = bond
            cyx_residues.add(res_id1)
            cyx_residues.add(res_id2)
        
        with open(input_pdb, 'r', encoding='utf-8') as file:
            pdb_lines = file.readlines()

        new_pdb = []
        modifications_made = 0
        
        for line in pdb_lines:
            if line.startswith(("ATOM", "HETATM")):
                res_name = line[17:20].strip()
                res_id = int(line[22:26].strip())
                
                # Change CYS to CYX if in disulfide bond
                if res_name == "CYS" and res_id in cyx_residues:
                    line = line[:17] + "CYX" + line[20:]
                    modifications_made += 1
            
            new_pdb.append(line)

        # Write modified PDB file
        with open(output_pdb, 'w', encoding='utf-8') as file:
            file.writelines(new_pdb)
        
        # Return the actual number of disulfide bonds, not the number of record changes
        num_disulfide_bonds = len(disulfide_bonds)
        logger.info(f"Applied {num_disulfide_bonds} disulfide bonds: {modifications_made} PDB records changed from CYS to CYX")
        logger.info(f"Modified PDB saved as {output_pdb}")
        
        return num_disulfide_bonds
    
    def run_pdb4amber_with_cap_fix(
        self,
        input_pdb: str,
        output_pdb: str,
        fix_caps: bool = True,
        pdb4amber_options: Optional[Dict[str, Any]] = None
    ) -> Dict[str, Any]:
        """
        Run pdb4amber on a PDB file with optional ACE/NME cap HETATM fix.
        
        pdb4amber adds hydrogens and prepares structures for AMBER MD simulations.
        When protein capping (ACE/NME) is present, pdb4amber converts cap ATOM 
        records to HETATM, which causes issues with downstream tools like 
        packmol-memgen. This method optionally fixes that issue.
        
        Args:
            input_pdb: Path to input PDB file
            output_pdb: Path for output PDB file
            fix_caps: If True, automatically convert ACE/NME HETATM back to ATOM 
                     after pdb4amber processing (default: True)
            pdb4amber_options: Optional dictionary of pdb4amber command-line options
                              Example: {'reduce': True, 'dry': False}
        
        Returns:
            Dictionary containing:
                - 'success': bool - Whether pdb4amber succeeded
                - 'output_file': str - Path to output file
                - 'hetatm_fixed': int - Number of HETATM records fixed (if fix_caps=True)
                - 'stdout': str - pdb4amber stdout
                - 'stderr': str - pdb4amber stderr
        
        Raises:
            FileNotFoundError: If input file doesn't exist
            PreparationError: If pdb4amber execution fails
            
        Example:
            >>> analyzer = PreparationManager()
            >>> result = analyzer.run_pdb4amber_with_cap_fix(
            ...     "protein_capped.pdb",
            ...     "protein_prepared.pdb",
            ...     fix_caps=True
            ... )
            >>> print(f"Fixed {result['hetatm_fixed']} cap records")
        """
        if not os.path.isfile(input_pdb):
            raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")
        
        # Prepare pdb4amber command
        cmd = ["pdb4amber", "-i", input_pdb, "-o", output_pdb]
        
        # Add optional arguments
        if pdb4amber_options:
            for key, value in pdb4amber_options.items():
                if isinstance(value, bool):
                    if value:
                        cmd.append(f"--{key}")
                else:
                    cmd.extend([f"--{key}", str(value)])
        
        logger.info(f"Running pdb4amber: {' '.join(cmd)}")
        
        try:
            # Run pdb4amber
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            
            logger.info(f"pdb4amber completed successfully")
            
            # Initialize result dictionary
            result_dict = {
                'success': True,
                'output_file': output_pdb,
                'hetatm_fixed': 0,
                'stdout': result.stdout,
                'stderr': result.stderr
            }
            
            # Apply ACE/NME cap fix if requested
            if fix_caps:
                logger.info("Applying ACE/NME HETATM fix...")
                hetatm_fixed = self._fix_cap_hetatm_records(output_pdb)
                result_dict['hetatm_fixed'] = hetatm_fixed
                
                if hetatm_fixed > 0:
                    logger.info(f"Fixed {hetatm_fixed} HETATM records for ACE/NME caps")
                else:
                    logger.debug("No ACE/NME caps found requiring HETATM fix")
            
            return result_dict
            
        except subprocess.CalledProcessError as e:
            error_msg = f"pdb4amber failed: {e.stderr}"
            logger.error(error_msg)
            raise PreparationError(error_msg)
        except Exception as e:
            error_msg = f"Error running pdb4amber: {str(e)}"
            logger.error(error_msg)
            raise PreparationError(error_msg)
    
    def _fix_cap_hetatm_records(self, pdb_file: str) -> int:
        """
        Fix HETATM records for ACE/NME caps after pdb4amber processing.
        
        pdb4amber converts ACE/NME caps from ATOM to HETATM records, which causes
        issues with downstream tools like packmol-memgen. This method converts
        HETATM records back to ATOM records for ACE and NME residues only.
        
        Args:
            pdb_file: Path to the PDB file to fix (modified in-place)
            
        Returns:
            Number of HETATM records converted to ATOM records
            
        Raises:
            IOError: If file cannot be read or written
            
        Note:
            This method modifies the file in-place. Other HETATM records 
            (water, ions, ligands) are preserved unchanged.
        """
        try:
            fixed_lines = []
            hetatm_fixed_count = 0
            
            with open(pdb_file, 'r') as f:
                for line in f:
                    # Check if line is a HETATM record for ACE or NME
                    if line.startswith('HETATM'):
                        # Extract residue name (columns 18-20 in PDB format)
                        if len(line) >= 20:
                            res_name = line[17:20].strip()
                            if res_name in ['ACE', 'NME']:
                                # Convert HETATM to ATOM (preserve spacing)
                                fixed_line = 'ATOM  ' + line[6:]
                                fixed_lines.append(fixed_line)
                                hetatm_fixed_count += 1
                                continue
                    
                    # Keep line as-is if not a ACE/NME HETATM record
                    fixed_lines.append(line)
            
            # Write fixed content back to file
            with open(pdb_file, 'w') as f:
                f.writelines(fixed_lines)
            
            logger.debug(f"Fixed {hetatm_fixed_count} HETATM records for ACE/NME caps in {pdb_file}")
            return hetatm_fixed_count
            
        except IOError as e:
            logger.error(f"Error fixing HETATM records in {pdb_file}: {e}")
            raise

# Convenience functions for backward compatibility
def run_propka(entrada_file: str, propka_version: str = '3') -> None:
    """
    Execute Propka to predict pKa values of a molecule.
    
    Args:
        entrada_file: Path to the input PDB file
        propka_version: Version of Propka to use (default: '3')
        
    Raises:
        PreparationError: If Propka execution fails
        FileNotFoundError: If input file doesn't exist
    """
    analyzer = PreparationManager(propka_version)
    analyzer.run_analysis(entrada_file)

def extract_summary_section(input_file: str, output_file: str) -> str:
    """
    Extract the 'SUMMARY OF THIS PREDICTION' section from a PropKa output file.
    
    Args:
        input_file: Path to the PropKa output file (.pka)
        output_file: Path where to save the extracted summary
        
    Returns:
        Path to the output file containing the summary
    """
    analyzer = PreparationManager()
    return analyzer.extract_summary(input_file, output_file)

def parse_summary_section(summary_file: str) -> List[Dict[str, Any]]:
    """
    Read PropKa output and extract ionizable residues and their pKa values.
    
    Args:
        summary_file: Path to the summary file containing PropKa predictions
        
    Returns:
        List of dictionaries containing residue information
    """
    analyzer = PreparationManager()
    return analyzer.parse_summary(summary_file)

def modify_pdb_based_on_summary(
    pdb_file: str, 
    residues: List[Dict[str, Any]], 
    ph: float, 
    output_pdb: str,
    custom_states: Optional[Dict[str, str]] = None
) -> None:
    """
    Modify PDB file according to pKa values from summary section and desired pH.
    
    Args:
        pdb_file: Path to input PDB file
        residues: List of residue dictionaries from parse_summary_section
        ph: Target pH value
        output_pdb: Path for output PDB file
        custom_states: Optional dictionary mapping residue IDs to custom protonation states
    """
    analyzer = PreparationManager()
    analyzer.protonable_residues = residues
    analyzer.apply_protonation_states(pdb_file, output_pdb, ph, custom_states, residues)

def run_pdb4amber(
    env_path: str, 
    pdb4amber_path: str, 
    input_pdb: str, 
    output_pdb: str
) -> None:
    """
    Run pdb4amber to prepare PDB file for Amber.
    
    Args:
        env_path: Path to conda/virtual environment
        pdb4amber_path: Path to pdb4amber executable
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
        
    Raises:
        PreparationError: If pdb4amber execution fails
    """
    logger.info(f"Running pdb4amber on {input_pdb}")
    
    command = [pdb4amber_path, "-i", input_pdb, "-o", output_pdb]
    
    try:
        # Set environment variables if needed
        env = os.environ.copy()
        if env_path:
            env["PATH"] = f"{env_path}/bin:{env.get('PATH', '')}"
        
        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True,
            env=env
        )
        
        logger.info("pdb4amber executed successfully")
        logger.debug(f"pdb4amber output: {result.stdout}")
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Error executing pdb4amber: {e.stderr if e.stderr else str(e)}"
        logger.error(error_msg)
        raise PreparationError(error_msg)
        
    except FileNotFoundError:
        error_msg = f"pdb4amber not found at {pdb4amber_path}"
        logger.error(error_msg)
        raise PreparationError(error_msg)

def convert_to_namd(
    top_file: str, 
    crd_file: str, 
    psf_output: str = 'system.psf', 
    pdb_output: str = 'system.pdb'
) -> None:
    """
    Convert Amber files to NAMD format using ParmEd.
    
    Args:
        top_file: Amber topology file (.top)
        crd_file: Amber coordinate file (.crd)  
        psf_output: Output PSF file name
        pdb_output: Output PDB file name
        
    Raises:
        PreparationError: If conversion fails or ParmEd not available
    """
    try:
        import parmed as pmd
        
        logger.info(f"Converting {top_file} and {crd_file} to NAMD format")
        
        # Load Amber files
        system = pmd.load_file(top_file, crd_file)
        
        # Save in NAMD format
        system.save(psf_output, overwrite=True)
        system.save(pdb_output, overwrite=True)
        
        logger.info(f"Files converted and saved as {psf_output} and {pdb_output}")
        
    except ImportError:
        raise PreparationError(
            "ParmEd is required for NAMD conversion. "
            "Install with: pip install parmed"
        )
    except Exception as e:
        raise PreparationError(f"Error converting to NAMD format: {e}")

def convert_to_amber(
    top_file: str, 
    crd_file: str, 
    prmtop_output: str = 'system.prmtop', 
    inpcrd_output: str = 'system.inpcrd'
) -> None:
    """
    Convert Amber top/crd files to prmtop/inpcrd format using ParmEd.
    
    Args:
        top_file: Amber topology file (.top)
        crd_file: Amber coordinate file (.crd)  
        prmtop_output: Output parameter topology file name (.prmtop)
        inpcrd_output: Output restart file name (.inpcrd)
        
    Raises:
        PreparationError: If conversion fails
    """
    try:
        import parmed as pmd
        
        logger.info(f"Converting {top_file} and {crd_file} to Amber prmtop/inpcrd format")
        
        # Load Amber files
        system = pmd.load_file(top_file, crd_file)
        
        # Save in Amber prmtop/inpcrd format
        system.save(prmtop_output, overwrite=True)
        system.save(inpcrd_output, overwrite=True)
        
        logger.info(f"Files converted and saved as {prmtop_output} and {inpcrd_output}")
        
    except ImportError:
        raise PreparationError(
            "ParmEd is required for Amber conversion. "
            "Install with: pip install parmed"
        )
    except Exception as e:
        raise PreparationError(f"Error converting to Amber format: {e}")