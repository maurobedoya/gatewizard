# gatewizard/utils/protein_capping.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Protein capping utility for adding ACE and NME groups to protein termini.

This module provides functionality to remove hydrogens and add ACE (N-terminal)
and NME (C-terminal) capping groups to protein structures before Propka analysis.

Based on the script by Mohd Ibrahim (Technical University of Munich).
"""

import numpy as np
import warnings
from pathlib import Path
from typing import Optional, Tuple, List, Union
from tempfile import NamedTemporaryFile

import MDAnalysis as mda
from MDAnalysis.core.universe import Universe as MDAUniverse

from gatewizard.utils.logger import get_logger

# Suppress MDAnalysis warnings
warnings.filterwarnings("ignore")

logger = get_logger(__name__)

class ProteinCappingError(Exception):
    """Custom exception for protein capping errors."""
    pass

class ProteinCapper:
    """
    Class for adding ACE and NME capping groups to protein termini.
    
    This class provides methods to remove hydrogens and add terminal capping
    groups to prepare proteins for Propka analysis.
    """
    
    def __init__(self):
        """Initialize the protein capper."""
        pass
    
    def remove_hydrogens_and_cap(
        self, 
        input_file: Union[str, Path], 
        output_file: Optional[Union[str, Path]] = None,
        target_dir: Optional[Union[str, Path]] = None
    ) -> Tuple[str, dict]:
        """
        Remove hydrogens from protein and add ACE/NME capping groups.
        
        Args:
            input_file: Path to input PDB file
            output_file: Path to output PDB file (optional, auto-generated if None)
            target_dir: Directory for intermediate files (optional, uses output_file dir if None)
            
        Returns:
            Tuple of (output_file_path, residue_mapping)
            residue_mapping: Dict mapping original residue IDs to new IDs
            Format: {(resname, chain, old_resid): (resname, chain, new_resid)}
            
        Raises:
            ProteinCappingError: If capping process fails
            FileNotFoundError: If input file doesn't exist
        """
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        if output_file is None:
            output_file = input_path.with_name(f"{input_path.stem}_capped{input_path.suffix}")
        
        # Determine target directory for intermediate files
        if target_dir:
            mapping_dir = Path(target_dir)
        else:
            mapping_dir = Path(output_file).parent
        
        logger.info(f"Starting protein capping process for {input_file}")
        
        try:
            # Step 1: Remove hydrogens
            temp_no_h_file = self._remove_hydrogens(input_file)
            
            # Step 2: Add capping groups and create residue mapping
            residue_mapping = self._add_capping_groups(temp_no_h_file, output_file, input_file)
            
            # Step 3: Generate GateWizard mapping file with custom filename based on final output
            output_basename = Path(output_file).stem
            mapping_filename = f"{output_basename}_gatewizard_residue_mapping.txt"
            mapping_file_path = self._generate_gatewizard_mapping_file(
                mapping_dir, residue_mapping, input_file, mapping_filename
            )
            
            # Clean up temporary file
            Path(temp_no_h_file).unlink(missing_ok=True)
            
            logger.info(f"Protein capping completed. Output: {output_file}")
            logger.info(f"GateWizard mapping file created: {mapping_file_path}")
            logger.debug(f"Residue mapping created with {len(residue_mapping)} entries")
            return str(output_file), residue_mapping
            
        except Exception as e:
            logger.error(f"Error during protein capping: {e}", exc_info=True)
            raise ProteinCappingError(f"Protein capping failed: {str(e)}")
    
    def _remove_hydrogens(self, input_file: Union[str, Path]) -> str:
        """
        Remove hydrogen atoms from the protein structure.
        
        Args:
            input_file: Path to input PDB file
            
        Returns:
            Path to temporary file without hydrogens
        """
        logger.debug("Removing hydrogen atoms")
        
        try:
            u = mda.Universe(str(input_file))
            
            # Select protein atoms without hydrogens
            protein_no_h = u.select_atoms("protein and not name H* 1H* 2H* 3H*")
            
            if len(protein_no_h) == 0:
                raise ProteinCappingError("No protein atoms found or all atoms are hydrogens")
            
            # Create temporary file
            with NamedTemporaryFile(mode='w', suffix='.pdb', delete=False) as temp_file:
                temp_filename = temp_file.name
            
            protein_no_h.write(temp_filename)
            logger.debug(f"Hydrogens removed, temporary file: {temp_filename}")
            
            return temp_filename
            
        except Exception as e:
            raise ProteinCappingError(f"Failed to remove hydrogens: {str(e)}")
    
    def _extract_original_residues(self, original_file: Union[str, Path]) -> dict:
        """
        Extract original residue information from the input PDB file.
        
        Args:
            original_file: Path to original PDB file
            
        Returns:
            Dict mapping (chain, old_resid) to resname for tracking
        """
        original_residues = {}
        
        with open(original_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    resname = line[17:20].strip()
                    chain = line[21:22].strip()
                    resid = int(line[22:26].strip())
                    key = (chain, resid)
                    if key not in original_residues:
                        original_residues[key] = resname
        
        return original_residues
    
    def _add_capping_groups(self, input_file: str, output_file: Union[str, Path], original_file: Union[str, Path]) -> dict:
        """
        Add ACE and NME capping groups to protein termini and create residue mapping.
        
        Args:
            input_file: Path to input PDB file (without hydrogens)
            output_file: Path to output PDB file
            original_file: Path to original PDB file (for residue mapping)
            
        Returns:
            Dict mapping original residue IDs to new IDs
            Format: {(resname, chain, old_resid): (resname, chain, new_resid)}
        """
        logger.debug("Adding ACE and NME capping groups")
        
        try:
            # Ensure output directory exists
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            # Create mapping from original residue numbers
            original_residue_mapping = self._extract_original_residues(original_file)
            
            # Group original residues by chain
            chains_in_original = {}
            for (chain, resid), resname in original_residue_mapping.items():
                if chain not in chains_in_original:
                    chains_in_original[chain] = []
                chains_in_original[chain].append(resid)
            
            u = mda.Universe(input_file)
            segment_universes = []
            res_start = 0
            residue_mapping = {}
            
            # Get list of chains to process in order
            chain_list = sorted(chains_in_original.keys())
            seg_to_chain_map = {}
            
            # Map segments to chains based on order
            for i, seg in enumerate(u.segments):  # type: ignore[union-attr]
                if i < len(chain_list):
                    seg_to_chain_map[seg.segid] = chain_list[i]
                else:
                    logger.warning(f"More segments than original chains, using segment {seg.segid} as chain ID")
                    seg_to_chain_map[seg.segid] = seg.segid
            
            for seg in u.segments:  # type: ignore[union-attr]
                chain = u.select_atoms(f"segid {seg.segid}")
                original_chain_id = seg_to_chain_map.get(seg.segid, seg.segid)
                
                # Add ACE capping group to N-terminus
                resid_n = chain.residues.resids[0]
                n_term_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_n}")
                ace_positions = self._get_ace_positions(n_term_residue)
                
                ace_universe = self._create_universe(
                    n_atoms=len(ace_positions),
                    names=["C", "CH3", "O"],
                    resnames=["ACE"] * len(ace_positions),
                    positions=ace_positions,
                    resids=[resid_n] * len(ace_positions),
                    segid=chain.segids[0]
                )
                
                # Add NME capping group to C-terminus
                resid_c = chain.residues.resids[-1]
                c_term_residue = u.select_atoms(f"segid {seg.segid} and resid {resid_c}")
                nme_positions = self._get_nme_positions(c_term_residue)
                
                nme_universe = self._create_universe(
                    n_atoms=len(nme_positions),
                    names=["N", "C"],
                    resnames=["NME"] * len(nme_positions),
                    positions=nme_positions,
                    resids=[resid_c + 2] * len(nme_positions),
                    segid=chain.segids[0]
                )
                
                # Handle OXT removal if present
                if "OXT" in c_term_residue.names:
                    oxt_index = np.where(c_term_residue.names == "OXT")[0][0]
                    oxt_atom = c_term_residue[oxt_index]
                    chain_atoms = u.select_atoms(f"segid {seg.segid} and not index {oxt_atom.index}")
                else:
                    chain_atoms = u.select_atoms(f"segid {seg.segid}")
                
                # Merge ACE, protein chain, and NME
                capped_chain = mda.Merge(ace_universe.atoms, chain_atoms, nme_universe.atoms)
                
                # Renumber residues and create mapping
                n_ace_atoms = len(ace_positions)
                n_chain_atoms = len(chain_atoms.residues)
                n_nme_atoms = len(nme_positions)
                
                resids_ace = [res_start + 1] * n_ace_atoms
                resids_chain = list(range(res_start + 2, res_start + 2 + n_chain_atoms))
                resids_nme = [res_start + 2 + n_chain_atoms] * n_nme_atoms
                
                all_resids = resids_ace + resids_chain + resids_nme
                capped_chain.atoms.residues.resids = np.array(all_resids)  # type: ignore[union-attr]
                
                # Create residue mapping for this chain
                chain_id = original_chain_id  # Use the mapped chain ID
                original_resids = sorted(chain.residues.resids)
                
                # Map original residues to new residues (after ACE cap)
                for i, original_resid in enumerate(original_resids):
                    new_resid = res_start + 2 + i  # +2 because ACE takes position +1
                    original_key = (chain_id, original_resid)
                    
                    if original_key in original_residue_mapping:
                        original_resname = original_residue_mapping[original_key]
                        # Map: (original_resname, chain, original_resid) -> (resname, chain, new_resid)
                        mapping_key = (original_resname, chain_id, original_resid)
                        mapping_value = (original_resname, chain_id, new_resid)
                        residue_mapping[mapping_key] = mapping_value
                    else:
                        logger.warning(f"No original mapping found for {chain_id} {original_resid}")
                
                res_start = max(all_resids)
                segment_universes.append(capped_chain)
            
            # Merge all segments
            final_universe = mda.Merge(*(seg.atoms for seg in segment_universes))
            final_universe.atoms.write(str(output_file))  # type: ignore[union-attr]
            
            return residue_mapping
            
        except Exception as e:
            raise ProteinCappingError(f"Failed to add capping groups: {str(e)}")
    
    def _generate_gatewizard_mapping_file(self, output_dir: Union[str, Path], 
                                        residue_mapping: dict, 
                                        original_file: Union[str, Path],
                                        mapping_filename: Optional[str] = None) -> str:
        """
        Generate a GateWizard mapping file similar to tleap's system_for_tleap_renum.txt.
        This file shows the mapping from original residue numbering to final numbering
        after capping, including ACE/NME caps.
        
        Args:
            output_dir: Directory where to save the mapping file
            residue_mapping: Dictionary with residue mappings from _add_capping_groups
            original_file: Path to original PDB file
            mapping_filename: Custom filename for the mapping file (optional)
            
        Returns:
            Path to the generated mapping file
        """
        try:
            output_path = Path(output_dir)
            
            # Ensure the output directory exists
            output_path.mkdir(parents=True, exist_ok=True)
            
            if mapping_filename:
                mapping_file = output_path / mapping_filename
            else:
                mapping_file = output_path / "gatewizard_residue_mapping.txt"
            
            # Extract original residues to get complete structure info
            original_residue_mapping = self._extract_original_residues(original_file)
            
            # Group by chain for ordered output
            chains_data = {}
            
            # First, collect all mapped residues
            for (original_resname, chain, original_resid), (final_resname, final_chain, final_resid) in residue_mapping.items():
                if chain not in chains_data:
                    chains_data[chain] = []
                chains_data[chain].append({
                    'original_resname': original_resname,
                    'original_resid': original_resid,
                    'final_resname': final_resname,
                    'final_resid': final_resid,
                    'is_cap': False
                })
            
            # Add ACE and NME caps to the mapping
            for chain in sorted(chains_data.keys()):
                chain_residues = sorted(chains_data[chain], key=lambda x: x['original_resid'])
                
                if chain_residues:
                    # Add ACE cap (N-terminal) - it gets resid 1 less than first residue's final_resid
                    first_final_resid = min(r['final_resid'] for r in chain_residues)
                    ace_resid = first_final_resid - 1
                    
                    chains_data[chain].insert(0, {
                        'original_resname': '-',  # ACE didn't exist originally
                        'original_resid': '-',    # No original residue ID
                        'final_resname': 'ACE',
                        'final_resid': ace_resid,
                        'is_cap': True
                    })
                    
                    # Add NME cap (C-terminal) - it gets resid 1 more than last residue's final_resid
                    last_final_resid = max(r['final_resid'] for r in chain_residues if not r['is_cap'])
                    nme_resid = last_final_resid + 1
                    
                    chains_data[chain].append({
                        'original_resname': '-',  # NME didn't exist originally
                        'original_resid': '-',    # No original residue ID
                        'final_resname': 'NME',
                        'final_resid': nme_resid,
                        'is_cap': True
                    })
            
            # Write the mapping file
            with open(mapping_file, 'w') as f:
                f.write("# GateWizard Residue Mapping File\n")
                f.write("# Generated when protein capping is enabled\n")
                f.write("# Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID\n")
                f.write("# '-' indicates caps (ACE/NME) that did not exist in original structure\n")
                f.write("#\n")
                
                # Sort chains and write data
                for chain in sorted(chains_data.keys()):
                    # Sort by final residue ID to match tleap output format
                    chain_residues = sorted(chains_data[chain], key=lambda x: x['final_resid'])
                    
                    for residue in chain_residues:
                        line = f"{residue['original_resname']:>3} {chain}  {str(residue['original_resid']):>5}    {residue['final_resname']:>3}   {residue['final_resid']:>3}\n"
                        f.write(line)
            
            logger.info(f"Generated GateWizard mapping file: {mapping_file}")
            return str(mapping_file)
            
        except Exception as e:
            logger.error(f"Failed to generate GateWizard mapping file: {str(e)}")
            raise ProteinCappingError(f"Failed to generate mapping file: {str(e)}")
    
    def _create_universe(
        self, 
        n_atoms: int, 
        names: List[str], 
        resnames: List[str], 
        positions: List[np.ndarray], 
        resids: List[int], 
        segid: str
    ) -> MDAUniverse:
        """
        Create a new MDAnalysis universe with specified atoms.
        
        Args:
            n_atoms: Number of atoms
            names: Atom names
            resnames: Residue names
            positions: Atom positions
            resids: Residue IDs
            segid: Segment ID
            
        Returns:
            MDAnalysis Universe object
        """
        u_new = mda.Universe.empty(
            n_atoms=n_atoms,
            n_residues=n_atoms,
            atom_resindex=np.arange(n_atoms),
            residue_segindex=np.arange(n_atoms),
            n_segments=n_atoms,
            trajectory=True
        )
        
        u_new.add_TopologyAttr('name', names)
        u_new.add_TopologyAttr('resid', resids)
        u_new.add_TopologyAttr('resname', resnames)
        u_new.atoms.positions = np.array(positions)  # type: ignore[union-attr]
        u_new.add_TopologyAttr('segid', [segid] * n_atoms)
        u_new.add_TopologyAttr('chainID', [segid] * n_atoms)
        
        return u_new
    
    def _get_ace_positions(self, n_term_residue) -> List[np.ndarray]:
        """
        Calculate positions for ACE capping group atoms.
        
        Args:
            n_term_residue: N-terminal residue selection
            
        Returns:
            List of positions for ACE atoms [C, CH3, O]
        """
        try:
            c_idx = np.where(n_term_residue.names == "C")[0][0]
            ca_idx = np.where(n_term_residue.names == "CA")[0][0]
            n_idx = np.where(n_term_residue.names == "N")[0][0]
            
            pos_c = n_term_residue.positions[c_idx]
            pos_ca = n_term_residue.positions[ca_idx]
            pos_n = n_term_residue.positions[n_idx]
            
            # Calculate ACE positions using geometry
            # C-CA-N-C dihedral has two minima at -60,60
            c_position = self._calc_coordinate(pos_c, pos_ca, pos_n, 1.34, 120, -60)
            ch3_position = self._calc_coordinate(pos_ca, pos_n, c_position, 1.52, 120, 180)
            o_position = self._calc_coordinate(pos_ca, pos_n, c_position, 1.23, 120, 0)
            
            return [c_position, ch3_position, o_position]
            
        except Exception as e:
            raise ProteinCappingError(f"Failed to calculate ACE positions: {str(e)}")
    
    def _get_nme_positions(self, c_term_residue) -> List[np.ndarray]:
        """
        Calculate positions for NME capping group atoms.
        
        Args:
            c_term_residue: C-terminal residue selection
            
        Returns:
            List of positions for NME atoms [N, C]
        """
        try:
            o_idx = np.where(c_term_residue.names == "O")[0][0]
            ca_idx = np.where(c_term_residue.names == "CA")[0][0]
            c_idx = np.where(c_term_residue.names == "C")[0][0]
            
            pos_o = c_term_residue.positions[o_idx]
            pos_ca = c_term_residue.positions[ca_idx]
            pos_c = c_term_residue.positions[c_idx]
            
            # Calculate bisector of O-C-CA angle
            v1 = pos_o - pos_c
            v1 /= np.linalg.norm(v1)
            v2 = pos_ca - pos_c
            v2 /= np.linalg.norm(v2)
            bisector = v1 + v2
            bisector /= np.linalg.norm(bisector)
            
            # Place NME N atom
            bond_length = 1.34
            n_position = bond_length * (-bisector) + pos_c
            
            # Place NME C atom using dihedral of 0 degrees
            bond_length = 1.45
            c_position = self._calc_coordinate(pos_o, pos_c, n_position, bond_length, 120, 0)
            
            return [n_position, c_position]
            
        except Exception as e:
            raise ProteinCappingError(f"Failed to calculate NME positions: {str(e)}")
    
    def _calc_coordinate(
        self, 
        a: np.ndarray, 
        b: np.ndarray, 
        c: np.ndarray, 
        bond_len: float, 
        theta: float, 
        dihedral: float
    ) -> np.ndarray:
        """
        Calculate position of fourth atom 'd' based on atoms a, b, c and geometry.
        
        Args:
            a, b, c: Positions of three reference atoms
            bond_len: c-d bond length (angstrom)
            theta: b-c-d angle (degrees)
            dihedral: a-b-c-d dihedral angle (degrees)
            
        Returns:
            Position of fourth atom
        """
        dihedral_rad = np.deg2rad(dihedral)
        
        u = c - b
        x = a - b
        v = x - (np.dot(x, u) / np.dot(u, u)) * u
        
        w = np.cross(u, x)
        
        q = (v / np.linalg.norm(v)) * np.cos(dihedral_rad)
        e = (w / np.linalg.norm(w)) * np.sin(dihedral_rad)
        
        pos_temp2 = b + (q + e)
        
        u1 = b - c
        y1 = pos_temp2 - c
        
        mag_y1 = np.linalg.norm(y1)
        mag_u1 = np.linalg.norm(u1)
        
        theta_bcd = np.arccos(np.dot(u1, y1) / (mag_u1 * mag_y1))
        rotate = np.deg2rad(theta) - theta_bcd
        
        z = np.cross(u1, y1)
        n = z / np.linalg.norm(z)
        
        pos_ini = (c + y1 * np.cos(rotate) + 
                   np.cross(n, y1) * np.sin(rotate) +
                   n * np.dot(n, y1) * (1 - np.cos(rotate)))
        
        d = ((pos_ini - c) * (bond_len / np.linalg.norm(pos_ini - c))) + c
        
        return d

def cap_protein(
    input_file: Union[str, Path], 
    output_file: Optional[Union[str, Path]] = None
) -> Tuple[str, dict]:
    """
    Convenience function to cap a protein with ACE and NME groups.
    
    Args:
        input_file: Path to input PDB file
        output_file: Path to output PDB file (optional)
        
    Returns:
        Tuple of (capped_file_path, residue_mapping)
        residue_mapping: Dict mapping original to new residue IDs
        
    Raises:
        ProteinCappingError: If capping process fails
    """
    capper = ProteinCapper()
    return capper.remove_hydrogens_and_cap(input_file, output_file)