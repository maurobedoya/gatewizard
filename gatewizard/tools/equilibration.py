# gatewizard/tools/equilibration.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Equilibration tools for molecular dynamics simulations.

This module provides tools for generating equilibration protocols and
input files for various molecular dynamics engines.
"""

import os
import subprocess
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import json
import tempfile

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class EquilibrationProtocol:
    """Base class for equilibration protocols."""
    
    def __init__(self, name: str, description: str = ""):
        self.name = name
        self.description = description
        self.stages = []
    
    def add_stage(self, stage: Dict[str, Any]):
        """Add an equilibration stage to the protocol."""
        self.stages.append(stage)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert protocol to dictionary."""
        return {
            "name": self.name,
            "description": self.description,
            "stages": self.stages
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "EquilibrationProtocol":
        """Create protocol from dictionary."""
        protocol = cls(data["name"], data.get("description", ""))
        protocol.stages = data.get("stages", [])
        return protocol

class NAMDEquilibrationManager:
    """Manager for NAMD equilibration simulations."""
    
    def __init__(self, working_dir: Path, namd_executable: str = "namd3"):
        self.working_dir = Path(working_dir)
        self.namd_executable = namd_executable
        self.logger = get_logger(self.__class__.__name__)
        
        # Path to CHARMM-GUI templates
        self.charmm_gui_templates_dir = Path(__file__).parent.parent.parent / "equilibration" / "namd" / "charmm_gui_templates"
    
    def find_system_files(self) -> Optional[Dict[str, str]]:
        """
        Automatically detect system files in working directory.
        
        Looks for standard AMBER system files and bilayer PDB with CRYST1 record.
        
        Returns:
            Dictionary with system file paths, or None if required files not found:
            {
                'prmtop': Path to system.prmtop,
                'inpcrd': Path to system.inpcrd,
                'pdb': Path to system.pdb,
                'bilayer_pdb': Path to bilayer PDB with CRYST1
            }
        
        Example:
            >>> manager = NAMDEquilibrationManager(Path("work_dir"))
            >>> system_files = manager.find_system_files()
            >>> if system_files:
            ...     result = manager.setup_namd_equilibration(
            ...         system_files=system_files,
            ...         stage_params_list=stages
            ...     )
        """
        system_files = {}
        
        # Find AMBER topology file
        prmtop_files = list(self.working_dir.glob("*.prmtop"))
        if not prmtop_files:
            self.logger.error("No .prmtop file found in working directory")
            return None
        system_files['prmtop'] = str(prmtop_files[0])
        self.logger.info(f"Found topology: {prmtop_files[0].name}")
        
        # Find AMBER coordinate file
        inpcrd_files = list(self.working_dir.glob("*.inpcrd"))
        if not inpcrd_files:
            # Try alternative extensions
            inpcrd_files = list(self.working_dir.glob("*.crd"))
            if not inpcrd_files:
                inpcrd_files = list(self.working_dir.glob("*.rst"))
        
        if not inpcrd_files:
            self.logger.error("No .inpcrd/.crd/.rst file found in working directory")
            return None
        system_files['inpcrd'] = str(inpcrd_files[0])
        self.logger.info(f"Found coordinates: {inpcrd_files[0].name}")
        
        # Find system PDB file
        system_pdb = self.working_dir / "system.pdb"
        if not system_pdb.exists():
            # Try to find any PDB that's not a bilayer file
            pdb_files = [f for f in self.working_dir.glob("*.pdb") 
                        if "bilayer" not in f.name.lower()]
            if pdb_files:
                system_pdb = pdb_files[0]
            else:
                self.logger.error("No system.pdb file found in working directory")
                return None
        system_files['pdb'] = str(system_pdb)
        self.logger.info(f"Found system PDB: {system_pdb.name}")
        
        # Find bilayer PDB with CRYST1 record
        bilayer_pdb = self._find_bilayer_pdb_with_cryst1()
        if not bilayer_pdb:
            self.logger.error("No bilayer PDB with CRYST1 record found in working directory")
            self.logger.error("Required: bilayer*_lipid.pdb file from packmol-memgen --parametrize")
            return None
        system_files['bilayer_pdb'] = str(bilayer_pdb)
        self.logger.info(f"Found bilayer PDB with CRYST1: {bilayer_pdb.name}")
        
        return system_files
    
    def _get_config_name(self, stage_name: str, stage_index: int) -> str:
        """
        Convert GUI display names to valid config file names.
        
        Maps names like 'Equilibration 1', 'Equilibration 2', etc. to 'step1', 'step2', etc.
        For custom names, uses stage_index to generate sequential step names.
        
        Args:
            stage_name: Display name from GUI (e.g., "Equilibration 1" or "Initial Equilibration")
            stage_index: Zero-based index of the stage in the protocol
            
        Returns:
            Valid config file name (e.g., "step1", "step2", etc.)
        """
        # Handle the standard naming convention with spaces "Equilibration N"
        if stage_name.startswith("Equilibration "):
            try:
                stage_num = stage_name.split()[1]
                return f"step{stage_num}"
            except (IndexError, ValueError):
                pass
        
        # Handle Production stage
        if stage_name == "Production":
            return "step7_production"
        
        # Handle legacy names (in case they exist) - convert to new convention
        legacy_mapping = {
            "equilibration_1": "step1",
            "equilibration_2": "step2", 
            "equilibration_3": "step3",
            "equilibration_4": "step4",
            "equilibration_5": "step5",
            "equilibration_6": "step6",
            "eq1": "step1",  # Convert old naming
            "eq2": "step2",
            "eq3": "step3",
            "eq4": "step4",
            "eq5": "step5",
            "eq6": "step6",
            "production": "step7_production"
        }
        
        # Check legacy mapping first
        if stage_name.lower() in legacy_mapping:
            return legacy_mapping[stage_name.lower()]
        
        # For custom stage names, use stage_index to generate sequential step names
        # This ensures proper ordering: step1, step2, step3, etc.
        return f"step{stage_index + 1}"
    
    def _read_box_dimensions(self, pdb_file: Path) -> Tuple[float, float, float]:
        """
        Read box dimensions from PDB file.
        
        First tries to read CRYST1 record, then estimates from coordinates.
        Example CRYST1 line: CRYST1   70.335   70.833   85.067  90.00  90.00  90.00 P 1           1
        
        Args:
            pdb_file: Path to PDB file
            
        Returns:
            Tuple of (a, b, c) box dimensions in Angstroms
        """
        try:
            with open(pdb_file, 'r') as f:
                for line in f:
                    if line.startswith('CRYST1'):
                        # CRYST1 record contains unit cell parameters
                        # Format: CRYST1    a       b       c     alpha  beta  gamma sgroup    z
                        #         CRYST1 70.335  70.833  85.067  90.00  90.00  90.00 P 1        1
                        a = float(line[6:15].strip())
                        b = float(line[15:24].strip()) 
                        c = float(line[24:33].strip())
                        self.logger.info(f"Read box dimensions from CRYST1 record: {a:.2f} x {b:.2f} x {c:.2f} Å")
                        return a, b, c
        except (FileNotFoundError, ValueError, IndexError) as e:
            self.logger.warning(f"Could not read CRYST1 from {pdb_file}: {e}")
        
        # If no CRYST1 record found, estimate from coordinates
        try:
            self.logger.info(f"No CRYST1 record found in {pdb_file.name}, estimating box dimensions from coordinates")
            return self._estimate_box_from_coordinates(pdb_file)
        except Exception as e:
            self.logger.warning(f"Could not estimate box dimensions from {pdb_file}: {e}")
        
        # Return default dimensions for membrane systems (approximate)
        self.logger.info("Using default box dimensions (100x100x100 Å)")
        return 100.0, 100.0, 100.0
    
    def _estimate_box_from_coordinates(self, pdb_file: Path) -> Tuple[float, float, float]:
        """
        Estimate box dimensions from coordinate extremes with padding.
        
        Args:
            pdb_file: Path to PDB file
            
        Returns:
            Tuple of (a, b, c) box dimensions in Angstroms
        """
        x_coords, y_coords, z_coords = [], [], []
        
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        x_coords.append(x)
                        y_coords.append(y)
                        z_coords.append(z)
                    except (ValueError, IndexError):
                        continue
        
        if not x_coords:
            raise ValueError("No valid coordinates found in PDB file")
        
        # Calculate ranges and add padding (10 Å on each side)
        padding = 10.0
        x_range = max(x_coords) - min(x_coords) + 2 * padding
        y_range = max(y_coords) - min(y_coords) + 2 * padding
        z_range = max(z_coords) - min(z_coords) + 2 * padding
        
        self.logger.info(f"Estimated box dimensions: {x_range:.2f} x {y_range:.2f} x {z_range:.2f} Å")
        return x_range, y_range, z_range
    
    def _calculate_pme_grid_size(self, box_dimension: float, cutoff: float = 9.0) -> int:
        """
        Calculate optimal PME grid size for a given box dimension.
        
        The PME grid size should be:
        1. At least 2x the cutoff distance in grid points
        2. Efficiently factorable (preferably powers of 2, 3, 5)
        3. About 1 Å per grid point for good accuracy
        
        Args:
            box_dimension: Box dimension in Angstroms
            cutoff: Electrostatic cutoff in Angstroms
            
        Returns:
            Optimal PME grid size
        """
        # Rule of thumb: ~1 Å per grid point, but at least 2x cutoff
        min_grid_size = max(int(box_dimension), int(2 * cutoff))
        
        # Find the next efficiently factorable number
        # NAMD/FFTW work best with numbers that factor into small primes (2, 3, 5)
        grid_size = self._find_efficient_grid_size(min_grid_size)
        
        self.logger.debug(f"Box dimension: {box_dimension:.2f} Å, "
                         f"Min grid size: {min_grid_size}, "
                         f"Optimal grid size: {grid_size}")
        
        return grid_size
    
    def _find_efficient_grid_size(self, min_size: int) -> int:
        """
        Find the smallest efficiently factorable number >= min_size.
        
        Efficient numbers for FFT are those that factor into small primes (2, 3, 5).
        
        Args:
            min_size: Minimum required grid size
            
        Returns:
            Efficient grid size >= min_size
        """
        if min_size <= 1:
            return 1
        
        # Generate efficient numbers by multiplying powers of 2, 3, 5
        efficient_numbers = []
        
        # Generate numbers up to a reasonable limit (min_size * 2)
        limit = min_size * 2
        
        # Powers of 2
        power_2 = 1
        while power_2 <= limit:
            # Powers of 3
            power_3 = power_2
            while power_3 <= limit:
                # Powers of 5
                power_5 = power_3
                while power_5 <= limit:
                    efficient_numbers.append(power_5)
                    power_5 *= 5
                power_3 *= 3
            power_2 *= 2
        
        # Sort and find the first number >= min_size
        efficient_numbers.sort()
        
        for num in efficient_numbers:
            if num >= min_size:
                return num
        
        # Fallback: if no efficient number found, use next power of 2
        power_of_2 = 1
        while power_of_2 < min_size:
            power_of_2 *= 2
        
        return power_of_2
    
    def _find_bilayer_pdb_with_cryst1(self) -> Optional[Path]:
        """
        Find bilayer PDB file that contains CRYST1 record for box dimensions.
        
        Prioritizes bilayer*_lipid.pdb files generated by packmol-memgen --parametrize,
        which contain essential CRYST1 box information for MD simulations.
        
        Returns:
            Path to bilayer PDB with CRYST1 record, or None if not found
        """
        # First priority: Look for final prepared files with pattern bilayer*_lipid.pdb
        final_pattern = "bilayer*_lipid.pdb"
        final_files = list(self.working_dir.glob(final_pattern))
        
        for pdb_file in final_files:
            if self._is_final_prepared_pdb(pdb_file):
                self.logger.info(f"Found final prepared bilayer PDB with CRYST1: {pdb_file}")
                return pdb_file
        
        # Second priority: Look for other bilayer files with CRYST1
        other_patterns = [
            "bilayer_*.pdb",
            "*_bilayer.pdb", 
            "*membrane*.pdb"
        ]
        
        for pattern in other_patterns:
            bilayer_files = list(self.working_dir.glob(pattern))
            for pdb_file in bilayer_files:
                # Skip lipid-only files unless it's the final prepared pattern
                if ("lipid" in pdb_file.name.lower() and 
                    not pdb_file.name.endswith("_lipid.pdb")):
                    continue
                
                # Check if this file has CRYST1 record and is not intermediate
                if self._is_final_prepared_pdb(pdb_file):
                    self.logger.info(f"Found bilayer PDB with CRYST1: {pdb_file}")
                    return pdb_file
        
        return None
    
    def _find_original_bilayer_pdb(self) -> Optional[Path]:
        """
        Find the final prepared bilayer PDB file (from packmol-memgen).
        
        The correct file should:
        1. Have pattern bilayer*_lipid.pdb (final prepared file)
        2. Contain CRYST1 header (properly formatted)
        3. NOT start with "REMARK   Packmol generated pdb file" (intermediate file)
        
        Returns:
            Path to the final prepared bilayer PDB file, or None if not found
        """
        # First priority: Look for final prepared files with pattern bilayer*_lipid.pdb
        final_pattern = "bilayer*_lipid.pdb"
        final_files = list(self.working_dir.glob(final_pattern))
        
        for pdb_file in final_files:
            # Verify this is the correct final file (has CRYST1, not intermediate)
            if self._is_final_prepared_pdb(pdb_file):
                self.logger.info(f"Found final prepared bilayer PDB: {pdb_file}")
                return pdb_file
        
        # Second priority: Look for other bilayer files with CRYST1 header
        other_patterns = [
            "bilayer_*.pdb",
            "*_bilayer.pdb"
        ]
        
        for pattern in other_patterns:
            bilayer_files = list(self.working_dir.glob(pattern))
            for pdb_file in bilayer_files:
                # Skip if it's an intermediate file or doesn't have CRYST1
                if not self._is_final_prepared_pdb(pdb_file):
                    continue
                self.logger.info(f"Found bilayer PDB with CRYST1: {pdb_file}")
                return pdb_file
        
        # Fallback: Any bilayer file (warn user about potential issues)
        for pattern in ["bilayer_*.pdb", "*_bilayer.pdb"]:
            bilayer_files = list(self.working_dir.glob(pattern))
            if bilayer_files:
                self.logger.warning(f"Using fallback bilayer PDB (may not have CRYST1): {bilayer_files[0]}")
                return bilayer_files[0]
        
        return None
    
    def _is_final_prepared_pdb(self, pdb_file: Path) -> bool:
        """
        Check if a PDB file is the final prepared file (not an intermediate).
        
        The final prepared file should:
        1. Have CRYST1 header line
        2. NOT start with intermediate file markers like "REMARK   Packmol generated pdb file"
        
        Args:
            pdb_file: Path to PDB file to check
            
        Returns:
            True if this is the final prepared file, False otherwise
        """
        try:
            with open(pdb_file, 'r') as f:
                lines = f.readlines()
            
            if not lines:
                return False
            
            has_cryst1 = False
            is_intermediate = False
            
            # Check first few lines for indicators
            for i, line in enumerate(lines[:10]):  # Check first 10 lines
                line = line.strip()
                
                # Check for CRYST1 header (good indicator of final file)
                if line.startswith('CRYST1'):
                    has_cryst1 = True
                    
                # Check for intermediate file markers (bad indicators)
                if ('Packmol generated pdb file' in line and 
                    'Packmol Memgen estimated parameters' in line):
                    is_intermediate = True
                    break
                    
                if 'charmmlipid2amber.py transformed file' in line:
                    is_intermediate = True
                    break
            
            # Final file should have CRYST1 and NOT be intermediate
            result = has_cryst1 and not is_intermediate
            
            if result:
                self.logger.debug(f"✅ {pdb_file.name} is final prepared file (has CRYST1, not intermediate)")
            else:
                if is_intermediate:
                    self.logger.debug(f"❌ {pdb_file.name} is intermediate file (packmol-memgen generated)")
                elif not has_cryst1:
                    self.logger.debug(f"❌ {pdb_file.name} missing CRYST1 header")
            
            return result
            
        except Exception as e:
            self.logger.warning(f"Could not check PDB file {pdb_file}: {e}")
            return False
    
    def _read_amber_box_dimensions(self, coord_file: Path) -> Tuple[float, float, float]:
        """
        Read box dimensions from AMBER coordinate (.inpcrd, .crd or .rst) file.
        
        Args:
            coord_file: Path to AMBER coordinate file
            
        Returns:
            Tuple of (a, b, c) box dimensions in Angstroms
        """
        try:
            with open(coord_file, 'r') as f:
                lines = f.readlines()
                
            # Box information is in the last line of AMBER coordinate files
            if len(lines) >= 2:
                last_line = lines[-1].strip()
                box_values = last_line.split()
                
                if len(box_values) >= 3:
                    a = float(box_values[0])
                    b = float(box_values[1])
                    c = float(box_values[2])
                    self.logger.info(f"Read box dimensions from AMBER file: {a:.2f} x {b:.2f} x {c:.2f} Å")
                    return a, b, c
                    
        except (FileNotFoundError, ValueError, IndexError) as e:
            self.logger.warning(f"Could not read box dimensions from AMBER file {coord_file}: {e}")
        
        # Fallback to default
        self.logger.info("Using default box dimensions for AMBER system (100x100x100 Å)")
        return 100.0, 100.0, 100.0
    
    def generate_config_file(self, stage_name: str, stage_params: Dict[str, Any], 
                           stage_index: int, system_files: Dict[str, str],
                           previous_stage_name: Optional[str] = None) -> str:
        """
        Generate NAMD configuration file for a specific equilibration stage using AMBER force field.
        
        Args:
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths (should include 'prmtop' and 'inpcrd')
            previous_stage_name: Name of the previous stage for restart files (optional)
        
        Returns:
            Content of the NAMD configuration file
        """
        
        config_lines = []
        
        # Header
        config_lines.extend([
            "#############################################################",
            f"## NAMD Configuration File for {stage_params.get('name', stage_name)}",
            f"## Generated by Gatewizard",
            f"## Stage {stage_index + 1}: {stage_params.get('description', '')}",
            f"## Force Field: AMBER",
            "#############################################################",
            ""
        ])
        
        # Input files - AMBER format (now using local copies)
        config_lines.extend([
            "# Input files - AMBER format",
            "amber              on",
            f"parmfile           system.prmtop",
            f"ambercoor          system.inpcrd"
        ])
        
        # Restart files for stages after the first
        if stage_index > 0:
            # Use passed previous stage name or build from previous stage index
            if previous_stage_name:
                prev_stage = self._generate_output_name(previous_stage_name, stage_index - 1)
            else:
                # Fallback: use the input name generation method
                prev_stage = self._generate_input_name(stage_index, previous_stage_name)
            
            config_lines.extend([
                f"bincoordinates     {prev_stage}.coor",
                f"binvelocities      {prev_stage}.vel", 
                f"extendedsystem     {prev_stage}.xsc"
            ])
        
        config_lines.append("")
        
        # AMBER Force field settings
        config_lines.extend([
            "# AMBER Force field settings",
            "exclude            scaled1-4",
            "oneFourScaling         0.833333333",  # = 1/1.2 (SCEE=1.2 in AMBER)
            "scnb               2.0",          # VDW 1-4 scaling factor
            "readexclusions     yes",         # Read exclusions from PARM file
            "switching          off",         # Turn off switching (AMBER default)
            "LJcorrection       on",          # Apply analytical tail correction
            "zeromomentum       on",          # Remove COM drift (netfrc=1 in AMBER)
            ""
        ])
        
        # Output files
        output_name = self._generate_output_name(stage_name, stage_index)
        dcd_freq = stage_params.get('dcd_freq', 5000)
        
        config_lines.extend([
            "# Output files",
            f"outputName         {output_name}",
            f"dcdfile            {output_name}.dcd",
            f"dcdfreq            {dcd_freq}",
            f"outputEnergies     {dcd_freq}",
            f"outputPressure     {dcd_freq}",
            f"outputTiming       {dcd_freq}",
            f"xstFreq            {dcd_freq}",
            ""
        ])
        
        # Basic simulation parameters
        timestep = stage_params.get('timestep', 2.0)  # Default 2 fs (NAMD uses femtoseconds)
        steps = stage_params.get('steps', 125000)
        
        config_lines.extend([
            "# Simulation parameters",
            f"timestep           {timestep}",
            f"nonbondedFreq      1",
            f"fullElectFrequency 1",            # AMBER default
            f"stepspercycle      10",
            "# Note: numsteps parameter is intentionally omitted",
            "# Steps are controlled by the 'run' command at the end",
            "",
            "# AMBER-compatible force field settings",
            "rigidBonds         all",          # SHAKE all bonds (ntc=2, ntf=2)
            "rigidTolerance     1.0e-8",       # SHAKE tolerance (tol in AMBER)
            "rigidIterations    100",
            "useSettle          on",           # Use SETTLE for water (jfastw=0)
            f"cutoff             {stage_params.get('cutoff', 9.0)}",  # Default AMBER cutoff
            f"pairlistdist       {stage_params.get('cutoff', 9.0) + 2.0}",  # cutoff + 2.0
            "watermodel         tip3",         # Default water model
            ""
        ])
        
        # Temperature control - Langevin thermostat (corresponds to ntt=3 in AMBER)
        temperature = stage_params.get('temperature', 310.15)
        ensemble = stage_params.get('ensemble', 'NPT')
        
        if ensemble in ['NVT', 'NPT', 'NPAT', 'NPgT']:
            damping = stage_params.get('langevin_damping', 5.0)  # gamma_ln in AMBER
            
            # For first stage, set initial temperature
            # For subsequent stages, only set thermostat parameters (velocities come from restart)
            if stage_index == 0:
                config_lines.extend([
                    "# Temperature control - Langevin thermostat (first stage)",
                    f"temperature        {temperature}",      # tempi in AMBER
                    "langevin           on",                  # ntt=3 in AMBER
                    f"langevinTemp       {temperature}",      # temp0 in AMBER
                    f"langevinDamping    {damping}",          # gamma_ln in AMBER (ps^-1)
                    "langevinHydrogen   off",                 # AMBER default
                    ""
                ])
            else:
                config_lines.extend([
                    "# Temperature control - Langevin thermostat (restart stage)",
                    "langevin           on",                  # ntt=3 in AMBER
                    f"langevinTemp       {temperature}",      # temp0 in AMBER
                    f"langevinDamping    {damping}",          # gamma_ln in AMBER (ps^-1)
                    "langevinHydrogen   off",                 # AMBER default
                    ""
                ])
        
        # Pressure control - Berendsen barostat (corresponds to ntp=1 in AMBER)
        if ensemble in ['NPT', 'NPAT', 'NPgT']:
            pressure = stage_params.get('pressure', 1.01325)  # AMBER default pressure
            surface_tension = stage_params.get('surface_tension', 0.0)
            
            # Use Berendsen pressure control for AMBER compatibility
            config_lines.extend([
                "# Pressure control - Berendsen barostat",
                "BerendsenPressure     on",               # ntp=1 in AMBER
                f"BerendsenPressureTarget {pressure}",    # pres0 in AMBER
                "BerendsenPressureCompressibility  4.57e-5",  # compressibility in AMBER (1/bar)
                "BerendsenPressureRelaxationTime 100.0", # taup in AMBER (fs in NAMD, ps in AMBER)
                "useGroupPressure      yes",             # needed for rigidBonds
            ])
            
            # Configure pressure scaling based on ensemble
            if ensemble == 'NPAT':
                # Semi-isotropic pressure control for membrane systems
                config_lines.extend([
                    "useFlexibleCell       yes",             # allow cell shape changes
                    "useConstantRatio      yes",             # keep XY ratio constant
                ])
                
                # Add surface tension if specified
                if surface_tension > 0.0:
                    config_lines.extend([
                        f"# Surface tension control for NPAT ensemble",
                        f"surfaceTensionTarget  {surface_tension}",  # dyn/cm
                    ])
                    
            elif ensemble == 'NPgT':
                # Constant surface tension ensemble
                config_lines.extend([
                    "useFlexibleCell       yes",             # allow cell shape changes
                    "useConstantRatio      yes",             # keep XY ratio constant
                ])
                
                # Surface tension is required for NPgT
                if surface_tension == 0.0:
                    surface_tension = 0.0  # Default value in dyn/cm
                    self.logger.warning(f"NPgT ensemble requires surface tension. Using default: {surface_tension} dyn/cm")
                
                config_lines.extend([
                    f"# Surface tension control for NPgT ensemble",
                    f"surfaceTensionTarget  {surface_tension}",  # dyn/cm
                ])
                
            else:  # NPT
                # Isotropic pressure control
                config_lines.extend([
                    "useFlexibleCell       no",              # isotropic scaling
                    "useConstantArea       no",
                ])
            
            config_lines.append("")
        
        # Restraints/Constraints (if needed)
        constraints = stage_params.get('constraints', {})
        has_restraints = any(float(v) > 0 for v in constraints.values())
        
        if has_restraints:
            # Calculate constraint scaling based on the maximum restraint force
            max_force = max(float(v) for v in constraints.values() if float(v) > 0)
            constraint_scaling = min(10.0, max(1.0, max_force))  # Scale between 1.0 and 10.0
            
            # Use stage-specific restraints file if available, otherwise general file
            config_name = self._get_config_name(stage_name, stage_index)
            if config_name == "step7_production":
                stage_restraints_file = f"restraints/{config_name}_restraints.pdb"
            else:
                stage_restraints_file = f"restraints/{config_name}_equilibration_restraints.pdb"
            general_restraints_file = "restraints.pdb"
            
            config_lines.extend([
                "# Harmonic restraints",
                "constraints        on",
                "consexp            2",                    # harmonic restraints
                f"consref            {stage_restraints_file}",  # reference coordinates
                f"conskfile          {stage_restraints_file}",  # force constant file
                "conskcol           B",                    # use B-factor column
                f"constraintScaling  {constraint_scaling}",  # overall scaling factor
                f"# Note: If {stage_restraints_file} not found, use {general_restraints_file}",
                ""
            ])
        
        # Non-bonded interactions - AMBER settings are already included above
        # The cutoff, pairlistdist, switching, exclude, etc. are set in the simulation parameters
        
        # Get box dimensions for PME grid size calculation
        # PRIORITY: Use bilayer_pdb from system_files for CRYST1 record (most accurate)
        bilayer_pdb_path = self.working_dir / system_files.get('bilayer_pdb', '')
        inpcrd_file = self.working_dir / system_files.get('inpcrd', 'system.inpcrd')
        
        # Try to read from bilayer PDB with CRYST1 first (highest priority)
        if bilayer_pdb_path.exists():
            box_a, box_b, box_c = self._read_box_dimensions(bilayer_pdb_path)
            self.logger.info(f"Using box dimensions from bilayer PDB for cell basis: {bilayer_pdb_path.name}")
        elif inpcrd_file.exists():
            box_a, box_b, box_c = self._read_amber_box_dimensions(inpcrd_file)
            self.logger.info(f"Using box dimensions from AMBER inpcrd file: {inpcrd_file.name}")
        else:
            self.logger.warning("No coordinate file found, using default box dimensions")
            box_a, box_b, box_c = 100.0, 100.0, 100.0
        
        # PME electrostatics - AMBER compatible settings
        config_lines.extend([
            "# PME electrostatics - AMBER compatible",
            "PME                yes",
            "PMETolerance       1.0e-6",      # dsum_tol in AMBER
            "PMEInterpOrder     4",           # order=4 in AMBER (cubic spline)
            "PMEGridSpacing     1.0",         # Let NAMD automatically calculate grid sizes
            ""
        ])
        
        # Periodic boundary conditions
        # Box dimensions already calculated above for PME grid sizing
        config_lines.extend([
            "# Periodic boundary conditions",
            f"# Box dimensions: {box_a:.2f} x {box_b:.2f} x {box_c:.2f} Å",
            f"cellBasisVector1   {box_a:.6f}   0.000000   0.000000",
            f"cellBasisVector2   0.000000   {box_b:.6f}   0.000000", 
            f"cellBasisVector3   0.000000   0.000000   {box_c:.6f}",
            "cellOrigin         0.0   0.0   0.0",
            ""
        ])
        
        # Wrap output
        config_lines.extend([
            "# Wrap output",
            "wrapWater          on",
            "wrapAll            on",
            ""
        ])
        
        # Minimization specific settings
        integrator = stage_params.get('integrator', '')
        if 'minimization' in stage_name.lower() or integrator == 'conjugate_gradient':
            minimize_steps = min(5000, steps)
            config_lines.extend([
                "# Energy minimization",
                f"minimize           {minimize_steps}",
                ""
            ])
        
        # Run command
        if 'minimization' not in stage_name.lower():
            config_lines.extend([
                "# Run the simulation", 
                f"run               {steps}"
            ])
        
        return "\n".join(config_lines)
    
    def generate_restraints_file(self, system_pdb: Path, constraints: Dict[str, float], 
                               output_file: Path, stage_name: str = "") -> None:
        """
        Generate restraints PDB file with B-factors for constraint forces.
        
        Uses the final system.pdb file for generating restraints to ensure
        consistency with the parametrized system used in simulations.
        
        Args:
            system_pdb: Path to the system.pdb file  
            constraints: Dictionary of constraint types and forces
            output_file: Path for output restraints file
            stage_name: Name of the equilibration stage (for documentation)
        """
        
        # Use only the final system.pdb file for restraints
        if not system_pdb.exists():
            self.logger.error(f"System PDB file not found: {system_pdb}")
            return
            
        self.logger.info(f"Using system.pdb for restraints: {system_pdb.name}")
        
        # Read system PDB file
        with open(system_pdb, 'r') as f:
            lines = f.readlines()
        
        # Process each line and assign restraint forces
        restraint_lines = []
        atom_count = 0
        restraint_stats = {
            'protein_backbone': 0,
            'protein_sidechain': 0,
            'lipid_head': 0,
            'lipid_tail': 0,
            'water': 0,
            'ions': 0,
            'other': 0
        }
        
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                # Parse atom information
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                chain_id = line[21].strip()
                
                # Determine restraint force based on atom type
                restraint_force, atom_type = self._get_restraint_force_detailed(
                    atom_name, residue_name, chain_id, constraints
                )
                
                # Update statistics
                if atom_type in restraint_stats:
                    restraint_stats[atom_type] += 1
                
                # Replace B-factor with restraint force (columns 61-66)
                new_line = line[:60] + f"{restraint_force:6.2f}" + line[66:]
                restraint_lines.append(new_line)
                atom_count += 1
            else:
                # Keep non-atom lines as is
                restraint_lines.append(line)
        
        # Write restraints file
        with open(output_file, 'w') as f:
            f.writelines(restraint_lines)
        
        # Log statistics
        self.logger.info(f"Generated restraints file: {output_file}")
        self.logger.info(f"Source PDB: {system_pdb.name}")
        self.logger.info(f"Stage: {stage_name}")
        self.logger.info(f"Total atoms processed: {atom_count}")
        for atom_type, count in restraint_stats.items():
            if count > 0:
                force = constraints.get(atom_type, 0.0)
                self.logger.info(f"  {atom_type}: {count} atoms, force = {force} kcal/mol/Å²")
    
    def _get_restraint_force_detailed(self, atom_name: str, residue_name: str, 
                                    chain_id: str, constraints: Dict[str, float]) -> tuple:
        """
        Determine the appropriate restraint force for an atom with detailed classification.
        
        Args:
            atom_name: Name of the atom
            residue_name: Name of the residue
            chain_id: Chain identifier
            constraints: Dictionary of constraint forces
        
        Returns:
            Tuple of (restraint_force, atom_type)
        """
        
        # Standard protein residues (including protonation states from AMBER)
        protein_residues = {
            # Standard amino acids
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
            'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
            'THR', 'TRP', 'TYR', 'VAL', 'HSE', 'HSD', 'HSP', 'CYX',
            # Protonation states from AMBER/propka (based on PROTONATION_STATES dict)
            'ASH',  # Protonated aspartic acid
            'GLH',  # Protonated glutamic acid
            'HIE',  # Histidine with proton on epsilon nitrogen
            'HID',  # Histidine with proton on delta nitrogen
            'HIP',  # Histidine with both nitrogens protonated
            'LYN',  # Neutral lysine (deprotonated)
            'TYM',  # Deprotonated tyrosine
            'CYM',  # Deprotonated cysteine
            'SEP',  # Phosphorylated serine
            'T2P',  # Phosphorylated threonine
            # Terminal caps
            'ACE', 'NHE', 'NME', 'COO'
        }
        
        # Lipid residues (include common AMBER/CHARMM lipid names)
        lipid_residues = {
            # Standard lipid names
            'POPC', 'POPE', 'POPS', 'DPPC', 'DMPC', 'DOPC', 'DSPC',
            'CHOL', 'CHOLEST', 'PALM', 'OLEO', 'STEROL',
            # AMBER-style lipid names (common in packmol-memgen)
            'PC', 'PE', 'PS', 'PA', 'PG', 'PI', 'SM', 'CHL', 'CHOL',
            'OL', 'LA', 'MY', 'PA', 'ST', 'AR',  # AMBER lipid residue codes
            # Additional lipid variants
            'OLE', 'PAL', 'STE', 'LIN'
            # Note: LYN removed - it's neutral lysine (protein), not a lipid
        }
        
        # Water residues (include various naming conventions)
        water_residues = {'TIP3', 'HOH', 'WAT', 'SOL', 'TIP4', 'SPC', 'T3P', 'T4P'}
        
        # Ion residues (include various naming conventions)
        ion_residues = {
            'NA', 'CL', 'K', 'CA', 'MG', 'ZN', 'FE', 'CU',
            'Na+', 'Cl-', 'K+', 'Ca2+', 'Mg2+', 'Zn2+', 'Fe2+', 'Fe3+',
            'SOD', 'CLA', 'POT', 'CAL', 'MAG', 'ZIN', 'IRN', 'COP'
        }
        
        # Protein atoms
        if residue_name in protein_residues:
            # Backbone atoms (including hydrogens)
            backbone_atoms = {
                'N', 'CA', 'C', 'O', 'OXT',
                'H', 'HN', 'HA', 'HT1', 'HT2', 'HT3'
            }
            
            if atom_name in backbone_atoms:
                return constraints.get('protein_backbone', 0.0), 'protein_backbone'
            else:
                return constraints.get('protein_sidechain', 0.0), 'protein_sidechain'
        
        # Lipid molecules
        elif residue_name in lipid_residues:
            # Head group atoms (phosphate, choline, ethanolamine, etc.)
            head_atoms = {
                # Phosphate groups
                'P', 'O11', 'O12', 'O13', 'O14', 'O21', 'O22', 'O31', 'O32', 'O33', 'O34',
                'O1P', 'O2P', 'O3P', 'O4P', 'OP1', 'OP2', 'OP3', 'OP4',
                # Choline and ethanolamine heads (specific to head group only)
                'N', 'C11', 'C12', 'C13', 'C14',  # Removed C15, C16 as they can be tail carbons
                'N31', 'C32', 'C33', 'C34', 'C35',
                # Glycerol backbone (connects head to tails)
                'C1', 'C2', 'C3', 'O21', 'O31',
                # Common head group patterns
                'HN1', 'HN2', 'HN3', 'HO2', 'HO3', 'HS'
            }
            
            # Check atom name patterns for head vs tail classification
            if (atom_name in head_atoms or 
                atom_name.startswith(('P', 'O1', 'O2', 'N')) or
                'P' in atom_name or 'N3' in atom_name or
                (atom_name.startswith('C') and len(atom_name) <= 2 and atom_name in ['C1', 'C2', 'C3'])):
                return constraints.get('lipid_head', 0.0), 'lipid_head'
            else:
                # Tail carbons and other atoms (including C15, C16, etc.)
                return constraints.get('lipid_tail', 0.0), 'lipid_tail'
        
        # Water molecules
        elif residue_name in water_residues:
            return constraints.get('water', 0.0), 'water'
        
        # Ions
        elif residue_name in ion_residues:
            return constraints.get('ions', 0.0), 'ions'
        
        # Other molecules (ligands, cofactors, etc.)
        else:
            return constraints.get('other', 0.0), 'other'
    
    def setup_namd_equilibration(
        self,
        system_files: Optional[Dict[str, str]] = None,
        stage_params_list: Optional[List[Dict[str, Any]]] = None,
        output_name: str = "equilibration",
        scheme_type: Optional[str] = None,
        namd_executable: str = "namd3"
    ) -> Dict[str, Any]:
        """
        Complete NAMD equilibration setup - replicates GUI workflow.
        
        This simplified method performs all steps needed for equilibration:
        1. Auto-detects system files (if not provided)
        2. Creates output directory structure (equilibration/namd/)
        3. Copies system files to NAMD directory
        4. Generates NAMD configuration files for all stages
        5. Generates restraint files for each stage
        6. Creates run script
        7. Creates protocol summary
        
        Args:
            system_files: Dictionary with source file paths (optional - will auto-detect if None):
                - 'prmtop': AMBER topology file (.prmtop)
                - 'inpcrd': AMBER coordinate file (.inpcrd)
                - 'pdb': System PDB file
                - 'bilayer_pdb': (REQUIRED) Bilayer PDB file with CRYST1 record for box dimensions
                If None, will search working_dir for standard file names.
            stage_params_list: List of stage dictionaries with equilibration parameters.
                Each stage must include 'ensemble' key (NVT, NPT, NPAT, or NPgT).
            output_name: Output directory name (default: "equilibration")
            scheme_type: Equilibration scheme (optional - auto-detected from stages).
                If None, will be extracted from the 'ensemble' field of the first stage.
                Can be explicitly set if needed (NVT, NPT, NPAT, or NPgT).
            namd_executable: NAMD executable path (default: "namd3")
        
        Returns:
            Dictionary with paths to generated files:
            {
                'output_dir': Path to main output directory,
                'namd_dir': Path to NAMD directory,
                'config_files': List of config file paths,
                'restraints_dir': Path to restraints directory,
                'run_script': Path to run script,
                'summary_file': Path to protocol summary
            }
        
        Example 1 (Auto-detect files):
            >>> manager = NAMDEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_namd_equilibration(
            ...     stage_params_list=[
            ...         {'name': 'Equilibration 1', 'time_ns': 0.125, ...}
            ...     ],
            ...     scheme_type="NPT"
            ... )
        
        Example 2 (Explicit file paths):
            >>> manager = NAMDEquilibrationManager(Path("/work/dir"))
            >>> result = manager.setup_namd_equilibration(
            ...     system_files={
            ...         'prmtop': '/path/to/system.prmtop',
            ...         'inpcrd': '/path/to/system.inpcrd',
            ...         'pdb': '/path/to/system.pdb',
            ...         'bilayer_pdb': '/path/to/bilayer_lipid.pdb'
            ...     },
            ...     stage_params_list=[
            ...         {'name': 'Equilibration 1', 'time_ns': 0.125, ...}
            ...     ],
            ...     scheme_type="NPT"
            ... )
            >>> # Run: cd result['namd_dir'] && ./run_equilibration.sh
        """
        import shutil
        import json
        
        self.logger.info("=== Setting up NAMD equilibration ===")
        
        # Auto-detect scheme_type from stages if not provided
        if scheme_type is None:
            if stage_params_list and len(stage_params_list) > 0:
                # Extract ensemble from first stage
                first_ensemble = stage_params_list[0].get('ensemble', 'NPT')
                scheme_type = first_ensemble
                self.logger.info(f"Auto-detected scheme_type from stages: {scheme_type}")
            else:
                # Default fallback
                scheme_type = "NPT"
                self.logger.info(f"No stages provided, using default scheme_type: {scheme_type}")
        
        # Validate scheme_type
        valid_schemes = ['NVT', 'NPT', 'NPAT', 'NPgT']
        if scheme_type not in valid_schemes:
            raise ValueError(f"Invalid scheme_type '{scheme_type}'. Must be one of {valid_schemes}")
        
        # Auto-detect system files if not provided
        if system_files is None:
            self.logger.info("Auto-detecting system files in working directory...")
            system_files = self.find_system_files()
            if system_files is None:
                raise FileNotFoundError(
                    "Could not auto-detect required system files. "
                    "Please provide system_files dictionary explicitly."
                )
        
        # Validate required files exist
        for file_type in ['prmtop', 'inpcrd', 'pdb', 'bilayer_pdb']:
            if file_type not in system_files:
                raise ValueError(f"Missing required key '{file_type}' in system_files")
            file_path = Path(system_files[file_type])
            if not file_path.exists():
                raise FileNotFoundError(f"File not found: {file_path}")
        
        # Create output directory structure
        output_dir = Path(self.working_dir) / output_name
        namd_dir = output_dir / "namd"
        restraints_dir = namd_dir / "restraints"
        
        for directory in [output_dir, namd_dir, restraints_dir]:
            directory.mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Created directory: {directory}")
        
        # Step 1: Copy system files to NAMD directory
        self.logger.info("Copying system files to NAMD directory...")
        copied_files = {}
        
        for file_type, source_path in system_files.items():
            source_path_obj = Path(source_path)
            
            # Special handling for bilayer_pdb - keep original name, only used for CRYST1
            if file_type == 'bilayer_pdb':
                target_filename = source_path_obj.name  # Keep original filename
                target_path = namd_dir / target_filename
                shutil.copy2(source_path, target_path)
                copied_files[file_type] = target_filename
                self.logger.info(f"  Copied {file_type}: {source_path_obj.name} (for CRYST1 box info only)")
            else:
                # Standard naming: system.prmtop, system.inpcrd, system.pdb
                target_filename = f"system{source_path_obj.suffix}"
                target_path = namd_dir / target_filename
                shutil.copy2(source_path, target_path)
                copied_files[file_type] = target_filename
                self.logger.info(f"  Copied {file_type}: {source_path_obj.name} -> {target_filename}")
        
        # Step 2: Generate NAMD configuration files for each stage
        self.logger.info("Generating NAMD configuration files...")
        config_files = []
        protocols_dict = {}
        
        # Ensure stage_params_list is not None
        if stage_params_list is None:
            raise ValueError("stage_params_list cannot be None")
        
        # Convert stage list to protocols dictionary
        for i, stage_params in enumerate(stage_params_list):
            stage_name = stage_params.get('name', f'Equilibration {i+1}')
            protocols_dict[stage_name] = stage_params
        
        previous_stage_name = None
        for i, (stage_name, stage_params) in enumerate(protocols_dict.items()):
            # Generate config using CHARMM-GUI template system
            config_content = self.generate_charmm_gui_config_file(
                stage_name=stage_name,
                stage_params=stage_params,
                stage_index=i,
                system_files=copied_files,  # Use relative names
                scheme_type=scheme_type,
                previous_stage_name=previous_stage_name,
                all_stage_settings=protocols_dict
            )
            
            # Write configuration file
            config_name = self._get_config_name(stage_name, i)
            if config_name == "step7_production":
                config_file = namd_dir / f"{config_name}.conf"
            else:
                config_file = namd_dir / f"{config_name}_equilibration.conf"
            
            config_file.write_text(config_content)
            config_files.append(config_file)
            self.logger.info(f"  Generated: {config_file.name}")
            
            previous_stage_name = stage_name
        
        # Step 3: Generate restraint files for each stage
        self.logger.info("Generating restraint files...")
        system_pdb = namd_dir / copied_files.get('pdb', 'system.pdb')
        
        if system_pdb.exists():
            for i, (stage_name, stage_params) in enumerate(protocols_dict.items()):
                constraints = stage_params.get('constraints', {})
                has_restraints = any(float(v) > 0 for v in constraints.values())
                
                if has_restraints:
                    config_name = self._get_config_name(stage_name, i)
                    if config_name == "step7_production":
                        restraint_file = restraints_dir / f"{config_name}_restraints.pdb"
                    else:
                        restraint_file = restraints_dir / f"{config_name}_equilibration_restraints.pdb"
                    
                    self.generate_restraints_file(
                        system_pdb=system_pdb,
                        constraints=constraints,
                        output_file=restraint_file,
                        stage_name=stage_params.get('name', stage_name)
                    )
                    self.logger.info(f"  Generated: {restraint_file.name}")
        else:
            self.logger.warning(f"System PDB not found: {system_pdb}, skipping restraints")
        
        # Step 4: Generate run script
        self.logger.info("Generating run script...")
        run_script_content = self.generate_run_script(protocols_dict, namd_executable)
        run_script = namd_dir / "run_equilibration.sh"
        run_script.write_text(run_script_content)
        run_script.chmod(0o755)
        self.logger.info(f"  Generated: {run_script.name}")
        
        # Step 5: Create protocol summary
        protocol_summary = {
            "protocol_name": f"{scheme_type} Equilibration Protocol",
            "total_stages": len(protocols_dict),
            "scheme_type": scheme_type,
            "stages": protocols_dict,
            "namd_executable": namd_executable,
            "force_field": "AMBER"
        }
        
        summary_file = namd_dir / "protocol_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(protocol_summary, f, indent=2)
        self.logger.info(f"  Generated: {summary_file.name}")
        
        self.logger.info("=== Setup complete ===")
        self.logger.info(f"Output directory: {namd_dir}")
        self.logger.info(f"To run: cd {namd_dir} && ./run_equilibration.sh")
        
        return {
            'output_dir': output_dir,
            'namd_dir': namd_dir,
            'config_files': config_files,
            'restraints_dir': restraints_dir,
            'run_script': run_script,
            'summary_file': summary_file
        }
    
    def generate_colvars_file(self, system_pdb: Path, output_file: Path, 
                            stage_params: Optional[Dict[str, Any]] = None) -> None:
        """
        Generate NAMD colvars configuration file for bilayer thickness restraint using bilayer_utils.
        
        Args:
            system_pdb: Path to the system PDB file to analyze for phosphate atoms
            output_file: Path for output colvars configuration file (should be bilayer_thickness.col)
            stage_params: Stage parameters including bilayer thickness and force constant
        """
        try:
            if stage_params and stage_params.get('bilayer_thickness') is not None:
                self._generate_bilayer_thickness_restraint(system_pdb, output_file, stage_params)
            else:
                self.logger.info("Bilayer thickness restraint not enabled for this stage")
                
        except Exception as e:
            self.logger.error(f"Error generating colvars file: {e}")
            self.logger.info("Falling back to basic approach")
            # Ensure stage_params is not None for fallback
            if stage_params is None:
                stage_params = {}
            self._generate_bilayer_thickness_restraint_fallback(system_pdb, output_file, stage_params)
    
    def _generate_bilayer_thickness_restraint(self, system_pdb: Path, output_file: Path, stage_params: Optional[Dict[str, Any]]) -> None:
        """
        Generate colvars configuration for bilayer thickness harmonic restraint using bilayer_utils.
        This creates a simple harmonic restraint to maintain bilayer thickness.
        """
        try:
            # Import bilayer utilities
            import sys
            import os
            # Add the gatewizard utils path to sys.path
            gatewizard_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
            utils_path = os.path.join(gatewizard_root, 'utils')
            if utils_path not in sys.path:
                sys.path.insert(0, utils_path)
            
            from bilayer_utils import BilayerAnalyzer  # type: ignore
            
            # Initialize bilayer analyzer with all phosphate patterns (including P31)
            analyzer = BilayerAnalyzer()
            
            # Analyze bilayer to find ALL phosphate atoms
            upper_bilayer, lower_bilayer = analyzer.analyze_bilayer_from_pdb(str(system_pdb))
            
            if not upper_bilayer or not lower_bilayer:
                self.logger.warning("Could not find phosphate atoms in both bilayers for thickness restraint")
                return
            
            # Get statistics for logging
            stats = analyzer.get_bilayer_statistics(upper_bilayer, lower_bilayer)
            self.logger.info(f"Found {stats['total_phosphorus_atoms']} total phosphorus atoms for bilayer thickness restraint")
            self.logger.info(f"Upper bilayer: {stats['upper_bilayer_count']} atoms")
            self.logger.info(f"Lower bilayer: {stats['lower_bilayer_count']} atoms") 
            self.logger.info(f"Atom types found: {', '.join(stats['atom_types'])}")
            self.logger.info(f"Residue types: {', '.join(stats['residue_types'])}")
            
            # Get parameters from stage_params or use defaults
            if stage_params is None:
                stage_params = {}
            target_thickness = float(stage_params.get('bilayer_thickness', 39.1))  # Å
            force_constant = float(stage_params.get('force_constant', 10.0))  # kcal/mol/Å²
            
            # Generate simple harmonic restraint configuration
            colvars_content = self._create_harmonic_thickness_config(upper_bilayer, lower_bilayer, target_thickness, force_constant)
            
            # Write colvars file
            with open(output_file, 'w') as f:
                f.write(colvars_content)
            
            self.logger.info(f"Generated bilayer thickness restraint file with {stats['total_phosphorus_atoms']} phosphate atoms: {output_file}")
            self.logger.info(f"Target thickness: {target_thickness} Å, Force constant: {force_constant} kcal/mol/Å²")
            
        except ImportError as e:
            self.logger.error(f"Could not import bilayer_utils: {e}")
            self.logger.info("Falling back to basic approach")
            if stage_params is None:
                stage_params = {}
            self._generate_bilayer_thickness_restraint_fallback(system_pdb, output_file, stage_params)
        except Exception as e:
            self.logger.error(f"Error in bilayer analysis: {e}")
            self.logger.info("Falling back to basic approach")
            if stage_params is None:
                stage_params = {}
            self._generate_bilayer_thickness_restraint_fallback(system_pdb, output_file, stage_params)
    
    def _create_harmonic_thickness_config(self, upper_bilayer, lower_bilayer, target_thickness: float, force_constant: float) -> str:
        """
        Create simple harmonic restraint colvars configuration for bilayer thickness.
        
        Args:
            upper_bilayer: List of PhosphorusAtom objects for upper bilayer
            lower_bilayer: List of PhosphorusAtom objects for lower bilayer
            target_thickness: Target bilayer thickness in Å
            force_constant: Force constant in kcal/mol/Å²
            
        Returns:
            Complete colvars configuration string with harmonic restraint
        """
        # Extract ALL NAMD indices (0-based) - no limit, include all atoms
        upper_indices = [str(atom.namd_index) for atom in upper_bilayer]
        lower_indices = [str(atom.namd_index) for atom in lower_bilayer]
        
        # Format indices in readable chunks (10 per line for readability)
        def format_indices_multiline(indices: List[str], indent: str = "            ") -> str:
            if not indices:
                return ""
            
            lines = []
            for i in range(0, len(indices), 10):
                chunk = indices[i:i+10]
                lines.append(indent + " ".join(chunk))
            return "\n".join(lines)
        
        upper_formatted = format_indices_multiline(upper_indices)
        lower_formatted = format_indices_multiline(lower_indices)
        
        # Create simple harmonic restraint configuration
        config = f"""# NAMD Colvars Configuration for Bilayer Thickness Restraint
# Generated by Gatewizard using bilayer analysis
# Upper bilayer atoms: {len(upper_bilayer)} phosphates
# Lower bilayer atoms: {len(lower_bilayer)} phosphates
# Total atoms used: {len(upper_bilayer) + len(lower_bilayer)} phosphates

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNumbers {' '.join(upper_indices)}
        }}
        group2 {{
            atomNumbers {' '.join(lower_indices)}
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant {force_constant}  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
        
        return config
    
    def _generate_bilayer_thickness_restraint_fallback(self, system_pdb: Path, output_file: Path, stage_params: Dict[str, Any]) -> None:
        """
        Fallback method using basic phosphate detection for bilayer thickness restraint.
        Used when bilayer_utils is not available.
        """
        self.logger.warning("Using fallback phosphate detection method for bilayer thickness restraint")
        
        # Find phosphate atoms in the system using basic method
        phosphate_atoms = self._find_phosphate_atoms_basic(system_pdb)
        
        if len(phosphate_atoms) < 2:
            self.logger.warning("Insufficient phosphate atoms found for bilayer thickness restraint")
            return
        
        # Separate into upper and lower leaflets (simple z-coordinate based)
        z_coords = [atom['z'] for atom in phosphate_atoms]
        z_center = sum(z_coords) / len(z_coords)
        
        upper_leaflet = [atom for atom in phosphate_atoms if atom['z'] > z_center]
        lower_leaflet = [atom for atom in phosphate_atoms if atom['z'] < z_center]
        
        if len(upper_leaflet) == 0 or len(lower_leaflet) == 0:
            self.logger.warning("Could not separate phosphates into leaflets")
            return
        
        # Get parameters from stage_params or use defaults
        target_thickness = float(stage_params.get('bilayer_thickness', 39.1)) if stage_params else 39.1  # Å
        force_constant = float(stage_params.get('force_constant', 10.0)) if stage_params else 10.0  # kcal/mol/Å²
        
        # Generate basic harmonic restraint configuration
        colvars_content = self._create_basic_harmonic_thickness_config(upper_leaflet, lower_leaflet, target_thickness, force_constant)
        
        # Write colvars file
        with open(output_file, 'w') as f:
            f.write(colvars_content)
        
        self.logger.info(f"Generated fallback bilayer thickness restraint file: {output_file}")
        self.logger.info(f"Found {len(upper_leaflet)} upper and {len(lower_leaflet)} lower leaflet phosphates")
        self.logger.info(f"Target thickness: {target_thickness} Å, Force constant: {force_constant} kcal/mol/Å²")
    
    def _generate_phosphate_distance_colvars_fallback(self, system_pdb: Path, output_file: Path) -> None:
        """
        Fallback method using the original simple phosphate detection.
        Used when bilayer_utils is not available.
        """
        self.logger.warning("Using fallback phosphate detection method")
        
        # Find phosphate atoms in the system using basic method
        phosphate_atoms = self._find_phosphate_atoms_basic(system_pdb)
        
        if len(phosphate_atoms) < 2:
            self.logger.warning("Insufficient phosphate atoms found for distance measurement")
            return
        
        # Separate into upper and lower leaflets (simple z-coordinate based)
        z_coords = [atom['z'] for atom in phosphate_atoms]
        z_center = sum(z_coords) / len(z_coords)
        
        upper_leaflet = [atom for atom in phosphate_atoms if atom['z'] > z_center]
        lower_leaflet = [atom for atom in phosphate_atoms if atom['z'] < z_center]
        
        if len(upper_leaflet) == 0 or len(lower_leaflet) == 0:
            self.logger.warning("Could not separate phosphates into leaflets")
            return
        
        # Generate basic colvars configuration
        colvars_content = self._create_basic_phosphate_distance_config(upper_leaflet, lower_leaflet)
        
        # Write colvars file
        with open(output_file, 'w') as f:
            f.write(colvars_content)
        
        self.logger.info(f"Generated fallback colvars file: {output_file}")
        self.logger.info(f"Found {len(upper_leaflet)} upper and {len(lower_leaflet)} lower leaflet phosphates")
    
    def _find_phosphate_atoms_basic(self, system_pdb: Path) -> List[Dict]:
        """Find phosphate atoms (P atoms) in lipid molecules."""
        phosphate_atoms = []
        
        with open(system_pdb, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith(('ATOM', 'HETATM')):
                    atom_name = line[12:16].strip()
                    residue_name = line[17:20].strip()
                    
                    # Look for phosphorus atoms in lipid residues (3-character names)
                    # Support both 'P' and 'P31' atom names
                    if atom_name in ['P', 'P31'] and residue_name in ['PC', 'PA', 'PE', 'PS', 'PG', 'PI']:
                        try:
                            atom_id = int(line[6:11].strip())
                            x = float(line[30:38].strip())
                            y = float(line[38:46].strip())
                            z = float(line[46:54].strip())
                            
                            phosphate_atoms.append({
                                'id': atom_id,
                                'name': atom_name,
                                'residue': residue_name,
                                'x': x, 'y': y, 'z': z,
                                'line_num': line_num
                            })
                        except (ValueError, IndexError) as e:
                            self.logger.warning(f"Error parsing line {line_num}: {e}")
        
        return phosphate_atoms
    
    def _create_basic_phosphate_distance_config(self, upper_leaflet: List[Dict], 
                                        lower_leaflet: List[Dict]) -> str:
        """Create basic colvars configuration for phosphate distance measurement (fallback method)."""
        
        # Use ALL atoms from each leaflet (no artificial limit)
        upper_atoms = " ".join([str(atom['id']) for atom in upper_leaflet])
        lower_atoms = " ".join([str(atom['id']) for atom in lower_leaflet])
        
        colvars_config = f"""# NAMD Colvars Configuration for Phosphate Distance Measurement (Fallback)
# Generated for ABF simulation of bilayer thickness
# Upper leaflet atoms: {len(upper_leaflet)} phosphates
# Lower leaflet atoms: {len(lower_leaflet)} phosphates

colvarsTrajFrequency    1000
colvarsRestartFrequency 1000

# Define collective variable for phosphate-phosphate distance
colvar {{
    name phosphate_distance
    
    # Distance between center of mass of upper and lower leaflet phosphates
    distance {{
        group1 {{
            atomNumbers {upper_atoms}
        }}
        group2 {{
            atomNumbers {lower_atoms}
        }}
    }}
}}

# ABF bias for phosphate distance
abf {{
    colvars          phosphate_distance
    fullSamples      200
    historyFreq      1000
    inputPrefix      ""
    outputPrefix     "phosphate_distances"
}}

# Metadynamics hills for enhanced sampling (optional)
metadynamics {{
    colvars         phosphate_distance
    hillWeight      0.1
    hillWidth       0.5
    newHillFreq     1000
    writeHillsFreq  1000
    outputPrefix    "phosphate_hills"
}}
"""
        return colvars_config
    
    def _create_basic_harmonic_thickness_config(self, upper_leaflet: List[Dict], 
                                              lower_leaflet: List[Dict], 
                                              target_thickness: float, 
                                              force_constant: float) -> str:
        """Create basic harmonic restraint configuration for bilayer thickness (fallback method)."""
        
        # Use ALL atoms from each leaflet (no artificial limit)
        upper_atoms = " ".join([str(atom['id']) for atom in upper_leaflet])
        lower_atoms = " ".join([str(atom['id']) for atom in lower_leaflet])
        
        colvars_config = f"""# NAMD Colvars Configuration for Bilayer Thickness Restraint (Fallback)
# Generated for harmonic restraint of bilayer thickness
# Upper leaflet atoms: {len(upper_leaflet)} phosphates
# Lower leaflet atoms: {len(lower_leaflet)} phosphates

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNumbers {upper_atoms}
        }}
        group2 {{
            atomNumbers {lower_atoms}
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant {force_constant}  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
        return colvars_config
    
    def generate_run_script(self, protocols: Dict[str, Dict], 
                          namd_executable: Optional[str] = None) -> str:
        """
        Generate bash script to run all equilibration stages.
        Each stage uses its own CPU/GPU settings from the protocol configuration.
        
        Args:
            protocols: Dictionary of equilibration protocols with per-stage resource settings
            namd_executable: Path to NAMD executable
        
        Returns:
            Content of the run script
        """
        
        namd_exe = namd_executable or self.namd_executable
        
        script_lines = [
            "#!/bin/bash",
            "#############################################################",
            "## NAMD Equilibration Run Script",
            "## Generated by Gatewizard",
            "#############################################################",
            "",
            "# Set NAMD executable",
            f'NAMD="{namd_exe}"',
            "",
            "# Check if NAMD is available",
            'if ! command -v $NAMD &> /dev/null; then',
            '    echo "Error: NAMD executable not found: $NAMD"',
            '    exit 1',
            'fi',
            "",
            'echo "Starting NAMD equilibration protocol..."',
            'echo "Each stage uses individual CPU/GPU settings"',
            'echo ""',
            ""
        ]
        
        # Add commands for each stage
        for i, (stage_key, stage_data) in enumerate(protocols.items()):
            stage_num = i + 1
            stage_name = stage_data.get('name', stage_key)
            steps = stage_data.get('steps', 'N/A')
            timestep = stage_data.get('timestep', 'N/A')
            use_gpu = stage_data.get('use_gpu', False)
            cpu_cores = stage_data.get('cpu_cores', 1)
            gpu_id = stage_data.get('gpu_id', 0)
            num_gpus = stage_data.get('num_gpus', 1)
            
            # Build NAMD command with appropriate flags
            namd_cmd = f'$NAMD'
            
            # Add processor specification
            namd_cmd += f' +p{cpu_cores}'
            
            # Add GPU flags if enabled
            if use_gpu:
                if num_gpus == 1:
                    namd_cmd += f' +devices {gpu_id}'  # Single GPU device
                else:
                    # Multiple GPUs: create device list starting from gpu_id
                    device_list = ','.join(str(gpu_id + i) for i in range(num_gpus))
                    namd_cmd += f' +devices {device_list}'
            
            # Complete command - use config-safe names for file names with new step naming
            config_name = self._get_config_name(stage_key, i)
            if config_name == "step7_production":
                namd_cmd += f' step7_production.conf > step7_production.log 2>&1'
            else:
                namd_cmd += f' {config_name}_equilibration.conf > {config_name}_equilibration.log 2>&1'
            
            # Create detailed resource information
            gpu_info = "No"
            if use_gpu:
                if num_gpus == 1:
                    gpu_info = f"Yes (GPU {gpu_id})"
                else:
                    gpu_list = ','.join(str(gpu_id + i) for i in range(num_gpus))
                    gpu_info = f"Yes ({num_gpus} GPUs: {gpu_list})"
            
            script_lines.extend([
                f'# Stage {stage_num}: {stage_name}',
                f'echo "Running Stage {stage_num}: {stage_name}"',
                f'echo "Steps: {steps}, Timestep: {timestep} ps"',
                f'echo "Resources: {cpu_cores} CPU cores, GPU: {gpu_info}"',
                namd_cmd,
                "",
                'if [ $? -ne 0 ]; then',
                f'    echo "Error in Stage {stage_num}: {stage_name}"',
                '    exit 1',
                'fi',
                f'echo "Stage {stage_num} completed successfully"',
                'echo ""',
                ""
            ])
        
        script_lines.extend([
            'echo "All equilibration stages completed successfully!"',
            'echo "Check the log files for detailed output"'
        ])
        
        return "\n".join(script_lines)
    
    def run_equilibration(self, config_files: List[Path], 
                         num_processors: int = 4) -> subprocess.Popen:
        """
        Run NAMD equilibration simulation.
        
        Args:
            config_files: List of NAMD configuration files to run
            num_processors: Number of processors to use
        
        Returns:
            Process object for the running simulation
        """
        
        # Create run script
        script_content = self._create_run_script(config_files, num_processors)
        
        # Write script to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as f:
            f.write(script_content)
            script_path = f.name
        
        # Make script executable
        os.chmod(script_path, 0o755)
        
        # Run script
        process = subprocess.Popen(
            ['bash', script_path],
            cwd=self.working_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        self.logger.info(f"Started NAMD equilibration with PID: {process.pid}")
        return process
    
    def _create_run_script(self, config_files: List[Path], 
                          num_processors: int) -> str:
        """Create a run script for the given configuration files."""
        
        script_lines = [
            "#!/bin/bash",
            f"# NAMD Equilibration Script",
            f"# Generated by Gatewizard",
            "",
            f'NAMD="{self.namd_executable}"',
            f'NPROCS={num_processors}',
            ""
        ]
        
        for i, config_file in enumerate(config_files):
            stage_num = i + 1
            script_lines.extend([
                f'echo "Running stage {stage_num}: {config_file.name}"',
                f'$NAMD +p$NPROCS {config_file.name}',
                f'if [ $? -ne 0 ]; then echo "Error in stage {stage_num}"; exit 1; fi',
                ""
            ])
        
        script_lines.append('echo "Equilibration completed successfully!"')
        
        return "\n".join(script_lines)
    
    def load_charmm_gui_template(self, scheme_type: str, stage_number: int, 
                                 system_files: Dict[str, str], 
                                 target_thickness: Optional[float] = None) -> str:
        """
        Load and customize CHARMM-GUI template file for a specific scheme and stage.
        
        Args:
            scheme_type: Type of scheme (NVT, NPT, NPAT, NPgT)
            stage_number: Stage number (1-12 for equilibration, 13 for production)
            system_files: Dictionary of system file paths
            target_thickness: Target bilayer thickness (deprecated, not used anymore)
            
        Returns:
            Customized NAMD configuration content
        """
        # Map scheme types to folder names
        scheme_folders = {
            "NVT": "01_NVT",
            "NPT": "02_NPT", 
            "NPAT": "03_NPAT",
            "NPgT": "04_NPgT"
        }
        
        if scheme_type not in scheme_folders:
            raise ValueError(f"Unknown scheme type: {scheme_type}")
        
        # Build template file path
        scheme_folder = scheme_folders[scheme_type]
        if stage_number <= 12:
            # Use stage 6 template for all equilibration stages 7-12
            if stage_number <= 6:
                template_file = f"step6.{stage_number}_equilibration.inp"
            else:
                template_file = "step6.6_equilibration.inp"  # Reuse last equilibration template
        elif stage_number == 13:
            template_file = "step7_production.inp"
        else:
            raise ValueError(f"Invalid stage number: {stage_number}")
            
        template_path = self.charmm_gui_templates_dir / scheme_folder / template_file
        
        if not template_path.exists():
            raise FileNotFoundError(f"CHARMM-GUI template not found: {template_path}")
        
        # Read template content
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # Customize template for Gatewizard
        customized_content = self._customize_charmm_gui_template_old(
            template_content, system_files, target_thickness
        )
        
        return customized_content
    
    def _customize_charmm_gui_template_old(self, template_content: str, 
                                      system_files: Dict[str, str],
                                      target_thickness: Optional[float] = None) -> str:
        """
        (DEPRECATED - kept for backward compatibility)
        Customize CHARMM-GUI template content for Gatewizard.
        
        Args:
            template_content: Original template content
            system_files: Dictionary of system file paths
            target_thickness: Target bilayer thickness for restraints
            
        Returns:
            Customized template content
        """
        lines = template_content.split('\n')
        customized_lines = []
        
        skip_restraint_section = False
        
        for line in lines:
            # Replace system file paths (use paths as-is from system_files, no hardcoded ../../)
            # Match both CHARMM-GUI format (step5_input.*) and Gatewizard format (system.*)
            if line.strip().startswith('parmfile') and ('step5_input.parm7' in line or 'system.prmtop' in line):
                customized_lines.append(f"parmfile                {system_files.get('prmtop', 'system.prmtop')}")
            elif line.strip().startswith('ambercoor') and ('step5_input.rst7' in line or 'system.inpcrd' in line):
                customized_lines.append(f"ambercoor               {system_files.get('inpcrd', 'system.inpcrd')}")
            
            # Skip planar and dihedral restraints as requested
            elif 'planar restraint' in line.lower() or 'dihedral restraint' in line.lower():
                skip_restraint_section = True
                customized_lines.append(f"# {line.strip()} - DISABLED BY GATEWIZARD")
                continue
            elif skip_restraint_section and (line.strip().startswith('#') or line.strip() == ''):
                customized_lines.append(line)
                continue  
            elif skip_restraint_section:
                if any(keyword in line.lower() for keyword in ['colvars', 'extrabonds', 'exec sed']):
                    customized_lines.append(f"# {line.strip()} - DISABLED BY GATEWIZARD")
                    continue
                else:
                    skip_restraint_section = False
                    
            # Bilayer thickness restraint is now handled in _generate_restraints_block
            else:
                customized_lines.append(line)
        
        return '\n'.join(customized_lines)
    
    def generate_bilayer_thickness_colvar(self, target_thickness: float, 
                                         output_path: Path, 
                                         pdb_file: Optional[Path] = None) -> None:
        """
        Generate collective variable file for bilayer thickness restraint.
        
        Args:
            target_thickness: Target bilayer thickness in Angstroms
            output_path: Path where to save the colvar file
            pdb_file: Optional PDB file for automatic phosphorus atom detection
        """
        # Try to use the new bilayer analyzer if PDB file is provided
        if pdb_file and pdb_file.exists():
            try:
                # Import the bilayer utilities
                import sys
                import os
                # Add the gatewizard utils path to sys.path
                gatewizard_root = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
                utils_path = os.path.join(gatewizard_root, 'utils')
                if utils_path not in sys.path:
                    sys.path.insert(0, utils_path)
                
                from bilayer_utils import generate_bilayer_thickness_colvar_from_pdb  # type: ignore
                
                # Generate colvar configuration using automatic analysis
                colvar_config = generate_bilayer_thickness_colvar_from_pdb(
                    str(pdb_file),
                    colvar_name="bilayer_thickness"
                )
                
                # Add header comment and restraint configuration
                colvar_content = f"""# Collective variable for bilayer thickness control
# Generated by Gatewizard using automatic bilayer analysis
# PDB file: {pdb_file}

{colvar_config}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant 10.0  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
                
                self.logger.info(f"Generated bilayer thickness colvar using automatic analysis from {pdb_file}")
                
            except (ImportError, Exception) as e:
                self.logger.warning(f"Could not use automatic bilayer analysis: {e}")
                self.logger.info("Falling back to generic colvar configuration")
                
                # Fall back to the old generic method
                colvar_content = f"""# Collective variable for bilayer thickness control
# Generated by Gatewizard (generic configuration)

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNameResidueRange P 1-999999
        }}
        group2 {{
            atomNameResidueRange P 1-999999  
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant 10.0  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
        else:
            # No PDB file provided or doesn't exist, use generic method
            self.logger.info("No PDB file provided for automatic analysis, using generic colvar configuration")
            colvar_content = f"""# Collective variable for bilayer thickness control
# Generated by Gatewizard (generic configuration)

colvar {{
    name bilayer_thickness
    
    # Distance between upper and lower leaflet phosphate groups
    distance {{
        group1 {{
            atomNameResidueRange P 1-999999
            atomNameResidueRange P 1-999999
        }}
        group2 {{
            atomNameResidueRange P 1-999999  
            atomNameResidueRange P 1-999999
        }}
    }}
}}

# Harmonic restraint to target thickness
harmonic {{
    colvars bilayer_thickness
    centers {target_thickness}
    forceConstant 10.0  # kcal/mol/Å²
}}

# Output settings
colvarsTrajFrequency 500
colvarsRestartFrequency 5000
"""
        
        # Ensure output directory exists
        from gatewizard.utils.helpers import create_directory_robust
        create_directory_robust(output_path.parent)
        
        # Write colvar file
        with open(output_path, 'w') as f:
            f.write(colvar_content)
        
        self.logger.info(f"Generated bilayer thickness colvar file: {output_path}")
    
    def generate_charmm_gui_config_file(self, stage_name: str, stage_params: Dict[str, Any],
                                      stage_index: int, system_files: Dict[str, str],
                                      scheme_type: str,
                                      previous_stage_name: Optional[str] = None,
                                      all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None,
                                      force_scheme_type: bool = False) -> str:
        """
        Generate NAMD configuration file using CHARMM-GUI templates with GateWizard customizations.
        
        Args:
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths
            scheme_type: CHARMM-GUI scheme type (NVT, NPT, NPAT, NPgT) - default ensemble for protocol
            previous_stage_name: Name of the previous stage for restart files
            force_scheme_type: If True, always use scheme_type for template selection, 
                             ignoring stage-specific ensemble values (GUI mode)
            
        Returns:
            Content of the NAMD configuration file
        """
        # Skip minimization stage - it's now incorporated into the first equilibration
        if stage_name == "minimization":
            self.logger.info("Skipping separate minimization stage - now included in eq1_equilibration")
            return ""
        
        # Map stage names to template files using config name mapping with stage_index
        config_name = self._get_config_name(stage_name, stage_index)
        
        # Determine which ensemble/scheme to use for template selection
        # Priority: 1) custom_template key, 2) stage's ensemble key (unless forced), 3) global scheme_type
        if force_scheme_type:
            # GUI mode: always use the selected scheme_type for all stages
            stage_ensemble = scheme_type
        else:
            # API mode: allow per-stage ensemble customization
            stage_ensemble = stage_params.get('ensemble', scheme_type)
        
        # Check if user explicitly specified a custom template
        custom_template = stage_params.get('custom_template', None)
        
        if custom_template:
            # User explicitly specified which template to use (e.g., 'step6.3_equilibration.inp')
            template_filename = custom_template
            template_scheme = stage_ensemble  # Use stage's ensemble for folder selection
            self.logger.info(f"Using custom template for {stage_name}: {custom_template} from {stage_ensemble} ensemble")
        else:
            # Auto-select template based on stage index and ensemble
            # Define template mapping based on config names (6 equilibration stages + production)
            template_mapping = {
                "step1": "step6.1_equilibration.inp",
                "step2": "step6.2_equilibration.inp", 
                "step3": "step6.3_equilibration.inp",
                "step4": "step6.4_equilibration.inp",
                "step5": "step6.5_equilibration.inp",
                "step6": "step6.6_equilibration.inp",
                "step7_production": "step7_production.inp"
            }
            
            # Get template filename based on stage position
            template_filename = template_mapping.get(config_name, "step6.1_equilibration.inp")
            
            # Determine which ensemble scheme to use for template folder selection
            if stage_ensemble != scheme_type:
                # Stage uses different ensemble than the protocol default
                self.logger.warning(
                    f"Stage {stage_index + 1} ({stage_name}) uses ensemble '{stage_ensemble}' "
                    f"but protocol default is '{scheme_type}'. Using '{stage_ensemble}' template."
                )
                template_scheme = stage_ensemble
            else:
                template_scheme = scheme_type
        
        # Load and customize template
        return self._load_and_customize_charmm_gui_template(
            template_scheme, template_filename, stage_name, stage_params, 
            stage_index, system_files, previous_stage_name, all_stage_settings
        )
    
    
    def _load_and_customize_charmm_gui_template(self, scheme_type: str, template_filename: str,
                                               stage_name: str, stage_params: Dict[str, Any],
                                               stage_index: int, system_files: Dict[str, str],
                                               previous_stage_name: Optional[str] = None,
                                               all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None) -> str:
        """
        Load CHARMM-GUI template and customize with GateWizard parameters.
        
        Args:
            scheme_type: CHARMM-GUI scheme type (NVT, NPT, NPAT, NPgT)
            template_filename: Template file to load
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths
            previous_stage_name: Name of the previous stage for restart files
            
        Returns:
            Customized NAMD configuration content
        """
        # Map scheme types to template directories
        scheme_mapping = {
            "NVT": "01_NVT",
            "NPT": "02_NPT", 
            "NPAT": "03_NPAT",
            "NPgT": "04_NPgT"
        }
        
        scheme_folder = scheme_mapping.get(scheme_type, "01_NVT")
        template_path = self.charmm_gui_templates_dir / scheme_folder / template_filename
        
        if not template_path.exists():
            raise FileNotFoundError(f"CHARMM-GUI template not found: {template_path}")
        
        # Read template content
        with open(template_path, 'r') as f:
            template_content = f.read()
        
        # Customize template with GateWizard parameters
        return self._customize_charmm_gui_template(
            template_content, stage_name, stage_params, stage_index, 
            system_files, previous_stage_name, all_stage_settings
        )
    
    def _customize_charmm_gui_template(self, template_content: str, stage_name: str,
                                      stage_params: Dict[str, Any], stage_index: int,
                                      system_files: Dict[str, str],
                                      previous_stage_name: Optional[str] = None,
                                      all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None) -> str:
        """
        Customize CHARMM-GUI template content with GateWizard parameters.
        
        Args:
            template_content: Raw template content
            stage_name: Name of the equilibration stage
            stage_params: Parameters for this stage
            stage_index: Index of the stage (0-based)
            system_files: Dictionary of system file paths
            previous_stage_name: Name of the previous stage for restart files
            
        Returns:
            Customized configuration content
        """
        # Temperature
        temperature = stage_params.get("temperature", 303.15)
        
        # DCD frequency for trajectory output
        dcd_freq = stage_params.get("dcd_freq", 5000)
        
        # Margin parameter for NPAT simulations
        margin = stage_params.get("margin", 5.0)
        
        # Time in nanoseconds and timestep
        time_ns = stage_params.get("time_ns", 0.125)  # Default 125 ps = 0.125 ns
        timestep = stage_params.get("timestep", 1.0)  # Default 1 fs (NAMD uses femtoseconds)
        
        # Calculate steps for display (but template uses the equation)
        calculated_steps = int(time_ns * 1e6 / timestep)
        
        # Minimization steps (only for first stage)
        minimize_steps = stage_params.get("minimize_steps", 10000) if stage_index == 0 else 0
        
        # Calculate firsttimestep based on previous stages
        first_timestep = self._calculate_first_timestep(stage_index, stage_params, all_stage_settings)
        
        # Handle cell basis vectors
        cell_basis_block = self._generate_cell_basis_block(stage_index)
        
        # Handle PME settings
        pme_block = self._generate_pme_block()
        
        # Handle restraints
        restraints_block = self._generate_restraints_block(stage_name, stage_params, stage_index)
        
        # Handle production steps for step13
        production_steps = stage_params.get("steps", 50000000)  # Default 50M steps
        
        # Generate input/output names for NAMD TCL variables
        output_name = self._generate_output_name(stage_name, stage_index)
        input_name = self._generate_input_name(stage_index, previous_stage_name)
        
        # Generate initial temperature directive (only for first stage)
        if stage_index == 0:  # First stage gets initial temperature assignment
            initial_temp_directive = f"temperature            $temp               # Initial temperature assignment for first stage"
        else:  # Subsequent stages don't set initial temperature (read from restart files)
            initial_temp_directive = "# No initial temperature assignment - reading velocities from restart file"
        
        # Perform replacements
        customized_content = template_content.replace("{TEMPERATURE}", str(temperature))
        customized_content = customized_content.replace("{DCD_FREQ}", str(dcd_freq))
        customized_content = customized_content.replace("{OUTPUT_ENERGIES}", str(dcd_freq))
        customized_content = customized_content.replace("{XST_FREQ}", str(dcd_freq))
        customized_content = customized_content.replace("{OUTPUT_TIMING}", str(dcd_freq))
        customized_content = customized_content.replace("{MARGIN}", str(margin))
        customized_content = customized_content.replace("{TIME_NS}", str(time_ns))
        customized_content = customized_content.replace("{TIMESTEP}", str(timestep))
        customized_content = customized_content.replace("{CELL_BASIS_VECTORS}", cell_basis_block)
        customized_content = customized_content.replace("{PME_SETTINGS}", pme_block)
        customized_content = customized_content.replace("{RESTRAINTS_BLOCK}", restraints_block)
        customized_content = customized_content.replace("{PRODUCTION_STEPS}", str(production_steps))
        customized_content = customized_content.replace("{RUN_STEPS}", str(calculated_steps))
        customized_content = customized_content.replace("{MINIMIZE_STEPS}", str(minimize_steps))
        customized_content = customized_content.replace("{FIRST_TIMESTEP}", str(first_timestep))
        customized_content = customized_content.replace("{INITIAL_TEMPERATURE_DIRECTIVE}", initial_temp_directive)
        
        # Replace system file paths (parmfile and ambercoor)
        # This allows using either relative paths (when files are copied) or absolute paths
        import re
        if 'prmtop' in system_files:
            prmtop_path = system_files['prmtop']
            customized_content = re.sub(
                r'parmfile\s+[\w/.]+\.prmtop',
                f'parmfile                {prmtop_path}',
                customized_content
            )
            customized_content = re.sub(
                r'parmfile\s+[\w/.]+\.parm7',
                f'parmfile                {prmtop_path}',
                customized_content
            )
            customized_content = re.sub(
                r'parmfile\s+[\w/.]+\.top',
                f'parmfile                {prmtop_path}',
                customized_content
            )
        
        if 'inpcrd' in system_files:
            inpcrd_path = system_files['inpcrd']
            customized_content = re.sub(
                r'ambercoor\s+[\w/.]+\.inpcrd',
                f'ambercoor               {inpcrd_path}',
                customized_content
            )
            customized_content = re.sub(
                r'ambercoor\s+[\w/.]+\.rst7',
                f'ambercoor               {inpcrd_path}',
                customized_content
            )
            customized_content = re.sub(
                r'ambercoor\s+[\w/.]+\.rst',
                f'ambercoor               {inpcrd_path}',
                customized_content
            )
        
        # Replace NAMD TCL variable names if they exist as placeholders
        customized_content = customized_content.replace("{OUTPUT_NAME}", output_name)
        customized_content = customized_content.replace("{INPUT_NAME}", input_name)
        
        # Replace hardcoded outputname and inputname in templates using regex
        
        # Replace set outputname lines (e.g., "set outputname eq6_equilibration;" -> "set outputname eq7_equilibration;")
        customized_content = re.sub(
            r'set\s+outputname\s+\w+;',
            f'set outputname          {output_name};',
            customized_content
        )
        
        # Replace set inputname lines (e.g., "set inputname eq5_equilibration;" -> "set inputname eq6_equilibration;")
        # Only replace if we have an input name (not first stage)
        if input_name:
            customized_content = re.sub(
                r'set\s+inputname\s+\w+;',
                f'set inputname           {input_name};',
                customized_content
            )
        else:
            # For first stage, remove or comment out the inputname line AND the restart file directives
            customized_content = re.sub(
                r'set\s+inputname\s+\w+;',
                '# set inputname           (not needed for first stage);',
                customized_content
            )
            # Also comment out the restart file directives that depend on inputname
            customized_content = re.sub(
                r'binCoordinates\s+\$inputname\.coor;',
                '# binCoordinates          $inputname.coor;    # (not needed for first stage)',
                customized_content
            )
            customized_content = re.sub(
                r'binVelocities\s+\$inputname\.vel;',
                '# binVelocities           $inputname.vel;     # (not needed for first stage)',
                customized_content
            )
            customized_content = re.sub(
                r'extendedSystem\s+\$inputname\.xsc;',
                '# extendedSystem          $inputname.xsc;     # (not needed for first stage)',
                customized_content
            )
        
        return customized_content
    
    def _calculate_first_timestep(self, stage_index: int, stage_params: Dict[str, Any], 
                                 all_stage_settings: Optional[Dict[str, Dict[str, Any]]] = None) -> int:
        """Calculate the first timestep for a stage based on previous stages."""
        if stage_index == 0:
            return 0  # First stage always starts from 0
        
        if all_stage_settings is None:
            # Fallback to old behavior if all_stage_settings not provided
            current_stage_steps = stage_params.get("steps", 125000)
            cumulative_steps = stage_index * current_stage_steps
            return cumulative_steps
        
        # Build list of actual stage keys from all_stage_settings in order
        # The keys might be in format "Equilibration 1", "Equilibration 2", ..., "Production"
        # or "equilibration_1", "equilibration_2", ..., "production"
        stage_keys = list(all_stage_settings.keys())
        
        cumulative_steps = 0
        for i in range(stage_index):
            if i < len(stage_keys):
                stage_key = stage_keys[i]
                stage_config = all_stage_settings[stage_key]
                
                # Get run steps based on time_ns and timestep
                time_ns = stage_config.get('time_ns', 0.125)
                timestep = stage_config.get('timestep', 1.0)
                run_steps = int(time_ns * 1e6 / timestep)
                
                # Add minimize steps for first stage only
                if i == 0:
                    minimize_steps = stage_config.get('minimize_steps', 10000)
                    cumulative_steps += minimize_steps
                
                cumulative_steps += run_steps
            else:
                # Default steps if stage not found
                cumulative_steps += 125000
        
        return cumulative_steps

    def _generate_output_name(self, stage_name: str, stage_index: int) -> str:
        """Generate output name for NAMD configuration."""
        # Use the new config name mapping function with stage_index
        config_name = self._get_config_name(stage_name, stage_index)
        # Production uses step7_production, equilibration stages use step{N}_equilibration
        if config_name == "step7_production":
            return config_name
        else:
            return f"{config_name}_equilibration"
    
    def _generate_input_name(self, stage_index: int, previous_stage_name: Optional[str] = None) -> str:
        """Generate input name for NAMD restart files."""
        if stage_index == 0:
            # First stage doesn't need input name
            return ""
        
        if previous_stage_name:
            # Use previous stage output name - convert previous stage name to config name
            # Use stage_index - 1 for the previous stage
            prev_config_name = self._get_config_name(previous_stage_name, stage_index - 1)
            # Production uses step7_production, equilibration stages use step{N}_equilibration
            if prev_config_name == "step7_production":
                return prev_config_name
            else:
                return f"{prev_config_name}_equilibration"
        
        # Fallback: generate based on stage index using direct step naming
        # For stage_index 1 (2nd stage) -> need step1_equilibration as input
        # For stage_index 2 (3rd stage) -> need step2_equilibration as input, etc.
        if stage_index > 0:
            prev_step_num = stage_index  # This gives us the previous step number
            return f"step{prev_step_num}_equilibration"
        
        # Final fallback (should never reach here)
        return ""
    
    def _generate_cell_basis_block(self, stage_index: int) -> str:
        """Generate cell basis vectors block - only for first stage, others use .xsc files."""
        if stage_index > 0:
            # After first stage, box dimensions come from .xsc files
            return "# Cell dimensions read from .xsc file"
        
        # For first stage, get box dimensions with proper priority order
        try:
            # PRIORITY 1: Read from bilayer*_lipid.pdb files with CRYST1 records (most accurate)
            bilayer_pdb = self._find_bilayer_pdb_with_cryst1()
            if bilayer_pdb and bilayer_pdb.exists():
                a, b, c = self._read_box_dimensions(bilayer_pdb)
                self.logger.info(f"Using box dimensions from bilayer*_lipid.pdb for cell basis: {bilayer_pdb.name}")
            else:
                # PRIORITY 2: Try system.pdb as fallback
                system_pdb = self.working_dir / "system.pdb"
                if system_pdb.exists():
                    a, b, c = self._read_box_dimensions(system_pdb)
                    self.logger.info(f"Using box dimensions from system.pdb for cell basis: {system_pdb.name}")
                else:
                    # PRIORITY 3: Default dimensions
                    a, b, c = 100.0, 100.0, 100.0
                    self.logger.warning("No PDB files found, using default box dimensions for cell basis")
            
            zcen = 0.0
        except Exception as e:
            self.logger.warning(f"Error reading box dimensions for cell basis: {e}, using defaults")
            a, b, c = 100.0, 100.0, 100.0
            zcen = 0.0
        
        return f"""cellBasisVector1     {a:.3f}   0.0   0.0
cellBasisVector2     0.0   {b:.3f}   0.0
cellBasisVector3     0.0   0.0   {c:.3f}
cellOrigin           0.0   0.0   {zcen:.3f}"""
    
    def _generate_pme_block(self) -> str:
        """Generate PME settings block without hardcoded grid sizes."""
        # PME grid sizes should be calculated automatically by NAMD
        # based on the system size and PMEGridSpacing
        return f"""PME                     yes
PMEInterpOrder          6
PMEGridSpacing          1.0"""
    
    def _generate_restraints_block(self, stage_name: str, stage_params: Dict[str, Any], stage_index: int) -> str:
        """Generate restraints block using GateWizard's restraint system."""
        constraints = stage_params.get("constraints", {})
        
        # Check if any constraints are defined
        has_restraints = any(float(v) > 0 for v in constraints.values())
        
        restraints_lines = []
        
        # Add position restraints if defined
        if has_restraints:
            # Map stage names to restraint file names using the correct naming scheme
            # Based on the display name to config name conversion
            config_name = self._get_config_name(stage_name, stage_index)
            
            if config_name == "step7_production":
                restraint_file = f"{config_name}_restraints.pdb"  # step7_production_restraints.pdb
            else:
                # For equilibration stages: step1_equilibration_restraints.pdb, step2_equilibration_restraints.pdb, etc.
                restraint_file = f"{config_name}_equilibration_restraints.pdb"
            
            # Use constraintScaling = 1.0 and keep GUI force values
            restraints_lines.extend([
                "# Position restraints",
                "constraints             on",
                "consexp                 2",
                f"consref                 restraints/{restraint_file}",
                f"conskfile               restraints/{restraint_file}",
                "conskcol                B",
                "constraintScaling       1.0"
            ])
        
        # If no restraints at all, add a comment
        if not restraints_lines:
            return "# No restraints defined for this stage"
        
        return "\n".join(restraints_lines)

class GromacsEquilibrationManager:
    """Manager for GROMACS equilibration simulations (placeholder)."""
    
    def __init__(self, working_dir: Path):
        self.working_dir = Path(working_dir)
        self.logger = get_logger(self.__class__.__name__)
    
    def generate_mdp_file(self, stage_name: str, stage_params: Dict[str, Any]) -> str:
        """Generate GROMACS MDP file (placeholder)."""
        self.logger.info("GROMACS equilibration not yet implemented")
        return ""

class AmberEquilibrationManager:
    """Manager for AMBER equilibration simulations (placeholder)."""
    
    def __init__(self, working_dir: Path):
        self.working_dir = Path(working_dir)
        self.logger = get_logger(self.__class__.__name__)
    
    def generate_input_file(self, stage_name: str, stage_params: Dict[str, Any]) -> str:
        """Generate AMBER input file (placeholder)."""
        self.logger.info("AMBER equilibration not yet implemented")
        return ""

class EquilibrationAnalyzer:
    """Analyzer for equilibration simulation results."""
    
    def __init__(self, working_dir: Path):
        self.working_dir = Path(working_dir)
        self.logger = get_logger(self.__class__.__name__)
    
    def analyze_energy_convergence(self, log_files: List[Path]) -> Dict[str, Any]:
        """
        Analyze energy convergence from simulation log files.
        
        Args:
            log_files: List of simulation log files
        
        Returns:
            Dictionary with convergence analysis results
        """
        
        results = {
            "converged": False,
            "total_energy": [],
            "temperature": [],
            "pressure": [],
            "volume": []
        }
        
        # Placeholder implementation
        self.logger.info("Energy convergence analysis not yet implemented")
        
        return results
    
    def generate_plots(self, analysis_results: Dict[str, Any], 
                      output_dir: Path) -> List[Path]:
        """
        Generate plots for equilibration analysis.
        
        Args:
            analysis_results: Results from analysis
            output_dir: Directory to save plots
        
        Returns:
            List of generated plot files
        """
        
        plots = []
        
        # Placeholder for plot generation
        self.logger.info("Plot generation not yet implemented")
        
        return plots
