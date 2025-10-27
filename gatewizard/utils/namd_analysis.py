"""
NAMD log analysis utilities for extracting timing and performance information.

This module provides functions to parse NAMD log files and extract useful
information like simulation progress, performance metrics, and timing data.
"""

import re
from pathlib import Path
from typing import Dict, List, Optional, NamedTuple, TYPE_CHECKING, Union
from dataclasses import dataclass
from .logger import get_logger

if TYPE_CHECKING:
    import numpy as np

logger = get_logger(__name__)

@dataclass
class NAMDTiming:
    """Container for NAMD timing information."""
    steps_completed: int = 0
    total_steps: int = 0
    simulated_time_ns: float = 0.0
    real_time_hours: float = 0.0
    ns_per_day: float = 0.0
    sec_per_step: float = 0.0
    processors: int = 0
    gpus: int = 0
    atoms: int = 0
    timestep_fs: float = 0.0
    first_timestep: int = 0
    hostname: str = ""

@dataclass
class NAMDProgress:
    """Container for NAMD progress information."""
    stage_name: str = ""
    status: str = "not_started"  # not_started, running, completed, error
    progress_percent: float = 0.0
    timing: Optional[NAMDTiming] = None
    log_file: Optional[Path] = None
    last_updated: float = 0.0

def parse_namd_log(log_file_path: Path) -> NAMDTiming:
    """
    Parse a NAMD log file to extract timing and performance information.
    
    Based on the namd_timing script functionality.
    """
    timing = NAMDTiming()
    
    if not log_file_path.exists():
        logger.debug(f"Log file does not exist: {log_file_path}")
        return timing
    
    try:
        with open(log_file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        logger.debug(f"Parsing log file {log_file_path.name}, size: {len(content)} chars")
        
        # Extract basic system information
        proc_match = re.search(r'Running on (\d+) processors', content)
        if proc_match:
            timing.processors = int(proc_match.group(1))
        
        atoms_match = re.search(r'Info: (\d+) ATOMS', content)
        if atoms_match:
            timing.atoms = int(atoms_match.group(1))
        
        timestep_match = re.search(r'Info: TIMESTEP\s+(\d+(?:\.\d+)?)', content)
        if timestep_match:
            timing.timestep_fs = float(timestep_match.group(1))
        
        first_ts_match = re.search(r'FIRST TIMESTEP\s+(\d+)', content)
        if first_ts_match:
            timing.first_timestep = int(first_ts_match.group(1))
        
        # Extract hostname
        host_match = re.search(r'Info: \d+ NAMD.*?(\S+)', content)
        if host_match:
            timing.hostname = host_match.group(1)
        
        # Count CUDA devices
        cuda_matches = re.findall(r'CUDA device \d+', content)
        timing.gpus = len(cuda_matches)
        
        # Extract TIMING lines for performance analysis
        timing_lines = re.findall(r'^TIMING:\s+(\d+)\s+[\d\.\-\+e]+\s+[\d\.\-\+e]+\s+[\d\.\-\+e]+\s+[\d\.\-\+e]+\s+[\d\.\-\+e]+\s+[\d\.\-\+e]+\s+([\d\.\-\+e]+)', content, re.MULTILINE)
        
        # Try alternative TIMING patterns if the first doesn't work
        if not timing_lines:
            # Try simpler TIMING pattern
            timing_lines = re.findall(r'^TIMING:\s+(\d+)\s+.*?(\d+\.\d+)', content, re.MULTILINE)
            logger.debug(f" Using alternative TIMING pattern")
        
        # Also try ENERGY lines as an alternative
        if not timing_lines:
            energy_lines = re.findall(r'^ENERGY:\s+(\d+)', content, re.MULTILINE)
            if energy_lines:
                logger.debug(f" Found {len(energy_lines)} ENERGY lines, using last one")
                last_step = energy_lines[-1]
                timing_lines = [(last_step, "0.0")]  # Fake time since we don't have it from ENERGY
        
        # Check for final completion step patterns when simulation finishes
        # Look for patterns like "WRITING ... TO OUTPUT FILE AT STEP <number>"
        final_step_matches = re.findall(r'WRITING.*?TO OUTPUT FILE AT STEP (\d+)', content)
        if final_step_matches:
            final_step = int(final_step_matches[-1])
            logger.debug(f" Found final output step: {final_step}")
            # If this final step is higher than our last timing step, use it
            if timing_lines:
                last_timing_step = int(timing_lines[-1][0])
                if final_step > last_timing_step:
                    logger.debug(f" Using final step {final_step} instead of last timing step {last_timing_step}")
                    timing_lines.append((str(final_step), "0.0"))
            else:
                timing_lines = [(str(final_step), "0.0")]
        
        logger.debug(f" Found {len(timing_lines)} TIMING/ENERGY lines")
        
        if timing_lines:
            # Get the last timing line for current progress
            last_step, last_time = timing_lines[-1]
            timing.steps_completed = int(last_step) - timing.first_timestep
            logger.debug(f" Last step: {last_step}, first: {timing.first_timestep}, completed: {timing.steps_completed}")
            
            # Calculate average time per step
            total_time = 0.0
            step_count = 0
            for step_str, time_str in timing_lines:
                total_time += float(time_str)
                step_count += 1
            
            if step_count > 0:
                timing.sec_per_step = total_time / step_count
                timing.real_time_hours = float(last_time) / 3600.0
                
                # Calculate simulated time in nanoseconds
                if timing.timestep_fs > 0:
                    timing.simulated_time_ns = (timing.steps_completed * timing.timestep_fs) / 1000000.0
                    
                    # Calculate ns/day performance
                    if timing.real_time_hours > 0:
                        timing.ns_per_day = timing.simulated_time_ns / (timing.real_time_hours / 24.0)
        
        # Try to get total steps from the configuration
        # Primary pattern: TCL commands (NAMD 3.0 format)
        tcl_run_match = re.search(r'TCL: Running for (\d+) steps', content, re.IGNORECASE)
        tcl_minimize_match = re.search(r'TCL: Minimizing for (\d+) steps', content, re.IGNORECASE)
        
        # Fallback patterns for older NAMD versions
        run_match = re.search(r'run\s+(\d+)', content)
        minimize_match = re.search(r'minimize\s+(\d+)', content)
        numsteps_match = re.search(r'numsteps\s+(\d+)', content)
        
        # Initialize with 0
        run_steps = 0
        minimize_steps = 0
        
        # Get run steps
        if tcl_run_match:
            run_steps = int(tcl_run_match.group(1))
            logger.debug(f" Found 'TCL: Running' pattern: {run_steps} steps")
        elif run_match:
            run_steps = int(run_match.group(1))
            logger.debug(f" Found 'run' pattern: {run_steps} steps")
        elif numsteps_match:
            run_steps = int(numsteps_match.group(1))
            logger.debug(f" Found 'numsteps' pattern: {run_steps} steps")
        
        # Get minimize steps
        if tcl_minimize_match:
            minimize_steps = int(tcl_minimize_match.group(1))
            logger.debug(f" Found 'TCL: Minimizing' pattern: {minimize_steps} steps")
        elif minimize_match:
            minimize_steps = int(minimize_match.group(1))
            logger.debug(f" Found 'minimize' pattern: {minimize_steps} steps")
        
        # Total steps is the sum of minimization and run steps
        if run_steps > 0 or minimize_steps > 0:
            timing.total_steps = run_steps + minimize_steps
            logger.debug(f" Total steps calculation: {minimize_steps} (minimize) + {run_steps} (run) = {timing.total_steps}")
        
        # If we only found minimize steps but no run steps yet (simulation in progress),
        # try to read the input file to get the expected total steps
        if minimize_steps > 0 and run_steps == 0:
            inp_file_steps = _get_expected_steps_from_inp_file(log_file_path)
            if inp_file_steps > minimize_steps:
                timing.total_steps = inp_file_steps
                logger.debug(f" Using expected steps from input file: {inp_file_steps}")
        
        # Also try to find steps from NAMD output messages if we couldn't find run/minimize commands
        if timing.total_steps == 0:
            # Look for patterns like "Info: STEPS <number>"
            steps_info_match = re.search(r'Info:.*?STEPS?\s+(\d+)', content, re.IGNORECASE)
            if steps_info_match:
                timing.total_steps = int(steps_info_match.group(1))
            
            # Look for patterns in lines like "ETITLE:" or similar
            etitle_match = re.search(r'ETITLE:.*?(\d+)', content)
            if etitle_match and not timing.total_steps:
                # This might be less reliable, use as last resort
                pass
    
    except Exception as e:
        # Log error but return partial results
        pass
    
    return timing

def _get_expected_steps_from_inp_file(log_file_path: Path) -> int:
    """
    Try to read the corresponding input file to get expected total steps.
    
    Args:
        log_file_path: Path to the log file (e.g., eq1_equilibration.log)
        
    Returns:
        Expected total steps from input file, or 0 if not found
    """
    try:
        # Try to find corresponding .inp file
        inp_file = log_file_path.with_suffix('.inp')
        if not inp_file.exists():
            # Try common naming patterns
            base_name = log_file_path.stem.replace('_equilibration', '').replace('.log', '')
            possible_inp_files = [
                log_file_path.parent / f"{base_name}_equilibration.inp",
                log_file_path.parent / f"{base_name}.inp",
                log_file_path.parent / f"step6.{base_name.replace('eq', '')}_equilibration.inp"
            ]
            
            for possible_inp in possible_inp_files:
                if possible_inp.exists():
                    inp_file = possible_inp
                    break
            else:
                return 0
        
        with open(inp_file, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Look for minimize and run commands
        minimize_match = re.search(r'minimize\s+(\d+)', content)
        run_match = re.search(r'run\s+(\d+)', content)
        
        minimize_steps = int(minimize_match.group(1)) if minimize_match else 0
        run_steps = int(run_match.group(1)) if run_match else 0
        
        total_expected = minimize_steps + run_steps
        logger.debug(f" Expected steps from {inp_file.name}: {minimize_steps} (minimize) + {run_steps} (run) = {total_expected}")
        
        return total_expected
        
    except Exception as e:
        logger.debug(f" Error reading input file for {log_file_path}: {e}")
        return 0

def get_equilibration_progress(equilibration_dir: Path) -> Dict[str, NAMDProgress]:
    """
    Get progress information for all equilibration stages.
    
    Looks for log files in the equilibration/namd directory structure.
    """
    progress = {}
    
    # Standard stage names and their expected log file patterns
    # Updated to match current equilibration structure (6 equilibration stages + production)
    # step1* is for equilibration_1, step2* for equilibration_2, etc.
    stage_patterns = {
        'equilibration_1': ['step1_equilibration*.log', 'step1_*.log'],
        'equilibration_2': ['step2_equilibration*.log', 'step2_*.log'],
        'equilibration_3': ['step3_equilibration*.log', 'step3_*.log'],
        'equilibration_4': ['step4_equilibration*.log', 'step4_*.log'],
        'equilibration_5': ['step5_equilibration*.log', 'step5_*.log'],
        'equilibration_6': ['step6_equilibration*.log', 'step6_*.log'],
        'production': ['step7_production*.log', 'production*.log', 'prod*.log']
    }
    
    # Look for equilibration directory
    eq_namd_dir = equilibration_dir / 'namd'
    logger.debug(f"Checking primary path: {eq_namd_dir}")
    logger.debug(f"Primary path exists: {eq_namd_dir.exists()}")
    
    if not eq_namd_dir.exists():
        # Try alternative locations
        eq_namd_dir = equilibration_dir / 'equilibration' / 'namd'
        logger.debug(f"Trying fallback path: {eq_namd_dir}")
        logger.debug(f"Fallback path exists: {eq_namd_dir.exists()}")
        
        if not eq_namd_dir.exists():
            eq_namd_dir = equilibration_dir
            logger.debug(f"Using base directory: {eq_namd_dir}")
    
    logger.debug(f"Final search directory: {eq_namd_dir}")
    
    if eq_namd_dir.exists():
        all_files = list(eq_namd_dir.glob("*"))
        logger.debug(f" All files in directory: {[f.name for f in all_files]}")
        log_files = list(eq_namd_dir.glob("*.log"))
    
    for stage_name, patterns in stage_patterns.items():
        stage_progress = NAMDProgress(stage_name=stage_name)
        
        # Look for log files matching the patterns
        log_file = None
        for pattern in patterns:
            matches = list(eq_namd_dir.glob(pattern))
            if matches:
                # Use the most recent log file
                log_file = max(matches, key=lambda f: f.stat().st_mtime)
                break
        
        if log_file and log_file.exists():
            stage_progress.log_file = log_file
            stage_progress.last_updated = log_file.stat().st_mtime
            
            # Parse the log file
            timing = parse_namd_log(log_file)
            stage_progress.timing = timing
            
            # Determine status and progress
            if timing.steps_completed > 0:
                stage_progress.status = "running"
                if timing.total_steps > 0:
                    stage_progress.progress_percent = (timing.steps_completed / timing.total_steps) * 100.0
                    if timing.steps_completed >= timing.total_steps:
                        stage_progress.status = "completed"
                else:
                    # If we can't determine total steps, consider it running
                    stage_progress.progress_percent = 50.0  # Unknown progress
            else:
                stage_progress.status = "not_started"
        else:
            pass  # No log file found for stage
        
        progress[stage_name] = stage_progress
    
    return progress

def format_timing_info(timing: NAMDTiming) -> str:
    """Format timing information for display."""
    if not timing or timing.steps_completed == 0:
        return "No timing data available"
    
    lines = []
    
    # Basic info
    if timing.hostname:
        lines.append(f"Host: {timing.hostname}")
    
    if timing.processors > 0:
        gpu_info = f", {timing.gpus} GPUs" if timing.gpus > 0 else ""
        lines.append(f"Resources: {timing.processors} processors{gpu_info}")
    
    if timing.atoms > 0:
        lines.append(f"System: {timing.atoms:,} atoms")
    
    # Progress info
    if timing.total_steps > 0:
        lines.append(f"Steps: {timing.steps_completed:,} / {timing.total_steps:,}")
    else:
        lines.append(f"Steps completed: {timing.steps_completed:,}")
    
    # Performance info
    if timing.simulated_time_ns > 0:
        lines.append(f"Simulated time: {timing.simulated_time_ns:.3f} ns")
    
    if timing.real_time_hours > 0:
        lines.append(f"Real time: {timing.real_time_hours:.2f} hours")
    
    if timing.ns_per_day > 0:
        lines.append(f"Performance: {timing.ns_per_day:.4f} ns/day")
    
    if timing.sec_per_step > 0:
        lines.append(f"Time per step: {timing.sec_per_step:.3f} sec")
    
    return "\n".join(lines)

def format_progress_summary(progress_dict: Dict[str, NAMDProgress]) -> str:
    """Format a summary of all stage progress."""
    lines = ["Equilibration Progress Summary:"]
    
    for stage_name, progress in progress_dict.items():
        status_icon = {
            "not_started": "⏸️",
            "running": "🏃",
            "completed": "✅",
            "error": "❌"
        }.get(progress.status, "❓")
        
        stage_display = stage_name.replace('_', ' ').title()
        line = f"  {status_icon} {stage_display}: {progress.status}"
        
        if progress.progress_percent > 0:
            line += f" ({progress.progress_percent:.1f}%)"
        
        if progress.timing and progress.timing.ns_per_day > 0:
            line += f" - {progress.timing.ns_per_day:.4f} ns/day"
        
        lines.append(line)
    
    return "\n".join(lines)

# ============================================================================
# High-level Wrapper Classes for Easy Analysis
# ============================================================================

class EnergyAnalyzer:
    """
    Easy-to-use wrapper for NAMD energy analysis with built-in plotting and full GUI capabilities.
    
    Supports:
    - Single or multiple log files with custom time scaling
    - Selective energy property plotting
    - Multiple plots (same figure or separate)
    - Full customization (colors, grid, units, target values, etc.)
    
    Example:
        >>> # Single file - plot
        >>> analyzer = EnergyAnalyzer("step1_equilibration.log")
        >>> analyzer.plot_energy(save="energy.png")
        
        >>> # With custom target temperature and pressure
        >>> analyzer.plot_energy(
        ...     target_temperature=310.0,  # 310 K
        ...     target_pressure=1.01325,   # 1.01325 atm (1 bar)
        ...     save="energy.png"
        ... )
        
        >>> # Multiple files with custom time
        >>> analyzer = EnergyAnalyzer(
        ...     ["step1.log", "step2.log", "step3.log"],
        ...     file_times={"step1.log": 0.05, "step2.log": 0.05, "step3.log": 0.05}
        ... )
        
        >>> # Plot specific properties
        >>> analyzer.plot_properties(
        ...     properties=["Temperature", "Total Energy", "Pressure"],
        ...     save="specific_energies.png"
        ... )
        
        >>> # Plot each property separately
        >>> analyzer.plot_properties(
        ...     properties=["Temperature", "Pressure"],
        ...     separate_plots=True,
        ...     save_prefix="plot_"
        ... )
    """
    
    def __init__(self, 
                 log_file: Union[Path, str, List[Union[Path, str]]],
                 file_times: Optional[Dict[str, float]] = None):
        """
        Initialize energy analyzer with NAMD log file(s).
        
        Args:
            log_file: Path to NAMD log file, or list of paths for multi-file analysis
            file_times: Dict mapping filename (just the name, not full path) to duration in ns
                       Example: {"step1.log": 0.05, "step2.log": 0.05}
                       Time is the DURATION of each file, not cumulative
        """
        # Handle single file or list
        if isinstance(log_file, (str, Path)):
            self.log_files = [Path(log_file)]
        else:
            self.log_files = [Path(f) for f in log_file]
        
        # Store file times (duration of each file in ns)
        self.file_times = file_times or {}
        
        # Track file ranges for time calculation (must be before parsing!)
        self._file_ranges = {}  # {filepath: (start_idx, end_idx, min_ts, max_ts)}
        
        # Parse all log files
        self.data = self._parse_energy_data()
        self.timing = parse_namd_log(self.log_files[0]) if self.log_files else None
        
    def _parse_energy_data(self) -> Dict[str, List[float]]:
        """Parse ENERGY lines from NAMD log file(s)."""
        data = {
            'timestep': [],
            'bond': [],
            'angle': [],
            'dihedral': [],
            'improper': [],
            'elect': [],
            'vdw': [],
            'boundary': [],
            'misc': [],
            'kinetic': [],
            'total': [],
            'temp': [],
            'potential': [],
            'total3': [],
            'tempavg': [],
            'pressure': [],
            'gpressure': [],
            'volume': [],
            'pressavg': [],
            'gpressavg': []
        }
        
        for log_file in self.log_files:
            if not log_file.exists():
                logger.warning(f"Log file not found: {log_file}")
                continue
            
            start_idx = len(data['timestep'])
            
            try:
                with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                    content = f.read()
                
                # Find ETITLE line to get column labels
                etitle_match = re.search(r'ETITLE:\s+(.+)', content)
                if not etitle_match:
                    logger.warning(f"No ETITLE line found in {log_file}")
                    continue
                
                # Parse ENERGY lines
                energy_lines = re.findall(r'^ENERGY:\s+(.+)$', content, re.MULTILINE)
                
                min_ts, max_ts = None, None
                
                for line in energy_lines:
                    values = line.split()
                    if len(values) >= 14:  # Minimum expected columns
                        try:
                            ts = int(values[0])
                            if min_ts is None:
                                min_ts = ts
                            max_ts = ts
                            
                            data['timestep'].append(ts)
                            data['bond'].append(float(values[1]))
                            data['angle'].append(float(values[2]))
                            data['dihedral'].append(float(values[3]))
                            data['improper'].append(float(values[4]))
                            data['elect'].append(float(values[5]))
                            data['vdw'].append(float(values[6]))
                            data['boundary'].append(float(values[7]))
                            data['misc'].append(float(values[8]))
                            data['kinetic'].append(float(values[9]))
                            data['total'].append(float(values[10]))
                            data['temp'].append(float(values[11]))
                            data['potential'].append(float(values[12]))
                            data['total3'].append(float(values[13]))
                            
                            if len(values) >= 15:
                                data['tempavg'].append(float(values[14]))
                            if len(values) >= 16:
                                data['pressure'].append(float(values[15]))
                            if len(values) >= 17:
                                data['gpressure'].append(float(values[16]))
                            if len(values) >= 18:
                                data['volume'].append(float(values[17]))
                            if len(values) >= 19:
                                data['pressavg'].append(float(values[18]))
                            if len(values) >= 20:
                                data['gpressavg'].append(float(values[19]))
                        except (ValueError, IndexError) as e:
                            continue
                
                # Store file range information
                end_idx = len(data['timestep'])
                self._file_ranges[str(log_file)] = (start_idx, end_idx, min_ts or 0, max_ts or 0)
                        
            except Exception as e:
                logger.error(f"Error parsing {log_file}: {e}")
                
        return data
    
    def get_statistics(self) -> Dict[str, Dict[str, float]]:
        """
        Get statistical summary of energy data.
        
        Returns:
            Dictionary with statistics for each energy component
        """
        import numpy as np
        
        stats = {}
        for key, values in self.data.items():
            if values and key != 'timestep':
                arr = np.array(values)
                stats[key] = {
                    'mean': float(np.mean(arr)),
                    'std': float(np.std(arr)),
                    'min': float(np.min(arr)),
                    'max': float(np.max(arr)),
                    'initial': float(arr[0]) if len(arr) > 0 else 0.0,
                    'final': float(arr[-1]) if len(arr) > 0 else 0.0
                }
        
        return stats
    
    def _calculate_time_array(self) -> 'np.ndarray':
        """
        Calculate time array for all data points with proper file time scaling.
        
        Returns:
            Time array in nanoseconds
        """
        import numpy as np
        
        time_array_ns = []
        
        # Check if we have custom time assignments
        has_custom_times = bool(self.file_times) and any(t > 0 for t in self.file_times.values())
        
        if has_custom_times and self._file_ranges:
            # Use custom time assignments
            cumulative_time_ns = 0.0
            
            for log_file in self.log_files:
                filepath_str = str(log_file)
                if filepath_str not in self._file_ranges:
                    continue
                
                start_idx, end_idx, min_ts, max_ts = self._file_ranges[filepath_str]
                num_points = end_idx - start_idx
                
                if num_points == 0:
                    continue
                
                # Get assigned time for this file (duration in ns)
                filename = log_file.name
                assigned_time_ns = self.file_times.get(filename, 0.0)
                
                if assigned_time_ns <= 0:
                    # No time assigned, use timestep-based calculation
                    timestep_fs = self.timing.timestep_fs if self.timing and self.timing.timestep_fs > 0 else 2.0
                    for i in range(num_points):
                        time_array_ns.append(cumulative_time_ns + i * timestep_fs / 1_000_000.0)
                    cumulative_time_ns += num_points * timestep_fs / 1_000_000.0
                else:
                    # Distribute points evenly across assigned time
                    if num_points == 1:
                        time_array_ns.append(cumulative_time_ns)
                    else:
                        file_times = np.linspace(cumulative_time_ns, cumulative_time_ns + assigned_time_ns, num_points)
                        time_array_ns.extend(file_times.tolist())
                    
                    cumulative_time_ns += assigned_time_ns
            
            return np.array(time_array_ns)
        else:
            # Use timestep-based calculation
            timestep_fs = self.timing.timestep_fs if self.timing and self.timing.timestep_fs > 0 else 2.0
            timesteps = np.array(self.data['timestep'])
            return timesteps * timestep_fs / 1_000_000.0  # Convert fs to ns
    
    def _normalize_property_name(self, prop_name: str) -> Optional[str]:
        """
        Normalize property name to internal data key (case-insensitive).
        
        Handles various input formats:
        - Display names: "Total Energy", "Temperature", etc.
        - NAMD column names: "TOTAL", "TEMP", "PRESSURE", etc.
        - Short names: "total", "temp", "pressure", etc.
        - Any case variation: "Total", "TOTAL", "total", "ToTaL", etc.
        
        Args:
            prop_name: Property name in any format/case
            
        Returns:
            Normalized internal key (lowercase) or None if not recognized
        """
        # Master mapping of all possible property names to internal keys
        # Using lowercase keys for everything
        property_mappings = {
            # Full display names
            'total energy': 'total',
            'potential energy': 'potential',
            'kinetic energy': 'kinetic',
            'electrostatic energy': 'elect',
            'van der waals energy': 'vdw',
            'bond energy': 'bond',
            'angle energy': 'angle',
            'dihedral energy': 'dihedral',
            'improper energy': 'improper',
            'temperature': 'temp',
            'pressure': 'pressure',
            'volume': 'volume',
            # Short forms (NAMD column names and common usage)
            'total': 'total',
            'potential': 'potential',
            'kinetic': 'kinetic',
            'elect': 'elect',
            'electrostatic': 'elect',
            'vdw': 'vdw',
            'bond': 'bond',
            'angle': 'angle',
            'dihedral': 'dihedral',
            'improper': 'improper',
            'temp': 'temp',
            'pressure': 'pressure',
            'volume': 'volume',
            # Additional aliases
            'pot': 'potential',
            'kin': 'kinetic',
            'elec': 'elect',
            'press': 'pressure',
            'vol': 'volume'
        }
        
        # Normalize input to lowercase for case-insensitive matching
        normalized_input = prop_name.lower().strip()
        
        # Direct lookup
        if normalized_input in property_mappings:
            return property_mappings[normalized_input]
        
        return None
    
    def get_available_properties(self) -> List[str]:
        """
        Get list of available energy properties that can be plotted.
        
        Returns:
            List of property names with units
        """
        available = []
        property_map = {
            'total': 'Total Energy',
            'potential': 'Potential Energy',
            'kinetic': 'Kinetic Energy',
            'elect': 'Electrostatic Energy',
            'vdw': 'Van der Waals Energy',
            'bond': 'Bond Energy',
            'angle': 'Angle Energy',
            'dihedral': 'Dihedral Energy',
            'improper': 'Improper Energy',
            'temp': 'Temperature',
            'pressure': 'Pressure',
            'volume': 'Volume'
        }
        
        for key, name in property_map.items():
            if key in self.data and self.data[key]:
                available.append(name)
        
        return available
    
    def plot_properties(self,
                       properties: Optional[List[str]] = None,
                       separate_plots: bool = False,
                       energy_units: str = "kcal/mol",
                       time_units: str = "ns",
                       pressure_units: str = "atm",
                       temperature_units: str = "K",
                       volume_units: str = "Å³",
                       line_colors: Optional[List[str]] = None,
                       bg_color: str = "#2b2b2b",
                       fig_bg_color: str = "#212121",
                       text_color: str = "Auto",
                       grid_color: Optional[str] = None,
                       show_grid: bool = True,
                       xlim: Optional[tuple] = None,
                       ylim: Optional[tuple] = None,
                       title: Optional[str] = None,
                       xlabel: Optional[str] = None,
                       ylabel: Optional[str] = None,
                       save: Optional[str] = None,
                       save_prefix: Optional[str] = None,
                       show: bool = False,
                       figsize: tuple = (10, 6),
                       dpi: int = 300):
        """
        Plot selected energy properties with full GUI-level customization.
        
        Property names are case-insensitive and support multiple formats:
        - Full names: "Temperature", "Total Energy", "Pressure" (any case)
        - Short names: "TEMP", "TOTAL", "PRESSURE" (any case)
        - Aliases: "pot" (potential), "kin" (kinetic), "elec" (electrostatic)
        
        Args:
            properties: List of property names to plot. If None, plots 4-panel view.
                       Available: "Total Energy", "Potential Energy", "Kinetic Energy",
                                 "Electrostatic Energy", "Van der Waals Energy", "Bond Energy",
                                 "Angle Energy", "Dihedral Energy", "Improper Energy",
                                 "Temperature", "Pressure", "Volume"
                       Note: Property names are case-insensitive, so "TEMP", "temp", 
                             "Temperature" all work.
            separate_plots: If True, create separate plot files for each property
            energy_units: "kcal/mol" or "kJ/mol"
            time_units: "ps", "ns", or "µs"
            pressure_units: "atm", "bar", "Pa", "kPa", "MPa"
            temperature_units: "K", "°C", or "°F"
            volume_units: "Å³", "nm³", "mL", "L"
            line_colors: List of colors for each property line
            bg_color: Plot area background color (hex or "none" for transparent)
            fig_bg_color: Figure border background color (hex or "none")
            text_color: Text/axes color ("Auto", color name, or hex)
            grid_color: Grid line color (None to match text_color)
            show_grid: Show grid lines
            xlim: X-axis limits (min, max)
            ylim: Y-axis limits (min, max)
            title: Plot title (auto-generated if None)
            xlabel: X-axis label (auto-generated if None)
            ylabel: Y-axis label (auto-generated if None)
            save: Filename to save plot (only for single plot or combined)
            save_prefix: Prefix for filenames when separate_plots=True
            show: Display plot interactively
            figsize: Figure size (width, height) in inches
            dpi: Resolution for saved figure
        
        Example:
            >>> # Plot specific properties on same figure (case-insensitive)
            >>> analyzer.plot_properties(
            ...     properties=["TEMP", "total energy", "Pressure"],  # Any case works!
            ...     line_colors=["red", "blue", "green"],
            ...     save="combined.png"
            ... )
            
            >>> # Plot each property separately
            >>> analyzer.plot_properties(
            ...     properties=["Temperature", "Pressure"],
            ...     separate_plots=True,
            ...     save_prefix="energy_"
            ... )
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        if not self.data['timestep']:
            logger.warning("No energy data to plot")
            return
        
        # If no properties specified, use the 4-panel plot
        if properties is None:
            return self.plot_energy(
                energy_units=energy_units,
                time_units=time_units,
                bg_color=bg_color,
                fig_bg_color=fig_bg_color,
                text_color=text_color,
                show_grid=show_grid,
                title=title,
                save=save,
                show=show,
                figsize=figsize,
                dpi=dpi
            )
        
        # Calculate time array with proper file scaling
        time_ns = self._calculate_time_array()
        
        # Convert time units
        if time_units == "ps":
            plot_time = time_ns * 1000.0
        elif time_units == "µs":
            plot_time = time_ns / 1000.0
        else:  # ns
            plot_time = time_ns
        
        # Default line colors
        if line_colors is None:
            line_colors = ['blue', 'red', 'green', 'orange', 'purple', 'brown', 'pink', 'cyan']
        
        # Auto-determine text color
        if text_color == "Auto":
            text_color = self._auto_text_color(bg_color)
        
        # Grid color defaults to text color if not specified
        if grid_color is None:
            grid_color = text_color
        
        if separate_plots:
            # Create separate plot for each property
            for i, prop_name in enumerate(properties):
                self._plot_single_property(
                    prop_name, plot_time, time_units,
                    energy_units, pressure_units, temperature_units, volume_units,
                    line_colors[i % len(line_colors)],
                    bg_color, fig_bg_color, text_color, grid_color, show_grid,
                    xlim, ylim, title, xlabel, ylabel,
                    save_prefix, show, figsize, dpi
                )
        else:
            # Plot all properties on same figure
            self._plot_combined_properties(
                properties, plot_time, time_units,
                energy_units, pressure_units, temperature_units, volume_units,
                line_colors, bg_color, fig_bg_color, text_color, grid_color,
                show_grid, xlim, ylim, title, xlabel, ylabel,
                save, show, figsize, dpi
            )
    
    def _auto_text_color(self, bg_color: str) -> str:
        """Auto-determine text color based on background luminance."""
        if bg_color == "none":
            return 'black'
        try:
            hex_color = bg_color.lstrip('#')
            r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
            luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
            return 'black' if luminance > 0.5 else 'white'
        except:
            return 'white'
    
    def _convert_property_units(self, data_key: str, values: 'np.ndarray',
                                energy_units: str, pressure_units: str,
                                temperature_units: str, volume_units: str) -> tuple:
        """Convert property values to requested units and return (values, unit_label)."""
        import numpy as np
        
        # Energy properties
        if data_key in ['total', 'potential', 'kinetic', 'elect', 'vdw', 'bond', 'angle', 'dihedral', 'improper']:
            if energy_units == "kJ/mol":
                return values * 4.184, "kJ/mol"
            return values, "kcal/mol"
        
        # Temperature
        elif data_key == 'temp':
            if temperature_units == "°C":
                return values - 273.15, "°C"
            elif temperature_units == "°F":
                return (values - 273.15) * 9/5 + 32, "°F"
            return values, "K"
        
        # Pressure
        elif data_key == 'pressure':
            conversions = {
                'atm': 1.0,
                'bar': 1.01325,
                'Pa': 101325.0,
                'kPa': 101.325,
                'MPa': 0.101325
            }
            factor = conversions.get(pressure_units, 1.0)
            return values * factor, pressure_units
        
        # Volume
        elif data_key == 'volume':
            conversions = {
                'Å³': 1.0,
                'nm³': 0.001,
                'mL': 1.66054e-24,
                'L': 1.66054e-27
            }
            factor = conversions.get(volume_units, 1.0)
            return values * factor, volume_units
        
        return values, ""
    
    def _plot_single_property(self, prop_name, plot_time, time_units,
                              energy_units, pressure_units, temperature_units, volume_units,
                              line_color, bg_color, fig_bg_color, text_color, grid_color,
                              show_grid, xlim, ylim, title, xlabel, ylabel,
                              save_prefix, show, figsize, dpi):
        """Plot a single property in its own figure."""
        import matplotlib.pyplot as plt
        import numpy as np
        
        data_key = self._normalize_property_name(prop_name)
        if not data_key or data_key not in self.data or not self.data[data_key]:
            logger.warning(f"Property '{prop_name}' not available")
            return
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set backgrounds
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        if bg_color != "none":
            ax.set_facecolor(bg_color)
        
        # Get data and convert units
        y_data = np.array(self.data[data_key])
        y_data, unit_label = self._convert_property_units(
            data_key, y_data, energy_units, pressure_units, temperature_units, volume_units
        )
        
        # Plot
        ax.plot(plot_time, y_data, color=line_color, linewidth=1.5)
        
        # Styling
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        
        if show_grid:
            ax.grid(True, alpha=0.3, color=grid_color)
        
        # Labels
        ax.set_xlabel(xlabel or f"Time ({time_units})", color=text_color)
        ax.set_ylabel(ylabel or f"{prop_name} ({unit_label})", color=text_color)
        ax.set_title(title or prop_name, color=text_color, fontweight='bold')
        
        # Limits
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        plt.tight_layout()
        
        # Save
        if save_prefix:
            safe_name = prop_name.lower().replace(' ', '_')
            filename = f"{save_prefix}{safe_name}.png"
            plt.savefig(filename, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {filename}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def _plot_combined_properties(self, properties, plot_time, time_units,
                                  energy_units, pressure_units, temperature_units, volume_units,
                                  line_colors, bg_color, fig_bg_color, text_color, grid_color,
                                  show_grid, xlim, ylim, title, xlabel, ylabel,
                                  save, show, figsize, dpi):
        """Plot multiple properties on the same figure."""
        import matplotlib.pyplot as plt
        import numpy as np
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set backgrounds
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        if bg_color != "none":
            ax.set_facecolor(bg_color)
        
        # Plot each property
        lines = []
        labels = []
        
        for i, prop_name in enumerate(properties):
            data_key = self._normalize_property_name(prop_name)
            if not data_key or data_key not in self.data or not self.data[data_key]:
                logger.warning(f"Property '{prop_name}' not available or not recognized")
                continue
            
            y_data = np.array(self.data[data_key])
            y_data, unit_label = self._convert_property_units(
                data_key, y_data, energy_units, pressure_units, temperature_units, volume_units
            )
            
            color = line_colors[i % len(line_colors)]
            line = ax.plot(plot_time, y_data, color=color, linewidth=1.5,
                          marker='o', markersize=2, label=f"{prop_name} ({unit_label})")
            lines.extend(line)
            labels.append(f"{prop_name} ({unit_label})")
        
        # Styling
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        
        if show_grid:
            ax.grid(True, alpha=0.3, color=grid_color)
        
        # Labels
        ax.set_xlabel(xlabel or f"Time ({time_units})", color=text_color)
        ax.set_ylabel(ylabel or "Multiple Properties", color=text_color)
        ax.set_title(title or f"NAMD Analysis - {len(properties)} Properties",
                    color=text_color, fontweight='bold')
        
        # Legend
        if len(properties) > 1:
            legend = ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.setp(legend.get_texts(), color=text_color)
        
        # Limits
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        plt.tight_layout()
        
        # Save
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_energy(self, 
                    energy_units: str = "kcal/mol",
                    time_units: str = "ns",
                    bg_color: str = "#2b2b2b",
                    fig_bg_color: str = "#212121",
                    text_color: str = "Auto",
                    show_grid: bool = True,
                    title: Optional[str] = None,
                    target_temperature: Optional[float] = None,
                    target_pressure: Optional[float] = None,
                    save: Optional[str] = None, show: bool = False, 
                    figsize: tuple = (12, 10), dpi: int = 300):
        """
        Create energy analysis plot with full GUI customization.
        
        Args:
            energy_units: 'kcal/mol' or 'kJ/mol'
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            bg_color: Background color for plot area
            fig_bg_color: Background color for figure border
            text_color: Text/axes color ('Auto' or specific color)
            show_grid: Show grid lines on plots
            title: Main title for the figure (default: auto-generated)
            target_temperature: Target temperature in Kelvin (default: auto-calculated from last 50% of trajectory)
            target_pressure: Target pressure in atm (default: 1.0 atm)
            save: Filename to save plot (e.g., "energy.png")
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        if not self.data['timestep']:
            logger.warning("No energy data to plot")
            return
        
        # Calculate time array with proper file scaling
        time_ns = self._calculate_time_array()
        
        # Convert time units
        if time_units == "ps":
            plot_time = time_ns * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = time_ns / 1000.0  # Convert ns to µs
        else:  # ns
            plot_time = time_ns
        
        # Energy unit conversion factor
        if energy_units == "kJ/mol":
            energy_factor = 4.184  # kcal/mol to kJ/mol
        else:  # kcal/mol
            energy_factor = 1.0
        
        # Create figure with subplots
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        
        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        
        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white'
        
        # Set title
        if len(self.log_files) == 1:
            title_text = title or f'Energy Analysis - {self.log_files[0].name}'
        else:
            title_text = title or f'Energy Analysis - {len(self.log_files)} Files'
        fig.suptitle(title_text, fontsize=14, fontweight='bold', color=text_color)
        
        # Configure all subplots with common settings
        for ax in axes.flat:
            if bg_color != "none":
                ax.set_facecolor(bg_color)
            ax.tick_params(colors=text_color)
            for spine in ax.spines.values():
                spine.set_color(text_color)
        
        # Panel 1: Total Energy
        if self.data['total']:
            energy_data = np.array(self.data['total']) * energy_factor / 1000
            axes[0, 0].plot(plot_time, energy_data, 'b-', linewidth=0.8)
            axes[0, 0].set_xlabel(f'Time ({time_units})', color=text_color)
            axes[0, 0].set_ylabel(f'Total Energy (×10³ {energy_units})', color=text_color)
            axes[0, 0].set_title('Total Energy Convergence', color=text_color)
            if show_grid:
                axes[0, 0].grid(True, alpha=0.3, color=text_color)
        
        # Panel 2: Potential and Kinetic Energy
        if self.data['potential'] and self.data['kinetic']:
            pot_energy = np.array(self.data['potential']) * energy_factor / 1000
            kin_energy = np.array(self.data['kinetic']) * energy_factor / 1000
            axes[0, 1].plot(plot_time, pot_energy, 
                           'r-', linewidth=0.8, label='Potential', alpha=0.8)
            axes[0, 1].plot(plot_time, kin_energy, 
                           'g-', linewidth=0.8, label='Kinetic', alpha=0.8)
            axes[0, 1].set_xlabel(f'Time ({time_units})', color=text_color)
            axes[0, 1].set_ylabel(f'Energy (×10³ {energy_units})', color=text_color)
            axes[0, 1].set_title('Potential and Kinetic Energy', color=text_color)
            legend = axes[0, 1].legend()
            plt.setp(legend.get_texts(), color=text_color)
            if show_grid:
                axes[0, 1].grid(True, alpha=0.3, color=text_color)
        
        # Panel 3: Temperature
        if self.data['temp']:
            temp_array = np.array(self.data['temp'])
            # Use user-provided target or auto-calculate from last 50% of trajectory
            target_temp = target_temperature if target_temperature is not None else 300.0
            axes[1, 0].plot(plot_time, temp_array, 'orange', linewidth=0.8)
            axes[1, 0].axhline(y=target_temp, color=text_color, linestyle='--', 
                              linewidth=1, label=f'Target: {target_temp:.1f} K', alpha=0.7)
            axes[1, 0].set_xlabel(f'Time ({time_units})', color=text_color)
            axes[1, 0].set_ylabel('Temperature (K)', color=text_color)
            axes[1, 0].set_title('Temperature Stability', color=text_color)
            legend = axes[1, 0].legend()
            plt.setp(legend.get_texts(), color=text_color)
            if show_grid:
                axes[1, 0].grid(True, alpha=0.3, color=text_color)
        
        # Panel 4: Pressure (if available)
        if self.data['pressure']:
            # Use user-provided target or default to 1.0 atm
            target_press = target_pressure if target_pressure is not None else 1.0
            axes[1, 1].plot(plot_time, self.data['pressure'], 'purple', linewidth=0.8)
            axes[1, 1].axhline(y=target_press, color=text_color, linestyle='--', 
                              linewidth=1, label=f'Target: {target_press:.1f} atm', alpha=0.7)
            axes[1, 1].set_xlabel(f'Time ({time_units})', color=text_color)
            axes[1, 1].set_ylabel('Pressure (atm)', color=text_color)
            axes[1, 1].set_title('Pressure Fluctuations', color=text_color)
            legend = axes[1, 1].legend()
            plt.setp(legend.get_texts(), color=text_color)
            if show_grid:
                axes[1, 1].grid(True, alpha=0.3, color=text_color)
        else:
            # If no pressure, show energy components
            if self.data['elect'] and self.data['vdw']:
                elect_energy = np.array(self.data['elect']) * energy_factor / 1000
                vdw_energy = np.array(self.data['vdw']) * energy_factor / 1000
                axes[1, 1].plot(plot_time, elect_energy, 
                               'b-', linewidth=0.8, label='Electrostatic', alpha=0.7)
                axes[1, 1].plot(plot_time, vdw_energy, 
                               'r-', linewidth=0.8, label='van der Waals', alpha=0.7)
                axes[1, 1].set_xlabel(f'Time ({time_units})', color=text_color)
                axes[1, 1].set_ylabel(f'Energy (×10³ {energy_units})', color=text_color)
                axes[1, 1].set_title('Energy Components', color=text_color)
                legend = axes[1, 1].legend()
                plt.setp(legend.get_texts(), color=text_color)
                if show_grid:
                    axes[1, 1].grid(True, alpha=0.3, color=text_color)
        
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()

class TrajectoryAnalyzer:
    """
    Easy-to-use wrapper for MD trajectory analysis with built-in plotting.
    
    Requires MDAnalysis to be installed.
    
    Supports all GUI options including:
    - Multiple trajectory files with time scaling
    - Full plot customization (colors, units, labels, limits)
    - RMSF-specific X-axis formatting with residue labels
    - Unit conversions (Å/nm, ps/ns/µs, kcal/kJ)
    
    Example:
        >>> analyzer = TrajectoryAnalyzer("system.psf", "trajectory.dcd")
        >>> analyzer.plot_rmsd(selection="protein and backbone", save="rmsd.png")
        >>> 
        >>> # Multiple trajectories with time scaling
        >>> analyzer = TrajectoryAnalyzer(
        ...     "system.psf", 
        ...     ["eq1.dcd", "eq2.dcd", "prod.dcd"],
        ...     file_times={"eq1.dcd": 1.0, "eq2.dcd": 2.0, "prod.dcd": 10.0}  # durations in ns
        ... )
        >>> analyzer.plot_rmsd(selection="protein", time_units="ns", distance_units="Å")
        >>> 
        >>> # Customized plotting
        >>> analyzer.plot_rmsf(
        ...     selection="protein and name CA",
        ...     xaxis_type="residue_type_number",  # Show "ALA123" style labels
        ...     residue_name_format="triple",      # Use 3-letter codes
        ...     label_frequency="every_5",         # Label every 5th residue
        ...     line_color="#1f77b4",
        ...     bg_color="#2b2b2b",
        ...     show_grid=True,
        ...     save="rmsf.png"
        ... )
    """
    
    def __init__(self, topology: Path, 
                 trajectory: Path | List[Path], 
                 file_times: Optional[Dict[str, float]] = None):
        """
        Initialize trajectory analyzer.
        
        Args:
            topology: Path to topology file (PSF, PDB, PRMTOP, etc.)
            trajectory: Path(s) to trajectory file(s) (DCD, XTC, TRR, etc.)
                       Can be a single path or list of paths for concatenated analysis
            file_times: Optional dictionary mapping trajectory filenames to their
                       simulation durations in nanoseconds. Used for proper time scaling.
                       Example: {"eq1.dcd": 1.0, "eq2.dcd": 2.0, "prod.dcd": 10.0}
        """
        try:
            import MDAnalysis as mda
        except ImportError:
            raise ImportError(
                "MDAnalysis is required for trajectory analysis. "
                "Install with: conda install -c conda-forge mdanalysis"
            )
        
        self.topology = Path(topology)
        
        # Handle single or multiple trajectories
        if isinstance(trajectory, (str, Path)):
            self.trajectories = [Path(trajectory)]
        else:
            self.trajectories = [Path(t) for t in trajectory]
        
        # Store file times for proper time scaling
        self.file_times = file_times or {}
        
        # Load trajectories into MDAnalysis
        if len(self.trajectories) == 1:
            self.universe = mda.Universe(str(self.topology), str(self.trajectories[0]))
        else:
            # Concatenate multiple trajectories
            self.universe = mda.Universe(str(self.topology), 
                                        [str(t) for t in self.trajectories])
        
        logger.info(f"Loaded trajectory: {len(self.universe.trajectory)} frames "
                   f"from {len(self.trajectories)} file(s)")
    
    def _calculate_time_array(self) -> 'np.ndarray':
        """
        Calculate proper time array based on file_times dict.
        
        Returns:
            Array of time values in nanoseconds
        """
        try:
            import numpy as np
            import MDAnalysis as mda
        except ImportError:
            raise ImportError("numpy and MDAnalysis are required")
        
        if not self.file_times:
            # Fallback: use frame indices with default timestep (2 fs)
            timestep_ps = 0.002  # 2 fs in ps
            n_frames = len(self.universe.trajectory)
            return np.arange(n_frames) * timestep_ps / 1000.0  # Convert to ns
        
        # Calculate time array based on user-specified file durations
        time_array = []
        cumulative_time_ns = 0.0
        
        # Load each trajectory separately to get frame counts
        for traj_path in self.trajectories:
            # Get filename for lookup in file_times dict
            filename = str(traj_path.name)
            
            # Load this trajectory to count frames
            temp_universe = mda.Universe(str(self.topology), str(traj_path))
            n_frames = len(temp_universe.trajectory)
            
            # Get duration for this file (in ns)
            duration_ns = self.file_times.get(filename, 0.0)
            
            if duration_ns > 0 and n_frames > 1:
                # Create linearly spaced time points for this trajectory
                file_times = np.linspace(cumulative_time_ns, 
                                        cumulative_time_ns + duration_ns, 
                                        n_frames)
            else:
                # Fallback: use frame indices (assume 0.01 ns = 10 ps per frame)
                file_times = cumulative_time_ns + np.arange(n_frames) * 0.01
                if duration_ns > 0:
                    cumulative_time_ns += duration_ns
                else:
                    cumulative_time_ns += n_frames * 0.01
                continue
            
            time_array.extend(file_times)
            cumulative_time_ns += duration_ns
        
        return np.array(time_array)
    
    def calculate_rmsd(self, selection: str = "protein and backbone", 
                       reference_frame: int = 0,
                       align: bool = True) -> Dict[str, 'np.ndarray']:
        """
        Calculate RMSD for selected atoms.
        
        Args:
            selection: MDAnalysis selection string
            reference_frame: Frame to use as reference (0 = first frame)
            align: If True, perform alignment (rotation + translation) before RMSD.
                  If False, calculate raw coordinate RMSD without alignment.
            
        Returns:
            Dictionary with 'time' (ns) and 'rmsd' (Angstroms) arrays
        """
        try:
            import MDAnalysis as mda
            from MDAnalysis.analysis import rms
            from MDAnalysis.analysis import align as mda_align
            import numpy as np
        except ImportError:
            raise ImportError("MDAnalysis and numpy are required")
        
        # Select atoms
        atoms = self.universe.select_atoms(selection)
        
        # Set reference frame
        self.universe.trajectory[reference_frame]
        ref_coords = self.universe.select_atoms(selection).positions.copy()
        
        if align:
            # Perform alignment using MDAnalysis align module
            # This modifies the trajectory in-place
            aligner = mda_align.AlignTraj(self.universe, self.universe,
                                         select=selection,
                                         ref_frame=reference_frame,
                                         in_memory=True)
            aligner.run()
            logger.info("Alignment completed (rotation + translation applied)")
        
        # Calculate RMSD for each frame
        rmsd_values = []
        for ts in self.universe.trajectory:
            if align:
                # Structures already aligned, compute RMSD directly
                rmsd = rms.rmsd(atoms.positions, ref_coords, superposition=False)
            else:
                # Raw RMSD without any alignment
                diff = atoms.positions - ref_coords
                rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
            rmsd_values.append(rmsd)
        
        rmsd_array = np.array(rmsd_values)
        
        # Get proper time array (in nanoseconds)
        time_ns = self._calculate_time_array()
        
        return {
            'time': time_ns,
            'rmsd': rmsd_array  # RMSD values in Angstroms
        }
    
    def calculate_rmsf(self, selection: str = "protein and name CA") -> Dict[str, 'np.ndarray']:
        """
        Calculate RMSF for selected atoms.
        
        Args:
            selection: MDAnalysis selection string
            
        Returns:
            Dictionary with 'resids', 'rmsf' (Angstroms), 'resnames', and 'atom_indices' arrays
        """
        try:
            import MDAnalysis as mda
            from MDAnalysis.analysis import rms
            import numpy as np
        except ImportError:
            raise ImportError("MDAnalysis and numpy are required")
        
        # Select atoms
        atoms = self.universe.select_atoms(selection)
        
        # Calculate RMSF
        rmsf_analysis = rms.RMSF(atoms).run()
        
        return {
            'resids': atoms.resids,
            'rmsf': rmsf_analysis.results.rmsf,  # RMSF in Angstroms
            'resnames': atoms.resnames,  # Residue names (e.g., ALA, GLY)
            'atom_indices': atoms.indices  # Atom indices
        }
    
    def calculate_distances(self, selections: Dict[str, tuple]) -> Dict[str, Dict[str, 'np.ndarray']]:
        """
        Calculate distances between atom selections.
        
        Args:
            selections: Dictionary of {name: (selection1, selection2)}
            
        Example:
            selections = {
                "gate": ("resid 50-70 and name CA", "resid 150-170 and name CA"),
                "salt_bridge": ("resid 125 and name NH1", "resid 200 and name OD1 OD2")
            }
            
        Returns:
            Dictionary with distance data for each named selection.
            Each entry contains 'time' (ns) and 'distance' (Angstroms) arrays.
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("numpy is required")
        
        results = {}
        
        for name, (sel1, sel2) in selections.items():
            atoms1 = self.universe.select_atoms(sel1)
            atoms2 = self.universe.select_atoms(sel2)
            
            distances = []
            
            for ts in self.universe.trajectory:
                # Calculate center of mass distance in Angstroms
                dist = np.linalg.norm(
                    atoms1.center_of_mass() - atoms2.center_of_mass()
                )
                distances.append(dist)
            
            # Get proper time array (in nanoseconds)
            time_ns = self._calculate_time_array()
            
            results[name] = {
                'time': time_ns,
                'distance': np.array(distances)  # Distance in Angstroms
            }
        
        return results
    
    def calculate_radius_of_gyration(self, selection: str = "protein") -> Dict[str, 'np.ndarray']:
        """
        Calculate radius of gyration over trajectory.
        
        Args:
            selection: MDAnalysis selection string
            
        Returns:
            Dictionary with 'time' (ns) and 'rg' (Angstroms) arrays
        """
        try:
            import numpy as np
        except ImportError:
            raise ImportError("numpy is required")
        
        atoms = self.universe.select_atoms(selection)
        
        rg_values = []
        
        for ts in self.universe.trajectory:
            rg_values.append(atoms.radius_of_gyration())  # In Angstroms
        
        # Get proper time array (in nanoseconds)
        time_ns = self._calculate_time_array()
        
        return {
            'time': time_ns,
            'rg': np.array(rg_values)  # Rg in Angstroms
        }
    
    def plot_rmsd(self, selection: str = "protein and backbone", 
                  reference_frame: int = 0,
                  align: bool = True,
                  distance_units: str = "Å",
                  time_units: str = "ns",
                  line_color: str = "blue",
                  line_width: float = 1.2,
                  line_style: str = "-",
                  bg_color: str = "#2b2b2b",
                  fig_bg_color: str = "#212121",
                  text_color: str = "Auto",
                  show_grid: bool = True,
                  xlim: Optional[tuple] = None,
                  ylim: Optional[tuple] = None,
                  title: Optional[str] = None,
                  xlabel: Optional[str] = None,
                  ylabel: Optional[str] = None,
                  highlight_threshold: Optional[float] = None,
                  highlight_color: str = "orange",
                  highlight_alpha: float = 0.2,
                  show_convergence: bool = True,
                  convergence_color: str = "red",
                  convergence_style: str = "--",
                  convergence_width: float = 1.5,
                  hlines: Optional[List[float]] = None,
                  hline_colors: Optional[List[str]] = None,
                  hline_styles: Optional[List[str]] = None,
                  hline_widths: Optional[List[float]] = None,
                  vlines: Optional[List[float]] = None,
                  vline_colors: Optional[List[str]] = None,
                  vline_styles: Optional[List[str]] = None,
                  vline_widths: Optional[List[float]] = None,
                  save: Optional[str] = None, show: bool = False,
                  figsize: tuple = (10, 6), dpi: int = 300):
        """
        Plot RMSD with full GUI customization options.
        
        Args:
            selection: MDAnalysis selection string
            reference_frame: Frame to use as reference
            align: Perform alignment before RMSD calculation
            distance_units: 'Å' (angstrom) or 'nm' (nanometer)
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            line_color: Color for the plot line (matplotlib color string or hex)
            line_width: Width of the plot line (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area (hex or 'none' for transparent)
            fig_bg_color: Background color for figure border (hex or 'none')
            text_color: Text/axes color ('Auto', matplotlib color, or hex)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: auto-generated)
            xlabel: X-axis label (default: auto-generated with units)
            ylabel: Y-axis label (default: auto-generated with units)
            highlight_threshold: If set, highlight regions above this RMSD value
            highlight_color: Color for highlight region and line (default: 'orange')
            highlight_alpha: Alpha transparency for highlight fill (default: 0.2)
            show_convergence: Show convergence line (mean of last 20% of trajectory)
            convergence_color: Color for convergence line (default: 'red')
            convergence_style: Line style for convergence line (default: '--')
            convergence_width: Width of convergence line (default: 1.5)
            hlines: List of Y values for horizontal reference lines
            hline_colors: List of colors for horizontal lines (default: cycle through standard colors)
            hline_styles: List of line styles for horizontal lines (default: '--')
            hline_widths: List of line widths for horizontal lines (default: 1.0)
            vlines: List of X values for vertical reference lines
            vline_colors: List of colors for vertical lines (default: cycle through standard colors)
            vline_styles: List of line styles for vertical lines (default: '--')
            vline_widths: List of line widths for vertical lines (default: 1.0)
            save: Filename to save plot
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        data = self.calculate_rmsd(selection, reference_frame, align)
        
        # Convert units
        plot_time = data['time'].copy()  # Time is in ns
        plot_rmsd = data['rmsd'].copy()  # RMSD is in Angstroms
        
        # Convert distance units
        if distance_units == "nm":
            plot_rmsd = plot_rmsd / 10.0  # Convert Å to nm
        
        # Convert time units
        if time_units == "ps":
            plot_time = plot_time * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = plot_time / 1000.0  # Convert ns to µs
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        
        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)
        
        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                # Calculate luminance
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white'
        
        # Plot data
        ax.plot(plot_time, plot_rmsd, color=line_color, linewidth=line_width, 
                linestyle=line_style, alpha=0.7)
        
        # Add highlight threshold if specified
        if highlight_threshold is not None:
            threshold_display = highlight_threshold if distance_units == "Å" else highlight_threshold / 10.0
            ax.axhline(y=threshold_display, color=highlight_color, linestyle=':', 
                      linewidth=1.5, alpha=0.8, 
                      label=f'Threshold: {highlight_threshold} {distance_units}')
            # Fill region above threshold
            # Convert boolean array to list for type compatibility
            where_condition = (plot_rmsd >= threshold_display).tolist()  # type: ignore
            ax.fill_between(plot_time, threshold_display, plot_rmsd, 
                          where=where_condition,
                          alpha=highlight_alpha, color=highlight_color, label='Above threshold')
        
        # Add convergence line (last 20% of trajectory)
        if show_convergence:
            cutoff_idx = int(len(plot_rmsd) * 0.8)
            converged_value = float(np.mean(plot_rmsd[cutoff_idx:]))
            ax.axhline(y=converged_value, color=convergence_color, 
                      linestyle=convergence_style, linewidth=convergence_width,
                      label=f'Converged: {converged_value:.2f} {distance_units}')
        
        # Add custom horizontal reference lines
        if hlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, yval in enumerate(hlines):
                color = hline_colors[i] if hline_colors and i < len(hline_colors) else default_colors[i % len(default_colors)]
                style = hline_styles[i] if hline_styles and i < len(hline_styles) else '--'
                width = hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                ax.axhline(y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Add custom vertical reference lines
        if vlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, xval in enumerate(vlines):
                color = vline_colors[i] if vline_colors and i < len(vline_colors) else default_colors[i % len(default_colors)]
                style = vline_styles[i] if vline_styles and i < len(vline_styles) else '--'
                width = vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                ax.axvline(x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Set labels with appropriate color
        xlabel_text = xlabel or f'Time ({time_units})'
        ylabel_text = ylabel or f'RMSD ({distance_units})'
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)
        
        # Set title
        title_text = title or f'RMSD - {selection}'
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight='bold')
        
        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        
        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color)
        
        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        ax.legend()
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_rmsf(self, selection: str = "protein and name CA",
                  xaxis_type: str = "residue_number",
                  show_residue_labels: bool = True,
                  residue_name_format: str = "single",
                  label_frequency: str = "auto",
                  distance_units: str = "Å",
                  line_color: str = "blue",
                  line_width: float = 1.2,
                  line_style: str = "-",
                  bg_color: str = "#2b2b2b",
                  fig_bg_color: str = "#212121",
                  text_color: str = "Auto",
                  show_grid: bool = True,
                  xlim: Optional[tuple] = None,
                  ylim: Optional[tuple] = None,
                  title: Optional[str] = None,
                  xlabel: Optional[str] = None,
                  ylabel: Optional[str] = None,
                  highlight_threshold: Optional[float] = None,
                  highlight_color: str = "orange",
                  highlight_alpha: float = 0.2,
                  hlines: Optional[List[float]] = None,
                  hline_colors: Optional[List[str]] = None,
                  hline_styles: Optional[List[str]] = None,
                  hline_widths: Optional[List[float]] = None,
                  vlines: Optional[List[float]] = None,
                  vline_colors: Optional[List[str]] = None,
                  vline_styles: Optional[List[str]] = None,
                  vline_widths: Optional[List[float]] = None,
                  save: Optional[str] = None, show: bool = False,
                  figsize: tuple = (12, 6), dpi: int = 300):
        """
        Plot RMSF with full GUI customization including residue labeling.
        
        Args:
            selection: MDAnalysis selection string
            xaxis_type: X-axis type - 'residue_number', 'residue_type_number', or 'atom_index'
            show_residue_labels: Show residue labels on X-axis
            residue_name_format: 'single' (A, G, V) or 'triple' (ALA, GLY, VAL)
            label_frequency: 'all', 'auto', 'every_2', 'every_5', 'every_10', 'every_20'
            distance_units: 'Å' (angstrom) or 'nm' (nanometer) for Y-axis
            line_color: Color for the plot line
            line_width: Width of the plot line (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area
            fig_bg_color: Background color for figure border
            text_color: Text/axes color ('Auto' or specific color)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: auto-generated)
            xlabel: X-axis label (default: auto-generated)
            ylabel: Y-axis label (default: auto-generated with units)
            highlight_threshold: If set, highlight residues above this RMSF value
            highlight_color: Color for highlight region and line (default: 'orange')
            highlight_alpha: Alpha transparency for highlight fill (default: 0.2)
            hlines: List of Y values for horizontal reference lines
            hline_colors: List of colors for horizontal lines
            hline_styles: List of line styles for horizontal lines (default: '--')
            hline_widths: List of line widths for horizontal lines (default: 1.0)
            vlines: List of X values for vertical reference lines
            vline_colors: List of colors for vertical lines
            vline_styles: List of line styles for vertical lines (default: '--')
            vline_widths: List of line widths for vertical lines (default: 1.0)
            save: Filename to save plot
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        data = self.calculate_rmsf(selection)
        
        # Convert distance units
        plot_rmsf = data['rmsf'].copy()  # RMSF is in Angstroms
        if distance_units == "nm":
            plot_rmsf = plot_rmsf / 10.0  # Convert Å to nm
        
        # Prepare X-axis data and labels
        if xaxis_type == "residue_number":
            x_data = data['resids']
            xlabel_default = "Residue Number"
        elif xaxis_type == "residue_type_number":
            x_data = data['resids']
            xlabel_default = "Residue"
        elif xaxis_type == "atom_index":
            x_data = data['atom_indices']
            xlabel_default = "Atom Index"
        else:
            x_data = data['resids']
            xlabel_default = "Residue Number"
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        
        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)
        
        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white'
        
        # Plot data
        ax.plot(x_data, plot_rmsf, color=line_color, linewidth=line_width, linestyle=line_style)
        ax.fill_between(x_data, 0, plot_rmsf, alpha=0.3, color=line_color)
        
        # Highlight flexible regions if threshold specified
        if highlight_threshold is not None:
            threshold_display = highlight_threshold if distance_units == "Å" else highlight_threshold / 10.0
            ax.axhline(y=threshold_display, color=highlight_color, linestyle=':', 
                      linewidth=1.5, alpha=0.8, 
                      label=f'Threshold: {highlight_threshold} {distance_units}')
            # Fill region above threshold
            # Convert boolean array to list for type compatibility
            where_condition = (plot_rmsf >= threshold_display).tolist()  # type: ignore
            ax.fill_between(x_data, threshold_display, plot_rmsf, 
                          where=where_condition,
                          alpha=highlight_alpha, color=highlight_color, label='Above threshold')
        
        # Add custom horizontal reference lines
        if hlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, yval in enumerate(hlines):
                color = hline_colors[i] if hline_colors and i < len(hline_colors) else default_colors[i % len(default_colors)]
                style = hline_styles[i] if hline_styles and i < len(hline_styles) else '--'
                width = hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                ax.axhline(y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Add custom vertical reference lines
        if vlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, xval in enumerate(vlines):
                color = vline_colors[i] if vline_colors and i < len(vline_colors) else default_colors[i % len(default_colors)]
                style = vline_styles[i] if vline_styles and i < len(vline_styles) else '--'
                width = vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                ax.axvline(x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Set labels
        xlabel_text = xlabel or xlabel_default
        ylabel_text = ylabel or f'RMSF ({distance_units})'
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)
        
        # Set title
        title_text = title or f'RMSF - {selection}'
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight='bold')
        
        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        
        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color, axis='y')
        
        # Handle residue labels on X-axis
        if show_residue_labels and xaxis_type == "residue_type_number":
            # Create residue labels with names
            resnames = data['resnames']
            resids = data['resids']
            
            # Convert residue names if needed
            if residue_name_format == "single":
                # Use 1-letter amino acid codes
                aa_codes = {
                    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
                }
                labels = [f"{aa_codes.get(name, name)}{resid}" 
                         for name, resid in zip(resnames, resids)]
            else:  # triple
                labels = [f"{name}{resid}" for name, resid in zip(resnames, resids)]
            
            # Determine label frequency
            n_residues = len(resids)
            if label_frequency == "all":
                step = 1
            elif label_frequency == "auto":
                # Auto-determine based on number of residues
                if n_residues < 20:
                    step = 1
                elif n_residues < 50:
                    step = 2
                elif n_residues < 100:
                    step = 5
                elif n_residues < 200:
                    step = 10
                else:
                    step = 20
            elif label_frequency == "every_2":
                step = 2
            elif label_frequency == "every_5":
                step = 5
            elif label_frequency == "every_10":
                step = 10
            elif label_frequency == "every_20":
                step = 20
            else:
                step = 1
            
            # Set tick positions and labels
            tick_positions = resids[::step]
            tick_labels = [labels[i] for i in range(0, len(labels), step)]
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels, rotation=45, ha='right')
        
        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        if highlight_threshold is not None:
            ax.legend()
        
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_distances(self, selections: Dict[str, tuple],
                       distance_units: str = "Å",
                       time_units: str = "ns",
                       line_colors: Optional[List[str]] = None,
                       line_width: float = 1.2,
                       line_style: str = "-",
                       bg_color: str = "#2b2b2b",
                       fig_bg_color: str = "#212121",
                       text_color: str = "Auto",
                       show_grid: bool = True,
                       xlim: Optional[tuple] = None,
                       ylim: Optional[tuple] = None,
                       title: Optional[str] = None,
                       xlabel: Optional[str] = None,
                       ylabel: Optional[str] = None,
                       show_mean_lines: bool = True,
                       hlines: Optional[List[float]] = None,
                       hline_colors: Optional[List[str]] = None,
                       hline_styles: Optional[List[str]] = None,
                       hline_widths: Optional[List[float]] = None,
                       vlines: Optional[List[float]] = None,
                       vline_colors: Optional[List[str]] = None,
                       vline_styles: Optional[List[str]] = None,
                       vline_widths: Optional[List[float]] = None,
                       save: Optional[str] = None, show: bool = False,
                       figsize: tuple = (10, 6), dpi: int = 300):
        """
        Plot distances with full GUI customization.
        
        Args:
            selections: Dictionary of {name: (selection1, selection2)}
            distance_units: 'Å' (angstrom) or 'nm' (nanometer)
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            line_colors: List of colors for each distance pair (default: auto-cycle)
            line_width: Width of the plot lines (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area
            fig_bg_color: Background color for figure border
            text_color: Text/axes color ('Auto' or specific color)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: "Distance Analysis")
            xlabel: X-axis label (default: auto-generated with units)
            ylabel: Y-axis label (default: auto-generated with units)
            show_mean_lines: Show mean distance lines
            save: Filename to save plot
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        results = self.calculate_distances(selections)
        
        # Default colors if not provided
        if line_colors is None:
            line_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        
        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)
        
        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white'
        
        # Plot each distance pair
        for i, (name, data) in enumerate(results.items()):
            color = line_colors[i % len(line_colors)]
            
            # Convert units
            plot_time = data['time'].copy()  # Time is in ns
            plot_distance = data['distance'].copy()  # Distance is in Å
            
            # Convert distance units
            if distance_units == "nm":
                plot_distance = plot_distance / 10.0  # Convert Å to nm
            
            # Convert time units
            if time_units == "ps":
                plot_time = plot_time * 1000.0  # Convert ns to ps
            elif time_units == "µs":
                plot_time = plot_time / 1000.0  # Convert ns to µs
            
            ax.plot(plot_time, plot_distance, color=color, 
                   linewidth=line_width, linestyle=line_style, label=name, alpha=0.7)
            
            # Add mean line (calculated from second half of trajectory)
            if show_mean_lines:
                mean_dist = float(np.mean(plot_distance[int(len(plot_distance)*0.5):]))
                ax.axhline(y=mean_dist, color=color, linestyle='--', 
                          linewidth=1.0, alpha=0.5)
        
        # Set labels
        xlabel_text = xlabel or f'Time ({time_units})'
        ylabel_text = ylabel or f'Distance ({distance_units})'
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)
        
        # Set title
        title_text = title or 'Distance Analysis'
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight='bold')
        
        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        
        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color)
        
        # Add custom horizontal reference lines
        if hlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, yval in enumerate(hlines):
                color = hline_colors[i] if hline_colors and i < len(hline_colors) else default_colors[i % len(default_colors)]
                style = hline_styles[i] if hline_styles and i < len(hline_styles) else '--'
                width = hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                ax.axhline(y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Add custom vertical reference lines
        if vlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, xval in enumerate(vlines):
                color = vline_colors[i] if vline_colors and i < len(vline_colors) else default_colors[i % len(default_colors)]
                style = vline_styles[i] if vline_styles and i < len(vline_styles) else '--'
                width = vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                ax.axvline(x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        ax.legend()
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_radius_of_gyration(self, selection: str = "protein",
                                distance_units: str = "Å",
                                time_units: str = "ns",
                                line_color: str = "purple",
                                line_width: float = 1.2,
                                line_style: str = "-",
                                bg_color: str = "#2b2b2b",
                                fig_bg_color: str = "#212121",
                                text_color: str = "Auto",
                                show_grid: bool = True,
                                xlim: Optional[tuple] = None,
                                ylim: Optional[tuple] = None,
                                title: Optional[str] = None,
                                xlabel: Optional[str] = None,
                                ylabel: Optional[str] = None,
                                show_convergence: bool = True,
                                convergence_color: str = "red",
                                convergence_style: str = "--",
                                convergence_width: float = 1.5,
                                hlines: Optional[List[float]] = None,
                                hline_colors: Optional[List[str]] = None,
                                hline_styles: Optional[List[str]] = None,
                                hline_widths: Optional[List[float]] = None,
                                vlines: Optional[List[float]] = None,
                                vline_colors: Optional[List[str]] = None,
                                vline_styles: Optional[List[str]] = None,
                                vline_widths: Optional[List[float]] = None,
                                save: Optional[str] = None, show: bool = False,
                                figsize: tuple = (10, 6), dpi: int = 300):
        """
        Plot radius of gyration with full GUI customization options.
        
        Args:
            selection: MDAnalysis selection string (default: "protein")
            distance_units: 'Å' (angstrom) or 'nm' (nanometer)
            time_units: 'ps' (picoseconds), 'ns' (nanoseconds), or 'µs' (microseconds)
            line_color: Color for the plot line (matplotlib color string or hex)
            line_width: Width of the plot line (default: 1.2)
            line_style: Line style: '-' (solid), '--' (dashed), '-.' (dash-dot), ':' (dotted)
            bg_color: Background color for plot area (hex or 'none' for transparent)
            fig_bg_color: Background color for figure border (hex or 'none')
            text_color: Text/axes color ('Auto', matplotlib color, or hex)
            show_grid: Show grid lines on plot
            xlim: X-axis limits as (min, max) tuple
            ylim: Y-axis limits as (min, max) tuple
            title: Plot title (default: auto-generated)
            xlabel: X-axis label (default: auto-generated with units)
            ylabel: Y-axis label (default: auto-generated with units)
            show_convergence: Show convergence line (mean of last 20% of trajectory)
            convergence_color: Color for convergence line (default: 'red')
            convergence_style: Line style for convergence line (default: '--')
            convergence_width: Width of convergence line (default: 1.5)
            save: Filename to save plot
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        data = self.calculate_radius_of_gyration(selection)
        
        # Convert units
        plot_time = data['time'].copy()  # Time is in ns
        plot_rg = data['rg'].copy()  # Rg is in Angstroms
        
        # Convert distance units
        if distance_units == "nm":
            plot_rg = plot_rg / 10.0  # Convert Å to nm
        
        # Convert time units
        if time_units == "ps":
            plot_time = plot_time * 1000.0  # Convert ns to ps
        elif time_units == "µs":
            plot_time = plot_time / 1000.0  # Convert ns to µs
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Set figure background
        if fig_bg_color != "none":
            fig.patch.set_facecolor(fig_bg_color)
        
        # Set plot background
        if bg_color != "none":
            ax.set_facecolor(bg_color)
        
        # Auto-determine text color if needed
        if text_color == "Auto":
            if bg_color == "none":
                text_color = 'black'
            else:
                # Calculate luminance
                try:
                    hex_color = bg_color.lstrip('#')
                    r, g, b = int(hex_color[0:2], 16), int(hex_color[2:4], 16), int(hex_color[4:6], 16)
                    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255
                    text_color = 'black' if luminance > 0.5 else 'white'
                except:
                    text_color = 'white'
        
        # Plot data
        ax.plot(plot_time, plot_rg, color=line_color, linewidth=line_width, 
                linestyle=line_style, alpha=0.7)
        
        # Add convergence line (last 20% of trajectory)
        if show_convergence:
            cutoff_idx = int(len(plot_rg) * 0.8)
            converged_value = float(np.mean(plot_rg[cutoff_idx:]))
            ax.axhline(y=converged_value, color=convergence_color, 
                      linestyle=convergence_style, linewidth=convergence_width,
                      label=f'Converged: {converged_value:.2f} {distance_units}')
        
        # Add custom horizontal reference lines
        if hlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, yval in enumerate(hlines):
                color = hline_colors[i] if hline_colors and i < len(hline_colors) else default_colors[i % len(default_colors)]
                style = hline_styles[i] if hline_styles and i < len(hline_styles) else '--'
                width = hline_widths[i] if hline_widths and i < len(hline_widths) else 1.0
                ax.axhline(y=yval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Add custom vertical reference lines
        if vlines:
            default_colors = ['gray', 'darkgray', 'lightgray', 'silver']
            for i, xval in enumerate(vlines):
                color = vline_colors[i] if vline_colors and i < len(vline_colors) else default_colors[i % len(default_colors)]
                style = vline_styles[i] if vline_styles and i < len(vline_styles) else '--'
                width = vline_widths[i] if vline_widths and i < len(vline_widths) else 1.0
                ax.axvline(x=xval, color=color, linestyle=style, linewidth=width, alpha=0.7)
        
        # Set labels with appropriate color
        xlabel_text = xlabel or f'Time ({time_units})'
        ylabel_text = ylabel or f'Radius of Gyration ({distance_units})'
        ax.set_xlabel(xlabel_text, color=text_color, fontsize=12)
        ax.set_ylabel(ylabel_text, color=text_color, fontsize=12)
        
        # Set title
        title_text = title or f'Radius of Gyration - {selection}'
        ax.set_title(title_text, color=text_color, fontsize=14, fontweight='bold')
        
        # Configure axes colors
        ax.tick_params(colors=text_color)
        for spine in ax.spines.values():
            spine.set_color(text_color)
        
        # Grid
        if show_grid:
            ax.grid(True, alpha=0.3, color=text_color)
        
        # Set axis limits if specified
        if xlim:
            ax.set_xlim(xlim)
        if ylim:
            ax.set_ylim(ylim)
        
        if show_convergence:
            ax.legend()
        
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_summary(self, save: Optional[str] = None, show: bool = False,
                          figsize: tuple = (14, 10), dpi: int = 300):
        """
        Create summary analysis plot with RMSD, RMSF, and Rg.
        
        Args:
            save: Filename to save plot
            show: Whether to display plot interactively
            figsize: Figure size (width, height)
            dpi: Resolution for saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import numpy as np
        except ImportError:
            logger.error("matplotlib and numpy are required for plotting")
            return
        
        # Calculate all metrics
        rmsd_data = self.calculate_rmsd("protein and backbone")
        rmsf_data = self.calculate_rmsf("protein and name CA")
        rg_data = self.calculate_radius_of_gyration("protein")
        
        # Create figure
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        fig.suptitle('Trajectory Analysis Summary', 
                     fontsize=14, fontweight='bold')
        
        # RMSD
        axes[0, 0].plot(rmsd_data['time'], rmsd_data['rmsd'], 'b-', linewidth=1.2, alpha=0.7)
        cutoff = int(len(rmsd_data['rmsd']) * 0.8)
        converged = np.mean(rmsd_data['rmsd'][cutoff:])
        axes[0, 0].axhline(y=converged, color='r', linestyle='--', linewidth=1.5,
                          label=f'Converged: {converged:.2f} Å')
        axes[0, 0].set_xlabel('Time (ns)')
        axes[0, 0].set_ylabel('RMSD (Å)')
        axes[0, 0].set_title('RMSD - Protein Backbone')
        axes[0, 0].legend()
        axes[0, 0].grid(True, alpha=0.3)
        
        # RMSF
        axes[0, 1].plot(rmsf_data['resids'], rmsf_data['rmsf'], 'g-', linewidth=1.2)
        axes[0, 1].fill_between(rmsf_data['resids'], 0, rmsf_data['rmsf'], alpha=0.3, color='g')
        axes[0, 1].set_xlabel('Residue Number')
        axes[0, 1].set_ylabel('RMSF (Å)')
        axes[0, 1].set_title('RMSF - C-alpha Atoms')
        axes[0, 1].grid(True, alpha=0.3, axis='y')
        
        # Radius of Gyration
        axes[1, 0].plot(rg_data['time'], rg_data['rg'], 'purple', linewidth=1.2, alpha=0.7)
        cutoff = int(len(rg_data['rg']) * 0.8)
        converged_rg = np.mean(rg_data['rg'][cutoff:])
        axes[1, 0].axhline(y=converged_rg, color='r', linestyle='--', linewidth=1.5,
                          label=f'Converged: {converged_rg:.2f} Å')
        axes[1, 0].set_xlabel('Time (ns)')
        axes[1, 0].set_ylabel('Rg (Å)')
        axes[1, 0].set_title('Radius of Gyration')
        axes[1, 0].legend()
        axes[1, 0].grid(True, alpha=0.3)
        
        # Statistics summary
        axes[1, 1].axis('off')
        stats_text = f"""
Statistics Summary:

RMSD:
  Mean: {np.mean(rmsd_data['rmsd']):.2f} ± {np.std(rmsd_data['rmsd']):.2f} Å
  Converged: {converged:.2f} Å

RMSF:
  Mean: {np.mean(rmsf_data['rmsf']):.2f} ± {np.std(rmsf_data['rmsf']):.2f} Å
  Max: {np.max(rmsf_data['rmsf']):.2f} Å (Residue {rmsf_data['resids'][np.argmax(rmsf_data['rmsf'])]})

Radius of Gyration:
  Mean: {np.mean(rg_data['rg']):.2f} ± {np.std(rg_data['rg']):.2f} Å
  Converged: {converged_rg:.2f} Å
  
Trajectory Info:
  Frames: {len(self.universe.trajectory)}
  Time: {rmsd_data['time'][-1]:.2f} ns
"""
        axes[1, 1].text(0.1, 0.5, stats_text, fontsize=10, family='monospace',
                       verticalalignment='center')
        
        plt.tight_layout()
        
        if save:
            plt.savefig(save, dpi=dpi, bbox_inches='tight')
            logger.info(f"Plot saved: {save}")
        
        if show:
            plt.show()
        else:
            plt.close()

