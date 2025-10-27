# gatewizard/core/job_monitor.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Job monitoring module for tracking system preparation progress.

This module provides functionality to monitor running and completed
preparation jobs, track their status, and manage job data.
"""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
from enum import Enum

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)

class JobStatus(Enum):
    """Enumeration of possible job statuses."""
    RUNNING = "running"
    COMPLETED = "completed"
    ERROR = "error"
    UNKNOWN = "unknown"

class JobInfo:
    """Container for job information."""
    
    def __init__(self, job_id: str, job_dir: Path, status_data: Dict[str, Any]):
        self.job_id = job_id
        self.job_dir = job_dir
        self.status_data = status_data
        self._last_updated = time.time()
    
    @property
    def status(self) -> JobStatus:
        """Get current job status."""
        status_str = self.status_data.get("status", "unknown")
        try:
            return JobStatus(status_str)
        except ValueError:
            return JobStatus.UNKNOWN
    
    @property
    def pdb_file(self) -> str:
        """Get PDB file path."""
        return self.status_data.get("pdb_file", "")
    
    @property
    def local_pdb(self) -> str:
        """Get local PDB file path."""
        return self.status_data.get("local_pdb", "")
    
    @property
    def start_time(self) -> Optional[datetime]:
        """Get job start time."""
        start_str = self.status_data.get("start_time")
        if start_str:
            try:
                return datetime.fromisoformat(start_str)
            except ValueError:
                pass
        return None
    
    @property
    def end_time(self) -> Optional[datetime]:
        """Get job end time."""
        end_str = self.status_data.get("end_time")
        if end_str:
            try:
                return datetime.fromisoformat(end_str)
            except ValueError:
                pass
        return None
    
    @property
    def elapsed_time(self) -> Optional[float]:
        """Get elapsed time in seconds."""
        start = self.start_time
        if start:
            # For completed jobs, use end_time. For running jobs, use current time.
            if self.status == JobStatus.COMPLETED and self.end_time:
                end = self.end_time
            else:
                end = datetime.now()
            return (end - start).total_seconds()
        return None
    
    @property
    def current_step(self) -> int:
        """Get current step number."""
        return self.status_data.get("current_step", 0)
    
    @property
    def steps_completed(self) -> List[str]:
        """Get list of completed steps."""
        return self.status_data.get("steps_completed", [])
    
    @property
    def error_message(self) -> Optional[str]:
        """Get error message if job failed."""
        return self.status_data.get("error")
    
    @property
    def command(self) -> str:
        """Get the command that was executed."""
        return self.status_data.get("command", "")
    
    @property
    def progress(self) -> float:
        """Get job progress as a fraction (0.0 to 1.0)."""
        if self.status == JobStatus.COMPLETED:
            return 1.0
        elif self.status == JobStatus.ERROR:
            return 0.0
        else:
            # Dynamic progress calculation based on workflow configuration
            config = self.status_data.get('config', {})
            
            # Determine the expected number of steps based on configuration
            preoriented = config.get('preoriented', True)
            parametrize = config.get('parametrize', True)
            
            # Calculate expected steps
            expected_steps = []
            if not preoriented:
                expected_steps.append("MEMEMBED")
            expected_steps.append("Packmol")
            if parametrize:
                expected_steps.extend(["pdb4amber", "tleap"])
            
            total_steps = len(expected_steps)
            current = self.current_step
            
            # Calculate progress based on actual workflow
            if current >= total_steps:
                return 1.0
            elif current <= 0:
                return 0.05  # Just started
            else:
                # Linear progress based on step completion
                return min(0.95, (current / total_steps))  # Cap at 95% until truly complete
    
    @property
    def step_description(self) -> str:
        """Get a detailed description of the current step."""
        # First, check if we have completed steps to show most recent
        if self.steps_completed:
            last_step = self.steps_completed[-1]
            if last_step in ["Completed", "Done"]:
                return "Completed"
            return f"Last: {last_step}"
        
        # Otherwise use current_step
        current = self.current_step
        descriptions = {
            0: "Starting preparation...",
            1: "Running MEMEMBED orientation...",
            2: "Packing lipids with Packmol...",
            3: "Parametrizing system...",
            4: "Building Amber input files...",
            5: "Completed"
        }
        return descriptions.get(current, f"Step {current}")
    
    def update_status(self, new_status_data: Dict[str, Any]):
        """Update job status data."""
        self.status_data.update(new_status_data)
        self._last_updated = time.time()
    
    def is_stale(self, max_age_seconds: float = 300) -> bool:
        """Check if job data is stale (hasn't been updated recently)."""
        return time.time() - self._last_updated > max_age_seconds

class JobMonitor:
    """
    Monitor for tracking preparation jobs.
    
    This class scans for job directories and tracks their status
    by reading status.json files.
    """
    
    def __init__(self, working_directory: Optional[Union[str, Path]] = None):
        """
        Initialize the job monitor.
        
        Args:
            working_directory: Directory to monitor for jobs (Path or str)
        """
        if working_directory is None:
            self.working_directory = Path.cwd()
        else:
            self.working_directory = Path(working_directory)
        self.jobs: Dict[str, JobInfo] = {}
        self._last_scan = 0
        self._scan_interval = 3.0  # seconds - much slower to reduce I/O
    
    def set_working_directory(self, directory: Path):
        """Set the working directory to monitor."""
        self.working_directory = Path(directory)
        logger.info(f"Job monitor set to directory: {self.working_directory}")
        self.scan_for_jobs()
    
    def scan_for_jobs(self, force: bool = False):
        """
        Scan working directory for preparation jobs.
        
        Args:
            force: Force scan even if within scan interval
        """
        current_time = time.time()
        if not force and (current_time - self._last_scan) < self._scan_interval:
            return
        
        if not self.working_directory.exists():
            logger.warning(f"Working directory does not exist: {self.working_directory}")
            return
        
        logger.debug(f"Scanning for jobs in: {self.working_directory}")
        
        # Find all job directories (look for directories with status.json file)
        # This is more robust than hardcoding "membrane_*" pattern
        job_dirs = []
        for item in self.working_directory.iterdir():
            if item.is_dir() and (item / "status.json").exists():
                job_dirs.append(item)
        current_job_ids = set()
        
        for job_dir in job_dirs:
            status_file = job_dir / "status.json"
            job_id = str(job_dir)
            current_job_ids.add(job_id)
            
            try:
                # Check if file is empty or being written
                file_size = status_file.stat().st_size
                if file_size == 0:
                    logger.debug(f"Skipping empty status file (being written): {status_file}")
                    continue
                
                # Read status file with explicit encoding and fresh read
                with open(status_file, 'r', encoding='utf-8') as f:
                    content = f.read().strip()
                    if not content:
                        logger.debug(f"Skipping empty status file: {status_file}")
                        continue
                    status_data = json.loads(content)
                
                # Check if this is a new job or existing job
                if job_id in self.jobs:
                    # Update existing job - always refresh the data
                    old_status = self.jobs[job_id].status
                    old_progress = self.jobs[job_id].progress
                    self.jobs[job_id].update_status(status_data)
                    
                    # Log progress changes for debugging
                    if self.jobs[job_id].progress != old_progress:
                        logger.debug(f"Job {job_id}: progress {old_progress:.2f} -> {self.jobs[job_id].progress:.2f}")
                    
                    # Check for process status if job claims to be running
                    if (self.jobs[job_id].status == JobStatus.RUNNING and 
                        old_status == JobStatus.RUNNING):
                        self._verify_job_process(job_dir, self.jobs[job_id])
                else:
                    # Create new job
                    self.jobs[job_id] = JobInfo(job_id, job_dir, status_data)
                    logger.info(f"Found new job: {job_id}")
                
            except json.JSONDecodeError as e:
                # This is expected when file is being written - use debug level
                logger.debug(f"Temporary JSON parse error for {status_file} (file may be updating): {e}")
                continue
            except (IOError, OSError) as e:
                # File access errors might be more serious
                logger.debug(f"Error reading status file {status_file}: {e}")
                continue
        
        # Remove jobs that no longer exist
        removed_jobs = set(self.jobs.keys()) - current_job_ids
        for job_id in removed_jobs:
            logger.info(f"Job removed: {job_id}")
            del self.jobs[job_id]
        
        self._last_scan = current_time
        logger.debug(f"Scan completed. Found {len(self.jobs)} jobs.")
    
    def _verify_job_process(self, job_dir: Path, job_info: JobInfo):
        """
        Verify if a job process is actually still running.
        
        Args:
            job_dir: Job directory path
            job_info: Job information object
        """
        try:
            pid_file = job_dir / "process.pid"
            if pid_file.exists():
                with open(pid_file, 'r') as f:
                    pid = int(f.read().strip())
                
                # Check if process is still running (works on Unix-like systems)
                try:
                    import os
                    import signal
                    os.kill(pid, 0)  # This will raise OSError if process doesn't exist
                    logger.debug(f"Job process {pid} is still running")
                except (OSError, ProcessLookupError):
                    # Process is not running, but status says running
                    # This might indicate a stuck job or incomplete status update
                    logger.warning(f"Job {job_info.job_id} status is 'running' but process {pid} not found")
                    
                    # Check if job should be marked as completed based on output files
                    self._check_completion_by_output(job_dir, job_info)
                    
                except ImportError:
                    # On Windows or if os/signal not available, skip process check
                    logger.debug("Process verification not available on this platform")
        except (IOError, ValueError) as e:
            logger.debug(f"Could not verify process for job {job_info.job_id}: {e}")
    
    def _check_completion_by_output(self, job_dir: Path, job_info: JobInfo):
        """
        Check if job should be marked complete based on output files.
        
        Args:
            job_dir: Job directory path
            job_info: Job information object
        """
        # Look for typical output files that indicate completion
        output_patterns = [
            "*.prmtop", "*.inpcrd", "*.top", "*.crd", 
            "system.prmtop", "system.inpcrd"
        ]
        
        found_outputs = []
        for pattern in output_patterns:
            matches = list(job_dir.glob(pattern))
            found_outputs.extend(matches)
        
        if found_outputs:
            logger.info(f"Job {job_info.job_id} appears complete - found output files: "
                       f"{[f.name for f in found_outputs]}")
            
            # Update status to completed
            try:
                status_file = job_dir / "status.json"
                if status_file.exists():
                    with open(status_file, 'r') as f:
                        status_data = json.load(f)
                    
                    status_data['status'] = 'completed'
                    status_data['end_time'] = datetime.now().isoformat()
                    
                    with open(status_file, 'w') as f:
                        json.dump(status_data, f, indent=2)
                    
                    # Update job info
                    job_info.update_status(status_data)
                    logger.info(f"Automatically marked job {job_info.job_id} as completed")
                    
            except Exception as e:
                logger.error(f"Error updating completion status for {job_info.job_id}: {e}")
    
    def get_active_jobs(self) -> Dict[str, JobInfo]:
        """Get all currently active (running) jobs."""
        return {
            job_id: job for job_id, job in self.jobs.items()
            if job.status == JobStatus.RUNNING
        }
    
    def get_completed_jobs(self) -> Dict[str, JobInfo]:
        """Get all completed jobs."""
        return {
            job_id: job for job_id, job in self.jobs.items()
            if job.status in (JobStatus.COMPLETED, JobStatus.ERROR)
        }
    
    def get_job(self, job_id: str) -> Optional[JobInfo]:
        """Get specific job by ID."""
        return self.jobs.get(job_id)
    
    def refresh_job(self, job_id: str) -> bool:
        """
        Force refresh of a specific job's status.

        Args:
            job_id: Job ID to refresh

        Returns:
            True if job was refreshed successfully
        """
        job = self.jobs.get(job_id)
        if not job:
            return False

        status_file = job.job_dir / "status.json"
        if not status_file.exists():
            logger.debug(f"Status file not found for job {job_id}: {status_file}")
            return False

        # Check if status file is empty (being written) or very small (incomplete)
        try:
            file_size = status_file.stat().st_size
            if file_size == 0:
                logger.debug(f"Status file is empty (being written) for job {job_id}")
                return False
        except OSError:
            logger.debug(f"Could not check status file size for job {job_id}")
            return False

        try:
            # Store old values for comparison
            old_progress = job.progress
            old_steps = job.steps_completed.copy()
            old_current_step = job.current_step
            old_status = job.status

            # Force read fresh data from file (no caching)
            with open(status_file, 'r', encoding='utf-8') as f:
                content = f.read().strip()
                if not content:
                    logger.debug(f"Status file is empty for job {job_id}")
                    return False
                status_data = json.loads(content)

            # Force update the job's status
            job.update_status(status_data)

            # Log if there were actual changes
            if (job.progress != old_progress or 
                job.steps_completed != old_steps or 
                job.current_step != old_current_step or
                job.status != old_status):
                logger.info(f"Job {job_id} updated: progress={job.progress:.2f}, "
                           f"steps={job.steps_completed}, current_step={job.current_step}, "
                           f"status={job.status.value}")
                return True
                
            # Removed stuck job detection - long-running jobs are normal for membrane preparation

            return True

        except json.JSONDecodeError as e:
            # This is expected when file is being written - use debug level
            logger.debug(f"Temporary JSON parse error for job {job_id} (file may be updating): {e}")
            return False
        except IOError as e:
            # File access errors are more serious
            logger.warning(f"Error accessing status file for job {job_id}: {e}")
            return False
    
    def remove_job(self, job_id: str) -> bool:
        """
        Remove a job from monitoring (does not delete files).
        
        Args:
            job_id: Job ID to remove
            
        Returns:
            True if job was removed
        """
        if job_id in self.jobs:
            del self.jobs[job_id]
            logger.info(f"Removed job from monitoring: {job_id}")
            return True
        return False
    
    def cleanup_stale_jobs(self, max_age_seconds: float = 3600):
        """
        Remove stale jobs from monitoring.
        
        Args:
            max_age_seconds: Maximum age before job is considered stale
        """
        stale_jobs = [
            job_id for job_id, job in self.jobs.items()
            if job.is_stale(max_age_seconds)
        ]
        
        for job_id in stale_jobs:
            self.remove_job(job_id)
        
        if stale_jobs:
            logger.info(f"Cleaned up {len(stale_jobs)} stale jobs")
    
    def get_job_statistics(self) -> Dict[str, int]:
        """Get statistics about monitored jobs."""
        stats = {
            "total": len(self.jobs),
            "running": 0,
            "completed": 0,
            "error": 0,
            "unknown": 0
        }
        
        for job in self.jobs.values():
            if job.status == JobStatus.RUNNING:
                stats["running"] += 1
            elif job.status == JobStatus.COMPLETED:
                stats["completed"] += 1
            elif job.status == JobStatus.ERROR:
                stats["error"] += 1
            else:
                stats["unknown"] += 1
        
        return stats
    
    def export_job_summary(self, output_file: Path):
        """
        Export summary of all jobs to a JSON file.
        
        Args:
            output_file: Path to output JSON file
        """
        summary = {
            "export_time": datetime.now().isoformat(),
            "working_directory": str(self.working_directory),
            "statistics": self.get_job_statistics(),
            "jobs": []
        }
        
        for job_id, job in self.jobs.items():
            job_summary = {
                "job_id": job_id,
                "job_dir": str(job.job_dir),
                "status": job.status.value,
                "pdb_file": job.pdb_file,
                "start_time": job.start_time.isoformat() if job.start_time else None,
                "end_time": job.end_time.isoformat() if job.end_time else None,
                "elapsed_time": job.elapsed_time,
                "progress": job.progress,
                "current_step": job.current_step,
                "steps_completed": job.steps_completed,
                "command": job.command,
                "error_message": job.error_message
            }
            summary["jobs"].append(job_summary)
        
        with open(output_file, 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Job summary exported to: {output_file}")