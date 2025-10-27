# gatewizard/gui/widgets/progress_tracker.py
# SPDX-License-Identifier: MIT
# Copyright (c) 2025 Constanza González and Mauricio Bedoya

"""
Progress tracker widget for monitoring system preparation jobs.

This module provides a widget for displaying and monitoring the progress
of system preparation jobs in real-time.
"""

from typing import Optional, Dict, Any
from pathlib import Path
import threading
import time

try:
    import customtkinter as ctk
except ImportError:
    raise ImportError("CustomTkinter is required for GUI")

from gatewizard.gui.constants import (
    COLOR_SCHEME, FONTS, WIDGET_SIZES, LAYOUT, JOB_STATUS, PREPARATION_STEPS
)
from gatewizard.core.job_monitor import JobMonitor, JobInfo, JobStatus
from gatewizard.utils.logger import get_logger
from gatewizard.utils.helpers import format_time

logger = get_logger(__name__)

class ProgressTracker(ctk.CTkFrame):
    """
    Widget for tracking progress of system preparation jobs.
    
    This widget displays running and completed jobs with their current
    status, progress, and other relevant information.
    """
    
    def __init__(
        self,
        parent,
        working_directory: Optional[str] = None,
        refresh_interval: int = 1000  # milliseconds - faster refresh for better responsiveness
    ):
        """
        Initialize the progress tracker.
        
        Args:
            parent: Parent widget
            working_directory: Directory to monitor for jobs
            refresh_interval: How often to refresh job status (ms)
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.working_directory = working_directory
        self.refresh_interval = refresh_interval
        self.job_monitor = JobMonitor(Path(working_directory) if working_directory else None)
        self.job_widgets = {}
        self.refresh_after_id = None
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        
        # Start refresh cycle
        self._start_refresh_cycle()
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            if hasattr(self, 'refresh_after_id') and self.refresh_after_id:
                self.after_cancel(self.refresh_after_id)
                self.refresh_after_id = None
            
            if hasattr(self, 'auto_update_id') and self.auto_update_id:
                self.after_cancel(self.auto_update_id)
                self.auto_update_id = None
                
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")
    
    def _create_widgets(self):
        """Create progress tracker widgets."""
        # Header
        self.header_frame = ctk.CTkFrame(self, fg_color="transparent")
        
        self.title_label = ctk.CTkLabel(
            self.header_frame,
            text="Job Progress",
            font=FONTS['heading'],
            text_color=COLOR_SCHEME['text']
        )
        
        self.refresh_button = ctk.CTkButton(
            self.header_frame,
            text="Refresh",
            width=80,
            height=25,
            command=self.refresh_jobs
        )
        
        self.clear_button = ctk.CTkButton(
            self.header_frame,
            text="Clear Completed",
            width=120,
            height=25,
            command=self._clear_completed_jobs
        )
        
        # Jobs display area
        self.jobs_scroll = ctk.CTkScrollableFrame(
            self,
            height=WIDGET_SIZES['progress_height'],
            fg_color=COLOR_SCHEME['background']
        )
        
        # No jobs message
        self.no_jobs_label = ctk.CTkLabel(
            self.jobs_scroll,
            text="No active jobs",
            font=FONTS['body'],
            text_color=COLOR_SCHEME['inactive']
        )
        self.no_jobs_label.pack(pady=20)
    
    def _setup_layout(self):
        """Setup the layout of widgets."""
        # Header
        self.header_frame.pack(fill="x", padx=LAYOUT['padding_medium'], pady=(LAYOUT['padding_medium'], LAYOUT['padding_small']))
        
        self.title_label.pack(side="left")
        self.clear_button.pack(side="right", padx=(LAYOUT['padding_small'], 0))
        self.refresh_button.pack(side="right", padx=LAYOUT['padding_small'])
        
        # Jobs area
        self.jobs_scroll.pack(fill="both", expand=True, padx=LAYOUT['padding_medium'], pady=LAYOUT['padding_small'])
    
    def _start_refresh_cycle(self):
        """Start the automatic refresh cycle."""
        self._refresh_cycle()
    
    def _refresh_cycle(self):
        """Refresh cycle that runs periodically."""
        try:
            # Force refresh all jobs first to get latest status
            if self.job_monitor:
                self.job_monitor.scan_for_jobs(force=True)
                
                # Force refresh ALL running jobs for real-time updates
                active_jobs = self.job_monitor.get_active_jobs()
                for job_id, job_info in active_jobs.items():
                    if job_info.status == JobStatus.RUNNING:
                        self.job_monitor.refresh_job(job_id)

            # Update display after refreshing status
            self._update_jobs_display()

        except Exception as e:
            logger.error(f"Error in progress tracker refresh cycle: {e}")

        # Dynamic refresh rate based on active jobs
        refresh_interval = 2000  # 2 seconds default
        if self.job_monitor:
            active_jobs = self.job_monitor.get_active_jobs()
            if active_jobs:
                refresh_interval = 1000  # 1 second for active jobs for better responsiveness
            else:
                refresh_interval = 3000  # 3 seconds when no active jobs
        
        # Schedule next refresh
        self.refresh_after_id = self.after(refresh_interval, self._refresh_cycle)

    def set_working_directory(self, directory: str):
        """Set the working directory to monitor."""
        self.working_directory = directory
        self.job_monitor.set_working_directory(Path(directory))
        self.refresh_jobs()
    
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the progress tracker widget."""
        try:
            if hasattr(self, 'title_label'):
                self.title_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'refresh_button'):
                self.refresh_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'clear_button'):
                self.clear_button.configure(font=scaled_fonts['body'])
            if hasattr(self, 'no_jobs_label'):
                self.no_jobs_label.configure(font=scaled_fonts['body'])
            # Update job widgets if they have update_fonts
            for job_id, job_widget in getattr(self, 'job_widgets', {}).items():
                if hasattr(job_widget, 'update_fonts'):
                    job_widget.update_fonts(scaled_fonts)
        except Exception as e:
            logger.warning(f"Error updating fonts in ProgressTracker: {e}")
    
    def refresh_jobs(self):
        """Manually refresh job status."""
        if self.job_monitor:
            # Force scan with no cache
            self.job_monitor.scan_for_jobs(force=True)
            
            # Force refresh all running jobs
            active_jobs = self.job_monitor.get_active_jobs()
            for job_id, job_info in active_jobs.items():
                if job_info.status == JobStatus.RUNNING:
                    self.job_monitor.refresh_job(job_id)
            
            # Update display
            self._update_jobs_display()
            
            logger.info(f"Manual refresh completed. Found {len(active_jobs)} active jobs.")
    
    def _update_jobs_display(self):
        """Update the display of jobs."""
        if not self.job_monitor:
            return
        
        # Scan for jobs
        self.job_monitor.scan_for_jobs()
        
        # Get current jobs
        all_jobs = {**self.job_monitor.get_active_jobs(), **self.job_monitor.get_completed_jobs()}
        
        # Remove widgets for jobs that no longer exist
        current_job_ids = set(all_jobs.keys())
        widget_job_ids = set(self.job_widgets.keys())
        
        for job_id in widget_job_ids - current_job_ids:
            self._remove_job_widget(job_id)
        
        # Update or create widgets for current jobs
        for job_id, job_info in all_jobs.items():
            if job_id in self.job_widgets:
                self._update_job_widget(job_id, job_info)
            else:
                self._create_job_widget(job_id, job_info)
        
        # Show/hide no jobs message
        if all_jobs:
            self.no_jobs_label.pack_forget()
        else:
            self.no_jobs_label.pack(pady=20)
    
    def _create_job_widget(self, job_id: str, job_info: JobInfo):
        """Create a widget for displaying job information."""
        # Job frame
        job_frame = ctk.CTkFrame(self.jobs_scroll, fg_color=COLOR_SCHEME['canvas'])
        job_frame.pack(fill="x", pady=2, padx=5)
        
        # Job header
        header_frame = ctk.CTkFrame(job_frame, fg_color="transparent")
        header_frame.pack(fill="x", padx=10, pady=(10, 5))
        
        # Job name and status
        job_name = Path(job_info.job_dir).name
        job_label = ctk.CTkLabel(
            header_frame,
            text=job_name,
            font=FONTS['body'],
            anchor="w"
        )
        job_label.pack(side="left")
        
        status_label = ctk.CTkLabel(
            header_frame,
            text=job_info.status.value.upper(),
            font=FONTS['small'],
            text_color=self._get_status_color(job_info.status)
        )
        status_label.pack(side="right")
        
        # Progress bar
        progress_frame = ctk.CTkFrame(job_frame, fg_color="transparent")
        progress_frame.pack(fill="x", padx=10, pady=5)
        
        progress_bar = ctk.CTkProgressBar(
            progress_frame,
            height=15,
            progress_color=self._get_status_color(job_info.status)
        )
        progress_bar.pack(fill="x")
        progress_bar.set(job_info.progress)
        
        # Progress text
        progress_text = f"{job_info.progress * 100:.0f}% - {self._get_current_step_text(job_info)}"
        progress_label = ctk.CTkLabel(
            progress_frame,
            text=progress_text,
            font=FONTS['small'],
            text_color=COLOR_SCHEME['text']
        )
        progress_label.pack(pady=(2, 0))
        
        # Time information
        time_frame = ctk.CTkFrame(job_frame, fg_color="transparent")
        time_frame.pack(fill="x", padx=10, pady=(5, 10))
        
        elapsed_time = job_info.elapsed_time or 0
        time_text = f"Elapsed: {format_time(elapsed_time)}"
        
        if job_info.start_time:
            time_text += f" | Started: {job_info.start_time.strftime('%H:%M:%S')}"
        
        time_label = ctk.CTkLabel(
            time_frame,
            text=time_text,
            font=FONTS['small'],
            text_color=COLOR_SCHEME['inactive']
        )
        time_label.pack(side="left")
        
        # Action buttons frame
        action_frame = ctk.CTkFrame(time_frame, fg_color="transparent")
        action_frame.pack(side="right")
        
        # Details button
        details_button = ctk.CTkButton(
            action_frame,
            text="Details",
            width=60,
            height=20,
            command=lambda: self._show_job_details(job_info)
        )
        details_button.pack(side="left", padx=(0, 5))
        
        # Conditional action buttons
        if job_info.status == JobStatus.COMPLETED:
            open_button = ctk.CTkButton(
                action_frame,
                text="Open Folder",
                width=80,
                height=20,
                command=lambda: self._open_job_folder(str(job_info.job_dir))
            )
            open_button.pack(side="left")
        elif job_info.status == JobStatus.ERROR:
            retry_button = ctk.CTkButton(
                action_frame,
                text="Retry",
                width=60,
                height=20,
                command=lambda: self._retry_job(job_info)
            )
            retry_button.pack(side="left")
        
        # Store widget references
        self.job_widgets[job_id] = {
            'frame': job_frame,
            'job_label': job_label,
            'status_label': status_label,
            'progress_bar': progress_bar,
            'progress_label': progress_label,
            'time_label': time_label,
            'details_button': details_button
        }
    
    def _update_job_widget(self, job_id: str, job_info: JobInfo):
        """Update an existing job widget."""
        if job_id not in self.job_widgets:
            return
        
        widgets = self.job_widgets[job_id]
        
        # Update status
        status_color = self._get_status_color(job_info.status)
        widgets['status_label'].configure(
            text=job_info.status.value.upper(),
            text_color=status_color
        )
        
        # Update progress
        progress_value = job_info.progress
        widgets['progress_bar'].set(progress_value)
        widgets['progress_bar'].configure(
            progress_color=status_color
        )
        
        # Update progress text with more detailed information
        progress_text = f"{progress_value * 100:.0f}% - {self._get_current_step_text(job_info)}"
        widgets['progress_label'].configure(text=progress_text)
        
        # Update time
        elapsed_time = job_info.elapsed_time or 0
        time_text = f"Elapsed: {format_time(elapsed_time)}"
        
        if job_info.start_time:
            time_text += f" | Started: {job_info.start_time.strftime('%H:%M:%S')}"
        
        widgets['time_label'].configure(text=time_text)
        
        # Force update the widget display
        try:
            widgets['frame'].update_idletasks()
        except Exception:
            pass  # Ignore errors during update
    
    def _remove_job_widget(self, job_id: str):
        """Remove a job widget."""
        if job_id in self.job_widgets:
            self.job_widgets[job_id]['frame'].destroy()
            del self.job_widgets[job_id]
    
    def _get_status_color(self, status: JobStatus) -> str:
        """Get color for job status."""
        color_map = {
            JobStatus.RUNNING: COLOR_SCHEME['highlight'],
            JobStatus.COMPLETED: COLOR_SCHEME['active'],
            JobStatus.ERROR: "#dc3545",  # Red
            JobStatus.UNKNOWN: COLOR_SCHEME['inactive']
        }
        return color_map.get(status, COLOR_SCHEME['inactive'])
    
    def _get_current_step_text(self, job_info: JobInfo) -> str:
        """Get text description of current step."""
        if job_info.status == JobStatus.COMPLETED:
            return "Completed"
        elif job_info.status == JobStatus.ERROR:
            error_msg = job_info.error_message or 'Unknown error'
            # Truncate long error messages
            if len(error_msg) > 50:
                error_msg = error_msg[:47] + "..."
            return f"Error: {error_msg}"
        elif job_info.steps_completed:
            last_step = job_info.steps_completed[-1]
            # Provide more descriptive text for each step
            step_descriptions = {
                "Starting": "Initializing preparation",
                "Initializing": "Initializing preparation", 
                "MEMEMBED": "Running MEMEMBED orientation",
                "Running Packmol": "Packing lipids with Packmol",
                "Packmol": "Packing lipids with Packmol", 
                "Running": "Processing system",
                "Processing": "Processing system",
                "pdb4amber": "Preparing PDB with pdb4amber",
                "Parametrizing": "Parametrizing system",
                "Parameterization": "Parametrizing system",
                "tleap": "Running tleap parametrization",
                "Building Amber input": "Converting to Amber format",
                "Transforming to AMBER": "Converting to Amber format",
                "Building Input": "Creating output files",
                "Finalizing": "Finalizing output",
                "Completed": "Completed"
            }
            return step_descriptions.get(last_step, f"Running: {last_step}")
        else:
            # Check current_step number if no completed steps
            current_step = job_info.current_step
            
            # Get workflow configuration to determine expected steps
            config = getattr(job_info, 'status_data', {}).get('config', {})
            preoriented = config.get('preoriented', True)
            parametrize = config.get('parametrize', True)
            
            # Create step names based on actual workflow
            step_names = ["Starting"]
            if not preoriented:
                step_names.append("MEMEMBED")
            step_names.append("Packmol")
            if parametrize:
                step_names.extend(["pdb4amber", "tleap"])
            step_names.append("Completed")
            
            if 0 <= current_step < len(step_names):
                step_name = step_names[current_step]
                step_descriptions = {
                    "Starting": "Initializing preparation",
                    "MEMEMBED": "Running MEMEMBED orientation",
                    "Packmol": "Packing lipids with Packmol",
                    "pdb4amber": "Preparing PDB with pdb4amber",
                    "tleap": "Running tleap parametrization",
                    "Completed": "Completed"
                }
                return step_descriptions.get(step_name, step_name)
            
            # Check if job has been running for a while to provide better feedback
            if job_info.start_time and job_info.elapsed_time:
                elapsed = job_info.elapsed_time
                if elapsed > 60:  # More than 1 minute
                    return "Processing... (please wait)"
                elif elapsed > 30:  # More than 30 seconds
                    return "Starting preparation (initializing)..."
            
            return "Starting preparation..."
    
    def _show_job_details(self, job_info: JobInfo):
        """Show detailed job information."""
        try:
            # Ensure parent window is properly visible before creating dialog
            self.update_idletasks()
            
            dialog = JobDetailDialog(self, job_info)
            # Ensure dialog is brought to front
            dialog.lift()
            dialog.focus_force()
        except Exception as e:
            logger.error(f"Error showing job details: {e}")
            # Try to show a simple message dialog instead
            try:
                import tkinter.messagebox as messagebox
                messagebox.showerror(
                    "Error", 
                    f"Could not open job details dialog:\n{str(e)}"
                )
            except Exception:
                # If even that fails, just log it
                logger.error("Could not show error dialog either")
    
    def _open_job_folder(self, job_dir: str):
        """Open job folder in file explorer."""
        try:
            import os
            import subprocess
            import sys
            
            if sys.platform == "win32":
                os.startfile(job_dir)
            elif sys.platform == "darwin":  # macOS
                subprocess.run(["open", job_dir])
            else:  # Linux
                subprocess.run(["xdg-open", job_dir])
                
        except Exception as e:
            logger.error(f"Error opening job folder: {e}")
    
    def _retry_job(self, job_info: JobInfo):
        """Retry a failed job."""
        # This would need to be implemented based on the specific retry mechanism
        logger.info(f"Retry requested for job: {job_info.job_dir}")
        # Placeholder for retry functionality
    
    def _clear_completed_jobs(self):
        """Clear completed jobs from the display."""
        if not self.job_monitor:
            return
        
        completed_jobs = self.job_monitor.get_completed_jobs()
        
        if not completed_jobs:
            return
        
        # Remove completed job widgets
        for job_id in completed_jobs.keys():
            self._remove_job_widget(job_id)
            self.job_monitor.remove_job(job_id)
        
        logger.info(f"Cleared {len(completed_jobs)} completed jobs")
    
    def cleanup(self):
        """Cleanup resources when widget is destroyed."""
        # Cancel refresh cycle
        if self.refresh_after_id:
            self.after_cancel(self.refresh_after_id)
            self.refresh_after_id = None
        
        # Cleanup job monitor
        if self.job_monitor:
            self.job_monitor.cleanup_stale_jobs()

class JobDetailDialog(ctk.CTkToplevel):
    """
    Dialog for displaying detailed job information.
    
    Shows information about a specific job including
    full command, output logs, and detailed progress.
    """
    
    def __init__(self, parent, job_info: JobInfo):
        """
        Initialize the job detail dialog.
        
        Args:
            parent: Parent widget
            job_info: JobInfo object with job details
        """
        super().__init__(parent)
        
        self.job_info = job_info
        self.auto_update_id = None
        
        # Configure window - made bigger in y direction
        self.title(f"Job Details - {Path(job_info.job_dir).name}")
        self.geometry("700x900")
        self.resizable(True, True)
        self.transient(parent)
        
        # Set protocol for window close button
        self.protocol("WM_DELETE_WINDOW", self._close_dialog)
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        self._populate_data()
        
        # Center dialog and make it visible
        self._center_window()
        
        # Update the window to ensure it's fully rendered
        self.update_idletasks()
        
        # Set grab after window is fully initialized and visible
        self.after(10, self._set_grab)
        
        # Start auto-update for logs if job is running
        if self.job_info.status == JobStatus.RUNNING:
            self._start_auto_update()
    
    def cleanup_callbacks(self):
        """Cancel all scheduled callbacks to prevent errors during shutdown."""
        try:
            if hasattr(self, 'auto_update_id') and self.auto_update_id:
                self.after_cancel(self.auto_update_id)
                self.auto_update_id = None
                
            logger.debug(f"Cleaned up callbacks for {type(self).__name__}")
        except Exception as e:
            logger.debug(f"Error cleaning up callbacks in {type(self).__name__}: {e}")
    
    def _create_widgets(self):
        """Create dialog widgets."""
        # Main frame
        self.main_frame = ctk.CTkFrame(self)
        
        # Job information frame
        self.info_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.info_label = ctk.CTkLabel(
            self.info_frame,
            text="Job Information",
            font=FONTS['heading']
        )
        
        # Details text area
        self.details_text = ctk.CTkTextbox(
            self.info_frame,
            height=200,
            font=FONTS['small']
        )
        
        # Command frame
        self.command_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.command_label = ctk.CTkLabel(
            self.command_frame,
            text="Command",
            font=FONTS['heading']
        )
        
        self.command_text = ctk.CTkTextbox(
            self.command_frame,
            height=80,
            font=FONTS['code']
        )
        
        # Logs frame
        self.logs_frame = ctk.CTkFrame(self.main_frame, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.logs_label = ctk.CTkLabel(
            self.logs_frame,
            text="Recent Log Output",
            font=FONTS['heading']
        )
        
        self.logs_text = ctk.CTkTextbox(
            self.logs_frame,
            height=300,  # Increased height for better visibility
            font=FONTS['code']
        )
        
        # Buttons frame
        self.buttons_frame = ctk.CTkFrame(self.main_frame, fg_color="transparent")
        
        self.close_button = ctk.CTkButton(
            self.buttons_frame,
            text="Close",
            command=self._close_dialog
        )
        
        self.open_folder_button = ctk.CTkButton(
            self.buttons_frame,
            text="Open Folder",
            command=self._open_folder
        )
        
        self.copy_command_button = ctk.CTkButton(
            self.buttons_frame,
            text="Copy Command",
            command=self._copy_command
        )
    
    def _setup_layout(self):
        """Setup dialog layout."""
        self.main_frame.pack(fill="both", expand=True, padx=20, pady=20)
        
        # Info frame
        self.info_frame.pack(fill="x", pady=(0, 10))
        
        self.info_label.pack(anchor="w", padx=10, pady=(10, 5))
        self.details_text.pack(fill="x", padx=10, pady=(0, 10))
        
        # Command frame
        self.command_frame.pack(fill="x", pady=(0, 10))
        
        self.command_label.pack(anchor="w", padx=10, pady=(10, 5))
        self.command_text.pack(fill="x", padx=10, pady=(0, 10))
        
        # Logs frame
        self.logs_frame.pack(fill="both", expand=True, pady=(0, 10))
        
        self.logs_label.pack(anchor="w", padx=10, pady=(10, 5))
        self.logs_text.pack(fill="both", expand=True, padx=10, pady=(0, 10))
        
        # Buttons
        self.buttons_frame.pack(fill="x")
        
        self.close_button.pack(side="right", padx=(10, 0))
        self.copy_command_button.pack(side="right", padx=10)
        self.open_folder_button.pack(side="right")
    
    def _populate_data(self):
        """Populate dialog with job data."""
        # Job details
        details = f"""Job Directory: {self.job_info.job_dir}
Status: {self.job_info.status.value.upper()}
Start Time: {self.job_info.start_time.strftime('%Y-%m-%d %H:%M:%S') if self.job_info.start_time else 'Unknown'}
End Time: {self.job_info.end_time.strftime('%Y-%m-%d %H:%M:%S') if self.job_info.end_time else 'Running'}
Elapsed Time: {format_time(self.job_info.elapsed_time) if self.job_info.elapsed_time else 'Unknown'}
Progress: {self.job_info.progress * 100:.1f}%
Current Step: {self.job_info.current_step}

Completed Steps:
"""
        
        for step in self.job_info.steps_completed:
            details += f"  ✓ {step}\n"
        
        if self.job_info.error_message:
            details += f"\nError: {self.job_info.error_message}"
        
        self.details_text.insert("1.0", details)
        self.details_text.configure(state="disabled")
        
        # Command
        self.command_text.insert("1.0", self.job_info.command)
        self.command_text.configure(state="disabled")
        
        # Logs
        self._load_recent_logs()
    
    def _load_recent_logs(self):
        """Load recent log output."""
        try:
            log_file = Path(self.job_info.job_dir) / "logs" / "preparation.log"
            
            if log_file.exists():
                # Read last 50 lines of log file
                with open(log_file, 'r') as f:
                    lines = f.readlines()
                
                recent_lines = lines[-50:] if len(lines) > 50 else lines
                log_content = ''.join(recent_lines)
                
                self.logs_text.insert("1.0", log_content)
            else:
                self.logs_text.insert("1.0", "No log file found")
                
        except Exception as e:
            self.logs_text.insert("1.0", f"Error reading logs: {e}")
        
        self.logs_text.configure(state="disabled")
    
    def _start_auto_update(self):
        """Start auto-updating logs for running jobs."""
        self._auto_update_logs()
    
    def _auto_update_logs(self):
        """Auto-update log output every 2 seconds."""
        try:
            # Only update if job is still running
            if self.job_info.status == JobStatus.RUNNING:
                # Clear and reload logs
                self.logs_text.configure(state="normal")
                self.logs_text.delete("1.0", "end")
                
                log_file = Path(self.job_info.job_dir) / "logs" / "preparation.log"
                
                if log_file.exists():
                    # Read last 100 lines of log file for better context
                    with open(log_file, 'r') as f:
                        lines = f.readlines()
                    
                    recent_lines = lines[-100:] if len(lines) > 100 else lines
                    log_content = ''.join(recent_lines)
                    
                    self.logs_text.insert("1.0", log_content)
                    # Auto-scroll to bottom
                    self.logs_text.see("end")
                else:
                    self.logs_text.insert("1.0", "No log file found")
                
                self.logs_text.configure(state="disabled")
                
                # Schedule next update
                self.auto_update_id = self.after(2000, self._auto_update_logs)
        except Exception as e:
            logger.error(f"Error auto-updating logs: {e}")
            # Schedule next update anyway
            self.auto_update_id = self.after(2000, self._auto_update_logs)
    
    def _center_window(self):
        """Center the dialog window."""
        self.update_idletasks()
        
        # Get screen dimensions
        screen_width = self.winfo_screenwidth()
        screen_height = self.winfo_screenheight()
        
        # Get window dimensions
        window_width = self.winfo_width()
        window_height = self.winfo_height()
        
        # Calculate position
        x = (screen_width - window_width) // 2
        y = (screen_height - window_height) // 2
        
        self.geometry(f"+{x}+{y}")
    
    def _set_grab(self):
        """Set grab after window is visible to avoid 'grab failed' error."""
        try:
            self.grab_set()
        except Exception as e:
            logger.warning(f"Could not set grab on dialog: {e}")
    
    def _close_dialog(self):
        """Properly close the dialog."""
        # Stop auto-update if running
        if self.auto_update_id:
            self.after_cancel(self.auto_update_id)
            self.auto_update_id = None
        
        try:
            self.grab_release()
        except Exception:
            pass  # Ignore errors when releasing grab
        self.destroy()
    
    def _open_folder(self):
        """Open the job folder."""
        try:
            import os
            import subprocess
            import sys
            
            if sys.platform == "win32":
                os.startfile(self.job_info.job_dir)
            elif sys.platform == "darwin":  # macOS
                subprocess.run(["open", self.job_info.job_dir])
            else:  # Linux
                subprocess.run(["xdg-open", self.job_info.job_dir])
                
        except Exception as e:
            logger.error(f"Error opening job folder: {e}")
    
    def _copy_command(self):
        """Copy the command to clipboard."""
        try:
            self.clipboard_clear()
            self.clipboard_append(self.job_info.command)
            self.update()  # Required for clipboard to work
            
            # Show brief feedback
            original_text = self.copy_command_button.cget("text")
            self.copy_command_button.configure(text="Copied!")
            self.after(1000, lambda: self.copy_command_button.configure(text=original_text))
            
        except Exception as e:
            logger.error(f"Error copying command: {e}")

class SimpleProgressBar(ctk.CTkFrame):
    """
    Simple progress bar widget for individual operations.
    
    A lightweight progress indicator for showing progress of
    individual operations within the application.
    """
    
    def __init__(self, parent, width: int = 200, height: int = 20):
        """
        Initialize the progress bar.
        
        Args:
            parent: Parent widget
            width: Progress bar width
            height: Progress bar height
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['inactive'], width=width, height=height)
        
        self.width = width
        self.height = height
        self.progress = 0.0
        
        # Configure frame
        self.pack_propagate(False)
        
        # Progress fill
        self.progress_fill = ctk.CTkFrame(
            self,
            fg_color=COLOR_SCHEME['highlight'],
            width=0,
            height=height - 4
        )
        self.progress_fill.place(x=2, y=2)
        
        # Progress text
        self.progress_text = ctk.CTkLabel(
            self,
            text="0%",
            font=FONTS['small'],
            text_color="white"
        )
        self.progress_text.place(relx=0.5, rely=0.5, anchor="center")
    
    def set_progress(self, progress: float):
        """
        Set the progress value.
        
        Args:
            progress: Progress value between 0.0 and 1.0
        """
        self.progress = max(0.0, min(1.0, progress))
        
        # Update fill width
        fill_width = int((self.width - 4) * self.progress)
        self.progress_fill.configure(width=max(0, fill_width))
        
        # Update text
        self.progress_text.configure(text=f"{self.progress * 100:.0f}%")
    
    def set_color(self, color: str):
        """
        Set the progress bar color.
        
        Args:
            color: Color string (hex or name)
        """
        self.progress_fill.configure(fg_color=color)
    
    def set_text(self, text: str):
        """
        Set custom text for the progress bar.
        
        Args:
            text: Custom text to display
        """
        self.progress_text.configure(text=text)
    
    def pulse(self):
        """Start a pulsing animation for indeterminate progress."""
        # This could be implemented for indeterminate progress indication
        pass

class JobSummaryWidget(ctk.CTkFrame):
    def update_fonts(self, scaled_fonts):
        """Update all fonts in the job summary widget."""
        try:
            if hasattr(self, 'title_label'):
                self.title_label.configure(font=scaled_fonts['heading'])
            if hasattr(self, 'running_count'):
                self.running_count.configure(font=scaled_fonts['title'])
            if hasattr(self, 'running_label'):
                self.running_label.configure(font=scaled_fonts['small'])
            if hasattr(self, 'completed_count'):
                self.completed_count.configure(font=scaled_fonts['title'])
            if hasattr(self, 'completed_label'):
                self.completed_label.configure(font=scaled_fonts['small'])
            if hasattr(self, 'failed_count'):
                self.failed_count.configure(font=scaled_fonts['title'])
            if hasattr(self, 'failed_label'):
                self.failed_label.configure(font=scaled_fonts['small'])
        except Exception as e:
            logger.warning(f"Error updating fonts in JobSummaryWidget: {e}")
    """
    Widget for displaying a summary of all jobs.
    
    Shows statistics about running, completed, and failed jobs
    in a compact format.
    """
    
    def __init__(self, parent, job_monitor: JobMonitor):
        """
        Initialize the job summary widget.
        
        Args:
            parent: Parent widget
            job_monitor: JobMonitor instance to get statistics from
        """
        super().__init__(parent, fg_color=COLOR_SCHEME['content_inside_bg'])
        
        self.job_monitor = job_monitor
        
        # Create widgets
        self._create_widgets()
        self._setup_layout()
        self.update_summary()
    
    def _create_widgets(self):
        """Create summary widgets."""
        self.title_label = ctk.CTkLabel(
            self,
            text="Job Summary",
            font=FONTS['heading']
        )
        
        # Statistics frame
        self.stats_frame = ctk.CTkFrame(self, fg_color="transparent")
        
        # Running jobs
        self.running_frame = ctk.CTkFrame(self.stats_frame, fg_color=COLOR_SCHEME['highlight'])
        self.running_count = ctk.CTkLabel(
            self.running_frame,
            text="0",
            font=FONTS['title'],
            text_color="white"
        )
        self.running_label = ctk.CTkLabel(
            self.running_frame,
            text="Running",
            font=FONTS['small'],
            text_color="white"
        )
        
        # Completed jobs
        self.completed_frame = ctk.CTkFrame(self.stats_frame, fg_color=COLOR_SCHEME['active'])
        self.completed_count = ctk.CTkLabel(
            self.completed_frame,
            text="0",
            font=FONTS['title'],
            text_color="white"
        )
        self.completed_label = ctk.CTkLabel(
            self.completed_frame,
            text="Completed",
            font=FONTS['small'],
            text_color="white"
        )
        
        # Failed jobs
        self.failed_frame = ctk.CTkFrame(self.stats_frame, fg_color="#dc3545")
        self.failed_count = ctk.CTkLabel(
            self.failed_frame,
            text="0",
            font=FONTS['title'],
            text_color="white"
        )
        self.failed_label = ctk.CTkLabel(
            self.failed_frame,
            text="Failed",
            font=FONTS['small'],
            text_color="white"
        )
    
    def _setup_layout(self):
        """Setup the layout."""
        self.title_label.pack(pady=(10, 5))
        
        self.stats_frame.pack(fill="x", padx=10, pady=(0, 10))
        
        # Pack stat frames
        for frame in [self.running_frame, self.completed_frame, self.failed_frame]:
            frame.pack(side="left", fill="both", expand=True, padx=5)
        
        # Pack labels within each frame
        for count_label, text_label in [
            (self.running_count, self.running_label),
            (self.completed_count, self.completed_label),
            (self.failed_count, self.failed_label)
        ]:
            count_label.pack(pady=(10, 2))
            text_label.pack(pady=(0, 10))
    
    def update_summary(self):
        """Update the summary statistics."""
        if not self.job_monitor:
            return
        
        stats = self.job_monitor.get_job_statistics()
        
        self.running_count.configure(text=str(stats.get('running', 0)))
        self.completed_count.configure(text=str(stats.get('completed', 0)))
        self.failed_count.configure(text=str(stats.get('error', 0)))