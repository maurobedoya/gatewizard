from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path

# Create monitor for your working directory
monitor = JobMonitor(working_directory=Path("./systems"))

# Scan for jobs
monitor.scan_for_jobs(force=True)

# Get active jobs
active_jobs = monitor.get_active_jobs()

if active_jobs:
    print(f"✓ Found {len(active_jobs)} active job(s)")
    
    for job_id, job_info in active_jobs.items():
        print(f"Job: {job_info.job_dir.name}")
        print(f"  Status: {job_info.status.value}")
        print(f"  Progress: {job_info.progress:.1f}%")
        print(f"  Current step: {job_info.current_step}")
        print(f"  Elapsed: {job_info.elapsed_time:.1f}s")
        
        # Show completed steps
        if job_info.steps_completed:
            print(f"  Completed steps:")
            for step in job_info.steps_completed:
                print(f"    ✓ {step}")
else:
    print("No active jobs found")

# Check for completed jobs
completed_jobs = monitor.get_completed_jobs()
if completed_jobs:
    print(f"✓ Found {len(completed_jobs)} completed job(s)")
    for job_id, job_info in completed_jobs.items():
        print(f"  - {job_info.job_dir.name}: {job_info.status.value}")
