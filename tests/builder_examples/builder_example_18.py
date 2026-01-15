from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path

# Monitor all jobs in a directory
monitor = JobMonitor(working_directory=Path("./systems"))
monitor.scan_for_jobs(force=True)

# Check active jobs
active_jobs = monitor.get_active_jobs()
print(f"\n{'='*60}")
print(f"ACTIVE JOBS: {len(active_jobs)}")
print(f"{'='*60}")

for job_id, job_info in active_jobs.items():
    print(f"\n📁 {job_info.job_dir.name}")
    print(f"   Status: {job_info.status.value}")
    print(f"   Progress: {job_info.progress:.1f}%")
    print(f"   Current: {job_info.current_step}")
    print(f"   Runtime: {job_info.elapsed_time:.0f}s")

    if job_info.steps_completed:
        print(f"   Completed steps:")
        for step in job_info.steps_completed:
            print(f"     ✓ {step}")

# Check completed jobs
completed_jobs = monitor.get_completed_jobs()
print(f"\n{'='*60}")
print(f"COMPLETED JOBS: {len(completed_jobs)}")
print(f"{'='*60}")

for job_id, job_info in completed_jobs.items():
    status_icon = "✓" if job_info.status.value == "completed" else "✗"
    print(f"{status_icon} {job_info.job_dir.name}: {job_info.status.value} ({job_info.elapsed_time:.0f}s)")