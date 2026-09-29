from gatewizard.core.job_monitor import JobMonitor
from pathlib import Path
import time

# Start a system preparation (runs in background)
from gatewizard.core.builder import Builder

builder = Builder()
success, message, job_dir = builder.prepare_system(
    pdb_file="protein_protonated_prepared.pdb",
    working_dir="./systems",
    upper_lipids=["POPC"],
    lower_lipids=["POPC"],
    lipid_ratios="1//1",
)

if not success:
    print(f"✗ Preparation not started: {message}")
else:
    print(f"✓ Job started: {job_dir}")

    # Monitor progress; cap wait so a missing packmol-memgen cannot hang pytest.
    monitor = JobMonitor(working_directory=Path("./systems"))
    deadline = time.time() + 30

    while time.time() < deadline:
        monitor.scan_for_jobs(force=True)
        active_jobs = monitor.get_active_jobs()

        if not active_jobs:
            completed = monitor.get_completed_jobs()
            if completed:
                for job_id, job_info in completed.items():
                    if str(job_dir) in job_id:
                        print(f"\n✓ Job completed: {job_info.status.value}")
                        print(f"  Total time: {job_info.elapsed_time:.1f}s")
                        break
            break

        for job_id, job_info in active_jobs.items():
            if str(job_dir) in job_id:
                print(
                    f"\rProgress: {job_info.progress:.1f}% - {job_info.current_step}",
                    end="",
                    flush=True,
                )

        time.sleep(2)
    else:
        print("\nStopped waiting after 30s (job still running)")
