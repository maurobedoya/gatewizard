"""Slurm scheduler adapter."""

from __future__ import annotations

from typing import List, Optional

from gatewizard.utils.cluster.resources import parse_sbatch_output, parse_squeue_me
from gatewizard.utils.cluster.templates import render_batch_script
from gatewizard.utils.cluster.types import BatchScriptRequest, RemoteJobHandle


class SlurmAdapter:
    name = "slurm"

    def render_batch_script(self, req: BatchScriptRequest) -> str:
        return render_batch_script(req)

    def submit_command(self, script_name: str) -> List[str]:
        return ["sbatch", script_name]

    def parse_submit_output(self, text: str) -> str:
        return parse_sbatch_output(text)

    def status_command(self, job_id: Optional[str] = None) -> List[str]:
        # Pipe-separated: reason/node/cpus stay aligned when %R is empty.
        fmt = "%i|%T|%j|%P|%R|%N|%M|%C"
        if job_id:
            return ["squeue", "-j", str(job_id), "-h", "-o", fmt]
        return ["squeue", "--me", "-h", "-o", fmt]

    def parse_status(self, text: str, job_id: Optional[str] = None) -> List[RemoteJobHandle]:
        handles: List[RemoteJobHandle] = []
        for row in parse_squeue_me(text):
            if job_id and str(row.get("job_id")) != str(job_id):
                continue
            handles.append(
                RemoteJobHandle(
                    scheduler=self.name,
                    job_id=str(row.get("job_id") or ""),
                    state=str(row.get("state") or ""),
                    name=str(row.get("name") or ""),
                    partition=str(row.get("partition") or ""),
                    reason=str(row.get("reason") or ""),
                    node_list=str(row.get("node_list") or ""),
                    elapsed=str(row.get("elapsed") or ""),
                    cpus=int(row.get("cpus") or 0),
                    raw=str(row.get("raw") or ""),
                )
            )
        return handles

    def cancel_command(self, job_id: str) -> List[str]:
        return ["scancel", str(job_id)]

    def name_status_command(self, job_name: str) -> List[str]:
        fmt = "%i|%T|%j|%P|%R|%N|%M|%C"
        return ["squeue", "--me", "-n", str(job_name), "-h", "-o", fmt]

    def name_accounting_command(self, job_name: str) -> List[str]:
        return [
            "sacct",
            "--name",
            str(job_name),
            "-X",
            "-n",
            "-P",
            "-o",
            "JobID,JobName,State,End",
            "--starttime",
            "2024-01-01",
        ]

    def inventory_command(self) -> List[str]:
        return ["sinfo", "-o", "%P %a %D %c %G %m %l"]

    def accounting_command(self, job_id: str) -> List[str]:
        return [
            "sacct",
            "-j",
            str(job_id),
            "-n",
            "-o",
            "JobID,State,Elapsed,ExitCode",
            "-P",
        ]


def get_scheduler(name: str = "slurm"):
    key = (name or "slurm").lower()
    if key == "slurm":
        return SlurmAdapter()
    if key in {"pbs", "torque"}:
        from gatewizard.utils.cluster.scheduler.pbs import PbsAdapter

        return PbsAdapter()
    # Future: sge
    return SlurmAdapter()
