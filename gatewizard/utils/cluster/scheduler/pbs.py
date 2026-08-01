"""PBS/Torque scheduler stub (future).

Implements the same surface as SlurmAdapter so a second backend can be
wired without changing the GUI. Not used in production yet.
"""

from __future__ import annotations

from typing import List, Optional

from gatewizard.utils.cluster.templates import render_batch_script
from gatewizard.utils.cluster.types import BatchScriptRequest, RemoteJobHandle


class PbsAdapter:
    """Minimal PBS adapter stub — submit/status parsing TBD per site."""

    name = "pbs"

    def render_batch_script(self, req: BatchScriptRequest) -> str:
        # Reuse Slurm-oriented template for now; sites should use custom_template.
        return render_batch_script(req)

    def submit_command(self, script_name: str) -> List[str]:
        return ["qsub", script_name]

    def parse_submit_output(self, text: str) -> str:
        # Typical: "12345.hostname"
        token = (text or "").strip().split()[0] if (text or "").strip() else ""
        return token.split(".")[0] if token else ""

    def status_command(self, job_id: Optional[str] = None) -> List[str]:
        if job_id:
            return ["qstat", "-f", str(job_id)]
        return ["qstat", "-u", "$USER"]

    def parse_status(self, text: str, job_id: Optional[str] = None) -> List[RemoteJobHandle]:
        # Intentionally minimal until a real PBS site is wired.
        _ = text, job_id
        return []

    def cancel_command(self, job_id: str) -> List[str]:
        return ["qdel", str(job_id)]

    def inventory_command(self) -> List[str]:
        return ["pbsnodes", "-a"]
