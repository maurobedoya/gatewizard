"""Scheduler adapter protocol."""

from __future__ import annotations

from typing import List, Optional, Protocol

from gatewizard.utils.cluster.types import BatchScriptRequest, RemoteJobHandle


class SchedulerAdapter(Protocol):
    """Common interface for Slurm / PBS / SGE adapters."""

    name: str

    def render_batch_script(self, req: BatchScriptRequest) -> str:
        ...

    def submit_command(self, script_name: str) -> List[str]:
        ...

    def parse_submit_output(self, text: str) -> str:
        ...

    def status_command(self, job_id: Optional[str] = None) -> List[str]:
        ...

    def parse_status(self, text: str, job_id: Optional[str] = None) -> List[RemoteJobHandle]:
        ...

    def cancel_command(self, job_id: str) -> List[str]:
        ...

    def inventory_command(self) -> List[str]:
        ...
