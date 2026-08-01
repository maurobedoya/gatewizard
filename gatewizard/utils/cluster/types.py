"""Shared types for remote cluster / HPC equilibration."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional


WORKDIR_STRATEGIES = (
    "run_in_place",
    "scratch_job_id",
    "scratch_named",
    "tmpdir",
    "custom_template",
)

SCHEDULERS = ("slurm", "pbs", "sge")


@dataclass
class ModulePackage:
    """One Environment Modules / Lmod package entry."""

    name: str
    full_name: str
    category: str = ""
    version: str = ""
    features: List[str] = field(default_factory=list)
    is_default: bool = False
    engine: Optional[str] = None  # namd|gromacs|amber|openmm|cuda|other

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class PartitionInfo:
    name: str
    nodes: int = 0
    cpus_per_node: int = 0
    max_gpus: int = 0
    max_mem: str = ""
    max_time: str = ""
    avail: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class NodeInfo:
    """One compute node from ``sinfo -N``."""

    name: str
    partition: str = ""
    state: str = ""
    cpus: int = 0
    gpus: int = 0
    gres: str = ""
    features: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class ProbeResult:
    """Snapshot of remote software / paths / scheduler inventory."""

    hostname: str = ""
    home: str = ""
    data_dir: str = ""
    scratch_dir: str = ""
    modules: List[ModulePackage] = field(default_factory=list)
    engine_modules: Dict[str, List[ModulePackage]] = field(default_factory=dict)
    partitions: List[PartitionInfo] = field(default_factory=list)
    nodes: List[NodeInfo] = field(default_factory=list)
    raw_module_avail: str = ""
    raw_sinfo: str = ""
    raw_sinfo_nodes: str = ""
    probed_at: str = ""
    errors: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "hostname": self.hostname,
            "home": self.home,
            "data_dir": self.data_dir,
            "scratch_dir": self.scratch_dir,
            "modules": [m.to_dict() for m in self.modules],
            "engine_modules": {
                k: [m.to_dict() for m in v] for k, v in self.engine_modules.items()
            },
            "partitions": [p.to_dict() for p in self.partitions],
            "nodes": [n.to_dict() for n in self.nodes],
            "raw_module_avail": self.raw_module_avail,
            "raw_sinfo": self.raw_sinfo,
            "raw_sinfo_nodes": self.raw_sinfo_nodes,
            "probed_at": self.probed_at,
            "errors": list(self.errors),
        }

    @classmethod
    def from_dict(cls, data: Optional[Dict[str, Any]]) -> "ProbeResult":
        if not data:
            return cls()
        modules = [ModulePackage(**m) for m in data.get("modules") or [] if isinstance(m, dict)]
        engine_modules: Dict[str, List[ModulePackage]] = {}
        for key, vals in (data.get("engine_modules") or {}).items():
            engine_modules[key] = [
                ModulePackage(**m) for m in vals if isinstance(m, dict)
            ]
        partitions = [
            PartitionInfo(**{k: v for k, v in p.items() if k in PartitionInfo.__dataclass_fields__})
            for p in data.get("partitions") or []
            if isinstance(p, dict)
        ]
        nodes = [
            NodeInfo(**{k: v for k, v in n.items() if k in NodeInfo.__dataclass_fields__})
            for n in data.get("nodes") or []
            if isinstance(n, dict)
        ]
        return cls(
            hostname=str(data.get("hostname") or ""),
            home=str(data.get("home") or ""),
            data_dir=str(data.get("data_dir") or ""),
            scratch_dir=str(data.get("scratch_dir") or ""),
            modules=modules,
            engine_modules=engine_modules,
            partitions=partitions,
            nodes=nodes,
            raw_module_avail=str(data.get("raw_module_avail") or ""),
            raw_sinfo=str(data.get("raw_sinfo") or ""),
            raw_sinfo_nodes=str(data.get("raw_sinfo_nodes") or ""),
            probed_at=str(data.get("probed_at") or ""),
            errors=list(data.get("errors") or []),
        )


@dataclass
class ClusterProfile:
    """User-level reusable remote cluster definition (no passwords)."""

    id: str
    name: str
    host: str
    username: str
    port: int = 22
    identity_file: str = ""
    scheduler: str = "slurm"
    submit_root: str = "/data/$USER/gatewizard"
    scratch_root: str = "$SCRATCH_DIR"
    workdir_strategy: str = "scratch_job_id"
    purge_modules: bool = True
    mail_user: str = ""
    mail_type: str = "NONE"
    extra_sbatch_lines: List[str] = field(default_factory=list)
    batch_template: Optional[str] = None
    module_hints: Dict[str, List[str]] = field(default_factory=dict)
    last_probe: Optional[Dict[str, Any]] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ClusterProfile":
        strategy = str(data.get("workdir_strategy") or "scratch_job_id")
        if strategy not in WORKDIR_STRATEGIES:
            strategy = "scratch_job_id"
        scheduler = str(data.get("scheduler") or "slurm").lower()
        if scheduler not in SCHEDULERS:
            scheduler = "slurm"
        return cls(
            id=str(data.get("id") or ""),
            name=str(data.get("name") or data.get("id") or "Cluster"),
            host=str(data.get("host") or ""),
            username=str(data.get("username") or ""),
            port=int(data.get("port") or 22),
            identity_file=str(data.get("identity_file") or ""),
            scheduler=scheduler,
            submit_root=str(data.get("submit_root") or "/data/$USER/gatewizard"),
            scratch_root=str(data.get("scratch_root") or "$SCRATCH_DIR"),
            workdir_strategy=strategy,
            purge_modules=bool(data.get("purge_modules", True)),
            mail_user=str(data.get("mail_user") or ""),
            mail_type=str(data.get("mail_type") or "NONE"),
            extra_sbatch_lines=[
                str(x) for x in (data.get("extra_sbatch_lines") or []) if str(x).strip()
            ],
            batch_template=data.get("batch_template"),
            module_hints={
                str(k): [str(x) for x in (v or [])]
                for k, v in (data.get("module_hints") or {}).items()
            },
            last_probe=data.get("last_probe")
            if isinstance(data.get("last_probe"), dict)
            else None,
        )


@dataclass
class BatchScriptRequest:
    """Inputs for rendering a scheduler batch wrapper."""

    job_name: str
    cpus: int = 8
    gpus: int = 0
    mem: str = ""
    time_limit: str = "24:00:00"
    partition: str = ""
    modules: List[str] = field(default_factory=list)
    purge_modules: bool = True
    mail_user: str = ""
    mail_type: str = "NONE"
    extra_sbatch_lines: List[str] = field(default_factory=list)
    workdir_strategy: str = "scratch_job_id"
    scratch_root: str = "$SCRATCH_DIR"
    job_folder_name: str = ""
    run_command: str = "bash run_equilibration_cluster.sh"
    template: Optional[str] = None
    nodelist: str = ""
    constraint: str = ""


@dataclass
class RemoteJobHandle:
    scheduler: str
    job_id: str
    state: str = ""
    name: str = ""
    partition: str = ""
    reason: str = ""
    node_list: str = ""
    elapsed: str = ""
    cpus: int = 0
    raw: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)
