"""Tests for remote cluster utilities (parsers, templates, metadata)."""

from __future__ import annotations

from pathlib import Path

import pytest

from gatewizard.utils.cluster import (
    BatchScriptRequest,
    ClusterProfile,
    default_template_for_strategy,
    expand_remote_path,
    get_scheduler,
    group_engine_modules,
    join_remote,
    parse_module_avail,
    parse_sbatch_output,
    parse_sinfo,
    parse_squeue_me,
    prefer_gpu_modules,
    render_batch_script,
    update_execution_fields,
    write_execution_metadata,
)

FIXTURES = Path(__file__).parent / "cluster_fixtures"


def test_parse_module_avail_ohpc_finds_md_engines():
    text = (FIXTURES / "module_avail_ohpc.txt").read_text(encoding="utf-8")
    packages = parse_module_avail(text)
    names = {p.full_name for p in packages}
    assert "md/namd/3.0b6+cuda" in names
    assert "md/gromacs/2024.1+cuda" in names
    assert "cuda/12.3.2" in names
    grouped = group_engine_modules(packages)
    assert any(p.full_name.startswith("md/namd") for p in grouped["namd"])
    assert any("gromacs" in p.full_name for p in grouped["gromacs"])
    gpu_first = prefer_gpu_modules(grouped["namd"], want_gpu=True)
    assert "cuda" in gpu_first[0].full_name.lower()


def test_parse_module_avail_rejects_help_garbage():
    text = "No modules found\nUse spider to search\nany Key to search"
    assert parse_module_avail(text) == []


def test_prefer_partitions_gpu_first():
    from gatewizard.utils.cluster import prefer_partitions
    from gatewizard.utils.cluster.types import PartitionInfo

    parts = [
        PartitionInfo(name="normal", max_gpus=0),
        PartitionInfo(name="gpu", max_gpus=2),
        PartitionInfo(name="debug", max_gpus=0),
    ]
    ranked = prefer_partitions(parts, want_gpu=True)
    assert ranked[0].name == "gpu"
    ranked_cpu = prefer_partitions(parts, want_gpu=False)
    assert ranked_cpu[0].name == "normal"


def test_parse_sinfo_nodes_and_prefer_gpu():
    from gatewizard.utils.cluster.resources import parse_sinfo_nodes, prefer_nodes

    text = (FIXTURES / "sinfo_nodes_ohpc.txt").read_text(encoding="utf-8")
    nodes = parse_sinfo_nodes(text)
    names = {n.name for n in nodes}
    assert names >= {"wc01", "gpu01", "gpu02"}
    gpu01 = next(n for n in nodes if n.name == "gpu01")
    assert gpu01.gpus == 4
    assert gpu01.partition == "gpu"
    ranked = prefer_nodes(nodes, want_gpu=True, partition="gpu")
    assert ranked[0].name == "gpu01"
    assert all(n.name != "gpu03" for n in ranked)  # drain filtered
    cpuish = prefer_nodes(nodes, want_gpu=False, partition="normal")
    assert cpuish[0].name in {"wc01", "wc02"}


def test_batch_script_includes_nodelist():
    script = render_batch_script(
        BatchScriptRequest(
            job_name="eq",
            cpus=4,
            gpus=1,
            partition="gpu",
            nodelist="gpu01",
            modules=["md/namd/3.0.1+cuda"],
        )
    )
    assert "#SBATCH --nodelist=gpu01" in script
    assert "#SBATCH --partition=gpu" in script
    assert "#SBATCH --gpus=1" in script


def test_cluster_engine_executable_maps_wsl_paths():
    from gatewizard.utils.equilibration_cluster_script import (
        cluster_engine_executable,
        script_has_wsl_or_windows_path,
        stamp_cluster_run_script_header,
    )

    assert (
        cluster_engine_executable(
            "namd", "/mnt/c/software/namd/NAMD_3.0.2_Linux-x86_64-multicore-CUDA/namd3"
        )
        == "namd3"
    )
    assert cluster_engine_executable("gromacs", "/usr/local/gromacs/bin/gmx") == "gmx"
    assert cluster_engine_executable("openmm", "/home/user/miniconda/bin/python") == "python"
    assert cluster_engine_executable("amber", "pmemd.cuda") == "pmemd.cuda"
    assert cluster_engine_executable("amber", "/opt/amber/bin/pmemd.cuda") == "pmemd.cuda"
    # Cluster GPU submenu / protocol use_gpu must upgrade plain pmemd → pmemd.cuda
    assert cluster_engine_executable("amber", "pmemd", use_gpu=True) == "pmemd.cuda"
    assert (
        cluster_engine_executable("amber", "pmemd.MPI", use_gpu=True) == "pmemd.cuda.MPI"
    )
    assert cluster_engine_executable("amber", "pmemd.cuda", use_gpu=False) == "pmemd"
    stamped = stamp_cluster_run_script_header("#!/bin/bash\necho hi\n")
    assert "Cluster runner" in stamped
    assert stamped.startswith("#!/bin/bash")
    assert script_has_wsl_or_windows_path('NAMD="/mnt/c/software/namd3"')
    assert not script_has_wsl_or_windows_path('NAMD="namd3"')


def test_namd_writes_local_and_cluster_run_scripts(tmp_path: Path):
    from gatewizard.tools.equilibration import NAMDEquilibrationManager
    from gatewizard.utils.equilibration_cluster_script import CLUSTER_RUN_SCRIPT

    mgr = NAMDEquilibrationManager(tmp_path, namd_executable="namd3")
    protocols = {
        "eq1": {
            "name": "Equilibration 1",
            "steps": 100,
            "timestep": 1.0,
            "use_gpu": True,
            "cpu_cores": 4,
            "gpu_id": 0,
            "num_gpus": 1,
        }
    }
    local = tmp_path / "run_equilibration.sh"
    local.write_text(
        mgr.generate_run_script(
            protocols, "/mnt/c/software/namd/namd3"
        )
    )
    from gatewizard.utils.equilibration_cluster_script import (
        cluster_engine_executable,
        write_cluster_run_script,
    )

    write_cluster_run_script(
        tmp_path,
        mgr.generate_run_script(
            protocols, cluster_engine_executable("namd", "/mnt/c/software/namd/namd3")
        ),
    )
    cluster = tmp_path / CLUSTER_RUN_SCRIPT
    assert local.is_file() and cluster.is_file()
    assert "/mnt/c/" in local.read_text()
    assert 'NAMD="namd3"' in cluster.read_text()
    assert "Cluster runner" in cluster.read_text()
    from gatewizard.utils.equilibration_cluster_script import resolve_cluster_launch_script

    assert resolve_cluster_launch_script(tmp_path).name == CLUSTER_RUN_SCRIPT


def test_parse_sinfo_and_squeue():
    sinfo = parse_sinfo((FIXTURES / "sinfo_ohpc.txt").read_text(encoding="utf-8"))
    assert {p.name for p in sinfo} >= {"normal", "gpu", "debug"}
    normal = next(p for p in sinfo if p.name == "normal")
    assert normal.max_gpus == 2
    assert normal.cpus_per_node == 64

    jobs = parse_squeue_me((FIXTURES / "squeue_me.txt").read_text(encoding="utf-8"))
    assert jobs[0]["job_id"] == "238"
    assert jobs[0]["state"] == "RUNNING"
    assert jobs[1]["state"] == "PENDING"


def test_render_scratch_job_id_template():
    script = render_batch_script(
        BatchScriptRequest(
            job_name="eq_popc",
            cpus=8,
            gpus=1,
            time_limit="12:00:00",
            modules=["cuda/12.3.2", "md/namd/3.0b6+cuda"],
            purge_modules=True,
            mail_user="user@example.com",
            mail_type="ALL",
            workdir_strategy="scratch_job_id",
            scratch_root="$SCRATCH_DIR",
        )
    )
    assert "#SBATCH -J eq_popc" in script
    assert "#SBATCH --gpus=1" in script
    assert "ml purge" in script
    assert "module load md/namd/3.0b6+cuda" in script
    assert 'workdir="$SCRATCH_DIR/$SLURM_JOB_ID"' in script
    assert "bash run_equilibration_cluster.sh" in script
    assert "rsync -a" in script
    # Mid-run log sync so Watching can show progress while scratch is node-local
    assert "include='step*.log'" in script
    assert "_GW_SYNC_PID" in script


def test_render_run_in_place():
    script = render_batch_script(
        BatchScriptRequest(
            job_name="localish",
            cpus=4,
            workdir_strategy="run_in_place",
            modules=["md/gromacs/2024.1+cuda"],
            purge_modules=False,
        )
    )
    assert "SLURM_SUBMIT_DIR" in script
    assert "SCRATCH_DIR" not in script or "workdir=" not in script
    assert "ml purge" not in script


def test_custom_template_placeholders():
    custom = default_template_for_strategy("custom_template")
    assert "{{job_name}}" in custom
    script = render_batch_script(
        BatchScriptRequest(
            job_name="x",
            template="# {{job_name}}\n{{run_command}}\n",
            workdir_strategy="custom_template",
        )
    )
    assert script.startswith("# x\n")
    assert "bash run_equilibration_cluster.sh" in script


def test_paths_and_profile():
    assert expand_remote_path("/data/$USER/gw", username="alice") == "/data/alice/gw"
    assert join_remote("/data/alice", "jobs", "eq1") == "/data/alice/jobs/eq1"
    profile = ClusterProfile.from_dict(
        {
            "id": "demo",
            "name": "Demo cluster",
            "host": "hpc.example.edu",
            "username": "alice",
            "workdir_strategy": "scratch_job_id",
        }
    )
    assert profile.scheduler == "slurm"
    assert profile.port == 22


def test_slurm_adapter_submit_parse():
    adapter = get_scheduler("slurm")
    assert adapter.submit_command("run_equilibration.slurm") == [
        "sbatch",
        "run_equilibration.slurm",
    ]
    assert parse_sbatch_output("Submitted batch job 238\n") == "238"
    handles = adapter.parse_status(
        (FIXTURES / "squeue_me.txt").read_text(encoding="utf-8"), job_id="238"
    )
    assert len(handles) == 1
    assert handles[0].state == "RUNNING"


def test_parse_squeue_pipe_format():
    from gatewizard.utils.cluster import parse_squeue_me, summarize_node_gpu_label

    text = (
        "4947|RUNNING|namd_nvt2|normal|cn01|cn01|7:56|12\n"
        "241|PENDING|wait|normal|(Resources)||0:00|8\n"
    )
    rows = parse_squeue_me(text)
    assert rows[0]["node_list"] == "cn01"
    assert rows[0]["cpus"] == 12
    assert rows[1]["state"] == "PENDING"
    assert rows[1]["node_list"] == ""
    assert summarize_node_gpu_label("gpu:l4:2,gpu:2080ti:1") == "L4 / 2080 Ti"


def test_gpus_from_gres_sums_multi_type():
    from gatewizard.utils.cluster.resources import _gpus_from_gres, canonicalize_slurm_state

    assert _gpus_from_gres("gpu:2080ti:2,gpu:3090:1") == 3
    assert _gpus_from_gres("gpu:l4:2") == 2
    assert _gpus_from_gres("gpu:2") == 2
    assert canonicalize_slurm_state("CANCELLED by 1002") == "CANCELLED"
    assert canonicalize_slurm_state("COMPLETED+") == "COMPLETED"
    assert canonicalize_slurm_state("RUNNING") == "RUNNING"


def test_ensure_amber_cluster_runner_for_gpus(tmp_path: Path):
    from gatewizard.utils.equilibration_cluster_script import (
        CLUSTER_RUN_SCRIPT,
        ensure_amber_cluster_runner_for_gpus,
    )

    (tmp_path / "step0_minimization.mdin").write_text("&cntrl\n/\n", encoding="utf-8")
    (tmp_path / "step1_equilibration.mdin").write_text("&cntrl\n/\n", encoding="utf-8")
    (tmp_path / "run_equilibration.sh").write_text(
        '#!/bin/bash\nAMBER="pmemd"\nMINI_AMBER="pmemd"\n'
        'PRMTOP="system.prmtop"\nINPCRD="system.inpcrd"\n',
        encoding="utf-8",
    )
    (tmp_path / "equilibration_resources.json").write_text(
        '{"use_gpu": true, "num_gpus": 1, "cpu_cores": 4, "gpu_id": 0}',
        encoding="utf-8",
    )

    assert ensure_amber_cluster_runner_for_gpus(tmp_path, gpus=1)
    cluster = (tmp_path / CLUSTER_RUN_SCRIPT).read_text(encoding="utf-8")
    assert 'AMBER="pmemd.cuda"' in cluster
    assert 'MINI_AMBER="pmemd"' in cluster
    assert "CUDA_VISIBLE_DEVICES" in cluster
    assert "GPU: Yes" in cluster

    assert ensure_amber_cluster_runner_for_gpus(tmp_path, gpus=0)
    cluster_cpu = (tmp_path / CLUSTER_RUN_SCRIPT).read_text(encoding="utf-8")
    assert 'AMBER="pmemd"' in cluster_cpu
    assert "GPU: No" in cluster_cpu


def test_read_batch_script_resources(tmp_path: Path):
    from gatewizard.utils.cluster import enrich_execution_resources, read_batch_script_resources

    (tmp_path / "run_equilibration.slurm").write_text(
        "#!/bin/bash\n#SBATCH -c 12\n#SBATCH --gpus=1\n",
        encoding="utf-8",
    )
    assert read_batch_script_resources(tmp_path)["cpus"] == 12
    assert read_batch_script_resources(tmp_path)["gpus"] == 1
    enriched = enrich_execution_resources(tmp_path, {"mode": "remote"})
    assert enriched["allocated_cpus"] == 12
    assert enriched["resources"]["cpus"] == 12
    assert enriched["resources"]["gpus"] == 1


def test_execution_metadata_roundtrip(tmp_path: Path):
    write_execution_metadata(
        tmp_path,
        {
            "mode": "remote",
            "cluster_id": "demo",
            "scheduler": "slurm",
            "remote_path": "/data/u/gw/job",
        },
    )
    execution = update_execution_fields(
        tmp_path, scheduler_job_id="238", last_remote_state="RUNNING"
    )
    assert execution["scheduler_job_id"] == "238"
    assert execution["mode"] == "remote"
    text = (tmp_path / "equilibration_job.json").read_text(encoding="utf-8")
    assert '"scheduler_job_id": "238"' in text


def test_archive_previous_run_outputs_and_newer_than(tmp_path: Path):
    from gatewizard.utils.equilibration_failure import (
        archive_previous_run_outputs,
        find_equilibration_failure,
        remote_job_is_active,
    )
    from gatewizard.utils.namd_analysis import get_equilibration_progress
    import time

    old = tmp_path / "step1_equilibration.log"
    old.write_text(
        "Charm++> Running in Multicore mode\n"
        "FATAL ERROR: CUDA error cudaGetDeviceCount(&deviceCount) "
        "in file src/DeviceCUDA.C: CUDA driver is a stub library\n",
        encoding="utf-8",
    )
    (tmp_path / "namd_nvt.4931.out").write_text("", encoding="utf-8")
    msg, src = find_equilibration_failure(tmp_path)
    assert msg and src == "step1_equilibration.log"
    assert get_equilibration_progress(tmp_path)["equilibration_1"].status == "error"

    time.sleep(0.05)
    cutoff = time.time()
    assert find_equilibration_failure(tmp_path, newer_than=cutoff + 10) == (None, None)

    n = archive_previous_run_outputs(tmp_path)
    assert n >= 1
    assert not old.exists()
    assert find_equilibration_failure(tmp_path) == (None, None)
    assert remote_job_is_active({"mode": "remote", "last_remote_state": "RUNNING"})
    assert not remote_job_is_active({"mode": "remote", "last_remote_state": "FAILED"})


def test_archive_preserves_namd_resume_checkpoints(tmp_path: Path) -> None:
    from gatewizard.utils.equilibration_failure import archive_previous_run_outputs

    stem = "step1_equilibration"
    (tmp_path / f"{stem}.conf").write_text("steps 100")
    (tmp_path / f"{stem}.coor").write_text("coor")
    (tmp_path / f"{stem}.log").write_text("End of program\n")
    failed = tmp_path / "step2_equilibration.log"
    failed.write_text("FATAL ERROR: CUDA\n", encoding="utf-8")

    n = archive_previous_run_outputs(tmp_path, engine="namd")
    assert n >= 1
    assert (tmp_path / f"{stem}.coor").is_file()
    assert (tmp_path / f"{stem}.log").is_file()
    assert not failed.exists()


def test_resolve_compute_node_and_midrun_sync(monkeypatch):
    from gatewizard.utils.cluster import midrun as midrun_mod

    calls = []

    def fake_run(session_id, cmd, timeout=60):
        calls.append(cmd)
        if "squeue" in cmd:
            return 0, "cn01\n", ""
        if "GW_MIDRUN" in cmd or "tar czf" in cmd:
            assert "cn01" in cmd
            assert "/scratch/testuser/4944" in cmd
            assert "/home/testuser/job" in cmd
            return 0, "GW_MIDRUN_DONE\n", ""
        return 0, "", ""

    monkeypatch.setattr(midrun_mod, "run_remote", fake_run)
    assert midrun_mod.resolve_compute_node("s1", "4944") == "cn01"
    ok, msg = midrun_mod.sync_scratch_progress_to_submit(
        "s1",
        job_id="4944",
        node="cn01",
        scratch_root="/scratch/testuser",
        remote_submit_dir="/home/testuser/job",
    )
    assert ok
    assert "cn01:/scratch/testuser/4944" in msg
    assert any("tar czf" in c for c in calls)


def test_midrun_missing_scratch(monkeypatch):
    from gatewizard.utils.cluster import midrun as midrun_mod

    def fake_run(session_id, cmd, timeout=60):
        return 0, "GW_MIDRUN_MISSING\nGW_MIDRUN_DONE\n", ""

    monkeypatch.setattr(midrun_mod, "run_remote", fake_run)
    ok, msg = midrun_mod.sync_scratch_progress_to_submit(
        "s1",
        job_id="1",
        node="cn01",
        scratch_root="/scratch/testuser",
        remote_submit_dir="/home/testuser/job",
    )
    assert not ok
    assert "scratch not found" in msg


def test_parse_scratch_workdir_from_slurm(tmp_path: Path) -> None:
    from gatewizard.utils.cluster.midrun import (
        parse_scratch_workdir_from_slurm,
        resolve_scratch_job_dir,
    )

    slurm = tmp_path / "run_equilibration.slurm"
    slurm.write_text('workdir="/scratch/testuser/$SLURM_JOB_ID"\n')
    root, strategy = parse_scratch_workdir_from_slurm(tmp_path)
    assert root == "/scratch/testuser"
    assert strategy == "scratch_job_id"
    scratch = resolve_scratch_job_dir(
        job_id="4948",
        execution={},
        local_dir=tmp_path,
    )
    assert scratch == "/scratch/testuser/4948"


def test_stage_scratch_to_login(monkeypatch, tmp_path: Path) -> None:
    from gatewizard.utils.cluster import midrun as midrun_mod

    def fake_run(session_id, cmd, timeout=60):
        assert "cn01" in cmd
        assert "/scratch/testuser/4948" in cmd
        assert "/tmp/gw_pull_4948" in cmd
        return 0, "GW_PULL_DONE\n", ""

    monkeypatch.setattr(midrun_mod, "run_remote", fake_run)
    ok, msg = midrun_mod.stage_scratch_to_login(
        "s1",
        node="cn01",
        scratch_job_dir="/scratch/testuser/4948",
        staging_dir="/tmp/gw_pull_4948",
        full=True,
    )
    assert ok
    assert "cn01:/scratch/testuser/4948" in msg


def test_build_remote_submit_path() -> None:
    from gatewizard.utils.cluster.midrun import build_remote_submit_path

    assert (
        build_remote_submit_path(
            "/home/$USER/gw_jobs",
            "openmm_nvt",
            username="testuser",
        )
        == "/home/testuser/gw_jobs/openmm_nvt"
    )


def test_resolve_remote_job_dir_prefers_existing_sacct(monkeypatch) -> None:
    from gatewizard.utils.cluster import midrun as midrun_mod

    good = "/home/testuser/gw_jobs/openmm_nvt"
    bad = "/home/testuser/gw_jobs/openmm_nvt_wrong"

    monkeypatch.setattr(
        midrun_mod,
        "remote_path_is_dir",
        lambda _sid, path: path.rstrip("/") == good,
    )
    monkeypatch.setattr(
        midrun_mod,
        "resolve_slurm_workdir",
        lambda _sid, _jid: good,
    )

    path, source, tried = midrun_mod.resolve_remote_job_dir(
        "s1",
        stored_path=bad,
        job_id="4948",
        submit_root="/home/$USER/gw_jobs",
        username="testuser",
        job_folder="openmm_nvt",
    )
    assert path == good
    assert source == "sacct WorkDir"
    assert bad in tried
    assert good in tried


def test_parse_rsync_progress_line() -> None:
    from gatewizard.utils.cluster.ssh import parse_rsync_progress_line

    evt = parse_rsync_progress_line(
        "  12,345,678  45%  10.50MB/s    0:01:23 (xfr#5, to-chk=10/100)"
    )
    assert evt is not None
    assert evt["percent"] == 45
    assert evt["bytes"] == 12345678
    assert evt["speed"] == "10.50MB/s"
    assert evt["phase"] == "sync"
    assert parse_rsync_progress_line("random noise") is None


def test_format_byte_size_and_local_dir_byte_size(tmp_path: Path) -> None:
    from gatewizard.utils.cluster.ssh import format_byte_size, local_dir_byte_size

    assert format_byte_size(500) == "500 B"
    assert "KB" in format_byte_size(2048)
    (tmp_path / "a.log").write_bytes(b"x" * 1000)
    (tmp_path / "skip.pid").write_bytes(b"y" * 50)
    assert local_dir_byte_size(tmp_path, excludes=["*.pid"]) == 1000

