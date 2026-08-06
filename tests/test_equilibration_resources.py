"""Tests for equilibration resource inference (v2 per-stage)."""

import json
from pathlib import Path

from gatewizard.utils.equilibration_resources import (
    RESOURCES_VERSION,
    aggregate_slurm_resources,
    enrich_resources_display,
    infer_equilibration_resources,
    resolve_all_stage_resources,
    resolve_compute_resources_from_eq_dir,
    resolve_compute_resources_from_stages,
    resolve_stage_resources,
    slurm_resources_from_eq_dir,
    write_equilibration_resources,
)


def test_resolve_stage_resources_minimization_forces_cpu() -> None:
    resolved = resolve_stage_resources(
        {"name": "Minimization", "stage_kind": "minimization", "use_gpu": True},
        {"use_gpu": True, "num_gpus": 1},
        engine="amber",
    )
    assert resolved["use_gpu"] is False
    assert resolved["num_gpus"] == 0
    assert resolved["cpu_cores"] == 6


def test_engine_resource_profile_amber_cpu_eq_gpu_prod() -> None:
    from gatewizard.utils.equilibration_resources import engine_resource_profile

    profile = engine_resource_profile("amber")
    assert profile["equilibration"]["cpu_cores"] == 6
    assert profile["equilibration"]["use_gpu"] is False
    assert profile["production"]["cpu_cores"] == 1
    assert profile["production"]["use_gpu"] is True


def test_engine_resource_profile_gromacs_gpu_md() -> None:
    from gatewizard.utils.equilibration_resources import engine_resource_profile

    profile = engine_resource_profile("gromacs")
    assert profile["equilibration"]["cpu_cores"] == 6
    assert profile["equilibration"]["use_gpu"] is True
    assert profile["production"]["cpu_cores"] == 6


def test_aggregate_slurm_resources_amber_engine_defaults() -> None:
    stages = resolve_all_stage_resources(
        [
            {"name": "Minimization", "stage_kind": "minimization"},
            {"name": "Equilibration 1", "stage_kind": "equilibration", "resources_inherit": True},
            {"name": "Production", "stage_kind": "production"},
        ],
        engine="amber",
    )
    slurm = aggregate_slurm_resources(stages)
    assert slurm["cpu_cores"] == 6
    assert slurm["num_gpus"] == 1
    assert stages[0]["use_gpu"] is False
    assert stages[1]["use_gpu"] is False
    assert stages[2]["use_gpu"] is True


def test_aggregate_slurm_resources_max_cpu_and_gpu() -> None:
    stages = resolve_all_stage_resources(
        [
            {
                "name": "Minimization",
                "stage_kind": "minimization",
                "cpu_cores": 4,
                "use_gpu": False,
            },
            {
                "name": "Equilibration 1",
                "stage_kind": "equilibration",
                "resources_inherit": True,
            },
            {
                "name": "Production",
                "stage_kind": "production",
                "num_gpus": 1,
                "use_gpu": True,
            },
        ],
        {"cpu_cores": 1, "num_gpus": 1, "use_gpu": True},
        engine="amber",
    )
    slurm = aggregate_slurm_resources(stages)
    assert slurm["cpu_cores"] == 4
    assert slurm["num_gpus"] == 1
    assert slurm["use_gpu"] is True


def test_write_equilibration_resources_v2(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    stages = [
        {
            "name": "Minimization",
            "stage_kind": "minimization",
            "cpu_cores": 4,
            "use_gpu": False,
        },
        {
            "name": "Equilibration 1",
            "stage_kind": "equilibration",
            "resources_inherit": True,
        },
    ]
    path = write_equilibration_resources(
        eq,
        "amber",
        stages,
        compute_defaults={"cpu_cores": 1, "num_gpus": 1, "use_gpu": True},
        stems=["step0_minimization", "step1_equilibration"],
    )
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["version"] == RESOURCES_VERSION
    assert len(data["stages"]) == 2
    assert data["stages"][0]["stem"] == "step0_minimization"
    assert data["stages"][0]["use_gpu"] is False
    assert data["slurm"]["cpu_cores"] == 4


def test_infer_equilibration_resources_from_protocol_summary(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    (eq / "protocol_summary.json").write_text(
        json.dumps(
            {
                "stages": {
                    "step1": {
                        "cpu_cores": 4,
                        "gpu_id": 0,
                        "num_gpus": 1,
                        "use_gpu": True,
                    },
                    "step2": {
                        "cpu_cores": 8,
                        "gpu_id": 1,
                        "num_gpus": 2,
                        "use_gpu": True,
                    },
                }
            }
        )
    )

    resources = infer_equilibration_resources(eq, "namd")

    assert resources["cpu_cores_min"] == 4
    assert resources["cpu_cores_max"] == 8
    assert resources["gpu_id_min"] == 0
    assert resources["gpu_id_max"] == 1
    assert resources["num_gpus"] == 2
    assert resources["use_gpu"] is True


def test_write_equilibration_resources_roundtrip(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    stages = [
        {
            "name": "Equilibration 1",
            "stage_kind": "equilibration",
            "cpu_cores": 2,
            "gpu_id": 0,
            "num_gpus": 1,
            "use_gpu": False,
        }
    ]
    write_equilibration_resources(eq, "namd", stages)

    resources = infer_equilibration_resources(eq, "namd")

    assert resources["cpu_cores_min"] == 2
    assert resources["use_gpu"] is False
    assert resources["version"] == RESOURCES_VERSION


def test_resolve_compute_resources_from_stages() -> None:
    flat = resolve_compute_resources_from_stages(
        [
            {
                "name": "Minimization",
                "stage_kind": "minimization",
                "cpu_cores": 4,
                "use_gpu": False,
            },
            {
                "name": "Equilibration 1",
                "stage_kind": "equilibration",
                "cpu_cores": 1,
                "gpu_id": 0,
                "num_gpus": 2,
                "use_gpu": True,
            },
        ],
        {"cpu_cores": 1, "num_gpus": 1, "use_gpu": True},
        engine="amber",
    )
    assert flat["cpu_cores"] == 4
    assert flat["use_gpu"] is True
    assert flat["gpu_id"] == 0
    assert flat["num_gpus"] == 2


def test_resolve_compute_resources_from_eq_dir(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    write_equilibration_resources(
        eq,
        "gromacs",
        [
            {
                "name": "Minimization",
                "stage_kind": "minimization",
                "cpu_cores": 4,
                "use_gpu": False,
            },
            {
                "name": "Equilibration 1",
                "stage_kind": "equilibration",
                "resources_inherit": True,
            },
        ],
        compute_defaults={"cpu_cores": 1, "num_gpus": 1, "use_gpu": True},
        stems=["step0_minimization", "step1_equilibration"],
    )
    flat = resolve_compute_resources_from_eq_dir(eq)
    assert flat["cpu_cores"] == 4
    assert flat["use_gpu"] is True
    assert flat["num_gpus"] == 1
    assert len(flat["stages"]) == 2


def test_slurm_resources_from_eq_dir(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    write_equilibration_resources(
        eq,
        "amber",
        [
            {"name": "Minimization", "stage_kind": "minimization"},
            {"name": "Production", "stage_kind": "production", "num_gpus": 1, "use_gpu": True},
        ],
        stems=["step0_minimization", "step7_production"],
    )
    slurm = slurm_resources_from_eq_dir(eq)
    assert slurm["cpu_cores"] == 6
    assert slurm["num_gpus"] == 1


def test_enrich_resources_display_summary() -> None:
    payload = enrich_resources_display(
        {
            "version": 2,
            "engine": "amber",
            "stages": [
                {"name": "Minimization", "stage_kind": "minimization", "cpu_cores": 4, "use_gpu": False},
                {"name": "Equilibration 1", "stage_kind": "equilibration", "use_gpu": True, "num_gpus": 1},
            ],
            "slurm": {"cpu_cores": 4, "num_gpus": 1, "use_gpu": True},
        }
    )
    assert "Min CPU×4" in payload["summary"]
    assert payload["cpu_cores_max"] == 4
