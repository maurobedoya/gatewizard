"""Tests for equilibration resource inference."""

import json
from pathlib import Path

from gatewizard.utils.equilibration_resources import (
    infer_equilibration_resources,
    resolve_compute_resources_from_eq_dir,
    resolve_compute_resources_from_stages,
    write_equilibration_resources,
)


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
    stages = [{"cpu_cores": 2, "gpu_id": 0, "num_gpus": 1, "use_gpu": False}]
    write_equilibration_resources(eq, "namd", stages)

    resources = infer_equilibration_resources(eq, "namd")

    assert resources["cpu_cores_min"] == 2
    assert resources["use_gpu"] is False


def test_resolve_compute_resources_from_stages() -> None:
    flat = resolve_compute_resources_from_stages(
        [
            {"cpu_cores": 4, "gpu_id": 0, "num_gpus": 1, "use_gpu": True},
            {"cpu_cores": 8, "gpu_id": 0, "num_gpus": 2, "use_gpu": True},
        ]
    )
    assert flat["cpu_cores"] == 8
    assert flat["use_gpu"] is True
    assert flat["gpu_id"] == 0
    assert flat["num_gpus"] == 2


def test_resolve_compute_resources_from_eq_dir(tmp_path: Path) -> None:
    eq = tmp_path / "job"
    eq.mkdir()
    write_equilibration_resources(
        eq,
        "gromacs",
        [{"cpu_cores": 4, "gpu_id": 0, "num_gpus": 1, "use_gpu": True}],
    )
    flat = resolve_compute_resources_from_eq_dir(eq)
    assert flat["cpu_cores"] == 4
    assert flat["use_gpu"] is True
    assert flat["gpu_id"] == 0
    assert flat["num_gpus"] == 1
