"""Tests for equilibration job metadata read/write and inference."""

from __future__ import annotations

import json
from pathlib import Path

from gatewizard.utils.equilibration_job_metadata import (
    JOB_METADATA_FILE,
    infer_equilibration_job_metadata,
    write_equilibration_job_metadata,
)


def test_write_and_read_equilibration_job_metadata(tmp_path: Path) -> None:
    eq = tmp_path / "03_equilibration_02_build_demo"
    eq.mkdir()
    protocol = {
        "name": "Demo",
        "description": "Test protocol",
        "stages": [
            {
                "name": "Equilibration 1",
                "ensemble": "NPT",
                "constraints": [
                    {
                        "name": "Protein backbone",
                        "force_constant": 10.0,
                        "selection": "protein_backbone",
                    }
                ],
            }
        ],
    }
    write_equilibration_job_metadata(
        eq,
        input_dir=str(tmp_path / "02_build_demo"),
        ensemble="npt",
        protocol=protocol,
        engine="namd",
    )

    meta = infer_equilibration_job_metadata(eq, working_dir=tmp_path)
    assert meta["input_dir"] == str((tmp_path / "02_build_demo").resolve())
    assert meta["ensemble"] == "NPT"
    assert meta["protocol"]["name"] == "Demo"
    assert (eq / JOB_METADATA_FILE).is_file()


def test_infer_input_dir_from_folder_name(tmp_path: Path) -> None:
    build = tmp_path / "02_build_demo"
    build.mkdir()
    (build / "system.prmtop").write_text("x", encoding="utf-8")

    eq = tmp_path / "03_equilibration_02_build_demo_openmm"
    eq.mkdir()

    meta = infer_equilibration_job_metadata(eq, working_dir=tmp_path)
    assert meta["input_dir"] == str(build.resolve())
    assert meta["ensemble"] is None
    assert meta["protocol"] is None


def test_infer_protocol_from_namd_summary(tmp_path: Path) -> None:
    eq = tmp_path / "03_equilibration_demo"
    eq.mkdir()
    (eq / "protocol_summary.json").write_text(
        json.dumps(
            {
                "protocol_name": "NPT Equilibration Protocol",
                "scheme_type": "NPT",
                "stages": {
                    "Equilibration 1": {
                        "name": "Equilibration 1",
                        "ensemble": "NPT",
                        "time_ns": 0.125,
                        "steps": 125000,
                        "timestep": 1.0,
                        "temperature": 303.15,
                        "constraints": {"protein_backbone": 10.0, "water": 0.0},
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    meta = infer_equilibration_job_metadata(eq)
    assert meta["ensemble"] == "NPT"
    assert meta["protocol"] is not None
    stage = meta["protocol"]["stages"][0]
    assert stage["name"] == "Equilibration 1"
    assert len(stage["constraints"]) == 2
    assert stage["constraints"][0]["selection"] == "protein_backbone"
