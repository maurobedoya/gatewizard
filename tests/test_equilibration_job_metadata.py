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


def test_infer_input_dir_by_topology_match(tmp_path: Path) -> None:
    build = tmp_path / "01_build_test"
    build.mkdir()
    prmtop = "same-topology-bytes"
    (build / "system.prmtop").write_text(prmtop, encoding="utf-8")
    (build / "system.inpcrd").write_text("coords", encoding="utf-8")

    eq = tmp_path / "namd_nvt2"
    eq.mkdir()
    (eq / "run_equilibration.sh").write_text("#!/bin/bash\nnamd3\n", encoding="utf-8")
    (eq / "system.prmtop").write_text(prmtop, encoding="utf-8")
    (eq / "system.inpcrd").write_text("coords", encoding="utf-8")

    meta = infer_equilibration_job_metadata(eq, working_dir=tmp_path, heal=False)
    assert meta["input_dir"] == str(build.resolve())


def test_infer_input_dir_from_gromacs_gro_atom_count(tmp_path: Path) -> None:
    build = tmp_path / "01_build_test"
    build.mkdir()
    (build / "system.prmtop").write_text("top", encoding="utf-8")
    (build / "system.inpcrd").write_text("title\n3\n", encoding="utf-8")

    eq = tmp_path / "gromacs_nvt"
    eq.mkdir()
    (eq / "run_equilibration.sh").write_text("#!/bin/bash\ngmx mdrun\n", encoding="utf-8")
    (eq / "system.gro").write_text("title\n3\n", encoding="utf-8")

    meta = infer_equilibration_job_metadata(eq, working_dir=tmp_path, heal=False)
    assert meta["input_dir"] == str(build.resolve())


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
    assert "protein_backbone" in meta["protocol"]["selections"]


def test_infer_protocol_from_gromacs_mdps(tmp_path: Path) -> None:
    eq = tmp_path / "gromacs_nvt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps({"execution": {"mode": "remote"}}), encoding="utf-8"
    )
    (eq / "step0_minimization.mdp").write_text(
        "integrator = steep\nnsteps = 10000\n", encoding="utf-8"
    )
    (eq / "step1_equilibration.mdp").write_text(
        "integrator = md\nnsteps = 125000\ndt = 0.001\nref_t = 303.15\npcoupl = no\n",
        encoding="utf-8",
    )
    (eq / "step7_production.mdp").write_text(
        "integrator = md\nnsteps = 100000000\ndt = 0.002\nref_t = 303.15\npcoupl = no\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq)
    assert meta["protocol"] is not None
    names = [s["name"] for s in meta["protocol"]["stages"]]
    assert names[0] == "Minimization"
    assert "Equilibration 1" in names
    assert names[-1] == "Production"
    assert meta["ensemble"] == "NVT"
    prod = meta["protocol"]["stages"][-1]
    assert prod["time_ns"] == 200.0
    assert prod["steps"] == 100000000
    assert prod["timestep"] == 2.0


def test_infer_protocol_production_time_from_amber_mdin(tmp_path: Path) -> None:
    eq = tmp_path / "a_npt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps({"execution": {"mode": "remote"}, "engine": "amber"}),
        encoding="utf-8",
    )
    (eq / "step7_production.mdin").write_text(
        "&cntrl\n  imin=0, nstlim=100000000, dt=0.002, temp0=303.15, ntp=1,\n/\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    assert meta["protocol"] is not None
    prod = next(s for s in meta["protocol"]["stages"] if s["name"] == "Production")
    assert prod["time_ns"] == 200.0
    assert prod["steps"] == 100000000


def test_infer_protocol_production_time_from_namd_conf(tmp_path: Path) -> None:
    eq = tmp_path / "namd_nvt2"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps({"execution": {"mode": "remote"}, "engine": "namd"}),
        encoding="utf-8",
    )
    (eq / "step7_production.conf").write_text(
        "set time 200.0;\nset tstep 2.0;\ntimestep 2.0\nrun [expr int($time * 1e6 / $tstep)]\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    assert meta["protocol"] is not None
    prod = next(s for s in meta["protocol"]["stages"] if s["name"] == "Production")
    assert prod["time_ns"] == 200.0
    assert prod["steps"] == 100000000


def test_patch_production_time_when_summary_omits_it(tmp_path: Path) -> None:
    eq = tmp_path / "g_npgt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps(
            {
                "engine": "gromacs",
                "protocol": {
                    "name": "NVT",
                    "stages": [
                        {"name": "Equilibration 1", "time_ns": 0.125, "steps": 125000, "timestep": 1.0},
                        {"name": "Production", "time_ns": 0, "steps": 0, "timestep": 2.0},
                    ],
                },
            }
        ),
        encoding="utf-8",
    )
    (eq / "step7_production.mdp").write_text(
        "nsteps = 100000000\ndt = 0.002\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    prod = next(s for s in meta["protocol"]["stages"] if s["name"] == "Production")
    assert prod["time_ns"] == 200.0
    assert prod["steps"] == 100000000


def test_patch_wrong_production_time_from_disk(tmp_path: Path) -> None:
    eq = tmp_path / "gromacs_nvt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps(
            {
                "engine": "gromacs",
                "protocol": {
                    "name": "NVT",
                    "stages": [
                        {"name": "Equilibration 1", "time_ns": 0.125, "steps": 125000, "timestep": 1.0},
                        {"name": "Production", "time_ns": 1.0, "steps": 500000, "timestep": 2.0},
                    ],
                },
            }
        ),
        encoding="utf-8",
    )
    (eq / "step7_production.mdp").write_text(
        "nsteps = 100000000\ndt = 0.002\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    prod = next(s for s in meta["protocol"]["stages"] if s["name"] == "Production")
    assert prod["time_ns"] == 200.0
    assert prod["steps"] == 100000000


def test_infer_protocol_from_openmm_inps(tmp_path: Path) -> None:
    eq = tmp_path / "openmm_nvt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps({"execution": {"mode": "remote"}}), encoding="utf-8"
    )
    (eq / "run_equilibration.sh").write_text(
        "#!/bin/bash\npython openmm_run.py\n", encoding="utf-8"
    )
    (eq / "step1_equilibration.inp").write_text(
        "nstep = 125000\ndt = 0.001\npcouple = no\n", encoding="utf-8"
    )
    (eq / "step7_production.inp").write_text(
        "nstep = 100000000\ndt = 0.002\npcouple = no\n", encoding="utf-8"
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    assert meta["engine"] == "openmm"
    assert meta["protocol"] is not None
    prod = next(s for s in meta["protocol"]["stages"] if s["name"] == "Production")
    assert prod["time_ns"] == 200.0
    assert prod["steps"] == 100000000


def test_infer_fills_protocol_when_job_json_is_execution_only(tmp_path: Path) -> None:
    """Cluster Watching may leave equilibration_job.json with only execution."""
    eq = tmp_path / "gromacs_nvt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps(
            {
                "execution": {
                    "mode": "remote",
                    "scheduler_job_id": "4945",
                    "last_remote_state": "RUNNING",
                }
            }
        ),
        encoding="utf-8",
    )
    (eq / "protocol_summary.json").write_text(
        json.dumps(
            {
                "protocol_name": "NVT Equilibration Protocol",
                "scheme_type": "NVT",
                "stages": {
                    "Equilibration 1": {
                        "name": "Equilibration 1",
                        "ensemble": "NVT",
                        "constraints": {"protein_backbone": 4.0},
                    }
                },
            }
        ),
        encoding="utf-8",
    )

    meta = infer_equilibration_job_metadata(eq)
    assert meta["ensemble"] == "NVT"
    assert meta["protocol"] is not None
    assert meta["protocol"]["stages"][0]["name"] == "Equilibration 1"
    # Healed back onto equilibration_job.json for subsequent Use in form.
    healed = json.loads((eq / JOB_METADATA_FILE).read_text(encoding="utf-8"))
    assert healed.get("execution", {}).get("scheduler_job_id") == "4945"
    assert healed.get("ensemble") == "NVT"
    assert isinstance(healed.get("protocol"), dict)


def test_infer_gromacs_positional_restraints_from_mdp(tmp_path: Path) -> None:
    eq = tmp_path / "gromacs_nvt"
    eq.mkdir()
    (eq / JOB_METADATA_FILE).write_text(
        json.dumps(
            {
                "engine": "gromacs",
                "protocol": {
                    "name": "NVT",
                    "stages": [
                        {
                            "name": "Equilibration 1",
                            "time_ns": 0.125,
                            "steps": 125000,
                            "timestep": 1.0,
                            "constraints": [],
                        }
                    ],
                },
            }
        ),
        encoding="utf-8",
    )
    (eq / "step1_equilibration.mdp").write_text(
        "define = -DPOSRES -DPOSRES_FC_BB=4184.0 -DPOSRES_FC_SC=2092.0 "
        "-DPOSRES_FC_LIPID=1046.0 -DPOSRES_FC_ION=4184.0\n"
        "integrator = md\nnsteps = 125000\ndt = 0.001\npcoupl = no\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    stage = meta["protocol"]["stages"][0]
    by_sel = {c["selection"]: c["force_constant"] for c in stage["constraints"]}
    assert by_sel["protein_backbone"] == 10.0
    assert by_sel["protein_sidechain"] == 5.0
    assert by_sel["lipid_head"] == 2.5
    assert by_sel["lipid_tail"] == 2.5
    assert by_sel["ions"] == 10.0


def test_infer_openmm_positional_restraints_from_inp(tmp_path: Path) -> None:
    eq = tmp_path / "openmm_nvt"
    eq.mkdir()
    (eq / "step1_equilibration.inp").write_text(
        "rest = yes\nnstep = 125000\ndt = 0.001\n"
        "fc_bb = 4184.0\nfc_sc = 2092.0\nfc_lpos = 1046.0\npcouple = no\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    stage = meta["protocol"]["stages"][0]
    by_sel = {c["selection"]: c["force_constant"] for c in stage["constraints"]}
    assert by_sel["protein_backbone"] == 10.0
    assert by_sel["protein_sidechain"] == 5.0
    assert by_sel["lipid_head"] == 2.5


def test_infer_amber_positional_restraints_from_mdin(tmp_path: Path) -> None:
    eq = tmp_path / "amber_nvt"
    eq.mkdir()
    (eq / "step1_equilibration.mdin").write_text(
        "&cntrl\n  ntr=1, nstlim=125000, dt=0.001,\n/\n"
        "protein backbone\n10.0\nATOM 1\nEND\n"
        "protein sidechain\n5.0\nATOM 2\nEND\n"
        "lipid head\n2.5\nATOM 3\nEND\n"
        "lipid tail\n2.5\nATOM 4\nEND\n"
        "ions\n10.0\nATOM 5\nEND\nEND\n",
        encoding="utf-8",
    )
    meta = infer_equilibration_job_metadata(eq, heal=False)
    stage = meta["protocol"]["stages"][0]
    by_sel = {c["selection"]: c["force_constant"] for c in stage["constraints"]}
    assert by_sel["protein_backbone"] == 10.0
    assert by_sel["protein_sidechain"] == 5.0
    assert by_sel["ions"] == 10.0
