#!/usr/bin/env python3
"""OpenMM equilibration inputs for the 7AAQ build (same files as the GUI job).

Closest API pieces:
  gatewizard/tools/equilibration.py          OpenMMEquilibrationManager
  tests/equilibration_examples/openmm_example_01.py   one stage
  tests/study_cases/prepare_systems.py       run_protocol (GUI membrane schedule)

Input:
  02_build_7aaq_protein_oriented_protonated/{system.prmtop,system.inpcrd,system.pdb,
                                   bilayer_7aaq_protein_oriented_protonated_lipid.pdb}

This writes inputs only. It does not run OpenMM.

Run:
  conda activate gatewizard
  python equilibrate_7aaq.py
"""

from __future__ import annotations

from pathlib import Path

from gatewizard.tools.equilibration import EquilibrationStage, OpenMMEquilibrationManager
from gatewizard.utils.equilibration_job_metadata import write_equilibration_job_metadata
from gatewizard.utils.equilibration_resources import write_equilibration_resources

WORK = Path(__file__).resolve().parent
BUILD_DIR = WORK / "02_build_7aaq_protein_oriented_protonated"
EQ_DIR = WORK / "03_equilibration_02_build_7aaq_protein_oriented_protonated"
TEMPERATURE = 303.15
PRODUCTION_NS = 100.0
SCHEME = "NPT"
# Desktop OpenMM start: 1 CPU, 1 GPU, compute target Auto.
# Auto leaves PLATFORM empty so OpenMM picks the platform when the job starts.
COMPUTE_DEFAULTS = {
    "cpu_cores": 1,
    "gpu_id": 0,
    "num_gpus": 1,
    "use_gpu": True,
    "compute_target": "auto",
}

EXPECTED = (
    "restraints",
    "bilayer_7aaq_protein_oriented_protonated_lipid.pdb",
    "equilibration_job.json",
    "equilibration_resources.json",
    "openmm_run.py",
    "run_equilibration.sh",
    "run_equilibration_cluster.sh",
    "step1_equilibration.inp",
    "step2_equilibration.inp",
    "step3_equilibration.inp",
    "step4_equilibration.inp",
    "step5_equilibration.inp",
    "step6_equilibration.inp",
    "step7_production.inp",
    "system.inpcrd",
    "system.pdb",
    "system.prmtop",
)


def stage_dicts() -> list[dict]:
    stages = OpenMMEquilibrationManager.get_default_stage_params(
        scheme_type=SCHEME,
        temperature=TEMPERATURE,
        include_production=True,
    )
    # The desktop protocol leads with this card. OpenMM writes the steps into
    # step1_equilibration.inp and does not emit a separate minimization file.
    out = [
        {
            "name": "Minimization",
            "stage_kind": "minimization",
            "ensemble": "NVT",
            "time_ns": 0.0,
            "steps": 0,
            "minimize_steps": 10000,
            "temperature": TEMPERATURE,
            "timestep": 1.0,
            "cpu_cores": COMPUTE_DEFAULTS["cpu_cores"],
            "gpu_id": COMPUTE_DEFAULTS["gpu_id"],
            "use_gpu": True,
            "num_gpus": COMPUTE_DEFAULTS["num_gpus"],
            "compute_target": COMPUTE_DEFAULTS["compute_target"],
        }
    ]
    for stage in stages:
        data = stage.to_dict() if isinstance(stage, EquilibrationStage) else dict(stage)
        if data.get("stage_kind") == "production" or data.get("name") == "Production":
            data["time_ns"] = PRODUCTION_NS
            data["steps"] = int(round(PRODUCTION_NS * 1_000_000 / float(data["timestep"])))
        data["cpu_cores"] = COMPUTE_DEFAULTS["cpu_cores"]
        data["gpu_id"] = COMPUTE_DEFAULTS["gpu_id"]
        data["use_gpu"] = True
        data["num_gpus"] = COMPUTE_DEFAULTS["num_gpus"]
        data["compute_target"] = COMPUTE_DEFAULTS["compute_target"]
        out.append(data)
    return out


def main() -> None:
    required = (
        BUILD_DIR / "system.prmtop",
        BUILD_DIR / "system.inpcrd",
        BUILD_DIR / "system.pdb",
        BUILD_DIR / "bilayer_7aaq_protein_oriented_protonated_lipid.pdb",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing build files:\n" + "\n".join(missing))

    stages = stage_dicts()
    manager = OpenMMEquilibrationManager(BUILD_DIR)
    result = manager.setup_openmm_equilibration(
        system_files={
            "prmtop": str(BUILD_DIR / "system.prmtop"),
            "inpcrd": str(BUILD_DIR / "system.inpcrd"),
            "pdb": str(BUILD_DIR / "system.pdb"),
            "bilayer_pdb": str(BUILD_DIR / "bilayer_7aaq_protein_oriented_protonated_lipid.pdb"),
        },
        stage_params_list=stages,
        output_name=str(EQ_DIR),
        scheme_type=SCHEME,
        add_com_restraint=False,
        add_rotation_restraint=False,
    )
    openmm_dir = Path(result["openmm_dir"])

    write_equilibration_job_metadata(
        openmm_dir,
        input_dir=str(BUILD_DIR),
        ensemble=SCHEME,
        protocol={"name": "Default Protocol", "stages": stages},
        engine="openmm",
    )
    write_equilibration_resources(
        openmm_dir,
        "openmm",
        stages,
        compute_defaults=COMPUTE_DEFAULTS,
    )

    print(f"Output: {openmm_dir}")
    for name in EXPECTED:
        path = openmm_dir / name
        print(f"  {'OK' if path.exists() else 'MISSING'}  {name}")
    print("This only generated inputs. Run OpenMM with:")
    print(f"  cd {openmm_dir} && bash run_equilibration.sh")


if __name__ == "__main__":
    main()
