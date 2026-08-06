"""Tests for energetic per-file stride helpers and OpenMM parse stride."""

from pathlib import Path

from gatewizard.utils.energy_stride import (
    apply_energy_stride_to_result,
    energy_keep_indices,
    lookup_file_map,
)
from gatewizard.utils.openmm_analysis import OpenMMLogAnalyzer, run_openmm_energetic_analysis


def test_lookup_file_map_basename():
    path = Path("/tmp/prod.log")
    assert lookup_file_map({"prod.log": 10}, path) == 10
    assert lookup_file_map({"PROD.LOG": 5}, path) == 5
    assert lookup_file_map(None, path) is None


def test_energy_keep_indices_uniform_stride():
    files = [Path("a.log"), Path("b.log")]
    ranges = {"a.log": (0, 10), "b.log": (10, 20)}
    # keys may be full path — also test with Path str
    ranges = {str(files[0]): (0, 10), str(files[1]): (10, 20)}
    keep = energy_keep_indices(20, files, ranges, {"a.log": 2, "b.log": 5})
    assert keep == [0, 2, 4, 6, 8, 10, 15]


def test_apply_energy_stride_to_result():
    result = {
        "x": list(range(10)),
        "series": [{"name": "T", "key": "temp", "unit": "K", "y": [float(i) for i in range(10)]}],
        "statistics": {
            "temp": {
                "mean": 4.5,
                "std": 1.0,
                "min": 0.0,
                "max": 9.0,
                "initial": 0.0,
                "final": 9.0,
            }
        },
    }
    files = [Path("prod.log")]
    ranges = {str(files[0]): (0, 10)}
    out = apply_energy_stride_to_result(result, files, ranges, {"prod.log": 3})
    assert out["x"] == [0, 3, 6, 9]
    assert out["series"][0]["y"] == [0.0, 3.0, 6.0, 9.0]
    assert out["n_points_after_stride"] == 4


def test_openmm_analyzer_honors_file_strides(tmp_path: Path):
    log = tmp_path / "state.log"
    lines = ['#"Step"\t"Temperature (K)"\n']
    for step in range(0, 100, 10):
        lines.append(f"{step}\t{300.0 + step * 0.01}\n")
    log.write_text("".join(lines), encoding="utf-8")

    full = OpenMMLogAnalyzer([log])
    assert len(full.data["step"]) == 10

    strided = OpenMMLogAnalyzer([log], file_strides={"state.log": 2})
    assert len(strided.data["step"]) == 5
    assert strided.data["step"] == [0.0, 20.0, 40.0, 60.0, 80.0]

    result = run_openmm_energetic_analysis(
        [str(log)],
        properties=["Temperature"],
        file_strides={"state.log": 5},
    )
    assert len(result["x"]) == 2


def test_list_openmm_energy_properties_header_only(tmp_path: Path):
    from gatewizard.utils.openmm_analysis import list_openmm_energy_properties

    log = tmp_path / "huge_header.log"
    # Header + one row is enough; do not require full-file parse.
    log.write_text(
        '#"Step"\t"Potential Energy (kJ/mole)"\t"Temperature (K)"\n'
        "0\t-1000.0\t300.0\n",
        encoding="utf-8",
    )
    props = list_openmm_energy_properties([str(log)])
    assert "Potential Energy" in props
    assert "Temperature" in props
