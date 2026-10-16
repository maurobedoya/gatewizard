"""One trajectory frame is written as a single-model PDB."""

from pathlib import Path

import MDAnalysis as mda

from gatewizard.utils.pdb_frame import write_current_frame_pdb


def _multimodel(path: Path, n_models: int) -> None:
    lines = []
    for model in range(n_models):
        x = float(model)
        lines.append(f"MODEL     {model + 1:4d}\n")
        lines.append(
            f"ATOM      1  N   MET A   1     {x:7.3f}   0.000   0.000  1.00  0.00           N\n"
        )
        lines.append("ENDMDL\n")
    lines.append("END\n")
    path.write_text("".join(lines), encoding="utf-8")


def test_write_current_frame_pdb_keeps_one_model(tmp_path: Path):
    src = tmp_path / "traj.pdb"
    _multimodel(src, 4)
    universe = mda.Universe(str(src))
    universe.trajectory[2]
    dest = tmp_path / "frame.pdb"
    written = write_current_frame_pdb(universe, dest)
    text = Path(written).read_text(encoding="utf-8")
    assert text.count("MODEL") <= 1
    assert text.count("\nATOM") == 1 or text.count("ATOM") == 1
    again = mda.Universe(written)
    assert len(again.trajectory) == 1
    assert abs(float(again.atoms.positions[0, 0]) - 2.0) < 1e-3
