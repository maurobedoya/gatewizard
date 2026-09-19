"""Tests for ion recognition and tleap neutralization ion selection."""

from __future__ import annotations

from pathlib import Path

import pytest

from gatewizard.core.builder import Builder
from gatewizard.utils.ions import (
    DEFAULT_ANION,
    DEFAULT_CATION,
    ion_mda_selection,
    is_ion_resname,
    resolve_amber_anion,
    resolve_amber_cation,
    tleap_neutralization_lines,
)


def test_is_ion_resname_amber_mixed_case():
    """Amber PDB writes Na+/K+/Cl-; recognition must be case-insensitive."""
    for name in ("Na+", "NA+", "na+", "K+", "k+", "Cl-", "CL-", "Cl", "SOD", "POT", "CLA"):
        assert is_ion_resname(name), name
    assert not is_ion_resname("HOH")
    assert not is_ion_resname("ALA")
    assert not is_ion_resname("POPC")
    assert not is_ion_resname("ETA")


def test_ion_selection_does_not_include_protein(tmp_path: Path):
    """PDB potassium (resname K) must not pull in protein atoms via ``ion``."""
    pytest.importorskip("MDAnalysis")
    import MDAnalysis as mda

    pdb = tmp_path / "prot_k.pdb"
    pdb.write_text(
        "ATOM      1  N   MET A   1       0.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      2  CA  MET A   1       1.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      3  C   MET A   1       2.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    4  K     K A 302      10.000   0.000   0.000  1.00  0.00           K\n"
        "END\n",
        encoding="utf-8",
    )
    u = mda.Universe(str(pdb))
    by_resname = u.select_atoms(ion_mda_selection())
    by_mask = u.atoms[[is_ion_resname(str(rn)) for rn in u.atoms.resnames]]
    assert len(by_resname) == 1
    assert list(by_resname.resnames) == ["K"]
    assert len(by_mask) == 1
    assert list(by_mask.resnames) == ["K"]
    assert "MET" not in set(by_resname.resnames)
    assert "MET" not in set(by_mask.resnames)


def test_resolve_amber_ions_defaults_and_aliases():
    assert resolve_amber_cation(None) == DEFAULT_CATION
    assert resolve_amber_cation("") == DEFAULT_CATION
    assert resolve_amber_cation("na+") == "Na+"
    assert resolve_amber_cation("K+") == "K+"
    assert resolve_amber_anion(None) == DEFAULT_ANION
    assert resolve_amber_anion("cl-") == "Cl-"
    assert resolve_amber_anion("Br-") == "Br-"


def test_tleap_neutralization_uses_selected_ions_not_hardcoded_na():
    lines = tleap_neutralization_lines("K+", "Cl-")
    assert "addIonsRand system K+ 0" in lines
    assert "addIonsRand system Cl- 0" in lines
    assert "Na+" not in lines

    lines_na = tleap_neutralization_lines("Na+", "Cl-")
    assert "addIonsRand system Na+ 0" in lines_na
    assert "addIonsRand system Cl- 0" in lines_na


def test_create_tleap_input_uses_builder_cation_anion(tmp_path: Path):
    pdb = tmp_path / "sys.pdb"
    pdb.write_text("END\n", encoding="utf-8")
    builder = Builder()
    content = builder._create_tleap_input(
        str(pdb),
        {**builder.config, "cation": "K+", "anion": "Cl-"},
    )
    assert "addIonsRand system K+ 0" in content
    assert "addIonsRand system Cl- 0" in content
    assert "addIonsRand system Na+ 0" not in content

    content_na = builder._create_tleap_input(
        str(pdb),
        {**builder.config, "cation": "Na+", "anion": "Br-"},
    )
    assert "addIonsRand system Na+ 0" in content_na
    assert "addIonsRand system Br- 0" in content_na


def test_generate_preparation_script_embeds_neutralization_ions(tmp_path: Path):
    builder = Builder()
    success, message, job_dir = builder.generate_preparation_inputs(
        pdb_file=None,
        working_dir=str(tmp_path),
        upper_lipids=["POPC"],
        lower_lipids=["POPC"],
        lipid_ratios="1//1",
        distxy_fix=100,
        water_model="tip3p",
        protein_ff="ff19SB",
        lipid_ff="lipid21",
        cation="K+",
        anion="Cl-",
        add_salt=True,
        salt_concentration=0.15,
        output_folder_name="02_build_ion_test",
    )
    assert success, message
    script = (Path(job_dir) / "run_preparation.sh").read_text(encoding="utf-8")
    assert "addIonsRand system K+ 0" in script
    assert "addIonsRand system Cl- 0" in script
    assert "addIonsRand system Na+ 0" not in script


def test_structure_manager_groups_amber_ions(tmp_path: Path):
    pytest.importorskip("MDAnalysis")
    from gatewizard.core.structure_manager import StructureManager

    pdb = tmp_path / "ions.pdb"
    pdb.write_text(
        "HETATM    1 Na+  Na+ A   1       0.000   0.000   0.000  1.00  0.00          Na\n"
        "HETATM    2  K+   K+ A   2       1.000   0.000   0.000  1.00  0.00           K\n"
        "HETATM    3 Cl-  Cl- A   3       2.000   0.000   0.000  1.00  0.00          Cl\n"
        "ATOM      4  N   ALA A   4       3.000   0.000   0.000  1.00  0.00           N\n"
        "ATOM      5  CA  ALA A   4       4.000   0.000   0.000  1.00  0.00           C\n"
        "END\n",
        encoding="utf-8",
    )
    sm = StructureManager()
    sm.load_structure(str(pdb))
    sels = sm.auto_detect_molecules()
    names = [s.name for s in sels]
    assert "Ions" in names
    assert "Na+" not in names
    assert "K+" not in names
    assert "Cl-" not in names
    ions = next(s for s in sels if s.name == "Ions")
    assert len(ions.atom_indices) == 3
