"""Tests for MDAnalysis → cpptraj mask conversion used by Fix PBC."""

from gatewizard.utils.trajectory_tools import mda_selection_to_cpptraj_mask


def test_protein_mask():
    assert mda_selection_to_cpptraj_mask("protein").startswith(":ALA")


def test_resname_space_separated():
    assert mda_selection_to_cpptraj_mask("resname PA PC OL") == ":PA,PC,OL"


def test_resname_comma_separated():
    assert mda_selection_to_cpptraj_mask("resname PA,PC,OL") == ":PA,PC,OL"


def test_atom_name():
    assert mda_selection_to_cpptraj_mask("name P31") == "@P31"
    assert mda_selection_to_cpptraj_mask("name P31 P32") == "@P31,P32"


def test_already_amber():
    assert mda_selection_to_cpptraj_mask(":PA,PC,OL") == ":PA,PC,OL"
    assert mda_selection_to_cpptraj_mask("@P31") == "@P31"


def test_protein_or_membrane():
    mask = mda_selection_to_cpptraj_mask("protein or resname PA PC OL")
    assert mask.startswith("(")
    assert ":PA,PC,OL" in mask
    assert ":ALA" in mask


def test_bare_residue_list():
    assert mda_selection_to_cpptraj_mask("PA PC OL") == ":PA,PC,OL"
