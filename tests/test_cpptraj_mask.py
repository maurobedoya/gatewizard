"""Tests for MDAnalysis → cpptraj mask conversion used by Fix PBC."""

from gatewizard.utils.trajectory_tools import (
    build_cpptraj_autoimage_attempts,
    compact_autoimage_anchor,
    mda_selection_to_cpptraj_mask,
)


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


def test_compact_anchor_prefers_protein_half():
    mask = mda_selection_to_cpptraj_mask("protein or resname PA PC OL")
    anchor = compact_autoimage_anchor(mask)
    assert anchor is not None
    assert "ALA" in anchor
    assert "PA" not in anchor


def test_compact_anchor_omits_multi_lipid():
    assert compact_autoimage_anchor(":PA,PC,OL") is None
    assert compact_autoimage_anchor("@P31") is None


def test_membrane_autoimage_attempts_no_origin():
    mask = mda_selection_to_cpptraj_mask("resname PA PC OL")
    attempts = build_cpptraj_autoimage_attempts(mask, membrane_like=True)
    assert attempts
    joined = "\n".join(attempts)
    assert "origin" not in joined.lower()
    assert "mode byvec" in joined
    assert "moveanchor" in joined
    assert f"center {mask} mass" in joined


def test_protein_membrane_autoimage_uses_protein_anchor():
    mask = mda_selection_to_cpptraj_mask("protein or resname PA PC OL")
    attempts = build_cpptraj_autoimage_attempts(mask, membrane_like=True)
    assert attempts[0].startswith("autoimage anchor")
    assert "mode byvec moveanchor" in attempts[0]
    assert "ALA" in attempts[0].splitlines()[0]
    assert "origin" not in "\n".join(attempts).lower()
