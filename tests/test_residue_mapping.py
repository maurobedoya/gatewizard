"""Tests for GateWizard residue mapping file parsing / remapping."""

from pathlib import Path

import pytest

from gatewizard.utils.residue_mapping import (
    align_mapping_keys_to_pdb,
    compose_amber_to_original,
    parse_residue_mapping_file,
    parse_residue_mapping_text,
    parse_residue_mapping_text_by_chain,
    remap_residue_type_labels,
    remap_resids_to_original,
    should_apply_residue_mapping,
)

SAMPLE = """\
# GateWizard Residue Mapping File
# Generated when protein capping is enabled
# Format: ORIGINAL_RESNAME CHAIN ORIGINAL_ID FINAL_RESNAME FINAL_ID
# '-' indicates caps (ACE/NME) that did not exist in original structure
#
  - A      -    ACE   1
VAL A     21    VAL   2
ALA A     22    ALA   3
GLY A     23    GLY   4
  - A      -    NME   5
"""

PROTONATED_RENUM = """\
ACE A     1    ACE     1
VAL A     2    VAL    21
ALA A     3    ALA    22
GLY A     4    GLY    23
"""


def test_parse_skips_caps_and_comments():
    m = parse_residue_mapping_text(SAMPLE)
    assert m == {2: 21, 3: 22, 4: 23}
    assert 1 not in m
    assert 5 not in m


def test_parse_protonated_renum_final_then_original():
    m = parse_residue_mapping_text(
        PROTONATED_RENUM, path="protein_protonated_renum.txt"
    )
    assert m[2] == 21
    assert m[3] == 22
    assert m[1] == 1


def test_parse_pdb4amber_renum_original_then_final():
    text = "ALA A    10    ALA     1\nGLY A    20    GLY     2\nLIG A   900    LIG     3\n"
    m = parse_residue_mapping_text(text, path="prepared_renum.txt")
    assert m == {1: 10, 2: 20, 3: 900}


def test_parse_four_column_renum_blank_chain():
    text = "MET       1    MET     1\nLYS       2    LYS     2\nY01 A   301    Y01   523\n"
    by_chain = parse_residue_mapping_text_by_chain(text, path="9g9v_protonated_renum.txt")
    assert by_chain[(" ", 1)] == 1
    assert by_chain[(" ", 2)] == 2
    assert by_chain[("A", 301)] == 523 or by_chain[("A", 523)] == 301


def test_align_mapping_flips_swapped_columns():
    swapped = {("A", 10): 1, ("A", 20): 2}
    aligned = align_mapping_keys_to_pdb(swapped, [("A", 1), ("A", 2)])
    assert aligned == {("A", 1): 10, ("A", 2): 20}


def test_should_apply_gatewizard_not_charmm():
    m = parse_residue_mapping_text(SAMPLE)
    assert should_apply_residue_mapping([2, 3, 4], m) is True
    assert should_apply_residue_mapping([21, 22, 23, 100], m) is False
    assert should_apply_residue_mapping([], m) is False


def test_remap_resids():
    m = parse_residue_mapping_text(SAMPLE)
    assert remap_resids_to_original([2, 3, 4], m) == [21, 22, 23]
    # Charmm-like: left unchanged
    assert remap_resids_to_original([21, 22, 100], m) == [21, 22, 100]
    # force remaps topology finals even when mixed with originals
    assert remap_resids_to_original([2, 3], m, force=True) == [21, 22]
    # over-remapped full-length axis rebuilt from sorted finals
    assert remap_resids_to_original([40, 41, 42], m) == [21, 22, 23]


def test_remap_labels():
    m = parse_residue_mapping_text(SAMPLE)
    assert remap_residue_type_labels(
        ["VAL2", "ALA3", "GLY4"], [2, 3, 4], m
    ) == ["VAL21", "ALA22", "GLY23"]
    assert remap_residue_type_labels(
        ["VAL21", "ALA22"], [21, 22], m
    ) == ["VAL21", "ALA22"]


def test_parse_file(tmp_path: Path):
    p = tmp_path / "prot_gatewizard_residue_mapping.txt"
    p.write_text(SAMPLE, encoding="utf-8")
    assert parse_residue_mapping_file(p)[2] == 21


def test_parse_by_chain_keeps_caps_and_distinguishes_chains():
    two_chain = SAMPLE + "VAL B     21    VAL   2\n"
    m = parse_residue_mapping_text_by_chain(two_chain)
    assert m[("A", 1)] is None
    assert m[("A", 2)] == 21
    assert m[("A", 5)] is None
    assert m[("B", 2)] == 21


def test_compose_amber_to_original_through_cap_map():
    remum = {("A", 1): 1, ("A", 2): 2, ("A", 3): 3}
    cap = {("A", 1): None, ("A", 2): 200, ("A", 3): 205}
    composed = compose_amber_to_original(remum, cap)
    assert composed[("A", 1)] is None
    assert composed[("A", 2)] == 200
    assert composed[("A", 3)] == 205
