"""Tests for GateWizard residue mapping file parsing / remapping."""

from pathlib import Path

import pytest

from gatewizard.utils.residue_mapping import (
    parse_residue_mapping_text,
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
    from gatewizard.utils.residue_mapping import parse_residue_mapping_file

    assert parse_residue_mapping_file(p)[2] == 21
