"""Tests for NAMD OPC water helpers."""

import pytest

from gatewizard.tools.namd_water import (
    namd_water_model_config_block,
    needs_namd_flexible_water_tleap,
    tleap_flexible_water_lines,
)


class TestNamdWaterHelpers:
    def test_flexible_water_only_namd_opc(self):
        assert needs_namd_flexible_water_tleap("namd", "opc") is True
        assert needs_namd_flexible_water_tleap("namd", "tip3p") is False
        assert needs_namd_flexible_water_tleap("gromacs", "opc") is False
        assert needs_namd_flexible_water_tleap(None, "opc") is False

    def test_tleap_lines_opc_namd(self):
        lines = tleap_flexible_water_lines("namd", "opc")
        assert "FlexibleWater on" in lines
        assert tleap_flexible_water_lines("namd", "tip3p") == ""

    def test_namd_config_opc(self):
        block = namd_water_model_config_block("opc")
        assert "waterModel              tip4" in block
        assert "useSettle               on" in block
        assert "rigidTolerance          1.0e-8" in block

    def test_namd_config_tip3p_empty(self):
        block = namd_water_model_config_block("tip3p")
        assert "waterModel" not in block
        assert block.strip().startswith("#")
