"""
NAMD-specific water model helpers.

OPC (and other 4-site Amber waters) need FlexibleWater in LEaP when the target
MD engine is NAMD, plus ``waterModel tip4`` in NAMD config files.
"""

from __future__ import annotations

from typing import Optional

# 4-site waters that use NAMD ``waterModel tip4`` + FlexibleWater prmtop
NAMD_FOUR_SITE_WATER_MODELS = frozenset({"opc", "tip4pd", "tip4pew"})


def normalize_water_model(water_model: Optional[str]) -> str:
    return (water_model or "tip3p").lower().strip()


def normalize_md_engine(md_engine: Optional[str]) -> str:
    return (md_engine or "").lower().strip()


def needs_namd_flexible_water_tleap(
    md_engine: Optional[str], water_model: Optional[str]
) -> bool:
    """True when tleap should use ``set default FlexibleWater on``."""
    return normalize_md_engine(md_engine) == "namd" and normalize_water_model(
        water_model
    ) in NAMD_FOUR_SITE_WATER_MODELS


def tleap_flexible_water_lines(
    md_engine: Optional[str], water_model: Optional[str]
) -> str:
    """LEaP directives inserted after sourcing the water leaprc."""
    if not needs_namd_flexible_water_tleap(md_engine, water_model):
        return ""
    return (
        "# NAMD 4-site water: write H-O-H angle into prmtop (required for waterModel tip4)\n"
        "set default FlexibleWater on\n"
    )


def namd_water_model_config_block(water_model: Optional[str]) -> str:
    """
    Extra NAMD keywords inserted after ``rigidBonds all`` in equilibration templates.

    Empty (comment only) for TIP3P and other 3-site models.
    """
    wm = normalize_water_model(water_model)
    if wm in NAMD_FOUR_SITE_WATER_MODELS:
        return (
            "waterModel              tip4\n"
            "useSettle               on\n"
            "rigidTolerance          1.0e-8\n"
        )
    return "# TIP3P / default Amber water (no extra NAMD water keywords)\n"


def read_water_model_from_builder_status(working_dir) -> Optional[str]:
    """Read ``water_model`` from a builder job ``status.json`` if present."""
    from pathlib import Path
    import json

    status_file = Path(working_dir) / "status.json"
    if not status_file.is_file():
        return None
    try:
        data = json.loads(status_file.read_text(encoding="utf-8"))
        wm = data.get("config", {}).get("water_model")
        return normalize_water_model(wm) if wm else None
    except (OSError, json.JSONDecodeError, AttributeError):
        return None
