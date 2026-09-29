"""Execute one documentation example script exactly once.

Numbered scripts under ``tests/<module>_examples/`` are the public
"copy this and it works" contract. Pytest execs each file once after
changing into that folder so relative paths resolve.
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import sys
import time
from pathlib import Path

import pytest


def discover_example_scripts(examples_dir: Path, *patterns: str) -> list[Path]:
    """Return unique scripts matching any of *patterns*, sorted by name."""
    if not examples_dir.is_dir():
        return []
    scripts: list[Path] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(examples_dir.glob(pattern)):
            resolved = path.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            scripts.append(path)
    return scripts


def example_id(script: Path) -> str:
    return script.name


def parametrize_example_scripts(examples_dir: Path, *patterns: str):
    """Return ``pytest.mark.parametrize`` over discovered example scripts."""
    scripts = discover_example_scripts(examples_dir, *patterns)
    if not scripts:
        return pytest.mark.parametrize(
            "script",
            [
                pytest.param(
                    None,
                    marks=pytest.mark.skip(
                        reason=f"No example scripts in {examples_dir}"
                    ),
                )
            ],
            ids=["missing"],
        )
    return pytest.mark.parametrize(
        "script",
        scripts,
        ids=[example_id(path) for path in scripts],
    )


def run_example_script(script: Path | None) -> None:
    """importlib-exec *script* with cwd set to its parent directory."""
    if script is None:
        pytest.skip("No example script")

    script = Path(script).resolve()
    if not script.is_file():
        pytest.fail(f"Example script not found: {script}")

    spec = importlib.util.spec_from_file_location(f"gw_example_{script.stem}", script)
    if spec is None or spec.loader is None:
        pytest.fail(f"Could not load example: {script}")

    module = importlib.util.module_from_spec(spec)
    original_dir = os.getcwd()
    os.chdir(script.parent)
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is None:
            continue
        try:
            reconfigure(encoding="utf-8", errors="replace")
        except Exception:
            pass
    try:
        spec.loader.exec_module(module)
    except SystemExit as exc:
        code = exc.code
        if code in (0, None):
            pytest.skip("Example script exited 0")
        if isinstance(code, str) and code.startswith("Skip"):
            pytest.skip(code)
        pytest.fail(f"Example exited with {code}")
    finally:
        os.chdir(original_dir)


def remove_tree(path: Path, attempts: int = 3) -> None:
    """Best-effort rmtree for example output dirs (cloud-sync locks)."""
    if not path.exists():
        return
    for _ in range(attempts):
        try:
            shutil.rmtree(path)
            return
        except (PermissionError, OSError):
            time.sleep(1)
