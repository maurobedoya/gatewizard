"""
Optional dependency handling utilities for Gatewizard.

This module provides utilities for handling optional dependencies
in a consistent way across the application.
"""

import importlib
import json
import os
import re
import shutil
import subprocess
from importlib import metadata
from pathlib import Path
from typing import Optional, Dict, Any, Tuple, List
from functools import wraps

from gatewizard.utils.logger import get_logger

logger = get_logger(__name__)


class OptionalDependencyError(Exception):
    """Exception raised when an optional dependency is missing."""

    pass


# import_name, distribution_name, required, description
DependencySpec = Tuple[str, str, bool, str]

# candidates, description
ExternalToolSpec = Tuple[Tuple[str, ...], str]

DEPENDENCY_REGISTRY: Dict[str, DependencySpec] = {
    "gatewizard": (
        "gatewizard",
        "gatewizard",
        True,
        "GateWizard core package",
    ),
    "numpy": ("numpy", "numpy", True, "Numerical arrays and linear algebra"),
    "matplotlib": (
        "matplotlib",
        "matplotlib",
        True,
        "Plotting and visualization",
    ),
    "requests": ("requests", "requests", True, "HTTP client library"),
    "MDAnalysis": (
        "MDAnalysis",
        "MDAnalysis",
        True,
        "Molecular dynamics trajectory analysis",
    ),
    "lipyphilic": (
        "lipyphilic",
        "lipyphilic",
        True,
        "Lipid bilayer analysis",
    ),
    "propka": ("propka", "propka", True, "Protonation state prediction"),
    "Pillow": ("PIL", "Pillow", True, "Image handling for ligand previews"),
    "rdkit": ("rdkit", "rdkit", True, "Cheminformatics and ligand parametrization"),
    "psique": ("psique", "psique", True, "Structure visualization support"),
    "parmed": (
        "parmed",
        "parmed",
        False,
        "ParmEd for NAMD/Amber conversion (install: pip install gatewizard[md])",
    ),
    "openmm": (
        "openmm",
        "openmm",
        False,
        "OpenMM molecular dynamics engine (install: pip install gatewizard[md])",
    ),
    "mempro": (
        "mempro",
        "mempro",
        False,
        "Membrane protein orientation (install separately: pip install git+https://github.com/pstansfeld/MemPrO.git)",
    ),
}

# Pip install groups for packages not in the core pyproject dependencies list.
INSTALL_GROUPS: Dict[str, str] = {
    "parmed": "md",
    "openmm": "md",
    "mempro": "orientation",
}

EXTERNAL_TOOL_REGISTRY: Dict[str, ExternalToolSpec] = {
    "mempro": (("mempro",), "Membrane protein orientation CLI (MemPrO)"),
    "packmol": (("packmol",), "Molecule packing (AmberTools)"),
    "packmol-memgen": (
        ("packmol-memgen",),
        "Membrane system builder (AmberTools)",
    ),
    "ambertools": (
        ("antechamber", "tleap", "pdb4amber"),
        "AmberTools suite (conda-forge; provides tleap, antechamber, packmol-memgen)",
    ),
    "namd": (("namd3", "namd2"), "NAMD molecular dynamics"),
    "gromacs": (("gmx",), "GROMACS molecular dynamics"),
}

_VERSION_PROBE_ARGS: Dict[str, List[List[str]]] = {
    "namd": [["+version"], ["-version"], ["--version"]],
    "gromacs": [["--version"]],
    "packmol": [["--version"]],
    "packmol-memgen": [["--version"], ["-h"], ["--help"]],
    "mempro": [["--version"]],
    "tleap": [["-h"]],
    "antechamber": [["-L"], ["-h"]],
    "pdb4amber": [["-h"]],
}

_VERSION_PATTERNS: Dict[str, List[re.Pattern[str]]] = {
    "packmol": [re.compile(r"Version\s+([\d.]+)", re.I)],
    "packmol-memgen": [
        re.compile(r"packmol[- ]memgen[^\d]*([\d.]+)", re.I),
        re.compile(r"Version\s+([\d.]+)", re.I),
    ],
    "gromacs": [
        re.compile(r"GROMACS.*?,\s*([\d.]+)", re.I),
        re.compile(r"GROMACS version[:\s]+([\d.]+)", re.I),
    ],
    "namd": [
        re.compile(r"NAMD\s+(?:version\s+)?([0-9]+(?:\.[0-9]+)+)", re.I),
        re.compile(r"Info:\s*NAMD\s+([0-9]+(?:\.[0-9]+)+)", re.I),
    ],
    "mempro": [
        re.compile(r"mempro[^\d]*([0-9]+\.[0-9]+\.[0-9]+)", re.I),
        re.compile(r"version\s+([0-9]+\.[0-9]+\.[0-9]+)", re.I),
    ],
    "tleap": [
        re.compile(r"AmberTools\s+([\d.]+)", re.I),
        re.compile(r"Leap:\s*.*Release\s+([\d.]+)", re.I),
    ],
    "antechamber": [
        re.compile(r"AmberTools\s+([\d.]+)", re.I),
        re.compile(r"Antechamber\s+([0-9.]+)", re.I),
    ],
    "ambertools": [
        re.compile(r"AmberTools\s+([\d.]+)", re.I),
        re.compile(r"^([\d.]+)\s*$"),
    ],
}

_NOISE_LINE_PATTERNS = [
    re.compile(r"^WARNING:", re.I),
    re.compile(r"^INFO:", re.I),
    re.compile(r"^DEBUG:", re.I),
    re.compile(r"^ERROR:", re.I),
    re.compile(r"Charm\+\+", re.I),
    re.compile(r"Converse\+\+", re.I),
    re.compile(r"^-I:\s", re.I),
    re.compile(r"invalid option", re.I),
    re.compile(r"teLeap:", re.I),
    re.compile(r"PMEMD not found", re.I),
    re.compile(r"Adding .+ to search path", re.I),
    re.compile(r"No provisioning arguments", re.I),
    re.compile(r"Falling back to cpu", re.I),
    re.compile(r"^Usage:", re.I),
]

_TOOLS_WITHOUT_GENERIC_FALLBACK = frozenset({"namd", "ambertools", "packmol-memgen"})

_CANONICAL_TOOL_NAMES = {
    "namd2": "namd",
    "namd3": "namd",
    "gmx": "gromacs",
}


def _canonical_tool_name(tool_name: str) -> str:
    return _CANONICAL_TOOL_NAMES.get(tool_name.lower(), tool_name.lower())


def require_optional_dependency(package_name: str, install_command: str = None):
    """
    Decorator to check for optional dependencies before function execution.

    Args:
        package_name: Name of the package to check
        install_command: Command to install the package (if different from pip install <package>)

    Raises:
        OptionalDependencyError: If the dependency is not available
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if not is_package_available(package_name):
                install_cmd = install_command or f"pip install {package_name}"
                raise OptionalDependencyError(
                    f"{package_name} is required for this operation. "
                    f"Install with: {install_cmd}"
                )
            return func(*args, **kwargs)

        return wrapper

    return decorator


def is_package_available(package_name: str) -> bool:
    """
    Check if a package is available for import.

    Args:
        package_name: Name of the package to check

    Returns:
        True if package is available, False otherwise
    """
    try:
        importlib.import_module(package_name)
        return True
    except ImportError:
        return False


def safe_import(package_name: str, alternative_name: str = None) -> Optional[Any]:
    """
    Safely import a package, returning None if not available.

    Args:
        package_name: Name of the package to import
        alternative_name: Alternative import name (e.g., 'parmed' for package, 'pmd' for import)

    Returns:
        The imported module or None if not available
    """
    try:
        module = importlib.import_module(package_name)
        logger.debug(f"Successfully imported {package_name}")
        return module
    except ImportError as e:
        logger.debug(f"Optional dependency {package_name} not available: {e}")
        return None


def get_package_version(
    import_name: str,
    *,
    distribution_name: Optional[str] = None,
) -> Optional[str]:
    """
    Resolve an installed package version.

    Tries importlib.metadata first, then module.__version__.
    """
    distribution = distribution_name or import_name
    try:
        return metadata.version(distribution)
    except metadata.PackageNotFoundError:
        pass

    module = safe_import(import_name)
    if module is not None:
        version = getattr(module, "__version__", None)
        if version is not None:
            return str(version)
    return None


def resolve_executable(candidates: Tuple[str, ...]) -> Optional[str]:
    """Return the first executable found on PATH from a list of candidate names."""
    for name in candidates:
        path = shutil.which(name)
        if path:
            return path
    return None


def _is_noise_line(line: str) -> bool:
    cleaned = line.strip()
    if not cleaned:
        return True
    return any(pattern.search(cleaned) for pattern in _NOISE_LINE_PATTERNS)


def _is_plausible_version(version: Optional[str]) -> bool:
    if not version:
        return False
    cleaned = version.strip()
    if len(cleaned) > 80:
        return False
    lower = cleaned.lower()
    if any(
        token in lower
        for token in ("invalid", "error", "option", "usage", "teLeap", "not found")
    ):
        return False
    if "/" in cleaned or "\\" in cleaned:
        return False
    if re.fullmatch(r"\d", cleaned):
        return False
    return bool(re.search(r"\d", cleaned))


def parse_tool_version(text: str, tool_name: str) -> Optional[str]:
    """Extract a concise version string from command output."""
    if not text:
        return None

    tool = _canonical_tool_name(tool_name)
    if tool == "namd":
        for line in text.splitlines():
            if _is_noise_line(line):
                continue
            for pattern in _VERSION_PATTERNS.get("namd", []):
                match = pattern.search(line)
                if match and _is_plausible_version(match.group(1)):
                    return match.group(1)
        return None

    for pattern in _VERSION_PATTERNS.get(tool, []):
        match = pattern.search(text)
        if match and _is_plausible_version(match.group(1)):
            return match.group(1)

    if tool in _TOOLS_WITHOUT_GENERIC_FALLBACK:
        return None

    for line in text.splitlines():
        if _is_noise_line(line):
            continue
        for pattern in _VERSION_PATTERNS.get(tool, []):
            match = pattern.search(line)
            if match and _is_plausible_version(match.group(1)):
                return match.group(1)
        generic = re.search(r"\b([0-9]+\.[0-9]+(?:\.[0-9]+)?)\b", line)
        if generic and not re.search(r"^\d{4}-\d{2}-\d{2}", line):
            return generic.group(1)

    return None


def _extract_version_line(text: str, tool_name: str = "") -> Optional[str]:
    return parse_tool_version(text, tool_name)


def probe_executable_version(executable: str, tool_name: str) -> Optional[str]:
    """Run lightweight version/help probes for an external executable."""
    canonical = _canonical_tool_name(tool_name)
    for args in _VERSION_PROBE_ARGS.get(canonical, _VERSION_PROBE_ARGS.get(tool_name, [["--version"]])):
        try:
            proc = subprocess.run(
                [executable] + args,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=4,
                check=False,
            )
            text = (proc.stdout or "").strip()
            if text:
                parsed = _extract_version_line(text, canonical)
                if _is_plausible_version(parsed):
                    return parsed
        except Exception:
            continue
    return None


def _version_from_install_path(executable: str, tool_name: str) -> Optional[str]:
    canonical = _canonical_tool_name(tool_name)
    if canonical == "namd":
        match = re.search(r"NAMD[_-]?(\d+(?:\.\d+)+)", executable, re.I)
        if match and _is_plausible_version(match.group(1)):
            return match.group(1)
    return None


def _ambertools_distribution_version() -> Optional[str]:
    for distribution in ("ambertools", "AmberTools"):
        try:
            return metadata.version(distribution)
        except metadata.PackageNotFoundError:
            continue
    return None


def _ambertools_version_from_conda_meta() -> Optional[str]:
    prefix = os.environ.get("CONDA_PREFIX")
    if not prefix:
        return None
    meta_dir = Path(prefix) / "conda-meta"
    if not meta_dir.is_dir():
        return None
    for path in sorted(meta_dir.glob("ambertools-*.json")):
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        version = data.get("version")
        if _is_plausible_version(version):
            return str(version)
    return None


def _ambertools_version_from_env() -> Optional[str]:
    roots = []
    for env_var in ("AMBERHOME", "CONDA_PREFIX"):
        value = os.environ.get(env_var)
        if value:
            roots.append(Path(value))

    for base in roots:
        for rel in (
            "VERSION",
            "share/ambertools/VERSION",
            "dat/leap/cmd/leap_version",
        ):
            path = base / rel
            if not path.is_file():
                continue
            text = path.read_text(encoding="utf-8", errors="ignore")
            parsed = parse_tool_version(text, "ambertools")
            if _is_plausible_version(parsed):
                return parsed
    return None


def _probe_ambertools_cli_version() -> Optional[str]:
    for candidate in ("antechamber", "tleap", "pdb4amber"):
        executable = resolve_executable((candidate,))
        if not executable:
            continue
        version = probe_executable_version(executable, candidate)
        if _is_plausible_version(version):
            return version
    return None


def _resolve_ambertools_version() -> Optional[str]:
    for resolver in (
        _ambertools_distribution_version,
        _ambertools_version_from_conda_meta,
        _ambertools_version_from_env,
        _probe_ambertools_cli_version,
    ):
        version = resolver()
        if _is_plausible_version(version):
            return version
    return None


def _external_tool_entry(name: str, spec: ExternalToolSpec) -> Dict[str, Any]:
    candidates, description = spec
    path = resolve_executable(candidates)
    version = None

    if name == "mempro":
        version = get_package_version("mempro", distribution_name="mempro")
        if not version and path:
            version = probe_executable_version(path, "mempro")
    elif name == "ambertools":
        version = _resolve_ambertools_version()
    elif name == "packmol-memgen":
        if path:
            version = probe_executable_version(path, "packmol-memgen")
        if not _is_plausible_version(version):
            version = _resolve_ambertools_version()
    elif name == "namd":
        if path:
            version = probe_executable_version(path, "namd")
            if not _is_plausible_version(version):
                version = _version_from_install_path(path, "namd")
    elif path:
        version = probe_executable_version(path, candidates[0])

    return {
        "name": name,
        "path": path,
        "version": version,
        "available": path is not None or (name == "ambertools" and version is not None),
        "description": description,
    }


def get_external_tool_versions() -> List[Dict[str, Any]]:
    """Return version information for external CLI tools used by GateWizard."""
    return [_external_tool_entry(name, spec) for name, spec in EXTERNAL_TOOL_REGISTRY.items()]


def _dependency_entry(name: str, spec: DependencySpec) -> Dict[str, Any]:
    import_name, distribution_name, required, description = spec
    available = is_package_available(import_name)
    version = get_package_version(import_name, distribution_name=distribution_name)
    if available and version is None and import_name != distribution_name:
        version = get_package_version(distribution_name, distribution_name=distribution_name)

    return {
        "available": available or version is not None,
        "required": required,
        "install_group": INSTALL_GROUPS.get(name, "core"),
        "description": description,
        "version": version,
    }


def get_dependency_versions(
    *,
    include_optional: bool = True,
    include_platform: bool = False,
    include_external_tools: bool = False,
) -> Dict[str, Any]:
    """
    Report installed versions for GateWizard dependencies.

    Returns:
        Dictionary with ``dependencies`` mapping package names to status info,
        optionally ``platform`` metadata, and ``executables`` for CLI tools.
    """
    dependencies = {}
    for name, spec in DEPENDENCY_REGISTRY.items():
        _import_name, _distribution, required, _description = spec
        if not include_optional and not required:
            continue
        dependencies[name] = _dependency_entry(name, spec)

    result: Dict[str, Any] = {"dependencies": dependencies}
    if include_platform:
        from gatewizard.utils.helpers import get_platform_info

        result["platform"] = get_platform_info()
    if include_external_tools:
        result["executables"] = get_external_tool_versions()
    return result


def get_optional_dependencies_status() -> Dict[str, Dict[str, Any]]:
    """
    Get the status of all optional dependencies.

    Returns:
        Dictionary mapping package names to availability status and version.
    """
    status = {}
    for name, spec in DEPENDENCY_REGISTRY.items():
        _import_name, _distribution, required, _description = spec
        if required:
            continue
        status[name] = _dependency_entry(name, spec)
    return status


def check_and_warn_missing_dependencies():
    """Log warnings for missing optional dependencies."""
    status = get_optional_dependencies_status()

    missing = [pkg for pkg, info in status.items() if not info["available"]]

    if missing:
        logger.warning(f"Missing optional dependencies: {', '.join(missing)}")
        logger.info("Some features may be limited. Install missing packages if needed.")
    else:
        logger.info("All optional dependencies are available")
