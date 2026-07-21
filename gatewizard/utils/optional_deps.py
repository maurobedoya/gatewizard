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
        "Membrane protein orientation (install: pip install -r requirements-orientation.txt)",
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
    # NAMD puts its version banner on "Info:" lines — do not treat those as noise.
    if re.search(r"^Info:\s*NAMD\s+[0-9]", cleaned, re.I):
        return False
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


def _probe_binary_output(executable: str, engine: str) -> str:
    """Return combined stdout from a short version probe (empty if all probes fail)."""
    args_map = {
        "namd": [["-version"], ["+version"], ["--version"]],
        "gromacs": [["--version"]],
        "openmm": [["--version"]],
    }
    timeout = 12 if engine == "namd" else 6
    for args in args_map.get(engine, [["--version"]]):
        try:
            proc = subprocess.run(
                [executable, *args],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                timeout=timeout,
                check=False,
            )
            text = (proc.stdout or "").strip()
            if text:
                return text
        except Exception:
            continue
    return ""


def _probe_binary_version(executable: str, engine: str) -> Optional[str]:
    """Run a short version probe and parse a concise version string."""
    text = _probe_binary_output(executable, engine)
    if text:
        parsed = parse_tool_version(text, engine)
        if parsed:
            return parsed
    if engine == "namd":
        return _version_from_install_path(executable, "namd")
    return None


def parse_engine_variant(
    text: str, engine: str, executable: str = ""
) -> Optional[str]:
    """Classify an MD engine install as CPU / CUDA / OpenCL / … when possible.

    Returns a short label suitable for UI pickers (``CUDA``, ``CPU``, ``OpenCL``,
    ``SYCL``, ``HIP``, ``Metal``), or ``None`` if unknown.
    """
    engine = (engine or "").strip().lower()
    blob = f"{text or ''}\n{executable or ''}"
    lower = blob.lower()

    if engine == "gromacs":
        # Official ``gmx --version`` line, e.g. "GPU support:             CUDA"
        m = re.search(r"GPU support:\s*([^\n\r]+)", text or "", re.I)
        if m:
            support = m.group(1).strip().lower()
            if not support or support in ("disabled", "no", "none", "false", "off"):
                return "CPU"
            if "cuda" in support:
                return "CUDA"
            if "opencl" in support:
                return "OpenCL"
            if "sycl" in support:
                return "SYCL"
            if "hip" in support or "rocm" in support:
                return "HIP"
            if "metal" in support:
                return "Metal"
            return support.split()[0].upper() if support else "CPU"
        # Conda build strings / paths: nompi_cuda, cuda12, …
        if re.search(r"nompi_cuda|cuda\d+|_cuda(?:_|\b)|/cuda", lower):
            return "CUDA"
        if "cuda" in lower and "disabled" not in lower:
            return "CUDA"
        return "CPU"

    if engine == "namd":
        # Install folder names: NAMD_3.0_Linux-x86_64-multicore-CUDA
        if re.search(r"\bcuda\b", lower) or re.search(r"multicore-cuda", lower):
            return "CUDA"
        if re.search(r"\bopencl\b", lower):
            return "OpenCL"
        if "multicore" in lower or "mpi" in lower:
            return "CPU"
        # Banner rarely states CUDA; path is the reliable signal.
        return "CPU" if executable else None

    if engine == "openmm":
        try:
            import openmm  # type: ignore

            names = [
                openmm.Platform.getPlatform(i).getName()
                for i in range(openmm.Platform.getNumPlatforms())
            ]
        except Exception:
            return None
        for preferred in ("CUDA", "OpenCL", "Metal", "HIP"):
            if preferred in names:
                return preferred
        if "CPU" in names:
            return "CPU"
        return None

    return None


def _probe_engine_variant(executable: str, engine: str, version_text: str = "") -> Optional[str]:
    text = version_text or _probe_binary_output(executable, engine)
    return parse_engine_variant(text, engine, executable)


def _format_engine_label(
    name: str,
    version: Optional[str],
    variant: Optional[str],
    source: str,
    executable: str,
) -> str:
    bits = [name]
    if version:
        bits[0] = f"{name} {version}"
    if variant:
        bits.append(f"· {variant}")
    bits.append(f"({source})")
    return " ".join(bits) + f" — {executable}"


def _discover_gmxrc_near(gmx_path: str) -> Optional[str]:
    """Return GMXRC next to a GROMACS install, if present."""
    p = Path(gmx_path)
    candidates = [
        p.parent / "GMXRC",
        p.parent.parent / "bin" / "GMXRC",
        p.parent.parent / "GMXRC",
    ]
    for c in candidates:
        if c.is_file():
            return str(c.resolve())
    return None


def _scan_gromacs_prefix(prefix: Path) -> List[Tuple[str, Optional[str]]]:
    """Return ``[(gmx_path, gmxrc_path|None), ...]`` under a GROMACS prefix."""
    found: List[Tuple[str, Optional[str]]] = []
    gmxrc = prefix / "bin" / "GMXRC"
    if not gmxrc.is_file():
        gmxrc = prefix / "GMXRC"
    gmx = prefix / "bin" / "gmx"
    if not gmx.is_file():
        gmx = prefix / "bin" / "gmx_mpi"
    if gmx.is_file() and os.access(gmx, os.X_OK):
        found.append(
            (str(gmx.resolve()), str(gmxrc.resolve()) if gmxrc.is_file() else None)
        )
    return found


def list_md_engine_candidates(engine: str) -> List[Dict[str, Any]]:
    """Discover installed MD engine binaries / packages for a picker UI.

    Engines:
      - ``gromacs``: conda ``gmx``, PATH hits, common prefixes (+ GMXRC when found)
      - ``namd``: ``namd3`` / ``namd2`` on PATH and common install dirs
      - ``openmm``: current Python OpenMM import (platforms remain separate)

    Returns a list of dicts with keys:
      ``id``, ``label``, ``executable``, ``version``, ``variant`` (CPU/CUDA/…),
      ``source``, ``gmxrc`` (optional).
    """
    engine = (engine or "").strip().lower()
    results: List[Dict[str, Any]] = []

    if engine == "openmm":
        version = get_package_version("openmm")
        py = shutil.which("python3") or shutil.which("python") or "python"
        variant = parse_engine_variant("", "openmm", py)
        label = (
            f"OpenMM {version}"
            if version
            else "OpenMM (not importable in current Python)"
        )
        if version and variant:
            label = f"OpenMM {version} · {variant} (current Python)"
        elif version:
            label = f"OpenMM {version} (current Python)"
        results.append(
            {
                "id": "openmm-current",
                "label": label,
                "executable": py,
                "version": version,
                "variant": variant,
                "source": "python",
                "available": version is not None,
            }
        )
        return results

    if engine == "gromacs":
        candidates: List[Tuple[str, Optional[str], str]] = []
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            conda_gmx = Path(conda_prefix) / "bin" / "gmx"
            if conda_gmx.is_file():
                candidates.append(
                    (str(conda_gmx.resolve()), None, "conda")
                )
        for name in ("gmx", "gmx_mpi"):
            which = shutil.which(name)
            if which:
                candidates.append((which, _discover_gmxrc_near(which), "path"))

        search_roots = [
            Path("/usr/local/gromacs"),
            Path("/opt/gromacs"),
            Path.home() / "gromacs",
            Path.home() / "local" / "gromacs",
        ]
        # Versioned prefixes e.g. /usr/local/gromacs-2024.2
        for parent in (Path("/usr/local"), Path("/opt"), Path.home()):
            if not parent.is_dir():
                continue
            try:
                for child in parent.iterdir():
                    if child.is_dir() and "gromacs" in child.name.lower():
                        search_roots.append(child)
            except OSError:
                pass

        for root in search_roots:
            if not root.is_dir():
                continue
            for gmx_path, gmxrc in _scan_gromacs_prefix(root):
                candidates.append((gmx_path, gmxrc, "gmxrc" if gmxrc else "prefix"))

        seen: set[str] = set()
        for exe, gmxrc, source in candidates:
            key = str(Path(exe).resolve()) if os.path.isfile(exe) else exe
            if key in seen:
                continue
            if not os.path.isfile(exe):
                continue
            seen.add(key)
            raw = _probe_binary_output(exe, "gromacs")
            version = parse_tool_version(raw, "gromacs") if raw else None
            variant = parse_engine_variant(raw, "gromacs", exe)
            results.append(
                {
                    "id": f"gmx-{len(results)}",
                    "label": _format_engine_label(
                        "GROMACS", version, variant, source, exe
                    ),
                    "executable": exe,
                    "version": version,
                    "variant": variant,
                    "source": source,
                    "gmxrc": gmxrc,
                    "available": True,
                }
            )
        return results

    if engine == "namd":
        candidates: List[Tuple[str, str]] = []
        for name in ("namd3", "namd2", "namd"):
            which = shutil.which(name)
            if which:
                candidates.append((which, "path"))
        search_roots = [
            Path("/usr/local"),
            Path("/opt"),
            Path.home(),
            Path.home() / "NAMD",
            Path.home() / "namd",
        ]
        for root in search_roots:
            if not root.is_dir():
                continue
            try:
                for child in root.iterdir():
                    if not child.is_dir():
                        continue
                    low = child.name.lower()
                    if "namd" not in low:
                        continue
                    for name in ("namd3", "namd2", "namd"):
                        for sub in (child, child / "bin"):
                            exe = sub / name
                            if exe.is_file() and os.access(exe, os.X_OK):
                                candidates.append((str(exe.resolve()), "prefix"))
            except OSError:
                pass

        seen: set[str] = set()
        for exe, source in candidates:
            key = str(Path(exe).resolve()) if os.path.isfile(exe) else exe
            if key in seen or not os.path.isfile(exe):
                continue
            seen.add(key)
            raw = _probe_binary_output(exe, "namd")
            version = parse_tool_version(raw, "namd") if raw else None
            if not _is_plausible_version(version):
                version = _version_from_install_path(exe, "namd")
            variant = parse_engine_variant(raw, "namd", exe)
            results.append(
                {
                    "id": f"namd-{len(results)}",
                    "label": _format_engine_label(
                        "NAMD", version, variant, source, exe
                    ),
                    "executable": exe,
                    "version": version,
                    "variant": variant,
                    "source": source,
                    "available": True,
                }
            )
        return results

    raise ValueError(f"Unsupported engine for discovery: {engine}")
