"""Parse Environment Modules / Lmod ``module avail`` output."""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from gatewizard.utils.cluster.types import ModulePackage

_ENGINE_PATTERNS: List[Tuple[str, re.Pattern[str]]] = [
    ("namd", re.compile(r"(^|/)namd(/|$)", re.I)),
    ("gromacs", re.compile(r"(^|/)(gromacs|gmx)(/|$)", re.I)),
    ("amber", re.compile(r"(^|/)(amber|pmemd|sander)(/|$)", re.I)),
    ("openmm", re.compile(r"(^|/)openmm(/|$)", re.I)),
    ("cuda", re.compile(r"(^|/)cuda(/|$)", re.I)),
]

_MODULE_TOKEN_RE = re.compile(
    r"(?P<full>(?:(?P<category>[A-Za-z0-9_.+-]+)/)?"
    r"(?P<body>[A-Za-z0-9_.+-]+(?:/[A-Za-z0-9_.+-]+)*))"
    r"(?:\s+\((?P<flags>[DL,\s]+)\))?"
)

# Help / error tokens that look like modules when ``module avail`` fails.
_GARBAGE_TOKENS = frozenset(
    {
        "no",
        "any",
        "modules",
        "module",
        "found",
        "available",
        "avail",
        "key",
        "to",
        "search",
        "use",
        "load",
        "unload",
        "list",
        "spider",
        "help",
        "where",
        "default",
        "currently",
        "loaded",
        "error",
        "command",
        "not",
        "see",
        "for",
        "more",
        "information",
        "info",
        "please",
        "try",
        "again",
        "unable",
        "unknown",
        "warning",
        "lmod",
        "tcl",
        "environment",
    }
)


def _is_plausible_module(full_name: str) -> bool:
    """Reject help-text words; require ``name/version`` or a known engine token."""
    name = (full_name or "").strip()
    if not name or len(name) < 3:
        return False
    lower = name.lower()
    if lower in _GARBAGE_TOKENS:
        return False
    if "/" in name:
        left, right = name.split("/", 1)
        if left.lower() in _GARBAGE_TOKENS or not right:
            return False
        return True
    # Bare names only if they look like an MD/CUDA engine package.
    return _infer_engine(name) is not None


def _split_version_features(body: str) -> Tuple[str, str, List[str]]:
    """Split ``namd/2.14+cuda`` → name, version, features."""
    parts = body.split("/")
    if len(parts) == 1:
        name = parts[0]
        version = ""
        features: List[str] = []
        if "+" in name:
            base, *feats = name.split("+")
            return base, "", feats
        return name, version, features
    name = parts[0]
    ver_feat = "/".join(parts[1:])
    features = []
    if "+" in ver_feat:
        version, *features = ver_feat.split("+")
    else:
        version = ver_feat
    return name, version, features


def _infer_engine(full_name: str) -> Optional[str]:
    for engine, pattern in _ENGINE_PATTERNS:
        if pattern.search(full_name):
            return engine
    return None


def parse_module_avail(text: str) -> List[ModulePackage]:
    """Parse ``module avail`` / ``ml av`` textual output into packages."""
    packages: List[ModulePackage] = []
    seen: set = set()
    for raw_line in text.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("---") or line.startswith("Where:"):
            continue
        if line.startswith("If the avail") or line.startswith("Use "):
            continue
        lower_line = line.lower()
        if "no module" in lower_line or "modules not found" in lower_line:
            continue
        # Drop leading path headers already handled by --- lines
        for token in re.split(r"\s{2,}", line):
            token = token.strip().rstrip(",")
            if not token or token.startswith("/"):
                continue
            match = _MODULE_TOKEN_RE.match(token)
            if not match:
                continue
            full = match.group("full")
            if full in seen:
                continue
            if not _is_plausible_module(full):
                continue
            seen.add(full)
            category = match.group("category") or ""
            body = match.group("body") or full
            # If category was parsed as first segment of md/namd/... keep full path
            if category and "/" not in body:
                # e.g. category=md body=namd/2.14+cuda from md/namd/2.14+cuda
                full_name = f"{category}/{body}" if not full.startswith(category) else full
            else:
                full_name = full
            # Prefer full match group
            full_name = full
            # Extract category prefix when present (md/namd/...)
            cat = ""
            name_body = full_name
            if full_name.count("/") >= 2:
                cat, name_body = full_name.split("/", 1)
            elif full_name.count("/") == 1:
                # Could be category/name or name/version — prefer no category for cuda/12.3.2
                left, right = full_name.split("/", 1)
                if left.lower() in {
                    "md",
                    "qm",
                    "mm",
                    "lang",
                    "compilers",
                    "libs",
                    "tools",
                }:
                    cat, name_body = left, right
                else:
                    name_body = full_name
            name, version, features = _split_version_features(name_body)
            flags = match.group("flags") or ""
            is_default = "D" in flags.replace(" ", "")
            packages.append(
                ModulePackage(
                    name=name,
                    full_name=full_name,
                    category=cat,
                    version=version,
                    features=features,
                    is_default=is_default,
                    engine=_infer_engine(full_name),
                )
            )
    return packages


def group_engine_modules(
    packages: List[ModulePackage],
    *,
    hints: Optional[Dict[str, List[str]]] = None,
) -> Dict[str, List[ModulePackage]]:
    """Group modules by inferred MD engine (+ cuda)."""
    grouped: Dict[str, List[ModulePackage]] = {
        "namd": [],
        "gromacs": [],
        "amber": [],
        "openmm": [],
        "cuda": [],
    }
    hint_map = hints or {}
    for pkg in packages:
        engine = pkg.engine
        if engine is None:
            for key, prefixes in hint_map.items():
                for prefix in prefixes:
                    if pkg.full_name.startswith(prefix) or prefix in pkg.full_name:
                        engine = key
                        break
                if engine:
                    break
        if engine and engine in grouped:
            grouped[engine].append(pkg)
    return grouped


def prefer_gpu_modules(
    packages: List[ModulePackage], *, want_gpu: bool
) -> List[ModulePackage]:
    """Sort packages putting GPU (+cuda) variants first when want_gpu."""
    def score(pkg: ModulePackage) -> Tuple[int, int, str]:
        has_cuda = any("cuda" in f.lower() or "gpu" in f.lower() for f in pkg.features)
        has_cuda = has_cuda or "+cuda" in pkg.full_name.lower()
        if want_gpu:
            return (0 if has_cuda else 1, 0 if pkg.is_default else 1, pkg.full_name)
        return (0 if not has_cuda else 1, 0 if pkg.is_default else 1, pkg.full_name)

    return sorted(packages, key=score)
