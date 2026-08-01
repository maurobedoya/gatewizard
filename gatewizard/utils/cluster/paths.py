"""Remote path helpers for cluster profiles."""

from __future__ import annotations

import re
from typing import Optional


_ENV_VAR_RE = re.compile(r"\$([A-Za-z_][A-Za-z0-9_]*)|\$\{([A-Za-z_][A-Za-z0-9_]*)\}")


def expand_remote_path(
    template: str,
    *,
    username: str = "",
    home: str = "",
    data_dir: str = "",
    scratch_dir: str = "",
) -> str:
    """Expand ``$USER``, ``$HOME``, ``$DATA_DIR``, ``$SCRATCH_DIR`` in a path template."""
    mapping = {
        "USER": username,
        "HOME": home or (f"/home/{username}" if username else ""),
        "DATA_DIR": data_dir or (f"/data/{username}" if username else ""),
        "SCRATCH_DIR": scratch_dir or (f"/scratch/{username}" if username else ""),
    }

    def repl(match: re.Match) -> str:
        key = match.group(1) or match.group(2)
        return mapping.get(key, match.group(0))

    return _ENV_VAR_RE.sub(repl, template or "")


def join_remote(root: str, *parts: str) -> str:
    """POSIX join for remote absolute/relative paths."""
    base = (root or "").rstrip("/")
    chunks = [p.strip("/") for p in parts if p and str(p).strip("/")]
    if not base:
        return "/".join(chunks)
    if not chunks:
        return base
    return base + "/" + "/".join(chunks)


def suggest_submit_root(
    *,
    data_dir: str = "",
    home: str = "",
    username: str = "",
    subdir: str = "gatewizard",
) -> str:
    if data_dir:
        return join_remote(data_dir, subdir)
    if home:
        return join_remote(home, subdir)
    if username:
        return f"/data/{username}/{subdir}"
    return f"$DATA_DIR/{subdir}"


def suggest_scratch_root(*, scratch_dir: str = "", username: str = "") -> str:
    if scratch_dir:
        return scratch_dir
    if username:
        return f"/scratch/{username}"
    return "$SCRATCH_DIR"
