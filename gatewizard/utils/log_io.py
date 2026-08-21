"""Read large simulation logs without loading the entire file into memory."""

from __future__ import annotations

from pathlib import Path

DEFAULT_HEAD_BYTES = 96 * 1024
DEFAULT_TAIL_BYTES = 4 * 1024 * 1024


def read_text_head_tail(
    path: Path,
    *,
    head_bytes: int = DEFAULT_HEAD_BYTES,
    tail_bytes: int = DEFAULT_TAIL_BYTES,
) -> str:
    """Return file text, sampling head+tail when the file is large.

    Status polling only needs header parameters and the recent tail
    (last steps, timings, completion banners). Full ENERGY/mdout dumps
    in the middle are skipped.
    """
    try:
        size = path.stat().st_size
    except OSError:
        return ""
    if size <= 0:
        return ""
    if size <= head_bytes + tail_bytes:
        try:
            return path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return ""
    try:
        with path.open("rb") as handle:
            head = handle.read(head_bytes)
            handle.seek(max(0, size - tail_bytes))
            tail = handle.read(tail_bytes)
    except OSError:
        return ""
    head_text = head.decode("utf-8", errors="replace")
    tail_text = tail.decode("utf-8", errors="replace")
    nl = tail_text.find("\n")
    if nl >= 0:
        tail_text = tail_text[nl + 1 :]
    return head_text + "\n" + tail_text


def find_first_line_containing(
    path: Path,
    needles: tuple[bytes, ...],
    *,
    max_scan_bytes: int = 128 * 1024 * 1024,
    chunk_size: int = 256 * 1024,
) -> str:
    """Return the first line that contains any *needles*, scanning from the start.

    Used when a marker (e.g. GROMACS ``Started mdrun``) sits after a large
    topology dump that is neither in the log head nor the tail.
    """
    if not needles:
        return ""
    overlap = max(64, max(len(n) for n in needles))
    try:
        size = path.stat().st_size
    except OSError:
        return ""
    if size <= 0:
        return ""
    limit = min(size, max_scan_bytes)
    try:
        with path.open("rb") as handle:
            prev = b""
            scanned = 0
            while scanned < limit:
                data = handle.read(min(chunk_size, limit - scanned))
                if not data:
                    break
                buf = prev + data
                for needle in needles:
                    idx = buf.find(needle)
                    if idx < 0:
                        continue
                    end = buf.find(b"\n", idx)
                    line = buf[idx : end if end >= 0 else len(buf)]
                    return line.decode("utf-8", errors="replace").rstrip("\r")
                scanned += len(data)
                prev = data[-overlap:]
    except OSError:
        return ""
    return ""
