"""SSH / rsync helpers for remote cluster sessions.

Prefer system ``ssh`` / ``rsync`` (keys + ControlMaster). Password auth uses
optional Paramiko and keeps the password only in process memory.
"""

from __future__ import annotations

import os
import re
import shlex
import subprocess
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

ProgressCallback = Callable[[Dict[str, object]], None]

# rsync --info=progress2 lines look like:
#   12,345,678  45%  10.50MB/s    0:01:23 (xfr#5, to-chk=10/100)
_RSYNC_PROGRESS_RE = re.compile(
    r"(?P<bytes>[\d,]+)\s+(?P<pct>\d+)%\s+(?P<speed>[\d.]+[kKMGT]?B/s)\s+(?P<eta>\d+:\d+(?::\d+)?)"
)


class ClusterSSHError(RuntimeError):
    """Remote connection / command failure."""


def parse_rsync_progress_line(line: str) -> Optional[Dict[str, object]]:
    """Parse one rsync ``--info=progress2`` stderr line into a progress dict."""
    text = (line or "").strip().replace("\r", "")
    if not text:
        return None
    match = _RSYNC_PROGRESS_RE.search(text)
    if not match:
        return None
    try:
        transferred = int(match.group("bytes").replace(",", ""))
    except ValueError:
        transferred = 0
    try:
        percent = max(0, min(100, int(match.group("pct"))))
    except ValueError:
        percent = 0
    return {
        "phase": "sync",
        "percent": percent,
        "bytes": transferred,
        "speed": match.group("speed"),
        "eta": match.group("eta"),
        "message": f"Downloading… {percent}% ({match.group('speed')})",
    }


def format_byte_size(num: int) -> str:
    """Human-readable byte size for progress messages."""
    n = float(max(0, int(num)))
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024.0 or unit == "TB":
            if unit == "B":
                return f"{int(n)} {unit}"
            return f"{n:.1f} {unit}"
        n /= 1024.0
    return f"{int(num)} B"


def compute_rsync_timeout(
    expected_bytes: Optional[int] = None,
    *,
    base: int = 600,
    max_wall: int = 86400,
    bytes_per_second_floor: int = 256 * 1024,
) -> int:
    """Wall-clock limit for rsync, scaled when the remote payload size is known.

    Uses a conservative floor transfer rate (256 KiB/s by default) so large
    trajectory pulls are not cut off at the small-job default (600 s).
    """
    if not expected_bytes or expected_bytes <= 0:
        return base
    scaled = base + int(expected_bytes / max(1, bytes_per_second_floor))
    return max(base, min(max_wall, scaled))


def local_dir_byte_size(
    local_dir: str | Path, *, excludes: Optional[List[str]] = None
) -> int:
    """Sum file sizes under *local_dir*, honoring simple exclude names/globs."""
    root = Path(local_dir)
    if not root.is_dir():
        return 0
    total = 0
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        try:
            rel = path.relative_to(root).as_posix()
        except ValueError:
            continue
        if _path_excluded(rel, excludes) or _path_excluded(path.name, excludes):
            continue
        try:
            total += path.stat().st_size
        except OSError:
            continue
    return total


def remote_dir_byte_size(session_id: str, remote_dir: str) -> int:
    """Return total bytes under a remote directory (``du -sb``), or 0 on failure."""
    path = (remote_dir or "").strip().rstrip("/")
    if not path:
        return 0
    path_q = shlex.quote(path)
    try:
        _rc, out, _err = run_remote(
            session_id,
            f"du -sb {path_q} 2>/dev/null | awk '{{print $1}}'",
            timeout=90,
        )
    except ClusterSSHError:
        return 0
    text = (out or "").strip().splitlines()
    if not text:
        return 0
    try:
        return max(0, int(text[-1].strip()))
    except ValueError:
        return 0


@dataclass
class SSHSession:
    """In-memory SSH session (never persisted)."""

    session_id: str
    host: str
    username: str
    port: int = 22
    identity_file: str = ""
    control_path: Optional[str] = None
    password: Optional[str] = None  # memory only
    use_paramiko: bool = False
    created_at: float = field(default_factory=time.time)
    last_used: float = field(default_factory=time.time)


_SESSIONS: Dict[str, SSHSession] = {}
_LOCK = threading.Lock()


def _expand_identity(path: str) -> str:
    if not path:
        return ""
    return str(Path(path).expanduser())


def _ssh_base_args(session: SSHSession) -> List[str]:
    args = [
        "ssh",
        "-p",
        str(session.port),
        "-o",
        "BatchMode=yes" if not session.password else "BatchMode=no",
        "-o",
        "StrictHostKeyChecking=accept-new",
        "-o",
        "ConnectTimeout=20",
    ]
    identity = _expand_identity(session.identity_file)
    if identity:
        args.extend(["-i", identity])
    if session.control_path:
        args.extend(
            [
                "-o",
                "ControlMaster=auto",
                "-o",
                f"ControlPath={session.control_path}",
                "-o",
                "ControlPersist=10m",
            ]
        )
    args.append(f"{session.username}@{session.host}")
    return args


def connect_ssh(
    *,
    host: str,
    username: str,
    port: int = 22,
    identity_file: str = "",
    password: Optional[str] = None,
) -> SSHSession:
    """Open a reusable SSH session. Password is never written to disk."""
    if not host or not username:
        raise ClusterSSHError("host and username are required")

    session_id = uuid.uuid4().hex
    control_dir = Path(tempfile.gettempdir()) / "gatewizard-ssh"
    control_dir.mkdir(mode=0o700, exist_ok=True)
    control_path = str(control_dir / f"cm-{session_id}")

    session = SSHSession(
        session_id=session_id,
        host=host,
        username=username,
        port=int(port or 22),
        identity_file=identity_file or "",
        control_path=control_path,
        password=password,
        use_paramiko=bool(password),
    )

    if password:
        _connect_paramiko(session)
    else:
        # Establish ControlMaster with a trivial command
        result = _run_system_ssh(session, ["echo", "ok"], timeout=30)
        if result[0] != 0:
            raise ClusterSSHError(result[2] or result[1] or "SSH connection failed")

    with _LOCK:
        _SESSIONS[session_id] = session
    return session


def get_session(session_id: str) -> SSHSession:
    with _LOCK:
        session = _SESSIONS.get(session_id)
    if not session:
        raise ClusterSSHError("SSH session expired or not found; connect again")
    session.last_used = time.time()
    return session


def close_session(session_id: str) -> None:
    with _LOCK:
        session = _SESSIONS.pop(session_id, None)
    if not session:
        return
    session.password = None
    if session.control_path and not session.use_paramiko:
        try:
            subprocess.run(
                [
                    "ssh",
                    "-O",
                    "exit",
                    "-o",
                    f"ControlPath={session.control_path}",
                    f"{session.username}@{session.host}",
                ],
                capture_output=True,
                timeout=5,
                check=False,
            )
        except Exception:
            pass
        try:
            Path(session.control_path).unlink(missing_ok=True)
        except Exception:
            pass


def run_remote(
    session_id: str,
    command: str,
    *,
    timeout: int = 120,
) -> Tuple[int, str, str]:
    """Run a remote shell command; returns (rc, stdout, stderr).

    Uses a login bash (``bash -lc``) so Environment Modules / Lmod from
    ``/etc/profile`` are available — including Paramiko password sessions.

    The remote command is always passed as **one** shell-quoted string so SSH
    does not re-split spaces (which breaks ``for`` loops / pipelines).
    """
    session = get_session(session_id)
    # -l: login (profile). Avoid -i: it adds ioctl noise without helping here.
    wrapped = f"bash -lc {shlex.quote(command)}"
    if session.use_paramiko:
        return _run_paramiko(session, wrapped, timeout=timeout)
    # Single argv: OpenSSH joins remote args with spaces then re-parses via the
    # user shell; multiple argv would break commands that contain spaces.
    return _run_system_ssh(session, [wrapped], timeout=timeout)


def rsync_to_remote(
    session_id: str,
    local_dir: str,
    remote_dir: str,
    *,
    delete: bool = False,
    excludes: Optional[List[str]] = None,
    timeout: int = 600,
) -> Tuple[int, str, str]:
    session = get_session(session_id)
    local = str(Path(local_dir))
    if not local.endswith(os.sep):
        local = local + os.sep
    # Ensure remote parent exists
    run_remote(session_id, f"mkdir -p {shlex.quote(remote_dir)}", timeout=60)
    if session.use_paramiko:
        return _rsync_via_paramiko_sftp(session, local_dir, remote_dir, excludes=excludes)

    ssh_cmd = "ssh"
    ssh_opts = [
        "-p",
        str(session.port),
        "-o",
        "StrictHostKeyChecking=accept-new",
    ]
    identity = _expand_identity(session.identity_file)
    if identity:
        ssh_opts.extend(["-i", identity])
    if session.control_path:
        ssh_opts.extend(
            [
                "-o",
                f"ControlPath={session.control_path}",
                "-o",
                "ControlMaster=auto",
                "-o",
                "ControlPersist=10m",
            ]
        )
    remote = f"{session.username}@{session.host}:{remote_dir}/"
    cmd = ["rsync", "-az", "-e", " ".join([ssh_cmd] + ssh_opts)]
    if delete:
        cmd.append("--delete")
    for ex in excludes or []:
        cmd.extend(["--exclude", ex])
    cmd.extend([local, remote])
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def rsync_from_remote(
    session_id: str,
    remote_dir: str,
    local_dir: str,
    *,
    excludes: Optional[List[str]] = None,
    includes: Optional[List[str]] = None,
    timeout: int = 600,
    idle_timeout: int = 600,
    ignore_times: bool = False,
    on_progress: Optional[ProgressCallback] = None,
    expected_bytes: Optional[int] = None,
) -> Tuple[int, str, str]:
    session = get_session(session_id)
    Path(local_dir).mkdir(parents=True, exist_ok=True)
    local = str(Path(local_dir))
    if not local.endswith(os.sep):
        local = local + os.sep
    if session.use_paramiko:
        return _rsync_from_paramiko_sftp(
            session,
            remote_dir,
            local_dir,
            excludes=excludes,
            includes=includes,
            on_progress=on_progress,
        )

    ssh_opts = [
        "ssh",
        "-p",
        str(session.port),
        "-o",
        "StrictHostKeyChecking=accept-new",
    ]
    identity = _expand_identity(session.identity_file)
    if identity:
        ssh_opts.extend(["-i", identity])
    if session.control_path:
        ssh_opts.extend(
            [
                "-o",
                f"ControlPath={session.control_path}",
                "-o",
                "ControlMaster=auto",
                "-o",
                "ControlPersist=10m",
            ]
        )
    remote = f"{session.username}@{session.host}:{remote_dir}/"
    # Omit -z during progress pulls: compression delays / mutes progress2 on pipes.
    cmd = ["rsync", "-a", "--info=progress2", "--outbuf=N", "-e", " ".join(ssh_opts)]
    # WSL/OneDrive mounts often have unreliable mtimes; growing stage logs then
    # look "up to date" locally while the cluster has advanced several steps.
    if ignore_times:
        cmd.append("--ignore-times")
    # Includes must come before excludes; trailing exclude '*' makes an allow-list.
    for inc in includes or []:
        cmd.extend(["--include", inc])
    ex_list = list(excludes or [])
    if includes and "*" not in ex_list:
        ex_list.append("*")
    for ex in ex_list:
        cmd.extend(["--exclude", ex])
    cmd.extend([remote, local])

    wall_timeout = compute_rsync_timeout(expected_bytes, base=timeout)
    idle_limit = max(60, int(idle_timeout))

    if on_progress is None:
        proc = subprocess.run(
            cmd, capture_output=True, text=True, timeout=wall_timeout, check=False
        )
        return proc.returncode, proc.stdout or "", proc.stderr or ""

    total_bytes = int(expected_bytes or 0)
    if total_bytes <= 0:
        total_bytes = remote_dir_byte_size(session_id, remote_dir)
    on_progress(
        {
            "phase": "sync",
            "percent": 0,
            "bytes": 0,
            "total_bytes": total_bytes or None,
            "message": (
                f"Downloading… 0 / {format_byte_size(total_bytes)}"
                if total_bytes > 0
                else "Starting download…"
            ),
        }
    )
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            # Binary pipes stay unbuffered enough for live size polling + progress2.
            text=False,
            bufsize=0,
        )
    except FileNotFoundError as exc:
        raise ClusterSSHError("rsync not found on PATH") from exc

    stderr_chunks: List[str] = []
    stdout_chunks: List[str] = []
    stop = threading.Event()
    last_pct = -1
    lock = threading.Lock()
    progress2_seen = threading.Event()
    start_local = local_dir_byte_size(local_dir, excludes=ex_list)
    start_time = time.time()
    activity = {"time": start_time, "bytes": start_local}

    def _mark_activity(bytes_hint: Optional[int] = None) -> None:
        activity["time"] = time.time()
        if isinstance(bytes_hint, int) and bytes_hint >= 0:
            activity["bytes"] = max(activity["bytes"], bytes_hint)

    def _emit(evt: Dict[str, object]) -> None:
        nonlocal last_pct
        transferred = evt.get("bytes")
        if isinstance(transferred, int):
            _mark_activity(transferred)
        else:
            _mark_activity()
        with lock:
            pct_raw = evt.get("percent")
            pct = int(pct_raw) if isinstance(pct_raw, (int, float)) else -1
            # Throttle duplicate percentages (size poller can fire often).
            if pct == last_pct and pct not in {0, 100}:
                return
            if pct >= 0:
                last_pct = pct
            on_progress(evt)

    def _progress_from_transferred(transferred: int) -> Dict[str, object]:
        if total_bytes > 0:
            pct = max(0, min(99, int(100.0 * transferred / total_bytes)))
            return {
                "phase": "sync",
                "percent": pct,
                "bytes": transferred,
                "total_bytes": total_bytes,
                "message": (
                    f"Downloading… {format_byte_size(transferred)} / "
                    f"{format_byte_size(total_bytes)} ({pct}%)"
                ),
            }
        return {
            "phase": "sync",
            "percent": 0,
            "bytes": transferred,
            "message": f"Downloading… {format_byte_size(transferred)}",
        }

    def _consume_progress_line(line: str) -> None:
        evt = parse_rsync_progress_line(line)
        if evt is None:
            return
        progress2_seen.set()
        transferred = evt.get("bytes")
        if isinstance(transferred, int) and transferred >= 0:
            base = _progress_from_transferred(transferred)
            if evt.get("speed"):
                base["speed"] = evt["speed"]
                base["message"] = f"{base['message']} · {evt['speed']}"
            if evt.get("eta"):
                base["eta"] = evt["eta"]
            _emit(base)
        else:
            _emit(evt)

    def _read_stdout() -> None:
        assert proc.stdout is not None
        # Some rsync builds emit progress2 on stdout when not a TTY.
        buf = ""
        try:
            while True:
                raw = proc.stdout.read(1)
                if not raw:
                    break
                ch = raw.decode("utf-8", errors="replace")
                if ch in {"\r", "\n"}:
                    if buf:
                        stdout_chunks.append(buf + "\n")
                        _consume_progress_line(buf)
                        buf = ""
                    continue
                buf += ch
            if buf:
                stdout_chunks.append(buf)
                _consume_progress_line(buf)
        except Exception:
            pass

    def _read_stderr() -> None:
        assert proc.stderr is not None
        buf = ""
        try:
            while True:
                raw = proc.stderr.read(1)
                if not raw:
                    break
                ch = raw.decode("utf-8", errors="replace")
                if ch in {"\r", "\n"}:
                    if buf:
                        stderr_chunks.append(buf + "\n")
                        _consume_progress_line(buf)
                        buf = ""
                    continue
                buf += ch
            if buf:
                stderr_chunks.append(buf)
                _consume_progress_line(buf)
        except Exception:
            pass

    def _poll_local_size() -> None:
        """Fallback progress when progress2 is silent: local size vs remote total.

        When local is already nearly full (typical ignore-times re-pull), size
        stays flat — emit a status heartbeat until progress2 reports bytes.
        """
        if total_bytes <= 0:
            return
        nearly_full = start_local >= int(total_bytes * 0.95)
        while not stop.wait(0.4):
            if progress2_seen.is_set():
                continue
            if nearly_full:
                _emit(
                    {
                        "phase": "sync",
                        "percent": max(1, last_pct) if last_pct > 0 else 1,
                        "bytes": start_local,
                        "total_bytes": total_bytes,
                        "message": (
                            f"Re-syncing… {format_byte_size(start_local)} on disk / "
                            f"{format_byte_size(total_bytes)} remote "
                            "(waiting for transfer counters)"
                        ),
                    }
                )
                continue
            cur = local_dir_byte_size(local_dir, excludes=ex_list)
            if cur > activity["bytes"]:
                _mark_activity(cur)
            pct = max(0, min(99, int(100.0 * cur / total_bytes)))
            _emit(
                {
                    "phase": "sync",
                    "percent": pct,
                    "bytes": cur,
                    "total_bytes": total_bytes,
                    "message": (
                        f"Downloading… {format_byte_size(cur)} / "
                        f"{format_byte_size(total_bytes)} ({pct}%)"
                    ),
                }
            )

    out_thread = threading.Thread(target=_read_stdout, daemon=True)
    err_thread = threading.Thread(target=_read_stderr, daemon=True)
    size_thread = threading.Thread(target=_poll_local_size, daemon=True)
    out_thread.start()
    err_thread.start()
    size_thread.start()

    while proc.poll() is None:
        now = time.time()
        if now - activity["time"] > idle_limit:
            proc.kill()
            stop.set()
            transferred = max(0, local_dir_byte_size(local_dir, excludes=ex_list) - start_local)
            raise ClusterSSHError(
                f"rsync stalled (no progress for {idle_limit}s; "
                f"{format_byte_size(transferred)} transferred)"
            )
        if now - start_time > wall_timeout:
            proc.kill()
            stop.set()
            transferred = max(0, local_dir_byte_size(local_dir, excludes=ex_list) - start_local)
            raise ClusterSSHError(
                f"rsync timed out after {wall_timeout}s "
                f"({format_byte_size(transferred)} transferred)"
            )
        time.sleep(0.2)

    stop.set()
    out_thread.join(timeout=5)
    err_thread.join(timeout=5)
    size_thread.join(timeout=2)
    rc = proc.returncode if proc.returncode is not None else 0
    if rc == 0:
        cur = local_dir_byte_size(local_dir, excludes=ex_list)
        _emit(
            {
                "phase": "sync",
                "percent": 100,
                "bytes": cur,
                "total_bytes": total_bytes or cur,
                "message": (
                    f"Download complete ({format_byte_size(cur)})"
                    if cur
                    else "Download complete"
                ),
            }
        )
    return rc, "".join(stdout_chunks), "".join(stderr_chunks)


def _run_system_ssh(
    session: SSHSession, remote_argv: List[str], *, timeout: int
) -> Tuple[int, str, str]:
    cmd = _ssh_base_args(session) + remote_argv
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise ClusterSSHError(f"SSH timed out after {timeout}s") from exc
    except FileNotFoundError as exc:
        raise ClusterSSHError("ssh client not found on PATH") from exc
    return proc.returncode, proc.stdout or "", proc.stderr or ""


def _connect_paramiko(session: SSHSession) -> None:
    try:
        import paramiko
    except ImportError as exc:
        raise ClusterSSHError(
            "Password authentication requires the 'paramiko' package. "
            "Install it or use an SSH key / ssh-agent instead."
        ) from exc
    client = paramiko.SSHClient()
    client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    connect_kwargs = {
        "hostname": session.host,
        "port": session.port,
        "username": session.username,
        "password": session.password,
        "timeout": 20,
        "allow_agent": True,
        "look_for_keys": True,
    }
    identity = _expand_identity(session.identity_file)
    if identity:
        connect_kwargs["key_filename"] = identity
    try:
        client.connect(**connect_kwargs)
    except Exception as exc:
        raise ClusterSSHError(f"SSH connection failed: {exc}") from exc
    # Stash client on session via private attribute
    session._client = client  # type: ignore[attr-defined]


def _run_paramiko(
    session: SSHSession, command: str, *, timeout: int
) -> Tuple[int, str, str]:
    client = getattr(session, "_client", None)
    if client is None:
        _connect_paramiko(session)
        client = getattr(session, "_client", None)
    assert client is not None
    try:
        _stdin, stdout, stderr = client.exec_command(command, timeout=timeout)
        out = stdout.read().decode("utf-8", errors="replace")
        err = stderr.read().decode("utf-8", errors="replace")
        rc = stdout.channel.recv_exit_status()
        return rc, out, err
    except Exception as exc:
        raise ClusterSSHError(f"Remote command failed: {exc}") from exc


def _path_excluded(rel: str, excludes: Optional[List[str]]) -> bool:
    """Return True if relative path matches an exclude pattern."""
    if not excludes:
        return False
    name = Path(rel).name
    for pattern in excludes:
        pat = str(pattern)
        if not pat:
            continue
        if pat.startswith("*.") and name.endswith(pat[1:]):
            return True
        if pat in {rel, name}:
            return True
        if "/" in pat and (rel == pat or rel.endswith("/" + pat) or rel.startswith(pat.rstrip("*"))):
            return True
    return False


def _rsync_via_paramiko_sftp(
    session: SSHSession,
    local_dir: str,
    remote_dir: str,
    *,
    excludes: Optional[List[str]] = None,
) -> Tuple[int, str, str]:
    """Recursive upload via SFTP when rsync+password is unavailable."""
    try:
        import paramiko  # noqa: F401
    except ImportError as exc:
        raise ClusterSSHError("paramiko required for password-based file transfer") from exc
    client = getattr(session, "_client", None)
    if client is None:
        _connect_paramiko(session)
        client = getattr(session, "_client")
    local_root = Path(local_dir)
    if not local_root.is_dir():
        raise ClusterSSHError(f"Local directory not found: {local_dir}")

    expected_files = [
        p
        for p in local_root.rglob("*")
        if p.is_file() and not _path_excluded(p.relative_to(local_root).as_posix(), excludes)
    ]
    if not expected_files:
        raise ClusterSSHError(
            f"No files to upload from {local_dir} (directory empty or everything excluded)"
        )

    sftp = client.open_sftp()
    uploaded = 0
    errors: List[str] = []
    try:
        _mkdir_p_sftp(sftp, remote_dir)
        for path in expected_files:
            rel = path.relative_to(local_root).as_posix()
            remote_path = f"{remote_dir.rstrip('/')}/{rel}"
            parent = str(Path(remote_path).parent).replace("\\", "/")
            try:
                _mkdir_p_sftp(sftp, parent)
                sftp.put(str(path), remote_path)
                uploaded += 1
            except Exception as exc:
                errors.append(f"{rel}: {exc}")
    finally:
        sftp.close()

    if errors:
        raise ClusterSSHError(
            f"SFTP upload incomplete ({uploaded}/{len(expected_files)} files). "
            + "; ".join(errors[:5])
        )
    if uploaded == 0:
        raise ClusterSSHError(
            f"SFTP uploaded 0 files from {local_dir} (expected {len(expected_files)})"
        )
    return 0, f"uploaded {uploaded} files via SFTP", ""


def verify_remote_files(
    session_id: str,
    remote_dir: str,
    filenames: List[str],
) -> Tuple[bool, str]:
    """Check that required files exist on the remote host. Returns (ok, message)."""
    missing: List[str] = []
    for name in filenames:
        remote_path = f"{remote_dir.rstrip('/')}/{name}"
        rc, out, err = run_remote(
            session_id,
            f"test -f {shlex.quote(remote_path)} && echo OK || echo MISSING",
            timeout=30,
        )
        text = (out or err or "").strip()
        if rc != 0 or "OK" not in text:
            missing.append(name)
    if missing:
        rc2, listing, _ = run_remote(
            session_id,
            f"ls -la {shlex.quote(remote_dir)} 2>&1 | head -40",
            timeout=30,
        )
        return (
            False,
            f"Missing on remote after upload: {', '.join(missing)}. "
            f"Remote listing:\n{listing or '(empty)'}",
        )
    return True, "ok"


def remote_file_count(session_id: str, remote_dir: str) -> int:
    rc, out, _ = run_remote(
        session_id,
        f"find {shlex.quote(remote_dir)} -type f 2>/dev/null | wc -l",
        timeout=60,
    )
    try:
        return int((out or "0").strip().split()[0])
    except Exception:
        return 0


def _rsync_from_paramiko_sftp(
    session: SSHSession,
    remote_dir: str,
    local_dir: str,
    *,
    excludes: Optional[List[str]] = None,
    includes: Optional[List[str]] = None,
    on_progress: Optional[ProgressCallback] = None,
) -> Tuple[int, str, str]:
    client = getattr(session, "_client", None)
    if client is None:
        _connect_paramiko(session)
        client = getattr(session, "_client")
    sftp = client.open_sftp()
    downloaded = 0
    try:
        if on_progress:
            on_progress(
                {
                    "phase": "sync",
                    "percent": 0,
                    "message": "Downloading via SFTP…",
                }
            )
        downloaded = _download_dir_sftp(
            sftp,
            remote_dir,
            local_dir,
            excludes=excludes or [],
            includes=includes,
            on_progress=on_progress,
        )
    finally:
        sftp.close()
    if on_progress:
        on_progress(
            {
                "phase": "sync",
                "percent": 100,
                "message": f"Downloaded {downloaded} files",
            }
        )
    return 0, f"downloaded files into {local_dir} ({downloaded} local files)", ""


def _mkdir_p_sftp(sftp, remote_path: str) -> None:
    parts = remote_path.strip("/").split("/")
    cur = ""
    for part in parts:
        cur = f"{cur}/{part}"
        try:
            sftp.stat(cur)
        except IOError:
            try:
                sftp.mkdir(cur)
            except IOError:
                pass


def _name_matches_includes(name: str, includes: Optional[List[str]]) -> bool:
    """True when *name* matches an allow-list glob (ignores directory ``*/`` entries)."""
    if not includes:
        return True
    import fnmatch

    for pat in includes:
        if pat in {"*/", "*"}:
            continue
        if fnmatch.fnmatch(name, pat):
            return True
    return False


def _count_remote_files_sftp(
    sftp, remote_dir: str, *, excludes: List[str], includes: Optional[List[str]] = None
) -> int:
    import stat as statmod

    total = 0
    try:
        entries = sftp.listdir_attr(remote_dir)
    except IOError:
        return 0
    for entry in entries:
        name = entry.filename
        if name in {".", ".."}:
            continue
        if _path_excluded(name, excludes) or name in excludes:
            continue
        rpath = f"{remote_dir.rstrip('/')}/{name}"
        if statmod.S_ISDIR(entry.st_mode):
            total += _count_remote_files_sftp(
                sftp, rpath, excludes=excludes, includes=includes
            )
        elif _name_matches_includes(name, includes):
            total += 1
    return total


def _download_dir_sftp(
    sftp,
    remote_dir: str,
    local_dir: str,
    *,
    excludes: List[str],
    includes: Optional[List[str]] = None,
    on_progress: Optional[ProgressCallback] = None,
    _state: Optional[Dict[str, int]] = None,
) -> int:
    Path(local_dir).mkdir(parents=True, exist_ok=True)
    import stat as statmod

    if _state is None:
        total = max(
            1,
            _count_remote_files_sftp(
                sftp, remote_dir, excludes=excludes, includes=includes
            ),
        )
        _state = {"done": 0, "total": total}

    for entry in sftp.listdir_attr(remote_dir):
        name = entry.filename
        if name in {".", ".."}:
            continue
        if _path_excluded(name, excludes) or name in excludes:
            continue
        rpath = f"{remote_dir.rstrip('/')}/{name}"
        lpath = str(Path(local_dir) / name)
        if statmod.S_ISDIR(entry.st_mode):
            _download_dir_sftp(
                sftp,
                rpath,
                lpath,
                excludes=excludes,
                includes=includes,
                on_progress=on_progress,
                _state=_state,
            )
        elif _name_matches_includes(name, includes):
            sftp.get(rpath, lpath)
            _state["done"] += 1
            if on_progress:
                done = _state["done"]
                total = max(1, _state["total"])
                pct = max(0, min(99, int(100 * done / total)))
                on_progress(
                    {
                        "phase": "sync",
                        "percent": pct,
                        "message": f"Downloading via SFTP… {done}/{total}",
                        "files_done": done,
                        "files_total": total,
                    }
                )
    return int(_state["done"])
