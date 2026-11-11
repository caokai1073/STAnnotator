"""What differs between Windows and macOS / Linux: file locks, whether a process still runs, ending one, and starting
helper processes (ssh, git) without a console window popping up from the windowed app.

On Windows os.kill(pid, 0) does not test a process: it ends it (TerminateProcess with exit code 0), so pid_alive()
asks the system instead.
"""
import os, subprocess, threading, time, uuid
from contextlib import contextmanager

WINDOWS = os.name == "nt"
NO_WINDOW = subprocess.CREATE_NO_WINDOW if WINDOWS else 0  # creationflags of helper processes


def _unique():
    """A temporary file's name part of its own: the process, the thread and a random part (two threads, or two calls
    of one, writing the same file must not share it: one rename would find it gone)."""
    return f"{os.getpid()}.{threading.get_ident()}.{uuid.uuid4().hex[:8]}"


def dump_json(obj, path, **kw):
    """json.dump into path through a temporary file next to it, then one rename: another process (the agent's tools,
    the app) reading it meanwhile sees the old file or the new one, never half of one. On Windows a rename onto a
    file another process has open can fail for a moment: it is tried again, then written in place."""
    import json
    from pathlib import Path
    path = Path(path)
    tmp = path.with_name(f".{path.name}.{_unique()}.tmp")
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(obj, f, **kw)
    for i in range(20):
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.05 * (i + 1))
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, **kw)
    try:
        os.remove(tmp)
    except OSError:
        pass


def release_memory():
    """Hand memory Python has freed back to the system (glibc on Linux and Apple's malloc keep freed blocks for
    reuse, so the app would still look as large to the system after its data is gone)."""
    import ctypes, gc, sys
    gc.collect()
    try:
        if sys.platform.startswith("linux"):
            ctypes.CDLL("libc.so.6").malloc_trim(0)
        elif sys.platform == "darwin":
            ctypes.CDLL("/usr/lib/libSystem.B.dylib").malloc_zone_pressure_relief(None, 0)
    except (OSError, AttributeError):
        pass


def save_counts(path, X):
    """A sparse matrix as scipy.sparse.save_npz writes it (load_npz reads it, older versions of the app too), with the
    fastest zip compression: about 4x faster to write than save_npz's default for ~25% more bytes (a 177,000 x 5,000
    panel: 10 s -> 2.4 s, 79 -> 98 MB). Written at a temporary name, then put in place."""
    import zipfile
    import numpy as np
    X = X.tocsr() if X.format not in ("csr", "csc") else X

    def write(tmp):
        arrays = {"indices": X.indices, "indptr": X.indptr, "format": np.array(X.format.encode("ascii")),
                  "shape": np.array(X.shape), "data": X.data}
        with zipfile.ZipFile(tmp, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=1) as z:
            for k, a in arrays.items():
                with z.open(k + ".npy", "w", force_zip64=True) as fh:
                    np.lib.format.write_array(fh, np.asanyarray(a), allow_pickle=False)
    replace_atomic(path, write)


def replace_atomic(path, write):
    """write(tmp_path) makes the file at a temporary name next to path, then one rename puts it in place (a reader or
    a stopped run never leaves half a file that later reads would trust). The temporary name keeps path's suffix
    (numpy adds .npy / .npz otherwise)."""
    from pathlib import Path
    path = Path(path)
    tmp = path.with_name(f".{_unique()}.{path.name}")
    try:
        write(tmp)
        for i in range(20):
            try:
                os.replace(tmp, path)
                return
            except PermissionError:  # Windows: a reader has it open for a moment
                time.sleep(0.05 * (i + 1))
        os.replace(tmp, path)
    finally:
        if tmp.exists():
            try:
                os.remove(tmp)
            except OSError:
                pass


def read_text(path):
    """A text file the user gave (a marker list, a CSV): UTF-8 (a BOM too), else the system's encoding (a file saved
    by Windows Notepad or Excel in GBK), never failing on a stray byte."""
    from pathlib import Path
    b = Path(path).expanduser().read_bytes()
    try:
        return b.decode("utf-8-sig")
    except UnicodeDecodeError:
        import locale
        return b.decode(locale.getpreferredencoding(False), errors="replace")


@contextmanager
def locked(f):
    """An exclusive lock on the open file f for the duration of the block (several processes write the same file:
    the agent, the viewer, reviewers)."""
    if WINDOWS:
        import msvcrt
        f.flush()
        here = f.tell()
        f.seek(0)
        while True:
            try:
                msvcrt.locking(f.fileno(), msvcrt.LK_NBLCK, 1)
                break
            except OSError:
                time.sleep(0.05)
        f.seek(here)
        try:
            yield
        finally:
            f.flush()
            here = f.tell()
            f.seek(0)
            msvcrt.locking(f.fileno(), msvcrt.LK_UNLCK, 1)
            f.seek(here)
    else:
        import fcntl
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def pid_alive(pid):
    """True while process pid runs (never signals it)."""
    try:
        pid = int(pid)
    except (TypeError, ValueError):
        return False
    if pid <= 0:
        return False
    if WINDOWS:
        import ctypes
        k = ctypes.windll.kernel32
        h = k.OpenProcess(0x1000, False, pid)  # PROCESS_QUERY_LIMITED_INFORMATION
        if not h:
            return False
        try:
            code = ctypes.c_ulong()
            return bool(k.GetExitCodeProcess(h, ctypes.byref(code))) and code.value == 259  # STILL_ACTIVE
        finally:
            k.CloseHandle(h)
    try:
        os.kill(pid, 0)
        return True
    except PermissionError:  # runs, under another user
        return True
    except OSError:
        return False


def terminate(pid, tree=False):
    """Ask process pid to end (Windows: taskkill, with its child processes when tree)."""
    if WINDOWS:
        subprocess.run(["taskkill", "/PID", str(int(pid)), *(["/T"] if tree else []), "/F"], capture_output=True,
                       creationflags=NO_WINDOW, check=False)
        return
    import signal
    try:
        os.kill(int(pid), signal.SIGTERM)
    except OSError:
        pass


def open_folder(path):
    """Show a folder in the system's file manager (Explorer, Finder, the desktop's)."""
    path = str(path)
    if WINDOWS:
        os.startfile(path)  # noqa: S606 (a folder, opened by Explorer)
    elif __import__("sys").platform == "darwin":
        subprocess.Popen(["open", path])
    else:
        subprocess.Popen(["xdg-open", path], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def https_context():
    """TLS for the app's own requests to the internet (updates, bug reports): the system's certificates and certifi's.
    The built app's Python may not find the system's (macOS keeps them in the Keychain, out of OpenSSL's reach);
    certifi comes with the app."""
    import ssl
    ctx = ssl.create_default_context()
    try:
        import certifi
        ctx.load_verify_locations(certifi.where())
    except Exception:
        pass
    return ctx
