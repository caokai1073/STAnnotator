"""What differs between Windows and macOS / Linux: file locks, whether a process still runs, ending one, and starting
helper processes (ssh, git) without a console window popping up from the windowed app.

On Windows os.kill(pid, 0) does not test a process: it ends it (TerminateProcess with exit code 0), so pid_alive()
asks the system instead.
"""
import os, subprocess, time
from contextlib import contextmanager

WINDOWS = os.name == "nt"
NO_WINDOW = subprocess.CREATE_NO_WINDOW if WINDOWS else 0  # creationflags of helper processes


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
