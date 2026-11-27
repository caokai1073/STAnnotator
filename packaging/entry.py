"""Entry point of the standalone app: no arguments start the desktop app; otherwise the stagent subcommands
(mcp for Claude, viewer, import, ...)."""
import multiprocessing
import os
import sys
from pathlib import Path


LOG_MB = 5  # app.log is started anew past this size


def _local_dir():
    """The app's own files on this computer (not the data): compiled code, the log of the windowed app."""
    if os.name == "nt":
        return Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local") / "STAnnotator"
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "STAnnotator"
    return Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "STAnnotator"


def _numba_cache():
    """numba keeps the code it compiles (pynndescent's neighbor search, UMAP) next to the modules, which in the app
    are read-only (Program Files, a signed .app): without a cache folder of its own it compiles again on every start."""
    os.environ.setdefault("NUMBA_CACHE_DIR", str(_local_dir() / "numba"))


def _mpl_cache():
    """PyInstaller's run-time hook gives matplotlib a new temporary config folder at every start (for one-file apps,
    whose files move), so it built its font cache again each time ("Matplotlib is building the font cache" in every
    tester's log). This app is a folder that stays where it is installed: a folder of our own keeps the cache, one per
    install place (a cache of fonts in a moved or older app would point at files that are not there)."""
    import hashlib
    where = hashlib.sha1(str(Path(sys.executable).resolve().parent).encode("utf-8")).hexdigest()[:10]
    d = _local_dir() / "matplotlib" / where
    try:
        d.mkdir(parents=True, exist_ok=True)
        os.environ["MPLCONFIGDIR"] = str(d)
    except OSError:  # (not writable: the hook's temporary folder stays)
        pass


def _log_streams(to_file=False):
    """The app's output goes to app.log in _local_dir(): always for the windowed Windows program (STAnnotator.exe has
    no console: sys.stdout and sys.stderr are None, which the web server's logging (uvicorn) and every print fail on),
    and with to_file for the macOS app opened from Finder or the Dock, whose output the system drops. There the file
    also takes the process's own output (C libraries, child processes) and, through faulthandler, the Python stack of
    a crash. A Report bug takes the end of this log."""
    if not to_file and sys.stdout is not None and sys.stderr is not None:
        return
    log = _local_dir() / "app.log"
    try:
        log.parent.mkdir(parents=True, exist_ok=True)
        if log.exists() and log.stat().st_size > LOG_MB * 1e6:
            log.replace(log.with_suffix(".old.log"))
    except OSError:  # another instance holds it (Windows): keep writing to it
        pass
    try:
        f = open(log, "a", encoding="utf-8", buffering=1)
    except OSError:
        f = open(os.devnull, "w")
    if to_file:
        try:
            for fd in (1, 2):
                os.dup2(f.fileno(), fd)
        except OSError:
            pass
    if to_file or sys.stdout is None:
        sys.stdout = f
    if to_file or sys.stderr is None:
        sys.stderr = f
    try:
        import faulthandler
        faulthandler.enable(file=f)
    except Exception:
        pass
    import time
    print(f"== STAnnotator started {time.strftime('%Y-%m-%d %H:%M:%S')}", flush=True)


RELEASES = "https://github.com/caokai1073/STAnnotator/releases/latest"


def _message(text):
    """A message box without the app's libraries (they may be what is missing): Windows' own, or AppleScript's."""
    try:
        if os.name == "nt":
            import ctypes
            ctypes.windll.user32.MessageBoxW(None, text, "STAnnotator", 0x40 | 0x10000)  # (information, in front)
        elif sys.platform == "darwin":
            import subprocess
            t = text.replace("\\", "\\\\").replace('"', '\\"')
            subprocess.run(["osascript", "-e", f'display dialog "{t}" buttons {{"OK"}} default button 1 with title "STAnnotator"'],
                           timeout=600)
    except Exception:
        pass
    print(text, flush=True)


def _update_installing():
    """The Windows installer of an update is replacing the app's files (stagent/updater.py marks it, the installer's
    script removes the mark and opens the app when it is done): a start meanwhile would fail on files not there yet.
    A mark older than its limit is from an installation that did not end: not a reason to wait."""
    import tempfile, time
    f = Path(tempfile.gettempdir()) / "STAnnotator-update" / "installing.json"
    try:
        if not (f.exists() and time.time() - f.stat().st_mtime < 20 * 60):
            return False
    except OSError:
        return False
    if os.name == "nt":  # the installer (STAnnotator-Setup.exe and the .tmp it runs) still there? Else it stopped early
        try:
            import subprocess
            out = subprocess.run(["tasklist", "/FO", "CSV", "/NH"], capture_output=True, text=True, timeout=10,
                                 creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0)).stdout
            return "stannotator-setup" in out.lower()
        except Exception:
            return True
    return True


def _installation_whole():
    """The files a started app needs are in its folder (an update stopped half way leaves some out)."""
    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))  # (the app's folder, or the repository)
    return (base / "stagent" / "static" / "app.html").exists()


if __name__ == "__main__":
    multiprocessing.freeze_support()
    args = [a for a in sys.argv[1:] if not a.startswith("-psn_")]  # Finder may pass a process serial number
    # the macOS app opened from Finder or the Dock (no subcommand, no terminal); a subcommand keeps its output (Claude
    # Code talks to "mcp" over stdin / stdout)
    _log_streams(to_file=sys.platform == "darwin" and not args and not (sys.stdout and sys.stdout.isatty()))
    _numba_cache()
    frozen, gui = getattr(sys, "frozen", False), not args or args == ["app"]
    if frozen and gui and _update_installing():
        _message("An update of STAnnotator is being installed. The app opens by itself when the installation has "
                 "finished (a minute or two).")
        sys.exit(0)
    if frozen and gui and not _installation_whole():
        _message("STAnnotator's installation is incomplete (an update may have stopped half way). Please download "
                 f"the installer again and run it:\n{RELEASES}\nYour data stays where it is.")
        sys.exit(1)
    if frozen:
        _mpl_cache()
    try:
        from stagent.cli import main
    except Exception:  # a module or one of its files missing: an installation not whole, said instead of no window
        if not (frozen and gui):
            raise
        import traceback
        traceback.print_exc()
        _message("STAnnotator could not start: part of its installation is missing (an update may have stopped half "
                 f"way). Please download the installer again and run it:\n{RELEASES}\nYour data stays where it is.")
        sys.exit(1)
    main(args or ["app"])
