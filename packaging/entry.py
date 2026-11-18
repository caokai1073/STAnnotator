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


if __name__ == "__main__":
    multiprocessing.freeze_support()
    args = [a for a in sys.argv[1:] if not a.startswith("-psn_")]  # Finder may pass a process serial number
    # the macOS app opened from Finder or the Dock (no subcommand, no terminal); a subcommand keeps its output (Claude
    # Code talks to "mcp" over stdin / stdout)
    _log_streams(to_file=sys.platform == "darwin" and not args and not (sys.stdout and sys.stdout.isatty()))
    _numba_cache()
    if getattr(sys, "frozen", False):
        _mpl_cache()
    from stagent.cli import main
    main(args or ["app"])
