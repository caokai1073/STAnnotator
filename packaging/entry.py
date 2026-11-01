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


def _log_streams():
    """The windowed Windows program (STAnnotator.exe) has no console: sys.stdout and sys.stderr are None, which the
    web server's logging (uvicorn) and every print fail on. They go to app.log in _local_dir() instead."""
    if sys.stdout is not None and sys.stderr is not None:
        return
    log = _local_dir() / "app.log"
    try:
        log.parent.mkdir(parents=True, exist_ok=True)
        if log.exists() and log.stat().st_size > LOG_MB * 1e6:
            log.replace(log.with_suffix(".old.log"))
        f = open(log, "a", encoding="utf-8", buffering=1)
    except OSError:
        f = open(os.devnull, "w")
    if sys.stdout is None:
        sys.stdout = f
    if sys.stderr is None:
        sys.stderr = f


if __name__ == "__main__":
    multiprocessing.freeze_support()
    _log_streams()
    _numba_cache()
    args = [a for a in sys.argv[1:] if not a.startswith("-psn_")]  # Finder may pass a process serial number
    from stagent.cli import main
    main(args or ["app"])
