"""Entry point of the standalone app: no arguments start the desktop app; otherwise the stagent subcommands
(mcp for Claude, viewer, import, ...)."""
import multiprocessing
import os
import sys
from pathlib import Path


def _numba_cache():
    """numba keeps the code it compiles (pynndescent's neighbor search, UMAP) next to the modules, which in the app
    are read-only (Program Files, a signed .app): without a cache folder of its own it compiles again on every start."""
    if os.name == "nt":
        base = Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local") / "STAnnotator"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Caches" / "STAnnotator"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "STAnnotator"
    os.environ.setdefault("NUMBA_CACHE_DIR", str(base / "numba"))


if __name__ == "__main__":
    multiprocessing.freeze_support()
    _numba_cache()
    args = [a for a in sys.argv[1:] if not a.startswith("-psn_")]  # Finder may pass a process serial number
    from stagent.cli import main
    main(args or ["app"])
