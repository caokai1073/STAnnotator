"""Paths and devices.

All data lives under one home directory: $STAGENT_HOME, default ~/stagent (the older names still work:
$CELLTYPEAGENT_HOME, and an existing ~/celltypeagent folder when there is no ~/stagent)
  samples/<sample_id>/   imported samples (stagent import ...)
  sessions/<session_id>/ annotations, versions, chats
  marker_lists/          marker lists saved from the viewer
  trash/                 samples reset in the viewer
"""
import os, sys
from pathlib import Path

PACKAGE = Path(__file__).parent
def _home():
    for var in ("STAGENT_HOME", "CELLTYPEAGENT_HOME"):
        if os.environ.get(var):
            return Path(os.environ[var]).expanduser()
    new, old = Path.home() / "stagent", Path.home() / "celltypeagent"
    return old if old.is_dir() and not new.exists() else new


HOME = _home()
WORK = HOME  # older name kept for the analysis scripts
SAMPLES = HOME / "samples"
SESSIONS = HOME / "sessions"
TRASH = HOME / "trash"
USER_MARKERS = HOME / "marker_lists"
CHAT_MCP = HOME / "chat_mcp.json"
# working directory of the Claude Code CLI behind the viewer's chat box (only used when `claude` is installed)
ROOT = Path(os.environ.get("STAGENT_ROOT", PACKAGE.parent))
# H&E pixel size used when a sample's meta.json does not record one (bundles made before import existed)
HE_UM_PER_PX = 0.69
# a standalone app (PyInstaller): there is no "python -m", the app's executable takes the subcommands itself
FROZEN = bool(getattr(sys, "frozen", False))
if FROZEN:  # numba (pynndescent) caches compiled code; the app bundle is read-only
    os.environ.setdefault("NUMBA_CACHE_DIR", str(HOME / "cache" / "numba"))


def self_command(*args):
    """Command line that runs one of this installation's subcommands (mcp, viewer, ...)."""
    return [sys.executable, *args] if FROZEN else [sys.executable, "-m", "stagent.cli", *args]


def release_memory():
    """Hand memory Python has freed back to the system (the allocator keeps it otherwise, so the app looks as large
    as at its peak): malloc_zone_pressure_relief on macOS, malloc_trim on Linux."""
    import ctypes, gc
    gc.collect()
    try:
        if sys.platform == "darwin":
            ctypes.CDLL("libc.dylib").malloc_zone_pressure_relief(None, 0)
        elif sys.platform.startswith("linux"):
            ctypes.CDLL("libc.so.6").malloc_trim(0)
    except (OSError, AttributeError):
        pass


def ensure_dirs():
    for d in (SAMPLES, SESSIONS, USER_MARKERS):
        d.mkdir(parents=True, exist_ok=True)


def device():
    """'cuda' (NVIDIA) or 'mps' (Apple silicon GPU) when torch is installed (optional: it only speeds up the kNN
    graph of the marker step), else 'cpu'."""
    try:
        import torch
    except ImportError:
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return "mps"
    return "cpu"


def build():
    """Commit and date of this code: stamped into data/build.txt by git archive (export-subst), else read from git."""
    text = (PACKAGE / "data" / "build.txt").read_text().strip() if (PACKAGE / "data" / "build.txt").exists() else ""
    if text and not text.startswith("$Format"):
        return text
    try:
        import subprocess
        return subprocess.run(["git", "-C", str(PACKAGE), "log", "-1", "--format=%h %cs"], capture_output=True,
                              text=True, timeout=5).stdout.strip() or "dev"
    except Exception:
        return "dev"
