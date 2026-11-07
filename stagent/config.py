"""Paths and devices.

All data lives under one home directory, the data folder: $STAGENT_HOME when it is set, else the folder chosen in the
app (Change data folder, or the Windows installer), kept in settings_file(), else ~/stagent (the older names still
work: $CELLTYPEAGENT_HOME, and an existing ~/celltypeagent folder when there is no ~/stagent)
  samples/<sample_id>/   imported samples (stagent import ...)
  importing/<sample_id>/ a sample while it is being imported (moved to samples/ once complete)
  sessions/<session_id>/ annotations, versions, chats
  marker_lists/          marker lists saved from the viewer
  trash/                 samples reset in the viewer
"""
import json, os, sys
from pathlib import Path

PACKAGE = Path(__file__).parent


def settings_file():
    """The app's own settings (the data folder chosen in the app), outside the data folder: %APPDATA%\\STAnnotator on
    Windows (the installer writes it too), ~/Library/Application Support/STAnnotator on macOS, ~/.config/stagent."""
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA") or Path.home() / "AppData" / "Roaming") / "STAnnotator"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support" / "STAnnotator"
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME") or Path.home() / ".config") / "stagent"
    return base / "settings.json"


def read_settings():
    try:
        return json.load(open(settings_file(), encoding="utf-8-sig"))  # (the Windows installer writes a BOM)
    except (OSError, ValueError):
        return {}


def write_settings(**kw):
    """Written whole and then put in place: another process (Claude Code starting the agent's tools) reading it at
    that moment would otherwise get half a file, take it for none and use the default data folder."""
    f = settings_file()
    f.parent.mkdir(parents=True, exist_ok=True)
    tmp = f.with_name(f"{f.name}.{os.getpid()}.tmp")
    with open(tmp, "w", encoding="utf-8") as out:
        json.dump(read_settings() | kw, out, indent=1)
    os.replace(tmp, f)


def _home():
    """(data folder, what chose it: "environment", "settings" or "default")."""
    for var in ("STAGENT_HOME", "CELLTYPEAGENT_HOME"):
        if os.environ.get(var):
            return Path(os.environ[var]).expanduser(), "environment"
    chosen = read_settings().get("data_folder")
    if chosen:
        return Path(chosen).expanduser(), "settings"
    new, old = Path.home() / "stagent", Path.home() / "celltypeagent"
    return (old if old.is_dir() and not new.exists() else new), "default"


HOME, HOME_SOURCE = _home()
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
    """Command line that runs one of this installation's subcommands (mcp, viewer, ...). The Windows app is two
    programs: STAnnotator.exe (the window, no console) and stagent.exe beside it (a console program: Claude Code talks
    to its MCP server over stdin / stdout)."""
    if not FROZEN:
        return [sys.executable, "-m", "stagent.cli", *args]
    exe = Path(sys.executable)
    if os.name == "nt" and (exe.parent / "stagent.exe").exists():
        exe = exe.parent / "stagent.exe"
    return [str(exe), *args]


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
                              text=True, timeout=5, creationflags=subprocess.CREATE_NO_WINDOW if os.name == "nt" else 0
                              ).stdout.strip() or "dev"
    except Exception:
        return "dev"
