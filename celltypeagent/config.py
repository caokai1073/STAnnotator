"""Paths and devices.

All data lives under one home directory: $CELLTYPEAGENT_HOME, default ~/celltypeagent
  samples/<sample_id>/   imported samples (celltypeagent import ...)
  sessions/<session_id>/ annotations, versions, chats
  marker_lists/          marker lists saved from the viewer
  trash/                 samples reset in the viewer
"""
import os
from pathlib import Path

PACKAGE = Path(__file__).parent
HOME = Path(os.environ.get("CELLTYPEAGENT_HOME", Path.home() / "celltypeagent")).expanduser()
WORK = HOME  # older name kept for the analysis scripts
SAMPLES = HOME / "samples"
SESSIONS = HOME / "sessions"
TRASH = HOME / "trash"
USER_MARKERS = HOME / "marker_lists"
CHAT_MCP = HOME / "chat_mcp.json"
DEFAULT_MARKERS = PACKAGE / "data" / "ibd_xenium480_draft.txt"
# working directory of the Claude Code CLI behind the viewer's chat box (only used when `claude` is installed)
ROOT = Path(os.environ.get("CELLTYPEAGENT_ROOT", PACKAGE.parent))
# H&E pixel size used when a sample's meta.json does not record one (bundles made before import existed)
HE_UM_PER_PX = 0.69


def ensure_dirs():
    for d in (SAMPLES, SESSIONS, USER_MARKERS):
        d.mkdir(parents=True, exist_ok=True)


def device():
    """'cuda' (NVIDIA), 'mps' (Apple silicon GPU) or 'cpu'."""
    import torch
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return "mps"
    return "cpu"
