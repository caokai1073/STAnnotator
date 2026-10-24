"""Connecting the app to the user's own Claude: Claude Desktop (MCP server entry) and Claude Code (a workspace whose
MCP config points at this installation, opened in Terminal). The app never signs in to Claude itself; each user
talks to Claude in Anthropic's apps with their own account.

The viewer and the MCP server (started by Claude) are separate processes. They share two small files in the data
folder: viewer.json (what the viewer shows: port, sample, session) and viewer_request.json (a session the agent
asks the viewer to show).
"""
import json, os, platform, shutil, subprocess, sys, time
from pathlib import Path
from . import config

# where Claude Code usually lives; a Dock-launched app has a minimal PATH
CLAUDE_PATHS = ["~/.local/bin", "~/.claude/local", "/opt/homebrew/bin", "/usr/local/bin", "~/.npm-global/bin",
                "~/.bun/bin", "~/.volta/bin"]


def desktop_config_path():
    system = platform.system()
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "Claude" / "claude_desktop_config.json"
    if system == "Windows":
        return Path(os.environ.get("APPDATA", Path.home())) / "Claude" / "claude_desktop_config.json"
    return Path.home() / ".config" / "Claude" / "claude_desktop_config.json"


def server_entry():
    return {"command": sys.executable, "args": ["-m", "celltypeagent.mcp_server"],
            "env": {"CELLTYPEAGENT_HOME": str(config.HOME)}}


def register_desktop(config_path=None, name="celltypeagent"):
    """Add this installation's MCP server to Claude Desktop (other entries are kept; a backup is written)."""
    path = Path(config_path) if config_path else desktop_config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(path.read_text()) if path.exists() and path.read_text().strip() else {}
    if path.exists():
        shutil.copy(path, path.with_suffix(f".backup-{time.strftime('%Y%m%d-%H%M%S')}.json"))
    cfg.setdefault("mcpServers", {})[name] = server_entry()
    path.write_text(json.dumps(cfg, indent=2))
    return path


def claude_code_path():
    extra = os.pathsep.join(str(Path(p).expanduser()) for p in CLAUDE_PATHS)
    return shutil.which("claude", path=os.environ.get("PATH", "") + os.pathsep + extra)


def status():
    """What is installed and connected (shown in the app's Claude dialog)."""
    path = desktop_config_path()
    entry = None
    try:
        entry = json.loads(path.read_text()).get("mcpServers", {}).get("celltypeagent") if path.exists() else None
    except (ValueError, OSError):
        pass
    apps = [Path("/Applications/Claude.app"), Path.home() / "Applications" / "Claude.app"]
    return {"platform": platform.system(),
            "desktop_installed": any(p.exists() for p in apps) if platform.system() == "Darwin" else None,
            "desktop_registered": bool(entry),
            "desktop_current": bool(entry) and entry == server_entry(),
            "desktop_config": str(path),
            "claude_code": claude_code_path()}


# ---------------- Claude Code workspace ----------------
CLAUDE_MD = """# Spatial Annotator workspace

This folder connects Claude Code to the Spatial Annotator app (celltypeagent tools, data folder {home}).
Use only the celltypeagent tools. The user looks at the sample in the app; call use_viewer_session to load what
they are looking at before acting, and open_viewer after starting a new annotation (it shows it in the app). Follow the instructions of
the celltypeagent MCP server. Reply in the user's language.
"""


def workspace():
    """A folder for Claude Code with this installation's MCP server and its tools allowed."""
    ws = config.HOME / "claude_workspace"
    (ws / ".claude").mkdir(parents=True, exist_ok=True)
    (ws / "celltypeagent_mcp.json").write_text(json.dumps({"mcpServers": {"celltypeagent": server_entry()}}, indent=2))
    (ws / ".claude" / "settings.json").write_text(json.dumps({"permissions": {"allow": ["mcp__celltypeagent"]}}, indent=2))
    (ws / "CLAUDE.md").write_text(CLAUDE_MD.format(home=config.HOME))
    return ws


def claude_code_command(prompt=None):
    ws = workspace()
    cmd = f'cd "{ws}" && claude --mcp-config "{ws / "celltypeagent_mcp.json"}"'
    return cmd + (f" '{prompt}'" if prompt else ""), ws


def open_claude_code(sample=None, session=None, request=None):
    """Open the user's own Claude Code in Terminal (macOS), connected to this app; elsewhere return the command.
    request: the user's first message (otherwise Claude loads the session and waits)."""
    prompt = None
    if session:
        prompt = (f"I am looking at sample {sample} (session {session}) in the Spatial Annotator app. "
                  "Call use_viewer_session to load it, then "
                  + (f"do this: {request}" if request else "wait for my request."))
    elif request:
        prompt = request
    if prompt:
        prompt = prompt.replace("'", "\u2019")  # it is passed in single quotes to the shell
    cmd, ws = claude_code_command(prompt)
    if platform.system() != "Darwin":
        return {"opened": False, "command": cmd}
    script = ws / "start_claude.command"
    exe = claude_code_path()
    body = ["#!/bin/zsh -l", f'cd "{ws}"']
    if exe:
        body.append(f'exec "{exe}" --mcp-config "{ws / "celltypeagent_mcp.json"}"' + (f" '{prompt}'" if prompt else ""))
    else:
        body += ['if ! command -v claude >/dev/null; then',
                 '  echo "Claude Code is not installed. Install it (https://claude.com/claude-code), log in with your',
                 '  echo "own Claude account, then press the button in the app again."; read; exit 1; fi',
                 f'exec claude --mcp-config "{ws / "celltypeagent_mcp.json"}"' + (f" '{prompt}'" if prompt else "")]
    script.write_text("\n".join(body) + "\n")
    script.chmod(0o755)
    subprocess.run(["open", "-a", "Terminal", str(script)], check=False)
    return {"opened": True, "command": cmd, "claude_code": exe}


# ---------------- viewer <-> agent ----------------
def _vfile(name):
    return config.HOME / name


def write_viewer_state(**kw):
    p = _vfile("viewer.json")
    try:
        cur = json.loads(p.read_text()) if p.exists() else {}
    except ValueError:
        cur = {}
    cur.update(kw, updated=time.time())
    p.write_text(json.dumps(cur))


def viewer_state():
    p = _vfile("viewer.json")
    try:
        return json.loads(p.read_text()) if p.exists() else {}
    except ValueError:
        return {}


def viewer_alive(state=None):
    """The viewer recorded in viewer.json answers on its port."""
    import urllib.request
    st = state or viewer_state()
    if not st.get("port"):
        return False
    try:
        urllib.request.urlopen(f"http://127.0.0.1:{st['port']}/api/config", timeout=2)
        return True
    except Exception:
        return False


def request_show(sample, session):
    _vfile("viewer_request.json").write_text(json.dumps({"sample": sample, "session": session, "time": time.time()}))


def viewer_request():
    p = _vfile("viewer_request.json")
    try:
        return json.loads(p.read_text()) if p.exists() else {}
    except ValueError:
        return {}
