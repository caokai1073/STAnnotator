"""Connecting the app to the user's own Claude Code: a workspace whose MCP config points at this installation, opened
in Terminal. The app never signs in to Claude itself; each user talks to Claude in Claude Code with their own account.

The viewer and the MCP server (started by Claude) are separate processes. They share two small files in the data
folder: viewer.json (what the viewer shows: port, sample, session) and viewer_request.json (a session the agent
asks the viewer to show).
"""
import glob, json, os, platform, re, shutil, subprocess, time
from pathlib import Path
from . import config

# where Claude Code usually lives; a Dock-launched app has a minimal PATH
CLAUDE_PATHS = ["~/.local/bin", "~/.claude/local", "/opt/homebrew/bin", "/usr/local/bin", "~/.npm-global/bin",
                "~/.bun/bin", "~/.volta/bin"]


def server_entry():
    cmd = config.self_command("mcp")
    return {"command": cmd[0], "args": cmd[1:], "env": {"CELLTYPEAGENT_HOME": str(config.HOME)}}


_shell_claude = None


def claude_code_path():
    """Claude Code's executable: the app's PATH, the usual install folders (nvm included), and on macOS the user's
    interactive login shell, where nvm, conda or Homebrew set PATH up (~/.zshrc)."""
    extra = [str(Path(p).expanduser()) for p in CLAUDE_PATHS]
    extra += sorted(glob.glob(str(Path.home() / ".nvm" / "versions" / "node" / "*" / "bin")), reverse=True)
    hit = shutil.which("claude", path=os.pathsep.join([os.environ.get("PATH", "")] + extra))
    if hit or platform.system() != "Darwin":
        return hit
    global _shell_claude
    if _shell_claude is None:
        try:
            out = subprocess.run([os.environ.get("SHELL") or "/bin/zsh", "-ilc",
                                  "whence -p claude 2>/dev/null || type -P claude 2>/dev/null || command -v claude"],
                                 capture_output=True, text=True, timeout=15, stdin=subprocess.DEVNULL).stdout
            paths = [l.strip() for l in out.splitlines() if l.strip().startswith("/")]
            _shell_claude = paths[-1] if paths and os.access(paths[-1], os.X_OK) else ""
        except (OSError, subprocess.SubprocessError):
            _shell_claude = ""
    return _shell_claude or None


def mcp_seen():
    """When Claude last used the tools (the MCP server notes it), and whether that Claude Code is still running."""
    p = config.HOME / "mcp_seen.json"
    if not p.exists():
        return {"seen": False}
    d = json.load(open(p))
    try:
        os.kill(int(d["pid"]), 0)
        alive = True
    except (OSError, ValueError, KeyError):
        alive = False
    return {"seen": True, "last": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(d["time"])), "via": d.get("via"),
            "running": alive}


def status():
    """Whether Claude Code is installed (shown when the app's Claude Code button cannot open it)."""
    return {"platform": platform.system(), "claude_code": claude_code_path()}


# ---------------- Claude Code workspace ----------------
CLAUDE_MD = """# STAnnotator workspace

This folder connects Claude Code to the STAnnotator app (celltypeagent tools, data folder {home}).
Use only the celltypeagent tools. The user looks at the sample in the app; call use_viewer_session to load what
they are looking at before acting, and open_viewer after starting a new annotation (it shows it in the app). Follow the instructions of
the celltypeagent MCP server. Reply in the user's language.
"""


def workspace(remote=None):
    """A folder for Claude Code with this installation's MCP server and its tools allowed. remote: a server connection
    (remote.py); the tools then run on the server through ssh, with the server's data folder."""
    if remote:
        from . import remote as rm
        ws = config.HOME / f"claude_workspace_{re.sub(r'[^A-Za-z0-9_.-]+', '_', remote['host'])}"
        entry, home = rm.mcp_entry(remote), f"{remote.get('home') or '~/celltypeagent'} on {remote['host']}"
    else:
        ws, entry, home = config.HOME / "claude_workspace", server_entry(), config.HOME
    (ws / ".claude").mkdir(parents=True, exist_ok=True)
    (ws / "celltypeagent_mcp.json").write_text(json.dumps({"mcpServers": {"celltypeagent": entry}}, indent=2))
    (ws / ".claude" / "settings.json").write_text(json.dumps({"permissions": {"allow": ["mcp__celltypeagent"]}}, indent=2))
    (ws / "CLAUDE.md").write_text(CLAUDE_MD.format(home=home))
    return ws


def claude_code_command(prompt=None, remote=None):
    ws = workspace(remote)
    # the prompt goes first: --mcp-config takes several values and would swallow it
    cmd = f'cd "{ws}" && claude ' + (f"'{prompt}' " if prompt else "") + f'--mcp-config "{ws / "celltypeagent_mcp.json"}"'
    return cmd, ws


TAB_TITLE = "STAnnotator - Claude Code"
# Bring the Terminal tab that runs this app's Claude Code to the front (found by the tty recorded when it was opened,
# still running claude) and type the given message into it; otherwise open it in a new window (the message goes in as
# the first prompt). Prints "existing <tty>" or "new <tty>".
# (Claude Code changes the tab's title while it runs, so the title cannot identify the tab.)
FOCUS_SCRIPT = """on run argv
  set theTTY to item 1 of argv
  set theScript to item 2 of argv
  set theTitle to item 3 of argv
  set theText to item 4 of argv
  tell application "Terminal"
    if theTTY is not "" then
      repeat with w in windows
        repeat with t in tabs of w
          try
            if (tty of t) is theTTY and (busy of t) then
              set ps to processes of t
              if ps contains "claude" or ps contains "node" then
                set selected of t to true
                set index of w to 1
                activate
                if theText is not "" then do script theText in t
                return "existing " & theTTY
              end if
            end if
          end try
        end repeat
      end repeat
    end if
    set t to do script quoted form of theScript
    try
      set custom title of t to theTitle
    end try
    activate
    return "new " & (tty of t)
  end tell
end run
"""


def open_claude_code(sample=None, session=None, request=None, send_if_open=False, remote=None):
    """Open the user's own Claude Code in Terminal (macOS), connected to this app; elsewhere return the command.
    request: the first message (otherwise Claude loads the session and waits); with send_if_open it is also typed into
    a Claude Code window that is already open. On a server the app reaches over SSH, Claude Code runs on the user's
    computer: the request goes back to the app there (which calls this with remote=its connection)."""
    if remote is None and os.environ.get("CELLTYPEAGENT_REMOTE") == "1":
        return {"remote": True, "args": {"sample": sample, "session": session, "request": request,
                                         "send_if_open": send_if_open}}
    prompt = None
    if session:
        prompt = (f"I am looking at sample {sample} (session {session}) in the STAnnotator app. "
                  "Call use_viewer_session to load it, then "
                  + (f"do this: {request}" if request else "wait for my request."))
    elif request:
        prompt = request
    if prompt:
        prompt = prompt.replace("'", "\u2019")  # it is passed in single quotes to the shell
    send = prompt if (send_if_open and request) else None
    cmd, ws = claude_code_command(prompt, remote)
    if platform.system() != "Darwin":
        return {"opened": False, "command": cmd}
    script = ws / "start_claude.command"
    exe = claude_code_path()
    # the prompt goes first: --mcp-config takes several values and would swallow it
    args = (f"'{prompt}' " if prompt else "") + f'--mcp-config "{ws / "celltypeagent_mcp.json"}"'
    body = ["#!/bin/zsh -il", f'cd "{ws}"']  # interactive login shell: the same PATH as the user's Terminal
    if exe:
        body.append(f'exec "{exe}" {args}')
    else:
        body += ['if ! command -v claude >/dev/null 2>&1; then',
                 '  echo "Claude Code was not found. Install it (https://docs.claude.com/en/docs/claude-code/setup),"',
                 '  echo "log in once with your own Claude account, then press the button in the app again."',
                 '  read -r; exit 1',
                 'fi',
                 f'exec claude {args}']
    script.write_text("\n".join(body) + "\n")
    script.chmod(0o755)
    focus, tab = ws / "focus_claude.applescript", ws / "claude_tab.txt"
    focus.write_text(FOCUS_SCRIPT)
    tty = tab.read_text().strip() if tab.exists() else ""
    r = subprocess.run(["osascript", str(focus), tty, str(script), TAB_TITLE, send or ""], capture_output=True,
                       text=True, check=False)
    out = r.stdout.strip().split(" ", 1)
    if r.returncode != 0 or len(out) < 2:  # e.g. the user did not allow the app to control Terminal
        subprocess.run(["open", "-a", "Terminal", str(script)], check=False)
        return {"opened": True, "reused": False, "command": cmd, "claude_code": exe}
    tab.write_text(out[1])
    return {"opened": True, "reused": out[0] == "existing", "command": cmd, "claude_code": exe}


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
