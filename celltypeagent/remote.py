"""Working on a remote server over SSH: the data and every computation stay on the server; this computer shows the
viewer and runs Claude Code.

On the server (started by the app over SSH):
    celltypeagent remote-serve [--check]
starts the viewer in the background on the server's own localhost (or finds the one already running with the same
code), records it in <data folder>/remote_viewer.json and prints {"port", "build", "home", ...} as JSON.

On this computer (the desktop app): Connection runs remote-serve through `ssh`, forwards a local port to the
server's viewer (`ssh -N -L`), restarts the tunnel when it drops, and gives Claude Code an MCP server that runs on
the server through `ssh` (stdio passes through it). SSH settings come from ~/.ssh/config; the login must work without
a password prompt (a key or ssh-agent), since there is no terminal to type it in.
"""
import json, os, platform, re, shlex, signal, subprocess, sys, threading, time, urllib.request
from pathlib import Path
from . import config

STATE = "remote_viewer.json"


# ---------------- on the server ----------------
def _answers(port, timeout=2):
    try:
        urllib.request.urlopen(f"http://127.0.0.1:{port}/api/config", timeout=timeout)
        return True
    except Exception:
        return False


def serve(check=False):
    """Start (or reuse) the background viewer of this data folder; print where it listens as JSON."""
    info = {"build": config.build(), "home": str(config.HOME), "python": platform.python_version(),
            "system": f"{platform.system()} {platform.machine()}", "host": platform.node()}
    if check:
        print(json.dumps(info))
        return info
    f = config.HOME / STATE
    st = json.loads(f.read_text()) if f.exists() else {}
    alive = bool(st.get("port")) and _answers(st["port"])
    if alive and st.get("build") == info["build"]:
        info.update(port=st["port"], reused=True)
    else:
        if alive:  # the viewer of an older version of the code: replace it
            try:
                os.kill(int(st["pid"]), signal.SIGTERM)
                time.sleep(1)
            except (OSError, ValueError, KeyError):
                pass
        from .desktop import free_port
        port = free_port(8790)
        env = dict(os.environ, CELLTYPEAGENT_DESKTOP="1", CELLTYPEAGENT_REMOTE="1")
        log = open(config.HOME / "remote_viewer.log", "a")
        p = subprocess.Popen(config.self_command("viewer", "--port", str(port)), env=env, stdin=subprocess.DEVNULL,
                             stdout=log, stderr=log, start_new_session=True)  # outlives the SSH session
        for _ in range(600):
            if _answers(port, 1):
                break
            if p.poll() is not None:
                raise SystemExit(f"the viewer did not start; see {config.HOME / 'remote_viewer.log'}")
            time.sleep(0.2)
        f.write_text(json.dumps({"port": port, "pid": p.pid, "build": info["build"], "started": time.time()}))
        info.update(port=port, reused=False)
    print(json.dumps(info))
    return info


# ---------------- on this computer ----------------
def _ssh(*args):
    """ssh without prompts; CELLTYPEAGENT_SSH_CONFIG points at another ssh config (tests)."""
    extra = ["-F", os.environ["CELLTYPEAGENT_SSH_CONFIG"]] if os.environ.get("CELLTYPEAGENT_SSH_CONFIG") else []
    return ["ssh", *extra, "-o", "BatchMode=yes", "-o", "ConnectTimeout=15", "-o", "ServerAliveInterval=30", *args]


def ssh_hosts():
    """Host names of ~/.ssh/config (no patterns), for the connection dialog."""
    f = Path(os.environ.get("CELLTYPEAGENT_SSH_CONFIG") or Path.home() / ".ssh" / "config")
    if not f.exists():
        return []
    out = []
    for line in f.read_text(errors="replace").splitlines():
        m = re.match(r"\s*Host\s+(.+)", line, re.I)
        if m:
            out += [h for h in m.group(1).split() if not any(c in h for c in "*?!")]
    return sorted(set(out))


def _path(p):
    """A path for the remote shell: quoted, with ~ still expanded there."""
    p = p.strip()
    if p == "~":
        return '"$HOME"'
    if p.startswith("~/"):
        return '"$HOME"/' + shlex.quote(p[2:])
    return shlex.quote(p)


def remote_command(conn, *args):
    """The shell command that runs a subcommand of STAnnotator on the server, with its data folder."""
    env = f"CELLTYPEAGENT_HOME={_path(conn['home'])} " if conn.get("home", "").strip() else ""
    return f"env {env}{_path(conn.get('command') or 'celltypeagent')} " + " ".join(shlex.quote(a) for a in args)


def _explain(err, conn):
    e = err.strip()
    if "Permission denied" in e or "publickey" in e:
        return (f"SSH could not log in to {conn['host']} without a password. Set up a key (ssh-keygen, then "
                f"ssh-copy-id {conn['host']}) or add it to ssh-agent, then try again. ({e.splitlines()[-1]})")
    if "Could not resolve" in e or "Name or service not known" in e:
        return f"Unknown server {conn['host']!r}: use a name from ~/.ssh/config or user@address. ({e.splitlines()[-1]})"
    if "timed out" in e or "Connection refused" in e or "No route" in e:
        return f"Cannot reach {conn['host']}: {e.splitlines()[-1]}"
    if "not found" in e or "No such file" in e:
        return (f"STAnnotator was not found on the server at {conn.get('command') or 'celltypeagent'}. Give the full "
                f"path of its 'celltypeagent' command (e.g. /path/to/STAnnotator/.venv/bin/celltypeagent). ({e.splitlines()[-1]})")
    return e.splitlines()[-1] if e else "unknown error"


def run_remote(conn, *args, timeout=180):
    """Run a subcommand on the server; returns the JSON it printed last."""
    r = subprocess.run(_ssh(conn["host"], remote_command(conn, *args)), capture_output=True, text=True,
                       timeout=timeout, stdin=subprocess.DEVNULL)
    lines = [l for l in r.stdout.splitlines() if l.strip().startswith("{")]
    if r.returncode != 0 or not lines:
        raise RuntimeError(_explain(r.stderr or r.stdout, conn))
    return json.loads(lines[-1])


def check(conn):
    """Log in and ask STAnnotator on the server for its version and data folder (nothing is started)."""
    info = run_remote(conn, "remote-serve", "--check", timeout=60)
    info["local_build"] = config.build()
    return info


def mcp_entry(conn):
    """Claude Code's MCP server for this connection: STAnnotator's tools, run on the server through ssh."""
    cmd = _ssh(conn["host"], remote_command(conn, "mcp"))
    return {"command": cmd[0], "args": cmd[1:]}


class Connection:
    """The viewer on a server, reached through an SSH tunnel to a local port."""

    def __init__(self, conn):
        self.conn = dict(conn)
        self.tunnel = None
        self.port = None
        self.remote_port = None
        self.info = {}
        self.closed = False

    def start(self, progress=lambda m: None):
        progress(f"starting STAnnotator on {self.conn['host']}…")
        self.info = run_remote(self.conn, "remote-serve")
        self.remote_port = self.info["port"]
        from .desktop import free_port
        self.port = free_port(8800)
        progress("connecting…")
        self._open_tunnel()
        threading.Thread(target=self._watch, daemon=True).start()
        self.info["local_build"] = config.build()
        return self

    def _open_tunnel(self):
        self.tunnel = subprocess.Popen(
            _ssh("-N", "-o", "ExitOnForwardFailure=yes", "-L", f"127.0.0.1:{self.port}:127.0.0.1:{self.remote_port}",
                 self.conn["host"]), stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        for _ in range(150):
            if _answers(self.port, 1):
                return
            if self.tunnel.poll() is not None:
                raise RuntimeError(_explain(self.tunnel.stderr.read().decode(errors="replace"), self.conn))
            time.sleep(0.2)
        raise RuntimeError(f"no answer from the viewer on {self.conn['host']} through the SSH tunnel")

    def _watch(self):
        """Reopen the tunnel when it drops (network change, sleep); the server's viewer keeps running."""
        while not self.closed:
            time.sleep(3)
            if self.closed or self.tunnel.poll() is None:
                continue
            try:
                self._open_tunnel()
            except Exception:
                time.sleep(10)

    @property
    def url(self):
        return f"http://127.0.0.1:{self.port}/"

    def download(self, path, dest):
        """Save a file the server's viewer serves (an export) to a path on this computer."""
        with urllib.request.urlopen(self.url.rstrip("/") + path, timeout=3600) as r, open(dest, "wb") as f:
            n = 0
            while chunk := r.read(1 << 20):
                f.write(chunk)
                n += len(chunk)
        return n

    def stop(self):
        self.closed = True
        if self.tunnel and self.tunnel.poll() is None:
            self.tunnel.terminate()


# ---------------- remembered settings (this computer) ----------------
def settings_file():
    return config.HOME / "remote.json"


def load_settings():
    try:
        return json.loads(settings_file().read_text())
    except (OSError, ValueError):
        return {}


def save_settings(**kw):
    s = load_settings()
    s.update(kw)
    config.HOME.mkdir(parents=True, exist_ok=True)
    settings_file().write_text(json.dumps(s, indent=1))
    return s
