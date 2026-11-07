"""Updates of the desktop app: a newer release on GitHub is offered when the app starts (and from About), downloaded
in the background and installed with one click.

- check(): the latest release (GitHub's API, at most every few hours; the answer is kept in the settings), its
  version, notes and the installer for this system (STAnnotator.dmg on macOS, STAnnotator-Setup.exe on Windows).
- Download: the installer goes to a folder of the app's own (not the data folder), with progress.
- install(): macOS: a small script waits for the app to quit, puts the new STAnnotator.app from the DMG in the place
  of this one and opens it; where the app's folder cannot be written, the DMG opens in Finder to drag it over.
  Windows: the installer runs silently into the same place (Inno Setup keeps the folder and the per-user or
  all-users choice of the installation it updates; it may ask for permission), then the app starts again.
The data folder and the settings are not touched. A copy run from source is told of the release only.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

from . import __version__, config

REPO = "caokai1073/STAnnotator"
API = f"https://api.github.com/repos/{REPO}/releases/latest"
PAGE = f"https://github.com/{REPO}/releases/latest"
CHECK_EVERY = 6 * 3600


def _num(v):
    return tuple(int(x) for x in re.findall(r"\d+", str(v))[:4])


def newer(a, b):
    """Whether version a is newer than version b ("v0.2.47" > "0.2.46")."""
    return _num(a) > _num(b)


def _platform():
    return "mac" if sys.platform == "darwin" else "win" if os.name == "nt" else "other"


def can_install():
    """The built app on macOS or Windows can replace itself; a copy run from source cannot."""
    return bool(config.FROZEN) and _platform() in ("mac", "win")


def _fetch_latest(timeout=15):
    import urllib.request
    req = urllib.request.Request(API, headers={"Accept": "application/vnd.github+json", "User-Agent": "STAnnotator-app"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        rel = json.loads(r.read().decode("utf-8"))
    want = {"mac": ".dmg", "win": ".exe"}.get(_platform())
    asset = next((a for a in rel.get("assets", []) if want and a.get("name", "").lower().endswith(want)), None)
    return {"version": str(rel.get("tag_name", "")).lstrip("v"), "notes": rel.get("body") or "",
            "page": rel.get("html_url") or PAGE, "published": rel.get("published_at"),
            "asset": {"name": asset["name"], "url": asset["browser_download_url"], "bytes": asset.get("size")} if asset else None}


def check(force=False):
    """{current, latest, available, can_install, skipped, ...}; asks GitHub at most every CHECK_EVERY seconds."""
    st = config.read_settings()
    cached = st.get("update_latest") or {}
    if force or not cached or time.time() - st.get("update_checked", 0) > CHECK_EVERY:
        try:
            cached = _fetch_latest()
            config.write_settings(update_latest=cached, update_checked=time.time())
        except Exception as e:
            if not cached:
                return {"current": __version__, "available": False, "error": f"could not reach GitHub: {e}"}
    latest = cached.get("version", "")
    available = bool(latest) and newer(latest, __version__)
    return {"current": __version__, "latest": latest, "available": available, "platform": _platform(),
            "skipped": st.get("update_skip") == latest, "can_install": can_install() and bool(cached.get("asset")),
            "notes": cached.get("notes", ""), "page": cached.get("page", PAGE), "asset": cached.get("asset")}


def skip(version):
    config.write_settings(update_skip=version)


# ---------------- download ----------------
def _folder():
    d = Path(tempfile.gettempdir()) / "STAnnotator-update"
    d.mkdir(parents=True, exist_ok=True)
    return d


class Download:
    """The installer of a release, fetched in a thread: state(), then path when done."""

    def __init__(self, info):
        self.info, self.done, self.total, self.error, self.path = info, 0, info["asset"].get("bytes") or 0, None, None
        self.stop = False
        threading.Thread(target=self._run, daemon=True).start()

    def _run(self):
        import urllib.request
        a = self.info["asset"]
        dest = _folder() / a["name"]
        tmp = dest.with_name(dest.name + ".part")
        try:
            req = urllib.request.Request(a["url"], headers={"User-Agent": "STAnnotator-app"})
            with urllib.request.urlopen(req, timeout=60) as r, open(tmp, "wb") as f:
                self.total = int(r.headers.get("Content-Length") or self.total or 0)
                while not self.stop:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
                    self.done += len(chunk)
            if self.stop:
                tmp.unlink(missing_ok=True)
                self.error = "stopped"
                return
            if self.total and self.done != self.total:
                raise OSError(f"the download stopped at {self.done} of {self.total} bytes")
            os.replace(tmp, dest)
            self.path = dest
        except Exception as e:
            tmp.unlink(missing_ok=True)
            self.error = str(e)

    def state(self):
        return {"version": self.info["version"], "done": self.done, "total": self.total, "error": self.error,
                "ready": bool(self.path)}


# ---------------- install ----------------
def _mac_bundle():
    """The STAnnotator.app this process runs from."""
    p = Path(sys.executable).resolve()
    for q in p.parents:
        if q.suffix == ".app":
            return q
    raise RuntimeError("this app is not running from an .app bundle")


MAC_SCRIPT = r'''#!/bin/bash
# STAnnotator update: wait for the app to quit, replace it with the one in the DMG, open it.
PID="$1"; DMG="$2"; APP="$3"; LOG="$4"
exec >>"$LOG" 2>&1
echo "== update $(date)"
case "$APP" in /*.app) ;; *) echo "not an app bundle: $APP"; open "$DMG"; exit 1;; esac
for i in $(seq 1 120); do kill -0 "$PID" 2>/dev/null || break; sleep 0.5; done
MNT=$(mktemp -d /tmp/stannotator-update.XXXXXX)
if ! hdiutil attach -nobrowse -noautoopen -mountpoint "$MNT" "$DMG"; then open "$DMG"; exit 1; fi
NEW="$MNT/STAnnotator.app"
if [ ! -d "$NEW" ]; then hdiutil detach "$MNT" -quiet; open "$DMG"; exit 1; fi
STAGE="$(dirname "$APP")/.STAnnotator-update.app"
rm -rf "$STAGE"
if ditto "$NEW" "$STAGE" && rm -rf "$APP" && mv "$STAGE" "$APP"; then
  hdiutil detach "$MNT" -quiet; rmdir "$MNT" 2>/dev/null
  rm -f "$DMG"
  open "$APP"
  echo "updated $APP"
else
  rm -rf "$STAGE"
  hdiutil detach "$MNT" -quiet
  open "$DMG"   # could not write there: the user drags the app over
  echo "could not replace $APP; opened the DMG"
fi
'''


def install(path):
    """Start the installation of a downloaded installer; the caller then quits the app. Returns what will happen."""
    path = Path(path)
    log = _folder() / "update.log"
    if _platform() == "mac":
        app = _mac_bundle()
        if not os.access(app.parent, os.W_OK):
            subprocess.Popen(["open", str(path)])
            return {"mode": "finder", "note": f"{app.parent} cannot be written: drag STAnnotator to Applications "
                                              "in the window that opens, replacing the old one"}
        script = _folder() / "install_update.sh"
        script.write_text(MAC_SCRIPT)
        script.chmod(0o755)
        subprocess.Popen(["/bin/bash", str(script), str(os.getpid()), str(path), str(app), str(log)],
                         start_new_session=True, stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
        return {"mode": "replace", "app": str(app)}
    if _platform() == "win":
        import base64
        exe = Path(sys.executable).resolve()
        q = lambda x: "'" + str(x).replace("'", "''") + "'"
        # /SILENT shows only the progress; Inno Setup keeps the folder and the privileges of the installation it
        # updates and closes this app (CloseApplications); then the app starts again. PowerShell gets the command
        # as UTF-16, so folders with non-ASCII names (a Chinese user name) work, which a .cmd file in the code page
        # would not.
        ps = (f"Start-Process -FilePath {q(path)} -Wait -ArgumentList "
              f"'/SILENT','/SP-','/SUPPRESSMSGBOXES','/NOCANCEL','/CLOSEAPPLICATIONS',{q(f'/LOG={log}')}; "
              f"Start-Process -FilePath {q(exe)}")
        flags = getattr(subprocess, "CREATE_NO_WINDOW", 0) | getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        subprocess.Popen(["powershell", "-NoProfile", "-ExecutionPolicy", "Bypass", "-EncodedCommand",
                          base64.b64encode(ps.encode("utf-16-le")).decode()], creationflags=flags, close_fds=True)
        return {"mode": "installer", "app": str(exe)}
    raise RuntimeError("updates install themselves only in the macOS and Windows apps")


def cleanup():
    """Installers of earlier updates (kept until the next start)."""
    d = Path(tempfile.gettempdir()) / "STAnnotator-update"
    for f in d.glob("*") if d.is_dir() else []:
        if f.name != "update.log" and time.time() - f.stat().st_mtime > 3600:
            shutil.rmtree(f, ignore_errors=True) if f.is_dir() else f.unlink(missing_ok=True)
