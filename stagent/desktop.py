"""STAnnotator desktop app: the viewer in a native window, with its local server running inside the app.

    stagent app [--browser]          (or: python -m stagent.desktop)

The window is pywebview (WebKit on macOS); without pywebview, or with --browser, the viewer opens in the default
browser instead. Claude is used through the user's own Claude Code (MCP), not inside the app.
"""
import argparse, os, socket, sys, threading, time, webbrowser
from pathlib import Path
from . import config

APP_NAME = "STAnnotator"


def free_port(preferred=8765):
    for p in [preferred] + list(range(8766, 8800)):
        with socket.socket() as s:
            try:
                s.bind(("127.0.0.1", p))
                return p
            except OSError:
                continue
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def start_server(port):
    """Run the viewer's server in a background thread; return the uvicorn server once it accepts connections."""
    import uvicorn
    from .webapp import build_app
    server = uvicorn.Server(uvicorn.Config(build_app(port), host="127.0.0.1", port=port, log_level="warning"))
    threading.Thread(target=server.run, daemon=True).start()
    for _ in range(200):
        try:
            with socket.create_connection(("127.0.0.1", port), timeout=0.2):
                return server
        except OSError:
            time.sleep(0.1)
    raise RuntimeError("the app's local server did not start")



# The file filters of the Open dialogs. pywebview accepts "Description (*.ext;*.ext)" with a description of letters,
# digits, underscores and spaces only (webview.util.parse_file_type): any other character (e.g. "/") makes it refuse
# the dialog, so the Choose button does nothing (macOS and Windows alike).
FILE_TYPES = {"h5ad": ("AnnData or R (*.h5ad;*.rds;*.rda;*.RData)", "All files (*.*)"),
              "xenium": ("Xenium files (*.parquet;*.csv.gz;*.h5;*.xenium;*.mtx.gz;*.tsv.gz;*.ome.tif)", "All files (*.*)"),
              "xcells": ("Xenium cells (*.parquet;*.csv.gz)", "All files (*.*)"),
              "xcounts": ("Xenium counts (*.h5;*.mtx.gz;*.tsv.gz)", "All files (*.*)"),
              "xexp": ("Xenium experiment (*.xenium)", "All files (*.*)"),
              "ximage": ("Images (*.ome.tif;*.tif;*.tiff)", "All files (*.*)"),
              "image": ("Images (*.tif;*.tiff;*.btf)", "All files (*.*)"),
              "csv": ("CSV (*.csv)", "All files (*.*)"),
              "markers": ("Marker lists (*.csv;*.txt)", "All files (*.*)"),
              "bundle": ("Saved annotations (*.stannot;*.zip)", "All files (*.*)"),
              "flat": ("CosMx and MERSCOPE files (*.csv;*.csv.gz)", "All files (*.*)"),
              "geojson": ("GeoJSON (*.geojson;*.json)", "All files (*.*)"),
              "heonly": ("Slide images (*.tif;*.tiff;*.btf;*.svs;*.ndpi;*.png;*.jpg;*.jpeg)", "All files (*.*)")}
FILE_FILTER = r'^([\w ]+)\((\*(?:\.(?:\w+|\*))*(?:;\*(?:\.(?:\w+|\*))*)*)\)$'  # pywebview 6's, for the test


class Api:
    """Native functions for the page (window.pywebview.api): file dialogs, opening links, and working on a remote
    server (remote.py): the window then shows the server's viewer through an SSH tunnel."""

    def __init__(self, local_url=None):
        self._window = None  # private: pywebview exposes public attributes to the page
        self._local_url = local_url
        self._conn = None
        self._auto_tried = False

    def pick(self, kind, multiple=False):
        import webview
        fd = getattr(webview, "FileDialog", None)
        if kind == "folder":
            mode = fd.FOLDER if fd else webview.FOLDER_DIALOG
            res = self._window.create_file_dialog(mode)
        else:
            mode = fd.OPEN if fd else webview.OPEN_DIALOG
            types = FILE_TYPES.get(kind, ("All files (*.*)",))
            try:
                res = self._window.create_file_dialog(mode, allow_multiple=bool(multiple), file_types=types)
            except ValueError:  # a filter pywebview does not accept: the dialog without it rather than none
                res = self._window.create_file_dialog(mode, allow_multiple=bool(multiple))
        return [str(p) for p in res] if res else []

    def save_path(self, filename):
        """Native Save dialog; returns the chosen path or an empty string."""
        import webview
        fd = getattr(webview, "FileDialog", None)
        res = self._window.create_file_dialog(fd.SAVE if fd else webview.SAVE_DIALOG, save_filename=filename)
        if isinstance(res, (list, tuple)):
            res = res[0] if res else None
        return str(res) if res else ""

    def open_url(self, url):
        webbrowser.open(url)

    # ---- remote server ----
    def connection(self):
        """Where the app works now, the remembered server, the hosts of ~/.ssh/config; auto: connect to the server
        used last (once, when the app starts)."""
        from . import remote
        st = remote.load_settings()
        auto = st.get("mode") == "remote" and bool(st.get("last")) and not self._auto_tried and not self._conn
        self._auto_tried = True
        return {"connected": bool(self._conn), "host": self._conn.conn["host"] if self._conn else None,
                "info": self._conn.info if self._conn else None, "settings": st, "hosts": remote.ssh_hosts(),
                "auto": auto, "local_build": config.build()}

    def test_remote(self, conn):
        from . import remote
        try:
            return {"ok": True, **remote.check(conn)}
        except Exception as e:
            return {"error": str(e)}

    def connect_remote(self, conn):
        from . import remote
        try:
            c = remote.Connection(conn).start()
        except Exception as e:
            return {"error": str(e)}
        if self._conn:
            self._conn.stop()
        self._conn = c
        remote.save_settings(mode="remote", last={k: conn.get(k, "") for k in ("host", "command", "home")})
        self._window.load_url(c.url)
        return {"ok": True, **c.info}

    def use_local(self):
        from . import remote
        if self._conn:
            self._conn.stop()
            self._conn = None
        remote.save_settings(mode="local")
        self._window.load_url(self._local_url)
        return {"ok": True}

    def claude_status(self):
        from . import claudelink
        return claudelink.status()

    def open_claude(self, args):
        """Claude Code on this computer; with a server connection its tools run on the server."""
        from . import claudelink
        return claudelink.open_claude_code(**args, remote=self._conn.conn if self._conn else None)

    def download(self, path, dest):
        """Save an export of the server's viewer (a URL path) to a file on this computer."""
        if not self._conn:
            return {"error": "not connected to a server"}
        try:
            return {"saved": dest, "bytes": self._conn.download(path, dest)}
        except Exception as e:
            return {"error": str(e)}

    def report_add_local(self, path, days=7, parts=("conversations", "app_log")):
        """A problem report made on the server: add this computer's agent conversations and app.log to it."""
        from . import report
        try:
            return {"added": report.add_local(path, int(days), parts)}
        except Exception as e:
            return {"error": str(e)}

    def report_upload(self, path, name, days=7, parts=("conversations", "app_log"), note="", sender="", url=""):
        """With the app on a server: fetch the report made there, add this computer's conversations and app.log, and
        upload it to the developers."""
        import tempfile
        from . import report
        if not self._conn:
            return {"error": "not connected to a server"}
        tmp = Path(tempfile.mkdtemp()) / name
        try:
            self._conn.download(path, str(tmp))
            report.add_local(tmp, int(days), parts)
            return report.upload(tmp, name, {"version": config.build()}, note, sender, url)
        except Exception as e:
            return {"error": str(e)}
        finally:
            tmp.unlink(missing_ok=True)

    def write_text(self, path, text):
        try:
            Path(path).write_text(text)
            return {"saved": path}
        except OSError as e:
            return {"error": str(e)}

    def _close(self):
        if self._conn:
            self._conn.stop()


def _mac_identity():
    """Show the app's name (not "Python") in the macOS menu bar, best effort."""
    if sys.platform != "darwin":
        return
    try:
        from Foundation import NSBundle
        info = NSBundle.mainBundle().infoDictionary()
        info["CFBundleName"] = APP_NAME
        icon = config.PACKAGE / "data" / "app_icon.png"
        if icon.exists():
            from AppKit import NSApplication, NSImage
            NSApplication.sharedApplication().setApplicationIconImage_(NSImage.alloc().initWithContentsOfFile_(str(icon)))
    except Exception:
        pass


def main(argv=None):
    ap = argparse.ArgumentParser(prog="stagent app", description=__doc__.splitlines()[0])
    ap.add_argument("--browser", action="store_true", help="open the viewer in the default browser instead of a window")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args(argv)
    os.environ["STAGENT_DESKTOP"] = "1"  # the page then offers Claude Code instead of a chat box
    config.ensure_dirs()
    port = free_port(a.port)
    server = start_server(port)
    url = f"http://127.0.0.1:{port}/"
    try:
        import webview
    except ImportError:
        webview = None
        if not a.browser:
            print("pywebview is not installed (pip install pywebview); opening the browser instead")
    if a.browser or webview is None:
        webbrowser.open(url)
        print(f"{APP_NAME} is running at {url}  (Ctrl+C to quit)")
        try:
            while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            pass
    else:
        webview.settings["ALLOW_DOWNLOADS"] = True  # links that download files (the app also uses Save dialogs)
        api = Api(url)
        api._window = webview.create_window(APP_NAME, url, js_api=api, width=1500, height=950, min_size=(1000, 650))
        _mac_identity()
        webview.start()
        api._close()  # the SSH tunnel; the server's viewer keeps running for the next connection
        try:  # the Claude Code the app opened in Terminal quits with it
            from . import claudelink
            claudelink.close_claude_code()
        except Exception as e:
            print(f"could not close Claude Code: {e}")
    server.should_exit = True


if __name__ == "__main__":
    main()
