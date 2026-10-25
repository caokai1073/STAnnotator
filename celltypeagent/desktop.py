"""Spatial Annotator desktop app: the viewer in a native window, with its local server running inside the app.

    celltypeagent app [--browser]          (or: python -m celltypeagent.desktop)

The window is pywebview (WebKit on macOS); without pywebview, or with --browser, the viewer opens in the default
browser instead. Claude is used through the user's own Claude Code (MCP), not inside the app.
"""
import argparse, os, socket, sys, threading, time, webbrowser
from . import config

APP_NAME = "Spatial Annotator"


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


class Api:
    """Native functions for the page (window.pywebview.api): file dialogs and opening links."""

    def __init__(self):
        self._window = None  # private: pywebview exposes public attributes to the page

    def pick(self, kind, multiple=False):
        import webview
        fd = getattr(webview, "FileDialog", None)
        if kind == "folder":
            mode = fd.FOLDER if fd else webview.FOLDER_DIALOG
            res = self._window.create_file_dialog(mode)
        else:
            mode = fd.OPEN if fd else webview.OPEN_DIALOG
            types = {"h5ad": ("AnnData (*.h5ad)", "All files (*.*)"),
                     "image": ("Images (*.tif;*.tiff;*.btf)", "All files (*.*)"),
                     "csv": ("CSV (*.csv)", "All files (*.*)"),
                     "markers": ("Marker lists (*.csv;*.txt)", "All files (*.*)")}.get(kind, ("All files (*.*)",))
            res = self._window.create_file_dialog(mode, allow_multiple=bool(multiple), file_types=types)
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
    ap = argparse.ArgumentParser(prog="celltypeagent app", description=__doc__.splitlines()[0])
    ap.add_argument("--browser", action="store_true", help="open the viewer in the default browser instead of a window")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args(argv)
    os.environ["CELLTYPEAGENT_DESKTOP"] = "1"  # the page then offers Claude Code instead of a chat box
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
        api = Api()
        api._window = webview.create_window(APP_NAME, url, js_api=api, width=1500, height=950, min_size=(1000, 650))
        _mac_identity()
        webview.start()
    server.should_exit = True


if __name__ == "__main__":
    main()
