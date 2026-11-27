"""An update installing on Windows: the updater marks it before starting the installer and its script removes the mark
before opening the app; the app started meanwhile (packaging/entry.py) sees a fresh mark and waits, a stale one (an
installation that never ended) does not stop it; an installation missing its files is told. Runs in a temporary
folder (TMPDIR), no installer is started."""
import base64, importlib.util, json, os, sys, tempfile, time
from pathlib import Path
tmp = tempfile.mkdtemp(prefix="stagent_upd_")
os.environ["TMPDIR"] = tmp
tempfile.tempdir = None  # (gettempdir reads TMPDIR again)
ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("entry", ROOT / "packaging" / "entry.py")
entry = importlib.util.module_from_spec(spec)
spec.loader.exec_module(entry)
from stagent import updater
mark = Path(tempfile.gettempdir()) / "STAnnotator-update" / "installing.json"
assert not entry._update_installing()
# the updater on Windows: the mark, then the installer, then the mark removed, then the app
calls = []
updater._platform = lambda: "win"
updater.subprocess.Popen = lambda args, **kw: calls.append(args)
setup = Path(tmp) / "STAnnotator-Setup.exe"
setup.write_bytes(b"x")
r = updater.install(setup)
assert r["mode"] == "installer" and mark.exists(), r
ps = base64.b64decode(calls[0][-1]).decode("utf-16-le")
i_setup, i_rm, i_app = ps.index(str(setup)), ps.index("Remove-Item"), ps.rindex("Start-Process")
assert i_setup < i_rm < i_app and str(mark) in ps, ps
print("installer script:", ps[:60], "...", ps[i_rm:i_rm + 40])
assert entry._update_installing()  # the app started meanwhile waits
old = time.time() - 30 * 60
os.utime(mark, (old, old))
assert not entry._update_installing()  # an installation that never ended does not stop the app
mark.unlink()
# the installation's files: whole (this repository) or not (a folder without the page)
assert entry._installation_whole()
sys._MEIPASS = tmp
assert not entry._installation_whole()
del sys._MEIPASS
print("ok")
