"""A sample whose image cannot be read still opens: its cells on a blank background, and he_info says why
(image_error) instead of failing with 500. The image moved; cut short (a copy or sync still going: 4 bytes); not a
TIFF; no permission (macOS keeping the app out of Downloads): chmod where it takes the right to read away (not on
Windows, where it only makes the file read-only, nor for root), and everywhere a PermissionError raised by the opener.
Runs in a temporary data folder."""
import os, shutil, tempfile, json
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_img_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")  # (never the user's settings)
from stagent import importer, config
try:
    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    importer.import_sample("s1", h5, progress=lambda *a: None)
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    mf = config.SAMPLES / "s1" / "meta.json"
    from stagent import render
    import errno

    def denied(path):
        raise PermissionError(errno.EACCES, "Permission denied", path)
    real = render.Slide
    for case in ("missing", "cut short", "not a tiff", "no permission", "permission error"):
        img = os.path.join(home, f"{case.replace(' ', '_')}.tif")
        if case == "cut short":
            open(img, "wb").write(b"II*\x00")
        elif case == "not a tiff":
            open(img, "w").write("this is not an image\n" * 20)
        elif case == "no permission":
            open(img, "wb").write(b"II*\x00")
            os.chmod(img, 0)
            if os.access(img, os.R_OK):  # (Windows, root: still readable; the opener's PermissionError case covers it)
                print(case, "-> skipped (chmod does not take the right to read away here)")
                continue
        elif case == "permission error":
            open(img, "wb").write(b"II*\x00")
            render.Slide = denied
        m = json.load(open(mf))
        m.update(he_path=img, image_kind="H&E")
        json.dump(m, open(mf, "w"))
        webapp._views.pop("s1", None)  # (the app's Refresh does this)
        try:
            r = c.get("/api/sample/s1/he_info")
        finally:
            render.Slide = real
        assert r.status_code == 200, (case, r.status_code, r.text[:200])
        j = r.json()
        print(case, "->", j["image_kind"], "|", j["image_error"][:90])
        assert j["image_kind"] == "none" and img in j["image_error"], j
        if case == "missing":
            assert "no longer there" in j["image_error"]
        if case in ("cut short", "not a tiff"):
            assert "copying" in j["image_error"], j["image_error"]
        if case == "permission error":
            assert "Permission denied" in j["image_error"] or "macOS" in j["image_error"], j["image_error"]
        assert c.get("/api/sample/s1/cells.bin").status_code == 200, case  # (the cells, on the blank background)
    print("ok")
finally:
    for f in os.listdir(home):
        try:
            os.chmod(os.path.join(home, f), 0o700)
        except OSError:
            pass
    shutil.rmtree(home, ignore_errors=True)
