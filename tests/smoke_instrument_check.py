"""The instrument's own DAPI put under the cells by their coordinates (Steps → Image, a file of the Xenium output folder)
is checked once against its nuclei: Steps says so (check_pending) and the check, when run, clears it; a DAPI of these
cells reads "fine", another run's DAPI (the same size, its nuclei elsewhere) is far off. Synthetic data: the cells of
the MERFISH example as a Xenium source, nuclei drawn at their places. Runs in a temporary data folder."""
import json, os, shutil, tempfile, time
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_inst_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")
import pandas as pd, tifffile
from stagent import importer, config
try:
    importer.import_sample("s1", os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    cf = config.SAMPLES / "s1" / "cells.parquet"
    cells = pd.read_parquet(cf)
    cells["x"] -= cells.x.min() - 20  # (a Xenium's coordinates start at 0)
    cells["y"] -= cells.y.min() - 20
    cells.to_parquet(cf)
    PS = 1.0  # um per pixel of the synthetic DAPI
    W, H = int(cells.x.max() / PS) + 40, int(cells.y.max() / PS) + 40
    out = os.path.join(home, "xenium_out")
    os.makedirs(os.path.join(out, "morphology_focus"))

    def dapi(path, x, y):  # nuclei (a small bright disk) at x, y (um)
        img = np.zeros((H, W), np.uint8)
        yy, xx = np.mgrid[-3:4, -3:4]
        disk = (xx ** 2 + yy ** 2 <= 9)
        px, py = np.round(x / PS).astype(int), np.round(y / PS).astype(int)
        for dy, dx in zip(*np.nonzero(disk)):
            ok = (py + dy - 3 >= 0) & (py + dy - 3 < H) & (px + dx - 3 >= 0) & (px + dx - 3 < W)
            img[py[ok] + dy - 3, px[ok] + dx - 3] = np.random.default_rng(1).integers(150, 250, ok.sum())
        with tifffile.TiffWriter(path) as t:
            t.write(img, tile=(256, 256), subifds=3, compression="zlib", metadata={"axes": "YX"})
            for k in (2, 4, 8):
                t.write(img[::k, ::k], tile=(256, 256), subfiletype=1, compression="zlib")
    x, y = cells.x.to_numpy(), cells.y.to_numpy()
    own = os.path.join(out, "morphology_focus", "ch0000_dapi.ome.tif")
    dapi(own, x, y)
    other_dir = os.path.join(home, "other_run", "morphology_focus")
    os.makedirs(other_dir)
    other = os.path.join(other_dir, "ch0000_dapi.ome.tif")
    rng = np.random.default_rng(0)  # another run's tissue: nuclei elsewhere, the same image size
    dapi(other, rng.uniform(20, x.max(), len(x)), rng.uniform(20, y.max(), len(y)))
    mf = config.SAMPLES / "s1" / "meta.json"
    m = json.load(open(mf))
    m["source"] = {**(m.get("source") or {}), "kind": "xenium", "pixel_size_um": PS, "folder": out}
    json.dump(m, open(mf, "w"))
    from starlette.testclient import TestClient
    from stagent import webapp, session as sess
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    pid = sess.placeholder("s1") or sess.Session.create("s1", "none", method="none").id
    reg = lambda: next(x for x in c.get("/api/sample/s1/steps").json()["steps"] if x["key"] == "registration")

    def check():
        r = c.post("/api/sample/s1/registration/check", json={"session": pid})
        assert r.status_code == 200, r.text
        for _ in range(1500):
            st = c.get("/api/job/regcheck:s1").json()
            if st.get("done"):
                return st
            time.sleep(0.2)
        raise AssertionError("the check did not end")
    for path, paired in ((own, True), (other, False)):
        if not paired:  # (another run's DAPI copied into this output folder, as a user might)
            dst = os.path.join(out, "morphology_focus", "ch0001_dapi.ome.tif")
            shutil.copy(path, dst)
            path = dst
        r = c.post("/api/sample/s1/image", json={"path": path, "kind": "DAPI"})
        assert r.status_code == 200, r.text
        webapp._drop_view("s1")
        x_ = reg()
        assert x_.get("check_pending") == "instrument" and x_["image"] == "placed", x_
        st = check()
        assert not st.get("error"), st["error"]
        res = st["result"]
        print("own DAPI" if paired else "another run's DAPI", "->", res["verdict"], "median", res.get("median_offset_um"),
              "read", res.get("read_share"), "align again:", res.get("new_alignment"))
        assert not reg().get("check_pending")  # (asked once)
        if paired:
            assert res["verdict"] == "fine" and not res.get("new_alignment"), res
        else:
            assert res.get("new_alignment") or res["verdict"] != "fine", res
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
