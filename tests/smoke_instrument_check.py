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
    # a pending check is dropped by any new placement (e.g. an H&E aligned by the agent), and Steps shows none for an H&E
    importer.mark_check("s1")
    assert reg().get("check_pending")
    importer.place_cells("s1", [[1 / PS, 0, 0], [0, 1 / PS, 0]], own, "H&E", by="agent", progress=lambda *a: None)
    webapp._drop_view("s1")
    assert not reg().get("check_pending") and "check_pending" not in json.load(open(mf)), reg()
    # a stain of the instrument that is not the DAPI (boundary): not checked on nuclei
    stain = os.path.join(out, "morphology_focus", "ch0002_boundary.ome.tif")
    shutil.copy(own, stain)
    r = c.post("/api/sample/s1/image", json={"path": stain, "kind": "IF"})
    assert r.status_code == 200 and not json.load(open(mf)).get("check_pending"), (r.text, json.load(open(mf)).get("check_pending"))
    # two DAPI in one output folder: told apart by their file in the image list (Image box, Discard)
    names = [x["name"] for x in importer.images("s1")]
    assert len(names) == len(set(names)) and any("ch0001_dapi" in n for n in names), names
    ims = c.get("/api/sample/s1/he_info").json()["images"]
    assert all(x.get("path") for x in ims), ims
    # a fluorescence file of several channels made from the instrument's images (its size, elsewhere), imported as an
    # "H&E" with no registration: placed by the cells' coordinates as a DAPI (not aligned as an H&E), then checked
    multi = os.path.join(home, "elsewhere", "four_channels.ome.tif")
    os.makedirs(os.path.dirname(multi))
    base = tifffile.imread(own)
    with tifffile.TiffWriter(multi, ome=True) as t:
        chans = np.stack([base, base // 2])
        t.write(chans, tile=(256, 256), subifds=2, photometric="minisblack", metadata={"axes": "CYX", "Channel": {"Name": ["DAPI", "CD45"]}})
        for k in (2, 4):
            t.write(chans[:, ::k, ::k], tile=(256, 256), subfiletype=1, photometric="minisblack")
    mt = c.get("/api/sample/s1/images/match", params={"path": multi}).json()
    assert mt["instrument"] and mt["kind"] == "DAPI", mt
    r = c.post("/api/sample/s1/image", json={"path": multi, "kind": "H&E", "fresh": True})
    assert r.status_code == 200, r.text
    m = json.load(open(mf))
    assert m["image_kind"] == "DAPI" and not m.get("image_pending") and m.get("check_pending") and os.path.realpath(m["he_path"]) == os.path.realpath(multi), \
        {k: m.get(k) for k in ("image_kind", "image_pending", "check_pending", "he_path")}
    st = check()
    print("multi-channel file of the instrument's size ->", st["result"]["verdict"], "read", st["result"].get("read_share"))
    assert st["result"]["verdict"] == "fine" and not st["result"].get("new_alignment"), st["result"]
    # the instrument's DAPI shown by the cells' coordinates (as a Xenium import leaves it, checked): Steps offers
    # Import (an H&E), Add and Discard, not Import alone
    m = json.load(open(mf)); m.update(he_path=own, image_kind="DAPI", registration_by=None); m.pop("check_pending", None)
    json.dump(m, open(mf, "w")); webapp._drop_view("s1")
    assert reg()["image"] == "instrument" and reg()["state"] == "optional", reg()
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
