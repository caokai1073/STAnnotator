"""An image moved or deleted after import (Steps → Image: Locate…): the step says it is no longer there and gives its
path; Locate with the same file at its new place points the sample to it (the placement, the images' ids, the hidden
channels and the registration check kept); a different image (the cells off it) or the same place is refused; an
added image moved is found again the same way. Runs in a temporary data folder."""
import json, os, shutil, tempfile
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_relink_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")
import pandas as pd, tifffile
from stagent import importer, config, regcheck, render
try:
    importer.import_sample("s1", os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    W, H = 1100, 1500

    def pyramid(path, w=W, h=H, seed=0):
        rgb = np.random.default_rng(seed).integers(0, 255, (h, w, 3), dtype=np.uint8)
        with tifffile.TiffWriter(path) as t:
            t.write(rgb, tile=(128, 128), subifds=2, photometric="rgb", compression="zlib", metadata={"axes": "YXS"})
            for k in (2, 4):
                t.write(rgb[::k, ::k], tile=(128, 128), subfiletype=1, photometric="rgb", compression="zlib")
    os.makedirs(os.path.join(home, "old"))
    img = os.path.join(home, "old", "he.ome.tif")
    pyramid(img)
    xy = pd.read_parquet(config.SAMPLES / "s1" / "cells.parquet", columns=["x", "y"]).to_numpy()
    lo, hi = xy.min(0), xy.max(0)
    sx, sy = (W - 2) / (hi[0] - lo[0]), (H - 2) / (hi[1] - lo[1])
    A = [[sx, 0.0, -lo[0] * sx], [0.0, sy, -lo[1] * sy]]
    mf = config.SAMPLES / "s1" / "meta.json"
    m = json.load(open(mf))
    m.update(he_path=img, image_kind="H&E", placement_affine=A, registration_affine=A, registration_by="file")
    json.dump(m, open(mf, "w"))
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    reg = lambda: next(x for x in c.get("/api/sample/s1/steps").json()["steps"] if x["key"] == "registration")
    assert reg().get("image") == "placed", reg()
    # an added image (same pixels), one of its channels hidden
    img2 = os.path.join(home, "old", "stain.ome.tif")
    pyramid(img2, seed=1)
    r = c.post("/api/sample/s1/images/add", json={"path": img2, "kind": "IF", "same_as": "main"})
    assert r.status_code == 200, r.text
    added = r.json()["added"]
    m = json.load(open(mf)); m["images_hidden"] = [added + "#c1"]; json.dump(m, open(mf, "w"))
    # its registration was checked (the app's check: Steps → Image → Check)
    from stagent import session as sess
    import time
    pid = sess.placeholder("s1") or sess.Session.create("s1", "none", method="none").id

    def check():
        j = c.post("/api/sample/s1/registration/check", json={"session": pid}).json()
        assert "job" in j or "error" not in j, j
        for _ in range(900):
            st = c.get(f"/api/job/regcheck:s1").json()
            if st.get("done"):
                return st
            time.sleep(0.2)
        raise AssertionError("the check did not end")
    st = check()
    assert not st.get("error"), st["error"]
    assert json.load(open(regcheck._file("s1")))["image"] == img
    c.get("/api/sample/s1/he_info")  # (the view opens the image: its pixel size is recorded)
    assert importer.image_sizes("s1")[os.path.realpath(img)] == [W, H]
    # the folder is moved away (the app's open files closed first: Windows cannot move an open file)
    webapp._drop_view("s1"); render.close(img); render.close(img2)
    shutil.move(os.path.join(home, "old"), os.path.join(home, "new"))
    x = reg()
    print("moved ->", x["state"], "|", x["detail"])
    assert x["image"] == "error" and x["path"] == img and "no longer there" in x["detail"] and "Locate" in x["detail"], x
    # refused: the same place, another image (the cells off it), not an image of the sample
    small = os.path.join(home, "small.ome.tif")
    pyramid(small, 300, 300)
    other = os.path.join(home, "other.ome.tif")  # another slide: larger, the cells fall on it, same pixel size
    pyramid(other, W + 64, H + 64, seed=3)
    for old, new, why in ((img, img, "same file"), (img, small, "not the same image"), (img, other, f"{W + 64} x {H + 64} pixels"),
                          ("/nope.tif", small, "not an image of this sample")):
        r = c.post("/api/sample/s1/image_relink", json={"old": old, "path": new})
        assert r.status_code == 400 and why in r.json()["error"], (why, r.status_code, r.text)
    new_img, new_img2 = os.path.join(home, "new", "he.ome.tif"), os.path.join(home, "new", "stain.ome.tif")
    r = c.post("/api/sample/s1/image_relink", json={"old": img, "path": new_img})
    assert r.status_code == 200 and r.json()["inside"] >= 0.9, r.text
    m = json.load(open(mf))
    assert m["he_path"] == str(os.path.realpath(new_img)) and m["placement_affine"] == A and m["image_relinks"][-1]["from"] == img, m
    assert reg().get("image") == "placed", reg()
    j = c.get("/api/sample/s1/he_info").json()
    assert j["image_kind"] == "H&E" and not j.get("image_error"), j
    w = json.load(open(regcheck._file("s1")))
    assert w["key"] == regcheck._key("s1", m["he_path"]) and w["image"] == m["he_path"], (w["key"], w["image"])  # (carried over)
    assert img not in json.dumps(m["registration_checked"]), m["registration_checked"]
    st = check()  # the check again at the new place (it reads the image of its windows again)
    assert not st.get("error"), st["error"]
    print("check after Locate:", st["result"]["verdict"])
    # the added image, moved too, its size not recorded (imported before it was): the page is asked to confirm
    sf = config.SAMPLES / "s1" / "image_sizes.json"
    d = json.load(open(sf)); d.pop(os.path.realpath(img2), None); json.dump(d, open(sf, "w"))
    r = c.post("/api/sample/s1/image_relink", json={"old": img2, "path": new_img2})
    assert r.status_code == 409 and r.json().get("confirm") and "not recorded" in r.json()["error"], r.text
    r = c.post("/api/sample/s1/image_relink", json={"old": img2, "path": new_img2, "unknown_size_ok": True})
    assert r.status_code == 200, r.text
    m = json.load(open(mf))
    e = next(x for x in m["images"] if x["id"] == added)
    assert e["path"] == str(os.path.realpath(new_img2)) and m["images_hidden"] == [added + "#c1"], (e, m["images_hidden"])
    ids = [x["id"] for x in importer.images("s1")]
    assert ids[:1] == ["main"] and added in ids, ids
    # back to the first place and again: each relink keeps what it said (history not rewritten)
    webapp._drop_view("s1"); render.close(new_img); render.close(new_img2)
    shutil.move(os.path.join(home, "new"), os.path.join(home, "old"))
    r = c.post("/api/sample/s1/image_relink", json={"old": json.load(open(mf))["he_path"], "path": img})
    assert r.status_code == 200, r.text
    hist = json.load(open(mf))["image_relinks"]
    assert hist[0]["from"] == img and hist[0]["to"].endswith("new/he.ome.tif") and hist[-1]["to"] == str(os.path.realpath(img)), hist
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
