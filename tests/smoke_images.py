"""Several images per sample (importer.images, add_image, remove_image; the viewer's warped tiles): another image of the
sample, rotated and mirrored against the one the cells are placed on and added with its own affine, is shown in the
frame's place matching the frame; an image added in the same pixels as another; the instrument's morphology images of
a Xenium source listed by themselves (DAPI and stains, named from their files); the image shown before stays listed
when the cells move onto another one; an added image removed. Runs in a temporary data folder."""
import io, json, os, shutil, tempfile
import numpy as np, tifffile
from PIL import Image
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_images_")
os.environ["STAGENT_HOME"] = home; os.environ["STAGENT_SETTINGS_DIR"] = home + "/cfg"
from stagent import importer, config, render, webapp


def pyramid(path, img):
    rgb = img.ndim == 3
    with tifffile.TiffWriter(path) as w:
        w.write(img, tile=(256, 256), subifds=2, photometric="rgb" if rgb else "minisblack", metadata={"axes": "YXS" if rgb else "YX"})
        for k in (2, 4):
            w.write(img[::k, ::k], tile=(256, 256), subfiletype=1, photometric="rgb" if rgb else "minisblack")


try:
    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    importer.import_sample("s1", h5, progress=lambda *a: None)
    import pandas as pd
    c = pd.read_parquet(config.SAMPLES / "s1" / "cells.parquet", columns=["x", "y"])
    x0, y0 = c.x.min(), c.y.min()
    W, H = int((c.x.max() - x0) / 2) + 200, int((c.y.max() - y0) / 2) + 200
    # the frame: a textured image, cells at 2 um/px; the other: the same tissue rotated 90 degrees and mirrored
    rng = np.random.default_rng(0)
    base = (np.kron(rng.integers(40, 255, (H // 16 + 1, W // 16 + 1, 3)), np.ones((16, 16, 1)))[:H, :W]).astype(np.uint8)
    pyramid(os.path.join(home, "frame.ome.tif"), base)
    other = np.ascontiguousarray(base.transpose(1, 0, 2))  # (x, y) -> (y, x): a 90 degree turn and a mirror
    pyramid(os.path.join(home, "other.ome.tif"), other)
    A = np.array([[0.5, 0, 100 - 0.5 * x0], [0, 0.5, 100 - 0.5 * y0]])
    importer.set_image("s1", os.path.join(home, "frame.ome.tif"))
    importer.place_cells("s1", A, os.path.join(home, "frame.ome.tif"), "H&E", by="file", progress=lambda *a: None)
    reg = os.path.join(home, "other_affine.csv")
    B = A[[1, 0]]  # the other image's pixels: (py, px)
    pd.DataFrame([dict(zip(["a11", "a12", "a21", "a22", "tx", "ty"], [B[0, 0], B[0, 1], B[1, 0], B[1, 1], B[0, 2], B[1, 2]]))]).to_csv(reg, index=False)
    r = importer.add_image("s1", os.path.join(home, "other.ome.tif"), "IF", registration=reg, name="turned")
    ims = importer.images("s1")
    assert [x["id"] for x in ims] == ["main", r["added"]] and ims[1]["name"] == "turned", ims
    # shown in the frame's place: the same picture as the frame
    v = webapp.view("s1")
    sl, T, _ = v.layer(r["added"])
    bx = (300, 200, 900, 700)
    a = np.asarray(render.slide(os.path.join(home, "frame.ome.tif")).read(*bx, 300, 250).convert("L"), float)
    b = np.asarray(render.read_warped(sl, T, *bx, 300, 250).convert("L"), float)
    rr = np.corrcoef(a.ravel(), b.ravel())[0, 1]
    assert rr > 0.95, f"the other image does not match the frame in its place (r = {rr:.2f})"
    # through the viewer's tiles
    from starlette.testclient import TestClient
    cl = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    info = cl.get("/api/sample/s1/he_info").json()
    assert [x["id"] for x in info["images"]] == ["main", r["added"]], info["images"]
    t1 = cl.get("/api/sample/s1/he/8/0_0.jpg"); t2 = cl.get("/api/sample/s1/he/8/0_0.jpg", params={"img": r["added"]})
    assert t1.status_code == t2.status_code == 200 and t1.content != t2.content
    # the same pixels as another image of the sample
    pyramid(os.path.join(home, "stain.ome.tif"), other[..., 0])
    r2 = cl.post("/api/sample/s1/images/add", json={"path": os.path.join(home, "stain.ome.tif"), "kind": "IF", "same_as": r["added"]}).json()
    assert next(x for x in importer.images("s1") if x["id"] == r2["added"])["affine"] == ims[1]["affine"], r2
    # the cells moved onto another image: the one shown before stays listed (with the affine that placed them on it)
    pyramid(os.path.join(home, "he2.ome.tif"), base)
    importer.place_cells("s1", A, os.path.join(home, "he2.ome.tif"), "H&E", by="file", progress=lambda *a: None)
    ims = importer.images("s1")
    old = [x for x in ims if x["path"].endswith("frame.ome.tif")]
    assert ims[0]["path"].endswith("he2.ome.tif") and old and np.allclose(old[0]["affine"], A), ims
    # Undo: back on the first image, the second one listed
    importer.undo_placement("s1")
    ims = importer.images("s1")
    assert ims[0]["path"].endswith("frame.ome.tif") and any(x["path"].endswith("he2.ome.tif") for x in ims), ims
    # removed
    assert cl.post("/api/sample/s1/images/remove", json={"id": r2["added"]}).status_code == 200
    assert r2["added"] not in [x["id"] for x in importer.images("s1")]
    # an image added "in the same pixels as" the frame follows it when the frame is placed again (another affine)
    pyramid(os.path.join(home, "ifsame.ome.tif"), base[..., 1])
    r3 = cl.post("/api/sample/s1/images/add", json={"path": os.path.join(home, "ifsame.ome.tif"), "kind": "IF", "same_as": "main"}).json()
    A2 = A + np.array([[0, 0, 7], [0, 0, -5]])
    importer.place_cells("s1", A2, os.path.join(home, "frame.ome.tif"), "H&E", by="file", progress=lambda *a: None)
    x3 = next(x for x in importer.images("s1") if x["id"] == r3["added"])
    assert np.allclose(x3["affine"], A2), "an image in the frame's pixels did not follow the frame"
    # the same file added again with another registration: other tiles (rev), not the cached warp
    info1 = {x["id"]: x["rev"] for x in cl.get("/api/sample/s1/he_info").json()["images"]}
    B2 = B + np.array([[0, 0, 30], [0, 0, 0]])
    pd.DataFrame([dict(zip(["a11", "a12", "a21", "a22", "tx", "ty"], [B2[0, 0], B2[0, 1], B2[1, 0], B2[1, 1], B2[0, 2], B2[1, 2]]))]).to_csv(reg, index=False)
    cl.post("/api/sample/s1/images/add", json={"path": os.path.join(home, "other.ome.tif"), "kind": "IF", "registration": reg})
    info2 = {x["id"]: x["rev"] for x in cl.get("/api/sample/s1/he_info").json()["images"]}
    assert info2[r["added"]] != info1[r["added"]], "the new placement of the same file kept its rev"
    v = webapp.view("s1")
    T_new = v.layer(r["added"], info2[r["added"]])[1]
    assert not np.allclose(T_new, T), "the warp of the new registration is the old one"
    # an extra image that cannot be read: a blank tile not kept, and the frame image not blamed
    os.remove(os.path.join(home, "ifsame.ome.tif")); render.unload(os.path.join(home, "ifsame.ome.tif")); v.forget_layers()
    t = cl.get("/api/sample/s1/he/8/0_0.jpg", params={"img": r3["added"], "iv": "x"})
    assert t.status_code == 200 and "no-store" in t.headers.get("cache-control", ""), t.headers
    assert not cl.get("/api/sample/s1/he_info").json()["image_error"], "the frame image was blamed for another one"
    # the Image window's image (fresh) waits for its own registration or alignment, though the sample has one
    m = json.load(open(config.SAMPLES / "s1" / "meta.json")); m["registration_affine"] = A.tolist(); json.dump(m, open(config.SAMPLES / "s1" / "meta.json", "w"))
    pyramid(os.path.join(home, "he3.ome.tif"), base)
    rq = cl.post("/api/sample/s1/image", json={"path": os.path.join(home, "he3.ome.tif"), "kind": "H&E", "fresh": True}).json()
    m = json.load(open(config.SAMPLES / "s1" / "meta.json"))
    assert rq.get("pending") and m["image_pending"].endswith("he3.ome.tif") and not m["he_path"].endswith("he3.ome.tif"), (rq, m["he_path"])
    # the instrument's morphology images of a Xenium source
    mf = config.SAMPLES / "s1" / "meta.json"
    morph = os.path.join(home, "xen", "morphology_focus"); os.makedirs(morph)
    for n in ("ch0000_dapi", "ch0001_atp1a1_cd45_e-cadherin", "ch0002_18s"):
        pyramid(os.path.join(morph, n + ".ome.tif"), base[..., 0])
    m = json.load(open(mf)); m["source"] = {**(m.get("source") or {}), "kind": "xenium", "pixel_size_um": 0.5,
                                            "morphology": os.path.join(morph, "ch0000_dapi.ome.tif")}
    json.dump(m, open(mf, "w"))
    names = [(x["name"], x["kind"]) for x in importer.images("s1") if x["source"] == "instrument"]
    assert names == [("DAPI", "DAPI"), ("ATP1A1 / CD45 / E-CADHERIN", "IF"), ("18S", "IF")], names
    print("smoke_images ok", f"(warped r = {rr:.3f})")
finally:
    shutil.rmtree(home, ignore_errors=True)
