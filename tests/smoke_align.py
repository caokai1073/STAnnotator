"""Import the cells alone, add the H&E, align it without a registration file, compare with the registration file;
then Align automatically (the app's button: align.auto through /api/sample/<s>/align/auto, no agent) on the image
waiting again, the same comparison, and its Undo. Runs in a temporary data folder. Needs the SCP2927 example files (paths below)."""
import os, shutil, sys, tempfile, time
import numpy as np
H5AD = os.environ.get("SMOKE_H5AD", "/data2/celltypeagent/share/validation/scp2927/scp2927_2sections.h5ad")
HE = os.environ.get("SMOKE_HE", "/data/kai/endothelial/data/SCP2927/HE/human_slideid-0050423-108731-147328.ome.tif")
REG = os.path.join(os.path.dirname(__file__), "..", "examples", "scp2927", "registration.csv")
home = tempfile.mkdtemp(prefix="stagent_align_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(os.environ["STAGENT_HOME"], "cfg")  # (never the user's settings)
from stagent import importer, align
try:
    t = time.time()
    importer.import_sample("t", [H5AD], section_column="donor_id", section="147328", progress=lambda *a: None)
    print("cells only:", importer._meta("t")["image_kind"], f"{time.time() - t:.0f}s")
    print("image:", importer.set_image("t", HE))
    s = align.Sample("t")
    st = align.search(s)
    print("candidates:", [(c["id"], c["orientation_name"], c["um_per_px"], c["score"]) for c in st["candidates"][:3]])
    print("refine:", align.refine(s, "K1")["rounds"][-1])
    print("apply:", align.apply(s, "K1"))
    A = importer.placement("t")
    R = importer.load_registration(REG, "147328")
    xy = align.Sample("t").cells[["x", "y"]].to_numpy()
    d = np.hypot(*((xy @ A[:, :2].T + A[:, 2]) - (xy @ R[:, :2].T + R[:, 2])).T) * 0.69
    print(f"vs the registration file: median {np.median(d):.1f} um (the file itself is ~13 um off the nuclei)")
    print("undo:", importer.undo_placement("t"), "| registration file:", importer.set_registration("t", REG, "147328"))
    # a wrong pixel size in the image file (0.2125 for 0.69): the wide search, then around its best size in fine steps
    real_size = align.image_um_per_px
    align.image_um_per_px = lambda p: 0.2125
    try:
        st = align.search(align.Sample("t"))
    finally:
        align.image_um_per_px = real_size
    best = st["candidates"][0]
    print("wrong size in the file:", best["orientation_name"], best["um_per_px"], best["score"], st.get("note", "")[:60])
    assert st.get("searched_beyond_file_um_per_px") and abs(np.log(best["um_per_px"] / 0.691)) < 0.025, best

    # Align automatically: a sample whose H&E waits, aligned by the app alone
    importer.import_sample("u", [H5AD], section_column="donor_id", section="147328", progress=lambda *a: None)
    importer.set_image("u", HE)
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    t = time.time()
    r = c.post("/api/sample/u/align/auto", json={})
    assert r.status_code == 200, r.text
    job = r.json()["job"]
    assert c.post("/api/sample/u/align/auto", json={}).status_code == 409, "a second run while one runs"
    while not (j := c.get(f"/api/job/{job}").json())["done"]:
        time.sleep(2)
    assert not j.get("error"), j
    res = j["result"]
    print("auto:", res, f"{time.time() - t:.0f}s")
    assert res["applied"] and res["verdict"] in ("fine", "local shifts"), res
    m = importer._meta("u")
    assert not m.get("image_pending") and m.get("registration_by") == "auto", m.get("registration_by")
    A = importer.placement("u")
    xy = align.Sample("u").cells[["x", "y"]].to_numpy()
    d = np.hypot(*((xy @ A[:, :2].T + A[:, 2]) - (xy @ R[:, :2].T + R[:, 2])).T) * 0.69
    print(f"auto vs the registration file: median {np.median(d):.1f} um")
    assert np.median(d) < 40, "the automatic alignment is not where the registration file puts the cells"
    # a Check that advises aligning again: the Image step is to do (Align again), not done
    from stagent import regcheck, config
    from stagent.osutil import dump_json
    mf = config.SAMPLES / "u" / "meta.json"
    mm = importer._meta("u")
    mm["registration_checked"] = {"key": regcheck._key("u", importer.image_of("u")), "time": "2099-01-01T00:00:00",
                                  "verdict": "fine", "median_um": 22.2, "new_alignment": True}
    dump_json(mm, mf, indent=1)
    step = next(x for x in c.get("/api/sample/u/steps").json()["steps"] if x["key"] == "registration")
    assert step["state"] == "todo" and step.get("advise") == "align", step
    mm["registration_checked"]["new_alignment"] = False
    dump_json(mm, mf, indent=1)
    step = next(x for x in c.get("/api/sample/u/steps").json()["steps"] if x["key"] == "registration")
    assert step["state"] == "done", step
    print("Image step after a check advising to align again:", "todo, Align again")
    assert c.post("/api/sample/u/placement/undo", json={}).status_code == 200
    assert importer._meta("u").get("image_pending"), "Undo did not put the image back to waiting"
finally:
    shutil.rmtree(home)
