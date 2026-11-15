"""Steps: Priors apart from Regions. Priors is a step to do (not optional): the organ's tissue priors checked in the
Prior window and accepted (POST /api/sample/<s>/priors/accept); done while they stay as accepted, to do again when
they change; Regions stays optional (outlines drawn). Runs in a temporary data folder."""
import os, shutil, tempfile
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_steps_")
os.environ["STAGENT_HOME"] = home; os.environ["STAGENT_SETTINGS_DIR"] = home + "/cfg"
from stagent import config, importer, hetissue
from stagent.session import Session
from starlette.testclient import TestClient
from stagent import webapp

quiet = lambda *a: None
try:
    h5 = os.path.join(home, "one.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    importer.import_sample("s1", h5, progress=quiet)
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    step = lambda r, k: next(x for x in r["steps"] if x["key"] == k)
    r = c.get("/api/sample/s1/steps").json()
    keys = [x["key"] for x in r["steps"]]
    assert keys.index("annotate") < keys.index("priors") < keys.index("regions"), keys
    assert step(r, "priors")["state"] == "todo", step(r, "priors")
    # the organ chosen in the Prior window (kept for the sample)
    assert c.post("/api/sample/s1/priors", json={"organ": "colon", "removed": [], "added": []}).status_code == 200
    p = step(c.get("/api/sample/s1/steps").json(), "priors")
    organ = p["organ"]
    assert organ and p["state"] == "todo" and "check them" in p["detail"], p
    # Accept: done; refused without an organ or for no sample
    assert c.post("/api/sample/s1/priors/accept", json={}).status_code == 400
    assert c.post("/api/sample/s1/priors/accept", json=None).status_code == 400
    assert c.post("/api/sample/nope/priors/accept", json={"organ": organ}).status_code == 404
    assert c.post("/api/sample/s1/priors/accept", json={"organ": organ}).status_code == 200
    p = step(c.get("/api/sample/s1/steps").json(), "priors")
    assert p["state"] == "done" and "accepted" in p["detail"], p
    # the priors changed (one removed): to do again
    first = hetissue.organ_priors(organ)["priors"][0]["id"]
    assert c.post("/api/sample/s1/priors", json={"organ": organ, "removed": [first], "added": []}).status_code == 200
    p = step(c.get("/api/sample/s1/steps").json(), "priors")
    assert p["state"] == "todo" and "changed since" in p["detail"], p
    c.post("/api/sample/s1/priors/accept", json={"organ": organ})
    assert step(c.get("/api/sample/s1/steps").json(), "priors")["state"] == "done"
    # with an annotation open: Regions optional (no outline), Priors as before
    s = Session.create("s1", "none", method="none")
    r = c.get(f"/api/s/{s.id}/steps").json()
    assert step(r, "regions")["state"] == "optional" and step(r, "priors")["state"] == "done", r["steps"]
    hetissue.set_priors(organ, [], [])  # (the user's changes of this organ undone: the shared settings of the test)
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
