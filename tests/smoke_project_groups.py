"""Groups of a project's slices (Data → Group): a grouping (a metadata column, one value per slice) saved through the
endpoint is kept in the project and each slice's meta.json, Compare lists it (Compare by) when its values differ
between slices, and DE genes compare its groups; names that are empty, repeated or the app's own are refused.
Runs in a temporary data folder."""
import json, os, shutil, tempfile, time
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_groups_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")
import h5py
from stagent import importer, joint, projects, config
try:
    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    with h5py.File(h5, "r+") as h:  # three slices: thirds of the tissue
        x = h["obs"]["x"][:]
        g = h["obs"].create_group("third")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = np.digitize(x, np.quantile(x, [1 / 3, 2 / 3])).astype(np.int8)
        g.create_dataset("categories", data=np.array(["a", "b", "c"], dtype=object), dtype=h5py.string_dtype())
    importer.import_project("demo", "h5ad", [h5], sections=["a", "b", "c"], section_column="third", progress=lambda *a: None)
    joint.annotate("demo", open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read())
    sids = list(joint.latest_sessions("demo").values())
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    pid = projects.project_of(json.load(open(config.SESSIONS / sids[0] / "session.json"))["sample_id"])
    g0 = c.get(f"/api/project/{pid}/groups").json()
    ids = [x["sample_id"] for x in g0["slices"]]
    assert len(ids) == 3, g0
    cond = {ids[0]: "UC", ids[1]: "UC", ids[2]: "control"}
    r = c.post(f"/api/project/{pid}/groups", json={"columns": {"condition": cond, "batch": {ids[0]: "1", ids[1]: "1", ids[2]: "1"}}})
    assert r.status_code == 200 and r.json()["columns"]["condition"] == cond, r.text
    assert json.load(open(config.SAMPLES / ids[2] / "meta.json"))["metadata"] == {"condition": "control", "batch": "1"}
    opts = c.post("/api/analysis/options", json={"session": sids[0], "scope": "project"}).json()
    assert opts["groupings"].get("condition") and set(opts["groupings"]["condition"]) == {"UC", "control"}, opts["groupings"].keys()
    assert "batch" not in opts["groupings"]  # (one value on every slice: nothing to compare)
    j = c.post("/api/analysis/run", json={"session": sids[0], "scope": "project", "kind": "de", "by": "condition", "group": "UC", "top": 5}).json()
    for _ in range(600):
        st = c.get(f"/api/job/{j['job']}").json()
        if st.get("done"):
            break
        time.sleep(0.2)
    assert not st.get("error") and st["result"]["n_group"] > 0 and st["result"]["up"], st.get("error")
    print("DE UC vs control:", [x["gene"] for x in st["result"]["up"]][:5])
    for bad in ({"": cond}, {"cell type": cond}, {"Slice": cond}, {"x" * 41: cond}):
        assert c.post(f"/api/project/{pid}/groups", json={"columns": bad}).status_code == 400, bad
    r = c.post(f"/api/project/{pid}/groups", json={"columns": {}}).json()  # all removed
    assert r["columns"] == {} and json.load(open(config.SAMPLES / ids[0] / "meta.json"))["metadata"] == {}
    assert "condition" not in c.post("/api/analysis/options", json={"session": sids[0], "scope": "project"}).json()["groupings"]
    print("smoke_project_groups ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
