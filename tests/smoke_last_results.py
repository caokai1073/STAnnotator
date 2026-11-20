"""Last result of analyses the agent runs and of analyses of all slices: the agent's trajectory, DE genes, top genes
and composition are kept as the open annotation's Last result (with their CSV), the app's suggestions poll carries a
results signature that changes with them (the app reloads Last result), and an analysis of all slices run from the
app is kept on every slice of the project. Runs in a temporary data folder."""
import os, shutil, tempfile, time
import h5py
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_last_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")  # (never the user's settings)
from stagent import importer, projects
from stagent.session import Session
try:
    h5 = os.path.join(home, "two.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    with h5py.File(h5, "r+") as h:
        x = h["obs"]["x"][:]
        g = h["obs"].create_group("half")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = (x >= np.median(x)).astype(np.int8)
        g.create_dataset("categories", data=np.array(["left", "right"], dtype=object), dtype=h5py.string_dtype())
    importer.import_project("demo", "h5ad", [h5], sections=["left", "right"], section_column="half", progress=lambda *a: None)
    p = projects.list_projects()[0]
    MK = open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read()
    sids = [Session.create(s, "none", method="markers", markers_text=MK).id for s in p["samples"]]
    from starlette.testclient import TestClient
    from stagent import webapp, mcp_server
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    sig0 = c.get(f"/api/s/{sids[0]}/plausibility/suggestions").json()["results"]
    # the agent's analyses on the open annotation
    mcp_server._active["session"] = s0 = Session(sids[0])
    s0.actor = "agent"
    mcp_server.trajectory("spatial", start={"by": "cell type", "values": ["Transit-amplifying cells"]})
    mcp_server.de_genes("cell type", "Transit-amplifying cells")
    mcp_server.group_markers("lineage")
    mcp_server.composition("lineage", "cell type")
    res = c.get(f"/api/s/{sids[0]}/results").json()
    assert {"trajectory_spatial", "de", "markers", "composition"} <= set(res), list(res)
    assert res["trajectory_spatial"]["request"]["kind"] == "trajectory" and res["de"]["request"]["kind"] == "de"
    for k in ("trajectory_spatial", "de", "markers", "composition"):
        r = c.get(f"/api/analysis/{res[k]['result']['id']}.csv")
        assert r.status_code == 200 and len(r.text.splitlines()) > 1, (k, r.status_code)
    sig1 = c.get(f"/api/s/{sids[0]}/plausibility/suggestions").json()["results"]
    assert sig1 != sig0 and len(sig1) == 4, sig1
    # an analysis of all slices, from the app: kept on both slices
    j = c.post("/api/analysis/run", json={"session": sids[0], "kind": "composition", "scope": "project",
                                         "rows": "slice", "cols": "lineage"}).json()
    for _ in range(300):
        st = c.get(f"/api/job/{j['job']}").json()
        if st.get("done"):
            break
        time.sleep(0.2)
    assert not st.get("error"), st
    other = c.get(f"/api/s/{sids[1]}/results").json()
    assert other.get("composition", {}).get("request", {}).get("scope") == "project", list(other)
    assert other["composition"]["result"]["id"] == st["result"]["id"]
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
