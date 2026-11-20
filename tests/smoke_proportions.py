"""Cell type shares between conditions with the replicates as the samples (analysis.proportions, as propeller): the
null is calibrated and real changes are found on simulated slices; on a project of four slices (the example MERFISH
slice cut in quarters, two per condition) the composition by condition carries the test, Compare in the app runs it
and keeps its own CSV, one slice per condition is refused with a reason, and the agent's composition tool returns it.
Runs in a temporary data folder."""
import os, shutil, tempfile, time, types
import h5py
import numpy as np, pandas as pd
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_prop_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")  # (never the user's settings)
from stagent import analysis as A, importer, projects
from stagent.session import Session
rng = np.random.default_rng(1)


def fake(n_per=(3, 3), effect=None, K=15, cells=3000, conc=200):
    base = rng.dirichlet(np.full(K, 2.0))
    rows = []
    for gi, n in enumerate(n_per):
        for r in range(n):
            p = base.copy()
            if effect and gi == 1:
                p[0] *= effect
                p /= p.sum()
            lab = rng.choice(K, cells, p=rng.dirichlet(p * conc))
            rows.append(pd.DataFrame({"slice": f"s{gi}{r}", "condition": f"c{gi}", "cell type": [f"T{x}" for x in lab]}))
    L = pd.concat(rows, ignore_index=True)
    return A.Cells([types.SimpleNamespace(panel=["g"])], L, L["slice"].to_numpy(), pid="p")


try:
    ps = np.array([v["p"] for _ in range(150) for v in A.proportions(fake(), "condition", "cell type")["values"]])
    print(f"null: p < 0.05 in {(ps < 0.05).mean():.3f}")
    assert 0.02 < (ps < 0.05).mean() < 0.08, (ps < 0.05).mean()
    hit = np.mean([next(v for v in A.proportions(fake(effect=3.0), "condition", "cell type")["values"]
                        if v["value"] == "T0")["fdr"] < 0.05 for _ in range(60)])
    print(f"3x on one type, 3 vs 3 slices: found in {hit:.2f}")
    assert hit > 0.6, hit
    r3 = A.proportions(fake(n_per=(2, 3, 3)), "condition", "cell type")
    assert r3["test"] == "F" and len(r3["groups"]) == 3 and "F" in r3["values"][0]

    # a project of four quarters, conditions A A B B
    h5 = os.path.join(home, "q.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    with h5py.File(h5, "r+") as h:
        x, y = h["obs"]["x"][:], h["obs"]["y"][:]
        g = h["obs"].create_group("quarter")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = ((x >= np.median(x)).astype(np.int8) * 2 + (y >= np.median(y)).astype(np.int8))
        g.create_dataset("categories", data=np.array(["q1", "q2", "q3", "q4"], dtype=object), dtype=h5py.string_dtype())
    importer.import_project("demo", "h5ad", [h5], sections=["q1", "q2", "q3", "q4"], section_column="quarter",
                            progress=lambda *a: None)
    p = projects.list_projects()[0]
    projects.save(p["project_id"], p["name"], [(s, p["slices"][s]["section"], {"condition": "AB"[i // 2]})
                                               for i, s in enumerate(p["samples"])])
    MK = open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read()
    sids = [Session.create(s, "none", method="markers", markers_text=MK).id for s in p["samples"]]
    c = A.cells(pid=p["project_id"])
    x = A.crosstab(c, "condition", "lineage")
    t = x["test"]
    assert t["groups"] == ["A", "B"] and t["n_replicates"] == {"A": 2, "B": 2} and t["test"] == "t", t
    assert all({"mean_share", "t", "p", "fdr", "log2_ratio"} <= set(v) for v in t["values"])
    print("lineage shares A vs B:", [(v["value"], v["mean_share"], v["fdr"]) for v in t["values"][:3]])
    assert "test" not in A.crosstab(c, "slice", "lineage")  # (rows = the replicates themselves: no test)
    # one slice per condition: refused with the reason
    projects.save(p["project_id"], p["name"], [(s, p["slices"][s]["section"], {"condition": "ABBB"[i]})
                                               for i, s in enumerate(p["samples"])])
    x1 = A.crosstab(A.cells(pid=p["project_id"]), "condition", "lineage")
    assert "test" not in x1 and "2 or more" in x1["test_note"], x1.get("test_note")
    projects.save(p["project_id"], p["name"], [(s, p["slices"][s]["section"], {"condition": "AB"[i // 2]})
                                               for i, s in enumerate(p["samples"])])
    # the app: Compare -> Composition, all slices, rows = condition
    from starlette.testclient import TestClient
    from stagent import webapp
    cl = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    j = cl.post("/api/analysis/run", json={"session": sids[0], "kind": "composition", "scope": "project",
                                           "rows": "condition", "cols": "cell type", "replicate": "slice"}).json()
    for _ in range(300):
        st = cl.get(f"/api/job/{j['job']}").json()
        if st.get("done"):
            break
        time.sleep(0.2)
    assert not st.get("error"), st
    r = st["result"]
    assert r["test"]["replicate"] == "slice" and r["test"]["id"], r.get("test_note")
    csv = cl.get(f"/api/analysis/{r['test']['id']}.csv").text.splitlines()
    assert csv[0].startswith("cell type,n_cells,mean_share_A,mean_share_B,log2_ratio,t,p,fdr"), csv[0]
    # the agent's tool
    from stagent import mcp_server
    mcp_server._active["session"] = Session(sids[0])
    out = mcp_server.composition("condition", "lineage", scope="project")
    assert '"test"' in out and '"fdr"' in out
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
