"""Spatially variable genes (svg.py) and gene program scores (programs.py): Moran's I equals the formula, is about 0
for noise and drops when a pattern only follows the cell types (beyond types); on the example MERFISH slice the app's
Compare runs both kinds, the program colors the map (program:<key>, legend, cache key), its CSV is kept, a program
with too few panel genes is refused, and the agent's tools return them. Runs in a temporary data folder."""
import os, shutil, tempfile, time
import numpy as np, scipy.sparse as sp
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_svg_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")  # (never the user's settings)
from stagent import importer, programs, svg
from stagent.session import Session
rng = np.random.default_rng(0)
try:
    n = 400
    xy = rng.uniform(0, 100, (n, 2))
    X = np.stack([rng.normal(size=n), np.sin(xy[:, 0] / 10) + 0.1 * rng.normal(size=n)], 1)
    W = svg.weights(xy, 6)
    I = svg.morans_i(sp.csr_matrix(X), W)
    Wd = W.toarray()
    ref = [n / Wd.sum() * ((x - x.mean()) @ Wd @ (x - x.mean())) / ((x - x.mean()) @ (x - x.mean())) for x in X.T]
    assert np.allclose(I, ref, atol=1e-4), (I, ref)
    assert abs(I[0]) < 0.1 and I[1] > 0.8, I
    t = (xy[:, 0] > 50).astype(int)
    y = (t * 2.0 + rng.normal(size=n))[:, None]
    a, b = svg.morans_i(sp.csr_matrix(y), W)[0], svg.morans_i(sp.csr_matrix(y), W, t)[0]
    print(f"a pattern of the types: I {a:.2f}, beyond the types {b:.2f}")
    assert a > 0.3 and abs(b) < 0.1

    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    importer.import_sample("s1", h5, progress=lambda *a: None)
    s = Session.create("s1", "none", method="markers", markers_text=open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read())
    s.labels()
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")

    def run(body):
        j = c.post("/api/analysis/run", json={"session": s.id, "scope": "slice", **body}).json()
        assert "job" in j, j
        for _ in range(600):
            st = c.get(f"/api/job/{j['job']}").json()
            if st.get("done"):
                return st
            time.sleep(0.2)
        raise AssertionError("timeout")
    st = run({"kind": "spatial", "beyond_types": True})
    r = st["result"]
    assert not st.get("error") and r["beyond_types"] and len(r["genes"]) == 60 and r["genes"][0]["morans_i"] >= r["genes"][-1]["morans_i"]
    print("spatial genes beyond types:", [(g["gene"], g["morans_i"]) for g in r["genes"][:5]])
    assert c.get(f"/api/analysis/{r['id']}.csv").text.startswith("gene,morans_i,z,fdr,pct_detected")
    st = run({"kind": "program", "program": "Proliferation", "by": "cell type"})
    r = st["result"]
    assert not st.get("error") and r["mode"] == "program:proliferation" and len(r["genes"]) >= 3, st
    print("proliferation by type:", [(g["group"], g["mean"]) for g in r["groups"][:3]], "genes", r["genes"])
    pg = c.get(f"/api/s/{s.id}/programs").json()
    assert "proliferation" in pg["scored"] and "Hypoxia" in pg["builtin"]
    v = s.versions()[-1]
    col = c.get(f"/api/s/{s.id}/colors.bin", params={"version": v, "mode": "program:proliferation", "only": ""})
    assert col.status_code == 200 and len(col.content) == 4 * len(s.cells), col.status_code
    leg = c.get(f"/api/s/{s.id}/legend", params={"version": v, "mode": "program:proliferation"})
    assert leg.status_code == 200 and sum(x["count"] for x in leg.json()["items"] if "count" in x) == len(s.cells), leg.text[:200]
    k1 = webapp.mode_key(s, "program:proliferation")
    st = run({"kind": "program", "genes": "Mki67 Cenpf Pcna Mcm2", "program": "Proliferation", "by": "cell type"})
    assert not st.get("error"), st
    assert webapp.mode_key(s, "program:proliferation") != k1, (k1, webapp.mode_key(s, "program:proliferation"))  # scored again: new colors
    st = run({"kind": "program", "genes": "Mki67 NOPE1 NOPE2", "program": "Tiny", "by": "cell type"})
    assert "only 1 of its 3 genes" in st["error"], st
    from stagent import mcp_server
    mcp_server._active["session"] = Session(s.id)
    assert '"morans_i"' in mcp_server.spatially_variable_genes(top=5)
    assert '"pct_above"' in mcp_server.program_scores("Hypoxia" if len(programs.resolve("Hypoxia", None, s.panel)[2]) >= 3 else "Proliferation")
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
