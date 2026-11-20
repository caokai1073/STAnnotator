"""Protein / IF stains as evidence (protein.py): a CosMx slide made from the example MERFISH slice (PanCK from its
epithelial genes, CD45 from Ptprc, DAPI flat) is imported with its stains kept (protein.parquet); an annotation with
stains (the default) reweighs the types (cells positive for a stain move to the types mostly positive for it; DAPI is
left out; rates taken at two levels only at a clear gap, membrane stains left out), one without them and annotations made before (no "stains" in the scheme) are as before, the agent's review
names the stains, the samples list shows them, and the notebook replays the annotation with them. Runs in a
temporary data folder."""
import os, shutil, tempfile, json
import numpy as np, pandas as pd
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_stains_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")  # (never the user's settings)
from stagent import importer, protein, notebook
from stagent.session import Session
rng = np.random.default_rng(0)
try:
    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    r = importer.read_h5ad(h5)
    X, genes = r["X"].tocsr(), list(r["genes"])
    n = X.shape[0]
    gi = {g: i for i, g in enumerate(genes)}
    epi = np.asarray(X[:, [gi[g] for g in ("Vil1", "Cdh17", "Epcam", "Krt8", "Krt20") if g in gi]].sum(1)).ravel()
    imm = np.asarray(X[:, [gi[g] for g in ("Ptprc",) if g in gi]].sum(1)).ravel()
    d = os.path.join(home, "cosmx")
    os.makedirs(d)
    fov, cid = np.ones(n, int), np.arange(1, n + 1)
    cm = pd.DataFrame(X.toarray().astype(int), columns=genes)
    cm.insert(0, "cell_ID", cid)
    cm.insert(0, "fov", fov)
    cm.to_csv(os.path.join(d, "demo_exprMat_file.csv"), index=False)
    md = pd.DataFrame({"fov": fov, "cell_ID": cid, "CenterX_global_px": r["xy"][:, 0] / importer.COSMX_PX_UM,
                       "CenterY_global_px": r["xy"][:, 1] / importer.COSMX_PX_UM,
                       "Mean.PanCK": np.where(epi > 0, 3000, 200) * rng.uniform(0.7, 1.3, n),
                       "Mean.CD45": np.where(imm > 0, 2500, 150) * rng.uniform(0.7, 1.3, n),
                       "Mean.DAPI": rng.uniform(900, 1100, n)})
    md.to_csv(os.path.join(d, "demo_metadata_file.csv"), index=False)
    importer.import_flat("cosmx", "cx", d, progress=lambda *a: None)
    st = protein.ensure("cx")
    assert st is not None and list(st.columns) == ["cell_id", "PanCK", "CD45", "DAPI"], None if st is None else list(st.columns)
    assert protein.available("cx") == ["PanCK", "CD45"]
    MK = open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read()
    on = Session.create("cx", "none", method="markers", markers_text=MK)
    off = Session.create("cx", "none", method="markers", markers_text=MK, stains=False)
    assert on.scheme["stains"] and not off.scheme["stains"] and on.vocabulary()["stains"] == ["PanCK", "CD45"]
    a, b = np.asarray(on.labels()), np.asarray(off.labels())
    lin = lambda s: np.array([s.lineage_of_type.get(t, t) for t in s.labels()])
    pck = md["Mean.PanCK"].to_numpy() > 1000
    ids = on.cells.cell_id.astype(str).to_numpy()
    pos = pd.Series(pck, index=[f"{f}_{c}" for f, c in zip(fov, cid)]).reindex(ids).to_numpy() if not set(ids) <= set(map(str, cid)) else pck
    e_on, e_off = (lin(on) == "Epithelial"), (lin(off) == "Epithelial")
    print(f"changed {(a != b).sum():,} of {len(a):,}; PanCK+ called epithelial: {e_off[pos].mean():.3f} -> {e_on[pos].mean():.3f}; "
          f"PanCK- called epithelial: {e_off[~pos].mean():.3f} -> {e_on[~pos].mean():.3f}")
    assert (a != b).sum() > 0 and e_on[pos].mean() >= e_off[pos].mean() and e_on[~pos].mean() <= e_off[~pos].mean()
    # an annotation made before (no "stains" in its scheme) stays as it was
    old = Session.create("cx", "none", method="markers", markers_text=MK)
    sj = json.load(open(old.dir / "session.json"))
    sj["scheme"].pop("stains")
    json.dump(sj, open(old.dir / "session.json", "w"))
    assert (np.asarray(Session(old.id).labels()) == b).all()
    # the agent's review, the samples list
    from stagent import mcp_server
    det = mcp_server._d_start_session({"sample_id": "cx", "markers_text": MK})["details"]
    assert any("PanCK, CD45" in x for x in det), det
    from starlette.testclient import TestClient
    from stagent import webapp
    cl = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    smp = next(x for x in cl.get("/api/samples").json() if x["sample_id"] == "cx")
    assert smp["stains"] == ["PanCK", "CD45"], smp["stains"]
    # the notebook replays it with the stains
    p = notebook.build(on)
    nb = json.load(open(p))
    code = ["".join(c["source"]) for c in nb["cells"] if c["cell_type"] == "code" and "%pip" not in "".join(c["source"])]
    assert any("sa.protein.apply" in c for c in code)
    g = {"__name__": "__nb__"}
    os.chdir(os.path.dirname(p))
    import io, contextlib
    out = io.StringIO()
    with contextlib.redirect_stdout(out):
        for c in code:
            exec(compile(c, "cell", "exec"), g)
            if "ORIGINAL" in g:
                break
    assert "labels identical to the original" in out.getvalue(), out.getvalue()[-300:]
    # two levels only at a clear gap (the Mac review's ASHUC1 rates): Membrane-like rates running on from type to
    # type (ILCs 0.67, CD8 T 0.79 ...) give every type the same rate; PanCK splits at its gap (Mast cells 0.55 low)
    # (ILCs, CD8 T, B, naive CD4 T, plasma, fibroblasts, macrophages, blood vessel EC, mast | goblet, cycling TA, unassigned)
    n = np.array([18954, 22147, 32822, 353, 258126, 106631, 26936, 28733, 55354, 50999, 36754, 63411])
    memb = np.array([0.67, 0.79, 0.77, 0.71, 0.65, 0.36, 0.51, 0.47, 0.48, 0.98, 0.98, 0.60])
    pck = np.array([0.24, 0.10, 0.07, 0.09, 0.25, 0.31, 0.36, 0.40, 0.55, 0.98, 0.98, 0.81])
    R = protein._two_levels(np.c_[memb, pck], n)
    assert np.ptp(R[:, 0]) == 0, R[:, 0]
    lo, hi = R[:9, 1], R[9:, 1]
    assert np.ptp(lo) == 0 and np.ptp(hi) == 0 and hi[0] > 0.8 > 0.4 > lo[0], R[:, 1]
    assert protein.skipped("Membrane") and protein.skipped("MembraneStain") and protein.skipped("DAPI") and not protein.skipped("PanCK")
    print("ok")
finally:
    os.chdir("/")
    shutil.rmtree(home, ignore_errors=True)
