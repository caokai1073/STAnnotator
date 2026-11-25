"""Program and gene colors on one continuous scale across the slices of a project (the map, its legend and the
project's overview): its top the highest of the slices' own (the 99th percentile of the cells above background, or
where the gene is detected), a slice scored again on its own joins it, one scored before the tops were kept is
read from its scores, a top can be given; genes color by normalized expression; a gene missing in a panel is refused. Runs in a temporary data folder."""
import os, shutil, tempfile, time
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_progov_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")  # (never the user's settings)
import h5py
from stagent import importer, joint, programs, projects
from stagent.session import Session
try:
    h5 = os.path.join(home, "a.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    with h5py.File(h5, "r+") as h:  # two slices: the left and the right half
        x = h["obs"]["x"][:]
        g = h["obs"].create_group("half")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = (x >= (x.min() + x.max()) / 2).astype(np.int8)
        g.create_dataset("categories", data=np.array(["left", "right"], dtype=object), dtype=h5py.string_dtype())
    importer.import_project("demo", "h5ad", [h5], sections=["left", "right"], section_column="half", progress=lambda *a: None)
    joint.annotate("demo", open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read())
    sids = list(joint.latest_sessions("demo").values())
    ss = [Session(x) for x in sids]
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")

    def run(body, sid):
        j = c.post("/api/analysis/run", json={"session": sid, **body}).json()
        assert "job" in j, j
        for _ in range(900):
            st = c.get(f"/api/job/{j['job']}").json()
            if st.get("done"):
                assert not st.get("error"), st
                return st["result"]
            time.sleep(0.2)
        raise AssertionError("timeout")

    def colors(s, mode):
        r = c.get(f"/api/s/{s.id}/colors.bin", params={"version": s.versions()[-1], "mode": mode})
        assert r.status_code == 200, r.text
        return np.frombuffer(r.content, np.uint8).reshape(-1, 4)

    def scale(s, mode):
        return c.get(f"/api/s/{s.id}/legend", params={"version": s.versions()[-1], "mode": mode}).json()["scale"]

    def ramp(v, top):  # the expected colors: at or below 0 faint gray, above it the color bar up to top
        out = np.tile(np.array(webapp.GENE_NONE, np.uint8), (len(v), 1))
        on = v > 0
        t = np.clip(v[on] / top, 0, 1)
        out[on, :3] = webapp.GENE_LUT[(t * 255).round().astype(int)]
        return out, on

    # programs: a continuous scale, its top the highest of the project's slices' tops (99th percentile above 0)
    run({"kind": "program", "program": "Proliferation", "by": "cell type", "scope": "project"}, sids[0])
    infos = [programs.saved(s)["proliferation"] for s in ss]
    v = [programs.values(s, "proliferation") for s in ss]
    own = [float(np.quantile(x[x > 0], 0.99)) for x in v]
    assert all(abs(i["slice_top"] - o) < 1e-3 for i, o in zip(infos, own)) and infos[0]["scope"] == "project", infos
    for s, vv in zip(ss, v):
        sc = scale(s, "program:proliferation")
        assert abs(sc["top"] - max(own)) < 1e-3 and sc["n_slices"] == 2 and sc["kind"] == "program", sc
        want, on = ramp(vv, sc["top"])
        col = colors(s, "program:proliferation")
        assert np.abs(col[:, :3].astype(int) - want[:, :3]).max() <= 1 and (col[~on] == webapp.GENE_NONE).all()
        leg = {x["key"]: x["count"] for x in c.get(f"/api/s/{s.id}/legend", params={"version": s.versions()[-1], "mode": "program:proliferation"}).json()["items"]}
        assert leg == {"above background": int(on.sum()), "at or below background": int((~on).sum())}, leg
    print("program tops per slice", [round(x, 4) for x in own], "project top", scale(ss[0], "program:proliferation")["top"])
    # one slice scored again on its own: its new top joins the project's scale (both slices, recolored)
    before = colors(ss[0], "program:proliferation").copy()
    run({"kind": "program", "program": "Proliferation", "by": "cell type", "scope": "slice"}, sids[1])
    i1 = programs.saved(ss[1])["proliferation"]
    assert i1["scope"] == "slice" and i1["run_id"] != infos[0]["run_id"]
    t1 = programs.top(programs.values(ss[1], "proliferation"))
    top = max(own[0], t1)
    assert abs(scale(ss[0], "program:proliferation")["top"] - top) < 1e-3 and abs(scale(ss[1], "program:proliferation")["top"] - top) < 1e-3
    if abs(top - max(own)) > 1e-3:
        assert (colors(ss[0], "program:proliferation") != before).any()
    # scored before the tops were kept: computed from its scores
    from stagent.osutil import dump_json
    allp = programs.saved(ss[1]); del allp["proliferation"]["slice_top"]; dump_json(allp, ss[1].dir / "programs" / "programs.json")
    assert abs(programs.project_top(ss[0], "proliferation")[0] - top) < 1e-3
    # a given top (the overview of all samples)
    want, _ = ramp(v[0], 2 * top)
    assert np.abs(colors(ss[0], f"program:proliferation~{2 * top}")[:, :3].astype(int) - want[:, :3]).max() <= 1
    # genes: normalized expression (as model.lognorm) on a continuous scale, its top the highest of the project's
    # slices' 99th percentiles of the cells it is detected in (the map and the overview alike), or the one given;
    # not detected faint gray
    from stagent import model
    gene = Session(sids[0]).panel[0]
    refs = [np.asarray(model.lognorm(s.counts)[:, s.panel.index(gene)].todense()).ravel() for s in ss]
    own = [float(np.quantile(r[r > 0], 0.99)) for r in refs]
    tops = []
    for s, ref, mine in zip(ss, refs, own):
        ver = s.versions()[-1]
        leg = c.get(f"/api/s/{s.id}/legend", params={"version": ver, "mode": f"gene:{gene.lower()}"}).json()
        sc = leg["scale"]
        assert abs(sc["slice_top"] - mine) < 1e-3 and abs(sc["top"] - max(own)) < 1e-3 and sc["n_slices"] == 2, (sc, own)
        assert {x["key"]: x["count"] for x in leg["items"]} == {"detected": int((ref > 0).sum()), "not detected": int((ref == 0).sum())}
        tops.append(sc["top"])
        r = c.get(f"/api/s/{s.id}/colors.bin", params={"version": ver, "mode": f"gene:{gene.lower()}"})
        assert r.status_code == 200 and len(r.content) == 4 * len(s.cells), r.status_code
        col = np.frombuffer(r.content, np.uint8).reshape(-1, 4)
        assert (col[ref == 0] == webapp.GENE_NONE).all()
        on = ref > 0
        t = np.clip(ref[on] / sc["top"], 0, 1)
        assert np.abs(col[on, :3].astype(int) - webapp.GENE_LUT[(t * 255).round().astype(int)]).max() <= 1
        assert (col[on, 3] >= 215).all() and col[on][t >= 1, 3].min() == 255
        hi = ref > np.quantile(ref[on], 0.5)  # the same values, the same colors on another scale's top
        r2 = c.get(f"/api/s/{s.id}/colors.bin", params={"version": ver, "mode": f"gene:{gene}~{2 * sc['top']}"})
        col2 = np.frombuffer(r2.content, np.uint8).reshape(-1, 4)
        t2 = np.clip(ref[on] / (2 * sc["top"]), 0, 1)
        assert np.abs(col2[on, :3].astype(int) - webapp.GENE_LUT[(t2 * 255).round().astype(int)]).max() <= 1, hi.sum()
        assert c.get(f"/api/s/{s.id}/legend", params={"version": ver, "mode": f"gene:{gene}~{2 * sc['top']}"}).json()["scale"]["top"] == 2 * sc["top"]
    assert tops[0] == tops[1], tops
    print("gene scale: slices' own tops", [round(x, 3) for x in own], "project top", tops[0])
    r = c.get(f"/api/s/{sids[0]}/colors.bin", params={"version": ss[0].versions()[-1], "mode": "gene:NOPE_GENE"})
    assert r.status_code == 400 and "panel" in r.json()["error"], r.text
    # every gene's top of a slice in one pass, kept next to its counts: the same as gene by gene
    from stagent import config
    assert (config.SAMPLES / ss[0].sample_id / webapp.GENE_TOPS_FILE).exists()
    for g in ss[0].panel[:20]:
        assert webapp.slice_top(ss[0].sample_id, g) == webapp.gene_top(webapp.gene_norm_of(ss[0].sample_id, g)), g
    assert webapp.project_top(ss[0].sample_id, "NOT_A_GENE") == (1.0, 1)
    # a program's scale comes from each slice's latest annotation: one annotated again without it no longer counts
    import json as _json
    from stagent.session import Session as _S
    newer = _S.create(ss[1].sample_id, method="markers", markers_text=open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read(), knn="exact")
    assert programs.sample_top(ss[1].sample_id, "proliferation") is None
    assert abs(programs.project_top(ss[0], "proliferation")[0] - programs.top(programs.values(ss[0], "proliferation"))) < 1e-3
    assert programs.project_top(ss[0], "proliferation")[1] == 1
    assert projects.list_projects()[0]["samples"]
    # a renamed slice: downloads are named by its new name (with its project's), not the id it was imported under
    assert c.post(f"/api/sample/{ss[0].sample_id}/rename", json={"name": "Left half"}).status_code == 200
    v0 = ss[0].versions()[-1]
    for url in (f"/api/s/{sids[0]}/export.csv?version={v0}", f"/api/sample/{ss[0].sample_id}/registration.csv"):  # (no regions here)
        r = c.get(url)
        cd = r.headers.get("content-disposition", "")
        assert r.status_code == 200 and "demo_Left_half" in cd and ss[0].sample_id not in cd, (url, r.status_code, cd)
    assert projects.file_stem(ss[1].sample_id) == "demo_right"
    print("smoke_program_overview ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
