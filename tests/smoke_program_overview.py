"""Program colors on the same scale across the slices of a project, and the colors the project's overview asks for:
a program scored on all slices together keeps the cutoffs of every slice's cells in each slice (the map, its legend
and the overview bin by them, so a color means the same score on every slice); a program scored on one slice keeps
its own; one scored before the cutoffs were kept bins by its own slice as before; genes color by normalized expression
on a continuous scale, one for all slices of a project (or one given); a gene missing in a panel is refused. Runs in a temporary data folder."""
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

    # all slices together: each slice keeps the cutoffs of every cell
    run({"kind": "program", "program": "Proliferation", "by": "cell type", "scope": "project"}, sids[0])
    infos = [programs.saved(s)["proliferation"] for s in ss]
    assert infos[0]["cutoffs"] == infos[1]["cutoffs"] and infos[0]["scope"] == "project" and infos[0]["n_slices"] == 2, infos
    v = [programs.values(s, "proliferation") for s in ss]
    allv = np.concatenate(v)
    assert np.allclose(infos[0]["cutoffs"], np.quantile(allv[allv > 0], [0.5, 0.8, 0.95]), atol=1e-5)
    b = [programs.cell_bins(s, "proliferation") for s in ss]
    both = np.concatenate(b)
    top = (both == 4).sum() / (allv > 0).sum()  # the top 5% bin: 5% of the cells above background over both slices
    print("top bin share over both slices", round(top, 4), "per slice", [round((x == 4).sum() / max(1, (y > 0).sum()), 4) for x, y in zip(b, v)])
    assert abs(top - 0.05) < 0.01, top
    for s, bb in zip(ss, b):  # the map, its legend and the overview's colors bin the same way
        ver = s.versions()[-1]
        leg = c.get(f"/api/s/{s.id}/legend", params={"version": ver, "mode": "program:proliferation"}).json()["items"]
        n = {x["key"]: x["count"] for x in leg}
        assert [n[k] for k, *_ in programs.BINS] == np.bincount(bb, minlength=5).tolist(), (n, np.bincount(bb))
        col = np.frombuffer(c.get(f"/api/s/{s.id}/colors.bin", params={"version": ver, "mode": "program:proliferation"}).content, np.uint8).reshape(-1, 4)
        assert (col[:, :3] == np.array([x[2] for x in programs.BINS], np.uint8)[bb]).all()
    # one slice alone: its own cutoffs (the overview tells this slice apart: another run)
    run({"kind": "program", "program": "Proliferation", "by": "cell type", "scope": "slice"}, sids[1])
    i1 = programs.saved(ss[1])["proliferation"]
    assert i1["scope"] == "slice" and i1["run_id"] != infos[0]["run_id"] and np.allclose(i1["cutoffs"], programs.cutoffs(programs.values(ss[1], "proliferation")))
    assert programs.saved(ss[0])["proliferation"]["run_id"] == infos[0]["run_id"]
    # scored before the cutoffs were kept: its own slice's quantiles, as before
    from stagent.osutil import dump_json
    allp = programs.saved(ss[1]); del allp["proliferation"]["cutoffs"]; dump_json(allp, ss[1].dir / "programs" / "programs.json")
    vv = programs.values(ss[1], "proliferation")
    assert (programs.cell_bins(ss[1], "proliferation") == programs.bins(vv)).all()
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
    assert projects.list_projects()[0]["samples"]
    print("smoke_program_overview ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
