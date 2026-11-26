"""Export > Reproduction guide for agents: the zip of an annotation (with regions, PCA / UMAP / Leiden, niches and
kept analyses) holds REPRODUCE.md with every step, manifest.json with the checks, data/ (cell table, marker list,
regions, registration, the analyses' tables, colors) and, when asked, the notebook. Runs in a temporary data folder."""
import io, json, os, shutil, tempfile, time, zipfile
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_guide_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(home, "cfg")
from stagent import importer, embedding, niches, __version__
from stagent.session import Session
try:
    importer.import_sample("demo", os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    s = Session.create("demo", "none", method="markers", knn="exact",
                       markers_text=open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read())
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    x0, x1, y0, y1 = np.quantile(px, 0.2), np.quantile(px, 0.8), np.quantile(py, 0.2), np.quantile(py, 0.8)
    s.add_region("inner", [[x0, y0], [x1, y0], [x1, y1], [x0, y1]])
    s = Session(s.id)
    embedding.run(s)
    niches.run(s, n_niches=4, identity="cell types")
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    v = s.versions()[-1]
    types = list(s.types)
    runs = [("/api/analysis/run", {"kind": "de", "by": "cell type", "group": types[0]}),
            ("/api/analysis/run", {"kind": "composition", "rows": "region", "cols": "cell type"}),
            ("/api/communication/run", {"kind": "run", "by": "cell type"}),
            ("/api/coloc/run", {"by": "cell type"}),
            ("/api/trajectory/run", {"method": "spatial", "by": "cell type",
                                     "start": {"by": "cell type", "values": ["Transit-amplifying cells"]},
                                     "end": {"by": "cell type", "values": ["Colonocytes"]}})]
    for url, req in runs:
        j = c.post(url, json={"session": s.id, "scope": "slice", "version": v, **req}).json()
        assert "job" in j, (url, j)
        for _ in range(600):
            st = c.get(f"/api/job/{j['job']}").json()
            if st.get("done"):
                break
            time.sleep(0.2)
        assert not st.get("error"), st.get("error")
    for nb in (False, True):
        r = c.get(f"/api/s/{s.id}/export_guide.zip?version={v}&notebook={int(nb)}")
        assert r.status_code == 200, r.text
        z = zipfile.ZipFile(io.BytesIO(r.content))
        names = z.namelist()
        root = names[0].split("/")[0]
        assert root == f"demo_v{v}_reproduce", root
        have = {n.split("/", 1)[1] for n in names}
        need = {"REPRODUCE.md", "manifest.json", "data/cells.csv", "data/marker_list.txt", "data/marker_list.csv",
                "data/regions.geojson", "data/colors.json", "data/state.json", "data/results/de.json", "data/results/de.csv",
                "data/results/composition.json", "data/analyses_log.json", "data/CHANGELOG.md",
                "data/results/run.json", "data/results/run.csv", "data/results/coloc.json", "data/results/coloc.csv",
                "data/results/trajectory_spatial.json", "data/results/trajectory_spatial.csv", "data/trajectory_spatial.csv"}
        assert need <= have, need - have
        assert ("analysis.ipynb" in have) == nb, have
        M = json.loads(z.read(f"{root}/manifest.json"))
        assert not M["skipped"], M["skipped"]
        guide = z.read(f"{root}/REPRODUCE.md").decode()
        for h in ("## How to use this guide", "Software", "Input data", "Cells and genes kept", "Annotation",
                  "PCA, UMAP and Leiden", "Niches", "Analyses kept", "Figures", "Determinism", "Appendix"):
            assert h in guide, h
        from stagent import reproguide
        assert f"/blob/{reproguide._ref()}/stagent/markers.py#L" in guide  # links into this build's code
        assert "Leiden resolution used at every level" in guide
        cells = z.read(f"{root}/data/cells.csv").decode().splitlines()
        assert len(cells) - 1 == len(s.cells), (len(cells), len(s.cells))
        steps = {x["step"] for x in M["checks"]}
        assert {"filter", "annotation", "embedding", "niches", "analysis de", "analysis composition", "analysis run",
                "analysis coloc", "analysis trajectory_spatial"} <= steps, steps
        for h in ("### Cell-cell communication (`run`)", "### Co-localization (neighborhood enrichment) (`coloc`)",
                  "### Trajectory (spatial axis) (`trajectory_spatial`)", "communication.py", "coloc.py", "trajectory.py#L"):
            assert h in guide, h
        tr = z.read(f"{root}/data/trajectory_spatial.csv").decode().splitlines()
        assert tr[0] == "cell_id,t,bin" and len(tr) - 1 == len(s.cells), tr[:2]
        ct = next(x for x in M["checks"] if x["step"] == "annotation")["expected"]
        assert sum(ct.values()) == len(s.cells)
        col = json.loads(z.read(f"{root}/data/colors.json"))
        assert set(col["cell_type"]) == set(types) and col["niche"] and col["leiden"], col.keys()
        print(f"guide (notebook {nb}): {len(names)} files, {len(r.content) // 1024} KB, {len(guide.splitlines())} lines")
    # two builds at once (a second click): both whole, neither breaks the other
    import threading
    from stagent import reproguide
    outs, errs = [], []
    def one():
        try:
            outs.append(reproguide.build(Session(s.id), v, False))
        except Exception as e:
            errs.append(e)
    th = [threading.Thread(target=one) for _ in range(2)]
    [t.start() for t in th]; [t.join() for t in th]
    assert not errs and len(set(outs)) == 2, (errs, outs)
    for o in outs:
        assert zipfile.ZipFile(o).testzip() is None and f"demo_v{v}_reproduce/REPRODUCE.md" in zipfile.ZipFile(o).namelist()
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
