"""Spatial information in the marker annotation, on the bundled mouse colon example:
1. a new annotation types the cells the marker step could not type with their neighbors too (scheme "neighbors"):
   only those cells can differ from the same annotation without it; the cells the marker step typed do not;
2. the neighbor rounds stop (no change, or a cell flipping between two types);
3. spillover (on by default): a gene in a cell below a tenth of its neighbors' mean count is not counted; the marker
   step runs with it and stores the option; the neighbors of a project's slices stay within each slice.
Run: python tests/smoke_spatial.py  (a temporary data folder removed at the end)."""
import os, shutil, tempfile
from pathlib import Path

HOME = Path(tempfile.mkdtemp(prefix="stagent_spatial_"))
os.environ["STAGENT_HOME"] = str(HOME)
import numpy as np, scipy.sparse as sp  # noqa: E402
from stagent import importer, markers  # noqa: E402
from stagent.session import Session  # noqa: E402

EX = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish"
try:
    importer.import_sample("M", str(EX / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    text = (EX / "mouse_colon_merfish_markers.txt").read_text()
    s = Session.create("M", method="markers", markers_text=text, knn="exact", spillover=False)
    assert s.scheme["neighbors"] is True and s.scheme["spillover"] is False
    lab = np.asarray(s.labels())
    res = s.marker_result()
    untyped = np.asarray(res["labels"]) == markers.UNASSIGNED
    st = s.state()
    logp = s.base(st)["logp"]
    tree = s._marker_tree([])
    import stagent.model as model  # the same restriction without the neighbors
    m = s._init_model(st, lambda *a: None)
    plain = markers.restrict(model.predict(m, s._model_counts()), res, tree, s.types, (), counts=s.counts,
                             panel=s.panel).argmax(1)
    withn = logp.argmax(1)
    changed = withn != plain
    assert changed.any(), "the neighbors changed nothing"
    assert not (changed & ~untyped).any(), "a cell the marker step typed changed"
    print(f"{untyped.sum():,} cells the marker step could not type; the neighbors changed {changed.sum():,} of them")

    # the rounds stop: a tiny case where two cells would flip forever
    # two pairs of cells, each cell leaning a little to the type its neighbor has: every round they swap
    lp = np.log(np.array([[0.51, 0.49], [0.49, 0.51], [0.51, 0.49], [0.49, 0.51]], np.float32))
    hist = markers.neighbor_rounds(lp, [0, 1, 2, 3], np.array([[0, 0], [1, 0], [10, 0], [11, 0]], float), k=1,
                                   rounds=50)
    assert hist[0] > 0 and len(hist) <= 3, hist
    print("neighbor rounds:", hist)

    # spillover: the detection matrix and the marker step with it
    xy = s.cells[["x", "y"]].to_numpy(float)
    D = markers.spill_detected(s.counts, xy)
    plain_det = (sp.csr_matrix(s.counts) > 0)
    assert D.nnz < plain_det.nnz and not (D > plain_det).nnz, "spillover must only remove detections"
    print(f"spillover leaves out {plain_det.nnz - D.nnz:,} of {plain_det.nnz:,} detections")
    g = np.repeat([0, 1], [len(xy) // 2, len(xy) - len(xy) // 2])
    nb = markers._spatial_neighbors(np.vstack([xy[: len(xy) // 2], xy[: len(xy) - len(xy) // 2]]), 6, g)
    assert (g[nb] == g[:, None]).all(), "neighbors crossed slices"
    s2 = Session.create("M", method="markers", markers_text=text, knn="exact")  # spillover on by default
    assert s2.scheme["spillover"] is True
    lab2 = np.asarray(s2.labels())
    print(f"with spillover: {(lab2 != lab).sum():,} cells labeled otherwise; Unassigned {(lab2 == markers.UNASSIGNED).sum()}")

    # an annotation made before 0.2.57 (no "neighbors" in its scheme) keeps the rules of then when it is edited:
    # a one-cell correction changes about that cell, not the thousands the new rules would re-type
    import json
    old = Session.create("M", method="markers", markers_text=text, knn="exact", spillover=False)
    info = json.load(open(old.dir / "session.json"))
    for k in ("neighbors", "spillover"):
        info["scheme"].pop(k, None)
    json.dump(info, open(old.dir / "session.json", "w"))
    old = Session(old.id)
    v1 = np.asarray(old.labels())
    un1 = int((v1 == markers.UNASSIGNED).sum())
    assert un1 > 1000, f"the old rules leave the untyped cells Unassigned ({un1})"
    target = s.types[0] if s.types[0] != markers.UNASSIGNED else s.types[1]
    v = old.label_cells(target, cells=[int(np.where(v1 != target)[0][0])])
    n_changed = old.diff(v - 1, v)["n_label_changes"]
    assert n_changed < 50, f"one correction of an old annotation changed {n_changed} labels"
    print(f"old annotation: {un1:,} Unassigned kept; a one-cell correction changed {n_changed} label(s)")

    # re-clustering one cluster: the types its level has elsewhere are not searched for as rare types in it
    rec = max((r for r in s.marker_result()["records"] if r["parent"] is None and r.get("assigned")),
              key=lambda r: r["n_cells"])
    assigned_elsewhere = {r.get("assigned") for r in s.marker_result()["records"] if r["parent"] is None and r["id"] != rec["id"]}
    res2 = markers.rerun(s.marker_result(), s.counts, s.panel, s._marker_tree([]), rec["id"], knn="exact")
    before = {r["id"] for r in s.marker_result()["records"]}
    new_rare = {r["assigned"] for r in res2["records"] if r.get("rare") and r["parent"] is None and r["id"] not in before}
    assert not (new_rare & assigned_elsewhere), f"re-clustering looked for types found elsewhere: {new_rare & assigned_elsewhere}"
    print(f"re-clustering {rec['id']} ({rec['assigned']}): no rare search for the level's other types")

    # cells of few transcripts are not drawn to a low-quality type whose markers are seldom detected (a 790,000-cell
    # Xenium section gave 59,000 of them to such a "Mast cells" cluster): detection scaled to each cell's depth
    rng = np.random.default_rng(0)
    T, K = 6, 8
    g = [f"g{t}_{k}" for t in range(T) for k in range(K)] + [f"m{k}" for k in range(4)]
    tr = markers.parse("\n".join(f"T{t}: " + " ".join(f"g{t}_{k}+" for k in range(K)) for t in range(T))
                       + "\nMast: " + " ".join(f"m{k}+" for k in range(4)), g)["tree"]
    tys = [f"T{t}" for t in range(T)] + ["Mast", markers.UNASSIGNED]
    def make(n, depth, own=None, p_own=0.5, p_bg=0.02, mast=0.0):
        P = np.full((n, len(g)), p_bg * depth / 400)
        if own is not None:
            P[:, own * K:(own + 1) * K] = p_own
        P[:, T * K:] = mast if mast else p_bg * depth / 400
        return rng.binomial(1, np.clip(P, 0, 1)) * 1.0, np.full(n, depth)
    blocks, lab_, dep, truth = [], [], [], []
    for t in range(T):
        x, d = make(1500, 400, t); blocks.append(x); lab_ += [t] * 1500; dep.append(d)
    x, d = make(3000, 60, None, mast=0.09); blocks.append(x); lab_ += [T] * 3000; dep.append(d)
    for t in range(T):
        x, d = make(500, 60, t, p_own=0.12, p_bg=0.0); keep = x.sum(1) > 0
        blocks.append(x[keep]); dep.append(d[keep]); lab_ += [-1] * int(keep.sum()); truth += [t] * int(keep.sum())
    X = np.vstack(blocks); dep = np.concatenate(dep)
    X = sp.csr_matrix(np.hstack([X, (dep - X.sum(1))[:, None]]))
    rows = np.where(np.array(lab_) == -1)[0]
    _, best = markers.marker_assign(X, g + ["filler"], tr, tys, rows, np.array(lab_))
    right, mast = np.mean(best == np.array(truth)), np.mean(best == T)
    assert mast < 0.05 and right > 0.9, (right, mast)
    print(f"low-depth cells: {right:.0%} to their own type, {mast:.0%} to the low-quality type")
    print("ok")
finally:
    shutil.rmtree(HOME, ignore_errors=True)
