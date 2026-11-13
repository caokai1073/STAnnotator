"""Spatial information in the marker annotation, on the bundled mouse colon example:
1. a new annotation types the cells the marker step could not type with their neighbors too (scheme "neighbors"):
   only those cells can differ from the same annotation without it; the cells the marker step typed do not;
2. the neighbor rounds stop (no change, or a cell flipping between two types);
3. spillover (an option): a gene in a cell below a tenth of its neighbors' mean count is not counted; the marker
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
    s = Session.create("M", method="markers", markers_text=text, knn="exact")
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
    s2 = Session.create("M", method="markers", markers_text=text, knn="exact", spillover=True)
    assert s2.scheme["spillover"] is True
    lab2 = np.asarray(s2.labels())
    print(f"with spillover: {(lab2 != lab).sum():,} cells labeled otherwise; Unassigned {(lab2 == markers.UNASSIGNED).sum()}")
    print("ok")
finally:
    shutil.rmtree(HOME, ignore_errors=True)
