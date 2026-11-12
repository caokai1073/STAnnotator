"""The speed-ups of 0.2.56 give the same numbers as before, and the app's exported cell positions work as a registration file:
1. fastpar.col_means = scipy's X[rows].mean(0) (and of expm1, and the detection rate), bit for bit;
2. fastpar.lr_edges = numpy's products and np.bincount, bit for bit; pmap keeps the order of its items;
3. markers.pca / knn: a second call on the same input returns the kept result, equal to the first;
4. an import whose registration is a cell positions file puts the cells it lists where it says (a correction of the
   placement: replayed by replay_placement, taken back by Undo correction), also when the image comes later.
Run: python tests/smoke_fastpar.py  (synthetic data, in a temporary data folder removed at the end)."""
import os, shutil, tempfile
from pathlib import Path

HOME = Path(tempfile.mkdtemp(prefix="stagent_fastpar_"))
os.environ["STAGENT_HOME"] = str(HOME)
import numpy as np, pandas as pd, scipy.sparse as sp  # noqa: E402
from stagent import fastpar, markers, importer  # noqa: E402

rng = np.random.default_rng(3)
try:
    # ---- 1. column means
    X = sp.random(5000, 300, density=0.08, format="csr", random_state=4, dtype=np.float32)
    X.data = np.log1p(X.data * 50).astype(np.float32)
    perm = rng.permutation(len(X.indices))  # unsorted column indices within rows, as lognorm's output can have
    for r in range(X.shape[0]):
        a, b = X.indptr[r], X.indptr[r + 1]
        X.indices[a:b], X.data[a:b] = X.indices[a:b][::-1].copy(), X.data[a:b][::-1].copy()
    X.has_sorted_indices = False
    E = np.expm1(X.data)
    for rows in (np.sort(rng.choice(5000, 37, replace=False)), np.arange(5000), np.flatnonzero(rng.random(5000) < 0.7)):
        sub = X[rows]
        e = sub.copy()
        e.data = np.expm1(e.data)
        ref = (np.asarray(sub.mean(0)).ravel(), np.asarray(e.mean(0)).ravel(), sub.getnnz(axis=0) / len(rows))
        got = fastpar.col_means(X, rows, E)
        assert all(np.array_equal(u, v) and u.dtype == v.dtype for u, v in zip(ref, got)), "col_means differs from scipy"
    print("col_means: the same as scipy, bit for bit")

    # ---- 2. ligand-receptor products and their sums per group pair; pmap
    n, m, K = 4000, 60000, 7
    src, dst = rng.integers(0, n, m).astype(np.int32), np.repeat(np.arange(n, dtype=np.int32), m // n)
    w = rng.random(m).astype(np.float32)
    grp = rng.integers(0, K, n)
    pair = grp[src] * K + grp[dst]
    L, R = rng.random(n).astype(np.float32) * (rng.random(n) < 0.4), rng.random(n).astype(np.float32)
    v, wv, sums = fastpar.lr_edges(src, dst, w, pair, L, R, K * K)
    v0 = L[src] * R[dst]
    assert np.array_equal(v, v0) and np.array_equal(wv, w * v0) and np.array_equal(sums, np.bincount(pair, w * v0, K * K))
    items = list(range(50))
    assert fastpar.pmap(lambda i: (i, float(np.sum(np.arange(i * 1000.0)))), items) == [(i, float(np.sum(np.arange(i * 1000.0)))) for i in items]
    print("lr_edges: the same as numpy and np.bincount; pmap: in order")

    # ---- 3. PCA and neighbors kept by the content of their input
    C = sp.random(3000, 80, density=0.2, format="csr", random_state=5, dtype=np.float32)
    Z1 = markers.pca(C, 10)
    Z1[:] = 0  # (the caller's copy: the kept one is untouched)
    Z2, Z3 = markers.pca(C, 10), markers._pca(sp.csr_matrix(C, dtype=np.float32), 10)
    assert np.array_equal(Z2, Z3), "the kept PCA differs from a new one"
    d1, i1 = markers.knn(Z2, 10, "exact")
    d2, i2 = markers.knn(Z2.copy(), 10, "exact")
    assert np.array_equal(d1, d2) and np.array_equal(i1, i2)
    print("pca / knn: kept results equal new ones")

    # ---- 4. cell positions as the registration file
    import tifffile
    h5 = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish" / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"
    r = importer.read_h5ad(str(h5), None, None, None, None, None, None)
    ids, xy = np.asarray(r["cell_id"]).astype(str), np.asarray(r["xy"], float)
    ncell = len(ids)
    img = HOME / "he.tif"
    tifffile.imwrite(img, np.full((1150, 1300, 3), 230, np.uint8), photometric="rgb")
    A = np.array([[0.5, 0.0, 600.0], [0.0, 0.5, 600.0]])  # (the cells' um span -1083..1298)
    p = xy @ A[:, :2].T + A[:, 2]
    moved = rng.random(ncell) < 0.3  # a third of the cells off the affine (as local corrections or a moved piece)
    p[moved] += rng.normal(0, 4, (moved.sum(), 2))
    listed = rng.random(ncell) < 0.9  # the file lists most cells (another sample's export, after a filter)
    pos = pd.DataFrame({"cell_id": ids, "x_um": xy[:, 0].round(3), "y_um": xy[:, 1].round(3), "px": p[:, 0].round(3),
                        "py": p[:, 1].round(3)})[listed]
    reg = HOME / "s_cell_positions.csv"
    pos.to_csv(reg, index=False)
    assert importer.read_positions(reg) is not None and importer.read_positions(h5) is None
    importer.import_sample("pos", str(h5), he=str(img), registration=str(reg), progress=lambda *a: None)
    m = importer._meta("pos")
    cells = pd.read_parquet(importer.config.SAMPLES / "pos" / "cells.parquet").set_index("cell_id")
    got = cells.loc[pos.cell_id, ["px", "py"]].to_numpy()
    assert np.abs(got - pos[["px", "py"]].to_numpy()).max() < 1e-9, "listed cells not where the file puts them"
    Afit = np.asarray(m["placement_affine"])
    unl = ids[~listed]
    expect = xy[~listed] @ Afit[:, :2].T + Afit[:, 2]
    assert np.abs(cells.loc[unl, ["px", "py"]].to_numpy() - expect).max() < 1e-6, "unlisted cells not on the fitted affine"
    assert any("positions" in f for f in m["registration_fixes"]) and "cell positions" in m["registration_by"]
    first, now = importer.replay_placement(m, cells[["x", "y"]].to_numpy(float), cells.index.to_numpy())
    assert np.abs(now - cells[["px", "py"]].to_numpy()).max() < 1e-6, "replay differs from the positions"
    from stagent import qcfilter  # a filter keeps every remaining cell where it was, and brings cells back in place
    qcfilter.apply(["pos"], {"min_counts": int(np.percentile(cells.total_counts, 30))})
    f1 = pd.read_parquet(importer.config.SAMPLES / "pos" / "cells.parquet").set_index("cell_id")
    assert len(f1) < len(cells) and np.abs(f1[["px", "py"]].to_numpy() - cells.loc[f1.index, ["px", "py"]].to_numpy()).max() < 1e-6
    qcfilter.apply(["pos"], {})
    f2 = pd.read_parquet(importer.config.SAMPLES / "pos" / "cells.parquet").set_index("cell_id")
    assert len(f2) == len(cells) and np.abs(f2.loc[cells.index, ["px", "py"]].to_numpy() - cells[["px", "py"]].to_numpy()).max() < 1e-6
    m = importer._meta("pos")
    importer.undo_registration_fix("pos")
    back = pd.read_parquet(importer.config.SAMPLES / "pos" / "cells.parquet").set_index("cell_id")
    assert np.abs(back[["px", "py"]].to_numpy() - (back[["x", "y"]].to_numpy() @ Afit[:, :2].T + Afit[:, 2])).max() < 1e-6
    print(f"cell positions as registration: {listed.sum()} cells placed as in the file, the rest on the fitted affine; "
          "kept through a filter and back, and undone")
    # the image later: the positions wait, then apply
    importer.import_sample("pos2", str(h5), registration=str(reg), progress=lambda *a: None)
    assert importer._meta("pos2").get("registration_positions")
    importer.set_image("pos2", str(img), progress=lambda *a: None)
    c2 = pd.read_parquet(importer.config.SAMPLES / "pos2" / "cells.parquet").set_index("cell_id")
    assert np.abs(c2.loc[pos.cell_id, ["px", "py"]].to_numpy() - pos[["px", "py"]].to_numpy()).max() < 1e-9
    print("cell positions given before the image: applied when the image is added")
    print("ok")
finally:
    shutil.rmtree(HOME, ignore_errors=True)
