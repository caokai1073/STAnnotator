"""Which step of the exact neighbor search on the Apple GPU (markers._topk_mlx) goes wrong, and on how many rows
(developer tool, macOS with MLX).

    STAGENT_HOME=<data folder> .venv/bin/python tools/diag_mlx_knn.py --sample <sample id> [--lineage-rows N]
    .venv/bin/python tools/diag_mlx_knn.py --npy Z.npy            (the PCA scores of a search: cells x components)

With --sample: the root level of an annotation (log-normalized counts, PCA), as markers._graph searches it. The rows
are centered and deduplicated as markers._exact does, then searched block by block as _topk_mlx does, and every GPU
step is compared with numpy on the same numbers:
  1. the distances (the product on the GPU against numpy's, float32)
  2. the minimum of each group of g columns (the reduction)
  3. the k groups with the smallest minima (argpartition on the GPU)
  4. the distances of the k*g columns of those groups (take_along_axis)
  5. the k smallest of them (argpartition on the GPU)
  6. the result: rows whose k-th distance is larger than the true k-th (over all columns)
Nothing is written. Prints the number of rows wrong at each step and a few examples."""
import argparse, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from stagent import markers  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--sample")
ap.add_argument("--npy")
ap.add_argument("--k", type=int, default=markers.N_NEIGHBORS + 1 + markers.KNN_EXTRA, help="candidates searched (kx)")
ap.add_argument("--max-rows", type=int, default=0, help="check only the first N rows (0: all)")
a = ap.parse_args()

import mlx.core as mx  # noqa: E402
print("mlx", mx.__version__, mx.default_device())
if a.npy:
    Z = np.load(a.npy)
else:
    import scipy.sparse as sp
    from stagent import config
    counts = sp.load_npz(config.SAMPLES / a.sample / "counts.npz")
    t = time.time()
    from stagent.model import lognorm
    Z = markers.pca(lognorm(counts))
    print(f"PCA {Z.shape} in {time.time() - t:.0f}s")
Z64 = np.asarray(Z, np.float64)
Zc = (Z64 - Z64.mean(0)).astype(np.float32)
Zu = np.unique(Zc, axis=0)
X = Zu if len(Zc) - len(Zu) > max(1000, len(Zc) // 100) else Zc
n, dim = X.shape
k = min(a.k, n)
g = max(1, int(round(np.sqrt(n / k))))
ng = -(-n // g)
print(f"rows searched {n} (of {len(Zc)}), dim {dim}, k {k}, g {g}, groups {ng}, sorted per row: {ng} then {k * g}")
if g == 1 or ng <= k:
    sys.exit("few rows: the direct path (whole rows sorted); nothing grouped to check")
w = -(-(dim + 1) // 16) * 16
A = np.zeros((ng * g, w), np.float32)
A[:n, :dim] = X
A[:n, dim] = (X ** 2).sum(1)
A[n:, dim] = 1e30
L = np.zeros((n, w), np.float32)
np.multiply(X, -2, out=L[:, :dim])
L[:, dim] = 1
rows = int(min(n, max(64, (256 * 2 ** 20) // (4 * len(A)))))
last = min(n, a.max_rows) if a.max_rows else n
bad = {s: [] for s in range(1, 7)}
maxdiff = 0.0
Lm, At = mx.array(L), mx.array(A).T
offs = mx.arange(g, dtype=mx.uint32)
t = time.time()
for i in range(0, last, rows):
    j = min(i + rows, last)
    P = Lm[i:j] @ At
    Pn = np.array(P)                                    # the GPU's distances, as numpy
    ref = L[i:j] @ A.T                                  # numpy's
    d = np.abs(Pn[:, :n] - ref[:, :n]).max(1)
    maxdiff = max(maxdiff, float(d.max()))
    G = P.reshape(-1, ng, g).min(axis=2)
    Gn, Gref = np.array(G), Pn.reshape(-1, ng, g).min(axis=2)
    bad[2] += list(i + np.flatnonzero((Gn != Gref).any(1)))
    grp = mx.argpartition(G, k - 1, axis=1)[:, :k].astype(mx.uint32)
    grpn = np.array(grp).astype(np.int64)
    kth = np.partition(Gref, k - 1, axis=1)[:, k - 1]
    got = np.take_along_axis(Gref, grpn, 1)
    bad[3] += list(i + np.flatnonzero((got.max(1) > kth) | (np.sort(grpn, 1)[:, 1:] == np.sort(grpn, 1)[:, :-1]).any(1)))
    cols = (grp[:, :, None] * g + offs).reshape(-1, k * g)
    v = mx.take_along_axis(P, cols, axis=1)
    coln, vn = np.array(cols).astype(np.int64), np.array(v)
    bad[4] += list(i + np.flatnonzero((vn != np.take_along_axis(Pn, coln, 1)).any(1)))
    sel = mx.argpartition(v, k - 1, axis=1)[:, :k]
    seln = np.array(sel).astype(np.int64)
    vk = np.partition(vn, k - 1, axis=1)[:, k - 1]
    bad[5] += list(i + np.flatnonzero((np.take_along_axis(vn, seln, 1).max(1) > vk)
                                      | (np.sort(seln, 1)[:, 1:] == np.sort(seln, 1)[:, :-1]).any(1)))
    c = np.take_along_axis(coln, seln, 1)
    true_k = np.partition(Pn[:, :n], k - 1, axis=1)[:, k - 1]
    bad[6] += list(i + np.flatnonzero(np.take_along_axis(Pn, c, 1).max(1) > true_k))
    print(f"\r{j}/{last} rows, {time.time() - t:.0f}s", end="", flush=True)
print()
print(f"1. distances: largest |GPU - numpy| {maxdiff:.3g} (float32 rounding: about 1e-5 of |z|^2)")
names = {2: "group minima (reduction)", 3: "k smallest groups (argpartition)", 4: "gathered distances (take_along_axis)",
         5: "k smallest columns (argpartition)", 6: "result: k-th neighbor farther than the true k-th"}
for s in range(2, 7):
    print(f"{s}. {names[s]}: {len(bad[s])} rows wrong" + (f", e.g. rows {bad[s][:8]}" if bad[s] else ""))
