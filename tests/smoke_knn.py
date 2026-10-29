"""The exact neighbor search gives the same neighbors whatever computes it (GPU, torch on the CPU, numpy), equal to
a float64 brute force; and the annotation time estimate covers every step.
Run: python tests/smoke_knn.py [sample_id] (default: synthetic cells)."""
import sys, time
import numpy as np, scipy.sparse as sp
from stagent import config, markers, model, timing

config.RECORD_TIMINGS = False
if len(sys.argv) > 1:
    X = sp.load_npz(config.SAMPLES / sys.argv[1] / "counts.npz").tocsr()
    X = X[np.random.default_rng(1).choice(X.shape[0], min(X.shape[0], 40000), replace=False)]
    Z = markers.pca(model.lognorm(X))
else:
    rng = np.random.default_rng(0)
    Z = (rng.standard_normal((20000, 30)) * rng.random(30)).astype(np.float32)
    Z[1:4] = Z[0]  # identical cells: ties go to the lower row number
k, kk = 15, 15 + 1 + markers.KNN_EXTRA


def rerank(cand):
    Z64, rows = Z.astype(np.float64), np.arange(len(Z))
    d = ((Z64[cand] - Z64[rows, None]) ** 2).sum(-1)
    d[cand == rows[:, None]] = np.inf
    return np.take_along_axis(cand, np.lexsort((cand, d), axis=1)[:, :k], 1)


backends = {"numpy": lambda: markers._topk_numpy(Z, kk)}
if markers.has_torch():
    backends["torch, CPU"] = lambda: markers._topk_torch(Z, kk, "cpu")
    if config.device() != "cpu":
        backends[f"torch, {config.device()}"] = lambda: markers._topk_torch(Z, kk, config.device())
got = {}
for name, f in backends.items():
    t = time.time()
    got[name] = rerank(f())
    print(f"{name}: {time.time() - t:.1f} s")
first = next(iter(got.values()))
for name, nb in got.items():
    assert (nb == first).all(), f"{name} differs in {(nb != first).any(1).sum()} rows"
_, nb = markers.knn(Z, k, "exact")
assert (nb == first).all()
r = np.arange(300)
D = ((Z[r, None].astype(np.float64) - Z[None].astype(np.float64)) ** 2).sum(-1)
D[r, r] = np.inf
true = np.lexsort((np.broadcast_to(np.arange(len(Z)), D.shape), D), axis=1)[:, :k]
assert (true == nb[:300]).all(), "exact neighbors differ from the float64 brute force"
print(f"exact neighbors identical on {len(got)} backends and to the brute force ({len(Z):,} cells)")

e = timing.estimate(len(Z), 480, None, "auto", "exact")
assert e["seconds"] > 0 and len(e["steps"]) == 6, e
print("estimate:", e)
