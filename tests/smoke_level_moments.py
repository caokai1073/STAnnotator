"""Whole-transcriptome levels (more than MAX_PCA_GENES genes): _level_moments (genes in groups straight from the
counts, the level's cells picked by binary search) gives the gene means and means of squares of _moments on the
level's log-normalized matrix to the bit, and _level_matrix the same PCA input as the old way (the most variable genes
of lognorm(counts)[idx]); every cell, a subset, cells without transcripts, small groups (many of them)."""
import numpy as np, scipy.sparse as sp
from stagent import markers, model, csrutil

rng = np.random.default_rng(3)
n, g = 3000, markers.MAX_PCA_GENES + 1500
X = sp.random(n, g, density=0.03, format="csr", random_state=3, dtype=np.float64)
X.data = np.ceil(X.data * 6)
X = X.astype(np.float32)
X[5] = 0  # a cell without transcripts
X.eliminate_zeros()
X.sum_duplicates()
assert X.has_canonical_format
for name, idx in (("all", np.arange(n)), ("subset", np.sort(rng.choice(n, 1900, replace=False))),
                  ("few", np.sort(rng.choice(n, 40, replace=False)))):
    every = len(idx) == n
    Xn = model.lognorm(X) if every else model.lognorm(X[idx])
    mu0, e0 = markers._moments(Xn)
    for group in (5_000, 200_000, markers.MOMENT_GROUP):
        mu1, e1 = markers._level_moments(X, idx, every, group=group)
        assert np.array_equal(mu0, mu1) and np.array_equal(e0, e1), (name, group)
    keep = np.sort(np.argsort(-(e0 - mu0 ** 2))[:markers.N_HVG])
    old = csrutil.columns_f32(Xn, keep)
    new = markers._level_matrix(X, idx)
    assert old.shape == new.shape and np.array_equal(old.indptr, new.indptr) and np.array_equal(old.indices, new.indices) \
        and np.array_equal(old.data, new.data), name
print("smoke_level_moments ok")
