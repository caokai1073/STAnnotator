"""PCA of cells without a single transcript (a cluster of empty cells clustered again): the covariance is 0, the
Lanczos iterations have no direction to start from (ARPACK error -9); the full eigen-decomposition takes over."""
import numpy as np, scipy.sparse as sp
from stagent import markers as M

Z = M._pca(sp.csr_matrix((336, 4624), dtype=np.float32))
assert Z.shape == (336, M.N_PCS) and np.isfinite(Z).all()
g, _ = M._graph(sp.csr_matrix((336, 4624), dtype=np.float32), "exact", None)
assert g.vcount() == 336
rng = np.random.default_rng(1)
Y = sp.csr_matrix(np.log1p(rng.poisson(0.3, (500, 400)).astype(np.float32)))
assert M._pca(Y).shape == (500, M.N_PCS)
print("ok")
