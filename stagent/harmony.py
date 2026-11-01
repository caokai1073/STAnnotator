"""Harmony (Korsunsky et al. 2019) in NumPy, used on Windows only: harmonypy 2.x is compiled and has no Windows
wheels, and the pure-Python harmonypy 0.0.10 takes many minutes on a few hundred thousand cells. macOS and Linux keep
harmonypy.

It follows harmonypy 2.1.0 (R harmony2) with its defaults, for one batch variable: soft k-means on the cosine-normalized
PCA with a diversity penalty (theta 2, sigma 0.1, min(N / 30, 100) clusters, 4 assignment passes in blocks of 5% of the
cells per round), then per cluster a ridge regression of the original PCA on the batches, with the penalty estimated as
alpha * E (alpha 0.2), whose batch terms are subtracted; at most 10 rounds, until the objective changes by less than 1%.
The random choices (initial centroids, block order) differ from harmonypy's, so results agree closely, not exactly.
"""
import numpy as np, scipy.sparse as sp


def _unit(X):
    """Rows of X scaled to length 1."""
    n = np.linalg.norm(X, axis=1, keepdims=True)
    return X / np.maximum(n, 1e-12)


def _softmax(L):
    L -= L.max(1, keepdims=True)
    np.exp(L, out=L)
    L /= L.sum(1, keepdims=True)
    return L


def _kmeans_init(X, K, rng):
    """K centroids from the cells as harmony2 picks them: for each, a random cell, replaced by a cell drawn with
    probability proportional to its cosine distance from it (no cell twice); then 10 Lloyd steps."""
    N = len(X)
    Y = X[rng.integers(0, N, K)].copy()
    chosen = np.zeros(N, bool)
    for i in range(K):
        w = np.abs(2 * (1 - X @ Y[i]))
        w[chosen] = 0
        s = w.sum()
        j = rng.choice(N, p=w / s) if s > 0 else rng.choice(np.flatnonzero(~chosen))
        chosen[j] = True
        Y[i] = X[j]
    for _ in range(10):
        a = np.argmin((Y * Y).sum(1)[None, :] - 2 * (X @ Y.T), axis=1)
        n = np.bincount(a, minlength=K)
        S = sp.csr_matrix((np.ones(N), (a, np.arange(N))), shape=(K, N)) @ X.astype(np.float64)
        keep = n > 0
        Y_new = Y.copy()
        Y_new[keep] = (S[keep] / n[keep, None]).astype(Y.dtype)
        if np.array_equal(Y_new, Y):
            break
        Y = Y_new
    return _unit(Y)


def harmonize(Z, batch, theta=2.0, sigma=0.1, nclust=None, alpha=0.2, block_size=0.05, max_iter=10,
              max_iter_kmeans=4, epsilon=1e-2, batch_prop_cutoff=1e-5, random_state=0, progress=None):
    """Z (cells x PCs) corrected for the batches in batch (one label per cell); returns float32 cells x PCs.
    progress(i, max_iter) is called after each round."""
    Z0 = np.ascontiguousarray(Z, dtype=np.float32)
    N, d = Z0.shape
    codes = np.unique(np.asarray(batch), return_inverse=True)[1].reshape(-1)
    n_b = np.bincount(codes)
    B = len(n_b)
    if B < 2:
        return Z0.copy()
    K = int(nclust or min(round(N / 30.0), 100))
    K = max(1, min(K, N))
    Pr = n_b / N
    theta = np.broadcast_to(np.asarray(theta, np.float64), (B,))
    sig = np.broadcast_to(np.asarray(sigma, np.float64), (K,)).astype(np.float32)
    members = [np.flatnonzero(codes == b) for b in range(B)]
    onehot = np.zeros((N, B), np.float32)
    onehot[np.arange(N), codes] = 1
    rng = np.random.default_rng(random_state)

    def assign(Zc, Y):
        """Cosine distances to the centroids, the assignments without the diversity penalty, and -dist / sigma."""
        dist = 2 * (1 - Zc @ Y.T)
        A = -dist / sig
        return dist, _softmax(A.copy()), A

    def tallies(R):
        O = (onehot.T @ R).T.astype(np.float64)        # K x B
        E = np.outer(R.sum(0, dtype=np.float64), Pr)  # K x B
        return O, E

    def objective(R, dist, O, E):
        err = float((R * dist).sum(dtype=np.float64))
        ent = float((sig * (R * np.log(np.maximum(R, 1e-30)))).sum(dtype=np.float64))
        cross = float((np.log((O + E + 1) / (2 * E + 1)) * theta[None, :] * sig[:, None].astype(np.float64) * O).sum())
        return (err + ent + cross) * 2000.0 / N

    Zc = _unit(Z0)
    Y = _kmeans_init(Zc, K, rng)
    dist, R, A = assign(Zc, Y)
    O, E = tallies(R)
    obj = [objective(R, dist, O, E)]
    n_blocks = int(np.ceil(1 / block_size))
    per_block = max(1, int(N * block_size))
    Zcorr = Z0
    for it in range(1, max_iter + 1):
        # clustering: assignments with the diversity penalty, cells updated a block at a time
        if len(obj) > 1:
            Zc = _unit(Zcorr)
            dist, R, A = assign(Zc, Y)
            O, E = tallies(R)
        for _ in range(max_iter_kmeans):
            order = rng.permutation(N)
            for i in range(n_blocks):
                idx = order[i * per_block: N if i == n_blocks - 1 else (i + 1) * per_block]
                if not len(idx):
                    break
                cb = codes[idx]
                Rb = R[idx]
                ob = onehot[idx]
                O -= (ob.T @ Rb).T
                E -= np.outer(Rb.sum(0, dtype=np.float64), Pr)
                log_div = (theta[:, None] * (np.log(2 * E.T + 1) - np.log(O.T + E.T + 1))).astype(np.float32)  # B x K
                Rb = A[idx]
                Rb += log_div[cb]
                Rb = _softmax(Rb)
                R[idx] = Rb
                O += (ob.T @ Rb).T
                E += np.outer(Rb.sum(0, dtype=np.float64), Pr)
        obj.append(objective(R, dist, O, E))

        # correction: per cluster, a ridge regression of the PCA on the batches (penalty alpha * E). The system
        # [[total, O'], [O, diag(O + alpha E)]] is solved by eliminating its diagonal block: the intercept W0 (the new
        # centroid) is the Schur complement's solution, each kept batch's term follows from it.
        kept =(O / n_b[None, :] > batch_prop_cutoff)
        kept &= kept.sum(1, keepdims=True) > 1
        RZ = np.stack([R[m].T.astype(np.float64) @ Z0[m].astype(np.float64) for m in members])  # B x K x d
        Ok = np.where(kept, O, 0.0)
        D = np.where(kept, O + alpha * E, 1.0)
        s = Ok.sum(1) - (Ok * Ok / D).sum(1)                                   # K
        rhs = np.where(kept.T[:, :, None], RZ, 0.0)                            # kept batches only
        rhs0 = rhs.sum(0) - ((Ok / D).T[:, :, None] * rhs).sum(0)              # K x d
        active = kept.any(1)
        W0 = np.where(active[:, None], rhs0 / np.where(active, s, 1.0)[:, None], 0.0)
        M = np.where(kept.T[:, :, None], (rhs - Ok.T[:, :, None] * W0[None]) / D.T[:, :, None], 0.0)  # B x K x d
        Zcorr = Z0.copy()
        for b, m in enumerate(members):
            Zcorr[m] -= (R[m] @ M[b].astype(np.float32))
        Y = Y.astype(np.float64)
        Y[active] = W0[active]
        Y = _unit(Y).astype(np.float32)

        if progress:
            progress(it, max_iter)
        old, new = obj[-2], obj[-1]
        if np.isfinite(old) and np.isfinite(new) and old != 0 and 0 <= (old - new) / abs(old) < epsilon:
            break
    return Zcorr
