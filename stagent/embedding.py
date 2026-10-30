"""PCA, UMAP and Leiden clusters of a sample's cells, from expression only (independent of any annotation), stored
with the sample (embedding/) and shared by all its annotations.

Fast path: the PCA works on the sparse matrix; the kNN graph comes from markers.knn (exact on a GPU, else
pynndescent) and is handed to UMAP as precomputed neighbors, so UMAP skips its own neighbor search; UMAP starts from
PC1/PC2 instead of a spectral layout and optimizes on all cores (n_jobs=-1, hence not bit-for-bit reproducible).
Leiden runs on the same kNN graph.
"""
import json, time
import numpy as np
from . import config, jobs, markers, model


FILES = ("umap.npy", "pca.npy", "leiden.npy", "meta.json")  # an embedding's files


def _dir(sample_id):
    d = config.SAMPLES / sample_id / "embedding"
    d.mkdir(exist_ok=True)
    return d


def run(s, n_pcs=30, n_neighbors=15, min_dist=0.3, resolution=1.0, progress=None):
    say = progress or (lambda f, m: None)
    t0 = time.time()
    say(0.05, "PCA")
    Z = markers.pca(model.lognorm(s.counts), n_pcs)[:, ::-1].copy()  # PC1 first
    say(0.2, f"{n_neighbors} nearest neighbors of {len(Z):,} cells")
    dist, nbr = markers.knn(Z, n_neighbors)
    say(0.4, "Leiden clusters")
    lab = markers._leiden(markers.knn_graph(dist, nbr), nbr, resolution, 0)
    say(0.5, "UMAP (the first run also compiles it, about half a minute)")
    import warnings
    with warnings.catch_warnings():  # "Tensorflow not installed; ParametricUMAP will be unavailable"
        warnings.simplefilter("ignore")
        import umap
    n = len(Z)
    idx = np.c_[np.arange(n), nbr].astype(np.int64)  # UMAP expects each cell as its own first neighbor
    dd = np.c_[np.zeros(n, np.float32), dist].astype(np.float32)
    init = Z[:, :2].astype(np.float64)
    init = 10 * (init - init.mean(0)) / max(float(np.abs(init - init.mean(0)).max()), 1e-9)
    Y = umap.UMAP(n_neighbors=idx.shape[1], min_dist=float(min_dist), precomputed_knn=(idx, dd, None), init=init,
                  n_jobs=-1, low_memory=True, verbose=False).fit_transform(Z)
    d = _dir(s.sample_id)
    meta = {"run_id": time.strftime("%Y%m%d-%H%M%S"), "n_pcs": int(n_pcs), "n_neighbors": int(n_neighbors),
            "min_dist": float(min_dist), "resolution": float(resolution), "n_clusters": int(lab.max()) + 1,
            "seconds": round(time.time() - t0, 1)}
    with jobs.replacing(*(d / f for f in FILES)):  # a stopped run puts the last embedding back
        np.save(d / "umap.npy", Y.astype(np.float32))
        np.save(d / "pca.npy", Z[:, :2].astype(np.float32))
        np.save(d / "leiden.npy", lab.astype(np.int32))
        json.dump(meta, open(d / "meta.json", "w"), indent=1)
    say(1.0, "done")
    return meta


def current(sample_id):
    d = config.SAMPLES / sample_id / "embedding"
    return json.load(open(d / "meta.json")) if (d / "meta.json").exists() else None


def coords(sample_id, kind="umap"):
    if kind not in ("umap", "pca"):
        raise ValueError("kind must be umap or pca")
    return np.load(config.SAMPLES / sample_id / "embedding" / f"{kind}.npy")


def leiden(sample_id):
    f = config.SAMPLES / sample_id / "embedding" / "leiden.npy"
    return np.load(f) if f.exists() else None
