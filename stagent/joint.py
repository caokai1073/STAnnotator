"""Analyzing the slices of a project together, so that labels, clusters and niches mean the same in every slice.

annotate(): the marker step and the expression model run once on the cells of all slices pooled (the same marker
  list, the same Leiden clusters, one model); every slice then gets an annotation (a Session) of its own built on
  that shared result, where regions, priors, corrections and reviews work as for a single slice.
embed(): PCA of the pooled cells (optionally corrected for slice-to-slice differences with Harmony), neighbors,
  Leiden and UMAP on all cells at once; each slice keeps its part (embedding/ of the sample).
niches(): every cell's identity (NMF programs fitted on the pooled cells, or its cell type) averaged over its
  neighbors within its own slice, then one k-means over all slices; each slice's annotation keeps its part.
Comparing the slices (composition, DE genes) is in analysis.py.
"""
import json, time
import numpy as np, pandas as pd, scipy.sparse as sp
from . import config, markers, model, projects


def _dir(pid):
    d = projects.PROJECTS / pid
    d.mkdir(parents=True, exist_ok=True)
    return d


def _pooled(pid):
    """Sample ids, row offsets and the pooled counts of a project's slices (they must share the gene panel)."""
    p = projects.load(pid)
    ids = [s for s in p["samples"] if (config.SAMPLES / s / "meta.json").exists()]
    panels = [(config.SAMPLES / s / "genes.txt").read_text().split() for s in ids]
    if any(g != panels[0] for g in panels[1:]):
        raise ValueError("the slices of this project do not share one gene panel")
    counts = [sp.load_npz(config.SAMPLES / s / "counts.npz").tocsr() for s in ids]
    off = np.cumsum([0] + [c.shape[0] for c in counts])
    return ids, off, sp.vstack(counts).tocsr(), panels[0]


def _slice_mres(res, o, n):
    """The part of a pooled marker-step result that belongs to rows o..o+n (row numbers made local)."""
    members = {}
    for k, rows in res["members"].items():
        r = np.asarray(rows)
        r = r[(r >= o) & (r < o + n)] - o
        if len(r):
            members[k] = r
    return {"labels": res["labels"][o:o + n], "cluster_of": res["cluster_of"][o:o + n], "members": members,
            "records": res["records"], "nodes": res["nodes"], "knn": res.get("knn")}


def latest_sessions(pid, chosen=None):
    """{sample_id: annotation} for the project's slices: the one chosen ({sample_id: session id}, the annotations
    open in the app) when it is an annotation of that slice, else the slice's most recent annotation made in the app."""
    from .session import list_sessions
    ss = list_sessions()
    out = {}
    for s in projects.load(pid)["samples"]:
        mine = ss[(ss.sample_id == s) & ss.method.isin(["markers", "scref", "existing"])] if len(ss) else ss
        if len(mine):
            want = (chosen or {}).get(s)
            out[s] = want if want in set(mine.session_id) else mine.sort_values("created_at").session_id.iloc[-1]
    return out


def annotate(pid, markers_text, resolution=None, carry=False, progress=None, knn=None):
    """One marker annotation of all the project's slices together; returns {sample_id: new session id}.
    knn: the neighbor search, "exact" (default) or "approx"."""
    from .session import Session, knn_mode
    from . import timing
    knn = knn_mode(knn)
    say = progress or (lambda f, m: None)
    ids, off, counts, panel = _pooled(pid)
    with timing.run(counts.shape[0]):
        return _annotate(pid, markers_text, resolution, carry, say, knn, ids, off, counts, panel)


def _annotate(pid, markers_text, resolution, carry, say, knn, ids, off, counts, panel):
    from .session import Session
    previous = latest_sessions(pid) if carry else {}
    from .session import filtered_genes
    parsed = markers.parse(markers_text, panel, filtered_genes(ids[0]))
    if not parsed["ok"]:
        raise ValueError("marker list has errors: " + "; ".join(parsed["errors"][:8]))
    say(0.02, f"marker step on {counts.shape[0]:,} cells of {len(ids)} slices")
    res = markers.annotate(counts, panel, parsed["tree"], lambda f, m: say(0.02 + 0.5 * f, m),
                           resolution=markers.AUTO if resolution in (None, "", markers.AUTO) else float(resolution),
                           knn=knn)
    types = parsed["leaves"]
    index = {t: i for i, t in enumerate(types)}
    lab = np.array([index.get(x, -1) for x in res["labels"]], np.int64)
    say(0.55, "training one expression model on all slices")
    m = model.train_pseudo(counts, lab, len(types), progress=lambda f, msg: say(0.55 + 0.25 * f, msg))
    run = time.strftime("%Y%m%d-%H%M%S")
    out = {}
    for i, sid in enumerate(ids):
        say(0.8 + 0.18 * i / len(ids), f"slice {i + 1} of {len(ids)}: {sid}")
        s = Session.create(sid, "none", method="markers", markers_text=markers_text, resolution=resolution,
                           carry_from=previous.get(sid), knn=knn)
        s._save_mres("init", _slice_mres(res, off[i], off[i + 1] - off[i]))
        model.save(m, s.dir / "base" / "init_model.npz")
        info = json.load(open(s.dir / "session.json"))
        info["scheme"]["joint"] = {"project": pid, "run": run, "slices": len(ids)}
        json.dump(info, open(s.dir / "session.json", "w"), indent=1)
        s = Session(s.id)
        s.result()
        out[sid] = s.id
    json.dump({"run": run, "sessions": out, "n_cells": int(counts.shape[0]), "time": time.strftime("%Y-%m-%dT%H:%M:%S")},
              open(_dir(pid) / "annotation.json", "w"), indent=1)
    say(1.0, "done")
    return out


def embed(pid, n_pcs=30, n_neighbors=15, min_dist=0.3, resolution=1.0, harmony=True, progress=None):
    """PCA (Harmony across slices when asked), neighbors, Leiden and UMAP of all slices' cells at once."""
    say = progress or (lambda f, m: None)
    ids, off, counts, panel = _pooled(pid)
    t0 = time.time()
    say(0.05, f"PCA of {counts.shape[0]:,} cells of {len(ids)} slices")
    Z = markers.pca(model.lognorm(counts), n_pcs)[:, ::-1].copy()
    if harmony and len(ids) > 1:
        say(0.15, "Harmony: removing slice-to-slice differences from the PCA")
        import harmonypy
        batch = pd.DataFrame({"slice": np.repeat(ids, np.diff(off))})
        ho = harmonypy.run_harmony(Z, batch, ["slice"], verbose=False, random_state=0)
        Z = np.ascontiguousarray(np.asarray(ho.Z_corr).T if np.asarray(ho.Z_corr).shape[0] == Z.shape[1] else ho.Z_corr,
                                 dtype=np.float32)
    say(0.3, f"{n_neighbors} nearest neighbors")
    dist, nbr = markers.knn(Z, n_neighbors)
    say(0.45, "Leiden clusters")
    lab = markers._leiden(markers.knn_graph(dist, nbr), nbr, resolution, 0)
    say(0.55, "UMAP")
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import umap
    n = len(Z)
    init = Z[:, :2].astype(np.float64)
    init = 10 * (init - init.mean(0)) / max(float(np.abs(init - init.mean(0)).max()), 1e-9)
    Y = umap.UMAP(n_neighbors=n_neighbors + 1, min_dist=float(min_dist), init=init, n_jobs=-1, low_memory=True,
                  precomputed_knn=(np.c_[np.arange(n), nbr].astype(np.int64),
                                   np.c_[np.zeros(n, np.float32), dist].astype(np.float32), None)).fit_transform(Z)
    run = time.strftime("%Y%m%d-%H%M%S")
    meta = {"run_id": run, "n_pcs": int(n_pcs), "n_neighbors": int(n_neighbors), "min_dist": float(min_dist),
            "resolution": float(resolution), "n_clusters": int(lab.max()) + 1, "seconds": round(time.time() - t0, 1),
            "project": pid, "slices": ids, "harmony": bool(harmony and len(ids) > 1)}
    for i, sid in enumerate(ids):
        d = config.SAMPLES / sid / "embedding"
        d.mkdir(exist_ok=True)
        np.save(d / "umap.npy", Y[off[i]:off[i + 1]].astype(np.float32))
        np.save(d / "pca.npy", Z[off[i]:off[i + 1], :2].astype(np.float32))
        np.save(d / "leiden.npy", lab[off[i]:off[i + 1]].astype(np.int32))
        json.dump(meta, open(d / "meta.json", "w"), indent=1)
    json.dump(meta, open(_dir(pid) / "embedding.json", "w"), indent=1)
    say(1.0, "done")
    return meta


def niches(pid, k=15, n_niches=8, identity="nmf", n_programs=20, progress=None, sessions=None):
    """Niches over all slices: identities averaged over each cell's neighbors in its own slice, one k-means. They
    go to the annotations in sessions ({sample_id: session id}, those open in the app), else the latest ones."""
    from scipy.spatial import cKDTree
    from sklearn.cluster import KMeans
    from .session import Session
    say = progress or (lambda f, m: None)
    sess = latest_sessions(pid, sessions)
    ids, off, counts, panel = _pooled(pid)
    missing = [s for s in ids if s not in sess]
    if missing:
        raise ValueError(f"annotate every slice first (missing: {', '.join(missing[:4])})")
    S = [Session(sess[s]) for s in ids]
    if identity == "nmf":
        from sklearn.decomposition import NMF
        say(0.05, f"NMF: {n_programs} programs fitted on the pooled cells")
        X = model.lognorm(counts).astype(np.float32).tocsr()
        fit = np.sort(np.random.default_rng(0).choice(X.shape[0], min(60000, X.shape[0]), replace=False))
        nmf = NMF(n_components=n_programs, init="nndsvda", max_iter=300, random_state=0, tol=1e-4).fit(X[fit])
        W = np.vstack([nmf.transform(X[i:i + 50000]) for i in range(0, X.shape[0], 50000)]).astype(np.float32)
        ident = W / np.maximum(W.sum(1, keepdims=True), 1e-9)
        names, H = [f"P{i + 1}" for i in range(n_programs)], nmf.components_.astype(np.float32)
    elif identity == "cell types":
        types = S[0].types
        if any(s.types != types for s in S[1:]):
            raise ValueError("the slices have different cell types: annotate the project jointly first")
        pred = np.concatenate([s.result()["pred"] for s in S])
        ident = np.zeros((len(pred), len(types)), np.float32)
        ident[np.arange(len(pred)), pred] = 1
        names, H = list(types), None
    else:
        raise ValueError("identity must be 'nmf' or 'cell types'")
    say(0.5, "averaging over neighbors within each slice")
    avg = np.zeros_like(ident)
    for i, s in enumerate(S):
        xy = s.cells[["x", "y"]].to_numpy()
        nn = cKDTree(xy).query(xy, k=int(k) + 1)[1]
        part = ident[off[i]:off[i + 1]]
        avg[off[i]:off[i + 1]] = np.vstack([part[nn[j:j + 20000]].mean(1) for j in range(0, len(xy), 20000)])
    say(0.8, f"clustering all slices into {n_niches} niches")
    km = KMeans(n_clusters=int(n_niches), n_init=4, random_state=0).fit(avg)
    order = np.argsort(-np.bincount(km.labels_, minlength=n_niches))
    lab = np.argsort(order)[km.labels_].astype(np.int16)
    run = time.strftime("%Y%m%d-%H%M%S")
    for i, s in enumerate(S):
        d = s.dir / "niches"
        d.mkdir(exist_ok=True)
        np.save(d / "labels.npy", lab[off[i]:off[i + 1]])
        np.save(d / "profile.npy", km.cluster_centers_[order].astype(np.float32))
        if H is not None:
            np.savez(d / f"nmf_{n_programs}.npz", W=W[off[i]:off[i + 1]], H=H)
        json.dump({"run_id": run, "k": int(k), "k_by_region": {}, "radius_um": None, "identity": identity,
                   "n_programs": int(n_programs) if identity == "nmf" else None, "n_niches": int(n_niches),
                   "identity_names": names, "version": s.state()["version"], "k_range": [int(k), int(k)],
                   "seconds": 0, "names": {}, "project": pid}, open(d / "run.json", "w"), indent=1)
    say(1.0, "done")
    return {"run": run, "sessions": sess}
