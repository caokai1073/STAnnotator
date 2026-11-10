"""Niches: places where cells of similar make-up live together.

1. Identity of every cell: NMF of the log-normalized expression into n_programs gene programs; a cell's identity is
   its program weights scaled to sum 1 (a soft composition of programs, independent of the annotation). Or
   identity="cell types": the cell type labels of the annotation (one-hot).
2. Neighborhood: the identities of the cell and its k nearest cells are averaged. k is one number, or one per drawn
   region (k_by_region; cells in no region use k), or set by density (radius_um: every cell within that distance,
   between K_MIN and K_MAX cells), so a neighborhood can cover a similar stretch of tissue in dense and sparse areas.
3. Niches: k-means of the averaged identities into n_niches clusters.
The latest run of a session is stored in niches/ (labels, parameters, composition per niche, names Claude gave them).
"""
import json, time
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from . import jobs, model
from .osutil import dump_json

K_MIN, K_MAX = 5, 80
NMF_FIT_CELLS = 60000   # NMF is fitted on a random subset, then every cell is projected onto the programs


def _dir(s):
    d = s.dir / "niches"
    d.mkdir(exist_ok=True)
    return d


def programs(s, n_programs=20, progress=None):
    """NMF of the sample's expression: cell weights (cells x n_programs) and program loadings (n_programs x genes);
    cached in the session."""
    say = progress or (lambda f, m: None)
    f = _dir(s) / f"nmf_{n_programs}.npz"
    if f.exists():
        z = np.load(f)
        return z["W"], z["H"]
    from sklearn.decomposition import NMF
    X = model.lognorm(s.counts).astype(np.float32).tocsr()
    rng = np.random.default_rng(0)
    fit = np.sort(rng.choice(X.shape[0], min(NMF_FIT_CELLS, X.shape[0]), replace=False))
    say(0.1, f"NMF: fitting {n_programs} programs on {len(fit):,} cells")
    nmf = NMF(n_components=n_programs, init="nndsvda", max_iter=300, random_state=0, tol=1e-4)
    nmf.fit(X[fit])
    say(0.4, "NMF: projecting all cells")
    W = np.vstack([nmf.transform(X[i:i + 50000]) for i in range(0, X.shape[0], 50000)]).astype(np.float32)
    H = nmf.components_.astype(np.float32)
    np.savez(f, W=W, H=H)
    return W, H


def _neighbors_k(s, k, k_by_region=None, radius_um=None, version=None):
    """Neighborhood size of every cell."""
    n = len(s.cells)
    if radius_um:
        xy = s.cells[["x", "y"]].to_numpy()
        cnt = cKDTree(xy).query_ball_point(xy, r=float(radius_um), return_length=True)
        return np.clip(np.asarray(cnt) - 1, K_MIN, K_MAX)
    kk = np.full(n, int(k))
    if k_by_region:
        reg = s.result(version)["region"]
        names = set(s.region_names(version))
        for name, kv in k_by_region.items():
            if name not in names:
                raise ValueError(f"no region {name!r}; regions: {sorted(names)}")
            kk[reg == name] = int(kv)
    if kk.min() < 1 or kk.max() > K_MAX:
        raise ValueError(f"k must be between 1 and {K_MAX}")
    return kk


def run(s, k=15, n_niches=8, identity="nmf", n_programs=20, k_by_region=None, radius_um=None, progress=None):
    """Compute niches and store them as the session's current niche result."""
    say = progress or (lambda f, m: None)
    t0 = time.time()
    if identity == "nmf":
        W, H = programs(s, n_programs, progress)
        ident = W / np.maximum(W.sum(1, keepdims=True), 1e-9)
        names = [f"P{i + 1}" for i in range(W.shape[1])]
    elif identity == "cell types":
        pred = s.result()["pred"]
        ident = np.zeros((len(pred), len(s.types)), np.float32)
        ident[np.arange(len(pred)), pred] = 1
        names, H = list(s.types), None
    else:
        raise ValueError("identity must be 'nmf' or 'cell types'")
    say(0.6, "averaging identities over each cell's neighbors")
    kk = _neighbors_k(s, k, k_by_region, radius_um)
    xy = s.cells[["x", "y"]].to_numpy()
    nn = cKDTree(xy).query(xy, k=int(kk.max()) + 1)[1]   # column 0 is the cell itself
    avg = np.zeros_like(ident)
    for kv in np.unique(kk):
        rows = np.where(kk == kv)[0]
        for i in range(0, len(rows), 20000):
            r = rows[i:i + 20000]
            avg[r] = ident[nn[r, :kv + 1]].mean(1)
    say(0.8, f"clustering into {n_niches} niches")
    from sklearn.cluster import KMeans
    km = KMeans(n_clusters=int(n_niches), n_init=4, random_state=0).fit(avg)
    lab = km.labels_.astype(np.int16)
    order = np.argsort(-np.bincount(lab, minlength=n_niches))  # N1 = the largest niche
    lab = np.argsort(order)[lab].astype(np.int16)
    d = _dir(s)
    meta = {"run_id": time.strftime("%Y%m%d-%H%M%S"), "k": int(k), "k_by_region": k_by_region or {},
            "radius_um": radius_um, "identity": identity, "n_programs": int(n_programs) if identity == "nmf" else None,
            "n_niches": int(n_niches), "identity_names": names, "version": s.state()["version"],
            "k_range": [int(kk.min()), int(kk.max())], "seconds": round(time.time() - t0, 1), "names": {}}
    with jobs.replacing(d / "labels.npy", d / "profile.npy", d / "run.json"):  # a stopped run puts the last back
        np.save(d / "labels.npy", lab)
        np.save(d / "profile.npy", km.cluster_centers_[order].astype(np.float32))
        dump_json(meta, d / "run.json", indent=1)
    say(1.0, "done")
    return summary(s)


def current(s):
    """(meta, labels) of the session's current niche result, or (None, None)."""
    d = s.dir / "niches"
    if not (d / "run.json").exists():
        return None, None
    return json.load(open(d / "run.json")), np.load(d / "labels.npy")


def niche_names(s):
    meta, lab = current(s)
    if meta is None:
        return []
    return [meta["names"].get(f"N{i + 1}") or f"N{i + 1}" for i in range(meta["n_niches"])]


def summary(s, version=None, top=5):
    """Per niche: size, cell types (from the annotation), identity profile (programs and their genes), regions,
    distance to the epithelium."""
    meta, lab = current(s)
    if meta is None:
        return None
    r = s.result(version)
    types = np.array(s.types, dtype=object)[r["pred"]]
    prof = np.load(s.dir / "niches" / "profile.npy")
    H = None
    if meta["identity"] == "nmf":
        H = np.load(s.dir / "niches" / f"nmf_{meta['n_programs']}.npz")["H"]
    depth = s.cells.depth_um.to_numpy()
    names = niche_names(s)
    out = []
    for i in range(meta["n_niches"]):
        m = lab == i
        item = {"niche": f"N{i + 1}", "name": names[i], "n_cells": int(m.sum()), "share": round(float(m.mean()), 3),
                "cell_types": pd.Series(types[m]).value_counts(normalize=True).head(top).round(3).to_dict(),
                "regions": pd.Series(r["region"][m]).value_counts(normalize=True).head(3).round(2).to_dict(),
                "depth_um_median": None if np.isnan(depth[m]).all() else round(float(np.nanmedian(depth[m])), 1)}
        o = np.argsort(-prof[i])[:3]
        item["identity"] = [{"name": meta["identity_names"][j], "weight": round(float(prof[i, j]), 3),
                             **({"genes": [s.panel[g] for g in np.argsort(-H[j])[:6]]} if H is not None else {})}
                            for j in o]
        out.append(item)
    return {"params": {k: meta[k] for k in ("identity", "n_programs", "k", "k_by_region", "radius_um", "n_niches",
                                            "k_range", "run_id")}, "niches": out}


def set_names(s, names):
    """Names for niches ({"N1": "crypt base", ...}); shown in the viewer's legend."""
    meta, _ = current(s)
    if meta is None:
        raise ValueError("no niches yet: run them first")
    for k_ in names:
        if k_ not in {f"N{i + 1}" for i in range(meta["n_niches"])}:
            raise ValueError(f"no niche {k_}")
    meta["names"].update({k_: str(v).strip()[:40] for k_, v in names.items()})
    dump_json(meta, s.dir / "niches" / "run.json", indent=1)
    return niche_names(s)


def context(s, version=None, ks=(10, 15, 20, 30, 50)):
    """For choosing k: per drawn region (and the cells in none), cell density, the median distance covered by k
    neighbors, and the main lineages."""
    r = s.result(version)
    xy = s.cells[["x", "y"]].to_numpy()
    d = cKDTree(xy).query(xy, k=max(ks) + 1)[0]
    lin = s.lineage_labels(version)
    out = []
    for name in list(dict.fromkeys(list(s.region_names(version)) + ["unknown"])):
        m = r["region"] == name
        if m.sum() < 50:
            continue
        # density from the distance to the 10th neighbor (cells per mm2 within that disc)
        rad = np.median(d[m, 10])
        dens = 10 / (np.pi * (rad / 1000) ** 2) if rad > 0 else None
        out.append({"region": "no region" if name == "unknown" else name, "n_cells": int(m.sum()),
                    "cells_per_mm2": None if dens is None else int(dens),
                    "radius_um_covered_by_k": {str(k_): round(float(np.median(d[m, k_])), 1) for k_ in ks},
                    "lineages": pd.Series(lin[m]).value_counts(normalize=True).head(4).round(2).to_dict()})
    return {"regions": out, "k_max": K_MAX,
            "note": "radius_um_covered_by_k = median distance from a cell to its k-th neighbor; similar radii across "
                    "regions mean neighborhoods of a similar tissue scale. Regions are the user's drawn regions."}
