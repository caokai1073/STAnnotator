"""Trajectories: where cells lie along a spatial axis or an expression pseudotime, and the genes that change along it.

Methods, on the cells of an analysis (one slice, or the pooled slices of a project), optionally only within a subset
(e.g. one lineage or region):
  spatial     the cells' position on an axis between two groups of cells of any grouping (cell types, regions,
              niches, ...): t = d_A / (d_A + d_B), with d the distance (um) to the nearest cell of the group in the same
              slice; 0 in A, 1 in B. With only a start group, the distance to it (um), shown scaled to 0..1.
  expression  diffusion pseudotime (Haghverdi et al. 2016, Nat Methods 13:845): a diffusion map of the kNN graph of
              the cells' PCA (adaptive Gaussian kernel, density-normalized), and t = the diffusion distance from the
              root group (the mean of its cells in the scaled diffusion components), scaled to 0..1.
  mixed       the same on a graph that also joins cells that are close in the tissue:
              W = (1 - alpha) W_expression + alpha W_space (alpha 0: expression only; 1: space only).
Up to MAX_GRAPH_CELLS cells the graph holds every cell (about 5 s for 100,000, 25 s for 350,000); above, a random sample (the same for
every root) and the other cells take the mean t of their neighbors that have one. Cells in a part of the graph that
the root does not reach get no t. The PCA, graphs and diffusion maps of the last cell sets are kept (_CACHE), so a new
root, weight of space or composition on the same cells takes about a second. Each method keeps its own last run.

Genes along t (gene_trends): each gene's correlation with the rank of t (rising and falling genes), and its mean in
NBINS bins of equal cell numbers along t: the share of a gene's variance between bins (eta squared) finds genes that
peak in the middle (transient). The cells are treated as independent: rank by effect, not by p.

The core functions take arrays only, so the exported notebook can carry this file and run them.
"""
import numpy as np, scipy.sparse as sp

K = 15                 # neighbors of the expression and the spatial graph
N_PCS = 30
N_DCS = 10             # diffusion components used for the pseudotime
MAX_GRAPH_CELLS = 1000000  # beyond: the graph of a random sample, the other cells placed by their neighbors
EIG_TOL = 1e-6  # of the diffusion map's eigenvectors: the order is the same as at machine precision (rank r = 1.0)
NBINS = 20
MAX_GENE_CELLS = 30000
MIN_DETECTED = 0.02    # genes detected in fewer cells than this share (of the cells on the trajectory) are left out
TOP = 25
SEED = 0


# ---------------- spatial axis ----------------
def distance_to(xy, slices, target):
    """Distance (um) of every cell to the nearest cell of target (a boolean mask) in its own slice; inf when its slice
    has none."""
    from scipy.spatial import cKDTree
    d = np.full(len(xy), np.inf)
    for sl in np.unique(slices):
        m = slices == sl
        tm = m & target
        if tm.any():
            d[m] = cKDTree(xy[tm]).query(xy[m])[0]
    return d


def spatial_axis(xy, slices, start, end=None):
    """(t, raw): with end, t = d_start / (d_start + d_end) in 0..1 (raw = t); without, raw = the distance to start (um)
    and t = raw scaled by its 99th percentile (clipped at 1). NaN where a slice lacks the group(s)."""
    dA = distance_to(xy, slices, start)
    if end is None:
        raw = np.where(np.isfinite(dA), dA, np.nan)
        hi = np.nanpercentile(raw, 99) if np.isfinite(raw).any() else 1.0
        return np.clip(raw / max(hi, 1e-9), 0, 1), raw
    dB = distance_to(xy, slices, end)
    with np.errstate(invalid="ignore", divide="ignore"):
        t = dA / (dA + dB)
    t[start] = 0.0
    t[end & ~start] = 1.0
    t[~(np.isfinite(dA) & np.isfinite(dB))] = np.nan
    return t, t.copy()


# ---------------- graphs and diffusion pseudotime ----------------
def _kernel(dist, nbr, n):
    """Symmetric sparse kernel from kNN lists: w_ij = exp(-d_ij^2 / (s_i s_j)), s_i the median distance of i's
    neighbors (the adaptive Gaussian of diffusion maps)."""
    s = np.maximum(np.median(dist, 1), 1e-9)
    rows = np.repeat(np.arange(n), dist.shape[1])
    w = np.exp(-(dist.ravel() ** 2) / (s[rows] * s[nbr.ravel()]))
    W = sp.csr_matrix((w, (rows, nbr.ravel())), shape=(n, n))
    return W.maximum(W.T).tocsr()


def spatial_graph(xy, slices, k=K):
    """Kernel of each cell's k nearest cells in the tissue (um), within its slice."""
    from scipy.spatial import cKDTree
    n = len(xy)
    s, parts = np.ones(n), []
    for sl in np.unique(slices):
        idx = np.where(slices == sl)[0]
        kk = min(k, len(idx) - 1)
        if kk < 1:  # a cell alone in its slice: no edges
            continue
        d, j = cKDTree(xy[idx]).query(xy[idx], kk + 1)
        d, j = d[:, 1:], idx[j[:, 1:]]
        s[idx] = np.maximum(np.median(d, 1), 1e-9)
        parts.append((np.repeat(idx, kk), j.ravel(), d.ravel()))
    if not parts:
        return sp.csr_matrix((n, n))
    r, c, dd = (np.concatenate(x) for x in zip(*parts))
    W = sp.csr_matrix((np.exp(-(dd ** 2) / (s[r] * s[c])), (r, c)), shape=(n, n))
    return W.maximum(W.T).tocsr()


def expression_graph(Z, k=K):
    from . import markers
    dist, nbr = markers.knn(Z, k, mode="auto")
    return _kernel(np.asarray(dist, float), np.asarray(nbr), len(Z))


def diffusion_components(W, n_dcs=N_DCS):
    """Scaled diffusion components (cells x n_dcs; lambda / (1 - lambda) times the eigenvectors after the first) of
    the density-normalized kernel W, and the eigenvalues."""
    from scipy.sparse.linalg import eigsh
    W = W.tocsr().astype(np.float64)
    q = np.asarray(W.sum(1)).ravel()
    q[q == 0] = 1
    Wn = sp.diags(1 / q) @ W @ sp.diags(1 / q)  # density normalization (alpha = 1)
    z = np.asarray(Wn.sum(1)).ravel()
    z[z == 0] = 1
    T = sp.diags(1 / np.sqrt(z)) @ Wn @ sp.diags(1 / np.sqrt(z))
    k = min(n_dcs + 1, W.shape[0] - 2)
    lam, U = eigsh(T, k=k, which="LA", v0=np.random.default_rng(SEED).random(W.shape[0]), tol=EIG_TOL)
    o = np.argsort(-lam)
    lam, U = lam[o], U[:, o]
    psi = U / np.sqrt(z)[:, None]  # right eigenvectors of the transition matrix
    lam = np.clip(lam[1:], 0, 1 - 1e-6)
    return psi[:, 1:] * (lam / (1 - lam))[None, :], lam


def pseudotime(Z, xy, slices, root, alpha=0.0, k=K, n_dcs=N_DCS, max_cells=MAX_GRAPH_CELLS, progress=None, cache=None):
    """Diffusion pseudotime from the root cells (mask), on the expression graph of Z (PCA) mixed with the spatial
    graph by alpha. Returns (t in 0..1 with NaN where the root does not reach, info). cache: a dict kept between runs
    on the same cells (Z, xy, slices): the graphs and diffusion maps in it are reused, so a new root or composition
    costs a second and a new alpha only its own graph mix."""
    say = progress or (lambda f, m: None)
    from scipy.sparse.csgraph import connected_components
    c = cache if cache is not None else {}
    n = len(Z)
    if not root.any():
        raise ValueError("the root group has no cells among the cells of the trajectory")

    def get(key, fn):
        if key not in c:
            c[key] = fn()
        return c[key]
    if n > max_cells:  # a random sample for the graph (the same for every root), the other cells placed after
        samp = get("samp", lambda: np.sort(np.random.default_rng(SEED).choice(n, max_cells, replace=False)))
        tag = "s"
        if not root[samp].any():  # a rare root: its cells join the sample (graphs of this run only)
            samp, tag, c = np.union1d(samp, np.where(root)[0]), "r", {}
    else:
        samp, tag = np.arange(n), "a"
    rs = root[samp]
    say(0.2, f"neighbor graph of {len(samp):,} cells" + (" (kept from the last run)" if (tag, "We") in c or (tag, "Ws") in c else ""))
    W = sp.csr_matrix((len(samp), len(samp)))
    if alpha < 1:
        W = W + (1 - alpha) * get((tag, "We"), lambda: expression_graph(Z[samp], k))
    if alpha > 0:
        W = W + alpha * get((tag, "Ws"), lambda: spatial_graph(xy[samp], slices[samp], k))
    ncomp, comp = get((tag, "cc", alpha), lambda: connected_components(W, directed=False))
    comps = tuple(int(x) for x in np.unique(comp[rs]))
    sub = np.where(np.isin(comp, comps))[0]
    say(0.45, "diffusion map")
    D, lam = get((tag, "dm", alpha, comps), lambda: diffusion_components(W[sub][:, sub], n_dcs))
    r0 = D[rs[sub]].mean(0)
    ts = np.full(len(samp), np.nan)
    ts[sub] = np.sqrt(((D - r0) ** 2).sum(1))
    hi = np.nanpercentile(ts, 99.5)
    ts = np.clip(ts / max(hi, 1e-12), 0, 1)
    t = np.full(n, np.nan)
    t[samp] = ts
    if len(samp) < n:  # the other cells: the mean t of their neighbors that have one (expression, and space when mixed)
        say(0.75, "placing the other cells")
        lists = []
        if alpha < 1:
            from . import markers
            lists.append((1 - alpha, get("nbr_all", lambda: np.asarray(markers.knn(Z, k, mode="auto")[1]))))
        if alpha > 0:
            lists.append((alpha, get("snbr_all", lambda: _spatial_neighbors(xy, slices, k))))
        rest = np.setdiff1d(np.arange(n), samp)
        for _ in range(8):  # cells whose neighbors are all unplaced take a value in a later round
            todo = rest[np.isnan(t[rest])]
            if not len(todo):
                break
            num, den = np.zeros(len(todo)), np.zeros(len(todo))
            for w, nb in lists:
                v = t[nb[todo]]
                f = np.isfinite(v)
                num += w * np.where(f, v, 0).sum(1) / np.maximum(f.sum(1), 1)
                den += w * (f.sum(1) > 0)
            new = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
            if not np.isfinite(new).any():
                break
            t[todo] = new
    info = {"n_graph_cells": int(len(samp)), "n_components": int(ncomp), "n_unreached": int(np.isnan(t).sum()),
            "eigenvalues": [round(float(x), 4) for x in lam]}
    return t, info


def _spatial_neighbors(xy, slices, k=K):
    """Each cell's k nearest cells in the tissue within its slice (itself where a slice has too few)."""
    from scipy.spatial import cKDTree
    out = np.tile(np.arange(len(xy))[:, None], (1, k))
    for sl in np.unique(slices):
        idx = np.where(slices == sl)[0]
        kk = min(k, len(idx) - 1)
        if kk >= 1:
            out[idx, :kk] = idx[cKDTree(xy[idx]).query(xy[idx], kk + 1)[1][:, 1:]]
    return out


_CACHE = {}  # graphs of the last few cell sets ordered (trajectory.run), most recent last
CACHE_SETS = 2


def _cache_for(key):
    if key in _CACHE:
        _CACHE[key] = _CACHE.pop(key)
    else:
        _CACHE[key] = {}
        while len(_CACHE) > CACHE_SETS:
            _CACHE.pop(next(iter(_CACHE)))
    return _CACHE[key]


# ---------------- genes and composition along t ----------------
def _bh(p):
    """Benjamini-Hochberg adjusted p-values."""
    p = np.asarray(p, float)
    n = len(p)
    if not n:
        return p
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.clip(q, 0, 1)
    return out



def _bins(t, nbins):
    """Bin of every value of t (equal numbers of cells), and the bins' centers (median t)."""
    edges = np.unique(np.quantile(t, np.linspace(0, 1, nbins + 1)))
    b = np.clip(np.searchsorted(edges, t, side="right") - 1, 0, len(edges) - 2)
    centers = np.array([np.median(t[b == i]) if (b == i).any() else np.nan for i in range(len(edges) - 1)])
    return b, centers


def gene_trends(X, t, genes, nbins=NBINS, top=TOP, max_cells=MAX_GENE_CELLS, min_detected=MIN_DETECTED):
    """Genes along t (cells with a finite t; X log-normalized, cells x genes, sparse): a table of every gene detected
    in min_detected of them (gene, r = correlation with the rank of t, p, fdr, eta2 = share of variance between bins,
    peak = t of its highest bin, mean, pct), the top rising / falling / transient genes, and the bin means of these
    genes for a heatmap (z-scored, genes ordered by their peak)."""
    import pandas as pd
    from scipy.stats import t as tdist, rankdata
    ok = np.where(np.isfinite(t))[0]
    if len(ok) < 30:
        raise ValueError("fewer than 30 cells on the trajectory")
    rng = np.random.default_rng(SEED)
    rows = np.sort(rng.choice(ok, max_cells, replace=False)) if len(ok) > max_cells else ok
    Xs = sp.csr_matrix(X[rows], dtype=np.float64)
    ts = t[rows]
    n = len(rows)
    pct = np.asarray((Xs > 0).mean(0)).ravel()
    mu = np.asarray(Xs.mean(0)).ravel()
    sd = np.sqrt(np.maximum(np.asarray(Xs.multiply(Xs).mean(0)).ravel() - mu ** 2, 0))
    r_ = rankdata(ts)
    rc = (r_ - r_.mean()) / (r_.std() + 1e-12)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (Xs.T @ rc) / n / sd
    r = np.where(sd > 0, r, 0.0)
    tt = r * np.sqrt(max(n - 2, 1) / np.maximum(1 - r ** 2, 1e-12))
    p = 2 * tdist.sf(np.abs(tt), max(n - 2, 1))
    b, centers = _bins(ts, nbins)
    nb = len(centers)
    B = sp.csr_matrix((np.ones(n), (b, np.arange(n))), shape=(nb, n))
    cnt = np.asarray(B.sum(1)).ravel()
    means = np.asarray((B @ Xs).todense()) / np.maximum(cnt, 1)[:, None]  # bins x genes
    between = (cnt[:, None] * (means - mu[None, :]) ** 2).sum(0) / n
    with np.errstate(invalid="ignore", divide="ignore"):
        eta2 = np.where(sd > 0, between / sd ** 2, 0.0)
    peak = centers[np.argmax(means, 0)]
    hi, lo = means.max(0), means.min(0)
    with np.errstate(invalid="ignore", divide="ignore"):  # how far it falls again on both sides of its peak (0..1)
        hump = np.where(hi > lo, (hi - np.maximum(means[0], means[-1])) / (hi - lo), 0.0)
    keep = pct >= min_detected
    tab = pd.DataFrame({"gene": genes, "r": r, "p": p, "eta2": eta2, "peak": peak, "hump": hump, "mean": mu,
                        "pct": pct})[keep]
    tab["fdr"] = _bh(tab.p.to_numpy())
    tab = tab.sort_values("r", ascending=False).reset_index(drop=True)
    rising = tab[tab.r > 0].head(top)
    falling = tab[tab.r < 0].sort_values("r").head(top)
    mid = (tab.peak > np.nanquantile(ts, 0.2)) & (tab.peak < np.nanquantile(ts, 0.8))
    trans = tab[mid & (tab.hump > 0.4) & (tab.eta2 > 2 * tab.r ** 2)].sort_values("eta2", ascending=False).head(top)
    shown = list(dict.fromkeys(list(rising.gene[:12]) + list(trans.gene[:8]) + list(falling.gene[:12])))
    gi = {g: i for i, g in enumerate(genes)}
    H = means[:, [gi[g] for g in shown]]
    Hz = (H - H.mean(0)) / (H.std(0) + 1e-9)
    order = np.argsort(centers[np.argmax(H, 0)], kind="stable")
    rec = lambda d: [{k: (round(float(v), 4) if isinstance(v, (float, np.floating)) else v) for k, v in x.items()}
                     for x in d[["gene", "r", "eta2", "peak", "pct", "fdr"]].to_dict("records")]
    return {"n_cells": int(len(ok)), "n_tested_cells": int(n), "bins": [round(float(c), 4) for c in centers],
            "bin_cells": cnt.astype(int).tolist(), "rising": rec(rising), "falling": rec(falling), "transient": rec(trans),
            "heatmap": {"genes": [shown[i] for i in order], "z": np.round(Hz[:, order].T, 3).tolist(),
                        "mean": np.round(H[:, order].T, 4).tolist()}}, tab


def composition(t, labels, nbins=NBINS, top=12):
    """Shares of the groups (labels) in each bin of t (equal numbers of cells): the largest `top` groups, the rest
    as "other"."""
    import pandas as pd
    ok = np.isfinite(t) & pd.notna(labels)
    tt, ll = t[ok], np.asarray(labels, dtype=object)[ok]
    if not len(tt):
        return {"groups": [], "shares": []}
    b, centers = _bins(tt, nbins)
    vc = pd.Series(ll).value_counts()
    keepg = list(vc.index[:top])
    ll = np.where(np.isin(ll, keepg), ll, "other")
    groups = keepg + (["other"] if (ll == "other").any() else [])
    ct = pd.crosstab(b, ll).reindex(columns=groups, fill_value=0)
    sh = ct.div(ct.sum(1), axis=0)
    return {"groups": [str(g) for g in groups], "bins": [round(float(c), 4) for c in centers],
            "shares": np.round(sh.to_numpy(), 4).tolist()}


# ---------------- runs on an analysis' cells (analysis.Cells) ----------------
def _xy(c):
    return np.vstack([s.cells[["x", "y"]].to_numpy(float) for s in c.sessions])


def run(c, method, within=None, start=None, end=None, root=None, alpha=0.3, by="cell type", progress=None):
    """A trajectory over the cells c (analysis.cells), only those inside within ({grouping: values}).
    spatial: start (and optionally end) = {"by": grouping, "values": [...]} (the groups searched among all cells);
    expression / mixed: root = {"by", "values"} among the cells of the trajectory; alpha: the weight of space (mixed).
    by: the grouping whose composition along t is shown. Returns (result, table of genes, t and raw of every cell,
    the masks {"within", "start", "end", "root"} that define the run, for its replay in the notebook)."""
    from . import analysis, markers, model
    say = progress or (lambda f, m: None)
    if method not in ("spatial", "expression", "mixed"):
        raise ValueError("method must be spatial, expression or mixed")
    sel = analysis._within(c, within)
    if sel.sum() < 30:
        raise ValueError(f"only {int(sel.sum())} cells within {within}: too few for a trajectory")
    xy, slices = _xy(c), np.asarray(c.slices)
    grp = lambda g: analysis._mask(c, g["by"], g["values"])
    note = []
    if method == "spatial":
        if not start or not start.get("values"):
            raise ValueError("choose where the axis starts")
        A = grp(start)
        Bm = grp(end) if end and end.get("values") else None
        say(0.2, "distances in the tissue")
        tt, raw = spatial_axis(xy, slices, A, Bm)
        masks = {"within": sel, "start": A, "end": Bm}
        t = np.where(sel, tt, np.nan)
        raw = np.where(sel, raw, np.nan)
        unit = "0..1" if Bm is not None else "um"
        info = {"axis": {"start": start, "end": end if Bm is not None else None, "unit": unit}}
        if Bm is None:
            info["axis"]["max_um"] = round(float(np.nanpercentile(raw, 99)), 1) if np.isfinite(raw).any() else None
        miss = sorted({str(x) for x in np.unique(slices[sel & ~np.isfinite(tt)])})
        if miss:
            note.append(f"no start or end cells in {', '.join(miss)}: its cells have no position")
    else:
        if not root or not root.get("values"):
            raise ValueError("choose the root (where the trajectory starts)")
        R = grp(root) & sel
        masks = {"within": sel, "root": R}
        if not R.any():
            raise ValueError(f"no {root['by']} {', '.join(map(str, root['values']))} among the cells of the trajectory")
        rows = np.where(sel)[0]
        a = float(alpha) if method == "mixed" else 0.0
        if not 0 <= a <= 1:
            raise ValueError("the weight of space must be between 0 and 1")
        import hashlib
        cache = _cache_for((tuple((s.id, getattr(s, "_cells_mtime", None)) for s in c.sessions),
                            hashlib.sha1(rows.tobytes()).hexdigest()))
        if "Z" not in cache:
            say(0.05, f"PCA of {len(rows):,} cells")
            cache["Z"] = markers.pca(c.X[rows], N_PCS)[:, ::-1].copy()
        ts, info = pseudotime(cache["Z"], xy[rows], slices[rows], R[rows], alpha=a, progress=say, cache=cache)
        t = np.full(len(sel), np.nan)
        t[rows] = ts
        raw = t.copy()
        info |= {"root": root, "alpha": a}
        if c.pid and len(set(slices[rows])) > 1:
            note.append("slices pooled without batch correction: differences between slices can shape the expression graph")
        if info["n_unreached"]:
            note.append(f"{info['n_unreached']:,} cells are in parts of the graph the root does not reach: they have no pseudotime")
    say(0.85, "genes along the trajectory")
    genes, tab = gene_trends(c.X, t, list(c.genes))
    comp = composition(t, c.labels[by].to_numpy() if by in c.labels.columns else np.full(len(t), None, dtype=object))
    res = {"method": method, "within": within or {}, "by": by, **info, **genes, "composition": comp,
           "n_cells_selected": int(sel.sum()), "n_on_trajectory": int(np.isfinite(t).sum()),
           "note": " ".join(note + ["Correlations treat cells as independent: rank genes by r and eta2, not by p."])}
    return res, tab, t, raw, masks


def label(res):
    """A short title of a run."""
    g = lambda x: f"{x['by']} {', '.join(map(str, x['values']))}" if x else ""
    if res["method"] == "spatial":
        a = res["axis"]
        return f"from {g(a['start'])} to {g(a['end'])}" if a.get("end") else f"distance to {g(a['start'])}"
    return f"pseudotime from {g(res['root'])}" + (f" (space {res['alpha']:g})" if res["method"] == "mixed" else "")


# ---------------- the last run of an annotation, for its colors on the map ----------------
_VIRIDIS = [(68, 1, 84), (72, 40, 120), (62, 74, 137), (49, 104, 142), (38, 130, 142), (31, 158, 137),
            (53, 183, 121), (109, 205, 89), (180, 222, 44), (253, 231, 37)]
TRAJ_BINS = [(f"t{i}", c) for i, c in enumerate(_VIRIDIS)]
NONE_KEY, NONE_RGB = "none", (200, 200, 200)


METHODS = ("spatial", "expression", "mixed")
METHOD_NAMES = {"spatial": "spatial axis", "expression": "pseudotime", "mixed": "expression + space"}


def save(s, t, raw, res, run_id, masks=None, request=None):
    """Keep a run's t of the cells of annotation s (its rows of the analysis) for the map (color "trajectory:<method>"),
    with the masks and the request that define it (masks: this annotation's rows; the notebook replays the run).
    Each method keeps its own last run (trajectory/<method>.npz and .json)."""
    import json, time
    d = s.dir / "trajectory"
    d.mkdir(exist_ok=True)
    m = {f"mask_{k}": np.asarray(v, bool) for k, v in (masks or {}).items() if v is not None}
    meth = res["method"]
    np.savez_compressed(d / f"{meth}.npz", t=np.asarray(t, np.float32), raw=np.asarray(raw, np.float32), **m)
    meta = {"run_id": run_id, "label": label(res), "method": meth, "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "unit": (res.get("axis") or {}).get("unit", "0..1"), "max_um": (res.get("axis") or {}).get("max_um"),
            "alpha": res.get("alpha"), "within": res.get("within"), "by": res.get("by"),
            "pooled": bool(request and request.get("scope") == "project"),
            "top": {k: [x["gene"] for x in res[k][:10]] for k in ("rising", "falling", "transient")}}
    json.dump(meta, open(d / f"{meth}.json", "w"), indent=1)


def _files(s, method):
    """The files of a method's last run (an annotation of an earlier build kept one run, last.*, of any method)."""
    import json
    d = s.dir / "trajectory"
    if (d / f"{method}.json").exists():
        return d / f"{method}.json", d / f"{method}.npz"
    if (d / "last.json").exists() and json.load(open(d / "last.json")).get("method") == method:
        return d / "last.json", d / "last.npz"
    return None, None


def last(s, method=None):
    """(meta, t, raw) of the annotation's last trajectory of a method (None: the most recent of any), or
    (None, None, None)."""
    import json
    if method is None:
        runs = [x for x in (all_last(s)).values()]
        if not runs:
            return None, None, None
        method = max(runs, key=lambda x: x.get("time") or "")["method"]
    j, z = _files(s, method)
    if j is None:
        return None, None, None
    z = np.load(z)
    if len(z["t"]) != len(s.cells):  # the cells changed (filter) since
        return None, None, None
    return json.load(open(j)), z["t"], z["raw"]


def all_last(s):
    """{method: meta} of the methods that have a last run."""
    import json
    out = {}
    for m in METHODS:
        j, z = _files(s, m)
        if j is not None:
            out[m] = json.load(open(j))
    return out


def mode_method(mode):
    """The method of a color mode "trajectory:<method>" ("trajectory" alone: the most recent)."""
    return mode.split(":", 1)[1] if ":" in mode else None


def bin_of(t):
    """Index into TRAJ_BINS of every cell, -1 without a position."""
    return np.where(np.isfinite(t), np.clip((np.nan_to_num(t) * 10).astype(int), 0, 9), -1)


def legend(meta):
    """Labels of the ten bins of a run."""
    if meta.get("unit") == "um" and meta.get("max_um"):
        m = meta["max_um"]
        return [f"{i * m / 10:.0f}–{(i + 1) * m / 10:.0f} µm" + (" and more" if i == 9 else "") for i in range(10)]
    return [f"{i / 10:.1f}–{(i + 1) / 10:.1f}" for i in range(10)]
