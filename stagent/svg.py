"""Spatially variable genes: genes whose expression is more alike in neighboring cells than chance, by Moran's I.

For each slice the cells' K nearest neighbors in space (row-standardized weights) give Moran's I of every gene's
log-normalized expression: I = (z . W z) / (z . z) with z the gene's centered values; 0 = no spatial pattern, near 1 =
smooth patches. Its z-score is the normality approximation (Cliff & Ord): E[I] = -1/(n-1), Var from S0, S1, S2 of the
weights; with tens of thousands of cells almost any pattern is significant, so genes are ranked by I.

Across all cells the strongest patterns are those of the cell types (epithelial genes in the epithelium...).
within (e.g. one cell type) looks inside a subset; beyond_types takes each gene's mean in every cell type out first,
so what is left is spatial variation that the cell type composition does not explain (e.g. a gradient along the
crypt axis within the epithelium, or an inflamed patch). A project's slices: each slice on its own, then the mean I
(weighted by cells) and in how many slices the gene is among the top.
"""
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy import stats
from scipy.spatial import cKDTree

K = 10              # spatial neighbors per cell
MIN_CELLS = 50      # a slice (or subset) with fewer cells is skipped
CHUNK = 64          # genes at once (dense cells x CHUNK)


def weights(xy, k=K):
    """Row-standardized K-nearest-neighbor weights (CSR, n x n)."""
    n = len(xy)
    k = min(k, n - 1)
    nn = cKDTree(xy).query(xy, k=k + 1)[1][:, 1:]
    W = sp.csr_matrix((np.full(n * k, 1.0 / k, np.float32), nn.ravel(), np.arange(0, n * k + 1, k)), shape=(n, n))
    return W


def _moments(W):
    """E[I], Var[I] under normality for weights W (S0, S1, S2 of Cliff & Ord)."""
    n = W.shape[0]
    S0 = float(W.sum())
    Wt = (W + W.T).tocsr()
    S1 = 0.5 * float(Wt.multiply(Wt).sum())
    S2 = float(((np.asarray(W.sum(1)).ravel() + np.asarray(W.sum(0)).ravel()) ** 2).sum())
    EI = -1.0 / (n - 1)
    var = (n * n * S1 - n * S2 + 3 * S0 * S0) / ((n * n - 1) * S0 * S0) - EI * EI
    return EI, max(var, 1e-12)


def morans_i(X, W, codes=None):
    """Moran's I of every column of X (cells x genes, sparse or dense) for weights W (row sums 1); codes (a group per
    cell): each group's mean taken out of every gene first. CHUNK genes at a time, so memory stays small."""
    n, G = X.shape
    out = np.zeros(G)
    groups = [np.flatnonzero(codes == c) for c in np.unique(codes)] if codes is not None else None
    for a in range(0, G, CHUNK):
        Z = X[:, a:a + CHUNK]
        Z = (Z.toarray() if sp.issparse(Z) else np.asarray(Z)).astype(np.float32)
        if groups is not None:
            for m in groups:
                Z[m] -= Z[m].mean(0)
        Z -= Z.mean(0)
        num = (Z * (W @ Z)).sum(0)
        den = (Z * Z).sum(0)
        out[a:a + CHUNK] = np.where(den > 0, num / np.maximum(den, 1e-12), 0.0) * n / float(W.sum())
    return out


def _bh(p):
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def one(X, xy, types=None, k=K):
    """{I, z, p, detected} per gene of one slice (rows of X at positions xy); types: each gene's mean per type taken
    out first (beyond the cell types)."""
    codes = pd.factorize(pd.Series(types).astype(str))[0] if types is not None else None
    W = weights(xy, k)
    I = morans_i(X, W, codes)
    EI, var = _moments(W)
    z = (I - EI) / np.sqrt(var)
    return {"I": I, "z": z, "p": stats.norm.sf(z), "n": X.shape[0]}


def spatial_genes(c, within=None, beyond_types=False, top=30, by_type="cell type", k=K, progress=None):
    """Spatially variable genes of the cells of `c` (analysis.Cells: one annotation or a project's slices), among the
    cells in `within`; beyond_types: after each gene's mean per value of `by_type` is taken out."""
    from . import analysis
    say = progress or (lambda f, m: None)
    w = analysis._within(c, within)
    genes = np.asarray(c.genes)
    X = c.X
    off = np.cumsum([0] + [len(s.cells) for s in c.sessions])
    per = []
    for i, s in enumerate(c.sessions):
        rows = np.arange(off[i], off[i + 1])
        rows = rows[w[rows]]
        if len(rows) < MIN_CELLS:
            continue
        say(0.1 + 0.8 * i / len(c.sessions), f"Moran's I of {len(genes)} genes in {len(rows):,} cells"
            + (f" of {s.sample_id}" if len(c.sessions) > 1 else ""))
        xy = s.cells[["x", "y"]].to_numpy(float)[rows - off[i]]
        t = analysis._column(c, by_type)[rows] if beyond_types else None
        r = one(X[rows], xy, t, k)
        det = np.asarray((X[rows] > 0).mean(0)).ravel()
        per.append({"slice": str(c.labels.slice.to_numpy()[rows[0]]) if "slice" in c.labels else s.sample_id,
                    "n": len(rows), "I": r["I"], "z": r["z"], "det": det})
    if not per:
        raise ValueError(f"fewer than {MIN_CELLS} cells" + (" in that subset" if within else ""))
    n = np.array([x["n"] for x in per], float)
    I = (np.stack([x["I"] for x in per]) * n[:, None]).sum(0) / n.sum()
    det = (np.stack([x["det"] for x in per]) * n[:, None]).sum(0) / n.sum()
    zc = np.stack([x["z"] for x in per]).sum(0) / np.sqrt(len(per))  # Stouffer over slices
    p = stats.norm.sf(zc)
    fdr = _bh(p)
    order = np.argsort(-I)
    topk = max(1, min(int(top), len(genes)))
    rank_in = [set(np.argsort(-x["I"])[:topk].tolist()) for x in per]
    table = pd.DataFrame({"gene": genes, "morans_i": I, "z": zc, "fdr": fdr, "pct_detected": det,
                          "n_slices_top": [sum(g in r for r in rank_in) for g in range(len(genes))]}).iloc[order]
    if len(per) > 1:
        for x in per:
            table[f"I_{x['slice']}"] = x["I"][order]
    rows_out = [{"gene": r.gene, "morans_i": round(float(r.morans_i), 3), "z": round(float(r.z), 1),
                 "fdr": float(f"{r.fdr:.2e}"), "pct_detected": round(float(r.pct_detected), 3),
                 **({"n_slices_top": int(r.n_slices_top)} if len(per) > 1 else {})}
                for r in table.head(topk).itertuples()]
    return {"scope": "project" if c.pid else "slice", "within": within or {}, "beyond_types": bool(beyond_types),
            "by_type": by_type if beyond_types else None, "k": k, "n_cells": int(n.sum()),
            "slices": [{"slice": x["slice"], "n_cells": x["n"]} for x in per], "genes": rows_out,
            "n_significant": int((fdr < 0.05).sum()), "table": table,
            "note": ("Moran's I of each gene's log-normalized expression over the " + str(k) + " nearest cells: 0 = no "
                     "spatial pattern, near 1 = smooth patches. Genes ranked by I; with many cells nearly every "
                     "pattern is significant. " + ("Each gene's mean per " + by_type + " was taken out first: the "
                     "patterns left are not explained by where the cell types are. " if beyond_types else
                     "Across all cell types the patterns are mostly where the types are: use Only within (one type) "
                     "or Beyond cell types for patterns inside them. ") +
                     ("Several slices: I is the cells-weighted mean, z combined (Stouffer), n_slices_top = in how "
                      "many slices the gene is in the top list." if len(per) > 1 else ""))}
