"""Spatial co-localization: which groups of cells sit next to each other more (or less) often than by chance.

Neighborhood enrichment (as squidpy's nhood_enrichment, at a fixed distance):
1. Neighbor pairs: every two cells of a slice whose centers lie within radius_um of each other (pairs never cross
   slices). Cells without a value of the grouping (or outside `within`) take no part.
2. Observed: the number of neighbor pairs between group A and group B (A next to A: pairs inside A).
3. Expected: the same count with the labels shuffled among the cells of each slice (the positions, the group sizes
   and the slices stay), n_perm times. z = (observed - mean) / sd of the shuffles; ratio = observed / mean; p:
   two-sided, from the shuffles; FDR (Benjamini-Hochberg) over the pairs A <= B.
The shuffles are seeded: the same cells, grouping and settings give the same numbers on every computer.

A pair above 0 (z > 0) sits together more than its group sizes imply (shared niche, contact, a layer); below 0 the
two avoid each other (separate compartments). A group next to itself (the diagonal) measures how clustered it is.
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from . import fastpar

RADIUS_UM = 30.0
MAX_RADIUS_UM = 100.0
N_PERM = 200
MAX_GROUPS = 60
SEED = 0
NOTE = ("Neighbors: cells whose centers lie within the radius. Expected: the labels shuffled among the cells of each "
        "slice. z > 0: the two groups sit together more often than their sizes imply; z < 0: they avoid each other; "
        "the diagonal is how clustered a group is. Co-localization suggests, not proves, interaction: a shared layer "
        "or compartment also brings cells together.")
_K = {}


def _kernel():
    if not _K:
        import os
        from numba import njit
        # numba's disk cache needs the module's source file: the built app carries only its bytecode (as fastpar)
        cache = bool(globals().get("__file__")) and os.path.isfile(globals()["__file__"])

        @njit(nogil=True, cache=cache)
        def count(a, b, lab, K):
            out = np.zeros(K * K, np.int64)
            for e in range(len(a)):
                i, j = lab[a[e]], lab[b[e]]
                if i <= j:
                    out[i * K + j] += 1
                else:
                    out[j * K + i] += 1
            return out
        _K["count"] = count
    return _K["count"]


MAX_PAIRS = 150_000_000  # neighbor pairs of one slice (16 bytes each, and copies): beyond, a smaller radius is asked for


def _check_pairs(tree, xy, radius):
    """Refuse a radius that would make more than MAX_PAIRS pairs (estimated from up to 2,000 cells): the pairs of a
    large radius on a million cells would not fit in memory."""
    n = len(xy)
    if n < 20000:
        return
    probe = xy[np.random.default_rng(0).choice(n, 2000, replace=False)]
    per_cell = float(np.mean([len(x) for x in tree.query_ball_point(probe, radius)])) - 1
    est = n * per_cell / 2
    if est > MAX_PAIRS:
        ok = radius * (MAX_PAIRS / est) ** 0.5
        raise ValueError(f"a radius of {radius:g} um makes about {est / 1e6:.0f} million neighbor pairs among {n:,} "
                         f"cells, too many for memory: use {max(5, int(ok)):d} um or less")


def _pairs(c, keep, radius):
    """Neighbor pairs (a, b), a < b, rows of the pooled table, slice by slice; and the row ranges of the slices."""
    A, B, ranges, off = [], [], [], 0
    for s in c.sessions:
        n = len(s.cells)
        rows = np.flatnonzero(keep[off:off + n])
        if len(rows) > 1:
            xy = s.cells[["x", "y"]].to_numpy(float)[rows]
            ok = np.isfinite(xy).all(1)
            rows, xy = rows[ok], xy[ok]
            tree = cKDTree(xy)
            _check_pairs(tree, xy, radius)
            p = tree.query_pairs(radius, output_type="ndarray")
            if len(p):
                A.append((rows[p[:, 0]] + off).astype(np.int64)); B.append((rows[p[:, 1]] + off).astype(np.int64))
        ranges.append((off, off + n))
        off += n
    if not A:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), ranges
    return np.concatenate(A), np.concatenate(B), ranges


def _bh(p):
    p = np.asarray(p, float)
    if not len(p):
        return p
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


def run(c, by="cell type", within=None, radius_um=RADIUS_UM, n_perm=N_PERM, progress=None):
    """Neighborhood enrichment of the groups of `by` among the cells of the table c (analysis.cells), within
    radius_um. Returns (result, table): the result for the app and the agent, the table of every pair for CSV."""
    from . import analysis
    say = progress or (lambda f, m: None)
    radius = float(radius_um)
    if not 0 < radius <= MAX_RADIUS_UM:
        raise ValueError(f"the radius is between 0 and {MAX_RADIUS_UM:g} µm")
    n_perm = int(n_perm)
    if not 20 <= n_perm <= 2000:
        raise ValueError("the number of shuffles is between 20 and 2000")
    col = analysis._column(c, by)
    keep = pd.notna(col) & (col != None)  # noqa: E711 (None labels: not part of this grouping)
    if within:
        keep &= analysis._within(c, within)
    names = sorted({str(x) for x in col[keep]}, key=str)
    if len(names) < 2:
        raise ValueError(f"fewer than two groups of {by} among these cells")
    if len(names) > MAX_GROUPS:
        raise ValueError(f"{len(names)} groups of {by}: at most {MAX_GROUPS} (choose a coarser grouping, or Only within)")
    K = len(names)
    index = {g: i for i, g in enumerate(names)}
    lab = np.full(len(col), -1, np.int64)
    lab[keep] = [index[str(x)] for x in col[keep]]
    say(0.05, "finding the neighbor pairs")
    a, b, ranges = _pairs(c, keep, radius)
    if not len(a):
        raise ValueError(f"no two cells lie within {radius:g} µm of each other: choose a larger radius")
    count = _kernel()
    obs = count(a, b, lab, K).reshape(K, K)
    # the shuffles: labels permuted among the kept cells of each slice, each shuffle seeded by its number
    parts = [np.flatnonzero(keep[s:e]) + s for s, e in ranges]

    def shuffle(k):
        rng = np.random.default_rng([SEED, k])
        perm = lab.copy()
        for rows in parts:
            if len(rows) > 1:
                perm[rows] = lab[rows][rng.permutation(len(rows))]
        return count(a, b, perm, K)

    done = lambda k: say(0.1 + 0.85 * k / n_perm, f"shuffling the labels ({k} of {n_perm})")
    null = np.array(fastpar.pmap(shuffle, range(n_perm), done=done), np.float64).reshape(n_perm, K, K)
    mean, sd = null.mean(0), null.std(0, ddof=1)
    z = np.where(sd > 0, (obs - mean) / np.where(sd > 0, sd, 1), 0.0)
    ratio = np.where(mean > 0, obs / np.where(mean > 0, mean, 1), np.nan)
    lr = np.full_like(ratio, np.nan)
    np.log2(ratio, out=lr, where=ratio > 0)
    p = (1 + (np.abs(null - mean) >= np.abs(obs - mean)[None]).sum(0)) / (n_perm + 1)  # (at least 1 / (n_perm + 1))
    iu = np.triu_indices(K)
    fdr = np.full((K, K), np.nan)
    fdr[iu] = _bh(p[iu])
    sym = lambda m: np.where(np.arange(K)[:, None] <= np.arange(K)[None, :], m, m.T)  # the lower half from the upper
    obs, mean, sd, z, ratio, lr, p, fdr = (sym(x) for x in (obs.astype(float), mean, sd, z, ratio, lr, p, fdr))
    sizes = {g: int((lab == i).sum()) for g, i in index.items()}
    rows = [{"group_a": names[i], "group_b": names[j], "n_a": sizes[names[i]], "n_b": sizes[names[j]],
             "observed": int(obs[i, j]), "expected": round(float(mean[i, j]), 1), "ratio": round(float(ratio[i, j]), 3),
             "log2_ratio": round(float(lr[i, j]), 3) if np.isfinite(lr[i, j]) else None,
             "z": round(float(z[i, j]), 2), "p": float(p[i, j]), "fdr": float(fdr[i, j])} for i, j in zip(*iu)]
    table = pd.DataFrame(rows).sort_values("z", ascending=False, kind="stable").reset_index(drop=True)
    r2 = lambda m: [[None if not np.isfinite(v) else round(float(v), 3) for v in row] for row in m]
    top = table[table.group_a != table.group_b]
    pick = lambda t: t[["group_a", "group_b", "observed", "expected", "ratio", "z", "fdr"]].to_dict("records")
    say(1.0, "done")
    result = {"by": by, "radius_um": radius, "n_perm": n_perm, "groups": names, "n_cells": sizes,
              "n_pairs": int(len(a)), "within": within or {}, "scope": "project" if len(c.sessions) > 1 else "slice",
              "z": r2(z), "log2_ratio": r2(lr), "observed": r2(obs),
              "expected": r2(mean), "together": pick(top[(top.z > 0) & (top.fdr < 0.05)].head(15)),
              "apart": pick(top[(top.z < 0) & (top.fdr < 0.05)].sort_values("z", kind="stable").head(15)),
              "note": NOTE}
    return result, table
