"""Annotation from a single-cell reference, in the spirit of CellTypist.

A multinomial logistic regression is trained on the labeled cells of a reference h5ad (at most CAP per label), on
the genes the reference shares with the sample's panel, and predicts every cell of the sample. Both sides are
normalized on the shared genes only (counts per cell scaled to the same total, then log1p), so a
whole-transcriptome reference and a panel of a few hundred genes are comparable. Raw counts are read from
layers["counts"], raw.X or X; an X that is already log-normalized is turned back with expm1. Gene names are matched,
case-insensitively, against the var index or the var column that shares the most genes with the panel (e.g.
feature_name when the index holds Ensembl ids). An optional group column (e.g. coarse cell classes) gives each label
a lineage. Optional majority voting (as in CellTypist): the sample is over-clustered and every cluster takes the
mean of its cells' probabilities, so a cluster agrees on one label.
"""
import numpy as np, pandas as pd, scipy.sparse as sp, h5py
from .importer import _index, _read_column, _rows_matrix, as_h5ad

CAP = 3000        # reference cells per label used for training
MIN_CELLS = 10    # labels with fewer reference cells are left out
MAX_LABELS = 300  # obs columns with more distinct values are not offered as label columns
MIN_SHARED = 20   # shared genes needed


def _matrix_source(h):
    """(matrix path, var group) holding raw counts if the file has them."""
    if "layers" in h and "counts" in h["layers"]:
        return "layers/counts", "var"
    if "raw" in h and "X" in h["raw"] and "var" in h["raw"]:
        return "raw/X", "raw/var"
    return "X", "var"


def _gene_names(var):
    """Candidate gene-name arrays: the var index and every var column of strings."""
    out = {"index": _index(var)}
    for c in var.keys():
        if c in ("_index", var.attrs.get("_index")) or c.startswith("__"):
            continue
        try:
            v = _read_column(var[c])
        except Exception:
            continue
        if len(v) and isinstance(v[0], str):
            out[c] = np.array([str(x) for x in v], dtype=object)
    return out


def _match_genes(var, panel):
    """Best naming of the reference genes: (source, panel indices, reference indices) of the shared genes."""
    up = {g.upper(): i for i, g in enumerate(panel)}
    best = ("index", {})
    for src, names in _gene_names(var).items():
        pos = {}
        for j, g in enumerate(names):
            i = up.get(str(g).upper())
            if i is not None and i not in pos:
                pos[i] = j
        if len(pos) > len(best[1]):
            best = (src, pos)
    src, pos = best
    pidx = np.array(sorted(pos), dtype=int)
    return src, pidx, np.array([pos[i] for i in pidx], dtype=int)


def _label_columns(obs):
    cols = {}
    for c in obs.keys():
        if c in ("_index", obs.attrs.get("_index")) or c.startswith("__"):
            continue
        g = obs[c]
        if isinstance(g, h5py.Group) and "categories" in g:
            cat = g["categories"]
            vals = [str(v) for v in (_read_column(cat) if isinstance(cat, h5py.Group) else
                                     (cat[:].astype(str) if cat.dtype.kind in "SOU" else cat[:]))]
        elif isinstance(g, h5py.Dataset) and g.dtype.kind in "SOU":
            vals = sorted({str(v) for v in _read_column(g)})
        else:
            continue
        if 2 <= len(vals) <= MAX_LABELS:
            cols[c] = {"n_values": len(vals), "examples": [str(v) for v in vals[:6]]}
    return cols


def inspect(path, panel):
    """For the Reference dialog: cells, label column candidates, gene naming and overlap with the panel."""
    with h5py.File(as_h5ad(path), "r") as h:
        mpath, vg = _matrix_source(h)
        src, pidx, _ = _match_genes(h[vg], panel)
        return {"path": str(path), "n_cells": int(len(_index(h["obs"]))), "matrix": mpath, "gene_names": src,
                "n_shared_genes": int(len(pidx)), "n_panel_genes": len(panel),
                "missing_examples": [panel[i] for i in sorted(set(range(len(panel))) - set(pidx.tolist()))[:8]],
                "label_columns": _label_columns(h["obs"])}


def scan_labels(path, label_column, group_column=None):
    """Labels kept for training (>= MIN_CELLS cells) with their cell counts, and each label's group (lineage)."""
    path = as_h5ad(path)
    with h5py.File(path, "r") as h:
        obs = h["obs"]
        if label_column not in obs:
            raise ValueError(f"the reference has no obs column {label_column!r}")
        labels = np.array([str(x) for x in _read_column(obs[label_column])], dtype=object)
        groups = np.array([str(x) for x in _read_column(obs[group_column])], dtype=object) if group_column else None
    ok = ~pd.Series(labels).isin(["", "nan", "None", "NaN"]).to_numpy()
    vc = pd.Series(labels[ok]).value_counts()
    vc = vc[vc >= MIN_CELLS]
    if len(vc) < 2:
        raise ValueError(f"column {label_column!r} has fewer than two labels with at least {MIN_CELLS} cells")
    types = sorted(vc.index.tolist())
    lineage = {t: t for t in types}
    if groups is not None:
        tab = pd.crosstab(labels[ok], groups[ok])
        lineage = {t: str(tab.loc[t].idxmax()) for t in types}
    return labels, types, {t: int(vc[t]) for t in types}, lineage


def load_training(path, label_column, panel, seed=0, progress=None):
    """Training cells of the reference: counts on the shared genes (cells x shared panel genes), label indices."""
    say = progress or (lambda f, m: None)
    labels, types, natural, _ = scan_labels(path, label_column)
    rng = np.random.default_rng(seed)
    rows = np.sort(np.concatenate([rng.choice(np.where(labels == t)[0], min(CAP, natural[t]), replace=False)
                                   for t in types]))
    say(0.3, f"reading {len(rows):,} reference cells")
    path = as_h5ad(path)
    with h5py.File(path, "r") as h:
        mpath, vg = _matrix_source(h)
        _, pidx, ridx = _match_genes(h[vg], panel)
        if len(pidx) < MIN_SHARED:
            raise ValueError(f"the reference shares only {len(pidx)} genes with the panel (need {MIN_SHARED})")
        n_vars = len(_index(h[vg]))
        X = sp.csr_matrix(_rows_matrix(h, mpath, rows, n_vars)[:, ridx], dtype=np.float32)
    if X.nnz and not np.allclose(X.data[:20000], np.round(X.data[:20000])):
        X.data = np.expm1(X.data)  # log-normalized values: back to (normalized) counts
    y = np.array([types.index(t) for t in labels[rows]], dtype=np.int64)
    return {"X": X, "y": y, "types": types, "natural": np.array([natural[t] for t in types], float),
            "panel_idx": pidx}


def over_resolution(n):
    """Leiden resolution of CellTypist-like over-clustering, by number of cells."""
    for limit, res in ((5000, 5), (20000, 10), (40000, 15), (100000, 20), (200000, 25)):
        if n < limit:
            return res
    return 30


def majority_vote(logp, clusters):
    """Every cell gets the mean probabilities of its over-cluster (log), so a cluster agrees on one label."""
    P = np.exp(logp.astype(np.float32))
    k = int(clusters.max()) + 1
    one = sp.csr_matrix((np.ones(len(clusters), np.float32), (clusters, np.arange(len(clusters)))), shape=(k, len(clusters)))
    mean = np.asarray(one @ P) / np.maximum(np.bincount(clusters, minlength=k), 1)[:, None]
    return np.log(np.maximum(mean[clusters], 1e-12)).astype(np.float16)
