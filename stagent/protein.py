"""Protein / immunofluorescence stains measured per cell (CosMx: the Mean.<stain> columns of the cell metadata, e.g.
PanCK, CD45, CD3, CD68; Xenium: the "Protein Expression" features of the cell feature matrix) as evidence of the
cell types, besides the transcripts.

No stain needs to be named for a type: each stain is made positive / negative per cell (Otsu's threshold on
log1p(intensity)), each cell type's positive rate is read from the transcript-based annotation itself (smoothed
toward the overall rate), and every cell's probabilities of the types are multiplied by the likelihood of its stains
under each type's rates, tempered by WEIGHT, then renormalized. The rates are taken at two levels per stain, split at
the widest gap between the types' rates, and only when that gap is clear (see _two_levels): a stain spilling into
neighbors gives immune or stromal types rates of 0.25-0.55 that would otherwise move cells between them (on ASHUC1
3,322 cells, 92% PanCK+, into "Mast cells", a type half PanCK+ from epithelial cells labeled so), and a stain whose
rates run on without a gap (a membrane stain, CD45 or CD3 on ASHUC1 and lung5) would split types alike in it by a
fixed line (ILCs at 0.67 vs CD8 T cells at 0.79 swapped 2,300 cells by Membrane alone). Nuclear and membrane stains
(DAPI, Membrane...) are left out by name; PanCK moves PanCK-negative cells away from the types that are mostly
PanCK-positive. Stains are read once from the import's files and kept as protein.parquet (cell_id +
one column per stain) in the sample's folder, so filters and corrections of the cells do not matter (matched by id).
"""
import json
import numpy as np, pandas as pd
from . import config

WEIGHT = 0.5          # the stains' log-likelihood is multiplied by this (spillover of a stain into neighbors)
SMOOTH = 20.0         # pseudo-cells of the overall rate in each type's rate
RATE_CLIP = (0.02, 0.98)
MIN_CELLS_POS = 0.01  # a stain positive in fewer (or more than 1 - this) of the cells is left out
SKIP = ("dapi", "membrane")  # nuclear and membrane stains (names starting so, any case): on in every cell type
MIN_GAP = 0.2         # a stain tells types apart only where their positive rates have a gap at least this wide
MIN_N = 200           # types with fewer cells (rates less sure; 1% of the cells when that is less, as in a small
                      # sample) do not make or fill a gap; they join the nearer side


def skipped(name):
    """A stain of no type information (nuclear, membrane)."""
    return str(name).strip().lower().startswith(SKIP)


def _file(sample_id):
    return config.SAMPLES / str(sample_id) / "protein.parquet"


def _from_cosmx(source):
    from .importer import _cosmx_id, read_csv
    f = (source.get("files") or {}).get("cells")
    if not f:
        return None
    try:
        head = read_csv(f, nrows=0).columns.tolist()
    except OSError:
        return None
    stains = [c for c in head if c.startswith("Mean.")]
    if not stains:
        return None
    low = {c.lower(): c for c in head}
    use = [low[c] for c in ("fov", "cell_id") if c in low]
    d = read_csv(f, usecols=use + stains)
    cid, ok = _cosmx_id(d)
    d = d[ok.to_numpy()]
    out = pd.DataFrame({s[5:]: d[s].to_numpy(float) for s in stains})
    out.insert(0, "cell_id", np.asarray(cid[ok.to_numpy()]).astype(str))
    return out


def _from_xenium(source):
    import h5py, scipy.sparse as sp
    from .importer import _strings
    files = source.get("files") or {}
    if "atera" in files and "matrix" not in files:  # (Atera: protein features of its AnnData zarr matrix, if any)
        from .importer import _atera_read
        obs, var, X = _atera_read(files["atera"], keep_types=("Protein Expression",))
        if "feature_type" not in var.columns or not len(var):
            return None
        names = (var["feature_name"] if "feature_name" in var.columns else var.index).astype(str).to_numpy()
        out = pd.DataFrame(X.toarray().astype(float), columns=[str(x) for x in names])
        out.insert(0, "cell_id", [str(b) for b in obs.index])
        return out
    f = files.get("matrix")
    if not f:
        return None
    try:
        with h5py.File(f, "r") as h:
            m = h["matrix"]
            names, ftype = _strings(m["features/name"][:]), _strings(m["features/feature_type"][:])
            keep = np.flatnonzero(ftype == "Protein Expression")
            if not len(keep):
                return None
            n_feat, n_cells = (int(v) for v in m["shape"][:])
            X = sp.csr_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(n_cells, n_feat))[:, keep]
            bc = _strings(m["barcodes"][:])
    except OSError:
        return None
    out = pd.DataFrame(X.toarray().astype(float), columns=[str(x) for x in names[keep]])
    out.insert(0, "cell_id", [str(b) for b in bc])
    return out


def read(source):
    """The stains of an import's files (source: the sample's meta["source"]: kind and files), or None."""
    source = source or {}
    kind = source.get("kind")
    return _from_cosmx(source) if kind == "cosmx" else _from_xenium(source) if kind == "xenium" else None


def ensure(sample_id):
    """The sample's stains (DataFrame: cell_id + one column per stain), read once from the import's files; None when
    it has none (or its files are gone)."""
    f = _file(sample_id)
    none = f.with_suffix(".none")
    if f.exists():
        return pd.read_parquet(f)
    if none.exists():
        return None
    meta = json.load(open(config.SAMPLES / str(sample_id) / "meta.json"))
    try:
        d = read(meta.get("source"))
    except Exception as e:  # unreadable: no stains (said once)
        print(f"stains of {sample_id} not read: {e}")
        d = None
    if d is None or d.shape[1] < 2:
        if (meta.get("source") or {}).get("files"):  # (files there but no stains: do not look again)
            none.write_text("")
        return None
    from .osutil import replace_atomic
    replace_atomic(f, lambda t: d.to_parquet(t, index=False))
    return d


def available(sample_id):
    """Names of the sample's stains that can be used (without the nuclear ones), [] if none."""
    try:
        d = ensure(sample_id)
    except Exception:
        return []
    return [] if d is None else [c for c in d.columns[1:] if not skipped(c)]


def _otsu(v):
    h, e = np.histogram(v, 128)
    c = (e[:-1] + e[1:]) / 2
    w0 = np.cumsum(h)
    w1 = w0[-1] - w0
    m0 = np.cumsum(h * c) / np.maximum(w0, 1)
    m1 = ((h * c).sum() - np.cumsum(h * c)) / np.maximum(w1, 1)
    return float(c[np.argmax(w0 * w1 * (m0 - m1) ** 2)])


def positives(sample_id, cell_ids, stains=None):
    """(stain names, cells x stains bool positives, finite mask) for the given cells (by id)."""
    return positives_df(ensure(sample_id), cell_ids, stains)


def positives_df(d, cell_ids, stains=None):
    """positives() of a stains table (protein.parquet or read())."""
    if d is None:
        return [], None, None
    d = d.set_index(d.cell_id.astype(str)).reindex(pd.Index(np.asarray(cell_ids).astype(str)))
    names = [c for c in d.columns[1:] if not skipped(c) and (stains is None or c in stains)]
    V = np.log1p(np.clip(d[names].to_numpy(float), 0, None))
    ok = np.isfinite(V).all(1)
    P = np.zeros(V.shape, bool)
    keep = []
    for j, n in enumerate(names):
        v = V[ok, j]
        if not len(v) or v.max() <= v.min():
            continue
        P[:, j] = V[:, j] > _otsu(v)
        frac = P[ok, j].mean()
        if MIN_CELLS_POS <= frac <= 1 - MIN_CELLS_POS:
            keep.append(j)
    return [names[j] for j in keep], P[:, keep], ok


def _two_levels(R, n):
    """Each stain separates the types it is on in from the others at the widest gap between the types' rates (types
    of at least MIN_N cells; the cut in the gap's middle), each side given its cells' mean rate. Differences within a
    side are mostly a stain spilling into neighbors (a PanCK rate of 0.25 or 0.55 in immune or stromal types, by how
    close they lie to the epithelium) and would move cells between those types; a stain without a gap of MIN_GAP
    (rates running on from type to type) tells nothing: any line through it splits types alike in the stain."""
    R = R.copy()
    have = n > 0
    big = n >= min(MIN_N, max(10.0, 0.01 * n.sum()))
    for j in range(R.shape[1]):
        r = np.sort(R[big, j])
        g = np.diff(r)
        if not len(g) or g.max() < MIN_GAP:
            R[:, j] = np.average(R[have, j], weights=n[have]) if have.any() else 0.5
            continue
        k = int(g.argmax())
        hi = R[:, j] > (r[k] + r[k + 1]) / 2
        for side in (hi, ~hi):
            m = side & have
            if m.any():
                R[side, j] = np.average(R[m, j], weights=n[m])
    return R


def apply(logp, labels, P, ok, weight=WEIGHT, two_level=True):
    """logp (cells x types) times the stains' likelihood under each type's positive rates (from labels)."""
    T = logp.shape[1]
    lp = logp.astype(np.float32).copy()
    if P is None or not P.shape[1]:
        return lp, {}
    lab = np.asarray(labels)
    m = ok & (lab >= 0)
    overall = P[m].mean(0)
    n = np.bincount(lab[m], minlength=T).astype(float)
    pos = np.stack([np.bincount(lab[m], weights=P[m, j], minlength=T) for j in range(P.shape[1])], 1)
    R = np.clip((pos + SMOOTH * overall) / (n[:, None] + SMOOTH), *RATE_CLIP)   # types x stains
    if two_level:
        R = _two_levels(R, n)
    ll = P[ok].astype(np.float32) @ np.log(R).T + (~P[ok]).astype(np.float32) @ np.log(1 - R).T
    lp[ok] += weight * (ll - ll.max(1, keepdims=True))
    lp -= lp.max(1, keepdims=True)
    lp -= np.log(np.exp(lp).sum(1, keepdims=True))
    return lp, {"rates": R, "overall": overall}
