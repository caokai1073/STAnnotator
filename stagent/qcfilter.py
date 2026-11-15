"""Filtering the cells and genes of a sample, or of all slices of a project together (they keep one gene panel).

Cells with too few or too many transcripts, too few genes detected, or too large a share of mitochondrial transcripts
(genes MT-...; imaging panels rarely have them), and genes detected in too few of the remaining cells are deleted from
the sample's data: cells.parquet, counts.npz and genes.txt hold only what is kept, and the app and the agent see
nothing else. The data as imported is kept once in original/ of the sample folder; every filter starts from it, so
new thresholds filter the imported cells again, and thresholds left empty bring everything back (original/ is then
removed). Annotations and the PCA / UMAP number the cells of a sample: the caller moves them to the trash first.

The cells' places on the image follow the current placement (meta placement_affine), also for the cells brought
back; the distance to the epithelium is estimated again on the cells kept.
"""
import json, re, shutil
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from . import config, importer, osutil, regcheck
from .osutil import dump_json

CELL_PARAMS = ("min_counts", "max_counts", "min_genes", "max_mito_pct")
PARAMS = CELL_PARAMS + ("min_cells",)
MITO = re.compile(r"^mt-", re.I)
_data = {}  # (sample_id, file time) -> the imported data with its per-cell numbers


def _src(sample_id):
    d = config.SAMPLES / sample_id
    return d / "original" if (d / "original" / "counts.npz").exists() else d


STATS = "qc_stats.npz"  # beside a counts file: its per-cell and per-gene numbers, so the dialog need not read it
_light = {}  # (sample_id, counts file time) -> those numbers, without the matrix


def _counts_key(src):
    st = (src / "counts.npz").stat()
    return np.array([st.st_mtime_ns, st.st_size], np.int64)


def _numbers(X, genes):
    """Per cell: transcripts, genes detected, mitochondrial share; per gene: the cells it is detected in (the
    matrix holds no zeros: a stored entry is a detection)."""
    mito = np.array([bool(MITO.match(g)) for g in genes])
    counts = np.asarray(X.sum(1)).ravel()
    return {"genes": genes, "mito": mito, "counts": counts, "n_genes": np.diff(X.indptr),
            "mito_pct": 100 * np.asarray(X[:, mito].sum(1)).ravel() / np.maximum(counts, 1) if mito.any()
            else np.zeros(X.shape[0]),
            "detected": np.bincount(X.indices, minlength=X.shape[1])}


def _load(sample_id):
    """The data as imported with its numbers (the matrix and the cells: what filtering needs). Saves the numbers
    beside the counts (STATS) for the dialog."""
    with _lock(sample_id):  # (the dialog's read in the background and a filter started meanwhile: read once)
        return _load_locked(sample_id)


_locks = {}


def _lock(sample_id):
    import threading
    return _locks.setdefault(sample_id, threading.Lock())


def _load_locked(sample_id):
    src = _src(sample_id)
    key = (sample_id, (src / "counts.npz").stat().st_mtime)
    if key not in _data:
        for k in [k for k in _data if k[0] == sample_id or len(_data) > 2]:
            _data.pop(k)
        from .session import shared_counts
        X = shared_counts(src)  # the copy the open annotations use (a CosMx 6K matrix is gigabytes)
        if X.nnz and not X.data.all():  # stored zeros: dropped on a copy (the shared matrix is never changed)
            X = X.copy()
            X.eliminate_zeros()
        num = _numbers(X, config.read_genes(src / "genes.txt"))
        _data[key] = num | {"X": X, "cells": pd.read_parquet(src / "cells.parquet")}
        try:  # (a folder that cannot be written to: the numbers are computed again next time)
            np.savez(src / STATS, key=_counts_key(src), **{k: num[k] for k in ("counts", "n_genes", "mito_pct", "detected")})
        except OSError:
            pass
    return _data[key]


def save_numbers(folder, X, genes):
    """At import: the numbers of the counts just written in folder, saved beside them (the filter dialog opens without
    reading the matrix)."""
    folder = Path(folder)
    X = sp.csr_matrix(X)
    if not X.has_canonical_format or (X.nnz and not X.data.all()):  # (as _load sees the matrix)
        X = X.copy()
        X.sum_duplicates()
        X.eliminate_zeros()
    num = _numbers(X, list(genes))
    try:
        np.savez(folder / STATS, key=_counts_key(folder), **{k: num[k] for k in ("counts", "n_genes", "mito_pct", "detected")})
    except OSError:
        pass


def _numbers_of(sample_id):
    """The per-cell and per-gene numbers of the data as imported, without reading the matrix when they were saved
    (STATS, checked against the counts file's time and size); else read with the data (_load)."""
    src = _src(sample_id)
    key = (sample_id, (src / "counts.npz").stat().st_mtime)
    if key in _data:
        return _data[key]
    if key not in _light:
        try:
            z = np.load(src / STATS)
            if not np.array_equal(z["key"], _counts_key(src)):
                raise ValueError("stale")
            genes = config.read_genes(src / "genes.txt")
            num = {k: z[k] for k in ("counts", "n_genes", "mito_pct", "detected")}
            if len(num["detected"]) != len(genes):
                raise ValueError("stale")
            _light[key] = num | {"genes": genes, "mito": np.array([bool(MITO.match(g)) for g in genes])}
        except (OSError, KeyError, ValueError):
            return _load(sample_id)
    return _light[key]


def prefetch(ids):
    """Read the data of these samples in the background (the dialog is open: the preview and the filter will need
    it), slices at once (unzipping lets go of the GIL)."""
    import threading
    for s in ids:
        threading.Thread(target=lambda x=s: _quiet(_load, x), daemon=True).start()


def _quiet(fn, *a):
    try:
        fn(*a)
    except Exception:
        pass


def params(p):
    """The thresholds that are set, as numbers ({} for none); empty or missing ones are not applied."""
    out = {}
    for k in PARAMS:
        v = (p or {}).get(k)
        if v is None or str(v).strip() == "":
            continue
        try:
            v = float(v)
        except ValueError:
            raise ValueError(f"{k}: {v!r} is not a number")
        if not np.isfinite(v) or v < 0:
            raise ValueError(f"{k} must be a number of at least 0")
        out[k] = int(v) if k != "max_mito_pct" else v
    if "max_mito_pct" in out and out["max_mito_pct"] > 100:
        raise ValueError("the mitochondrial share is a percentage (0 to 100)")
    if out.get("max_counts") is not None and out.get("min_counts", 0) > out["max_counts"]:
        raise ValueError("the minimum transcripts per cell is above the maximum")
    return out


def _cell_masks(d, prm):
    """Cells kept, and for each threshold the cells it removes (a cell can fail several)."""
    fail = {}
    if "min_counts" in prm:
        fail["min_counts"] = d["counts"] < prm["min_counts"]
    if "max_counts" in prm:
        fail["max_counts"] = d["counts"] > prm["max_counts"]
    if "min_genes" in prm:
        fail["min_genes"] = d["n_genes"] < prm["min_genes"]
    if "max_mito_pct" in prm and d["mito"].any():
        fail["max_mito_pct"] = d["mito_pct"] > prm["max_mito_pct"]
    keep = np.ones(len(d["counts"]), bool)
    for m in fail.values():
        keep &= ~m
    return keep, fail


def _plan(ids, prm, full=False):
    """What the thresholds keep. The matrix is read only when needed: full (to filter), or a minimum of cells per
    gene that some cells removed may change (the genes' detection counted over the cells kept)."""
    data = {s: _numbers_of(s) for s in ids}
    genes = data[ids[0]]["genes"]
    if any(data[s]["genes"] != genes for s in ids[1:]):
        raise ValueError("the slices of this project do not share one gene panel")
    keep = {s: _cell_masks(data[s], prm) for s in ids}
    if full or (prm.get("min_cells", 0) > 0 and not all(k[0].all() for k in keep.values())):
        data = _load_all(ids)
        detected = sum(_detected_in(data[s]["X"], keep[s][0]) for s in ids)
    else:
        detected = sum(data[s]["detected"] for s in ids)
    gkeep = detected >= prm.get("min_cells", 0)
    return data, genes, keep, gkeep, detected


def _detected_in(X, rows):
    """Per gene, the cells among rows (a mask) it is detected in."""
    return np.bincount(X.indices[np.repeat(rows, np.diff(X.indptr))], minlength=X.shape[1])


def _load_all(ids):
    """_load of every slice, at once (unzipping lets go of the GIL)."""
    if len(ids) == 1:
        return {ids[0]: _load(ids[0])}
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(min(4, len(ids))) as ex:
        return dict(zip(ids, ex.map(_load, ids)))


def preview(ids, p):
    """What a filter with these thresholds keeps and removes (from the data as imported)."""
    prm = params(p)
    data, genes, keep, gkeep, _ = _plan(ids, prm)
    return _summary(ids, prm, genes, keep, gkeep)


def _summary(ids, prm, genes, keep, gkeep):
    n = sum(len(k[0]) for k in keep.values())
    kept = sum(int(k[0].sum()) for k in keep.values())
    by = {k: int(sum(int(f[k].sum()) for _, f in keep.values() if k in f)) for k in CELL_PARAMS}
    return {"params": prm, "n_cells": n, "n_cells_kept": kept, "removed_by": {k: v for k, v in by.items() if k in prm},
            "n_genes": len(genes), "n_genes_kept": int(gkeep.sum()),
            "genes_removed": [g for g, k in zip(genes, gkeep) if not k],
            "slices": [{"sample_id": s, "n_cells": int(len(keep[s][0])), "n_cells_kept": int(keep[s][0].sum())} for s in ids]}


def _hist(v, bins):
    h, e = np.histogram(v, bins=bins)
    return {"edges": [round(float(x), 4) for x in e], "counts": h.tolist()}


def overview(ids):
    """The numbers the thresholds act on, over the cells as imported: histograms of transcripts and genes per cell
    (log scale), of the mitochondrial share, and of the cells each gene is detected in; the filter in use."""
    data = {s: _numbers_of(s) for s in ids}
    counts = np.concatenate([data[s]["counts"] for s in ids])
    ngenes = np.concatenate([data[s]["n_genes"] for s in ids])
    d0 = data[ids[0]]
    detected = sum(data[s]["detected"] for s in ids)
    def log(v, n):  # counts on a log scale, log10(1 + v), in bins that each hold whole numbers
        e = np.unique(np.floor(np.geomspace(1, max(float(np.max(v)), 9) + 2, n + 1))) - 1.5
        return _hist(np.log10(1 + np.asarray(v, float)), np.log10(1 + np.maximum(e, -0.5)))
    out = {"n_cells": int(len(counts)), "n_genes": len(d0["genes"]), "slices": ids,
           "mito_genes": [g for g, m in zip(d0["genes"], d0["mito"]) if m],
           "counts": log(counts, 60), "genes": log(ngenes, 60), "cells_per_gene": log(detected, 40),
           "median": {"counts": float(np.median(counts)), "genes": float(np.median(ngenes))},
           "filter": None}
    fs = [json.load(open(config.SAMPLES / s / "meta.json")).get("filter") for s in ids]
    if fs[0]:  # the filter in use, cells removed summed over the slices
        out["filter"] = fs[0] | {"n_cells_removed": sum(f["n_cells_removed"] for f in fs if f)}
    if d0["mito"].any():
        pct = np.concatenate([data[s]["mito_pct"] for s in ids])
        out["mito_pct"] = _hist(pct, np.linspace(0, max(float(pct.max()), 1), 51))
        out["median"]["mito_pct"] = float(np.median(pct))
    return out


def _placed(cells, meta, columns):
    """The cells' pixel positions under the sample's current placement (and the import placement when the sample's
    table has it: the registration was corrected)."""
    cells = cells.drop(columns=[c for c in ("px_import", "py_import") if c in cells])
    xy = cells[["x", "y"]].to_numpy(float)
    for cols, key in ((("px", "py"), "placement_affine"), (("px_import", "py_import"), "placement_affine_import")):
        A = meta.get(key)
        if A and (cols[0] in columns or key == "placement_affine"):
            A = np.asarray(A, float)
            cells[cols[0]], cells[cols[1]] = (xy @ A[:, :2].T + A[:, 2]).T
    return cells


def _corrected(cells, meta, now):
    """Positions of a sample whose registration was corrected: the cells shown now keep their place; cells the filter
    brings back get the corrections replayed (importer.replay_placement). Piece moves stored by row number are
    first stored by cell id (meta is updated), since the filter numbers the rows anew."""
    for f in meta["registration_fixes"]:
        pc = f.get("piece")
        if pc is not None and pc.get("cell_ids") is None:
            pc["cell_ids"] = now.cell_id.to_numpy()[importer.piece_rows(pc, now.cell_id.to_numpy())].tolist()
            pc.pop("cells", None)
    first, pos = importer.replay_placement(meta, cells[["x", "y"]].to_numpy(float), cells.cell_id.to_numpy())
    cells["px"], cells["py"] = pos[:, 0], pos[:, 1]
    cells["px_import"], cells["py_import"] = first[:, 0], first[:, 1]
    at = pd.Index(now.cell_id.astype(str)).get_indexer(cells.cell_id.astype(str))
    have = at >= 0
    for col in ("px", "py", "px_import", "py_import"):
        if col in now:
            cells.loc[have, col] = now[col].to_numpy()[at[have]]
    return cells


def apply(ids, p, progress=lambda *a: None, genes_removed=None):
    """Filter the samples ids with these thresholds (from the data as imported); no thresholds: every cell and gene
    back. genes_removed: the genes to leave out, as a saved filter did (bundle.py: a project's slices counted together
    for min_cells), instead of counting them here. Returns the preview of what was done."""
    prm = params(p)
    progress(0.02, "reading the cells and genes as imported")
    data, genes, keep, gkeep, _ = _plan(ids, prm, full=True)
    summary = _summary(ids, prm, genes, keep, gkeep)   # (from the plan above: the data is not read again)
    if genes_removed is not None:
        out = set(genes_removed)
        gkeep = np.array([g not in out for g in genes])
    if not any(k[0].any() for k in keep.values()):
        raise ValueError("these thresholds remove every cell")
    if not gkeep.any():
        raise ValueError("these thresholds remove every gene")
    if genes_removed is not None:
        summary["genes_removed"] = [g for g in genes if g in set(genes_removed)]
    n = len(ids)
    done = [0.0] * n   # each slice's share done: the slices are filtered at once (zlib, parquet and numpy let go of the GIL)

    def at(i, s, f, what):
        done[i] = f
        progress(0.3 + 0.7 * sum(done) / n, what + (f" ({s}, slice {i + 1} of {n})" if n > 1 else ""))

    def one(i, s):
        at(i, s, 0.0, "keeping the data as imported")
        d = config.SAMPLES / s
        meta = json.load(open(d / "meta.json"))
        meta.setdefault("placement_affine", importer.placement(s).tolist())
        prev = set(pd.read_parquet(d / "cells.parquet").columns)
        o = d / "original"
        if prm and not (o / "counts.npz").exists():  # the data as imported, once
            o.mkdir(exist_ok=True)
            for f in ("cells.parquet", "counts.npz", "genes.txt", STATS):
                if f != STATS or (d / f).exists():  # (the numbers stay valid: copy2 keeps the counts file's time)
                    shutil.copy2(d / f, o / f)
        at(i, s, 0.25, "filtering the cells and genes")
        D, km = data[s], keep[s][0]
        X = D["X"][km][:, gkeep].tocsr()
        cells = D["cells"].iloc[np.flatnonzero(km)].reset_index(drop=True)
        cells = _placed(cells, meta, prev)
        if meta.get("registration_fixes"):  # corrected: local shifts and moved pieces are no affine
            cells = _corrected(cells, meta, pd.read_parquet(d / "cells.parquet"))
        cells["cell"] = np.arange(len(cells))
        if "total_counts" in cells:
            cells["total_counts"] = np.asarray(X.sum(1)).ravel().astype(np.int64)
        g = [x for x, k in zip(genes, gkeep) if k]
        if "depth_um" in cells and prm:
            at(i, s, 0.45, "the distance of each cell to the epithelium")
            cells["depth_um"] = importer.epithelial_depth(X, g, cells[["x", "y"]].to_numpy(float))
        at(i, s, 0.85, "writing the filtered cells")
        image = importer.image_of(s)
        placed = regcheck._key(s, image)   # the placement's checks stay valid: no cell moves (regcheck.carry)
        cells.to_parquet(d / "cells.parquet")
        osutil.save_counts(d / "counts.npz", X)   # (the fastest zip compression: seconds, not tens of seconds)
        config.write_genes(d / "genes.txt", g)
        meta.update(n_cells=int(len(cells)), n_genes=len(g))
        if prm:
            sl = next(x for x in summary["slices"] if x["sample_id"] == s)
            meta["filter"] = {"params": prm, "n_cells_imported": sl["n_cells"], "n_cells_removed": sl["n_cells"] - sl["n_cells_kept"],
                              "n_genes_imported": len(genes), "genes_removed": summary["genes_removed"],
                              "slices": ids if len(ids) > 1 else None}
        else:
            meta.pop("filter", None)
        regcheck.carry(s, image, placed, meta)
        dump_json(meta, d / "meta.json", indent=1)
        if not prm and o.exists():
            shutil.rmtree(o)
        done[i] = 1.0

    if n == 1:
        one(0, ids[0])
    else:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(min(4, n)) as ex:
            for f in [ex.submit(one, i, s) for i, s in enumerate(ids)]:
                f.result()
    progress(1.0, "done")
    return summary
