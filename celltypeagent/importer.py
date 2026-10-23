"""Import one sample (tissue section) from h5ad files, a registered H&E image and a registration file.

    celltypeagent import --name 147328 --h5ad a.h5ad [b.h5ad ...] --he slide.ome.tif --registration reg.csv
                         [--section-column donor_id --section 147328] [--region-column histopathology]

Inputs
  h5ad          AnnData files with raw counts (layer "counts", or X) and cell coordinates in micrometers
                (obs x_centroid/y_centroid, x/y, or obsm["spatial"]). Several files (e.g. one per lineage) are merged;
                cell type annotations in them are ignored.
  H&E           pyramidal TIFF / OME-TIFF registered to the cells.
  registration  either an affine CSV with columns a11,a12,a21,a22,tx,ty mapping cell (x, y) in um to H&E level-0
                pixels (px = a11*x + a12*y + tx, py = a21*x + a22*y + ty; with several rows, the row whose other
                columns contain the section value is used), or a Xenium Explorer alignment matrix (3x3 CSV mapping
                H&E pixels to Xenium morphology pixels of --xenium-px-um).
Output: $CELLTYPEAGENT_HOME/samples/<name>/ with cells.parquet, counts.npz, genes.txt and meta.json.
"""
import json, time
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.spatial import cKDTree
from . import config

COORD_COLUMNS = [("x_centroid", "y_centroid"), ("x", "y"), ("center_x", "center_y"), ("x_location", "y_location")]
DEPTH_MARKERS, DEPTH_MIN_COUNTS, DEPTH_MIN_FRAC, DEPTH_K, DEPTH_CORE = ("EPCAM", "KRT8"), 3, 0.05, 30, 0.3


# ---------------- reading AnnData (h5ad) with h5py: only the selected rows are loaded ----------------
def _strings(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def _read_column(g):
    """An obs/var column: plain dataset or categorical group (codes + categories)."""
    if isinstance(g, h5py.Group):
        cats = _strings(g["categories"][:]) if g["categories"].dtype.kind in "SOU" else g["categories"][:]
        codes = g["codes"][:]
        out = np.empty(len(codes), dtype=object)
        out[codes >= 0] = np.asarray(cats, dtype=object)[codes[codes >= 0]]
        return out
    v = g[:]
    return _strings(v) if v.dtype.kind in "SOU" else v


def _index(g):
    key = g.attrs.get("_index", "_index")
    return _strings(g[key][:])


def _rows_csr(g, rows, n_cols):
    """Rows of a CSR group, read as contiguous runs."""
    indptr = g["indptr"][:]
    runs = np.split(rows, np.where(np.diff(rows) != 1)[0] + 1) if len(rows) else []
    data, ind, ptr = [], [], [0]
    for r in runs:
        a, b = indptr[r[0]], indptr[r[-1] + 1]
        data.append(g["data"][a:b]); ind.append(g["indices"][a:b])
        ptr.extend((indptr[r[0] + 1:r[-1] + 2] - a + ptr[-1]).tolist())
    m = sp.csr_matrix((np.concatenate(data) if data else [], np.concatenate(ind) if ind else [], ptr),
                      shape=(len(rows), n_cols))
    return m


def _rows_matrix(h, path, rows, n_cols):
    g = h[path]
    if isinstance(g, h5py.Group):
        enc = g.attrs.get("encoding-type", "")
        if "csc" in str(enc):
            full = sp.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g.attrs["shape"]))
            return full.tocsr()[rows]
        return _rows_csr(g, rows, n_cols)
    return sp.csr_matrix(g[np.sort(rows)] if len(rows) else np.zeros((0, n_cols)))


def read_h5ad(path, section_column=None, section=None, x_column=None, y_column=None, counts_layer=None,
              region_column=None):
    with h5py.File(path, "r") as h:
        names = _index(h["obs"])
        genes = _index(h["var"])
        rows = np.arange(len(names))
        if section_column:
            if section_column not in h["obs"]:
                raise ValueError(f"{path}: obs has no column {section_column!r}")
            col = _read_column(h["obs"][section_column])
            rows = np.where(np.array([str(v) for v in col]) == str(section))[0]
        if x_column and y_column:
            xy = np.c_[_read_column(h["obs"][x_column])[rows], _read_column(h["obs"][y_column])[rows]].astype(float)
        else:
            pair = next(((a, b) for a, b in COORD_COLUMNS if a in h["obs"] and b in h["obs"]), None)
            if pair:
                xy = np.c_[_read_column(h["obs"][pair[0]])[rows], _read_column(h["obs"][pair[1]])[rows]].astype(float)
            elif "obsm" in h and "spatial" in h["obsm"]:
                xy = h["obsm"]["spatial"][:][rows, :2].astype(float)
            else:
                raise ValueError(f"{path}: no cell coordinates found; pass --x-column / --y-column")
        layer = counts_layer or ("counts" if "layers" in h and "counts" in h["layers"] else None)
        X = _rows_matrix(h, f"layers/{layer}" if layer else "X", rows, len(genes)).astype(np.float32)
        regions = _read_column(h["obs"][region_column])[rows] if region_column else None
    if X.nnz and not np.allclose(X.data[:10000], np.round(X.data[:10000])):
        raise ValueError(f"{path}: the matrix is not raw counts; pass --counts-layer with the layer holding counts")
    return {"cell_id": names[rows], "xy": xy, "X": sp.csr_matrix(X), "genes": list(genes), "regions": regions}


# ---------------- registration ----------------
def load_registration(path, section=None, xenium_px_um=0.2125):
    """Return a 2x3 affine mapping cell (x, y) in um to H&E level-0 pixels."""
    text = Path(path).read_text().strip()
    first = text.splitlines()[0].lower()
    if "a11" in first:
        df = pd.read_csv(path)
        if len(df) > 1:
            if section is None:
                raise ValueError("the registration file has several rows; pass --section to choose one")
            hit = df[df.astype(str).apply(lambda col: col == str(section)).any(axis=1)]
            if len(hit) != 1:
                raise ValueError(f"no single registration row matches section {section!r}")
            df = hit
        r = df.iloc[0]
        return np.array([[r.a11, r.a12, r.tx], [r.a21, r.a22, r.ty]], float)
    M = np.loadtxt(path, delimiter="," if "," in first else None)
    if M.shape != (3, 3):
        raise ValueError("registration must be an a11..ty CSV or a 3x3 matrix")
    inv = np.linalg.inv(M)  # Xenium Explorer: H&E px -> morphology px; invert, and convert um to morphology px
    S = np.diag([1 / xenium_px_um, 1 / xenium_px_um, 1.0])
    return (inv @ S)[:2]


def _he_size(path):
    import tifffile
    with tifffile.TiffFile(path) as t:
        lv = t.series[0].levels[0].shape
    return (lv[2], lv[1]) if lv[0] == 3 else (lv[1], lv[0])


# ---------------- depth (distance to the epithelium) ----------------
def epithelial_depth(X, genes, xy):
    up = {g.upper(): i for i, g in enumerate(genes)}
    idx = [up[g] for g in DEPTH_MARKERS if g in up]
    if not idx:
        return np.full(len(xy), np.nan)
    tot = np.asarray(X.sum(1)).ravel()
    c = np.asarray(X[:, idx].sum(1)).ravel()
    pos = (c >= DEPTH_MIN_COUNTS) & (c / np.maximum(tot, 1) >= DEPTH_MIN_FRAC)
    tree = cKDTree(xy)
    _, nn = tree.query(xy, k=min(DEPTH_K + 1, len(xy)))
    core = pos[nn[:, 1:]].mean(1) >= DEPTH_CORE
    if not core.any():
        return np.full(len(xy), np.nan)
    return cKDTree(xy[core]).query(xy)[0]


# ---------------- import ----------------
def import_sample(name, h5ad, he, registration, section_column=None, section=None, x_column=None, y_column=None,
                  counts_layer=None, region_column=None, xenium_px_um=0.2125, description="", overwrite=False,
                  progress=print):
    name = str(name).strip()
    if not name or any(ch in name for ch in "/\\\\ "):
        raise ValueError("--name must be a short id without spaces or slashes")
    out = config.SAMPLES / name
    if out.exists() and not overwrite:
        raise ValueError(f"sample {name!r} already exists (use --overwrite to replace it)")
    h5ad = [h5ad] if isinstance(h5ad, (str, Path)) else list(h5ad)
    parts = []
    for p in h5ad:
        progress(f"reading {p}")
        parts.append(read_h5ad(p, section_column, section, x_column, y_column, counts_layer, region_column))
        if not len(parts[-1]["cell_id"]):
            progress(f"  no cells of section {section!r} in {p}")
    parts = [q for q in parts if len(q["cell_id"])]
    if not parts:
        raise ValueError("no cells selected; check --section-column / --section")
    genes = parts[0]["genes"]
    for q in parts[1:]:  # align gene order across files
        if q["genes"] != genes:
            pos = {g: i for i, g in enumerate(q["genes"])}
            missing = [g for g in genes if g not in pos]
            if missing:
                raise ValueError(f"h5ad files have different panels (missing {missing[:5]}...)")
            q["X"] = q["X"][:, [pos[g] for g in genes]]
    X = sp.vstack([q["X"] for q in parts]).tocsr()
    xy = np.vstack([q["xy"] for q in parts])
    ids = np.concatenate([q["cell_id"] for q in parts])
    if len(set(ids)) != len(ids):
        raise ValueError("cell ids are not unique across the h5ad files")
    progress(f"{len(ids):,} cells, {len(genes)} genes")

    A = load_registration(registration, section, xenium_px_um)
    pxy = xy @ A[:, :2].T + A[:, 2]
    W, H = _he_size(he)
    inside = ((pxy[:, 0] >= 0) & (pxy[:, 0] < W) & (pxy[:, 1] >= 0) & (pxy[:, 1] < H)).mean()
    if inside < 0.9:
        raise ValueError(f"only {inside:.0%} of the cells fall inside the H&E image ({W} x {H} px): check the "
                         "registration file, the section, and that coordinates are in micrometers")
    um_per_px = float(1 / np.sqrt(abs(np.linalg.det(A[:, :2]))))
    progress(f"registration: {inside:.1%} of cells inside the H&E; H&E {um_per_px:.3f} um/px")

    progress("estimating distance to the epithelium")
    depth = epithelial_depth(X, genes, xy)
    cells = pd.DataFrame({"cell": np.arange(len(ids)), "cell_id": ids, "x": xy[:, 0], "y": xy[:, 1],
                          "px": pxy[:, 0], "py": pxy[:, 1], "depth_um": depth,
                          "total_counts": np.asarray(X.sum(1)).ravel().astype(np.int64)})
    has_regions = False
    if region_column:
        cells["pathologist_region"] = np.concatenate([q["regions"] for q in parts]).astype(str)
        has_regions = bool((~cells.pathologist_region.isin(["", "nan", "None", "Unlabelled", "unknown"])).any())
    out.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(out / "cells.parquet")
    sp.save_npz(out / "counts.npz", X)
    (out / "genes.txt").write_text("\n".join(genes) + "\n")
    meta = {"sample_id": name, "description": description, "n_cells": int(len(ids)), "n_genes": len(genes),
            "he_path": str(Path(he).resolve()), "he_um_per_px": um_per_px, "bundle_dir": str(out),
            "has_pathologist_regions": has_regions, "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "source": {"h5ad": [str(Path(p).resolve()) for p in h5ad], "section_column": section_column,
                       "section": section, "registration": str(Path(registration).resolve()),
                       "region_column": region_column}}
    json.dump(meta, open(out / "meta.json", "w"), indent=1)
    progress(f"imported sample {name} into {out}")
    return meta
