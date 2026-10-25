"""Import one sample (tissue section) into the data folder, from either source:

  h5ad    celltypeagent import --name 147328 --h5ad a.h5ad [b.h5ad ...] --he slide.ome.tif --registration reg.csv
                               [--section-column donor_id --section 147328] [--region-column histopathology]
  Xenium  celltypeagent import-xenium --name ln1 --xenium /path/to/output-XETG... [--he he.ome.tif --registration
                               he_imagealignment.csv]

h5ad: AnnData files with raw counts (layer "counts", or X) and cell coordinates in micrometers (obs
  x_centroid/y_centroid, x/y, or obsm["spatial"]). Several files (e.g. one per lineage) are merged; cell type
  annotations in them are ignored.
Xenium: the output folder of the Xenium instrument / Xenium Ranger: cells.parquet or cells.csv.gz, and
  cell_feature_matrix.h5 or cell_feature_matrix/ (only "Gene Expression" features are kept). The background image is
  the H&E when given (with its Xenium Explorer alignment file), otherwise the DAPI morphology image of the folder,
  otherwise none.
H&E: pyramidal TIFF / OME-TIFF.
registration: either an affine CSV with columns a11,a12,a21,a22,tx,ty mapping cell (x, y) in um to H&E level-0 pixels
  (px = a11*x + a12*y + tx, py = a21*x + a22*y + ty; with several rows, the row whose other columns contain the
  section value is used), or a Xenium Explorer alignment matrix (3x3 CSV mapping H&E pixels to Xenium morphology
  pixels of --xenium-px-um).
Output: $CELLTYPEAGENT_HOME/samples/<name>/ with cells.parquet, counts.npz, genes.txt and meta.json.

The cells, the image and the registration can come one at a time: a sample needs only its cells; set_image and
set_registration add the others later (in any order), and place_cells puts the cells on the image with an affine
found otherwise (e.g. by align.py). Until an image is registered the cells are shown on a blank background (or the
Xenium DAPI image) and the image waits in meta["image_pending"]. Moving the cells to a new placement moves the
outlines drawn on the old one with them (regions, proposed areas).
"""
import json, time
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp
from scipy.spatial import cKDTree
from . import config

COORD_COLUMNS = [("x_centroid", "y_centroid"), ("x", "y"), ("center_x", "center_y"), ("x_location", "y_location")]
DEPTH_MARKERS, DEPTH_MIN_COUNTS, DEPTH_MIN_FRAC, DEPTH_K, DEPTH_CORE = ("EPCAM", "KRT8"), 3, 0.05, 30, 0.3
BLANK = np.array([[2.0, 0, 0], [0, 2.0, 0]])  # placement without an image: 0.5 um per pixel of a blank background
MIN_INSIDE = 0.9  # share of cells that must fall inside the image for a registration to be accepted


# ---------------- reading AnnData (h5ad) with h5py: only the selected rows are loaded ----------------
def _strings(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def _read_column(g):
    """An obs/var column: plain dataset, categorical group (codes + categories), or nullable group (values + mask,
    written by recent anndata for strings, integers and booleans)."""
    if isinstance(g, h5py.Group):
        if "codes" in g:
            cats = _strings(g["categories"][:]) if g["categories"].dtype.kind in "SOU" else g["categories"][:]
            codes = g["codes"][:]
            out = np.empty(len(codes), dtype=object)
            out[codes >= 0] = np.asarray(cats, dtype=object)[codes[codes >= 0]]
            return out
        v = g["values"][:]
        out = (_strings(v) if v.dtype.kind in "SOU" else v).astype(object)
        if "mask" in g:
            out[g["mask"][:].astype(bool)] = None
        return out
    v = g[:]
    return _strings(v) if v.dtype.kind in "SOU" else v


def _index(g):
    key = g.attrs.get("_index", "_index")
    return np.array([str(v) for v in _read_column(g[key])], dtype=object)


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


def inspect_h5ad(paths, max_values=60):
    """For the import dialog: cells per file, coordinates found, raw-count layer, and obs columns that could name
    the section or a region annotation (few distinct values) with their values."""
    files, cols = [], {}
    for p in paths:
        with h5py.File(p, "r") as h:
            obs = h["obs"]
            n = len(_index(obs))
            pair = next(([a, b] for a, b in COORD_COLUMNS if a in obs and b in obs), None)
            files.append({"path": str(p), "n_cells": int(n),
                          "coordinates": pair or ("obsm['spatial']" if "obsm" in h and "spatial" in h["obsm"] else None),
                          "counts_layer": "counts" if "layers" in h and "counts" in h["layers"] else "X",
                          "n_genes": int(len(_index(h["var"])))})
            for c in obs.keys():
                if c in ("_index", obs.attrs.get("_index")) or c.startswith("__"):
                    continue
                g = obs[c]
                if isinstance(g, h5py.Group) and "categories" in g:  # categorical: its categories are enough
                    cat = g["categories"]
                    vals = [str(v) for v in (_strings(cat[:]) if cat.dtype.kind in "SOU" else cat[:])]
                elif isinstance(g, h5py.Group) or g.dtype.kind in "SOUiu":  # strings or integer ids
                    vals = sorted({str(v) for v in _read_column(g)})
                else:
                    continue
                if 1 <= len(vals) <= max_values:  # a single value can still pick the registration row
                    cols.setdefault(c, set()).update(vals)
    return {"files": files, "n_cells": sum(f["n_cells"] for f in files),
            "columns": {c: sorted(v) for c, v in sorted(cols.items())}}


# ---------------- registration ----------------
def _rows_matching(df, value):
    return df[df.astype(str).apply(lambda col: col == str(value)).any(axis=1)]


def load_registration(path, section=None, xenium_px_um=0.2125, candidates=(), progress=None):
    """Return a 2x3 affine mapping cell (x, y) in um to H&E level-0 pixels. An affine CSV with one row per section
    needs the section; without one, a candidate value (e.g. the h5ad's donor id) that matches exactly one row is used."""
    text = Path(path).read_text().strip()
    first = text.splitlines()[0].lower()
    if "a11" in first:
        df = pd.read_csv(path)
        if len(df) > 1:
            if section is None:
                found = {str(v): _rows_matching(df, v) for v in candidates}
                found = {v: h for v, h in found.items() if len(h) == 1}
                if len(found) != 1:
                    raise ValueError(f"the registration file has {len(df)} rows (one per section): choose the section "
                                     "(Import dialog: Section column and value; command line: --section-column and "
                                     "--section)")
                section = next(iter(found))
                if progress:
                    progress(f"registration: using the row of section {section}")
            hit = _rows_matching(df, section)
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
    from .render import slide
    s = slide(str(path))
    return s.W0, s.H0


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
def _check_name(name, overwrite):
    name = str(name).strip()
    if not name or any(ch in name for ch in "/\\ "):
        raise ValueError("the sample name must be a short id without spaces or slashes")
    if (config.SAMPLES / name).exists() and not overwrite:
        raise ValueError(f"a sample named {name!r} already exists: choose another name, or remove all samples "
                         "(Remove all… next to Import) to start over")
    return name


def _write_bundle(name, ids, xy, X, genes, image, A, regions, description, source, progress, image_kind,
                  registration=None, image_pending=None):
    """Map cells onto the background image, check the registration, estimate depth and write the bundle.
    registration: the affine of a registration file (kept even before there is an image to use it on)."""
    out = config.SAMPLES / name
    pxy = xy @ A[:, :2].T + A[:, 2]
    if image:
        W, H = _he_size(image)
        inside = ((pxy[:, 0] >= 0) & (pxy[:, 0] < W) & (pxy[:, 1] >= 0) & (pxy[:, 1] < H)).mean()
        if inside < MIN_INSIDE:
            raise ValueError(f"only {inside:.0%} of the cells fall inside the image ({W} x {H} px): check the "
                             "registration file, the section, and that coordinates are in micrometers")
        progress(f"registration: {inside:.1%} of cells inside the image")
    um_per_px = float(1 / np.sqrt(abs(np.linalg.det(A[:, :2]))))
    progress(f"image: {um_per_px:.3f} um/px ({image_kind})")
    progress("estimating distance to the epithelium")
    depth = epithelial_depth(X, genes, xy)
    cells = pd.DataFrame({"cell": np.arange(len(ids)), "cell_id": ids, "x": xy[:, 0], "y": xy[:, 1],
                          "px": pxy[:, 0], "py": pxy[:, 1], "depth_um": depth,
                          "total_counts": np.asarray(X.sum(1)).ravel().astype(np.int64)})
    has_regions = False
    if regions is not None:
        cells["pathologist_region"] = np.asarray(regions).astype(str)
        has_regions = bool((~cells.pathologist_region.isin(["", "nan", "None", "Unlabelled", "unknown"])).any())
    out.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(out / "cells.parquet")
    sp.save_npz(out / "counts.npz", X)
    (out / "genes.txt").write_text("\n".join(genes) + "\n")
    meta = {"sample_id": name, "description": description, "n_cells": int(len(ids)), "n_genes": len(genes),
            "he_path": str(Path(image).resolve()) if image else None, "he_um_per_px": um_per_px,
            "image_kind": image_kind, "bundle_dir": str(out), "has_pathologist_regions": has_regions,
            "created_at": time.strftime("%Y-%m-%dT%H:%M:%S"), "source": source, "placement_affine": A.tolist(),
            "registration_affine": None if registration is None else np.asarray(registration).tolist(),
            "registration_by": "file" if registration is not None else None,
            "image_pending": str(Path(image_pending).resolve()) if image_pending else None}
    json.dump(meta, open(out / "meta.json", "w"), indent=1)
    progress(f"imported sample {name} into {out}")
    return meta


def import_sample(name, h5ad, he=None, registration=None, section_column=None, section=None, x_column=None, y_column=None,
                  counts_layer=None, region_column=None, xenium_px_um=0.2125, description="", overwrite=False,
                  progress=print):
    name = _check_name(name, overwrite)
    h5ad = [h5ad] if isinstance(h5ad, (str, Path)) else list(h5ad)
    parts = []
    for p in h5ad:
        progress(f"reading {p}")
        parts.append(read_h5ad(p, section_column, section, x_column, y_column, counts_layer, region_column))
        if not len(parts[-1]["cell_id"]):
            progress(f"  no cells of section {section!r} in {p}")
    parts = [q for q in parts if len(q["cell_id"])]
    if not parts:
        raise ValueError("no cells selected; check the section column and value")
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
    candidates = ()
    if section is None:  # values that are the same for every cell (e.g. the donor id) can pick the registration row
        candidates = [v[0] for v in inspect_h5ad(h5ad)["columns"].values() if len(v) == 1]
    R = load_registration(registration, section, xenium_px_um, candidates, progress) if registration else None
    regions = np.concatenate([q["regions"] for q in parts]) if region_column else None
    if he and R is None:
        progress("  the image waits for a registration (Import > Registration, or Align with Claude)")
    image, A, kind = (he, R, "H&E") if he and R is not None else (None, BLANK, "none")
    return _write_bundle(name, ids, xy, X, genes, image, A, regions, description,
                         {"kind": "h5ad", "h5ad": [str(Path(p).resolve()) for p in h5ad], "section_column": section_column,
                          "section": section, "registration": str(Path(registration).resolve()) if registration else None,
                          "region_column": region_column, "candidates": [str(c) for c in candidates]}, progress, kind,
                         registration=R, image_pending=he if he and R is None else None)


# ---------------- Xenium output folder ----------------
MORPHOLOGY = ["morphology_focus/ch0000_dapi.ome.tif", "morphology_focus.ome.tif", "morphology_mip.ome.tif"]


def _xenium_cells(d):
    if (d / "cells.parquet").exists():
        c = pd.read_parquet(d / "cells.parquet", columns=["cell_id", "x_centroid", "y_centroid"])
    elif (d / "cells.csv.gz").exists():
        c = pd.read_csv(d / "cells.csv.gz", usecols=["cell_id", "x_centroid", "y_centroid"])
    else:
        raise ValueError(f"{d}: no cells.parquet or cells.csv.gz (is this a Xenium output folder?)")
    c["cell_id"] = [x.decode() if isinstance(x, bytes) else str(x) for x in c.cell_id]
    return c


def _xenium_matrix(d):
    """(barcodes, genes, cells x genes CSR) with only "Gene Expression" features."""
    if (d / "cell_feature_matrix.h5").exists():
        with h5py.File(d / "cell_feature_matrix.h5", "r") as h:
            m = h["matrix"]
            n_feat, n_cells = (int(v) for v in m["shape"][:])
            X = sp.csr_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=(n_cells, n_feat))
            bc = _strings(m["barcodes"][:])
            names, ftype = _strings(m["features/name"][:]), _strings(m["features/feature_type"][:])
    elif (d / "cell_feature_matrix").is_dir():
        import scipy.io
        md = d / "cell_feature_matrix"
        X = sp.csr_matrix(scipy.io.mmread(md / "matrix.mtx.gz").T)
        bc = pd.read_csv(md / "barcodes.tsv.gz", header=None)[0].astype(str).to_numpy()
        f = pd.read_csv(md / "features.tsv.gz", header=None, sep="\t")
        names, ftype = f[1].astype(str).to_numpy(), f[2].astype(str).to_numpy()
    else:
        raise ValueError(f"{d}: no cell_feature_matrix.h5 or cell_feature_matrix/ folder")
    keep = np.where(ftype == "Gene Expression")[0]
    if not len(keep):
        raise ValueError("the matrix has no 'Gene Expression' features")
    return bc, [str(g) for g in names[keep]], X[:, keep].astype(np.float32).tocsr()


def inspect_xenium(xenium_dir):
    """What an import would use: counts, panel size, pixel size and the morphology image found."""
    d = Path(xenium_dir).expanduser()
    info = {"folder": str(d), "pixel_size_um": 0.2125, "morphology": None}
    if (d / "experiment.xenium").exists():
        e = json.load(open(d / "experiment.xenium"))
        info["pixel_size_um"] = float(e.get("pixel_size", 0.2125))
        info.update({k: e[k] for k in ("panel_name", "panel_num_targets_predesigned", "region_name", "run_name")
                     if k in e})
    info["morphology"] = next((str(d / m) for m in MORPHOLOGY if (d / m).exists()), None)
    c = _xenium_cells(d)
    info["n_cells"] = int(len(c))
    return info


def import_xenium(name, xenium_dir, he=None, registration=None, description="", overwrite=False, progress=print):
    name = _check_name(name, overwrite)
    d = Path(xenium_dir).expanduser()
    info = inspect_xenium(d)
    progress(f"reading {d}")
    cells = _xenium_cells(d)
    bc, genes, X = _xenium_matrix(d)
    pos = pd.Series(np.arange(len(bc)), index=bc)
    missing = ~cells.cell_id.isin(pos.index)
    if missing.any():
        progress(f"  {int(missing.sum()):,} cells have no counts and are left out")
        cells = cells[~missing]
    X = X[pos[cells.cell_id].to_numpy()]
    xy = cells[["x_centroid", "y_centroid"]].to_numpy(float)
    progress(f"{len(cells):,} cells, {len(genes)} genes")
    ps = info["pixel_size_um"]
    R = load_registration(registration, None, ps) if registration else None
    if he and R is not None:
        A, image, kind = R, he, "H&E"
    elif info["morphology"]:
        A, image, kind = np.array([[1 / ps, 0, 0], [0, 1 / ps, 0]]), info["morphology"], "DAPI"
    else:
        progress("  no H&E and no morphology image: the viewer shows the cells on a blank background")
        A, image, kind = BLANK, None, "none"
    if he and R is None:
        progress("  the H&E waits for its alignment (Import > Registration, or Align with Claude)")
    return _write_bundle(name, cells.cell_id.to_numpy(), xy, X, genes, image, A, None, description,
                         {"kind": "xenium", "folder": str(d), "registration": str(Path(registration).resolve())
                          if registration else None, **{k: v for k, v in info.items() if k != "folder"}},
                         progress, kind, registration=R, image_pending=he if he and R is None else None)


def apply_registration_fix(sample_id, correction):
    """Move every cell onto its nucleus on the image: p + A (p - origin) + b in H&E level-0 pixels, the correction
    that insights.registration fitted. The positions as imported are kept (px_import, py_import) for
    undo_registration_fix. Regions keep their cells (they store cell ids); only the cells' place on the image moves."""
    d = config.SAMPLES / sample_id
    cells = pd.read_parquet(d / "cells.parquet")
    if "px_import" not in cells:
        cells["px_import"], cells["py_import"] = cells.px.to_numpy(), cells.py.to_numpy()
    A, b, o = np.asarray(correction["A"], float), np.asarray(correction["b"], float), np.asarray(correction["origin"], float)
    p = cells[["px", "py"]].to_numpy()
    q = p + (p - o) @ A.T + b
    cells["px"], cells["py"] = q[:, 0], q[:, 1]
    cells.to_parquet(d / "cells.parquet")
    meta = json.load(open(d / "meta.json"))
    meta.setdefault("registration_fixes", []).append({"correction_px": correction,
                                                      "time": time.strftime("%Y-%m-%dT%H:%M:%S")})
    json.dump(meta, open(d / "meta.json", "w"), indent=1)
    return {"moved_um_median": round(float(np.median(np.hypot(*(q - p).T)) * meta.get("he_um_per_px", 1)), 1),
            "moved_um_max": round(float(np.hypot(*(q - p).T).max() * meta.get("he_um_per_px", 1)), 1)}


def undo_registration_fix(sample_id):
    """Put the cells back where the import placed them."""
    d = config.SAMPLES / sample_id
    cells = pd.read_parquet(d / "cells.parquet")
    if "px_import" not in cells:
        raise ValueError("the registration of this sample was never corrected")
    cells["px"], cells["py"] = cells.px_import.to_numpy(), cells.py_import.to_numpy()
    cells.drop(columns=["px_import", "py_import"]).to_parquet(d / "cells.parquet")
    meta = json.load(open(d / "meta.json"))
    meta.pop("registration_fixes", None)
    json.dump(meta, open(d / "meta.json", "w"), indent=1)
    return {"ok": True}


# ---------------- adding the image and the registration later ----------------
def _meta(sample_id):
    return json.load(open(config.SAMPLES / sample_id / "meta.json"))


def placement(sample_id):
    """The affine that currently maps cell (x, y) um to the shown image's level-0 pixels (fitted from the cells for
    samples imported before it was stored)."""
    m = _meta(sample_id)
    if m.get("placement_affine"):
        return np.asarray(m["placement_affine"], float)
    c = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["x", "y", "px", "py"])
    X = np.c_[c.x, c.y, np.ones(len(c))]
    return np.linalg.lstsq(X, c[["px", "py"]].to_numpy(), rcond=None)[0].T


def _move_outlines(sample_id, M):
    """Apply the 2x3 map M (old pixels -> new pixels) to every outline stored in pixels for this sample: drawn
    regions in every version of every session, proposed areas and compartments."""
    def mv(poly):
        P = np.asarray(poly, float)
        return (P @ M[:, :2].T + M[:, 2]).tolist() if len(P) else poly
    for d in config.SESSIONS.iterdir() if config.SESSIONS.exists() else []:
        sj = d / "session.json"
        if not sj.exists() or json.load(open(sj)).get("sample_id") != sample_id:
            continue
        for v in (d / "versions").glob("v*.json"):
            st = json.load(open(v))
            ops = st.get("regions", {}).get("ops", [])
            if any(o.get("polygon") for o in ops):
                for o in ops:
                    if o.get("polygon"):
                        o["polygon"] = mv(o["polygon"])
                json.dump(st, open(v, "w"), indent=1)
        f = d / "advice.json"
        if f.exists():
            a = json.load(open(f))
            for c in a.get("compartments", {}).values():
                c["polygons_px"] = [mv(p) for p in c["polygons_px"]]
            for it in a.get("items", []):
                if it.get("action") and it["action"].get("polygons_px"):
                    it["action"]["polygons_px"] = [mv(p) for p in it["action"]["polygons_px"]]
            json.dump(a, open(f, "w"), indent=1)
        f = d / "area_suggestions.json"
        if f.exists():
            a = json.load(open(f))
            for it in a:
                if it.get("outline_px"):
                    it["outline_px"] = mv(it["outline_px"])
            json.dump(a, open(f, "w"), indent=1)


def place_cells(sample_id, A, image, image_kind="H&E", by="file", progress=print):
    """Put the cells on `image` with the affine A (cell um -> image level-0 pixels): new pixel positions, outlines
    moved along, the previous placement kept for undo_placement. Refuses an A that puts fewer than MIN_INSIDE of
    the cells inside the image."""
    d = config.SAMPLES / sample_id
    A = np.asarray(A, float)
    cells = pd.read_parquet(d / "cells.parquet")
    xy = cells[["x", "y"]].to_numpy()
    pxy = xy @ A[:, :2].T + A[:, 2]
    if image:
        W, H = _he_size(image)
        inside = ((pxy[:, 0] >= 0) & (pxy[:, 0] < W) & (pxy[:, 1] >= 0) & (pxy[:, 1] < H)).mean()
        if inside < MIN_INSIDE:
            raise ValueError(f"only {inside:.0%} of the cells fall inside the image ({W} x {H} px) with this "
                             "registration: check the file, the section and that coordinates are in micrometers")
        progress(f"registration: {inside:.1%} of cells inside the image")
    m = _meta(sample_id)
    old = placement(sample_id)
    M = np.c_[A[:, :2] @ np.linalg.inv(old[:, :2]), A[:, 2] - A[:, :2] @ np.linalg.inv(old[:, :2]) @ old[:, 2]]
    _move_outlines(sample_id, M)
    m["placement_previous"] = {k: m.get(k) for k in ("placement_affine", "he_path", "image_kind", "he_um_per_px",
                                                     "image_pending", "registration_affine", "registration_by")}
    m["placement_previous"]["placement_affine"] = old.tolist()
    cells["px"], cells["py"] = pxy[:, 0], pxy[:, 1]
    cells.drop(columns=[c for c in ("px_import", "py_import") if c in cells]).to_parquet(d / "cells.parquet")
    m.pop("registration_fixes", None)
    m.update(placement_affine=A.tolist(), he_path=str(Path(image).resolve()) if image else None,
             image_kind=image_kind if image else "none", he_um_per_px=float(1 / np.sqrt(abs(np.linalg.det(A[:, :2])))),
             registration_by=by)
    if image and m.get("image_pending") and Path(m["image_pending"]).resolve() == Path(image).resolve():
        m["image_pending"] = None
    json.dump(m, open(d / "meta.json", "w"), indent=1)
    return {"placed_on": image, "by": by, "um_per_px": round(m["he_um_per_px"], 4)}


def undo_placement(sample_id):
    """Back to the placement before the last place_cells."""
    m = _meta(sample_id)
    prev = m.get("placement_previous")
    if not prev:
        raise ValueError("there is no earlier placement to go back to")
    r = place_cells(sample_id, prev["placement_affine"], prev.get("he_path"), prev.get("image_kind") or "H&E",
                    by=prev.get("registration_by") or "file", progress=lambda *a: None)
    m = _meta(sample_id)
    m.update(image_pending=prev.get("image_pending"), registration_affine=prev.get("registration_affine"),
             placement_previous=None)
    json.dump(m, open(config.SAMPLES / sample_id / "meta.json", "w"), indent=1)
    return r


def set_image(sample_id, path, kind="H&E", progress=print):
    """Add (or replace) the sample's image. With a registration already known the cells move onto it at once;
    otherwise it waits (image_pending) for set_registration or an alignment."""
    path = str(Path(path).expanduser().resolve())
    _he_size(path)  # readable?
    m = _meta(sample_id)
    if m.get("registration_affine"):
        return place_cells(sample_id, m["registration_affine"], path, kind, by=m.get("registration_by") or "file",
                           progress=progress)
    m["image_pending"] = path
    m["image_pending_kind"] = kind
    json.dump(m, open(config.SAMPLES / sample_id / "meta.json", "w"), indent=1)
    return {"pending": path}


def set_registration(sample_id, path, section=None, progress=print):
    """Add the registration file; with an image (pending or shown) the cells move onto it at once."""
    m = _meta(sample_id)
    before = {k: m.get(k) for k in ("registration_affine", "registration_by")}
    src = m.get("source") or {}
    R = load_registration(Path(path).expanduser(), section if section else src.get("section"),
                          float(src.get("pixel_size_um") or 0.2125), src.get("candidates") or (), progress)
    m["registration_affine"], m["registration_by"] = R.tolist(), "file"
    src["registration"] = str(Path(path).expanduser().resolve())
    m["source"] = src
    json.dump(m, open(config.SAMPLES / sample_id / "meta.json", "w"), indent=1)
    image = m.get("image_pending") or (m.get("he_path") if m.get("image_kind") == "H&E" else None)
    if not image:
        return {"stored": True, "note": "import the image to see the cells on it"}
    r = place_cells(sample_id, R, image, m.get("image_pending_kind") or "H&E", by="file", progress=progress)
    m = _meta(sample_id)  # undo goes back to before this file, registration included
    m["placement_previous"].update(before)
    json.dump(m, open(config.SAMPLES / sample_id / "meta.json", "w"), indent=1)
    return r
