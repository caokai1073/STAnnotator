"""Visium HD and Visium (Space Ranger output): bins (or spots) instead of cells.

Space Ranger writes, for Visium HD, binned_outputs/square_002um, square_008um, square_016um (and custom sizes); each
holds filtered_feature_bc_matrix.h5 (or the folder filtered_feature_bc_matrix/) with the bins under the tissue, and
spatial/: tissue_positions.parquet (barcode, in_tissue, array_row, array_col, pxl_row_in_fullres,
pxl_col_in_fullres), scalefactors_json.json (bin_size_um, microns_per_pixel, tissue_hires_scalef) and
tissue_hires_image.png. Visium (55 um spots) has the same files in outs/ (tissue_positions.csv or
tissue_positions_list.csv; the pixel size from the spot diameter).

The analyses run on bins as on cells, at BIN_UM: finer bins (2 um, or 4 um) are pooled into bins of BIN_UM by their
grid position (counts summed; the position is the mean of the pooled bins under the tissue); bins of BIN_UM or more
(16 um, Visium spots) are used as they are. Positions are pixels of the full-resolution image Space Ranger was given,
so that image needs no registration; without it the output's tissue_hires_image.png is shown.
"""
import json, re
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp

BIN_UM = 8          # the analyses work on bins of this size (finer bins are pooled into it)
SPOT_UM = 55.0      # a Visium spot's diameter
CHUNK = 1_000_000   # bins read at a time when pooling


KINDS = ("visium_hd", "visium")
NO_REGISTRATION_CHECK = ("the sample comes from Visium bins on Space Ranger's grid, placed on the image by Space Ranger's "
                         "own alignment (cells are built from those bins): the nuclei check does not apply")


def is_bins(meta):
    """True for a sample placed by its pipeline on a grid of spots: Space Ranger (bins, spots, or cells built from
    its bins) or Stereo-seq bins; the nuclei registration check does not apply to them."""
    src = meta.get("source") or {}
    return src.get("kind") in KINDS or (src.get("kind") == "stereoseq" and src.get("unit") == "bins")


def _size(d):
    m = re.fullmatch(r"square_(\d+)um", Path(d).name)
    return int(m.group(1)) if m else None


def _is_bins(d):
    d = Path(d)
    return d.is_dir() and (d / "spatial").is_dir() and ((d / "filtered_feature_bc_matrix.h5").exists()
                                                       or (d / "filtered_feature_bc_matrix").is_dir())


def folders(path):
    """The bin folders under a Space Ranger output (outs/, binned_outputs/ or one square_XXum folder):
    [(path, bin size in um or None for Visium spots)]."""
    p = Path(path).expanduser()
    if p.is_file():
        p = p.parent
    if p.name == "spatial":
        p = p.parent
    for cand in (p, p / "outs"):
        if (cand / "binned_outputs").is_dir():
            cand = cand / "binned_outputs"
        found = sorted(((d, _size(d)) for d in cand.glob("square_*um") if _is_bins(d)), key=lambda x: x[1])
        if found:
            return found
        if _is_bins(cand):
            return [(cand, _size(cand))]
    raise ValueError(f"{path}: no Space Ranger output found (binned_outputs/square_XXXum/ with "
                     "filtered_feature_bc_matrix.h5 and spatial/, or a Visium outs/ folder)")


def choose(found):
    """The bin folder to read: BIN_UM bins when there are, else the finest that pools into BIN_UM, else the finest.
    A Visium output (spots) has one."""
    sizes = {s: d for d, s in found}
    if BIN_UM in sizes:
        return sizes[BIN_UM], BIN_UM
    finer = [s for s in sizes if s and s < BIN_UM]
    if finer:
        s = max(finer)
        return sizes[s], s
    d, s = found[0]
    return d, s


def _scalefactors(d):
    f = Path(d) / "spatial" / "scalefactors_json.json"
    if not f.exists():
        raise ValueError(f"{f} is missing")
    return json.load(open(f))


def _positions(d):
    sd = Path(d) / "spatial"
    if (sd / "tissue_positions.parquet").exists():
        t = pd.read_parquet(sd / "tissue_positions.parquet")
    elif (sd / "tissue_positions.csv").exists():
        t = pd.read_csv(sd / "tissue_positions.csv")
    elif (sd / "tissue_positions_list.csv").exists():
        t = pd.read_csv(sd / "tissue_positions_list.csv", header=None,
                        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    else:
        raise ValueError(f"{sd}: no tissue_positions file")
    return t


def scale(d, s):
    """(bin size in um, um per full-resolution pixel)."""
    sf = _scalefactors(d)
    if s is None:  # Visium spots: the pixel size from the spot diameter
        return SPOT_UM, SPOT_UM / float(sf["spot_diameter_fullres"])
    return float(sf.get("bin_size_um", s)), float(sf["microns_per_pixel"])


def target(s):
    """The bin size the analyses use for bins of size s (um)."""
    if s >= BIN_UM:
        return s
    f = int(np.ceil(BIN_UM / s))
    return s * f


def _strings(a):
    return np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)


def inspect(path, segmentation=None):
    """What an import would read: the bin folders found, the one used, its size and the size analysed, bins and
    genes; the cell outlines (Space Ranger 4's, or `segmentation`) that make cells the unit."""
    found = folders(path)
    d, s = choose(found)
    size, mpp = scale(d, s)
    out = {"folder": str(d), "key": _key(d), "found_um": [x[1] for x in found], "bin_um": size,
           "analysis_um": target(size), "um_per_px": mpp, "visium_spots": s is None,
           "hires_image": str(Path(d) / "spatial" / "tissue_hires_image.png")
           if (Path(d) / "spatial" / "tissue_hires_image.png").exists() else None}
    m = Path(d) / "filtered_feature_bc_matrix.h5"
    if m.exists():
        with h5py.File(m, "r") as h:
            n_feat, n_bc = (int(v) for v in h["matrix/shape"][:])
            ft = _strings(h["matrix/features/feature_type"][:])
        out.update(n_bins=n_bc, n_genes=int((ft == "Gene Expression").sum()))
    if out["analysis_um"] != size:
        f = int(round(out["analysis_um"] / size))
        out["n_analysis_bins_about"] = int(out.get("n_bins", 0) / (f * f) * 1.1) if out.get("n_bins") else None
    seg = segmentation or find_segmentation(path)
    out["segmentation"] = str(seg) if seg else None
    out["fine_bins"] = next((str(fd) for fd, fs in sorted(found, key=lambda x: x[1] or 99) if fs and fs <= 4), None)
    if seg and Path(seg).exists():
        txt = open(seg).read(4_000_000)
        out["n_outlines_about"] = txt.count('"Feature"') if len(txt) < 4_000_000 else int(Path(seg).stat().st_size / len(txt) * txt.count('"Feature"'))
    return out


def _key(d):
    """A name for the sample from the folders (the output's folder above outs/ or binned_outputs/)."""
    parts = [x for x in Path(d).resolve().parts if x not in ("outs", "binned_outputs") and not re.fullmatch(r"square_\d+um", x)]
    return parts[-1] if parts else "visium"


def read(d, s, progress=print):
    """(bin ids, xy in full-resolution pixels (x, y), genes, bins x genes CSR, bin size analysed in um, um per
    full-resolution pixel) of a bin folder; finer bins than BIN_UM are pooled."""
    size, mpp = scale(d, s)
    want = target(size)
    f = int(round(want / size))
    pos = _positions(d).set_index("barcode")
    m = Path(d) / "filtered_feature_bc_matrix.h5"
    if not m.exists():
        return _read_mtx(d, pos, f, want, mpp, progress)
    with h5py.File(m, "r") as h:
        g = h["matrix"]
        n_feat, n_bc = (int(v) for v in g["shape"][:])
        ft = _strings(g["features/feature_type"][:])
        names = _strings(g["features/name"][:])
        keep = np.where(ft == "Gene Expression")[0]
        if not len(keep):
            raise ValueError("the matrix has no 'Gene Expression' features")
        progress(f"{n_bc:,} bins of {size:g} um under the tissue, {len(keep):,} genes")
        bc = _strings(g["barcodes"][:])
        at = pos.index.get_indexer(bc)
        if (at < 0).any():
            raise ValueError(f"{int((at < 0).sum()):,} bins of the matrix are not in the tissue positions file")
        row, col = pos.array_row.to_numpy()[at], pos.array_col.to_numpy()[at]
        px, py = pos.pxl_col_in_fullres.to_numpy(float)[at], pos.pxl_row_in_fullres.to_numpy(float)[at]
        del pos
        indptr = g["indptr"][:]
        gene_of = np.full(n_feat, -1, np.int64)
        gene_of[keep] = np.arange(len(keep))
        if f == 1:
            X = sp.csr_matrix((g["data"][:].astype(np.float32), g["indices"][:], indptr), shape=(n_bc, n_feat))[:, keep]
            ids = bc
        else:
            # pooled bin of every fine bin, by grid position
            key = (row // f).astype(np.int64) * (int(col.max()) // f + 1) + col // f
            uniq, pid = np.unique(key, return_inverse=True)
            n_pool = len(uniq)
            progress(f"pooling them into {n_pool:,} bins of {want:g} um")
            X = None
            for a in range(0, n_bc, CHUNK):
                b = min(n_bc, a + CHUNK)
                lo, hi = int(indptr[a]), int(indptr[b])
                part = sp.csr_matrix((g["data"][lo:hi].astype(np.float32), g["indices"][lo:hi], indptr[a:b + 1] - lo),
                                     shape=(b - a, n_feat))[:, keep]
                M = sp.csr_matrix((np.ones(b - a, np.float32), (pid[a:b], np.arange(b - a))), shape=(n_pool, b - a))
                X = M @ part if X is None else X + M @ part
                progress(f"  {b:,} of {n_bc:,} bins read")
            n_in = np.bincount(pid, minlength=n_pool)
            px, py = np.bincount(pid, px, n_pool) / n_in, np.bincount(pid, py, n_pool) / n_in
            ncol = int(col.max()) // f + 1
            ids = np.array([f"p{want:03g}um_{u // ncol:05d}_{u % ncol:05d}" for u in uniq], dtype=object)
    X = sp.csr_matrix(X, dtype=np.float32)
    tot = np.asarray(X.sum(1)).ravel()
    ok = tot > 0
    if not ok.all():
        progress(f"  {int((~ok).sum()):,} bins without counts left out")
    return ids[ok], np.c_[px, py][ok], [str(x) for x in names[keep]], X[ok], want, mpp


def _read_mtx(d, pos, f, want, mpp, progress):
    if f != 1:
        raise ValueError("pooling fine bins needs filtered_feature_bc_matrix.h5")
    import scipy.io
    md = Path(d) / "filtered_feature_bc_matrix"
    X = sp.csr_matrix(scipy.io.mmread(md / "matrix.mtx.gz").T, dtype=np.float32)
    bc = pd.read_csv(md / "barcodes.tsv.gz", header=None)[0].astype(str).to_numpy()
    ft = pd.read_csv(md / "features.tsv.gz", header=None, sep="\t")
    keep = np.where(ft[2].astype(str).to_numpy() == "Gene Expression")[0]
    at = pos.index.get_indexer(bc)
    if (at < 0).any():
        raise ValueError(f"{int((at < 0).sum()):,} spots of the matrix are not in the tissue positions file")
    X = X[:, keep].tocsr()
    xy = np.c_[pos.pxl_col_in_fullres.to_numpy(float)[at], pos.pxl_row_in_fullres.to_numpy(float)[at]]
    progress(f"{len(bc):,} spots under the tissue, {len(keep):,} genes")
    ok = np.asarray(X.sum(1)).ravel() > 0
    return bc[ok], xy[ok], [str(x) for x in ft[1].astype(str).to_numpy()[keep]], X[ok], want, mpp


# ---------------- background images: a pyramidal copy the viewer can read at every zoom ----------------
def pyramid(src, dst, um_per_px=None, progress=print, tile=512):
    """dst: a tiled, pyramidal OME-TIFF copy of the RGB image src (a PNG / JPEG, or a TIFF / BigTIFF with one level,
    as the full-resolution image Space Ranger is given often is); returns src itself when it is pyramidal already.
    Written strip by strip, the lower levels through temporary files next to dst, so a whole-slide image needs
    little memory."""
    import tempfile, tifffile, zarr
    src = str(src)
    if src.lower().endswith((".tif", ".tiff", ".btf")):
        with tifffile.TiffFile(src) as t:
            s = t.series[0]
            if len(s.levels) > 1:
                return src
            axes, shape = s.axes, s.shape
        z = zarr.open(tifffile.imread(src, aszarr=True), mode="r")
        H, W = shape[axes.index("Y")], shape[axes.index("X")]

        def strip(y0, y1):
            idx = tuple(slice(y0, y1) if c == "Y" else slice(None) if c in "XS" else slice(0, 3) if c == "C" else 0
                        for c in axes)
            a = np.asarray(z[idx])
            kept = [c for c in axes if c in "YXSC"]
            a = np.moveaxis(a, [kept.index("Y"), kept.index("X")], [0, 1])
            a = a.reshape(a.shape[0], a.shape[1], -1)
            if a.dtype != np.uint8:
                a = np.clip(a / max(float(a.max()), 1) * 255, 0, 255).astype(np.uint8)
            return a[..., :3] if a.shape[2] >= 3 else np.repeat(a[..., :1], 3, axis=2)
    else:
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None
        arr = np.asarray(Image.open(src).convert("RGB"))
        H, W = arr.shape[:2]
        strip = lambda y0, y1: arr[y0:y1]
    dims = [(H, W)]
    while max(dims[-1]) > 2 * tile:
        h, w = dims[-1]
        dims.append((h // 2, w // 2))
    progress(f"pyramidal copy of the image ({W:,} x {H:,} px, {len(dims)} levels)")
    tmpd = Path(tempfile.mkdtemp(prefix=".pyramid-", dir=Path(dst).parent))
    try:
        lv = [None] + [np.lib.format.open_memmap(tmpd / f"l{k}.npy", mode="w+", dtype=np.uint8, shape=(h, w, 3))
                       for k, (h, w) in enumerate(dims) if k]

        def tiles(get, k):
            h, w = dims[k]
            nxt = lv[k + 1] if k + 1 < len(dims) else None
            for y0 in range(0, h, tile):
                a = get(y0, min(h, y0 + tile))
                if nxt is not None:  # the next level's rows of this strip
                    hh, ww = a.shape[0] // 2, nxt.shape[1]
                    if hh:
                        b = a[:2 * hh, :2 * ww].reshape(hh, 2, ww, 2, 3).mean((1, 3))
                        r0 = y0 // 2
                        nxt[r0:r0 + hh] = b[:max(0, min(hh, nxt.shape[0] - r0))]
                for x0 in range(0, w, tile):
                    t = np.full((tile, tile, 3), 245, np.uint8)
                    blk = a[:, x0:x0 + tile]
                    t[:blk.shape[0], :blk.shape[1]] = blk
                    yield t
                if k == 0 and (y0 // tile) % 20 == 0:
                    progress(f"  image rows {y0:,} of {h:,}")
        part = Path(str(dst) + ".part")
        kw = dict(tile=(tile, tile), compression="jpeg", compressionargs={"level": 90}, photometric="rgb", dtype=np.uint8)
        meta = {"axes": "YXS"} | ({"PhysicalSizeX": um_per_px, "PhysicalSizeY": um_per_px} if um_per_px else {})
        with tifffile.TiffWriter(part, bigtiff=True, ome=True) as t:
            t.write(tiles(strip, 0), shape=(H, W, 3), subifds=len(dims) - 1, metadata=meta, **kw)
            for k in range(1, len(dims)):
                t.write(tiles(lambda a, b, k=k: lv[k][a:b], k), shape=(*dims[k], 3), subfiletype=1, **kw)
        part.replace(dst)
    finally:
        import shutil
        shutil.rmtree(tmpd, ignore_errors=True)
    return str(dst)


# ---------------- cells from a segmentation: the 2 um bins split between the cells by area ----------------
SUB = 4             # each fine bin is sampled at SUB x SUB points against the cell outlines (area shares in 1/16)
SEG_FILES = ("segmented_outputs/cell_segmentations.geojson",)  # Space Ranger 4


def find_segmentation(path):
    """Space Ranger 4's cell outlines next to the binned output (segmented_outputs/cell_segmentations.geojson)."""
    p = Path(path).expanduser()
    p = p.parent if p.is_file() else p
    for base in (p, *p.parents[:3]):
        for cand in (base, base / "outs"):
            for f in SEG_FILES:
                if (cand / f).exists():
                    return str((cand / f).resolve())
    return None


def fine_bins(found):
    """The finest bin folder (2 um), which cells are built from."""
    fine = [(d, s) for d, s in found if s and s <= 4]
    if not fine:
        raise ValueError("cells from a segmentation need the 2 um bins (binned_outputs/square_002um)")
    return min(fine, key=lambda x: x[1])


def read_outlines(path):
    """Cell outlines of a GeoJSON (Polygon / MultiPolygon features in full-resolution pixels, as Space Ranger 4,
    QuPath, StarDist or Cellpose write them): (ids, [list of rings per cell, each (n, 2) x, y])."""
    gj = json.load(open(path))
    feats = gj["features"] if isinstance(gj, dict) and "features" in gj else gj
    ids, rings = [], []
    for k, f in enumerate(feats):
        g = f.get("geometry") or {}
        polys = [g["coordinates"]] if g.get("type") == "Polygon" else g.get("coordinates", []) if g.get("type") == "MultiPolygon" else []
        outer = [np.asarray(p[0], float) for p in polys if p and len(p[0]) >= 3]
        if not outer:
            continue
        pr = f.get("properties") or {}
        ids.append(str(pr.get("cell_id", pr.get("id", pr.get("name", k + 1)))))
        rings.append(outer)
    if not ids:
        raise ValueError(f"{path}: no cell outlines (Polygon features) found")
    return ids, rings


def _area_centroid(ring):
    x, y = ring[:, 0], ring[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    c = x * y1 - x1 * y
    a = c.sum() / 2
    if abs(a) < 1e-9:
        return 0.0, x.mean(), y.mean()
    return abs(a), ((x + x1) * c).sum() / (6 * a), ((y + y1) * c).sum() / (6 * a)


def _grid(row, col, px, py):
    """The bin grid in image pixels: (origin, step along array_col, step along array_row) by least squares."""
    i = np.random.default_rng(0).choice(len(row), min(len(row), 200000), replace=False)
    M = np.c_[col[i], row[i], np.ones(len(i))]
    cx, *_ = np.linalg.lstsq(M, px[i], rcond=None)
    cy, *_ = np.linalg.lstsq(M, py[i], rcond=None)
    return np.array([cx[2], cy[2]]), np.array([cx[0], cy[0]]), np.array([cx[1], cy[1]])


def read_cells(d, s, seg, progress=print, seed=0):
    """Cells from a segmentation: every fine bin (2 um) is split between the cells its area overlaps, in proportion
    to the area (SUB x SUB points per bin against the outlines rasterized at a quarter bin); each count of the bin
    goes to a cell with the probability of its area share (a seeded binomial split: counts stay whole and their
    expectation is area-proportional), the share outside every cell is left out. Returns (cell ids, centroid xy in
    full-resolution pixels, genes, cells x genes CSR, info)."""
    from PIL import Image, ImageDraw
    size, mpp = scale(d, s)
    rng = np.random.default_rng(seed)
    progress(f"reading the cell outlines {Path(seg).name}")
    ids, rings = read_outlines(seg)
    n_cells = len(ids)
    geo = np.array([_area_centroid(max(r, key=len)) for r in rings])
    area_um2 = np.array([sum(_area_centroid(x)[0] for x in r) for r in rings]) * mpp ** 2
    pos = _positions(d).set_index("barcode")
    m = Path(d) / "filtered_feature_bc_matrix.h5"
    with h5py.File(m, "r") as h:
        g = h["matrix"]
        n_feat, n_bc = (int(v) for v in g["shape"][:])
        ft, names = _strings(g["features/feature_type"][:]), _strings(g["features/name"][:])
        keep = np.where(ft == "Gene Expression")[0]
        bc = _strings(g["barcodes"][:])
        at = pos.index.get_indexer(bc)
        if (at < 0).any():
            raise ValueError(f"{int((at < 0).sum()):,} bins of the matrix are not in the tissue positions file")
        row, col = pos.array_row.to_numpy()[at], pos.array_col.to_numpy()[at]
        px, py = pos.pxl_col_in_fullres.to_numpy(float)[at], pos.pxl_row_in_fullres.to_numpy(float)[at]
        del pos
        o, ec, er = _grid(row, col, px, py)
        # the outlines rasterized at a quarter bin (labels 1..n_cells; 0 = no cell)
        r_px = float(np.hypot(*ec)) / SUB
        x0, y0 = px.min() - 4 * r_px * SUB, py.min() - 4 * r_px * SUB
        W, H = int((px.max() - x0) / r_px) + 4 * SUB + 2, int((py.max() - y0) / r_px) + 4 * SUB + 2
        progress(f"{n_cells:,} cell outlines rasterized at {r_px * mpp:.2f} um ({W:,} x {H:,})")
        img = Image.new("I", (W, H), 0)
        dr = ImageDraw.Draw(img)
        for k, r in enumerate(rings):
            for ring in r:
                dr.polygon([((x - x0) / r_px, (y - y0) / r_px) for x, y in ring], fill=k + 1)
        lab = np.asarray(img, dtype=np.int32)
        del img, dr
        off = (np.arange(SUB) + 0.5) / SUB - 0.5
        du, dv = np.meshgrid(off, off)
        du, dv = du.ravel(), dv.ravel()
        indptr = g["indptr"][:]
        gene_of = np.full(n_feat, -1, np.int64)
        gene_of[keep] = np.arange(len(keep))
        X, n_split, tot_in, tot_all = None, 0, 0.0, 0.0
        progress(f"splitting {n_bc:,} bins of {size:g} um between the cells by area")
        for a in range(0, n_bc, CHUNK // 2):
            b = min(n_bc, a + CHUNK // 2)
            # labels at SUB x SUB points of each bin
            sx = px[a:b, None] + du[None] * ec[0] + dv[None] * er[0]
            sy = py[a:b, None] + du[None] * ec[1] + dv[None] * er[1]
            ix = np.clip(((sx - x0) / r_px).astype(np.int64), 0, W - 1)
            iy = np.clip(((sy - y0) / r_px).astype(np.int64), 0, H - 1)
            L = lab[iy, ix]                                   # (bins, SUB*SUB): the cell at each point, 0 = none
            # per bin: the distinct cells and their area shares
            key, n_pts = np.unique(np.arange(b - a)[:, None] * (n_cells + 1) + L, return_counts=True)
            r_, cell = key // (n_cells + 1), key % (n_cells + 1)
            share = n_pts / L.shape[1]
            ok = cell > 0
            r_, cell, share = r_[ok], cell[ok] - 1, share[ok]
            n_split += int((np.bincount(r_, minlength=b - a) > 1).sum())
            # the bins' counts
            lo, hi = int(indptr[a]), int(indptr[b])
            data = g["data"][lo:hi].astype(np.int64)
            gi = gene_of[g["indices"][lo:hi]]
            brow = np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))
            m_ = gi >= 0
            data, gi, brow = data[m_], gi[m_], brow[m_]
            tot_all += float(data.sum())
            # each (bin, gene) count split between the bin's cells: sequential binomials over its cells
            order = np.argsort(r_, kind="stable")
            r_, cell, share = r_[order], cell[order], share[order]
            start = np.searchsorted(r_, np.arange(b - a))
            nlab = np.bincount(r_, minlength=b - a)
            rem = data.copy()
            used = np.zeros(len(data))
            out_c, out_g, out_v = [], [], []
            for j in range(int(nlab.max()) if len(nlab) else 0):
                has = nlab[brow] > j
                if not has.any():
                    break
                k = start[brow[has]] + j
                p = share[k] / np.maximum(1.0 - used[has], 1e-12)
                take = rng.binomial(rem[has], np.clip(p, 0, 1))
                rem[has] -= take
                used[has] += share[k]
                nz = take > 0
                out_c.append(cell[k][nz]); out_g.append(gi[has][nz]); out_v.append(take[nz])
            if out_c:
                cc, gg, vv = np.concatenate(out_c), np.concatenate(out_g), np.concatenate(out_v)
                tot_in += float(vv.sum())
                part = sp.csr_matrix((vv.astype(np.float32), (cc, gg)), shape=(n_cells, len(keep)))
                X = part if X is None else X + part
            progress(f"  {b:,} of {n_bc:,} bins")
    X = sp.csr_matrix(X if X is not None else (n_cells, len(keep)), dtype=np.float32)
    tot = np.asarray(X.sum(1)).ravel()
    ok = tot > 0
    info = {"n_outlines": n_cells, "n_cells": int(ok.sum()), "bins_split": n_split, "bins": n_bc,
            "counts_in_cells": round(tot_in / max(tot_all, 1), 4), "median_area_um2": round(float(np.median(area_um2[ok])), 1),
            "median_counts": float(np.median(tot[ok])) if ok.any() else 0.0}
    progress(f"{info['n_cells']:,} cells with counts (of {n_cells:,} outlines); {n_split:,} bins split between "
             f"cells; {info['counts_in_cells']:.0%} of the counts fall in a cell")
    return np.array(ids, dtype=object)[ok], geo[ok, 1:], [str(x) for x in names[keep]], X[ok], info
