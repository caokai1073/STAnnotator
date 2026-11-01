"""Stereo-seq (STOmics SAW output): spots are DNBs 0.5 um apart, read from a GEM (text) or a cell-bin GEF (HDF5).

Like Visium HD (visium.py), the cells are the unit when the data has them, else bins of BIN_UM:
  cells  a cell-bin GEF (<SN>.cellbin.gef or .adjusted.cellbin.gef: cellBin/cell, cellExp, gene), or a GEM with a
         CellID column (cell-bin GEM; CellID 0 = outside the cells); a DNB is a point, so it belongs to one cell
  bins   a GEM of the DNBs (<SN>.tissue.gem.gz, bin1): pooled into BIN_UM bins (BIN DNBs a side, counts summed)
GEM: '#key=value' header lines (OffsetX / OffsetY), then a tab-separated table with geneID (and geneName), x, y,
MIDCount (MIDCounts / UMICount), optionally the cell (CellID, or label as in MOSTA). Genes are named by geneName when there is one; Ensembl ids that
share a name are summed. Positions are DNBs (+ the file's offset), which are the pixels of the registered image SAW
writes (<SN>_HE_regist.tif or <SN>_ssDNA_regist.tif / DAPI): no registration needed; another image (e.g. MOSTA's raw
ssDNA scan, not in the DNBs' coordinates) waits for an alignment.
"""
import json, re
from pathlib import Path
import h5py, numpy as np, pandas as pd, scipy.sparse as sp

DNB_UM = 0.5
BIN_UM = 8
BIN = int(BIN_UM / DNB_UM)  # DNBs a side
CHUNK = 5_000_000           # GEM rows read at a time
GEM_ENDS = (".gem", ".gem.gz", ".tsv", ".tsv.gz", ".txt", ".txt.gz")
CELL_COLS = ("cellid", "cell_id", "cell", "label")  # a cell-bin GEM's cell column (MOSTA: label; 0 / empty = no cell)


def _cell_col(cols):
    low = {c.lower(): c for c in cols}
    return next((low[c] for c in CELL_COLS if c in low), None)


def _image_kind(p):
    n = Path(p).name.lower()
    return "H&E" if re.search(r"(^|[_.-])he([_.-]|$)", n) or "h&e" in n else "DAPI"


def find(path):
    """The files of a SAW output (a folder, or files): {"cellbin_gef", "gem", "cell_gem", "image"} (each a path or
    None). A GEM with a CellID column is a cell-bin GEM."""
    paths = [Path(p).expanduser() for p in (path if isinstance(path, (list, tuple)) else [path]) if p]
    files = []
    for p in paths:
        if p.is_dir():
            files += [f for f in p.rglob("*") if f.is_file() and len(f.relative_to(p).parts) <= 3]
        elif p.exists():
            files.append(p)
        else:
            raise ValueError(f"{p} does not exist")
    out = {"cellbin_gef": None, "gem": None, "cell_gem": None, "image": None, "image_registered": False}
    gefs = [f for f in files if f.name.endswith(".cellbin.gef")]
    if gefs:  # the adjusted cell bins (cells expanded by SAW) when there are, else the cell bins
        out["cellbin_gef"] = str(sorted(gefs, key=lambda f: ("adjusted" not in f.name, f.name))[0])
    for f in sorted(files):
        if f.name.lower().endswith(GEM_ENDS) and ("gem" in f.name.lower() or "bin1" in f.name.lower()):
            cols = gem_columns(f)
            if cols and _cell_col(cols):
                out["cell_gem"] = out["cell_gem"] or str(f)
            elif cols:
                out["gem"] = out["gem"] or str(f)
    # SAW's registered images (_regist.tif) are in the DNBs' coordinates; another image (e.g. the raw ssDNA scan) waits
    # for an alignment (Align with Claude)
    imgs = [f for f in files if f.name.lower().endswith((".tif", ".tiff")) and "regist" in f.name.lower()]
    out["image_registered"] = bool(imgs)
    imgs = imgs or [f for f in files if f.name.lower().endswith((".tif", ".tiff")) and "mask" not in f.name.lower()]
    if imgs:  # an H&E before an ssDNA
        out["image"] = str(sorted(imgs, key=lambda f: (_image_kind(f) != "H&E", f.name))[0])
    if not (out["cellbin_gef"] or out["cell_gem"] or out["gem"]):
        raise ValueError(f"{path}: no Stereo-seq expression found (a .cellbin.gef, or a GEM: .gem.gz / .tsv.gz with "
                         "geneID, x, y, MIDCount)")
    return out


def _open(f):
    import gzip
    return gzip.open(f, "rt") if str(f).endswith(".gz") else open(f)


def _header(f):
    """('#key=value' lines as a dict, number of those lines)."""
    meta, n = {}, 0
    with _open(f) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            n += 1
            k, _, v = line[1:].strip().partition("=")
            meta[k.strip()] = v.strip()
    return meta, n


def gem_columns(f):
    try:
        meta, n = _header(f)
        with _open(f) as fh:
            for _ in range(n):
                fh.readline()
            cols = fh.readline().rstrip("\n").split("\t")
    except (OSError, UnicodeDecodeError, EOFError):
        return None
    low = [c.lower() for c in cols]
    return cols if "x" in low and "y" in low and any(c in low for c in ("geneid", "genename")) else None


def _offset(meta):
    return np.array([float(meta.get("OffsetX", 0) or 0), float(meta.get("OffsetY", 0) or 0)])


def inspect(path):
    """What an import would read: the unit (cells or bins), the files and the image."""
    f = find(path)
    out = dict(f, unit="cells" if (f["cellbin_gef"] or f["cell_gem"]) else "bins", bin_um=BIN_UM,
               image_kind=_image_kind(f["image"]) if f["image"] else None)
    if f["cellbin_gef"]:
        with h5py.File(f["cellbin_gef"], "r") as h:
            out["n_cells"] = int(h["cellBin/cell"].shape[0])
            out["n_genes"] = int(h["cellBin/gene"].shape[0])
    if f["image"]:
        import tifffile
        with tifffile.TiffFile(f["image"]) as t:
            s = t.series[0]
            out["image_size"] = [int(s.shape[s.axes.index("X")]), int(s.shape[s.axes.index("Y")])]
    return out


def _names(ids, names):
    """Gene names (the name when there is one, else the id) and the column of each gene after summing duplicates."""
    lab = np.array([n if n not in ("", "nan", "None") else i for i, n in zip(ids, names)], dtype=object)
    uniq, col = np.unique(lab, return_inverse=True)
    return [str(x) for x in uniq], col


def read_cellbin_gef(f, progress=print):
    """(cell ids, xy in DNBs, genes, cells x genes CSR) of a cell-bin GEF."""
    with h5py.File(f, "r") as h:
        off = np.array([float(h.attrs.get("offsetX", [0])[0]), float(h.attrs.get("offsetY", [0])[0])])
        cell = h["cellBin/cell"][:]
        gene = h["cellBin/gene"][:]
        exp = h["cellBin/cellExp"][:]
    s = lambda a: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in a], dtype=object)
    genes, col = _names(s(gene["geneID"]), s(gene["geneName"]) if "geneName" in gene.dtype.names else s(gene["geneID"]))
    indptr = np.r_[cell["offset"].astype(np.int64), cell["offset"][-1].astype(np.int64) + cell["geneCount"][-1]]
    if not np.array_equal(np.diff(indptr), cell["geneCount"].astype(np.int64)):
        raise ValueError(f"{f}: the cells' expression offsets are not in order")
    X = sp.csr_matrix((exp["count"].astype(np.float32), col[exp["geneID"]], indptr), shape=(len(cell), len(genes)))
    X.sum_duplicates()
    progress(f"{len(cell):,} cells, {len(genes):,} genes (cell-bin GEF)")
    ids = np.array([f"cell{int(i)}" for i in cell["id"]], dtype=object)
    return ids, np.c_[cell["x"], cell["y"]].astype(float) + off, genes, X


def read_gem(f, cells, progress=print):
    """(ids, xy in DNBs, genes, X) of a GEM: by CellID (cells=True; 0 = outside every cell, left out) or pooled into
    BIN x BIN DNB bins."""
    meta, n_head = _header(f)
    cols = gem_columns(f)
    low = {c.lower(): c for c in cols}
    cnt = next((low[k] for k in ("midcount", "midcounts", "umicount", "count") if k in low), None)
    if cnt is None:
        raise ValueError(f"{f}: no MIDCount column (columns: {', '.join(cols)})")
    cell_col = _cell_col(cols)
    use = [low["geneid"]] + ([low["genename"]] if "genename" in low else []) + [low["x"], low["y"], cnt] + \
        ([cell_col] if cells else [])
    gene_ix, gene_ids, gene_names = {}, [], []
    keys, X, sx, sy, nrow = {}, None, [], [], []
    rows = 0
    for ch in pd.read_csv(f, sep="\t", skiprows=n_head, usecols=use, chunksize=CHUNK, dtype={low["geneid"]: str}):
        gid = ch[low["geneid"]].astype(str).to_numpy()
        gnm = ch[low["genename"]].astype(str).to_numpy() if "genename" in low else gid
        cu, ci = np.unique(gid, return_index=True)
        for g, i in zip(cu, ci):
            if g not in gene_ix:
                gene_ix[g] = len(gene_ids); gene_ids.append(g); gene_names.append(gnm[i])
        gcol = pd.Index(cu).get_indexer(gid)
        gcol = np.array([gene_ix[g] for g in cu], np.int64)[gcol]
        x, y = ch[low["x"]].to_numpy(np.int64), ch[low["y"]].to_numpy(np.int64)
        v = ch[cnt].to_numpy(np.float32)
        if cells:
            k = pd.to_numeric(ch[cell_col], errors="coerce").fillna(0).to_numpy().astype(np.int64)
            ok = k > 0
            k, x, y, v, gcol = k[ok], x[ok], y[ok], v[ok], gcol[ok]
        else:
            k = (x // BIN) * (1 << 32) + (y // BIN)
        ku, kinv = np.unique(k, return_inverse=True)
        for u in ku[~pd.Index(ku).isin(list(keys))] if keys else ku:
            keys[u] = len(keys)
        kid = np.array([keys[u] for u in ku], np.int64)[kinv]
        part = sp.csr_matrix((v, (kid, gcol)), shape=(len(keys), len(gene_ids)))  # summed as it is read
        if X is None:
            X = part
        else:
            X.resize(part.shape)
            X = X + part
        sx.append(np.bincount(kid, x.astype(float))); sy.append(np.bincount(kid, y.astype(float)))
        nrow.append(np.bincount(kid))
        rows += len(ch)
        progress(f"  {rows:,} rows read")
    n, G = len(keys), len(gene_ids)
    X = sp.csr_matrix((n, G), dtype=np.float32) if X is None else X
    X.resize((n, G))
    pad = lambda a: np.pad(a, (0, n - len(a)))
    tx, ty, tn = (sum(pad(a) for a in z) for z in (sx, sy, nrow))
    order = np.empty(n, dtype=object)
    for u, i in keys.items():
        order[i] = u
    if cells:
        xy = np.c_[tx / tn, ty / tn]
        ids = np.array([f"cell{int(u)}" for u in order], dtype=object)
    else:  # the bin's centre
        bx, by = np.array([int(u) >> 32 for u in order]), np.array([int(u) & 0xFFFFFFFF for u in order])
        xy = np.c_[bx * BIN + BIN / 2, by * BIN + BIN / 2].astype(float)
        ids = np.array([f"b{BIN_UM}um_{a}_{b}" for a, b in zip(bx, by)], dtype=object)
    genes, col = _names(gene_ids, gene_names)
    X = sp.csr_matrix(X @ sp.csr_matrix((np.ones(G, np.float32), (np.arange(G), col)), shape=(G, len(genes))))
    progress(f"{n:,} {'cells' if cells else f'bins of {BIN_UM} um'}, {len(genes):,} genes (GEM, {rows:,} rows)")
    return ids, xy + _offset(meta), genes, X
