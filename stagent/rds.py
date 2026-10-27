"""Reading R data files (.rds, .rda / .RData) without R, and turning a Seurat or SpatialExperiment object into the
cells, genes, counts, coordinates and cell metadata that importer.py imports.

parse(): the R serialization format (XDR binary, versions 2 and 3, gzip / bzip2 / xz compressed), as written by
  saveRDS() and save(), read as a stream. Vectors become numpy arrays (NA: None in strings, NaN in doubles, INT_NA in
  integers and logicals), lists Python lists, S4 objects and attributes RObj, environments REnv. The values of slots
  and list entries named in SKIP (normalized data, scaled data, graphs, reductions, ...) are read past, not kept.
table(): from a Seurat object (Assay v3/v4 "counts" slot, or Assay5 "counts" layers, split per sample or not; cell
  coordinates from an FOV's centroids, a SlideSeq / Visium image, or meta.data columns) or a SpatialExperiment /
  SingleCellExperiment (assay "counts"; spatialCoords, or colData columns): {cell_id, genes, X (cells x genes CSR),
  obs (DataFrame), xy (or None)}.
"""
import bz2, gzip, io, lzma, struct
import numpy as np, pandas as pd, scipy.sparse as sp

INT_NA = -2 ** 31
NILVALUE, GLOBALENV, UNBOUNDVALUE, MISSINGARG, BASENAMESPACE = 254, 253, 252, 251, 250
NAMESPACE, PACKAGE, PERSIST, EMPTYENV, BASEENV, ATTRLANG, ATTRLIST, ALTREP, REF = 249, 248, 247, 242, 241, 240, 239, 238, 255
BCREPDEF, BCREPREF = 244, 243
NIL, SYM, LIST, CLO, ENV, PROM, LANG, SPECIAL, BUILTIN, CHAR = 0, 1, 2, 3, 4, 5, 6, 7, 8, 9
LGL, INT, REAL, CPLX, STR, DOT, ANY, VEC, EXPR, BCODE, EXTPTR, WEAKREF, RAW, S4 = 10, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25


class RObj:
    """An R value with attributes: kind (the SEXP type name), value, attrs {name: value}."""
    def __init__(self, kind, value=None, attrs=None):
        self.kind, self.value, self.attrs = kind, value, attrs or {}

    def cls(self):
        c = self.attrs.get("class")
        return [] if c is None else [str(x) for x in (c.value if isinstance(c, RObj) else c)]

    def __repr__(self):
        return f"RObj({self.kind}, class={self.cls()}, attrs={list(self.attrs)})"


class REnv:
    def __init__(self):
        self.frame = {}


class Sym(str):
    pass


SKIP = {"scale.data", "graphs", "neighbors", "reductions", "commands", "tools", "misc", "molecules", "segmentation"}
# "data" is skipped only after a "counts" slot (Seurat Assay); a SummarizedExperiment keeps its assays in "data"


def _open(path):
    f = open(path, "rb")
    magic = f.read(6)
    f.seek(0)
    if magic[:2] == b"\x1f\x8b":
        return io.BufferedReader(gzip.GzipFile(fileobj=f), 1 << 20)
    if magic[:3] == b"BZh":
        return io.BufferedReader(bz2.BZ2File(f), 1 << 20)
    if magic == b"\xfd7zXZ\x00":
        return io.BufferedReader(lzma.LZMAFile(f), 1 << 20)
    return io.BufferedReader(f, 1 << 20)


class _Reader:
    def __init__(self, f):
        self.f, self.refs, self.pending = f, [], None

    def bytes(self, n):
        v = self.f.read(n)
        if len(v) < n:
            raise ValueError("the R file ends early (truncated or not an R data file)")
        return v

    def int(self):
        return struct.unpack(">i", self.bytes(4))[0]

    def length(self):
        n = self.int()
        if n == -1:
            hi, lo = self.int(), self.int()
            n = (hi << 32) + (lo & 0xFFFFFFFF)
        return n

    def arr(self, dtype, n, skip=False):
        size = np.dtype(dtype).itemsize * n
        if skip:  # read past in chunks
            while size > 0:
                size -= len(self.bytes(min(size, 1 << 24)))
            return None
        return np.frombuffer(self.bytes(size), dtype=dtype)

    def header(self):
        fmt = self.bytes(2)
        if fmt != b"X\n":
            raise ValueError("only the binary (XDR) R serialization is supported; save with saveRDS(..., ascii = FALSE)")
        version = self.int()
        self.int(), self.int()
        if version == 3:
            self.bytes(self.int())
        elif version != 2:
            raise ValueError(f"unknown R serialization version {version}")

    def strvec(self):
        self.int()
        return [self.item() for _ in range(self.int())]

    def attrs(self, flags, skip=False):
        return _pairs(self.item(skip)) if flags & (1 << 9) else {}

    def item(self, skip=False):
        flags, self.pending = (self.pending, None) if self.pending is not None else (self.int(), None)
        t, has_attr, has_tag = flags & 0xFF, flags & (1 << 9), flags & (1 << 10)
        if t == NILVALUE or t == NIL:
            return None
        if t == REF:
            idx = flags >> 8
            return self.refs[(idx or self.int()) - 1]
        if t in (GLOBALENV, EMPTYENV, BASEENV, BASENAMESPACE, UNBOUNDVALUE, MISSINGARG):
            return RObj("special", t)
        if t in (NAMESPACE, PACKAGE, PERSIST):
            v = RObj("namespace", self.strvec())
            self.refs.append(v)
            return v
        if t == SYM:
            v = Sym(self.item() or "")
            self.refs.append(v)
            return v
        if t == ENV:
            self.int()  # locked
            env = REnv()
            self.refs.append(env)
            self.item()  # enclosure
            frame, hashtab = self.item(), self.item()
            self.item()  # attributes
            for part in [frame] + (list(hashtab.value) if isinstance(hashtab, RObj) and hashtab.kind == "list" else []):
                env.frame.update(_pairs(part))
            return env
        if t in (LIST, LANG, CLO, PROM, DOT, ATTRLANG, ATTRLIST):
            items = []  # a pairlist, read iteratively along its cdr
            attrs = {}
            while True:
                a = self.item(skip) if (has_attr or t in (ATTRLANG, ATTRLIST)) else None
                if a is not None and not attrs:
                    attrs = _pairs(a)
                tag = self.item() if has_tag else None
                drop = isinstance(tag, str) and (tag in SKIP or (tag == "data" and any(k == "counts" for k, _ in items)))
                items.append((tag, self.item(skip or drop)))
                nxt = self.int()
                nt = nxt & 0xFF
                if nt in (LIST, LANG, DOT, ATTRLANG, ATTRLIST) and t in (LIST, LANG, DOT, ATTRLANG, ATTRLIST):
                    t, has_attr, has_tag = nt, nxt & (1 << 9), nxt & (1 << 10)
                    continue
                self.pending = nxt
                last = self.item(skip)
                if last is not None:
                    items.append((None, last))
                break
            return RObj("pairlist", items, attrs)
        if t == EXTPTR:
            v = RObj("extptr")
            self.refs.append(v)
            self.item(), self.item()
            v.attrs = self.attrs(flags)
            return v
        if t == WEAKREF:
            v = RObj("weakref")
            self.refs.append(v)
            v.attrs = self.attrs(flags)
            return v
        if t in (SPECIAL, BUILTIN):
            return RObj("builtin", self.bytes(self.int()).decode())
        if t == CHAR:
            n = self.int()
            if n == -1:
                return None
            raw = self.bytes(n)
            try:
                return raw.decode("utf-8")
            except UnicodeDecodeError:
                return raw.decode("latin-1")
        if t in (LGL, INT):
            v = self.arr(">i4", self.length(), skip)
            return RObj("logical" if t == LGL else "integer", v, self.attrs(flags, skip))
        if t == REAL:
            return RObj("double", self.arr(">f8", self.length(), skip), self.attrs(flags, skip))
        if t == CPLX:
            v = self.arr(">f8", 2 * self.length(), skip)
            return RObj("complex", None if v is None else v[0::2] + 1j * v[1::2], self.attrs(flags, skip))
        if t == STR:
            n = self.length()
            return RObj("character", np.array([self.item() for _ in range(n)], dtype=object), self.attrs(flags, skip))
        if t in (VEC, EXPR):
            n = self.length()
            return RObj("list", [self.item(skip) for _ in range(n)], self.attrs(flags, skip))
        if t == RAW:
            v = self.arr("u1", self.length(), skip)
            return RObj("raw", None if v is None else v.tobytes(), self.attrs(flags, skip))
        if t == S4:
            return RObj("S4", None, self.attrs(flags, skip))
        if t == BCODE:
            reps = [None] * self.int()
            v = self._bc(reps)
            v.attrs = self.attrs(flags)
            return v
        if t == ALTREP:
            info, state, attr = self.item(), self.item(), self.item()
            v = _altrep(info, state)
            v.attrs = _pairs(attr)
            return v
        raise ValueError(f"unsupported R object type {t}")

    def _bc(self, reps):
        self.item()  # code
        for _ in range(self.int()):
            t = self.int()
            if t == BCODE:
                self._bc(reps)
            elif t in (LANG, LIST, BCREPDEF, BCREPREF, ATTRLANG, ATTRLIST):
                self._bclang(t, reps)
            else:
                self.item()
        return RObj("bytecode")

    def _bclang(self, t, reps):
        if t == BCREPREF:
            return reps[self.int()]
        if t in (BCREPDEF, LANG, LIST, ATTRLANG, ATTRLIST):
            pos = self.int() if t == BCREPDEF else -1
            if t == BCREPDEF:
                t = self.int()
            v = RObj("lang")
            if pos >= 0:
                reps[pos] = v
            if t in (ATTRLANG, ATTRLIST):
                self.item()
            self.item()  # tag
            self._bclang(self.int(), reps)
            self._bclang(self.int(), reps)
            return v
        return self.item()


def _altrep(info, state):
    """The few ALTREP classes R writes: compact integer / real sequences, deferred strings, wrappers."""
    name = str(info.value[0][1]) if isinstance(info, RObj) and info.kind == "pairlist" else ""
    if name == "compact_intseq":
        n, start, step = state.value[:3]
        return RObj("integer", (start + step * np.arange(int(n))).astype(np.int32))
    if name == "compact_realseq":
        n, start, step = state.value[:3]
        return RObj("double", start + step * np.arange(int(n)))
    if name == "deferred_string":
        arg = state.value[0][1] if state.kind == "pairlist" else state
        v = arg.value
        s = [None if (isinstance(x, (int, np.integer)) and x == INT_NA) or (isinstance(x, float) and np.isnan(x))
             else (str(int(x)) if float(x).is_integer() else repr(float(x))) for x in v]
        return RObj("character", np.array(s, dtype=object))
    if name.startswith("wrap_"):
        x = state.value[0] if state.kind == "list" else state
        return RObj(x.kind, x.value, x.attrs)
    raise ValueError(f"unsupported ALTREP class {name!r}")


def _pairs(p):
    if p is None:
        return {}
    if isinstance(p, RObj) and p.kind == "pairlist":
        return {str(k): v for k, v in p.value if k is not None}
    return {}


def parse(path):
    """The object of an .rds file, or {name: object} of an .rda / .RData file."""
    with _open(path) as f:
        rda = f.peek(5)[:5] in (b"RDX2\n", b"RDX3\n")
        if rda:
            f.read(5)
        r = _Reader(f)
        r.header()
        v = r.item()
    return _pairs(v) if rda else v


# ---------------- R values to Python ----------------
def _slot(o, name):
    return o.attrs.get(name) if isinstance(o, RObj) else None


def _names(o):
    n = _slot(o, "names")
    return [] if n is None else [str(x) for x in n.value]


def _named(o):
    """A named R list as {name: value}."""
    if not isinstance(o, RObj) or o.kind != "list":
        return {}
    return dict(zip(_names(o), o.value))


def _vec(o):
    """An atomic R vector as numpy (factors as their labels, NA as None / NaN)."""
    if o is None:
        return None
    v = o.value
    if isinstance(v, np.ndarray) and v.dtype.byteorder == ">":
        v = v.astype(v.dtype.newbyteorder("="))
    if o.kind == "integer" and "levels" in o.attrs:
        lev = np.asarray(o.attrs["levels"].value, dtype=object)
        out = np.full(len(v), None, dtype=object)
        ok = v != INT_NA
        out[ok] = lev[v[ok] - 1]
        return out
    if o.kind == "logical":
        out = v.astype(object)
        out[v == INT_NA] = None
        out[v == 1], out[v == 0] = True, False
        return out
    if o.kind == "integer":
        return np.where(v == INT_NA, np.nan, v) if (v == INT_NA).any() else v
    return v


def _dataframe(o):
    """A data.frame (or Bioconductor DataFrame) as pandas."""
    if "DFrame" in o.cls() or "DataFrame" in o.cls():
        cols = _named(_slot(o, "listData"))
        rn = _slot(o, "rownames")
        index = None if rn is None else [str(x) for x in rn.value]
    else:
        cols = dict(zip(_names(o), o.value))
        rn = o.attrs.get("row.names")
        index = [str(x) for x in rn.value] if rn is not None and rn.kind == "character" else None
    data = {}
    for k, v in cols.items():
        if isinstance(v, RObj) and v.kind in ("integer", "double", "logical", "character") and "dim" not in v.attrs:
            data[k] = _vec(v)
    return pd.DataFrame(data, index=index)


def _matrix(o):
    """A count matrix (dgCMatrix / dgTMatrix / dense matrix) as (scipy sparse features x cells, rownames, colnames)."""
    cls = o.cls()
    dimnames = _slot(o, "Dimnames") if any(c.startswith("dg") for c in cls) else o.attrs.get("dimnames")
    dn = [None, None]
    if isinstance(dimnames, RObj) and dimnames.kind == "list":
        dn = [None if d is None else [str(x) for x in d.value] for d in dimnames.value]
    if "dgCMatrix" in cls:
        nr, nc = (int(x) for x in _slot(o, "Dim").value)
        m = sp.csc_matrix((_slot(o, "x").value.astype(np.float32), _slot(o, "i").value.astype(np.int32),
                           _slot(o, "p").value.astype(np.int64)), shape=(nr, nc))
    elif "dgTMatrix" in cls:
        nr, nc = (int(x) for x in _slot(o, "Dim").value)
        m = sp.coo_matrix((_slot(o, "x").value.astype(np.float32), (_slot(o, "i").value.astype(np.int64),
                           _slot(o, "j").value.astype(np.int64))), shape=(nr, nc)).tocsc()
    elif o.kind in ("integer", "double") and "dim" in o.attrs:
        nr, nc = (int(x) for x in o.attrs["dim"].value)
        m = sp.csc_matrix(np.asarray(o.value, np.float32).reshape(nc, nr).T)
    else:
        raise ValueError(f"unsupported count matrix class {cls or o.kind} (save it as a dgCMatrix)")
    return m, dn[0], dn[1]


def _seurat(o):
    assays = _named(_slot(o, "assays"))
    active = _slot(o, "active.assay")
    name = str(active.value[0]) if active is not None else next(iter(assays))
    a = assays[name]
    meta = _dataframe(_slot(o, "meta.data"))
    cells = list(meta.index)
    if "Assay5" in a.cls():
        layers = _named(_slot(a, "layers"))
        feats = _slot(a, "features")
        cmap = _slot(a, "cells")
        fnames = [str(x) for x in feats.attrs["dimnames"].value[0].value]
        cnames = [str(x) for x in cmap.attrs["dimnames"].value[0].value]
        flayers = [str(x) for x in feats.attrs["dimnames"].value[1].value]
        clayers = [str(x) for x in cmap.attrs["dimnames"].value[1].value]
        fm = np.asarray(feats.value).reshape(len(flayers), len(fnames)).T.astype(bool)
        cm = np.asarray(cmap.value).reshape(len(clayers), len(cnames)).T.astype(bool)
        keys = [k for k in layers if k == "counts" or k.startswith("counts.")]
        if not keys:
            raise ValueError(f"assay {name!r} has no counts layer (layers: {', '.join(layers)})")
        cols, rows = [], []
        for k in keys:  # one layer per sample after split(): cells of each layer, possibly a subset of the features
            m, _, _ = _matrix(layers[k])
            f = np.where(fm[:, flayers.index(k)])[0]
            if len(f) != len(fnames):
                m = sp.csr_matrix((np.ones(len(f), np.float32), (f, np.arange(len(f)))), shape=(len(fnames), len(f))) @ m
            cols.append(sp.csr_matrix(m.T))
            rows.extend(np.where(cm[:, clayers.index(k)])[0])
        X = sp.vstack(cols).tocsr()
        order = np.argsort(rows)
        X, cell_order = X[order], [cnames[i] for i in np.asarray(rows)[order]]
        genes = fnames
    else:
        m, genes, cell_order = _matrix(_slot(a, "counts"))
        X = m.T.tocsr()
    if cell_order is not None and list(cell_order) != cells:
        idx = {c: i for i, c in enumerate(cell_order)}
        X = X[[idx[c] for c in cells]]
    return cells, genes, X, meta, _seurat_xy(o, cells, meta)


def _seurat_xy(o, cells, meta):
    """Cell coordinates of a Seurat object: the centroids of its FOV images, SlideSeq coordinates, or None."""
    xy = np.full((len(cells), 2), np.nan)
    pos = {c: i for i, c in enumerate(cells)}
    found = False
    for img in _named(_slot(o, "images")).values():
        cls = img.cls()
        if "FOV" in cls:
            cen = _named(_slot(img, "boundaries")).get("centroids")
            if cen is None:
                continue
            co = _slot(cen, "coords")
            n = int(co.attrs["dim"].value[0])
            c = np.asarray(co.value, float).reshape(2, n).T
            ids = [str(x) for x in _slot(cen, "cells").value]
        elif "SlideSeq" in cls or "VisiumV1" in cls:
            df = _dataframe(_slot(img, "coordinates"))
            xc, yc = ("x", "y") if "x" in df else ("imagecol", "imagerow")
            c, ids = df[[xc, yc]].to_numpy(float), list(df.index)
        else:
            continue
        for k, cid in enumerate(ids):
            if cid in pos:
                xy[pos[cid]] = c[k]
                found = True
    return xy if found and not np.isnan(xy).any() else None


def _sce(o):
    assays = _slot(o, "assays")
    data = _named(_slot(_slot(assays, "data"), "listData"))
    if "counts" not in data:
        raise ValueError(f"no 'counts' assay (assays: {', '.join(data)})")
    m, genes, cells = _matrix(data["counts"])
    names = _slot(o, "NAMES")
    if genes is None and names is not None:
        genes = [str(x) for x in names.value]
    meta = _dataframe(_slot(o, "colData"))
    cells = cells or list(meta.index)
    meta.index = cells
    xy = None
    sc = _named(_slot(_slot(o, "int_colData"), "listData")).get("spatialCoords")
    if sc is not None and "dim" in sc.attrs:
        n, k = (int(x) for x in sc.attrs["dim"].value)
        if n and k >= 2:
            xy = np.asarray(sc.value, float).reshape(k, n).T[:, :2]
    return cells, genes, m.T.tocsr(), meta, xy


def table(path):
    """Cells of a Seurat / SpatialExperiment / SingleCellExperiment object in an .rds (or the first such object of an
    .rda): {cell_id, genes, X, obs, xy}."""
    obj = parse(path)
    cands = list(obj.values()) if isinstance(obj, dict) else [obj]
    for o in cands:
        if not isinstance(o, RObj):
            continue
        cls = o.cls()
        if "Seurat" in cls:
            cells, genes, X, meta, xy = _seurat(o)
        elif any(c in cls for c in ("SpatialExperiment", "SingleCellExperiment", "SummarizedExperiment")):
            cells, genes, X, meta, xy = _sce(o)
        else:
            continue
        if genes is None:
            raise ValueError("the count matrix has no gene names")
        return {"cell_id": np.array(cells, dtype=object), "genes": list(genes), "X": X, "obs": meta, "xy": xy,
                "kind": "Seurat" if "Seurat" in cls else cls[0]}
    kinds = sorted({c for o in cands if isinstance(o, RObj) for c in o.cls()})
    raise ValueError("no Seurat, SpatialExperiment or SingleCellExperiment object in this file"
                     + (f" (found: {', '.join(kinds)})" if kinds else ""))
