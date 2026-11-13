"""Marker-gene based hierarchical annotation; needs no labeled reference.

Marker list format (plain text; indentation defines the hierarchy, any depth):
    # comment
    Epithelial: EPCAM+ KRT8+
      Absorptive: CA1|SLC26A3|AQP8+ SPINK4-
        BEST4+ colonocytes: BEST4+ OTOP2+
      Goblet cells: SPINK4+ RETNLB+
    T and NK cells: CD3E|NKG7+ EPCAM-
  - one type per line, "Name: markers"; children are indented under their parent (2 spaces or a tab per level)
  - GENE+ expressed, GENE- not expressed, no sign = +
  - A|B|C+ means any one of these genes (the strongest counts); markers on a line are averaged
  - "Parent > Child: markers" is accepted too
  - cell types are the leaves (most specific entries); a type with subtypes is a lineage / group, so its subtypes
    must cover all of its cells (add a catch-all subtype such as "Effector/memory CD4 T cells" if needed)
Annotation, top-down: at each level the parent's cells are re-clustered (PCA, kNN graph, Leiden; with resolution
"auto" several resolutions are tried and the one that lets the markers assign the most cells wins); every cluster gets a
marker score per child type (detection rate z-scored across the clusters of that parent, sign-adjusted,
averaged over the child's markers); a cluster goes to the best child only if the score and its margin over the
runner-up are high enough, otherwise its cells stay at the parent level. Children are then annotated the same
way inside their own cells. Cells no top-level type claims are "Unassigned".
A child type no cluster went to is then looked for (RARE): in finer clusterings of the same graph, with smaller
clusters, accepted only where the cells really show its markers; and, for cells too few or too scattered to form a
cluster, in the single cells that show its markers clearly. A type still without cells is listed in "not_found" with
the detection of its markers, and the expression model then gives it no cells (session._only_found).
"""
import difflib, re
import numpy as np, scipy.sparse as sp

UNASSIGNED = "Unassigned"
MIN_CELLS, N_PCS, N_NEIGHBORS, MIN_CLUSTER = 60, 30, 15, 30
RESOLUTION = "auto"  # default; a number fixes the Leiden resolution for every level
MIN_SCORE, MIN_MARGIN = 0.5, 0.25
TOKEN = re.compile(r"^([A-Za-z0-9._|-]+?)([+-]?)$")


# ---------------- parsing ----------------
CSV_COLUMNS = ("cell_type", "parent", "positive_markers", "negative_markers")


def is_csv(text):
    first = next((l for l in text.splitlines() if l.strip() and not l.lstrip().startswith("#")), "")
    return first.replace(" ", "").lower().startswith("cell_type,")


def csv_to_text(text):
    """Marker list CSV (cell_type, parent, positive_markers, negative_markers; one row per type, parent empty for
    lineages; markers separated by spaces or semicolons, A|B = any of A or B) -> the indented text format."""
    import csv, io
    rows = [r for r in csv.DictReader(io.StringIO("\n".join(l for l in text.splitlines() if not l.lstrip().startswith("#"))))]
    rows = [{(k or "").strip().lower(): (v or "").strip() for k, v in r.items()} for r in rows]
    if not rows or "cell_type" not in rows[0]:
        raise ValueError("a marker list CSV needs the columns " + ", ".join(CSV_COLUMNS))
    kids = {}
    for r in rows:
        if r["cell_type"]:
            kids.setdefault(r.get("parent", ""), []).append(r)
    out, seen = [], set()

    def toks(s, sign):
        return [f"{g}{sign}" for g in re.split(r"[\s;,]+", s) if g]

    def walk(parent, depth):
        for r in kids.get(parent, []):
            if r["cell_type"] in seen:
                continue
            seen.add(r["cell_type"])
            out.append("  " * depth + f"{r['cell_type']}: " + " ".join(toks(r.get("positive_markers", ""), "+")
                                                                       + toks(r.get("negative_markers", ""), "-")))
            walk(r["cell_type"], depth + 1)
    walk("", 0)
    missing = [r["cell_type"] for r in rows if r["cell_type"] and r["cell_type"] not in seen]
    if missing:  # parents that are not in the list: keep the rows so the check reports them
        out += [f"{kids_parent} > {r['cell_type']}: " + " ".join(toks(r.get("positive_markers", ""), "+") + toks(r.get("negative_markers", ""), "-"))
                for kids_parent, rs in kids.items() for r in rs if r["cell_type"] in missing]
    return "\n".join(out) + "\n"


def text_to_csv(text):
    """The indented text format -> marker list CSV (see csv_to_text)."""
    import csv, io
    p = parse(text, None)
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(CSV_COLUMNS)

    def walk(nodes, parent):
        for n in nodes:
            pos = " ".join("|".join(m["genes"]) for m in n["markers"] if m["sign"] > 0)
            neg = " ".join("|".join(m["genes"]) for m in n["markers"] if m["sign"] < 0)
            w.writerow([n["name"], parent or "", pos, neg])
            walk(n["children"], n["name"])
    walk(p["tree"], None)
    return buf.getvalue()


def parse(text, panel, removed=()):
    """Return {"ok", "errors", "warnings", "tree", "types", "parent", "top", "genes", "canonical"}.
    text: the indented text format or a marker list CSV. panel: the sample's genes (None: genes are not checked).
    removed: genes of the panel that the sample's filter deleted; markers naming them are left out with a warning."""
    if is_csv(text):
        try:
            text = csv_to_text(text)
        except ValueError as e:
            return {"ok": False, "errors": [str(e)], "warnings": [], "tree": [], "types": [UNASSIGNED],
                    "parent": {}, "top": {}, "genes": [], "canonical": "", "leaves": [UNASSIGNED], "n_types": 0,
                    "depth": 0}
    upper = {g.upper(): g for g in panel} if panel is not None else None
    removed = {g.upper() for g in removed or ()}
    errors, warnings, roots, stack, names = [], [], [], [], {}

    def gene(name, ln):
        if upper is None:
            return name
        g = upper.get(name.upper())
        if g is None and name.upper() in removed:
            warnings.append(f"line {ln}: gene {name} was removed by the sample's gene filter; it is left out")
            return False
        if g is None:
            sug = difflib.get_close_matches(name.upper(), list(upper), n=3, cutoff=0.6)
            errors.append(f"line {ln}: gene {name} is not in the panel" + (f"; did you mean {', '.join(upper[s] for s in sug)}?" if sug else ""))
        return g

    for ln, raw in enumerate(text.splitlines(), 1):
        line = raw.split("#", 1)[0].rstrip()
        if not line.strip():
            continue
        indent = len(raw.expandtabs(2)) - len(raw.expandtabs(2).lstrip())
        if ":" not in line:
            errors.append(f"line {ln}: expected 'Name: GENE+ GENE- ...'")
            continue
        name, mk = line.strip().split(":", 1)
        path = [p.strip() for p in name.split(">")]
        name = path[-1]
        markers = []
        for tok in re.split(r"[\s,;]+", mk.strip()):
            if not tok:
                continue
            m = TOKEN.match(tok)
            if not m:
                errors.append(f"line {ln}: cannot read marker '{tok}'")
                continue
            genes = [gene(g, ln) for g in m.group(1).split("|") if g]
            if None not in genes and any(genes):
                markers.append({"genes": [g for g in genes if g], "sign": -1 if m.group(2) == "-" else 1})
        if not name:
            errors.append(f"line {ln}: missing type name"); continue
        if name in names or name == UNASSIGNED:
            errors.append(f"line {ln}: type '{name}' is defined twice (names must be unique)"); continue
        if not markers:
            errors.append(f"line {ln}: '{name}' has no markers")
        elif all(m["sign"] < 0 for m in markers):
            warnings.append(f"line {ln}: '{name}' has only negative markers")
        node = {"name": name, "markers": markers, "children": [], "line": ln}
        if len(path) > 1:  # path syntax: parents must already exist
            parent = names.get(path[-2])
            if parent is None:
                errors.append(f"line {ln}: parent '{path[-2]}' is not defined above"); continue
            parent["children"].append(node)
        else:
            while stack and stack[-1][0] >= indent:
                stack.pop()
            (stack[-1][1]["children"] if stack else roots).append(node)
            stack.append((indent, node))
        names[name] = node
    types, parent, top = [UNASSIGNED], {UNASSIGNED: None}, {UNASSIGNED: UNASSIGNED}

    def walk(nodes, par, root):
        for n in nodes:
            types.append(n["name"]); parent[n["name"]] = par; top[n["name"]] = root or n["name"]
            walk(n["children"], n["name"], root or n["name"])
    walk(roots, None, None)
    if not roots and not errors:
        errors.append("the marker list is empty")
    for n in names.values():
        if len(n["children"]) == 1:
            warnings.append(f"'{n['name']}' has a single subtype '{n['children'][0]['name']}': every {n['name']} cell "
                            f"becomes {n['children'][0]['name']}; add another subtype if some are something else")
    genes = sorted({g for n in names.values() for m in n["markers"] for g in m["genes"]})
    return {"ok": not errors, "errors": errors, "warnings": warnings, "tree": roots, "types": types,
            "parent": parent, "top": top, "genes": genes, "canonical": to_text(roots),
            "leaves": [UNASSIGNED] + [n["name"] for n in names.values() if not n["children"]],
            "n_types": len(types) - 1, "depth": _depth(roots)}


def _depth(nodes):
    return 0 if not nodes else 1 + max(_depth(n["children"]) for n in nodes)


def to_text(nodes, level=0):
    out = []
    for n in nodes:
        mk = " ".join("|".join(m["genes"]) + ("-" if m["sign"] < 0 else "+") for m in n["markers"])
        out.append("  " * level + f"{n['name']}: {mk}")
        out += to_text(n["children"], level + 1).splitlines() if n["children"] else []
    return "\n".join(out)


def descendants(tree):
    """{node name: [node and all its descendants]}, used as type groups for priors."""
    out = {}

    def walk(n):
        names = [n["name"]]
        for c in n["children"]:
            names += walk(c)
        out[n["name"]] = names
        return names
    for n in tree:
        walk(n)
    return out


# ---------------- clustering ----------------
RESOLUTIONS = (0.5, 1.0, 2.0, 4.0)  # candidates tried per level when resolution is "auto"
AUTO = "auto"
# Rare-type search: a child type no cluster of its level went to is looked for in finer clusterings of the same
# neighbor graph (RARE_RESOLUTIONS, clusters down to RARE_MIN_CLUSTER cells). A fine cluster becomes that type when it
# scores best for it with the usual score and margin, and its cells really show the type's markers: at least one
# positive marker is detected in RARE_MIN_DET of its cells and RARE_FOLD times as often as in the level's cells (or,
# a common one, in RARE_COMMON_DET of them and RARE_COMMON_FOLD times as often), no
# positive marker is rarer there than in the level's cells (a list often holds weak or common markers too), and no
# negative marker is enriched there. Types still without cells are reported with the detection of their markers
# (res["not_found"]).
# A type found by no cluster even so (its cells too few or too scattered to form one: Paneth cells among stem cells,
# ICC between muscle cells) gets the single cells that show its markers clearly: with one specific positive marker (a
# gene detected in at most RARE_SPECIFIC of the level's cells), RARE_CELL_COUNTS of its transcripts; with several,
# at least RARE_PAIR_COUNTS of each of two of them, in RARE_FOLD times as many cells as the two would meet by chance;
# no negative marker detected; at least RARE_MIN_CELLS such cells.
RARE = True
RARE_RESOLUTIONS = (2.0, 4.0, 8.0)
RARE_MIN_CLUSTER = 10
RARE_MIN_DET, RARE_FOLD = 0.25, 3.0
RARE_COMMON_DET, RARE_COMMON_FOLD = 0.4, 1.5  # (a catch-all type's markers are common in its level: e.g. KLRB1 in CD4 T)
RARE_SPECIFIC, RARE_CELL_COUNTS, RARE_PAIR_COUNTS, RARE_MIN_CELLS = 0.05, 3, 2, 10


def _topk_torch(Z, k, dev):
    """Indices of the k nearest rows (the row itself included) with torch on `dev`, blockwise."""
    import torch
    Zt = torch.tensor(Z, device=dev)
    n = len(Z)
    chunk = int(np.clip(1.2e9 / n, 256, 8192))  # rows per distance block, about 5 GB at most
    nbr = [torch.topk(torch.cdist(Zt[i:i + chunk], Zt), k, largest=False)[1].cpu() for i in range(0, n, chunk)]
    del Zt
    return torch.cat(nbr).numpy().astype(np.int64)


MLX_BLOCK_BYTES = 512 * 2 ** 20  # distances of a block of rows to all rows on the Apple GPU (MLX): at most this much
_MLX = {}                        # whether MLX can use the Apple GPU here (checked once)
_MLX_LOCK = __import__("threading").Lock()  # one MLX search at a time (two jobs at once take turns)


def _mlx_gpu():
    """MLX can compute on the Apple GPU here: the macOS app on Apple silicon (MLX needs macOS 14 or later; an older
    macOS, or a failing import, leaves the search on the CPU). A product on the GPU checks that its kernels load and
    that it computes in full float32: MLX uses TF32 (about 10 bits) on the GPUs that have it unless MLX_ENABLE_TF32=0
    (set in stagent/__init__), which MLX reads once; if something multiplied on the GPU before it was set, or MLX
    changes its default, the search stays on the CPU."""
    if "ok" not in _MLX:
        import sys, platform
        ok = False
        if sys.platform == "darwin" and platform.machine() == "arm64":
            try:
                import mlx.core as mx
                ok = bool(mx.is_available(mx.gpu))
                if ok:
                    rng = np.random.default_rng(0)
                    a, b = rng.standard_normal((64, 32), np.float32), rng.standard_normal((32, 2048), np.float32)
                    with _MLX_LOCK, mx.stream(mx.gpu):
                        p = np.array(mx.array(a) @ mx.array(b))  # (8 rows or more: the matrix kernel the search uses)
                    err = float(np.abs(p - a.astype(np.float64) @ b).max() / (np.abs(a) @ np.abs(b)).max())
                    if err > 1e-5:  # float32 rounds by about 1e-7; TF32 by 1e-3
                        print(f"exact neighbors: the Apple GPU multiplies in reduced precision (relative error "
                              f"{err:.1e}; MLX_ENABLE_TF32 must be 0 before MLX's first product); on the CPU")
                        ok = False
            except Exception as e:
                print(f"exact neighbors: no Apple GPU through MLX ({type(e).__name__}: {e}); on the CPU")
                ok = False
        _MLX["ok"] = ok
    return _MLX["ok"]


def _topk_mlx(Z, k, dev=None):
    """Indices of the k nearest rows of every row by their float32 squared distance (any of equally near rows; the
    row itself is one of them unless k others are as near), with MLX on the Apple GPU (dev: an MLX device, the GPU by
    default). A block of rows at a time (MLX_BLOCK_BYTES of distances, at least 64 rows): |z_j|^2 - 2 z_i.z_j as one
    product. MLX sorts a whole row for a partition, which over all n columns would cost more than the distances, so
    the columns go in groups of g: the k groups with the smallest minimum hold k nearest rows (a group with a row
    nearer than the k-th has its minimum below the k-th distance, and at most k groups have), and only their k g
    columns are sorted. Rows with NaN or infinity are left to the CPU search (ValueError)."""
    import mlx.core as mx
    Z = np.asarray(Z, np.float32)
    n, dim = Z.shape
    if not np.isfinite(Z).all():
        raise ValueError("non-finite values")
    g = max(1, int(round(np.sqrt(n / k))))  # about as many groups per row (n / g) as columns sorted (k g)
    ng = -(-n // g)
    direct = g == 1 or ng <= k  # few rows: whole rows sorted
    if direct:
        g, ng = 1, n
    w = -(-(dim + 1) // 16) * 16  # (zeros up to a multiple of 16: the GPU's product runs on whole tiles)
    A = np.zeros((ng * g, w), np.float32)  # columns: [z_j, |z_j|^2]; the padding to whole groups is never near
    A[:n, :dim] = Z
    A[:n, dim] = (Z ** 2).sum(1)
    A[n:, dim] = 1e30
    L = np.zeros((n, w), np.float32)  # rows: [-2 z_i, 1]
    np.multiply(Z, -2, out=L[:, :dim])
    L[:, dim] = 1
    rows = int(min(n, max(64, MLX_BLOCK_BYTES // (4 * len(A)))))
    nbr = np.empty((n, k), np.int64)
    with _MLX_LOCK, mx.stream(dev if dev is not None else mx.gpu):
        try:
            Lm, At = mx.array(L), mx.array(A).T
            offs = mx.arange(g, dtype=mx.uint32)
            for i in range(0, n, rows):
                P = Lm[i:i + rows] @ At  # rows x columns: the squared distance less |z_i|^2
                if direct:
                    c = mx.argpartition(P, k - 1, axis=1)[:, :k]
                else:
                    grp = mx.argpartition(P.reshape(-1, ng, g).min(axis=2), k - 1, axis=1)[:, :k].astype(mx.uint32)
                    cols = (grp[:, :, None] * g + offs).reshape(-1, k * g)
                    v = mx.take_along_axis(P, cols, axis=1)
                    c = mx.take_along_axis(cols, mx.argpartition(v, k - 1, axis=1)[:, :k], axis=1)
                nbr[i:i + rows] = np.array(c)
            del Lm, At, P
        finally:
            try:
                mx.clear_cache()  # (the GPU memory MLX keeps for reuse)
            except Exception:
                pass
    if nbr.min() < 0 or nbr.max() >= n:  # (only with absurdly large values: their distances overflow)
        raise ValueError("neighbor out of range")
    return nbr


KNN_TILE = (64, 4096)  # rows x columns of a distance tile of the exact search on the CPU: 1 MB, it stays in the cache
KNN_BLOCK = 512        # rows a thread takes at once (their tiles share the columns read from memory)
KNN_SAMPLE = 32        # every KNN_SAMPLE-th row is searched first: it bounds each row's k-th distance from above


def cpu_threads():
    """Threads of the exact search on the CPU (numpy): KNN_THREADS, fewer on a computer with fewer cores."""
    import os
    return max(1, min(KNN_THREADS, os.cpu_count() or 1))


_BLAS = {"n": 0, "lim": None}       # searches under the one-thread limit now, and the limit to undo after the last
_BLAS_LOCK = __import__("threading").Lock()


@__import__("contextlib").contextmanager
def _one_blas_thread():
    """BLAS on one thread while the threads of the exact search each run their own products (OpenBLAS or MKL: a
    product this thin is several times slower spread over all cores, and several at once on all cores slower still).
    A BLAS threadpoolctl cannot set (Apple's Accelerate) keeps its own threads. The limit is process-wide: two
    searches at once (two jobs) share it, set by the first and lifted by the last (else the one ending last would put
    back the other's one thread for good)."""
    with _BLAS_LOCK:
        if _BLAS["n"] == 0:
            try:
                from threadpoolctl import threadpool_limits
                _BLAS["lim"] = threadpool_limits(1, user_api="blas")
            except Exception:
                _BLAS["lim"] = None
        _BLAS["n"] += 1
    try:
        yield
    finally:
        with _BLAS_LOCK:
            _BLAS["n"] -= 1
            if _BLAS["n"] == 0 and _BLAS["lim"] is not None:
                _BLAS["lim"].restore_original_limits()
                _BLAS["lim"] = None


def _blocks(n, rows, threads, fn):
    """fn(start, stop) for consecutive blocks of `rows` of range(n), on `threads` threads (numpy and BLAS let go of
    the GIL while they compute)."""
    from concurrent.futures import ThreadPoolExecutor
    import itertools, threading
    starts, lock = itertools.count(0, rows), threading.Lock()

    def work():
        while True:
            with lock:
                i = next(starts)
            if i >= n:
                return
            fn(i, min(i + rows, n))
    if threads == 1 or n <= rows:
        return work()
    with ThreadPoolExecutor(threads) as ex:
        for f in [ex.submit(work) for _ in range(threads)]:
            f.result()


class _Tiles:
    """Squared distances of rows of Z to all rows in float32, tile by tile: P[q, j] = |z_j|^2 - 2 z_q.z_j (the squared
    distance less |z_q|^2, the same along a row), one BLAS product per tile of KNN_TILE that stays in the cache; a
    block of KNN_BLOCK rows goes through the columns tile by tile, so each column tile is read from memory once per
    block. Columns are the rows of Z; with `step`, every step-th row comes first (self.order): the sample, self.m
    columns. A thread keeps its own buffers (a few MB)."""

    def __init__(self, Z, step=None):
        n, dim = Z.shape
        self.Z, self.n, self.m = Z, n, -(-n // step) if step else 0
        self.order = np.argsort(np.arange(n) % step, kind="stable") if step and step > 1 else np.arange(n)
        self.A = np.empty((n, dim + 1), np.float32)
        self.A[:, :dim] = Z[self.order]
        self.A[:, dim] = (self.A[:, :dim] ** 2).sum(1)
        self.local = __import__("threading").local()

    def _tile(self, L, c0, c1):
        P = self.local.P[:len(L), :c1 - c0]
        return np.matmul(L, self.A[c0:c1].T, out=P)

    def _below(self, P, t, c0):
        M = np.less_equal(P, t[:, None], out=self.local.M[:len(P), :P.shape[1]])
        idx = np.flatnonzero(M)
        r = idx // P.shape[1]
        return r, idx - r * P.shape[1] + c0, P.ravel()[idx]

    def scan(self, rows, t=None, k=None):
        """[(first row, (row, column, value))] of the entries at most t (one per row) of the rows `rows` (at most
        KNN_BLOCK), one per KNN_TILE[0] of them, sorted by row (columns ascending within a row). t=None: the k-th
        smallest of the sample (an upper bound for the k-th smallest of the row)."""
        h, w = KNN_TILE
        loc = self.local
        if getattr(loc, "P", None) is None:
            loc.P, loc.M = np.empty((h, max(w, self.m)), np.float32), np.empty((h, max(w, self.m)), bool)
        L = np.empty((len(rows), self.A.shape[1]), np.float32)
        np.multiply(self.Z[rows], -2, out=L[:, :-1])
        L[:, -1] = 1
        subs = range(0, len(rows), h)
        parts = {s: [] for s in subs}
        start = 0
        if t is None:
            t, start = np.empty(len(rows), np.float32), self.m
            for s in subs:
                S = self._tile(L[s:s + h], 0, self.m)
                t[s:s + h] = np.partition(S, k - 1, axis=1)[:, k - 1]
                parts[s].append(self._below(S, t[s:s + h], 0))
        for c0 in range(start, self.n, w):
            for s in subs:
                parts[s].append(self._below(self._tile(L[s:s + h], c0, min(c0 + w, self.n)), t[s:s + h], c0))
        out = []
        for s in subs:
            r, c, v = (np.concatenate(x) for x in zip(*parts[s]))
            o = np.argsort(r.astype(np.int16), kind="stable")  # (a radix sort)
            out.append((s, (r[o], c[o], v[o])))
        return out


def _topk_numpy(Z, k, threads=None):
    """Indices of the k nearest rows of every row by their float32 squared distance (any of equally near rows; the row
    itself is one of them unless k others are as near), through BLAS (a KD-tree is slow in 30 dimensions) on
    cpu_threads() threads, in tiles that stay in the cache (a few MB per thread, whatever the number of rows). Rows
    are first compared with every KNN_SAMPLE-th row: their k-th nearest there bounds the k-th distance; the other
    tiles only keep the few entries under that bound."""
    Z = np.asarray(Z, np.float32)
    n = len(Z)
    T = _Tiles(Z, max(1, min(KNN_SAMPLE, n // (4 * k))))  # the sample has at least 4 k rows (all when n is small)
    nbr = np.empty((n, k), np.int64)

    def block(i, j):
        for s, (r, c, v) in T.scan(np.arange(i, j), k=k):
            i0 = i + s
            b = min(KNN_TILE[0], j - i0)
            cnt = np.bincount(r, minlength=b)
            if cnt.min() >= k and cnt.max() <= 8192:
                pos = np.arange(len(r)) - np.repeat(np.cumsum(cnt) - cnt, cnt)
                V, C = np.full((b, int(cnt.max())), np.inf, np.float32), np.zeros((b, int(cnt.max())), np.int64)
                V[r, pos], C[r, pos] = v, c
                nbr[i0:i0 + b] = T.order[np.take_along_axis(C, np.argpartition(V, k - 1, axis=1)[:, :k], 1)]
                continue
            ends = np.cumsum(cnt)
            for q in range(b):  # many exact duplicates (or NaN): one row at a time
                cq, vq = c[ends[q] - cnt[q]:ends[q]], v[ends[q] - cnt[q]:ends[q]]
                if len(cq) < k:  # NaN: they go last
                    cq, vq = np.arange(n), T.A @ (-2 * np.append(Z[i0 + q], -0.5))
                nbr[i0 + q] = T.order[cq[np.argsort(vq, kind="stable")[:k]]]
    with _one_blas_thread():
        _blocks(n, KNN_BLOCK, threads or cpu_threads(), block)
    return nbr


def _within(Z, rows, bound, fn, threads=None):
    """fn(i, cols) for each row i of `rows`, cols: every row j of Z (ascending) with |z_j|^2 - 2 z_i.z_j <= bound[i]
    in float32, the squared distance less |z_i|^2 (tiles as in _topk_numpy)."""
    Z = np.asarray(Z, np.float32)
    T = _Tiles(Z)

    def block(i, j):
        sel = rows[i:j]
        t = np.nextafter(np.asarray(bound[sel], np.float32), np.float32(np.inf))  # rounded up
        for s, (r, c, _) in T.scan(sel, t):
            ends = np.cumsum(np.bincount(r, minlength=min(KNN_TILE[0], len(sel) - s)))
            for q in range(len(ends)):
                fn(sel[s + q], c[ends[q - 1] if q else 0:ends[q]])
    with _one_blas_thread():  # (KNN_TILE[0] rows at once: a row can have many columns under its bound)
        _blocks(len(rows), KNN_TILE[0], threads or cpu_threads(), block)


MAX_PCA_GENES, N_HVG = 6000, 3000  # above MAX_PCA_GENES genes (whole-transcriptome data such as Visium HD) the PCA
                                    # uses the N_HVG most variable ones; a panel (Xenium 5K too) uses all its genes
PCA_FIT_MAX = None  # cells an embedding's PCA is fitted on when the sample has more; None fits on every cell.
PCA_FIT_SEED = 0    # Measured (2026-11-12, 400,000 cells x 6,000 genes at the 100,000-of-5,000,000 ratio): sampling
                    # saves no peak memory once GRAM_SLAB is small -- the Gram pass is no longer what the peak is --
                    # and PC5/PC6 of the sample point elsewhere than those of a fit on all cells (|corr| 0.16, 0.12).
                    # It does cut the Gram time by about two thirds, so it stays available for speed.


def gene_columns(tree, panel):
    """Columns (in panel order) of the genes a marker tree names, positive or negative."""
    names = set()

    def walk(nodes):
        for n in nodes:
            for m in n.get("markers", []):
                names.update(g.upper() for g in m["genes"])
            walk(n.get("children", []))
    walk(tree or [])
    return [i for i, g in enumerate(panel) if g.upper() in names]


_RECENT = {}  # the last few PCA and neighbor results by the content of their input (_recall)
_RECENT_LOCK = __import__("threading").Lock()
RECENT_KEEP = 3  # of each kind


def _content(*arrays):
    """A key of the arrays' exact content (shape, type and bytes)."""
    import hashlib
    h = hashlib.blake2b(digest_size=20)
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(f"{a.shape}{a.dtype.str}".encode())
        h.update(memoryview(a).cast("B"))
    return h.hexdigest()


def _recall(kind, key, make):
    """make(), or its result kept from an earlier call on the same input (e.g. the PCA of a sample's cells computed for
    the annotation, then for PCA + UMAP and for a trajectory of all its cells): the same numbers, without computing
    them again. Copies go out, so a caller that changes its arrays changes no one else's."""
    copy = lambda r: tuple(np.array(x) for x in r) if isinstance(r, tuple) else np.array(r)
    with _RECENT_LOCK:
        hit = _RECENT.get((kind, key))
        if hit is not None:
            _RECENT[(kind, key)] = _RECENT.pop((kind, key))  # most recent last
            return copy(hit)
    r = make()
    with _RECENT_LOCK:
        _RECENT[(kind, key)] = r
        mine = [k for k in _RECENT if k[0] == kind]
        for k in mine[:max(0, len(mine) - RECENT_KEEP)]:
            del _RECENT[k]
    return copy(r)


def pca(Xn, n_pcs=None, fit_max=None):
    """PCA of standardized genes on the sparse matrix (no dense copy): cell scores (cells x n_pcs, float32). The last
    few results are kept by the content of the matrix (_recall).

    fit_max: with more cells than this, the components are fitted on that many of them drawn at random (seed
    PCA_FIT_SEED) and every cell is then projected onto them, which leaves out the heaviest part of the
    computation -- the gene means and the Gram matrix over all the cells. The scores are then no longer those of a
    fit on every cell, so the clusters and the UMAP change. None (the default) fits on every cell: the marker
    annotation must not sample, and passes nothing."""
    X = sp.csr_matrix(Xn, dtype=np.float32)
    key = (_content(X.data, X.indices, X.indptr, np.array(X.shape))
           + f":{n_pcs or N_PCS}:{MAX_PCA_GENES}:{N_HVG}:{int(fit_max or 0)}:{PCA_FIT_SEED}")
    return _recall("pca", key, lambda: _pca(X, n_pcs, fit_max))


def _pca(X, n_pcs=None, fit_max=None):
    n = X.shape[0]
    rows = None
    if fit_max and n > fit_max:  # the components from a sample of the cells, all of them projected onto those
        rows = np.sort(np.random.default_rng(PCA_FIT_SEED).choice(n, int(fit_max), replace=False))
    F = X if rows is None else X[rows]  # the cells the components are fitted on
    if X.shape[1] > MAX_PCA_GENES:  # the most variable genes (of the log-normalized values)
        mu, ex2 = _moments(F)
        from . import csrutil
        keep = np.sort(np.argsort(-(ex2 - mu ** 2))[:N_HVG])
        X = csrutil.columns_f32(X, keep)  # block by block: not a second copy of the whole matrix
        F = X if rows is None else csrutil.columns_f32(F, keep)
    C, mu, ex2 = gram(F)  # the gene means and means of squares come from the same pass as the Gram matrix
    m = F.shape[0]
    del F
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0) * m / max(m - 1, 1)) + 1e-3
    cov = (C / m - np.outer(mu, mu)) / np.outer(sd, sd)  # covariance of standardized genes
    del C
    k = min(n_pcs or N_PCS, len(cov))
    import scipy.linalg
    V = scipy.linalg.eigh(cov, subset_by_index=[len(cov) - k, len(cov) - 1], driver="evr")[1]  # the top k only
    V *= np.where(V[np.abs(V).argmax(0), np.arange(k)] < 0, -1.0, 1.0)  # each one's largest entry positive: a sign
    V /= sd[:, None]  # that does not depend on the computer's LAPACK. Ascending variance: the last column is PC1
    return _scores(X, V, mu)  # = standardized X times the top components


_GRAM = None
GRAM_SLAB = 256  # genes whose cell lists (the transposed half of the matrix) are held at once; see _slabs. Every
                 # size gives the same Gram matrix to the bit; 256 costs 3% more time than 1024 and holds a twelfth
                 # as much (measured at 400,000 cells x 6,000 genes: 0.02 vs 0.34 times the matrix).


def _slabs(X):
    """The matrix column block by column block as (first gene, CSC of those genes with cells ascending, their CSR).
    Only a slab of the transpose exists at a time: a whole X.tocsc() is a second copy of the matrix, which for a
    five-million-cell CosMx 6K panel is some 13 GB beside the 13 GB of the matrix itself."""
    R = sp.csr_matrix(X)
    if not R.has_sorted_indices:  # (the log-normalized matrix shares the counts' sorted indices: no copy here)
        R = R.sorted_indices()
    for g0 in range(0, R.shape[1], GRAM_SLAB):
        S = R[:, g0:g0 + GRAM_SLAB]
        Cc = S.tocsc()
        del S
        Cc.sort_indices()
        yield g0, Cc, R
        del Cc


def _slab_moments(Cc, mu, ex2, g0):
    """The gene means and means of squares of a slab, into mu and ex2. A slab holds every cell, so each is exactly
    the entry of X.mean(0) / X.multiply(X).mean(0) for that gene, to the bit: scipy's mean scales the values by
    1 / cells and sums each column over the cells ascending, which a column block does the same way. What the slab
    avoids is the copy of the whole matrix that both make (mean scales it; multiply allocates twice the nonzeros
    before pruning them back)."""
    g1 = g0 + Cc.shape[1]
    mu[g0:g1] = np.asarray(Cc.mean(0)).ravel()
    sq = Cc._with_data(Cc.data * Cc.data, copy=False)  # the same structure, the values squared (float32, as multiply)
    ex2[g0:g1] = np.asarray(sq.mean(0)).ravel()


def _moment_arrays(X):
    """Empty gene means / means of squares, in the dtype scipy's mean returns for this matrix (float32 for a
    float32 one): everything built from them -- sd, outer(mu, mu) -- is then computed in that same dtype."""
    dt = np.asarray(X.data[:1]).dtype
    return np.zeros(X.shape[1], dt), np.zeros(X.shape[1], dt)


def _moments(X):
    """The gene means and means of squares alone (for picking the most variable genes of whole-transcriptome data,
    before the Gram matrix of the kept ones)."""
    X = sp.csr_matrix(X)
    g, mu, ex2 = X.shape[1], *_moment_arrays(X)
    for g0, Cc, _ in _slabs(X):
        _slab_moments(Cc, mu, ex2, g0)
    return mu, ex2


def _mirror(C):
    """The upper half of C mirrored into the lower one, block by block: C + np.triu(C, 1).T needs three more
    arrays the size of C (at 6,000 genes, 275 MB each)."""
    g = len(C)
    for a in range(0, g, 1024):
        b = min(a + 1024, g)
        blk = C[a:b, a:b]
        C[a:b, a:b] = blk + np.triu(blk, 1).T  # the diagonal block
        if b < g:
            C[b:, a:b] = C[a:b, b:].T          # the block below it
    return C


def gram(X):
    """(X.T @ X, gene means, gene means of squares): the first as a dense float64 array, on up to 8 threads (numba).
    Each entry is summed over the cells in their order, so the same on every computer and for any number of threads,
    and exact in float64 (scipy's sparse product sums in float32 on one thread: 3 times slower on 8 cores for a
    5,000-gene panel, 21.6 -> 6.3 s on 188,644 cells, and off by 3e-5). The means come from the same pass."""
    global _GRAM
    import os
    import numba
    if _GRAM is None:
        from numba import njit, prange
        cache = bool(globals().get("__file__")) and os.path.isfile(globals()["__file__"])  # (not in a notebook's sa)

        @njit(parallel=True, cache=cache)
        def _upper(cp, ci, cv, rp, ri, rv, g0, C):
            for t in prange(len(cp) - 1):  # row g0 + t of the upper half from column g0 + t: cells ascending,
                i = g0 + t                 # one thread per row
                for a in range(cp[t], cp[t + 1]):
                    c, v = ci[a], np.float64(cv[a])
                    lo, hi = rp[c], rp[c + 1]
                    s = lo + np.searchsorted(ri[lo:hi], i)  # (the row's genes are sorted: from gene i on)
                    for b in range(s, hi):
                        C[i, ri[b]] += v * np.float64(rv[b])
        _GRAM = _upper
    Xc = sp.csr_matrix(X)
    g = Xc.shape[1]
    C, (mu, ex2) = np.zeros((g, g)), _moment_arrays(Xc)
    with numba_parallel():
        before = numba.get_num_threads()
        try:
            numba.set_num_threads(max(1, min(KNN_THREADS, int(numba.config.NUMBA_NUM_THREADS))))
            for g0, Cc, R in _slabs(X):
                _slab_moments(Cc, mu, ex2, g0)
                _GRAM(Cc.indptr, Cc.indices, Cc.data, R.indptr, R.indices, R.data, g0, C)
        finally:
            numba.set_num_threads(before)
    return _mirror(C), mu, ex2


SCORE_BLOCK = 50_000  # cells whose principal-component scores are computed at once; see _scores


def _scores(X, V, mu):
    """The cells' scores on the components, cell block by cell block: X @ V upcasts the whole matrix to float64
    (a second copy of the data) because V is float64. Each cell's scores are its own dot products, so the blocks
    give the same numbers."""
    X = sp.csr_matrix(X)
    Z = np.empty((X.shape[0], V.shape[1]), np.float32)
    off = mu @ V
    for a in range(0, X.shape[0], SCORE_BLOCK):
        b = min(a + SCORE_BLOCK, X.shape[0])
        lo, hi = X.indptr[a], X.indptr[b]
        blk = sp.csr_matrix((X.data[lo:hi], X.indices[lo:hi], X.indptr[a:b + 1] - lo),  # a view, not a copy
                            shape=(b - a, X.shape[1]))
        Z[a:b] = blk @ V - off
    return Z


KNN = "exact"  # neighbor search of the marker step: "exact" (default: the true neighbors, the same on every computer) or "approx" (see knn)
KNN_THREADS = 8  # the approximate search always runs on this many threads (its result depends on the count)
KNN_EXTRA = 8  # candidates beyond k that the exact search re-ranks in float64


def has_torch():
    try:
        import torch  # noqa: F401
        return True
    except ImportError:
        return False


def knn_backend(mode):
    """What knn(mode=...) computes with on this computer: (key, description). The key names the step in timing."""
    from .config import device
    import importlib.metadata as md
    if mode == "approx":
        return "knn_approx", f"approximate (pynndescent {md.version('pynndescent')}, numba {md.version('numba')})"
    dev = device()
    if dev != "cpu":
        return "knn_exact_gpu", f"exact, on the {'CUDA' if dev == 'cuda' else 'Apple'} GPU"
    if _mlx_gpu():  # the macOS app (no torch): MLX
        return "knn_exact_mlx", "exact, on the Apple GPU (MLX)"
    if has_torch():
        return "knn_exact_torch", "exact, on the CPU (torch)"
    return "knn_exact_numpy", "exact, on the CPU (numpy)"


def _expand_copies(Zu, cu, inv, cnt, kk):
    """Candidates of every row from those of the distinct points Zu (cu: kx of them per point, any order): the points
    in order of distance (float64; ties by point), each standing for its first kk rows by row number, the first kk
    rows of that list. So a row's candidates are its kk nearest rows, ties at one point by row number, as a search
    over all rows gives them (the re-rank and its check in _exact then make them exact)."""
    u = len(Zu)
    order = np.argsort(inv, kind="stable")  # rows grouped by their point, ascending within a point
    start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
    j = np.arange(kk)
    mem = np.where(j < cnt[:, None], order[np.minimum(start[:, None] + j, len(order) - 1)], -1)  # u x kk, -1 padded
    Zu64 = Zu.astype(np.float64)
    out = np.empty((u, kk), np.int64)
    for i in range(0, u, 4096):
        c = cu[i:i + 4096]
        d = ((Zu64[c] - Zu64[i:i + len(c), None]) ** 2).sum(-1)
        c = np.take_along_axis(c, np.lexsort((c, d), axis=1), 1)  # nearest point first (the point itself at 0)
        m = mem[c].reshape(len(c), -1)
        first = np.argsort(m < 0, axis=1, kind="stable")[:, :kk]  # the rows of the nearest points, in that order
        out[i:i + len(c)] = np.take_along_axis(m, first, 1)
    cand = out[inv]
    bad = cand < 0  # fewer than kk rows in all the points found (only with very few distinct points): any rows
    if bad.any():
        fill = np.arange(len(inv))
        for r in np.flatnonzero(bad.any(1)):
            have = cand[r][cand[r] >= 0]
            cand[r] = np.concatenate([have, np.setdiff1d(fill, have)[:kk - len(have)]])
    return cand


def _exact(Z, k):
    """Exact neighbors, the same on every computer: candidates from a GPU (CUDA or Apple MPS) when torch has one,
    else from the Apple GPU through MLX (the macOS app), else from the CPU (torch, or numpy without it), in float32, re-ranked by their float64 distance with ties going to
    the lower row number. A row whose k-th distance is not clearly below its farthest candidate's (by more than
    float32 rounding can move a distance: many exact duplicates, say) is searched again over all rows that rounding
    could put within its k-th distance, so rounding differences between the backends do not change the result."""
    import warnings
    from .config import device
    from . import timing
    n, dim = Z.shape
    if k <= 0 or n < 2:  # (a single row has no neighbors)
        return np.empty((n, 0), np.float32), np.empty((n, 0), np.int64)
    kk = min(k + 1 + KNN_EXTRA, n)
    Z64 = np.asarray(Z, np.float64)
    with warnings.catch_warnings():  # (NaN in the data)
        warnings.simplefilter("ignore")
        mu = np.nan_to_num(np.nanmean(Z64, 0), nan=0.0, posinf=0.0, neginf=0.0)
    Zc = (Z64 - mu).astype(np.float32)  # near the origin float32 rounds less
    key, _ = knn_backend("exact")
    # many identical rows (cells without transcripts, or with one, in unfiltered data: often thousands): the distinct
    # points are searched, and each row gets the members of its nearest points (its own first): copies cost nothing
    Zu, inv, cnt = np.unique(Zc, axis=0, return_inverse=True, return_counts=True)
    inv = np.asarray(inv).ravel()
    dedup = n - len(Zu) > max(1000, n // 100)
    X, kx = (Zu, min(kk, len(Zu))) if dedup else (Zc, kk)
    cand = None
    if key == "knn_exact_gpu":
        try:
            with timing.step(key, len(X) ** 2 / 1e9):
                cand = _topk_torch(X, kx, device())
        except RuntimeError:  # e.g. GPU memory taken by other jobs: the CPU gives the same neighbors
            key = "knn_exact_torch"
    elif key == "knn_exact_mlx":
        try:
            with timing.step(key, len(X) ** 2 / 1e9):
                cand = _topk_mlx(X, kx)
        except Exception as e:  # any failure of MLX or the GPU: the CPU gives the same neighbors
            if not isinstance(e, ValueError):  # (ValueError: NaN in the data, for the CPU search)
                print(f"exact neighbors: the Apple GPU (MLX) failed ({type(e).__name__}: {e}); on the CPU instead")
            key = "knn_exact_numpy"
    if cand is None:
        with timing.step(key, len(X) ** 2 / 1e9):
            cand = _topk_torch(X, kx, "cpu") if key == "knn_exact_torch" else _topk_numpy(X, kx)
    if dedup:
        cand = _expand_copies(Zu, cand, inv, cnt, kk)
    dist, nbr = np.empty((n, k), np.float32), np.empty((n, k), np.int64)
    dk, dmax = np.empty(n), np.empty(n)  # squared: the k-th neighbor's and the farthest candidate's

    def rerank(i, j):  # 4096 rows x 24 candidates x 30 dimensions in float64: about 25 MB per thread
        c = cand[i:j]
        x = Z64[c]
        x -= Z64[i:j, None]
        x *= x
        d = x.sum(-1)
        dmax[i:j] = d.max(1)
        d[c == np.arange(i, j)[:, None]] = np.inf  # not its own neighbor
        o = np.lexsort((c, d), axis=1)[:, :k]
        nbr[i:j] = np.take_along_axis(c, o, 1)
        d = np.take_along_axis(d, o, 1)
        dk[i:j] = d[:, -1]
        dist[i:j] = np.sqrt(d)

    def again(i, c):  # every row that float32 rounding could put within the k-th distance
        x = Z64[c]
        x -= Z64[i]
        x *= x
        d = x.sum(-1)
        d[c == i] = np.inf
        if np.isfinite(d).sum() >= k:  # (else NaN in the data)
            o = np.lexsort((c, d))[:k]
            nbr[i], dist[i] = c[o], np.sqrt(d[o])
    with timing.step("knn_rerank", n / 1e6):
        _blocks(n, 4096, cpu_threads(), rerank)
        if kk < n:
            a = np.sqrt((Zc.astype(np.float64) ** 2).sum(1))
            # float32 rounding moves each compared value by at most about 4 (dim + 4) 2^-24 (|z_i| + |z_j|)^2 (centered;
            # a row nearer than the k-th has |z_j| <= |z_i| + its distance): twice that for the two values compared,
            # times a safety factor 2 (torch computes them slightly differently)
            margin = 16 * (dim + 4) * 2.0 ** -24 * (a + np.sqrt(dmax)) ** 2
            redo = np.flatnonzero(~(dmax - dk > margin))
            if len(redo):
                # a row with at least k identical rows (cells without transcripts, say: often thousands): its k
                # nearest are the k lowest-numbered of them, at distance 0, without searching (a search over every
                # row under a zero bound would cost the square of their number)
                _, inv, cnt = np.unique(Z64, axis=0, return_inverse=True, return_counts=True)
                inv = np.asarray(inv).ravel()
                same = redo[cnt[inv[redo]] > k]
                if len(same):
                    order = np.argsort(inv, kind="stable")  # the rows grouped by their point, ascending in a group
                    start = np.concatenate([[0], np.cumsum(cnt)[:-1]])
                    for i in same:
                        g = order[start[inv[i]]:start[inv[i]] + k + 1]
                        nbr[i], dist[i] = g[g != i][:k], 0.0
                    redo = redo[cnt[inv[redo]] <= k]
            if len(redo):
                _within(Zc, redo, dk - a * a + margin, again)
    return dist, nbr


def knn(Z, k, mode="auto"):
    """k nearest neighbors of every row (itself excluded): (distances, indices).
    mode "exact": the true neighbors, the same on every computer whatever it computes with (see _exact).
    "approx": pynndescent, fast on a CPU, on KNN_THREADS threads (its result depends on the thread count): the same
    result on every computer for the same pynndescent and numba versions (x86 and ARM may still differ slightly).
    "auto" (embeddings): exact on a GPU and for small sets, else pynndescent."""
    from .config import device
    from . import timing
    n = len(Z)
    k = min(k, n - 1)
    exact = mode == "exact" or (mode == "auto" and (device() != "cpu" or n <= 20000))
    if exact:  # (the last few results are kept by the content of Z: _recall)
        return _recall("knn", f"{_content(Z)}:{k}:exact:{KNN_EXTRA}", lambda: _exact(Z, k))

    def approx():
        with timing.step("knn_approx", n / 1e6, warmup="approx_compile"):
            nbr, dist = _nndescent(Z, k + 1)
        return dist[:, 1:], nbr[:, 1:]
    return _recall("knn", f"{_content(Z)}:{k}:approx:{KNN_THREADS}", approx)


_NUMBA_LOCK = __import__("threading").Lock()


@__import__("contextlib").contextmanager
def numba_parallel():
    """One numba-parallel run at a time in this process (pynndescent, UMAP, the startup compile): with numba's
    "workqueue" threading layer (the one a Mac or Windows build may get) two threads running parallel code at once abort
    the whole app. Waits in short steps, so a Stop still reaches a run that is waiting."""
    while not _NUMBA_LOCK.acquire(timeout=0.25):
        pass
    try:
        yield
    finally:
        _NUMBA_LOCK.release()


def prewarm_approx():
    """Compile pynndescent's numba code now (a small search in this process: about 10 s here, longer the first time
    on a Mac), so that the first approximate annotation of this app session does not wait for it; from then on the
    time estimates leave the compile out. Run in a background thread when the app starts (one core, briefly more)."""
    import time
    from . import timing
    if "approx_compile" in timing._warm:
        return
    t = time.perf_counter()
    try:
        Z = np.random.default_rng(0).standard_normal((3000, 20)).astype(np.float32)
        _nndescent(Z, 16)
    except Exception as e:  # the first real run compiles instead
        print(f"approximate neighbors: could not compile ahead ({type(e).__name__}: {e})")
        return
    timing._warm.add("approx_compile")
    print(f"approximate neighbors: compiled ahead in {time.perf_counter() - t:.1f} s")


def knn_threads():
    """The threads the approximate search runs on: KNN_THREADS (stagent/__init__.py lets numba have that many on any
    computer); fewer only if numba was loaded with fewer allowed (then other computers may differ)."""
    import numba
    return min(KNN_THREADS, int(numba.config.NUMBA_NUM_THREADS))


def _nndescent(Z, k):
    """pynndescent's neighbors on a fixed number of threads (its result depends on how many), whatever the number of
    cores; the thread count of the calling thread is put back afterwards (UMAP and the rest keep theirs)."""
    import numba
    from pynndescent import NNDescent
    with numba_parallel():
        before = numba.get_num_threads()
        try:
            try:
                return NNDescent(Z, n_neighbors=k, random_state=0, n_jobs=knn_threads()).neighbor_graph
            except ValueError as e:  # a thread setting numba refuses: the default, rather than no result
                print(f"approximate neighbors: {e}; using numba's default threads")
                return NNDescent(Z, n_neighbors=k, random_state=0).neighbor_graph
        finally:
            numba.set_num_threads(before)


def knn_graph(dist, nbr):
    """Weighted undirected igraph graph from kNN lists (adaptive Gaussian weights)."""
    import igraph
    n, k = nbr.shape
    w = np.exp(-(dist / np.maximum(dist[:, -1], 1e-6)[:, None]) ** 2).ravel()
    A = sp.csr_matrix((w, (np.repeat(np.arange(n), k), nbr.ravel())), shape=(n, n))
    A = A.maximum(A.T).tocoo()
    m = A.row < A.col
    g = igraph.Graph(n=n, edges=np.c_[A.row[m], A.col[m]])
    g.es["weight"] = A.data[m]
    return g


def fingerprint(a):
    """A short checksum of an array (or a list of strings): two runs that give the same have the same."""
    import hashlib
    if isinstance(a, np.ndarray) and a.dtype != object:
        b = np.ascontiguousarray(a).tobytes()
    else:
        b = "\n".join(map(str, a)).encode("utf-8")
    return hashlib.sha256(b).hexdigest()[:12]


def _graph(Xn, knn_mode="auto", fp=None):
    """PCA on standardized genes and a kNN graph with adaptive Gaussian weights; built once per level. fp: a dict
    that gets the checksums of the PCA and of the neighbors (the root level, for checking that runs agree)."""
    from . import timing
    m = Xn.shape[0]
    with timing.step("pca", m * Xn.shape[1] / 1e9):
        Z = pca(Xn)
    dist, nbr = knn(Z, N_NEIGHBORS, knn_mode)
    if fp is not None:
        fp.update(pca=fingerprint(Z), neighbors=fingerprint(np.asarray(nbr, np.int64)))
    with timing.step("graph", m / 1e6):
        g = knn_graph(dist, nbr)
    return g, nbr


_LEIDEN_LOCK = __import__("threading").Lock()


def _leiden(g, nbr, resolution, seed, min_cluster=None):
    """igraph's built-in Leiden (modularity objective; same algorithm as leidenalg, several times faster).
    Clusters smaller than min_cluster (default MIN_CLUSTER) are merged into the cluster most of their neighbors
    belong to."""
    MIN = MIN_CLUSTER if min_cluster is None else min_cluster
    import random
    import igraph
    with _LEIDEN_LOCK:  # igraph draws from Python's global random: two runs at once (an annotation and a PCA / UMAP
        # in the app's threads) would take each other's numbers and depend on timing. igraph binds Python's random only
        # in the thread that imported it first; every other thread (each job of the app) used igraph's own unseeded
        # generator, so annotating twice gave other clusters: bound here, in this thread, before seeding
        igraph.set_random_number_generator(random)
        random.seed(seed)
        lab = np.array(g.community_leiden(objective_function="modularity", weights="weight", resolution=resolution,
                                          n_iterations=2).membership)
    sizes = np.bincount(lab)
    small = sizes[lab] < MIN
    if small.any() and (~small).any():
        for i in np.where(small)[0]:
            cand = lab[nbr[i]][sizes[lab[nbr[i]]] >= MIN]
            if len(cand):
                lab[i] = np.bincount(cand).argmax()
    return np.unique(lab, return_inverse=True)[1]


def _detected(counts):
    """Cells x genes, 1.0 where the gene is detected (count > 0), float32. A counts matrix in canonical form without
    stored zeros (as samples are read) shares its index arrays (no copy of them: a CosMx 6K or Visium HD matrix is
    gigabytes); else the same matrix computed from the comparison."""
    if counts.has_canonical_format and (counts.data > 0).all():
        return sp.csr_matrix((np.ones(counts.nnz, np.float32), counts.indices, counts.indptr), shape=counts.shape)
    return (counts > 0).astype(np.float32).tocsr()


def _all_rows(idx, n):
    """idx is every row in order (the root level): the matrix itself serves, without a copy of it."""
    idx = np.asarray(idx)
    return len(idx) == n and (n == 0 or (idx[0] == 0 and idx[-1] == n - 1 and bool((np.diff(idx) == 1).all())))


class _Annotator:
    """Top-down marker annotation over one sample; also re-runs single levels or clusters in place."""

    def __init__(self, counts, panel, progress=None, seed=0, knn_mode=KNN):
        from .model import lognorm
        from . import timing
        with timing.step("prepare", counts.shape[0] / 1e6):
            counts = sp.csr_matrix(counts)
            self.Xn, self.det = lognorm(counts), _detected(counts)
            self.counts = counts  # (the rare-type search reads transcript counts of a few marker genes)
        self.knn_mode = knn_mode
        self.gi = {g: i for i, g in enumerate(panel)}
        self.say = progress or (lambda f, m: None)
        self.seed, self.n = seed, counts.shape[0]
        self.done, self.total = 0, 1
        self.fp = {}  # checksums of the root level (see annotate)

    def _evaluate(self, idx, lab, children, context=None):
        """Score the clusters lab of the rows idx for each child type: the detection rate of each marker, z-scored
        across the clusters, sign-adjusted, averaged over the child's markers (A|B: the strongest). context: the rows
        of each of the level's other clusters, when only some of its cells are clustered again (rerun of one cluster):
        the markers are z-scored over all of them, as in the first run, not over the new sub-clusters alone (between
        near-identical sub-clusters the z-scores are noise)."""
        k = int(lab.max()) + 1
        onehot = sp.csr_matrix((np.ones(len(idx)), (lab, np.arange(len(idx)))), shape=(k, len(idx)))
        size = np.asarray(onehot.sum(1)).ravel()
        frac = np.asarray((onehot @ (self.det if _all_rows(idx, self.n) else self.det[idx])).todense()) / np.maximum(size, 1)[:, None]
        if context:  # the level's other clusters
            frac = np.vstack([frac] + [np.asarray(self.det[np.asarray(rows)].mean(0)).reshape(1, -1) for rows in context])
        S = np.zeros((k, len(children)))
        for j, ch in enumerate(children):
            parts = []
            for m in ch["markers"]:
                f = frac[:, [self.gi[g] for g in m["genes"]]]
                z = np.clip((f - f.mean(0)) / (f.std(0) + 0.02), -3, 3)[:k]
                parts.append(m["sign"] * (z.max(1) if m["sign"] > 0 else z.min(1)))
            S[:, j] = np.mean(parts, axis=0)
        order = np.argsort(-S, 1)
        best = S[np.arange(k), order[:, 0]]
        second = S[np.arange(k), order[:, 1]] if len(children) > 1 else np.full(k, -np.inf)
        ok = (best >= MIN_SCORE) & (best - second >= MIN_MARGIN)
        return {"k": k, "size": size, "S": S, "order": order, "best": best, "second": second, "ok": ok,
                "assigned_frac": float(size[ok].sum() / len(idx))}

    def cluster(self, parent, idx, children, resolution, context=None):
        """Leiden at one level. With resolution "auto", every candidate in RESOLUTIONS is tried on the same graph and
        the one assigning the largest share of cells to a child type wins (ties within 1% go to the lower one)."""
        label = parent or "all cells"
        self.say(self.done / self.total, f"building the neighbor graph of {label} ({len(idx):,} cells)")
        from . import timing
        root = parent is None and len(idx) == self.n
        g, nbr = _graph(self.Xn if _all_rows(idx, self.n) else self.Xn[idx], self.knn_mode, self.fp if root else None)
        self.graph = (g, nbr)  # (kept for the rare-type search of this level)
        cands = RESOLUTIONS if resolution == AUTO else (float(resolution),)
        best, tried = None, {}
        for r in cands:
            self.say(self.done / self.total, f"Leiden on {label}, resolution {r:g}")
            with timing.step("leiden", len(idx) / 1e6):
                lab = _leiden(g, nbr, r, self.seed)
                ev = self._evaluate(idx, lab, children, context)
            tried[f"{r:g}"] = round(ev["assigned_frac"], 3)
            if best is None or ev["assigned_frac"] > best[1]["assigned_frac"] + 0.01:
                best = (r, ev, lab)
        if root:
            self.fp.update(root_resolution=best[0], root_clusters=fingerprint(np.asarray(best[2], np.int64)))
        return best[0], best[1], best[2], tried

    def run(self, res, parent, idx, children, resolution, prefix="", context=None):
        """Annotate idx under `parent` with its `children` types, then recurse. Updates res in place. context: see
        _evaluate (this level only)."""
        if not children or len(idx) < MIN_CELLS:
            return
        r, ev, lab, tried = self.cluster(parent, idx, children, resolution, context)
        res["nodes"] = [x for x in res["nodes"] if not (x["node"] == (parent or "root") and x.get("prefix") == prefix)]
        res["nodes"].append({"node": parent or "root", "prefix": prefix, "n_cells": int(len(idx)), "resolution": r,
                             "tried": tried, "assigned_frac": round(ev["assigned_frac"], 3)})
        for c in range(ev["k"]):
            assigned = children[ev["order"][c, 0]]["name"] if ev["ok"][c] else None
            cid = f"{prefix}{parent or 'root'}/c{c}"
            rows = idx[lab == c]
            if assigned:
                res["labels"][rows] = assigned
            res["cluster_of"][rows] = cid
            res["members"][cid] = rows.astype(np.int32)
            rec = {"id": cid, "parent": parent, "n_cells": int(ev["size"][c]), "assigned": assigned, "resolution": r,
                   "scores": {ch["name"]: round(float(ev["S"][c, j]), 2) for j, ch in enumerate(children)},
                   "margin": round(float(ev["best"][c] - ev["second"][c]), 2) if len(children) > 1 else None}
            res["records"].append(rec)
        if RARE:
            self.rare(res, parent, idx, children, ev, r, prefix, context)
        self.done += 1
        for ch in children:
            self.run(res, ch["name"], idx[res["labels"][idx] == ch["name"]], ch["children"], resolution, prefix)

    def _groups_det(self, idx, rows_by_cluster, groups):
        """Detection rate (any gene of the group detected) of each marker group in each set of rows, and in idx."""
        det = self.det if _all_rows(idx, self.n) else self.det[idx]
        out = np.zeros((len(rows_by_cluster) + 1, len(groups)))
        for j, gr in enumerate(groups):
            hit = np.asarray(det[:, [self.gi[x] for x in gr["genes"]]].sum(1)).ravel() > 0
            out[-1, j] = hit.mean()
            for c, rows in enumerate(rows_by_cluster):
                out[c, j] = hit[rows].mean()
        return out

    def not_found(self, res, tree):
        """The cell types of the list no cell went to, each with the detection of its positive markers among the cells
        of its parent (the level it was looked for in): [{type, parent, n_parent_cells, markers: {genes: [cells with
        one or more transcripts, with two or more]}, reason}]."""
        out = []
        D = descendants(tree)
        have = set(res["labels"].tolist())

        def walk(nodes, parent):
            for nd in nodes:
                if nd["children"]:
                    walk(nd["children"], nd["name"])
                elif nd["name"] not in have:
                    pool = np.flatnonzero(np.isin(res["labels"], list(D.get(parent, [])) + [parent])) if parent \
                        else np.arange(self.n)
                    item = {"type": nd["name"], "parent": parent, "n_parent_cells": int(len(pool)), "markers": {}}
                    if not len(pool):
                        item["reason"] = f"its parent {parent} has no cells"
                        out.append(item)
                        continue
                    C = self.counts[pool]
                    for m in nd["markers"]:
                        if m["sign"] > 0:
                            v = np.asarray(C[:, [self.gi[x] for x in m["genes"]]].max(1).todense()).ravel()
                            item["markers"]["|".join(m["genes"])] = [int((v > 0).sum()), int((v >= 2).sum())]
                    best = max(item["markers"].items(), key=lambda kv: kv[1][1], default=None)
                    item["reason"] = ("no cluster went to it, also at finer resolutions, and no group of cells shows "
                                      "its markers together" + (
                        f"; among the {len(pool):,} {parent or 'cells'}: " + ", ".join(
                            f"{g} in {a:,} cells ({b:,} with 2+ transcripts)" for g, (a, b) in item["markers"].items())
                        if best else ""))
                    out.append(item)
        walk(tree, None)
        return out

    def _rare_cells(self, idx, ch):
        """(rows of the level's cells that show ch's markers clearly, how) or (None, how); see RARE_SPECIFIC."""
        pos = [m for m in ch["markers"] if m["sign"] > 0]
        neg = [m for m in ch["markers"] if m["sign"] < 0]
        if not pos:
            return None, {"accepted": False, "rule": "no positive marker"}
        C = self.counts if _all_rows(idx, self.n) else self.counts[idx]
        V = np.column_stack([np.asarray(C[:, [self.gi[x] for x in m["genes"]]].max(1).todense()).ravel() for m in pos])
        ok = np.ones(len(idx), bool)
        for m in neg:
            ok &= np.asarray(C[:, [self.gi[x] for x in m["genes"]]].sum(1)).ravel() == 0
        det = (V > 0).mean(0)
        if len(pos) == 1:
            j = 0
            if det[j] > RARE_SPECIFIC:
                return None, {"accepted": False, "rule": f"{'|'.join(pos[j]['genes'])} is not specific here "
                              f"(detected in {det[j]:.1%} of the cells)"}
            hit = ok & (V[:, j] >= RARE_CELL_COUNTS)
            rule = f"{'|'.join(pos[j]['genes'])} >= {RARE_CELL_COUNTS} transcripts"
            expect = None
        else:
            hi = V >= RARE_PAIR_COUNTS
            hit = ok & (hi.sum(1) >= 2)
            p = hi.mean(0)
            expect = float((np.sum(p) ** 2 - np.sum(p ** 2)) / 2 * len(idx))  # two markers at once by chance
            rule = f"two of {', '.join('|'.join(m['genes']) for m in pos)} with >= {RARE_PAIR_COUNTS} transcripts each"
        n = int(hit.sum())
        how = {"rule": rule, "n_cells": n, "expected": None if expect is None else round(expect, 1)}
        good = n >= RARE_MIN_CELLS and (expect is None or n >= RARE_FOLD * expect)
        how["accepted"] = bool(good)
        return (idx[hit] if good else None), how

    def rare(self, res, parent, idx, children, ev, r0, prefix, context=None):
        """The rare-type search of one level (see RARE): children no cluster went to, looked for in finer clusterings
        of the same graph; the first resolution that finds a type with real marker evidence gives it its clusters."""
        found = {children[ev["order"][c, 0]]["name"] for c in range(ev["k"]) if ev["ok"][c]}
        missing = [j for j, ch in enumerate(children) if ch["name"] not in found]
        if not missing or len(children) < 2 or len(idx) < 2 * RARE_MIN_CLUSTER:
            return
        g, nbr = self.graph
        label = parent or "all cells"
        tried = {}
        for rr in [x for x in RARE_RESOLUTIONS if x > r0]:
            if not missing:
                break
            self.say(self.done / self.total, f"looking for rare types of {label}: Leiden at resolution {rr:g}")
            lab = _leiden(g, nbr, rr, self.seed, RARE_MIN_CLUSTER)
            ev2 = self._evaluate(idx, lab, children, context)
            got = {}
            for c in range(ev2["k"]):
                j = int(ev2["order"][c, 0])
                if ev2["ok"][c] and j in missing:
                    got.setdefault(j, []).append(c)
            for j, cs in list(got.items()):
                ch = children[j]
                pos = [m for m in ch["markers"] if m["sign"] > 0]
                neg = [m for m in ch["markers"] if m["sign"] < 0]
                rows = [np.flatnonzero(lab == c) for c in cs]
                dp = self._groups_det(idx, rows, pos) if pos else np.zeros((len(cs) + 1, 0))
                dn = self._groups_det(idx, rows, neg) if neg else np.zeros((len(cs) + 1, 0))
                keep = []
                for k, c in enumerate(cs):
                    strong = (((dp[k] >= RARE_MIN_DET) & (dp[k] >= RARE_FOLD * dp[-1]))  # a specific marker, or a
                              | ((dp[k] >= RARE_COMMON_DET) & (dp[k] >= RARE_COMMON_FOLD * dp[-1])))  # common one, raised
                    not_less = dp[k] >= 0.8 * dp[-1]
                    okn = dn[k] <= np.maximum(1.5 * dn[-1], dn[-1] + 0.05)
                    good = bool(pos) and bool(strong.any()) and bool(not_less.all()) and bool(okn.all())
                    tried.setdefault(ch["name"], []).append(
                        {"resolution": rr, "n_cells": int(len(rows[k])), "accepted": good,
                         "positive": {"|".join(m["genes"]): [round(float(dp[k, i]), 3), round(float(dp[-1, i]), 3)]
                                      for i, m in enumerate(pos)},
                         "negative": {"|".join(m["genes"]): [round(float(dn[k, i]), 3), round(float(dn[-1, i]), 3)]
                                      for i, m in enumerate(neg)}})
                    if good:
                        keep.append(k)
                if not keep:
                    continue
                for k in keep:
                    c = cs[k]
                    cells = idx[rows[k]]
                    cid = f"{prefix}{parent or 'root'}/rare{rr:g}c{c}"
                    for old in set(res["cluster_of"][cells].tolist()):  # out of the clusters they were in
                        if old in res["members"]:
                            m = np.asarray(res["members"][old])
                            res["members"][old] = m[~np.isin(m, cells)].astype(np.int32)
                    res["labels"][cells] = ch["name"]
                    res["cluster_of"][cells] = cid
                    res["members"][cid] = cells.astype(np.int32)
                    res["records"].append({"id": cid, "parent": parent, "n_cells": int(len(cells)), "assigned": ch["name"],
                                           "resolution": rr, "rare": True,
                                           "scores": {x["name"]: round(float(ev2["S"][c, i]), 2) for i, x in enumerate(children)},
                                           "margin": round(float(ev2["best"][c] - ev2["second"][c]), 2)})
                missing.remove(j)
        for j in list(missing):  # still no cells: the single cells that show the type's markers clearly
            ch = children[j]
            cells, how = self._rare_cells(idx, ch)
            tried.setdefault(ch["name"], []).append({"resolution": None, "single_cells": True, **how})
            if cells is None:
                continue
            cid = f"{prefix}{parent or 'root'}/rare-cells:{ch['name']}"
            for old in set(res["cluster_of"][cells].tolist()):
                if old in res["members"]:
                    m = np.asarray(res["members"][old])
                    res["members"][old] = m[~np.isin(m, cells)].astype(np.int32)
            res["labels"][cells] = ch["name"]
            res["cluster_of"][cells] = cid
            res["members"][cid] = cells.astype(np.int32)
            res["records"].append({"id": cid, "parent": parent, "n_cells": int(len(cells)), "assigned": ch["name"],
                                   "resolution": None, "rare": "cells", "scores": {}, "margin": None,
                                   "evidence": how["rule"]})
            missing.remove(j)
        for rec in res["records"]:  # the clusters the rare cells left: their sizes now
            if rec["id"] in res["members"] and not rec.get("rare"):
                rec["n_cells"] = int(len(res["members"][rec["id"]]))
        node = next(x for x in res["nodes"] if x["node"] == (parent or "root") and x.get("prefix") == prefix)
        node["rare"] = tried


def training_labels(res, types):
    """The expression model's training label of every cell from a marker-step result: the index of its cell type, or
    -1 for a cell left at a lineage or parent level (the session and the exported notebook both use this)."""
    index = {t: i for i, t in enumerate(types)}
    return np.array([index.get(x, -1) for x in res["labels"]], np.int64)


def restrict(logp, res, tree, types, corrected=()):
    """The expression model's log-probabilities (cells x types) made to agree with the marker step, as float16 (the
    session and the exported notebook both use this): a cell left at a lineage or parent level gets a type below it
    (the others down-weighted by 1e-4); a type the model did not learn (no cell of it in the marker step) is not given,
    unless corrected holds it (type indices the user corrected cells to): its cells would be noise."""
    types = list(types)
    groups = descendants(tree)
    lp = np.asarray(logp, np.float32).copy()
    lab = res["labels"]
    for node in set(lab.tolist()) - set(types):  # left at a group: a type of that group
        rows = np.where(lab == node)[0]
        lp[np.ix_(rows, np.where(~np.isin(types, groups.get(node, [])))[0])] += np.float32(np.log(1e-4))
    found = {types[i] for i in np.unique(training_labels(res, types)) if i >= 0}
    found |= {types[int(i)] for i in np.asarray(corrected, np.int64)}
    off = np.array([t not in found for t in types])
    if off.any() and not off.all():
        lp[:, off] = np.float32(-1e4)
    lp -= lp.max(1, keepdims=True)
    return (lp - np.log(np.exp(lp).sum(1, keepdims=True))).astype(np.float16)


def _n_levels(tree):
    return 1 + sum(1 for nd in descendants(tree) if _find(tree, nd)["children"])


def annotate(counts, panel, tree, progress=None, seed=0, resolution=AUTO, knn=KNN):
    """Full top-down annotation. Returns {"labels" (type name per cell), "cluster_of" (deepest cluster id per cell),
    "members" (cluster id -> cell rows), "records" (one per cluster), "nodes" (resolution chosen per level)}.
    knn: the neighbor search, "exact" or "approx" (see knn)."""
    A = _Annotator(counts, panel, progress, seed, knn)
    A.total = _n_levels(tree)
    res = {"labels": np.full(A.n, UNASSIGNED, dtype=object), "cluster_of": np.full(A.n, "", dtype=object),
           "members": {}, "records": [], "nodes": []}
    A.run(res, None, np.arange(A.n), tree, resolution)
    res["not_found"] = A.not_found(res, tree)
    res["knn"] = knn_backend(knn)[1] if knn != "auto" else "automatic"
    res["fingerprints"] = run_fingerprints(A.fp, res["labels"], knn)
    A.say(1.0, "marker annotation done")
    return res


def run_fingerprints(fp, labels, knn):
    """What an annotation run computed, as checksums, and what it ran with: two runs of the same data and marker list
    should agree step by step; where they part shows which step differed (PCA, neighbors, root clusters, labels)."""
    import importlib.metadata as md, platform
    ver = {}
    for p in ("numpy", "scipy", "pynndescent", "numba", "igraph"):
        try:
            ver[p] = md.version(p)
        except md.PackageNotFoundError:
            pass
    return {**fp, "labels": fingerprint(labels), "knn": knn,
            "knn_threads": knn_threads() if knn in ("approx", "auto") else None,
            "machine": f"{platform.system()} {platform.machine()}", "cores": __import__("os").cpu_count(), "versions": ver}


def rerun(res, counts, panel, tree, target, resolution=AUTO, prefix="", progress=None, seed=0, knn=KNN):
    """Re-annotate part of a result: target "root", a type name (its cells are split again into its subtypes) or a
    cluster id (only that cluster's cells are re-clustered at its level). Returns a new result; res is not modified."""
    out = {"labels": res["labels"].copy(), "cluster_of": res["cluster_of"].copy(), "members": dict(res["members"]),
           "records": list(res["records"]), "nodes": list(res["nodes"]), "knn": res.get("knn"),
           "not_found": list(res.get("not_found") or [])}
    rec = {r["id"]: r for r in res["records"]}
    context = None
    if target == "root":
        parent, children = None, tree
        idx = np.arange(len(out["labels"]))
    elif target in rec:
        parent = rec[target]["parent"]
        node = _find(tree, parent) if parent else None
        children = node["children"] if node else tree
        idx = np.asarray(out["members"].get(target, np.zeros(0, np.int32)))
        # the level's other clusters, for z-scoring the markers as the first run did
        sib = [np.asarray(m) for r in res["records"] if r["parent"] == parent and r["id"] != target
               and len(m := res["members"].get(r["id"], ())) ]
        if sib:  # the level's other clusters: the markers are z-scored over them too, as in the first run
            context = sib
    else:
        node = _find(tree, target)
        if node is None:
            raise ValueError(f"unknown type or cluster {target!r}")
        if not node["children"]:
            raise ValueError(f"{target} has no subtypes to split into; add subtypes to the marker list first")
        parent, children = target, node["children"]
        idx = np.where(np.isin(out["labels"], descendants(tree)[target]))[0]
    if not len(idx):
        raise ValueError(f"no cells under {target}")
    if len(idx) < MIN_CELLS:  # run() would leave them at the parent level, unannotated
        raise ValueError(f"{target} has {len(idx)} cells, too few to cluster again (at least {MIN_CELLS}): label "
                         "them instead (label_cluster or label_cells)")
    # the clusters this rerun replaces: those of the target's level and below whose cells are all re-clustered now
    D = descendants(tree)
    affected = {None} | set(D) if parent is None else set(D.get(parent, [parent]))  # parents of the clusters concerned
    inside = np.zeros(len(out["labels"]), bool)
    inside[idx] = True
    gone = {r["id"] for r in res["records"] if r["parent"] in affected and
            (r["id"] == target or (r["id"] in out["members"] and inside[np.asarray(out["members"][r["id"]])].all()))}
    if target not in rec:  # a type or the root: its level is computed again whole
        out["nodes"] = [x for x in out["nodes"] if x["node"] not in {(p or "root") for p in affected}]
    for g in gone:
        out["members"].pop(g, None)
    kept = []
    for r in out["records"]:  # clusters partly re-clustered (a deeper level spanning the target): only their other cells
        if r["id"] in gone:
            continue
        m = out["members"].get(r["id"])
        if r["parent"] in affected and m is not None and inside[np.asarray(m)].any():
            m = np.asarray(m)[~inside[np.asarray(m)]]
            out["members"][r["id"]] = m
            r = dict(r, n_cells=int(len(m)))
        kept.append(r)
    out["records"] = kept
    out["labels"][idx] = parent or UNASSIGNED
    A = _Annotator(counts, panel, progress, seed, knn)
    A.total = _n_levels(children) if children else 1
    A.run(out, parent, idx, children, resolution, prefix, context)
    out["not_found"] = A.not_found(out, tree)
    A.say(1.0, "re-clustering done")
    return out


def _find(tree, name):
    for n in tree:
        if n["name"] == name:
            return n
        f = _find(n["children"], name)
        if f:
            return f
    return None
