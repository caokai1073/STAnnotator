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


def _topk_torch(Z, k, dev):
    """Indices of the k nearest rows (the row itself included) with torch on `dev`, blockwise."""
    import torch
    Zt = torch.tensor(Z, device=dev)
    n = len(Z)
    chunk = int(np.clip(1.2e9 / n, 256, 8192))  # rows per distance block, about 5 GB at most
    nbr = [torch.topk(torch.cdist(Zt[i:i + chunk], Zt), k, largest=False)[1].cpu() for i in range(0, n, chunk)]
    del Zt
    return torch.cat(nbr).numpy().astype(np.int64)


def _topk_numpy(Z, k):
    """Indices of the k nearest rows (the row itself included), blockwise through BLAS (a KD-tree is slow in 30
    dimensions)."""
    sq = (Z * Z).sum(1)
    nbr = np.empty((len(Z), k), np.int64)
    for i in range(0, len(Z), 4096):
        D = sq[i:i + 4096, None] + sq[None, :] - 2 * Z[i:i + 4096] @ Z.T
        nbr[i:i + 4096] = np.argpartition(D, k - 1, axis=1)[:, :k]
    return nbr


MAX_PCA_GENES, N_HVG = 6000, 3000  # above MAX_PCA_GENES genes (whole-transcriptome data such as Visium HD) the PCA
                                    # uses the N_HVG most variable ones; a panel (Xenium 5K too) uses all its genes


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


def pca(Xn, n_pcs=None):
    """PCA of standardized genes on the sparse matrix (no dense copy): cell scores (cells x n_pcs, float32)."""
    X = sp.csr_matrix(Xn, dtype=np.float32)
    n = X.shape[0]
    mu = np.asarray(X.mean(0)).ravel()
    ex2 = np.asarray(X.multiply(X).mean(0)).ravel()
    if X.shape[1] > MAX_PCA_GENES:  # the most variable genes (of the log-normalized values)
        keep = np.sort(np.argsort(-(ex2 - mu ** 2))[:N_HVG])
        X, mu, ex2 = X[:, keep], mu[keep], ex2[keep]
    sd = np.sqrt(np.maximum(ex2 - mu ** 2, 0) * n / max(n - 1, 1)) + 1e-3
    cov = ((X.T @ X).toarray() / n - np.outer(mu, mu)) / np.outer(sd, sd)  # covariance of standardized genes
    V = np.linalg.eigh(cov)[1][:, -(n_pcs or N_PCS):] / sd[:, None]  # ascending variance: the last column is PC1
    return np.asarray(X @ V - mu @ V, dtype=np.float32)  # = standardized X times the top components


KNN = "approx"  # neighbor search of the marker step: "approx" (default) or "exact" (see knn)
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
    if has_torch():
        return "knn_exact_torch", "exact, on the CPU (torch)"
    return "knn_exact_numpy", "exact, on the CPU (numpy)"


def _exact(Z, k):
    """Exact neighbors, the same on every computer: candidates from a GPU (CUDA or Apple MPS) when torch has one,
    else from the CPU (torch, or numpy without it), re-ranked by their float64 distance with ties going to the lower
    row number, so rounding differences between the backends do not change the result."""
    from .config import device
    from . import timing
    n = len(Z)
    kk = min(k + 1 + KNN_EXTRA, n)
    key, _ = knn_backend("exact")
    cand = None
    if key == "knn_exact_gpu":
        try:
            with timing.step(key, n * n / 1e9):
                cand = _topk_torch(Z, kk, device())
        except RuntimeError:  # e.g. GPU memory taken by other jobs: the CPU gives the same neighbors
            key = "knn_exact_torch"
    if cand is None:
        with timing.step(key, n * n / 1e9):
            cand = _topk_torch(Z, kk, "cpu") if key == "knn_exact_torch" else _topk_numpy(Z, kk)
    Z64 = Z.astype(np.float64)
    dist, nbr = np.empty((n, k), np.float32), np.empty((n, k), np.int64)
    with timing.step("knn_rerank", n / 1e6):
        for i in range(0, n, 32768):
            c = cand[i:i + 32768]
            rows = np.arange(i, i + len(c))
            d = ((Z64[c] - Z64[rows, None]) ** 2).sum(-1)
            d[c == rows[:, None]] = np.inf  # not its own neighbor
            o = np.lexsort((c, d), axis=1)[:, :k]
            nbr[i:i + len(c)] = np.take_along_axis(c, o, 1)
            dist[i:i + len(c)] = np.sqrt(np.take_along_axis(d, o, 1))
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
    if mode == "exact":
        return _exact(Z, k)
    if mode == "auto" and (device() != "cpu" or n <= 20000):
        return _exact(Z, k)
    with timing.step("knn_approx", n / 1e6, warmup="approx_compile"):
        nbr, dist = _nndescent(Z, k + 1)
    return dist[:, 1:], nbr[:, 1:]


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


def _leiden(g, nbr, resolution, seed):
    """igraph's built-in Leiden (modularity objective; same algorithm as leidenalg, several times faster).
    Clusters smaller than MIN_CLUSTER are merged into the cluster most of their neighbors belong to."""
    import random
    with _LEIDEN_LOCK:  # igraph draws from Python's global random: two runs at once (an annotation and a PCA / UMAP
        random.seed(seed)  # in the app's threads) would take each other's numbers and depend on timing
        lab = np.array(g.community_leiden(objective_function="modularity", weights="weight", resolution=resolution,
                                          n_iterations=2).membership)
    sizes = np.bincount(lab)
    small = sizes[lab] < MIN_CLUSTER
    if small.any() and (~small).any():
        for i in np.where(small)[0]:
            cand = lab[nbr[i]][sizes[lab[nbr[i]]] >= MIN_CLUSTER]
            if len(cand):
                lab[i] = np.bincount(cand).argmax()
    return np.unique(lab, return_inverse=True)[1]


def _marker_z(det_frac, gene_idx):
    f = det_frac[:, gene_idx]
    return np.clip((f - f.mean(0)) / (f.std(0) + 0.02), -3, 3)


class _Annotator:
    """Top-down marker annotation over one sample; also re-runs single levels or clusters in place."""

    def __init__(self, counts, panel, progress=None, seed=0, knn_mode=KNN):
        from .model import lognorm
        from . import timing
        with timing.step("prepare", counts.shape[0] / 1e6):
            counts = sp.csr_matrix(counts)
            self.Xn, self.det = lognorm(counts), (counts > 0).astype(np.float32).tocsr()
        self.knn_mode = knn_mode
        self.gi = {g: i for i, g in enumerate(panel)}
        self.say = progress or (lambda f, m: None)
        self.seed, self.n = seed, counts.shape[0]
        self.done, self.total = 0, 1
        self.fp = {}  # checksums of the root level (see annotate)

    def score(self, children, frac):
        S = np.zeros((frac.shape[0], len(children)))
        for j, ch in enumerate(children):
            parts = []
            for m in ch["markers"]:
                z = _marker_z(frac, [self.gi[g] for g in m["genes"]])
                parts.append(m["sign"] * (z.max(1) if m["sign"] > 0 else z.min(1)))
            S[:, j] = np.mean(parts, axis=0)
        return S

    def _evaluate(self, idx, lab, children, context=None):
        """context: detection rates (clusters x genes) of the level's other clusters, when only some of its cells are
        clustered again (rerun of one cluster): the markers are z-scored over all of them, as in the first run, not
        over the new sub-clusters alone (between near-identical sub-clusters the z-scores are noise)."""
        k = int(lab.max()) + 1
        onehot = sp.csr_matrix((np.ones(len(idx)), (lab, np.arange(len(idx)))), shape=(k, len(idx)))
        size = np.asarray(onehot.sum(1)).ravel()
        frac = np.asarray((onehot @ self.det[idx]).todense()) / np.maximum(size, 1)[:, None]
        S = self.score(children, frac if context is None or not len(context) else np.vstack([frac, context]))[:k]
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
        g, nbr = _graph(self.Xn[idx], self.knn_mode, self.fp if root else None)
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
            res["records"].append({"id": cid, "parent": parent, "n_cells": int(ev["size"][c]), "assigned": assigned,
                                   "resolution": r,
                                   "scores": {ch["name"]: round(float(ev["S"][c, j]), 2) for j, ch in enumerate(children)},
                                   "margin": round(float(ev["best"][c] - ev["second"][c]), 2) if len(children) > 1 else None})
        self.done += 1
        for ch in children:
            self.run(res, ch["name"], idx[res["labels"][idx] == ch["name"]], ch["children"], resolution, prefix)


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
           "records": list(res["records"]), "nodes": list(res["nodes"]), "knn": res.get("knn")}
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
        if sib:
            det = (sp.csr_matrix(counts) > 0).astype(np.float32).tocsr()
            context = np.vstack([np.asarray(det[m].mean(0)).ravel() for m in sib])
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
