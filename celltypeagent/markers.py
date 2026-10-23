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
def parse(text, panel):
    """Return {"ok", "errors", "warnings", "tree", "types", "parent", "top", "genes", "canonical"}."""
    upper = {g.upper(): g for g in panel}
    errors, warnings, roots, stack, names = [], [], [], [], {}

    def gene(name, ln):
        g = upper.get(name.upper())
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
            if all(genes):
                markers.append({"genes": genes, "sign": -1 if m.group(2) == "-" else 1})
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


def _graph(Xn):
    """PCA on standardized genes and a kNN graph with adaptive Gaussian weights; built once per level.
    On a GPU (CUDA or Apple MPS) neighbors are exact (blockwise distances); on a CPU they come from pynndescent."""
    import torch, igraph
    from .config import device
    dev = device()
    X = torch.tensor(Xn.toarray(), dtype=torch.float32, device=dev)
    X = (X - X.mean(0)) / (X.std(0) + 1e-3)
    cov = (X.T @ X / len(X)).cpu()  # 480 x 480: eigh on the CPU works everywhere (MPS has no eigh)
    Z = X @ torch.linalg.eigh(cov)[1][:, -N_PCS:].to(dev)
    del X
    k = min(N_NEIGHBORS, len(Z) - 1)
    n = len(Z)
    if dev == "cpu" and n > 20000:
        from pynndescent import NNDescent
        nbr, dist = NNDescent(Z.numpy(), n_neighbors=k + 1, random_state=0).neighbor_graph
        nbr, dist = nbr[:, 1:], dist[:, 1:]
    else:
        chunk = int(np.clip(1.2e9 / n, 256, 8192))  # rows per distance block, about 5 GB at most
        dist, nbr = [], []
        for i in range(0, n, chunk):
            d, j = torch.topk(torch.cdist(Z[i:i + chunk], Z), k + 1, largest=False)
            dist.append(d[:, 1:].cpu()); nbr.append(j[:, 1:].cpu())
        dist, nbr = torch.cat(dist).numpy(), torch.cat(nbr).numpy()
    w = np.exp(-(dist / np.maximum(dist[:, -1], 1e-6)[:, None]) ** 2).ravel()
    A = sp.csr_matrix((w, (np.repeat(np.arange(n), k), nbr.ravel())), shape=(n, n))
    A = A.maximum(A.T).tocoo()
    m = A.row < A.col
    g = igraph.Graph(n=n, edges=np.c_[A.row[m], A.col[m]])
    g.es["weight"] = A.data[m]
    return g, nbr


def _leiden(g, nbr, resolution, seed):
    """igraph's built-in Leiden (modularity objective; same algorithm as leidenalg, several times faster).
    Clusters smaller than MIN_CLUSTER are merged into the cluster most of their neighbors belong to."""
    import random
    random.seed(seed)
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

    def __init__(self, counts, panel, progress=None, seed=0):
        from .model import lognorm
        counts = sp.csr_matrix(counts)
        self.Xn, self.det = lognorm(counts), (counts > 0).astype(np.float32).tocsr()
        self.gi = {g: i for i, g in enumerate(panel)}
        self.say = progress or (lambda f, m: None)
        self.seed, self.n = seed, counts.shape[0]
        self.done, self.total = 0, 1

    def score(self, children, frac):
        S = np.zeros((frac.shape[0], len(children)))
        for j, ch in enumerate(children):
            parts = []
            for m in ch["markers"]:
                z = _marker_z(frac, [self.gi[g] for g in m["genes"]])
                parts.append(m["sign"] * (z.max(1) if m["sign"] > 0 else z.min(1)))
            S[:, j] = np.mean(parts, axis=0)
        return S

    def _evaluate(self, idx, lab, children):
        k = int(lab.max()) + 1
        onehot = sp.csr_matrix((np.ones(len(idx)), (lab, np.arange(len(idx)))), shape=(k, len(idx)))
        size = np.asarray(onehot.sum(1)).ravel()
        frac = np.asarray((onehot @ self.det[idx]).todense()) / np.maximum(size, 1)[:, None]
        S = self.score(children, frac)
        order = np.argsort(-S, 1)
        best = S[np.arange(k), order[:, 0]]
        second = S[np.arange(k), order[:, 1]] if len(children) > 1 else np.full(k, -np.inf)
        ok = (best >= MIN_SCORE) & (best - second >= MIN_MARGIN)
        return {"k": k, "size": size, "S": S, "order": order, "best": best, "second": second, "ok": ok,
                "assigned_frac": float(size[ok].sum() / len(idx))}

    def cluster(self, parent, idx, children, resolution):
        """Leiden at one level. With resolution "auto", every candidate in RESOLUTIONS is tried on the same graph and
        the one assigning the largest share of cells to a child type wins (ties within 1% go to the lower one)."""
        label = parent or "all cells"
        self.say(self.done / self.total, f"building the neighbor graph of {label} ({len(idx):,} cells)")
        g, nbr = _graph(self.Xn[idx])
        cands = RESOLUTIONS if resolution == AUTO else (float(resolution),)
        best, tried = None, {}
        for r in cands:
            self.say(self.done / self.total, f"Leiden on {label}, resolution {r:g}")
            lab = _leiden(g, nbr, r, self.seed)
            ev = self._evaluate(idx, lab, children)
            tried[f"{r:g}"] = round(ev["assigned_frac"], 3)
            if best is None or ev["assigned_frac"] > best[1]["assigned_frac"] + 0.01:
                best = (r, ev, lab)
        return best[0], best[1], best[2], tried

    def run(self, res, parent, idx, children, resolution, prefix=""):
        """Annotate idx under `parent` with its `children` types, then recurse. Updates res in place."""
        if not children or len(idx) < MIN_CELLS:
            return
        r, ev, lab, tried = self.cluster(parent, idx, children, resolution)
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


def annotate(counts, panel, tree, progress=None, seed=0, resolution=AUTO):
    """Full top-down annotation. Returns {"labels" (type name per cell), "cluster_of" (deepest cluster id per cell),
    "members" (cluster id -> cell rows), "records" (one per cluster), "nodes" (resolution chosen per level)}."""
    A = _Annotator(counts, panel, progress, seed)
    A.total = _n_levels(tree)
    res = {"labels": np.full(A.n, UNASSIGNED, dtype=object), "cluster_of": np.full(A.n, "", dtype=object),
           "members": {}, "records": [], "nodes": []}
    A.run(res, None, np.arange(A.n), tree, resolution)
    A.say(1.0, "marker annotation done")
    return res


def rerun(res, counts, panel, tree, target, resolution=AUTO, prefix="", progress=None, seed=0):
    """Re-annotate part of a result: target "root", a type name (its cells are split again into its subtypes) or a
    cluster id (only that cluster's cells are re-clustered at its level). Returns a new result; res is not modified."""
    out = {"labels": res["labels"].copy(), "cluster_of": res["cluster_of"].copy(), "members": dict(res["members"]),
           "records": list(res["records"]), "nodes": list(res["nodes"])}
    rec = {r["id"]: r for r in res["records"]}
    if target == "root":
        parent, children = None, tree
        idx = np.arange(len(out["labels"]))
    elif target in rec:
        parent = rec[target]["parent"]
        node = _find(tree, parent) if parent else None
        children = node["children"] if node else tree
        idx = out["members"][target]
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
    out["labels"][idx] = parent or UNASSIGNED
    A = _Annotator(counts, panel, progress, seed)
    A.total = _n_levels(children) if children else 1
    A.run(out, parent, idx, children, resolution, prefix)
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
