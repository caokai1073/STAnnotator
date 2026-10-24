"""Marker check of an annotation: how well the cells of every node of the marker tree support its markers.

For a node (a lineage, a group or a cell type) its cells are all cells whose label is one of its cell types; it is
compared with its siblings (the other children of its parent; the other lineages at the top). For each listed
marker the share of the node's cells in which it is detected (any gene of an A|B marker) is set against the
siblings: a positive marker is weak when detected in fewer than WEAK of the node's cells, not specific when the
siblings detect it at UNSPECIFIC or more of that share; a negative marker is contradicted when detected in NEG_HIGH
or more of the node's cells. The genes that best separate the node from its siblings in this sample are listed too,
so the agent can judge with biology whether the label fits, propose better markers from the panel, and recognize
cell types missing from the list (Unassigned cells are compared with all assigned cells).
"""
import numpy as np, scipy.sparse as sp
from . import markers

WEAK = 0.25
UNSPECIFIC = 0.7
NEG_HIGH = 0.3
MIN_CELLS = 20
TOP = 8


def _det(s):
    if not hasattr(s, "_det_X"):
        X = s.counts.tocsr()
        s._det_X = sp.csr_matrix((np.ones(X.nnz, np.float32), X.indices, X.indptr), shape=X.shape)
    return s._det_X


def _marker_text(m):
    return "|".join(m["genes"]) + ("-" if m["sign"] < 0 else "+")


def check(s, version=None):
    """Per node of the marker tree: marker support, flags and top distinguishing genes (cached per version)."""
    r = s.result(version)
    v = r["state"]["version"]
    cache = s.__dict__.setdefault("_markerqc", {})
    if v in cache:
        return cache[v]
    if s.scheme["kind"] != "markers":
        raise ValueError("the marker check needs a marker-based annotation")
    Xb = _det(s)
    lab = r["pred"]
    n, T = len(lab), len(s.types)
    oh = sp.csr_matrix((np.ones(n, np.float32), (lab, np.arange(n))), shape=(T, n))
    det = np.asarray((oh @ Xb).todense())                      # cells of each type with each gene detected
    nt = np.bincount(lab, minlength=T).astype(float)
    gi = {g: i for i, g in enumerate(s.panel)}
    leaves = {k: [s.types.index(t) for t in v] for k, v in s.type_groups.items()}
    for t in s.types:
        leaves.setdefault(t, [s.types.index(t)])
    any_cache = {}

    def any_rate(genes, idx):
        """Share of cells of the types idx in which any of genes is detected."""
        key = tuple(genes)
        if key not in any_cache:
            hit = np.asarray(Xb[:, [gi[g] for g in genes if g in gi]].sum(1)).ravel() > 0
            any_cache[key] = np.bincount(lab, weights=hit, minlength=T)
        k = nt[idx].sum()
        return float(any_cache[key][idx].sum() / k) if k else float("nan")

    def top_genes(idx, other):
        a, b = nt[idx].sum(), nt[other].sum()
        if not a or not b:
            return []
        ra, rb = det[idx].sum(0) / a, det[other].sum(0) / b
        score = (ra - rb) * np.log2((ra + 0.02) / (rb + 0.02))
        score[ra < 0.15] = -np.inf
        order = [i for i in np.argsort(-score)[:TOP] if np.isfinite(score[i]) and ra[i] > rb[i]]
        return [{"gene": s.panel[i], "in": round(float(ra[i]), 2), "siblings": round(float(rb[i]), 2)} for i in order]

    nodes, flags = [], []

    def walk(level_nodes, depth, parent):
        for nd in level_nodes:
            idx = leaves.get(nd["name"], [])
            sib = [i for o in level_nodes if o["name"] != nd["name"] for i in leaves.get(o["name"], [])]
            n_cells = int(nt[idx].sum())
            rec = {"node": nd["name"], "level": depth, "parent": parent, "is_cell_type": not nd["children"],
                   "n_cells": n_cells, "compared_with": [o["name"] for o in level_nodes if o["name"] != nd["name"]] or ["-"],
                   "markers": [], "top_genes": []}
            if n_cells >= MIN_CELLS:
                for m in nd["markers"]:
                    a, b = any_rate(m["genes"], idx), any_rate(m["genes"], sib) if sib else float("nan")
                    flag = None
                    if m["sign"] > 0 and a < WEAK:
                        flag = "weak: detected in few of these cells"
                    elif m["sign"] > 0 and sib and b >= UNSPECIFIC * a:
                        flag = "not specific: siblings detect it about as often"
                    elif m["sign"] < 0 and a >= NEG_HIGH:
                        flag = "negative marker often detected in these cells"
                    mr = {"marker": _marker_text(m), "in": round(a, 2), "siblings": None if np.isnan(b) else round(b, 2),
                          "flag": flag}
                    rec["markers"].append(mr)
                    if flag:
                        flags.append({"node": nd["name"], **mr})
                rec["top_genes"] = top_genes(idx, sib) if sib else []
            else:
                rec["note"] = f"fewer than {MIN_CELLS} cells"
            nodes.append(rec)
            walk(nd["children"], depth + 1, nd["name"])

    walk(s.scheme["tree"], 1, None)
    ua = s.types.index(markers.UNASSIGNED) if markers.UNASSIGNED in s.types else None
    unassigned = None
    if ua is not None and nt[ua] >= MIN_CELLS:
        rest = [i for i in range(T) if i != ua]
        unassigned = {"n_cells": int(nt[ua]), "top_genes": top_genes([ua], rest)}
    out = {"version": v, "nodes": nodes, "flags": flags, "unassigned": unassigned,
           "thresholds": {"weak_below": WEAK, "not_specific_ratio": UNSPECIFIC, "negative_high": NEG_HIGH},
           "note": "in / siblings = share of cells with the marker (any gene of A|B) detected among the node's cells "
                   "and its siblings' cells; top_genes = genes that best separate the node from its siblings here"}
    cache[v] = out
    return out


# ---------------- marker list proposals (written by Claude, opened in the viewer's Markers dialog) ----------------
def _pfile(s):
    return s.dir / "marker_proposal.json"


def proposal(s):
    import json
    p = _pfile(s)
    return json.load(open(p)) if p.exists() else None


def propose(s, text, summary=""):
    """Validate a revised marker list and store it as the session's open proposal; return its diff."""
    import difflib, json, time
    parsed = markers.parse(text, s.panel)
    if not parsed["ok"]:
        raise ValueError("marker list has errors: " + "; ".join(parsed["errors"][:8]))
    cur = s.scheme.get("text", "")
    diff = [l for l in difflib.unified_diff(cur.splitlines(), parsed["canonical"].splitlines(), lineterm="", n=0)
            if l[:1] in "+-" and not l.startswith(("+++", "---"))]
    rec = {"text": parsed["canonical"], "summary": summary.strip(), "diff": diff, "warnings": parsed["warnings"],
           "n_types": len(parsed["leaves"]), "status": "open", "created_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    json.dump(rec, open(_pfile(s), "w"), indent=1)
    return rec


def set_proposal_status(s, status):
    import json
    if status not in ("open", "opened", "dismissed"):
        raise ValueError("status must be open, opened or dismissed")
    rec = proposal(s)
    if rec is None:
        raise ValueError("no marker list proposal")
    rec["status"] = status
    json.dump(rec, open(_pfile(s), "w"), indent=1)
    return rec
