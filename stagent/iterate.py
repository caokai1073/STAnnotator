"""Refining a marker annotation with the agent, round by round (the app's Refine window).

Each round: the app annotates with the current marker list, report() sums up what may be wrong with the result
(computed, no agent), the agent reviews the whole annotation with the tissue's priors and that report and proposes
edits of the marker list one by one (propose_round), the user accepts some, dismisses the round or stops, and may
write to the agent at any time: the user's words come first. Accepted edits are applied to the list (apply_edits)
and the app annotates again, carrying the regions, priors and corrections; the next report compares with the round
before (changes()).

State: <annotation>/refine.json: {"organ", "rounds": [{"n", "session", "version", "marker_text", "report",
"messages": [{"from": "user" | "agent", "text", "time"}], "proposal": {...} | None, "decision": None | "accepted" |
"dismissed" | "stopped", "accepted": [edit ids]}]}. A round's next annotation is a new annotation (session) of the
same sample; refine.json follows the chain: each new session keeps the file of the round it came from.
"""
import difflib, json, time
import numpy as np, pandas as pd
from . import config, markers
from .osutil import dump_json

WEAK_FOLD = 1.5        # a type whose own positive markers are detected under this many times as often in its cells
                       # as in the other cells of its lineage has weak evidence
LOW_CONF = 0.5         # median confidence of a type under this: the model is unsure about it
BIG_CHANGE = 0.5       # a type whose cells change by this share (or appear / vanish) between rounds is pointed out
EDIT_KINDS = ("add_marker", "remove_marker", "add_type", "remove_type", "rename_type", "move_type", "set_markers")


# ---------------- the computed report of an annotation ----------------
def report(s, version=None):
    """What may be wrong with a marker annotation, computed: the organ, cell counts and confidence per type, types
    with weak marker evidence, markers that do not tell a type from its siblings, types the list has but no cell
    got, and the cells the marker step could not type (typed by their marker genes)."""
    from . import hetissue, markerqc
    if s.scheme["kind"] != "markers":
        raise ValueError("refining works on a marker-based annotation")
    r = s.result(version)
    st = r["state"]
    types = np.array(s.types, dtype=object)
    lab = types[r["pred"]]
    conf = np.asarray(r["conf"], float)
    n = len(lab)
    counts = pd.Series(lab).value_counts()
    med_conf = pd.Series(conf).groupby(lab).median()
    tree = s._marker_tree(st.get("marker_step", {}).get("ops", []))
    lineage = {}
    for top in tree:
        for d in markers.descendants(tree)[top["name"]]:
            lineage[d] = top["name"]
    lin = np.array([lineage.get(t, t) for t in lab], dtype=object)
    X, gi = markerqc.det_columns(s, markerqc.tree_genes(tree))  # (only the marker genes: fast on large samples)
    col = {}

    def det(g):
        if g not in col:
            col[g] = np.asarray(X[:, gi[g]].todense()).ravel() > 0
        return col[g]
    weak, per_type = [], []
    for t, k in counts.items():
        nd = markers._find(tree, t)
        if nd is None:
            continue
        m = lab == t
        rest = (lin == lineage.get(t)) & ~m
        pos = [mk for mk in nd["markers"] if mk["sign"] > 0]
        fin = [float(np.mean(np.any([det(g)[m] for g in mk["genes"] if g in gi], 0))) for mk in pos]
        fout = [float(np.mean(np.any([det(g)[rest] for g in mk["genes"] if g in gi], 0))) if rest.any() else np.nan
                for mk in pos]
        a, b = (float(np.mean(fin)), float(np.nanmean(fout)) if rest.any() else np.nan) if pos else (np.nan, np.nan)
        fold = None if not np.isfinite(b) else round(a / max(b, 1e-3), 2)
        row = {"type": t, "lineage": lineage.get(t), "n_cells": int(k), "share": round(k / n, 4),
               "median_confidence": round(float(med_conf.get(t, np.nan)), 2), "markers_in": round(a, 3),
               "markers_in_rest_of_lineage": None if not np.isfinite(b) else round(b, 3), "fold": fold}
        per_type.append(row)
        if fold is not None and fold < WEAK_FOLD:
            weak.append(row | {"why": f"its markers are detected in {a:.0%} of its cells and {b:.0%} of the other "
                                      f"{lineage.get(t)} cells"})
    qc = markerqc.check(s, version)
    flags = [f for f in qc["flags"] if f.get("flag")]
    ml = s.marker_result(st)
    untyped = np.asarray(ml["labels"]) == markers.UNASSIGNED
    tc = s.cells.total_counts.to_numpy()
    out = {"version": st["version"], "n_cells": n, "organ": hetissue.guess_organ(s),
           "types": sorted(per_type, key=lambda x: -x["n_cells"]),
           "weak_types": weak,
           "low_confidence_types": [x for x in per_type if x["median_confidence"] < LOW_CONF],
           "marker_flags": flags,
           "not_found": s.types_not_found(version),
           "untyped_by_markers": {
               "n_cells": int(untyped.sum()), "share": round(float(untyped.mean()), 4),
               "median_transcripts": int(np.median(tc[untyped])) if untyped.any() else None,
               "median_transcripts_others": int(np.median(tc[~untyped])) if (~untyped).any() else None,
               "typed_by_marker_genes_as": pd.Series(lab[untyped]).value_counts().head(8).to_dict(),
               "clusters": [{"n_cells": x["n_cells"], "best_scores": dict(sorted(x["scores"].items(),
                                                                                 key=lambda kv: -kv[1])[:3]),
                             "median_transcripts": int(np.median(tc[np.asarray(ml["members"][x["id"]])]))}
                            for x in sorted(ml["records"], key=lambda x: -x["n_cells"])
                            if x["parent"] is None and not x["assigned"] and len(ml["members"].get(x["id"], ())) > 0][:6]},
           "note": ("fold = how many times as often a type's own positive markers are detected in its cells as in "
                    "the other cells of its lineage; marker_flags: markers that are weak, not specific (siblings "
                    "detect them about as often) or negative markers detected; untyped_by_markers: cells no "
                    "cluster of the marker step went to, each typed by its marker genes (the found type whose cells most "
                    "often show the same marker genes detected and not detected, with its share of cells; no "
                    "marker gene detected: Unassigned)")}
    out["where_types_lie"] = where_types_lie(s, lab, lineage)
    out["rule_conflicts"] = rule_conflicts(s, version)
    return out


LIE_K = 15  # the nearest cells of another type looked at around each cell


def _lie_majority(s, lab, lineage, k=LIE_K):
    """Per cell: the index (in lineage names) of the lineage more than half of its k nearest cells of another type
    belong to (cells of no type left out), -1 mixed, -2 too few cells of another type around. Returns (that, the
    type names, each cell's type index, the lineage of each type, the lineage names)."""
    from scipy.spatial import cKDTree
    xy = s.cells[["x", "y"]].to_numpy()
    n = len(lab)
    lab = np.asarray(lab, dtype=object)
    names, code = np.unique(lab.astype(str), return_inverse=True)
    lin_of = np.array([lineage.get(t, t) for t in names], dtype=object)
    lnames, lcode = np.unique(lin_of.astype(str), return_inverse=True)
    none = np.isin(names, markers.NO_TYPE)
    maj = np.full(n, -2, np.int32)
    if n < 2:  # (no other cell to lie among)
        return maj, names, code, lin_of, lnames
    tree, kk = cKDTree(xy), min(3 * k + 1, n)
    for a in range(0, n, 200000):  # (in blocks: 46 neighbors of a million cells would take 700 MB)
        nb = tree.query(xy[a:a + 200000], k=kk, workers=-1)[1][:, 1:].astype(np.int32)
        own = code[a:a + 200000]
        other = (code[nb] != own[:, None]) & ~none[code[nb]]
        use = other & (np.cumsum(other, 1) <= k)
        cnt = np.zeros((len(nb), len(lnames)), np.int32)
        np.add.at(cnt, (np.repeat(np.arange(len(nb)), use.sum(1)), lcode[code[nb[use]]]), 1)
        tot = use.sum(1)
        m = np.where(cnt.max(1) * 2 > tot, cnt.argmax(1), -1)
        maj[a:a + 200000] = np.where(tot >= max(3, k // 3), m, -2)
    return maj, names, code, lin_of, lnames


def where_types_lie(s, lab, lineage, k=LIE_K):
    """Per cell type: among which lineage its cells lie. For each cell, its k nearest cells of another type (cells of
    no type left out); the lineage of more than half of them, else "mixed". Per type the share of its cells in each,
    the lineage most of its cells lie among ("mostly_among", when that is over half of them and not its own) and the
    rest. Computed only, the priors not applied: the agent compares it with them (a mislabeled block, e.g. epithelial
    cells called neurons, lies among the epithelium; interstitial cells among smooth muscle are expected)."""
    maj, names, code, lin_of, lnames = _lie_majority(s, lab, lineage, k)
    none = np.isin(names, markers.NO_TYPE)
    rows = []
    for i, t in enumerate(names):
        idx = np.where(code == i)[0]
        if none[i] or len(idx) < 20:
            continue
        v = maj[idx]
        sh = {str(lnames[j]): round(float((v == j).mean()), 3) for j in np.unique(v[v >= 0])}
        own_l = str(lin_of[i])
        among = dict(sorted(((l, x) for l, x in sh.items() if l != own_l), key=lambda kv: -kv[1])[:3])
        top = next(iter(among.items()), (None, 0.0))
        rows.append({"type": str(t), "lineage": own_l, "n_cells": int(len(idx)), "among": among,
                     "among_own_lineage": sh.get(own_l, 0.0), "mixed": round(float((v == -1).mean()), 3),
                     "mostly_among": top[0] if top[1] > 0.5 else None})
    return sorted(rows, key=lambda x: -x["n_cells"])


def _lineages(s, st):
    """{type: its top-level type (lineage)} in the marker list of an annotation state."""
    tree = s._marker_tree(st.get("marker_step", {}).get("ops", []))
    return {d: top["name"] for top in tree for d in markers.descendants(tree)[top["name"]]}


_LIE_CACHE = {}  # (session, version) -> _lie_majority of its labels


def finding_cells(s, f, version=None):
    """The cells a finding against the priors is about (rows), the cells of its type (rows) and a title: the cells of
    the finding's rule (rule_id, a rule of the user's that it goes against), else the cells of its type (type, or the
    longest type name in its cells text; a group: all its types) lying among its lineage (among), else all of them."""
    r = s.result(version)
    st = r["state"]
    types = np.array(s.types, dtype=object)
    lab = types[r["pred"]]
    lineage = _lineages(s, st)
    groups = markers.descendants(s._marker_tree(st.get("marker_step", {}).get("ops", [])))
    names = sorted({str(t) for t in lab if t not in markers.NO_TYPE} | set(groups), key=len, reverse=True)
    t = f.get("type") if f.get("type") in names else next((x for x in names if x.lower() in str(f.get("cells", "")).lower()), None)
    of_type = np.where(np.isin(lab, groups.get(t, [t])))[0] if t else np.zeros(0, int)
    cf = r.get("conflicts")
    if f.get("rule_id") and cf is not None and not cf.empty and (cf.constraint == f["rule_id"]).any():
        g = cf[cf.constraint == f["rule_id"]]
        if t:
            g = g[g.base_label.isin(groups.get(t, [t]))]
        rows = g.row.to_numpy(dtype=int)
        return rows, of_type, f"{t or 'cells'} against your rule {f['rule_id']}: {len(rows):,}"
    if not t:
        return None, None, None
    among = f.get("among")
    if among:
        ck = (s.id, st["version"])  # (several maps of one round: the neighbors computed once)
        if ck not in _LIE_CACHE:
            if len(_LIE_CACHE) > 3:
                _LIE_CACHE.clear()
            _LIE_CACHE[ck] = _lie_majority(s, lab, lineage)
        maj, _, _, _, lnames = _LIE_CACHE[ck]
        j = np.where(lnames == among)[0]
        if len(j):
            rows = of_type[maj[of_type] == j[0]]
            return rows, of_type, f"{t} among {among}: {len(rows):,} of {len(of_type):,}"
    return of_type, of_type, f"{t}: {len(of_type):,} cells"


def rule_conflicts(s, version=None):
    """The user's rules (priors drawn on regions) the annotation goes against: per rule, the cells of the types it
    names that lie where it does not allow them (confident cells), by how they lie (in a cluster, at the region's
    boundary, isolated) and whether the rule changed their label."""
    from . import constraints as cmod
    r = s.result(version)
    cf = r.get("conflicts")
    rules = {c["id"]: c for c in r["state"].get("constraints", [])}
    if cf is None or cf.empty:
        return []
    out = []
    for cid, g in cf.groupby("constraint"):
        c = rules.get(cid)
        out.append({"rule": cid, "text": cmod.describe(c) if c else str(cid), "n_cells": int(len(g)),
                    "by_how": {k: int(v) for k, v in g.category.value_counts().items()},
                    "label_changed": int((g.status == "overridden").sum()),
                    "types": {k: int(v) for k, v in g.base_label.value_counts().head(4).items()}})
    return sorted(out, key=lambda x: -x["n_cells"])


def changes(before, after):
    """Cell counts per type between two reports: the types that appeared, vanished or changed by BIG_CHANGE."""
    a = {x["type"]: x["n_cells"] for x in before["types"]}
    b = {x["type"]: x["n_cells"] for x in after["types"]}
    out = []
    for t in sorted(set(a) | set(b)):
        x, y = a.get(t, 0), b.get(t, 0)
        if x == y:
            continue
        big = (x == 0) != (y == 0) or abs(y - x) >= BIG_CHANGE * max(x, y)
        out.append({"type": t, "before": x, "after": y, "big": bool(big)})
    return sorted(out, key=lambda c: (not c["big"], -abs(c["after"] - c["before"])))


# ---------------- marker list edits ----------------
def _names(nodes):
    return [x for n in nodes for x in [n["name"], *_names(n.get("children", []))]]


def _node(tree, name):
    nd = markers._find(tree, name)
    if nd is None:
        names = _names(tree)
        near = difflib.get_close_matches(name, names, n=3, cutoff=0.5)
        raise ValueError(f"no type {name!r} in the marker list" + (f"; did you mean {', '.join(map(repr, near))}?" if near else ""))
    return nd


def _parent_list(tree, name):
    """The list of nodes that holds name (the tree itself at the top)."""
    if any(n["name"] == name for n in tree):
        return tree
    for n in tree:
        try:
            return _parent_list(n["children"], name)
        except ValueError:
            pass
    raise ValueError(f"no type {name!r} in the marker list")


def _tokens(text):
    """Markers written as in the list ("A|B+ C-") -> marker dicts."""
    out = []
    for tok in str(text or "").replace(",", " ").split():
        sign = -1 if tok.endswith("-") else 1
        genes = tok.rstrip("+-").split("|")
        if genes and all(genes):
            out.append({"genes": genes, "sign": sign})
    return out


def apply_edits(text, edits, panel=None):
    """The marker list text with the edits applied in order; each edit {"kind", "type", ...}:
    add_marker / remove_marker {type, marker ("GENE+", "GENE-" or "A|B+")}; set_markers {type, markers};
    add_type {type, parent (None: top level), markers}; remove_type {type} (its subtypes too);
    rename_type {type, new_name}; move_type {type, parent}. Returns (text, problems of the result)."""
    import copy
    tree = copy.deepcopy(markers.parse(text, None)["tree"])
    for e in edits:
        k, t = e.get("kind"), e.get("type")
        if k == "add_marker":
            nd = _node(tree, t)
            for m in _tokens(e.get("marker")):
                if not any(set(x["genes"]) == set(m["genes"]) and x["sign"] == m["sign"] for x in nd["markers"]):
                    nd["markers"].append(m)
        elif k == "remove_marker":
            nd = _node(tree, t)
            gone = {(g, m["sign"]) for m in _tokens(e.get("marker")) for g in m["genes"]}
            keep = []
            for m in nd["markers"]:
                genes = [g for g in m["genes"] if (g, m["sign"]) not in gone]
                if genes:
                    keep.append(dict(m, genes=genes))
            nd["markers"] = keep
        elif k == "set_markers":
            _node(tree, t)["markers"] = _tokens(e.get("markers"))
        elif k == "add_type":
            new = {"name": t, "markers": _tokens(e.get("markers")), "children": []}
            (tree if not e.get("parent") else _node(tree, e["parent"])["children"]).append(new)
        elif k == "remove_type":
            lst = _parent_list(tree, t)
            lst[:] = [n for n in lst if n["name"] != t]
        elif k == "rename_type":
            _node(tree, t)["name"] = e["new_name"]
        elif k == "move_type":
            lst = _parent_list(tree, t)
            nd = next(n for n in lst if n["name"] == t)
            lst[:] = [n for n in lst if n["name"] != t]
            (tree if not e.get("parent") else _node(tree, e["parent"])["children"]).append(nd)
        else:
            raise ValueError(f"unknown edit {k!r}: one of {', '.join(EDIT_KINDS)}")
    out = markers.to_text(tree) + "\n"
    p = markers.parse(out, panel)
    return out, (p["errors"] if not p["ok"] else []) + p.get("warnings", [])


# ---------------- the rounds ----------------
def _file(s):
    return s.dir / "refine.json"


def load(s):
    f = _file(s)
    return json.load(open(f, encoding="utf-8")) if f.exists() else {"organ": None, "rounds": []}


def save(s, st):
    dump_json(st, _file(s), indent=1, default=str)


def start_round(s, organ=None):
    """A new round on the annotation's latest version: its report, compared with the round before (if any)."""
    st = load(s)
    if organ:
        st["organ"] = organ
    rep = report(s)
    if st["rounds"] and st["rounds"][-1].get("report"):
        rep["changes"] = changes(st["rounds"][-1]["report"], rep)
    st["rounds"].append({"n": len(st["rounds"]) + 1, "session": s.id, "version": rep["version"],
                         "marker_text": s.scheme.get("text"), "report": rep, "messages": [], "proposal": None,
                         "decision": None, "accepted": [], "time": time.strftime("%Y-%m-%dT%H:%M:%S")})
    save(s, st)
    return st


def current(s):
    st = load(s)
    if not st["rounds"]:
        raise ValueError("no refine round yet: start one")
    return st, st["rounds"][-1]


def mark_asked(s):
    """The agent was asked to review the current round (Agent review, or a message of the user). A dismissed round is
    open again (the user asks for a new proposal); an accepted or stopped one stays decided."""
    st, rd = current(s)
    rd["asked"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    rd["asked_t"] = time.time()
    if rd.get("decision") == "dismissed":
        rd["decision"], rd["accepted"] = None, []
    save(s, st)


def awaiting(s):
    """Whether the current round waits for the agent's proposal: asked, not decided, nothing proposed since."""
    try:
        st, rd = current(s)
    except ValueError:
        return False
    p = rd.get("proposal") or {}
    if not rd.get("asked") or rd.get("decision"):
        return False
    if "asked_t" in rd:
        return float(p.get("t") or 0) <= float(rd["asked_t"])
    return (p.get("time") or "") < rd["asked"]


def edits_from_text(old, new):
    """A whole revised marker list as edits of the old one, type by type (the Refine window shows edits): add_type
    (parents first), move_type, set_markers where a type's markers changed, remove_type (only the topmost removed)."""
    pn = markers.parse(new or "", None)
    if not pn["ok"]:
        raise ValueError("the revised list has errors: " + "; ".join(pn["errors"][:5]))

    def flat(text):
        out = {}
        def walk(nodes, parent):
            for n in nodes:
                out[n["name"]] = (parent, n["markers"])
                walk(n["children"], n["name"])
        walk(markers.parse(text, None)["tree"], None)
        return out
    def mtext(ms):
        return " ".join("|".join(m["genes"]) + ("-" if m["sign"] < 0 else "+") for m in ms)
    def key(ms):
        return {(frozenset(m["genes"]), m["sign"]) for m in ms}
    a, b = flat(old or ""), flat(new or "")
    adds, moves, sets, removes = [], [], [], []
    for t, (par, ms) in b.items():  # in the new list's order: a parent before its children
        if t not in a:
            adds.append({"kind": "add_type", "type": t, "parent": par, "markers": mtext(ms)})
            continue
        if a[t][0] != par:
            moves.append({"kind": "move_type", "type": t, "parent": par})
        if key(a[t][1]) != key(ms):
            gone = key(a[t][1]) - key(ms)
            new_ = key(ms) - key(a[t][1])
            fmt = lambda ks: ", ".join("|".join(sorted(g)) + ("-" if s_ < 0 else "+") for g, s_ in sorted(ks, key=str))
            sets.append({"kind": "set_markers", "type": t, "markers": mtext(ms),
                         "change": "; ".join(x for x in (f"added {fmt(new_)}" if new_ else "",
                                                         f"removed {fmt(gone)}" if gone else "") if x)})
    for t, (par, _) in a.items():
        if t not in b and (par is None or par in b):
            removes.append({"kind": "remove_type", "type": t})
    return adds + moves + sets + removes


def add_message(s, who, text):
    st, rd = current(s)
    rd["messages"].append({"from": who, "text": str(text)[:4000], "time": time.strftime("%Y-%m-%dT%H:%M:%S")})
    save(s, st)
    return rd


def _group(edits, text):
    """Each edit's "group", the heading it is shown under in the Refine window: the problem it fixes as the agent named
    it ("problem"), else the lineage (top-level type) it changes in the round's list; groups keep the order in which
    they first come. The edits keep their order (they are applied in it: a rename before the edits of the new name)."""
    try:
        tree = markers.parse(text, None)["tree"]
    except Exception:
        tree = []
    top = {d: n["name"] for n in tree for d in markers.descendants(tree)[n["name"]]}
    renamed = {}
    for e in edits:
        prob = " ".join(str(e.get("problem") or "").split())[:160]
        if prob:
            e["problem"] = prob
        else:
            e.pop("problem", None)  # (a blank one: the edit goes under its lineage, on the page too)
        t = e.get("parent") if e.get("kind") == "add_type" and e.get("parent") else e.get("type")
        lin = top.get(t) or top.get(renamed.get(t)) or (e.get("type") if e.get("kind") == "add_type" else None)
        if e.get("kind") == "rename_type" and e.get("new_name"):
            renamed[e["new_name"]] = e.get("type")
        e["group"] = prob or (f"{lin}" if lin else "Other edits")
    return edits


JUDGEMENTS = ("likely mislabel", "may be real biology", "unsure")
VERDICTS = ("biology", "mislabel")


def finding_key(f):
    """The same finding in another round: its cells, rule, type and lineage, case and spacing ignored."""
    return " ".join("|".join(str(f.get(k) or "") for k in ("cells", "rule", "type", "among")).lower().split())


def _findings(found):
    """The agent's findings against the priors, cleaned: {key, cells, where, rule, judgement, evidence}."""
    out = []
    for i, f in enumerate(found or []):
        if not isinstance(f, dict) or not str(f.get("cells") or "").strip():
            raise ValueError(f"prior finding {i + 1}: give at least the cells it is about (cells)")
        g = {k: " ".join(str(f.get(k) or "").split())[:600] for k in ("cells", "where", "rule", "evidence")}
        g |= {k: " ".join(str(f[k]).split())[:200] for k in ("type", "among", "rule_id") if f.get(k)}
        j = str(f.get("judgement") or "unsure").strip().lower()
        g["judgement"] = j if j in JUDGEMENTS else "unsure"
        g["key"] = finding_key(g)
        out.append(g)
    return out


def set_verdict(s, key, verdict):
    """The user's verdict on a finding against the priors: "biology" (real: leave it), "mislabel" (fix the list) or
    None (cleared). Kept for all later rounds; the agent reads them (refine_report) and does not raise a finding the
    user called real biology again."""
    st, rd = current(s)
    v = st.setdefault("prior_verdicts", {})
    # (a finding of an earlier proposal, or one the user has given a verdict before: its proposal may be replaced)
    found = {k: x for k, x in v.items()} | {f["key"]: f for r in st["rounds"] for f in ((r.get("proposal") or {}).get("findings") or [])}
    if key not in found:
        raise ValueError("no such finding in this refinement")
    if verdict is None:
        v.pop(key, None)
    elif verdict in VERDICTS:
        f = found[key]
        v[key] = {"verdict": verdict, "cells": f["cells"], "where": f["where"], "rule": f["rule"], "round": rd["n"],
                  "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    else:
        raise ValueError(f"verdict: {' or '.join(VERDICTS)}")
    save(s, st)
    return v


def propose(s, summary, edits, organ=None, questions=None, findings=None):
    """The agent's proposal for the current round: a summary, the edits (each with an id, a reason and the evidence),
    questions for the user and findings against the priors (each the cells, where they lie, the rule, the agent's
    judgement and the evidence); the edits are checked by applying them all to the list."""
    st, rd = current(s)
    if rd.get("decision") in ("accepted", "stopped"):
        raise ValueError(f"round {rd['n']} is already {rd['decision']}")
    rd["decision"], rd["accepted"] = None, []  # (a dismissed round gets a new proposal)
    clean = []
    for i, e in enumerate(edits or []):
        if e.get("kind") not in EDIT_KINDS:
            raise ValueError(f"edit {i + 1}: kind must be one of {', '.join(EDIT_KINDS)}")
        clean.append({"id": f"E{i + 1}", **{k: e[k] for k in e if k != "id"}})
        o = e.get("optional")  # (shown unticked: the user chooses it; "false" written as a word is false)
        clean[-1]["optional"] = o.strip().lower() in ("true", "yes", "1") if isinstance(o, str) else bool(o)
        if e.get("question"):
            clean[-1]["question"] = " ".join(str(e["question"]).split())[:600]
    _group(clean, rd["marker_text"] or "")
    text, problems = apply_edits(rd["marker_text"] or "", clean, s.panel) if clean else (rd["marker_text"], [])
    found = _findings(findings)
    rd["proposal"] = {"summary": str(summary or "")[:6000], "edits": clean, "questions": list(questions or [])[:10],
                      "findings": found,
                      "organ": organ, "check": problems, "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
                      "t": time.time()}
    if organ and not st.get("organ"):
        st["organ"] = organ
    text_q = "".join(f"\n? {q}" for q in rd["proposal"]["questions"])
    rd["messages"].append({"from": "agent", "text": (rd["proposal"]["summary"] + text_q).strip(),
                           "time": rd["proposal"]["time"]})
    save(s, st)
    return rd["proposal"]


def accept(s, edit_ids, proposal_time=None):
    """The marker list of the next round: the current one with the chosen edits applied. proposal_time: the time of
    the proposal the user saw (refused when the agent has proposed again since)."""
    st, rd = current(s)
    if rd.get("decision") in ("accepted", "stopped"):
        raise ValueError(f"round {rd['n']} is already {rd['decision']}")
    if not rd.get("proposal"):
        raise ValueError("the agent has not proposed anything in this round yet")
    if proposal_time and proposal_time != rd["proposal"].get("time"):
        raise ValueError("the agent has proposed again since: look at the new edits and choose again")
    chosen = [e for e in rd["proposal"]["edits"] if e["id"] in set(edit_ids)]
    if not chosen:
        raise ValueError("choose at least one edit")
    text, problems = apply_edits(rd["marker_text"] or "", chosen, s.panel)
    errors = [p for p in problems if "error" in str(p).lower()]
    p = markers.parse(text, s.panel)
    if not p["ok"]:
        raise ValueError("the edited list has errors: " + "; ".join(p["errors"][:5]))
    rd["decision"], rd["accepted"] = "accepted", [e["id"] for e in chosen]
    save(s, st)
    return {"marker_text": text, "warnings": problems, "errors": errors}


def decide(s, decision):
    st, rd = current(s)
    if decision not in ("dismissed", "stopped"):
        raise ValueError("decision: dismissed or stopped")
    if rd.get("decision") == "accepted":
        raise ValueError(f"round {rd['n']} is already accepted")
    rd["decision"] = decision
    save(s, st)
    return rd


def carry(src, dst):
    """The rounds so far go with the annotation made from the accepted edits (the chain of rounds); only between
    annotations of the same sample."""
    if src.sample_id != dst.sample_id:
        return
    st = load(src)
    if st["rounds"]:
        save(dst, st)
