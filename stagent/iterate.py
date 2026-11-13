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
import json, time
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
    got, and the cells the marker step could not type (typed by their expression)."""
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
    X = markerqc._det(s).tocsc()
    gi = {g: i for i, g in enumerate(s.panel)}
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
               "typed_by_expression_as": pd.Series(lab[untyped]).value_counts().head(8).to_dict(),
               "clusters": [{"n_cells": x["n_cells"], "best_scores": dict(sorted(x["scores"].items(),
                                                                                 key=lambda kv: -kv[1])[:3]),
                             "median_transcripts": int(np.median(tc[np.asarray(ml["members"][x["id"]])]))}
                            for x in sorted(ml["records"], key=lambda x: -x["n_cells"])
                            if x["parent"] is None and not x["assigned"] and x["id"] in ml["members"]][:6]},
           "note": ("fold = how many times as often a type's own positive markers are detected in its cells as in "
                    "the other cells of its lineage; marker_flags: markers that are weak, not specific (siblings "
                    "detect them about as often) or negative markers detected; untyped_by_markers: cells no "
                    "cluster of the marker step went to, which the expression model typed by their expression")}
    return out


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
def _node(tree, name):
    nd = markers._find(tree, name)
    if nd is None:
        raise ValueError(f"no type {name!r} in the marker list")
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
            gone = {g for m in _tokens(e.get("marker")) for g in m["genes"]}
            keep = []
            for m in nd["markers"]:
                genes = [g for g in m["genes"] if g not in gone]
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


def add_message(s, who, text):
    st, rd = current(s)
    rd["messages"].append({"from": who, "text": str(text)[:4000], "time": time.strftime("%Y-%m-%dT%H:%M:%S")})
    save(s, st)
    return rd


def propose(s, summary, edits, organ=None, questions=None):
    """The agent's proposal for the current round: a summary, the edits (each with an id, a reason and the evidence)
    and questions for the user; the edits are checked by applying them all to the list."""
    st, rd = current(s)
    clean = []
    for i, e in enumerate(edits or []):
        if e.get("kind") not in EDIT_KINDS:
            raise ValueError(f"edit {i + 1}: kind must be one of {', '.join(EDIT_KINDS)}")
        clean.append({"id": f"E{i + 1}", **{k: e[k] for k in e if k != "id"}})
    text, problems = apply_edits(rd["marker_text"] or "", clean, s.panel) if clean else (rd["marker_text"], [])
    rd["proposal"] = {"summary": str(summary or "")[:4000], "edits": clean, "questions": list(questions or [])[:5],
                      "organ": organ, "check": problems, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    if organ and not st.get("organ"):
        st["organ"] = organ
    save(s, st)
    return rd["proposal"]


def accept(s, edit_ids):
    """The marker list of the next round: the current one with the chosen edits applied."""
    st, rd = current(s)
    if not rd.get("proposal"):
        raise ValueError("the agent has not proposed anything in this round yet")
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
    rd["decision"] = decision
    save(s, st)
    return rd


def carry(src, dst):
    """The rounds so far go with the annotation made from the accepted edits (the chain of rounds)."""
    st = load(src)
    if st["rounds"]:
        save(dst, st)
