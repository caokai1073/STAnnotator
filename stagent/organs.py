"""Marker lists of common organs (data/organs/), pruned to a sample's gene panel.

Each organ file lists its cell types in the marker list format, with more markers than a panel will hold, most
specific first (data/organs/README.md). draft() keeps the types the panel can resolve: a cell type needs two of its
markers in the panel, or one marked decisive (GENE!); it gets its first four. A group (a type with subtypes) keeps
its subtypes when two or more remain; with one left, that one is dropped and the group becomes a cell type, since
its only subtype would take all its cells. The result is a marker list CSV (cell_type, parent, positive_markers,
negative_markers) and a report of what was left out and why, for the user or the agent to complete.
"""
import csv, io, re
from . import config

DIR = config.PACKAGE / "data" / "organs"
SHARED = ("immune", "stroma", "enteric")
N_MARKERS, N_GROUP_MARKERS, N_NEGATIVE = 4, 5, 2


def organs():
    """[{"organ", "label"}] of the built-in lists."""
    out = []
    for f in sorted(DIR.glob("*.txt")):
        if f.stem in SHARED:
            continue
        first = f.read_text().splitlines()[0]
        out.append({"organ": f.stem, "label": first.lstrip("# ").strip() if first.startswith("#") else f.stem})
    return out


def _lines(name, seen=()):
    """The lines of an organ file with its @includes expanded."""
    f = DIR / f"{name}.txt"
    if not f.exists():
        raise ValueError(f"no built-in list for {name!r}; built-in: {', '.join(o['organ'] for o in organs())}")
    out = []
    for raw in f.read_text().splitlines():
        if raw.strip().startswith("@"):
            inc = raw.strip()[1:]
            if inc in seen:
                raise ValueError(f"{name}: @{inc} includes itself")
            out += _lines(inc, seen + (name,))
        elif raw.strip() and not raw.lstrip().startswith("#"):
            out.append(raw)
    return out


def tree(name):
    """The organ's types in file order: [{"name", "parent", "markers": [(gene, decisive)], "negative": [gene]}]."""
    nodes, stack = [], []
    for raw in _lines(name):
        head, _, body = raw.strip().partition(":")
        indent = (len(raw.expandtabs(2)) - len(raw.expandtabs(2).lstrip())) // 2
        if ">" in head:  # "Parent > Child": under a type defined earlier (e.g. in an included file)
            parent, child = (x.strip() for x in head.split(">", 1))
            if parent not in {n["name"] for n in nodes}:
                raise ValueError(f"{name}: no type {parent!r} for {child!r}")
            head, level = child, None
        else:
            head = head.strip()
            del stack[indent:]
            parent, level = (stack[-1] if stack else None), indent
        toks = body.split()
        node = {"name": head, "parent": parent,
                "markers": [(t.rstrip("!"), t.endswith("!")) for t in toks if not t.startswith("-")],
                "negative": [t[1:] for t in toks if t.startswith("-")]}
        nodes.append(node)
        if level is not None:
            stack.append(head)
    return nodes


def draft(organ, panel):
    """The organ's list pruned to the panel: {"csv", "n_types", "kept", "dropped", "merged", "organ"}."""
    upper = {g.upper(): g for g in panel}
    nodes = tree(organ)
    kids = {}
    for n in nodes:
        kids.setdefault(n["parent"], []).append(n)
    dropped, merged = [], []

    def present(n):
        pos, seen = [], set()
        for g, decisive in n["markers"]:
            h = upper.get(g.upper())
            if h and h not in seen:
                seen.add(h)
                pos.append((h, decisive))
        neg = [upper[g.upper()] for g in n["negative"] if g.upper() in upper]
        return pos, neg

    def build(n):
        """The rows of n and its kept subtypes, or [] when the panel cannot resolve it."""
        pos, neg = present(n)
        sub = [build(c) for c in kids.get(n["name"], [])]
        sub = [r for r in sub if r]
        if len(sub) == 1:  # a lone subtype would take every cell of the group: the group is the cell type
            merged.append({"type": sub[0][0]["cell_type"], "into": n["name"]})
            sub = []
        if sub:
            genes = [g for g, _ in pos][:N_GROUP_MARKERS] or [r[0]["positive"][0] for r in sub]
            row = {"cell_type": n["name"], "parent": n["parent"] or "", "positive": genes, "negative": neg[:N_NEGATIVE],
                   "group": True}
            return [row] + [x for r in sub for x in r]
        if len(pos) >= 2 or (pos and pos[0][1]):
            return [{"cell_type": n["name"], "parent": n["parent"] or "", "positive": [g for g, _ in pos][:N_MARKERS],
                     "negative": neg[:N_NEGATIVE], "group": False}]
        dropped.append({"type": n["name"], "in_panel": [g for g, _ in pos],
                        "markers": [g for g, _ in n["markers"]][:6]})
        return []

    rows = [x for n in kids.get(None, []) for x in build(n)]
    kept = {r["cell_type"] for r in rows}
    for r in rows:  # a kept type whose group was dropped moves up
        while r["parent"] and r["parent"] not in kept:
            r["parent"] = next(n["parent"] or "" for n in nodes if n["name"] == r["parent"])
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(["cell_type", "parent", "positive_markers", "negative_markers"])
    for r in rows:
        w.writerow([r["cell_type"], r["parent"], "|".join(r["positive"]) if r["group"] else " ".join(r["positive"]),
                    " ".join(r["negative"])])
    leaves = [r["cell_type"] for r in rows if not r["group"]]
    return {"organ": organ, "csv": buf.getvalue(), "n_types": len(leaves), "kept": [r["cell_type"] for r in rows],
            "dropped": dropped, "merged": merged}


def safe_name(organ):
    return re.sub(r"[^A-Za-z0-9_]+", "_", organ.strip().lower()).strip("_")
