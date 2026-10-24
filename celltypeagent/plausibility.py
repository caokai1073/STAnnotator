"""Spatial plausibility check: cells whose label does not fit where they lie.

A cell's same-lineage fraction is the share of its 15 nearest neighbors in its own lineage. Epithelium and muscle
form sheets; immune cells are scattered or cluster only in places (follicles, infiltrates), so being alone is normal
for them. A cell type counts as cohesive when its lineage and the type itself have a median same-lineage fraction of
at least COHESIVE and at most RARE of its cells are isolated (fraction below ISOLATED). Isolated cells of cohesive
types closer than LINK_UM to each other (same lineage) form a group. Groups whose cells the model is unsure about
are likely mislabels; confident isolated cells may be real (shed or displaced cells, deep crypts, fistulas) and need
a look at the H&E. The suggested type is the best type outside the cell's lineage by expression, weighted by the
lineages around it. Nothing is changed here: corrections go through Session.label_cells.

On three SCP2927 sections (dataset lineage split as truth, never shown to the app) 48-68% of "likely mislabel" cells
and 19-45% of confident isolated cells were in another lineage, against 1.6-6.5% of all cells of the checked types;
the suggested lineage was right for 70-78% of the wrong ones.
"""
import numpy as np, pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from . import markers

COHESIVE = 0.5   # minimum median same-lineage neighbor fraction of a cohesive lineage and type
ISOLATED = 0.2   # a cell is isolated below this fraction
RARE = 0.1       # at most this share of a cohesive type's cells is isolated
LOW_CONF = 0.6   # isolated cells below this confidence are likely mislabels
LINK_UM = 30.0   # isolated cells of one lineage closer than this form one group
LIKELY = "likely mislabel"
CHECK = "isolated, model confident"


def check(s, version=None):
    """Per-cell flags and ranked groups for one version (cached on the session)."""
    r = s.result(version)
    v = r["state"]["version"]
    cache = s.__dict__.setdefault("_plausibility", {})
    if v in cache:
        return cache[v]
    lut = np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)
    names, code = np.unique(lut[r["pred"]].astype(str), return_inverse=True)
    nn = s.knn()
    same = (code[nn] == code[:, None]).mean(1)
    types = np.array(s.types, dtype=object)
    by_lin = pd.Series(same).groupby(code).median()
    by_type = pd.DataFrame({"t": r["pred"], "same": same, "iso": same < ISOLATED}).groupby("t").agg(
        same=("same", "median"), iso=("iso", "mean"), lin=("t", lambda x: names.tolist().index(lut[x.iloc[0]])))
    ok = (by_type.same >= COHESIVE) & (by_type.iso <= RARE) & (by_lin.reindex(by_type.lin).to_numpy() >= COHESIVE)
    cohesive = [int(t) for t in by_type.index[ok] if lut[t] != markers.UNASSIGNED]
    cohesion = {str(n): float(by_lin[i]) for i, n in enumerate(names)}
    locked = np.zeros(len(code), bool)
    locked[r["locked_rows"]] = True
    flag = np.isin(r["pred"], cohesive) & (same < ISOLATED) & ~locked
    conf = r["conf"]
    status = np.full(len(code), "", dtype=object)
    status[flag] = np.where(conf[flag] < LOW_CONF, LIKELY, CHECK)

    rows = np.where(flag)[0]
    logp = s.base(r["state"])["logp"][rows].astype(np.float32)
    nl = len(names)
    nbr = np.stack([np.bincount(code[nn[row]], minlength=nl) for row in rows]) / nn.shape[1] if len(rows) else np.zeros((0, nl))
    pos = {n: i for i, n in enumerate(names)}
    type_lin = np.array([pos.get(l, -1) for l in lut])
    # best type outside the cell's lineage: expression log-probability + log of that lineage's share around the cell
    prior = np.log(np.where(type_lin >= 0, nbr[:, np.maximum(type_lin, 0)], 0.0) + 0.05)
    banned = (type_lin[None, :] == code[rows][:, None]) | (types == markers.UNASSIGNED)[None, :]
    score = np.where(banned, -np.inf, logp + prior)
    alt = types[score.argmax(1)] if len(rows) else np.zeros(0, dtype=object)
    other = nbr.copy()
    other[np.arange(len(rows)), code[rows]] = 0
    around = names[other.argmax(1)] if len(rows) else np.zeros(0, dtype=object)

    xy = s.cells[["x", "y"]].to_numpy()
    groups = []
    if len(rows):
        pairs = np.array(sorted(cKDTree(xy[rows]).query_pairs(LINK_UM)), dtype=int).reshape(-1, 2)
        pairs = pairs[code[rows[pairs[:, 0]]] == code[rows[pairs[:, 1]]]]
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(rows), len(rows)))
        _, comp = connected_components(g, directed=False)
        reg, px, py = r["region"], s.cells.px.to_numpy(), s.cells.py.to_numpy()
        depth, cell_id = s.cells.depth_um.to_numpy(), s.cells.cell.to_numpy()
        for c in np.unique(comp):
            k = np.where(comp == c)[0]
            gr = rows[k]
            low = conf[gr] < LOW_CONF
            groups.append({
                "n_cells": int(len(gr)), "lineage": str(names[code[gr[0]]]),
                "types": pd.Series(types[r["pred"][gr]]).value_counts().to_dict(),
                "surrounded_by": pd.Series(around[k]).value_counts().index[0],
                "suggested_type": pd.Series(alt[k]).value_counts().index[0],
                "median_confidence": round(float(np.median(conf[gr])), 3), "n_low_confidence": int(low.sum()),
                "category": LIKELY if low.mean() >= 0.5 else CHECK,
                "regions": pd.Series(reg[gr]).value_counts().to_dict(),
                "median_depth_um": None if np.isnan(depth[gr]).all() else round(float(np.nanmedian(depth[gr])), 1),
                "center_px": [float(px[gr].mean()), float(py[gr].mean())],
                "bbox_px": [float(px[gr].min()), float(py[gr].min()), float(px[gr].max()), float(py[gr].max())],
                "rows": gr, "cells": cell_id[gr].tolist()})
        groups.sort(key=lambda x: (x["category"] != LIKELY, -x["n_cells"] * (0.5 + x["n_low_confidence"] / x["n_cells"])))
        for i, x in enumerate(groups):
            x["id"] = f"P{i + 1}"
    out = {"version": v, "status": status, "groups": groups, "cohesion": cohesion,
           "cohesive_types": [str(types[t]) for t in cohesive],
           "n_flagged": int(flag.sum()), "n_likely": int((status == LIKELY).sum())}
    cache[v] = out
    return out


def lineage_by_region(s, version=None):
    """For each drawn region, how many cells of every lineage it holds (empty when no regions are drawn)."""
    r = s.result(version)
    if not s.region_names(version):
        return {}
    lut = np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)
    t = pd.crosstab(lut[r["pred"]], r["region"])
    return {reg: {k: int(n) for k, n in t[reg].items() if n} for reg in t.columns}


def summary(s, version=None, top=15):
    """What the agent needs to judge the groups (without cell lists)."""
    res = check(s, version)
    keep = ("id", "n_cells", "lineage", "types", "surrounded_by", "suggested_type", "median_confidence",
            "n_low_confidence", "category", "regions", "median_depth_um")
    return {"version": res["version"], "n_flagged_cells": res["n_flagged"], "n_likely_mislabels": res["n_likely"],
            "n_groups": len(res["groups"]),
            "lineage_cohesion": {k: round(c, 2) for k, c in res["cohesion"].items()},
            "cohesive_types": res["cohesive_types"],
            "groups": [{k: g[k] for k in keep} for g in res["groups"][:top]],
            "lineages_by_region": lineage_by_region(s, version),
            "note": f"only cohesive types are checked (they form sheets: median share of the 15 nearest cells in the "
                    f"same lineage >= {COHESIVE:.0%}, and <= {RARE:.0%} of their cells isolated); isolated = share below "
                    f"{ISOLATED:.0%}; likely mislabel = most cells of the group below confidence {LOW_CONF}; "
                    "suggested_type = best type outside the lineage by expression, weighted by the lineages around"}


def group(s, gid, version=None):
    g = next((x for x in check(s, version)["groups"] if x["id"] == gid), None)
    if g is None:
        raise ValueError(f"no plausibility group {gid} in version {check(s, version)['version']}")
    return g


def label_groups(s, gids, cell_type=None, version=None, text="", author=None):
    """Hand-label the cells of groups found on `version` (default: latest), each with cell_type or, if None, with its
    own suggested type. One new version; refused if the session moved on since the groups were computed."""
    res = check(s, version)
    gs = [group(s, g, res["version"]) for g in gids]
    items = [{"cells": g["cells"], "label": cell_type or g["suggested_type"],
              "text": text or f"plausibility {g['id']}: {g['n_cells']} isolated {g['lineage']} cells", "author": author}
             for g in gs]
    msg = "; ".join(f"{g['id']}: {g['n_cells']} cells -> {it['label']}" for g, it in zip(gs, items))
    return s.label_many(items, f"plausibility check, labeled {msg}", expected_parent=res["version"])
