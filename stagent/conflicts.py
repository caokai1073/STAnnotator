"""Conflict detection: cells whose base prediction is confident but violates a constraint.

A conflict is never resolved silently. Status:
  kept        expression evidence outweighed the prior; label unchanged
  overridden  the prior changed the label
Spatial category (same definitions as scripts/04_conflicts.py):
  boundary    within BOUNDARY_UM of a cell where the type is allowed by the constraint
  clustered   >= CLUSTER_FRAC of its K nearest neighbors share the base prediction group
  isolated    otherwise (often tissue edges, debris, segmentation artifacts)
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from .constraints import penalized_cells

CONF, BOUNDARY_UM, K, CLUSTER_FRAC = 0.9, 50.0, 15, 0.3


def find(cells, region, base_pred, base_conf, pred, constraints, types, ctx):
    rows = []
    xy = cells[["x", "y"]].to_numpy()
    tree = cKDTree(xy)
    for c in constraints:
        if c["mode"] == "enrich":
            continue
        ti = [types.index(t) for t in c["types_resolved"]]
        bad = penalized_cells(c, ctx)
        in_group = np.isin(base_pred, ti)
        viol = in_group & bad & (base_conf > CONF)
        if not viol.any():
            continue
        ok = in_group & ~bad  # cells of the group where the constraint allows them
        allowed = ~bad & ~np.isin(region, ["unknown"])
        d_allowed = cKDTree(xy[allowed]).query(xy[viol])[0] if allowed.any() else np.full(viol.sum(), np.inf)
        _, nn = tree.query(xy[viol], k=min(K + 1, len(xy)))
        frac = in_group[nn[:, 1:]].mean(1)
        cat = np.where(d_allowed < BOUNDARY_UM, "boundary", np.where(frac >= CLUSTER_FRAC, "clustered", "isolated"))
        sub = pd.DataFrame({
            "row": np.where(viol)[0], "constraint": c["id"], "region": region[viol],
            "base_label": np.array(types)[base_pred[viol]],
            "base_conf": base_conf[viol], "label": np.array(types)[pred[viol]],
            "category": cat,
            "d_allowed_um": d_allowed, "nbr_frac": frac})
        sub["status"] = np.where(pred[viol] == base_pred[viol], "kept", "overridden")
        sub["n_allowed_cells_of_group"] = int(ok.sum())
        rows.append(sub)
    cols = ["row", "constraint", "region", "base_label", "base_conf", "label",
            "category", "d_allowed_um", "nbr_frac", "status"]
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.concat(rows, ignore_index=True)


GROUP_KEYS = ("constraint", "region", "category", "status", "base_label")


def group_id(values):
    """Group ids are 'key=value;key=value' filters over GROUP_KEYS, e.g. 'constraint=C1;region=submucosa'."""
    return ";".join(f"{k}={v}" for k, v in values.items())


def select(conf, gid):
    m = np.ones(len(conf), bool)
    for part in filter(None, gid.split(";")):
        k, v = part.split("=", 1)
        if k not in GROUP_KEYS:
            raise ValueError(f"unknown group key {k!r}; valid: {GROUP_KEYS}")
        m &= conf[k].astype(str).to_numpy() == v
    return conf[m]


def summarize(conf, by=("constraint", "region")):
    if conf.empty:
        return pd.DataFrame(columns=list(by) + ["n", "kept_frac", "top_base_label", "group_id"])
    s = conf.groupby(list(by), dropna=False).agg(
        n=("row", "size"), kept_frac=("status", lambda v: round((v == "kept").mean(), 3)),
        top_base_label=("base_label", lambda v: v.value_counts().index[0]),
        top_category=("category", lambda v: v.value_counts().index[0]))
    s = s.sort_values("n", ascending=False).reset_index()
    s["group_id"] = [group_id(dict(zip(by, r))) for r in s[list(by)].itertuples(index=False)]
    return s
