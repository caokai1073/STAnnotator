"""Cells between the user's regions that belong to one of them, found for the agent to propose (never applied here).

candidates(): by position only (cell types are not used, so a region prior can never confirm itself):
  gaps      cells in no region within max_gap_um of two regions or more (between_only; else of any region, which
            also grows a region into tissue nobody outlined): the region whose cells are nearest, or, when most of
            the cell's in-region neighbors sit in one region, that one
  boundary  cells of a region whose neighbors (K nearest) are mostly (min_share) in another region: an outline drawn
            a little off
The cells are grouped into patches (G1, G2, ...: connected within LINK_UM; the cells of a gap between two regions
go to either side) and kept in the
annotation's folder (region_fill.json) with the version they were computed on; the agent looks at them on the H&E
(figure) and proposes the patches that hold up (mcp_server.assign_region_cells). The user's Accept applies them as
"cells_override" region operations (Session.assign_region_cells): the outlines stay as drawn, Undo takes it back.
"""
import json, time
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from . import render
from .osutil import dump_json

LINK_UM = 25.0   # candidate cells this close join one patch (each cell keeps its own target region)
MAX_PATCHES = 40


def _file(s):
    return s.dir / "region_fill.json"


def candidates(s, max_gap_um=30.0, min_share=0.7, regions=None, include_boundary=True, between_only=True, version=None):
    """Patches of cells to move into a region (see the module doc), stored for figure() and the proposal."""
    r = s.result(version)
    reg = np.asarray(r["region"], dtype=object)
    names = [n for n in s.region_names(version) if n != "unknown"]
    if not names:
        raise ValueError("there are no regions yet: draw or import them first")
    targets = [n for n in (regions or names) if n in names]
    if regions and len(targets) != len(regions):
        raise ValueError(f"unknown region(s): {[x for x in regions if x not in names]}; regions: {names}")
    xy = s.cells[["x", "y"]].to_numpy()
    n = len(xy)
    # distance of every cell to the nearest cell of each target region
    dist = np.full((n, len(targets)), np.inf)
    for j, name in enumerate(targets):
        m = reg == name
        if m.any():
            dist[:, j] = cKDTree(xy[m]).query(xy, distance_upper_bound=max(max_gap_um * 3, 60))[0]
    nb = s.knn()
    nreg = reg[nb]                                                  # neighbors' regions
    in_reg = nreg != "unknown"
    to = np.full(n, None, dtype=object)
    share = np.zeros(n)
    for j, name in enumerate(targets):  # the region most of the in-region neighbors sit in
        c = (nreg == name).sum(1)
        better = c > share * nb.shape[1]
        to[better & (c > 0)] = name
        share = np.maximum(share, c / nb.shape[1])
    near = np.array(targets, dtype=object)[np.argmin(dist, 1)]
    dmin = dist.min(1)
    n_in = in_reg.sum(1)
    major = np.where(n_in > 0, share * nb.shape[1] / np.maximum(n_in, 1), 0)   # share among in-region neighbors
    gap = (reg == "unknown") & (dmin <= max_gap_um)
    if between_only:  # between two regions, not at the edge of one region facing tissue nobody outlined
        gap &= (dist <= max_gap_um).sum(1) >= 2
    gap_to = np.where((major >= min_share) & (to != None), to, near)       # noqa: E711 (object array)
    move_to = np.full(n, None, dtype=object)
    move_to[gap] = gap_to[gap]
    why = np.full(n, "", dtype=object)
    why[gap] = "gap"
    if include_boundary:
        bnd = (reg != "unknown") & (to != None) & (to != reg) & (share >= min_share)   # noqa: E711
        bnd &= np.isin(to, targets)
        move_to[bnd] = to[bnd]
        why[bnd] = "boundary"
    rows = np.where(move_to != None)[0]                                    # noqa: E711
    patches = []
    if len(rows):  # a patch: candidate cells connected within LINK_UM, whatever region each goes to
        pairs = cKDTree(xy[rows]).query_pairs(LINK_UM, output_type="ndarray")
        g = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])) if len(pairs) else ([], ([], [])),
                       shape=(len(rows), len(rows)))
        _, comp = connected_components(g, directed=False)
        sizes = np.bincount(comp)
        types = np.array(s.types, dtype=object)[r["pred"]]
        px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
        ids = s.cells.cell.to_numpy()
        count = lambda a: {str(k): int(v) for k, v in pd.Series(a).value_counts().items()}
        for k in np.argsort(-sizes):
            m = rows[comp == k]
            near_m = [(targets[i], float(np.median(dist[m, i]))) for i in range(len(targets))
                      if np.median(dist[m, i]) <= max_gap_um * 2]
            patches.append({"n_cells": int(len(m)), "moves": count(move_to[m]),
                            "now": {("no region" if a == "unknown" else a): b for a, b in count(reg[m]).items()},
                            "kinds": count(why[m]),
                            "regions_near": [{"region": a, "distance_um": round(d, 1)} for a, d in sorted(near_m, key=lambda x: x[1])],
                            "neighbors_in_target": round(float(np.median(share[m])), 2),
                            "cell_types": {str(a): int(b) for a, b in pd.Series(types[m]).value_counts().head(4).items()},
                            "center_px": [round(float(px[m].mean())), round(float(py[m].mean()))],
                            "box_px": [round(float(px[m].min())), round(float(py[m].min())),
                                       round(float(px[m].max())), round(float(py[m].max()))],
                            "cell_ids": ids[m].tolist(), "cell_to": move_to[m].tolist()})
    small = [p for p in patches if p["n_cells"] < 3]
    patches = [p for p in patches if p["n_cells"] >= 3]
    for i, p in enumerate(patches):
        p["id"] = f"G{i + 1}"
    st = {"version": r["state"]["version"], "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
          "params": {"max_gap_um": max_gap_um, "min_share": min_share, "regions": targets, "include_boundary": include_boundary,
                     "between_only": between_only},
          "patches": patches}
    dump_json(st, _file(s))
    out = {k: v for k, v in st.items() if k != "patches"}
    out["n_patches"] = len(patches)
    out["n_cells"] = int(sum(p["n_cells"] for p in patches))
    out["left_out"] = f"{len(small)} single cells or pairs ({sum(p['n_cells'] for p in small)} cells)"
    into = {}
    for p in patches:
        for k2, v in p["moves"].items():
            into[k2] = into.get(k2, 0) + v
    out["cells_into"] = into
    out["patches"] = [{k: v for k, v in p.items() if k not in ("cell_ids", "cell_to")} for p in patches[:MAX_PATCHES]]
    out["note"] = ("found by position only (distance to the regions' cells and the regions of the K nearest neighbors); "
                   "cell types are listed for your judgment, but a move made because of them would let a region prior "
                   "confirm itself. Check the patches on the H&E (show_region_fill) before proposing them.")
    return out


def load(s, version=None):
    if not _file(s).exists():
        raise ValueError("no candidates yet: run preview_region_fill first")
    st = json.load(open(_file(s)))
    v = s.state(version)["version"]
    if st["version"] != v:
        raise ValueError(f"the candidates were computed on version {st['version']}, the annotation is at version {v}: "
                         "run preview_region_fill again")
    return st


def pick(s, patch_ids, to_region=None):
    """The moves of the chosen patches: [{"to", "cells", "patches", "from"}], one per target region (to_region: all
    of their cells into that one)."""
    st = load(s)
    by = {p["id"]: p for p in st["patches"]}
    want = [str(x).strip().upper() for x in patch_ids]
    bad = [x for x in want if x not in by]
    if bad or not want:
        raise ValueError(f"unknown patch(es) {bad}; patches: G1..G{len(by)}" if bad else "choose at least one patch")
    if to_region is not None and to_region not in s.region_names():
        raise ValueError(f"no region {to_region!r}; regions: {s.region_names()}")
    reg = dict(zip(s.cells.cell.to_numpy(), s.result()["region"]))
    moves = {}
    for x in want:
        for cid, t in zip(by[x]["cell_ids"], by[x]["cell_to"]):
            t = to_region or t
            m = moves.setdefault(t, {"to": t, "cells": [], "patches": [], "from": {}})
            m["cells"].append(cid)
            if x not in m["patches"]:
                m["patches"].append(x)
            f = reg.get(cid, "unknown")
            m["from"][f] = m["from"].get(f, 0) + 1
    return list(moves.values())


def figure(s, patch_ids=None, version=None):
    """Left: the sample's cells by region (gray: in none) with the patches outlined in their target's color.
    Right: up to four patches on the H&E (120 um or more): each cell of the patch ringed in its target region's
    color, the other cells as small dots in their region's color."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    st = load(s, version)
    r = s.result(version)
    reg = np.asarray(r["region"], dtype=object)
    names = [n for n in s.region_names(version) if n != "unknown"]
    cmap = plt.get_cmap("tab10")
    user = s.region_colors()
    col = {n: user.get(n) or matplotlib.colors.to_hex(cmap(i % 10)) for i, n in enumerate(names)}
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    by = {p["id"]: p for p in st["patches"]}
    show = [x.upper() for x in patch_ids] if patch_ids else [p["id"] for p in st["patches"][:4]]
    show = [x for x in show if x in by][:4]
    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_axes([0.01, 0.03, 0.47, 0.88])
    sub = np.random.default_rng(0).choice(len(px), min(len(px), 120000), replace=False)
    c = np.array([col.get(x, "#d0d0d0") for x in reg[sub]])
    ax.scatter(px[sub], py[sub], s=0.2, c=c, linewidths=0, rasterized=True)
    for p in st["patches"]:
        x0, y0, x1, y1 = p["box_px"]
        pad = 15 / s.he_um_per_px
        ax.add_patch(Rectangle((x0 - pad, y0 - pad), x1 - x0 + 2 * pad, y1 - y0 + 2 * pad, fill=False,
                               ec="k", lw=1.4 if p["id"] in show else 0.7))
        if p["id"] in show or p["n_cells"] >= 50:
            ax.text(x0 - pad, y0 - pad, p["id"], fontsize=7, color="k", va="bottom")
    inr = reg != "unknown"  # the regions and their surroundings, not the whole sample
    bx0, bx1, by0, by1 = px[inr].min(), px[inr].max(), py[inr].min(), py[inr].max()
    pad = 0.1 * max(bx1 - bx0, by1 - by0)
    ax.set_xlim(bx0 - pad, bx1 + pad); ax.set_ylim(by1 + pad, by0 - pad)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for n in names:
        ax.scatter([], [], s=20, c=col[n], label=n)
    ax.legend(loc="upper right", fontsize=7, markerscale=1, frameon=True)
    ax.set_title(f"{len(st['patches'])} patches of cells to put into a region (boxes; gray: no region)", fontsize=9)
    sl = render.slide(__import__("stagent.importer", fromlist=["image_of"]).image_of(s.sample_id))
    upp = s.he_um_per_px
    for k, pid in enumerate(show):
        p = by[pid]
        cx, cy = p["center_px"]
        half = max(60 / upp, 0.6 * max(p["box_px"][2] - p["box_px"][0], p["box_px"][3] - p["box_px"][1]))
        a2 = fig.add_axes([0.51 + 0.245 * (k % 2), 0.5 - 0.47 * (k // 2), 0.23, 0.42])
        a2.imshow(np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, 420, 420)),
                  extent=(cx - half, cx + half, cy + half, cy - half))
        m = (np.abs(px - cx) < half) & (np.abs(py - cy) < half)
        a2.scatter(px[m], py[m], s=4, c=[col.get(x, "#ffffff") for x in reg[m]], linewidths=0)
        to = dict(zip(p["cell_ids"], p["cell_to"]))
        here = s.cells.cell.to_numpy()[m]
        sel = np.isin(here, p["cell_ids"])
        a2.scatter(px[m][sel], py[m][sel], s=22, facecolors="none", edgecolors=[col.get(to[c], "k") for c in here[sel]],
                   linewidths=1.2)
        a2.set_xlim(cx - half, cx + half); a2.set_ylim(cy + half, cy - half); a2.set_xticks([]); a2.set_yticks([])
        a2.set_title(f"{pid}: " + ", ".join(f"{v} -> {k}" for k, v in p["moves"].items()) + f"; {2 * half * upp:.0f} um\n"
                     "ringed: the patch's cells in their target's color", fontsize=7)
    out = s.dir / "figures" / f"region_fill_v{st['version']}.jpg"
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=90, pil_kwargs={"quality": 85}); plt.close(fig)
    return str(out)
