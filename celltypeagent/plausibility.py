"""Spatial plausibility check at the level of tissue areas.

Evidence per cell: a cell type counts as cohesive (it forms sheets: epithelium, smooth muscle) when its lineage and
the type itself have a median same-lineage share of the 15 nearest cells of at least COHESIVE and at most RARE of
its cells are isolated (share below ISOLATED). Isolated cells of cohesive types are out of place; for each one the
best type outside its lineage (expression log-probability + log of that lineage's share around the cell) is noted.

Areas: the sample is cut into a grid of square tiles (named like a spreadsheet: rows A, B, ... top to bottom,
columns 1, 2, ... left to right, on the H&E as shown in the viewer). Tiles holding at least MIN_TILE out-of-place
cells of one lineage, and ENRICH times more than expected if those cells were spread like all cells, are flagged;
flagged tiles sharing an edge form a candidate area (for example "Epithelial cells
scattered in a smooth muscle area, tiles F3-H6"). Claude reviews the candidates with the whole-slide overview, the
H&E and biology, and posts area suggestions for the user (suggest_areas); the user draws or adopts the region, and
a region prior or a correction is applied to it. Nothing here changes the annotation.

On three SCP2927 sections (dataset lineage split as truth, never shown to the app) 48-68% of out-of-place cells
with confidence < 0.6 and 19-45% of confident ones were in another lineage, against 1.6-6.5% of all cells of the
checked types; the suggested lineage was right for 70-78% of the wrong ones.
"""
import json, math, string, time
import numpy as np, pandas as pd
from scipy import ndimage
from . import markers, render

COHESIVE = 0.5   # minimum median same-lineage neighbor share of a cohesive lineage and type
ISOLATED = 0.2   # a cell is isolated below this share
RARE = 0.1       # at most this share of a cohesive type's cells is isolated
LOW_CONF = 0.6   # out-of-place cells below this confidence are likely mislabels
MIN_TILE = 6     # out-of-place cells of one lineage that flag a tile, if also
ENRICH = 3.0     # ... this many times the count expected if they were spread like all cells
MAX_TILES = 32   # grid tiles along the longer side of the sample
LIKELY = "likely mislabel"
CHECK = "isolated, model confident"


# ---------------- per-cell evidence ----------------
def check(s, version=None):
    """Per-cell status ('' / LIKELY / CHECK) and, for out-of-place cells, the suggested type (cached per version)."""
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
    prior = np.log(np.where(type_lin >= 0, nbr[:, np.maximum(type_lin, 0)], 0.0) + 0.05)
    banned = (type_lin[None, :] == code[rows][:, None]) | (types == markers.UNASSIGNED)[None, :]
    alt = types[np.where(banned, -np.inf, logp + prior).argmax(1)] if len(rows) else np.zeros(0, dtype=object)
    suggested = np.full(len(code), "", dtype=object)
    suggested[rows] = alt
    out = {"version": v, "status": status, "suggested": suggested, "lineage": lut[r["pred"]],
           "cohesion": {str(n): float(by_lin[i]) for i, n in enumerate(names)},
           "cohesive_types": [str(types[t]) for t in cohesive],
           "n_flagged": int(flag.sum()), "n_likely": int((status == LIKELY).sum())}
    cache[v] = out
    return out


# ---------------- grid ----------------
def _row_name(i):
    s = ""
    i += 1
    while i:
        i, k = divmod(i - 1, 26)
        s = string.ascii_uppercase[k] + s
    return s


def grid(s):
    """Square tiles over the sample in H&E level-0 pixels; tile size a round number of um."""
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    x0, y0 = float(px.min()), float(py.min())
    long_um = max(px.max() - x0, py.max() - y0) * s.he_um_per_px
    tile_um = max(100, int(math.ceil(long_um / MAX_TILES / 50.0)) * 50)
    tile_px = tile_um / s.he_um_per_px
    nx, ny = int((px.max() - x0) // tile_px) + 1, int((py.max() - y0) // tile_px) + 1
    return {"x0": x0, "y0": y0, "tile_px": tile_px, "tile_um": tile_um, "nx": nx, "ny": ny,
            "rows": [_row_name(i) for i in range(ny)], "cols": [str(j + 1) for j in range(nx)]}


def tile_index(s, g=None):
    g = g or grid(s)
    j = ((s.cells.px.to_numpy() - g["x0"]) // g["tile_px"]).astype(int)
    i = ((s.cells.py.to_numpy() - g["y0"]) // g["tile_px"]).astype(int)
    return i, j


def tile_name(g, i, j):
    return f"{g['rows'][i]}{g['cols'][j]}"


def parse_tiles(g, spec):
    """'C3', 'C3:E6' (a block), or several separated by commas -> set of (row, col)."""
    def one(t):
        t = t.strip().upper()
        k = len(t) - len(t.lstrip(string.ascii_uppercase))
        if not k or not t[k:].isdigit() or t[:k] not in g["rows"] or not 1 <= int(t[k:]) <= g["nx"]:
            raise ValueError(f"unknown tile {t!r}: rows {g['rows'][0]}-{g['rows'][-1]}, columns 1-{g['nx']}")
        return g["rows"].index(t[:k]), int(t[k:]) - 1
    out = set()
    for part in str(spec).replace(";", ",").split(","):
        if not part.strip():
            continue
        if ":" in part:
            (a, b), (c, d) = (one(x) for x in part.split(":", 1))
            out |= {(i, j) for i in range(min(a, c), max(a, c) + 1) for j in range(min(b, d), max(b, d) + 1)}
        else:
            out.add(one(part))
    if not out:
        raise ValueError("no tiles given")
    return out


def tiles_text(g, tiles):
    """Compact tile list: one 'A1:B3' block per run of rows sharing the same column span."""
    by_row = {}
    for i, j in sorted(tiles):
        by_row.setdefault(i, []).append(j)
    parts = []
    for i, js in by_row.items():
        start = prev = js[0]
        for j in js[1:] + [None]:
            if j is not None and j == prev + 1:
                prev = j
                continue
            parts.append(tile_name(g, i, start) + (f":{tile_name(g, i, prev)}" if prev != start else ""))
            if j is not None:
                start = prev = j
    return ", ".join(parts)


def outline(g, tiles):
    """Outer boundary of a set of tiles as one polygon (H&E level-0 pixels); holes are ignored."""
    tiles = set(tiles)
    edges = {}
    for i, j in tiles:  # directed edges, interior on the left; shared edges cancel
        for a, b, ni, nj in (((j, i), (j + 1, i), i - 1, j), ((j + 1, i), (j + 1, i + 1), i, j + 1),
                             ((j + 1, i + 1), (j, i + 1), i + 1, j), ((j, i + 1), (j, i), i, j - 1)):
            if (ni, nj) not in tiles:
                edges.setdefault(a, []).append(b)
    loops = []
    while edges:
        start = next(iter(edges))
        loop, p = [start], start
        while True:
            q = edges[p].pop()
            if not edges[p]:
                del edges[p]
            if q == start:
                break
            loop.append(q)
            p = q
        loops.append(loop)
    loop = max(loops, key=len)
    pts = [loop[k] for k in range(len(loop))  # drop points in the middle of straight runs
           if not (loop[k - 1][0] == loop[k][0] == loop[(k + 1) % len(loop)][0]
                   or loop[k - 1][1] == loop[k][1] == loop[(k + 1) % len(loop)][1])]
    return [[g["x0"] + x * g["tile_px"], g["y0"] + y * g["tile_px"]] for x, y in pts]


# ---------------- candidate areas ----------------
def _area_stats(s, r, res, g, ti, tj, tiles, lineage=None):
    """Composition and out-of-place cells of a set of tiles."""
    key = ti * g["nx"] + tj
    inside = np.isin(key, [i * g["nx"] + j for i, j in tiles])
    lin = res["lineage"]
    comp = pd.Series(lin[inside]).value_counts()
    oop = inside & (res["status"] != "") & ((lin == lineage) if lineage else True)
    out = {"tiles": tiles_text(g, tiles), "n_cells": int(inside.sum()),
           "composition": {k: round(v / max(comp.sum(), 1), 3) for k, v in comp.head(6).items()},
           "n_out_of_place": int(oop.sum()), "n_likely_mislabel": int((res["status"][oop] == LIKELY).sum()),
           "out_of_place_types": pd.Series(np.array(s.types, dtype=object)[r["pred"][oop]]).value_counts().head(4).to_dict(),
           "suggested_types": pd.Series(res["suggested"][oop]).value_counts().head(4).to_dict(),
           "median_confidence": round(float(np.median(r["conf"][oop])), 3) if oop.any() else None}
    reg = pd.Series(r["region"][oop]).value_counts()
    out["regions"] = {k: int(v) for k, v in reg.items()}
    return out, inside, oop


def scan(s, version=None):
    """Candidate areas: touching tiles with >= MIN_TILE out-of-place cells of the same lineage (cached per version)."""
    res = check(s, version)
    cache = s.__dict__.setdefault("_areas", {})
    if res["version"] in cache:
        return cache[res["version"]]
    r = s.result(version)
    g = grid(s)
    ti, tj = tile_index(s, g)
    areas = []
    oop_all = res["status"] != ""
    density = np.zeros((g["ny"], g["nx"]))
    np.add.at(density, (ti, tj), 1)
    density /= density.sum()
    for lin in sorted(set(res["lineage"][oop_all])):
        m = oop_all & (res["lineage"] == lin)
        counts = np.zeros((g["ny"], g["nx"]), int)
        np.add.at(counts, (ti[m], tj[m]), 1)
        flag = (counts >= MIN_TILE) & (counts >= ENRICH * m.sum() * density)  # a local excess, not a diffuse spread
        lab, n = ndimage.label(flag)  # tiles sharing an edge; corners alone do not join
        for k in range(1, n + 1):
            tiles = {(int(i), int(j)) for i, j in zip(*np.where(lab == k))}
            st, inside, oop = _area_stats(s, r, res, g, ti, tj, tiles, lin)
            px, py = s.cells.px.to_numpy()[oop], s.cells.py.to_numpy()[oop]
            areas.append({"lineage": lin, **st,
                          "category": LIKELY if st["n_likely_mislabel"] * 2 >= st["n_out_of_place"] else CHECK,
                          "outline_px": outline(g, tiles), "tile_set": sorted(tiles),
                          "bbox_px": [float(px.min()), float(py.min()), float(px.max()), float(py.max())]})
    areas.sort(key=lambda a: -a["n_out_of_place"])
    for i, a in enumerate(areas):
        a["id"] = f"A{i + 1}"
    by_lin = {}
    for lin in sorted(set(res["lineage"][oop_all])):  # also report diffuse spreads, which form no area
        m = oop_all & (res["lineage"] == lin)
        in_areas = sum(a["n_out_of_place"] for a in areas if a["lineage"] == lin)
        by_lin[lin] = {"n_out_of_place": int(m.sum()), "in_candidate_areas": int(in_areas),
                       "by_region": {k: int(v) for k, v in pd.Series(r["region"][m]).value_counts().head(6).items()}}
    out = {"version": res["version"], "grid": g, "areas": areas, "by_lineage": by_lin}
    cache[res["version"]] = out
    return out


def lineage_by_region(s, version=None):
    """For each drawn region, how many cells of every lineage it holds (empty when no regions are drawn)."""
    r = s.result(version)
    if not s.region_names(version):
        return {}
    lut = np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)
    t = pd.crosstab(lut[r["pred"]], r["region"])
    return {reg: {k: int(n) for k, n in t[reg].items() if n} for reg in t.columns}


AREA_KEYS = ("id", "lineage", "category", "tiles", "n_out_of_place", "n_likely_mislabel", "out_of_place_types",
             "suggested_types", "median_confidence", "composition", "n_cells", "regions")


def summary(s, version=None, top=12):
    """What the agent needs to review the sample (no per-cell lists)."""
    res, sc = check(s, version), scan(s, version)
    g = sc["grid"]
    return {"version": res["version"],
            "grid": {"tile_um": g["tile_um"], "rows": f"{g['rows'][0]}-{g['rows'][-1]} (top to bottom)",
                     "columns": f"1-{g['nx']} (left to right)"},
            "cohesive_types": res["cohesive_types"],
            "n_out_of_place_cells": res["n_flagged"], "n_likely_mislabels": res["n_likely"],
            "n_candidate_areas": len(sc["areas"]),
            "out_of_place_by_lineage": sc["by_lineage"],
            "candidate_areas": [{k: a[k] for k in AREA_KEYS} for a in sc["areas"][:top]],
            "lineages_by_region": lineage_by_region(s, version),
            "note": f"out of place = cell of a type that forms sheets with < {ISOLATED:.0%} of its 15 nearest cells in "
                    f"its lineage; a candidate area = touching {g['tile_um']} um tiles with >= {MIN_TILE} such cells of "
                    f"one lineage and >= {ENRICH:g}x the count expected from cell density (a local excess; a diffuse "
                    "spread shows only in out_of_place_by_lineage: view it with show_tissue_overview(highlight)); composition = lineage shares of all cells in the area. Candidates are computed from "
                    "the annotation only: confirm them on the overview and H&E before suggesting anything."}


def area_info(s, spec, version=None):
    """A candidate area id ('A3') or a tile spec ('C3:E6, F4') -> tiles and statistics."""
    sc = scan(s, version)
    g = sc["grid"]
    a = next((x for x in sc["areas"] if x["id"] == str(spec).strip().upper()), None)
    if a is not None:
        return a
    tiles = parse_tiles(g, spec)
    r, res = s.result(version), check(s, version)
    ti, tj = tile_index(s, g)
    st, inside, oop = _area_stats(s, r, res, g, ti, tj, tiles)
    return {"id": None, "lineage": None, **st, "outline_px": outline(g, tiles), "tile_set": sorted(tiles)}


# ---------------- figures ----------------
def fig_overview(s, version=None, highlight=None):
    """Whole sample: H&E with the tile grid (left) and cells colored by lineage with candidate areas (right)."""
    res, sc = check(s, version), scan(s, version)
    lin = res["lineage"]
    keep = None
    if highlight:
        keep = set()
        for h in highlight:
            keep |= set(s.type_groups.get(h, [h]))
        keep = {s.lineage_of_type.get(t, t) for t in keep} | set(highlight)
    out = s.dir / "figures" / f"overview_v{res['version']}_{'-'.join(sorted(keep or [])) or 'all'}.jpg".replace(" ", "_")
    return render.tissue_overview(s.meta["he_path"], s.cells.px.to_numpy(), s.cells.py.to_numpy(), lin,
                                  s.lineage_colors(), sc["grid"], sc["areas"], out, highlight=keep,
                                  title=f"sample {s.sample_id} v{res['version']}")


def fig_area(s, spec, version=None, highlight=None):
    a = area_info(s, spec, version)
    res = check(s, version)
    g = scan(s, version)["grid"]
    tiles = a["tile_set"]
    i0, i1 = min(t[0] for t in tiles), max(t[0] for t in tiles) + 1
    j0, j1 = min(t[1] for t in tiles), max(t[1] for t in tiles) + 1
    pad = 0.25 * g["tile_px"]
    box = (g["x0"] + j0 * g["tile_px"] - pad, g["y0"] + i0 * g["tile_px"] - pad,
           g["x0"] + j1 * g["tile_px"] + pad, g["y0"] + i1 * g["tile_px"] + pad)
    highlight = highlight or ([a["lineage"]] if a.get("lineage") else None)  # a candidate: its lineage in color
    keep = None
    if highlight:
        keep = {s.lineage_of_type.get(t, t) for h in highlight for t in s.type_groups.get(h, [h])} | set(highlight)
    out = s.dir / "figures" / f"area_v{res['version']}_{a['id'] or a['tiles'].replace(':', '-').replace(', ', '_')}.jpg"
    path = render.area_view(s.meta["he_path"], s.cells.px.to_numpy(), s.cells.py.to_numpy(), res["lineage"],
                            s.lineage_colors(), res["status"] != "", box, a["outline_px"], out, g, highlight=keep,
                            um_per_px=s.he_um_per_px,
                            title=f"sample {s.sample_id} v{res['version']} " + (a["id"] or "") + f" tiles {a['tiles']}")
    return path, a


# ---------------- area suggestions (written by Claude, shown in the viewer) ----------------
def _sfile(s):
    return s.dir / "area_suggestions.json"


def suggestions(s):
    p = _sfile(s)
    return json.load(open(p)) if p.exists() else []


def suggest(s, areas, author="agent"):
    """Replace the open suggestions with a new list: [{"tiles" or "area", "title", "reason", "action",
    "region_name"}]. Dismissed and done ones are kept."""
    g = scan(s)["grid"]
    old = [x for x in suggestions(s) if x["status"] != "open"]
    used = {x["id"] for x in old}
    new = []
    for a in areas:
        info = area_info(s, a.get("area") or a.get("tiles") or "")
        k = 1
        while f"S{k}" in used:
            k += 1
        used.add(f"S{k}")
        new.append({"id": f"S{k}", "title": str(a.get("title", "")).strip() or info["tiles"],
                    "reason": str(a.get("reason", "")).strip(), "action": str(a.get("action", "")).strip(),
                    "region_name": str(a.get("region_name", "")).strip(), "tiles": info["tiles"],
                    "from_area": info.get("id"), "outline_px": info["outline_px"], "lineage": info.get("lineage"),
                    "status": "open", "author": author, "created_at": time.strftime("%Y-%m-%dT%H:%M:%S")})
    json.dump(old + new, open(_sfile(s), "w"), indent=1)
    return new


def set_suggestion_status(s, sid, status, region_name=None):
    if status not in ("open", "done", "dismissed"):
        raise ValueError("status must be open, done or dismissed")
    rows = suggestions(s)
    hit = [x for x in rows if x["id"] == sid]
    if not hit:
        raise ValueError(f"no area suggestion {sid}")
    hit[0]["status"] = status
    if region_name:
        hit[0]["region_name"] = region_name
    json.dump(rows, open(_sfile(s), "w"), indent=1)
    return hit[0]
