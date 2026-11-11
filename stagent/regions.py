"""Per-cell region assignment from regions the user drew on the map.

Region state (stored in every version):
  {"source": "none" | "pathologist", "ops": [op, ...]}
  source "none":         every cell starts as "unknown" (no region constraint applies to it)
  source "pathologist":  start from the dataset's pathologist annotation, mapped to the canonical vocabulary
  op "polygon": {"op": "polygon", "id": "R1", "name": "mucosa", "polygon": [[px, py], ...] (H&E level-0 pixels),
                 "cells": [cell ids inside, computed when drawn], "author": ...}
Ops are applied in order, so a later polygon wins where polygons overlap. Region names are free text.
"""
import numpy as np
from . import vocab

MAX_NAME = 40


def default_state(source):
    if source not in ("none", "pathologist"):
        raise ValueError("region source must be 'none' (draw your own) or 'pathologist'")
    return {"source": source, "ops": []}


def assign(cells, state):
    """Return an array of region names, one per cell ("unknown" where no region applies). Overlapping outlines: the
    smaller wins (a region inside another keeps its cells); cells put into a region without an outline come last."""
    n = len(cells)
    if state.get("source") == "pathologist":
        raw = cells.pathologist_region.astype(object)  # SCP2927 names map to the canonical ones, others stay as they are
        reg = raw.map(lambda v: vocab.PATHOLOGIST_TO_REGION.get(v, v)).to_numpy(dtype=object)
        reg[~np.array([isinstance(v, str) and v not in ("", "nan", "Unlabelled", "unknown") for v in raw])] = "unknown"
    else:  # "none"; sessions saved with the removed "inferred" source also start empty
        reg = np.full(n, "unknown", dtype=object)
    ids = cells.cell.to_numpy()
    ops = state.get("ops", [])
    # where outlines overlap, the smaller one wins: a region drawn inside another keeps its cells, whichever was
    # drawn first (polygons from the largest to the smallest; equal ones in drawing order, the later winning)
    polys = sorted((op for op in ops if op["op"] == "polygon"), key=lambda op: -_area(op.get("polygon")))
    for op in polys:
        reg[np.isin(ids, op["cells"])] = op["name"]
    for op in ops:  # cells put into a region without an outline: as they were put, after every outline
        if op["op"] == "cells_override":
            reg[np.isin(ids, op["cells"])] = op["name"]
    return reg


def _area(polygon):
    """Area of an outline (shoelace; 0 for none)."""
    if not polygon or len(polygon) < 3:
        return 0.0
    p = np.asarray(polygon, float)
    x, y = p[:, 0], p[:, 1]
    return float(abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2)


def names(state, region_array=None):
    """Region names that exist in this version (drawn regions first, in drawing order)."""
    out = []
    for op in state.get("ops", []):
        if op["op"] in ("polygon", "cells_override") and op["name"] not in out:
            out.append(op["name"])
    if region_array is not None:
        out += [r for r in dict.fromkeys(region_array.tolist()) if r not in out and r not in vocab.FREE_REGIONS]
    return out


def validate_name(name):
    name = " ".join(str(name).split())
    if not name:
        raise vocab.VocabError("give the region a name")
    if len(name) > MAX_NAME:
        raise vocab.VocabError(f"region names are limited to {MAX_NAME} characters")
    if name in vocab.FREE_REGIONS:
        raise vocab.VocabError(f"'{name}' is reserved for cells outside every region")
    return name


def read_geojson(path):
    """The text of a GeoJSON file: UTF-8 (the GeoJSON standard; not the Windows code page, e.g. GBK on Chinese
    Windows, which fails on any non-ASCII name), else the system's encoding."""
    from pathlib import Path
    b = Path(path).expanduser().read_bytes()
    try:
        return b.decode("utf-8-sig")
    except UnicodeDecodeError:
        import locale
        return b.decode(locale.getpreferredencoding(False), errors="replace")


def parse_geojson(text):
    """Polygons of a GeoJSON file (e.g. QuPath annotations exported in the image's full-resolution pixels):
    [{"name", "color" ("#rrggbb" or None), "outer": [[x, y], ...], "holes": [ring, ...]}] and the number of objects
    left out (points, lines, detections and cells). A region's name is the object's classification, else its name.
    An outline traced along pixels (QuPath's brush and wand) becomes a smooth curve ("smoothed": True), so that it
    is not a staircase when zoomed in; others are kept point for point.
    QuPath writes a FeatureCollection or a plain list of features."""
    import json
    try:
        d = json.loads(text)
    except ValueError as e:
        raise ValueError(f"not a GeoJSON file: {e}")
    feats = d if isinstance(d, list) else d.get("features") if isinstance(d, dict) and d.get("type") == "FeatureCollection" \
        else [d] if isinstance(d, dict) else []
    out, skipped = [], 0
    for i, f in enumerate(feats):
        geom, props = (f or {}).get("geometry") or {}, (f or {}).get("properties") or {}
        if str(props.get("objectType", "annotation")).lower() in ("detection", "cell", "tile"):
            skipped += 1
            continue
        cls = props.get("classification") or {}
        name = (cls.get("name") if isinstance(cls, dict) else cls) or props.get("name") or "annotation"
        col = cls.get("color") if isinstance(cls, dict) else None
        color = "#%02x%02x%02x" % tuple(int(c) for c in col[:3]) if isinstance(col, list) and len(col) >= 3 else None
        polys = [geom["coordinates"]] if geom.get("type") == "Polygon" else geom.get("coordinates", []) \
            if geom.get("type") == "MultiPolygon" else None
        if not polys:
            skipped += 1
            continue
        for rings in polys:
            if rings and len(rings[0]) >= 3:
                rs = [[[float(x), float(y)] for x, y, *_ in r] for r in rings if len(r) >= 3]
                pix = is_pixel_outline(rs[0])  # QuPath's brush or wand: a staircase of pixels, made a smooth curve
                if pix:
                    rs = [smooth_pixel_outline(r) if is_pixel_outline(r) else r for r in rs]
                out.append({"name": " ".join(str(name).split())[:MAX_NAME], "color": color, "outer": rs[0],
                            "holes": rs[1:], "smoothed": pix})
    if not out:
        raise ValueError("the file has no polygons (QuPath: export the annotations as GeoJSON)")
    return out, skipped


def to_geojson(state, colors=None, image_size=None):
    """The drawn outlines of a version as a GeoJSON FeatureCollection in the image's full-resolution pixels, as QuPath
    writes its annotations, so QuPath opens it (Drag & Drop, or File > Import objects) and parse_geojson reads it
    back: every feature has a UUID id, objectType "annotation", the region as its name and as its classification
    {"name", "color": [r, g, b]} (QuPath's PathClass needs the color). colors: {region: "#rrggbb"} (else the app's
    palette). Cells moved into regions without an outline (cells_override) have no shape and are not in it.
    image_size: kept for the caller; not written (QuPath reads only the features). Returns (text, {"n_outlines",
    "regions", "n_without_outline"})."""
    import json, uuid
    from .render import region_colors
    names = list(dict.fromkeys(op["name"] for op in state.get("ops", []) if op.get("op") == "polygon"))
    cols = region_colors(names, colors or {})
    feats, moved = [], 0
    for op in state.get("ops", []):
        if op.get("op") != "polygon" or len(op.get("polygon") or []) < 3:
            moved += op.get("op") == "cells_override"
            continue
        rings = [op["polygon"]] + [h for h in op.get("holes", []) if len(h) >= 3]
        rings = [[[round(float(x), 2), round(float(y), 2)] for x, y in r] for r in rings]
        rings = [r + [r[0]] if r[0] != r[-1] else r for r in rings]  # GeoJSON rings are closed
        c = cols[op["name"]]
        fid = str(uuid.uuid5(uuid.NAMESPACE_URL, f"stannotator:{op.get('id')}:{op['name']}:{rings[0][:3]}"))
        feats.append({"type": "Feature", "id": fid, "geometry": {"type": "Polygon", "coordinates": rings},
                      "properties": {"objectType": "annotation", "name": op["name"],
                                     "classification": {"name": op["name"], "color": [int(c[k:k + 2], 16) for k in (1, 3, 5)]}}})
    if not feats:
        raise ValueError("no region outlines drawn in this version")
    return json.dumps({"type": "FeatureCollection", "features": feats}, ensure_ascii=False, separators=(",", ":")), \
        {"n_outlines": len(feats), "regions": names, "n_without_outline": moved}


def cells_in_polygon(px, py, polygon):
    poly = np.asarray(polygon, float)
    if len(poly) < 3:
        raise vocab.VocabError("a region needs at least 3 points")
    x0, y0 = poly.min(0)
    x1, y1 = poly.max(0)
    cand = np.where((px >= x0) & (px <= x1) & (py >= y0) & (py <= y1))[0]
    return cand[points_in(np.asarray(px, float)[cand], np.asarray(py, float)[cand], poly)]


def cells_in_pieces(px, py, pieces):
    """How many cells pieces [{"outer", "holes"}] hold: inside an outer ring and outside its holes, counted once (as
    session.set_region_shape assigns them)."""
    inside = np.zeros(len(px), bool)
    for pc in pieces:
        rows = cells_in_polygon(px, py, pc["outer"])
        for h in pc.get("holes") or []:
            if len(h) >= 3:
                rows = np.setdiff1d(rows, cells_in_polygon(px, py, h))
        inside[rows] = True
    return int(inside.sum())


def points_in(px, py, poly, chunk=4_000_000):
    """Which points lie inside a closed outline (even-odd rule): for each point, the edges that cross the horizontal
    line through it on its right. The edges are matched to the points through the points sorted by y, so the work
    grows with the crossings, not with points x edges (an outline of 50,000 points from QuPath: milliseconds)."""
    x1, y1 = poly[:, 0], poly[:, 1]
    x2, y2 = np.roll(x1, -1), np.roll(y1, -1)
    m = y1 != y2
    x1, y1, x2, y2 = x1[m], y1[m], x2[m], y2[m]
    lo, hi = np.minimum(y1, y2), np.maximum(y1, y2)
    order = np.argsort(py, kind="stable")
    ys = py[order]
    s = np.searchsorted(ys, lo, "left")
    n = np.searchsorted(ys, hi, "left") - s  # the points with lo <= y < hi
    odd = np.zeros(len(px), bool)
    cs = np.cumsum(n)
    start = 0
    while start < len(n):  # edges in groups of about `chunk` (edge, point) pairs
        stop = max(start + 1, int(np.searchsorted(cs, (cs[start - 1] if start else 0) + chunk, "right")))
        nn = n[start:stop]
        tot = int(nn.sum())
        if tot:
            e = np.repeat(np.arange(start, stop), nn)
            c = order[s[e] + np.arange(tot) - np.repeat(np.cumsum(nn) - nn, nn)]
            xi = x1[e] + (py[c] - y1[e]) / (y2[e] - y1[e]) * (x2[e] - x1[e])
            odd ^= (np.bincount(c[px[c] < xi], minlength=len(px)) & 1).astype(bool)
        start = stop
    return odd


# ---------------- outlines on the pixel grid (QuPath's brush and wand) ----------------
def is_pixel_outline(ring):
    """An outline traced along image pixels: whole-pixel points, mostly horizontal and vertical steps. It looks smooth
    zoomed out and like a staircase zoomed in."""
    p = np.asarray(ring, float)
    if len(p) < 20:
        return False
    v = np.diff(p, axis=0)
    whole = np.mean(np.all(np.abs(p - np.round(p)) < 1e-6, axis=1))
    axis = np.mean((np.abs(v[:, 0]) < 1e-9) | (np.abs(v[:, 1]) < 1e-9))
    return whole > 0.9 and axis > 0.5


def _dp_closed(p, tol):
    """Douglas-Peucker on a closed ring: the points that keep it within tol."""
    n = len(p)
    if n <= 8:
        return p
    far = int(np.argmax(np.hypot(*(p - p[0]).T)))
    keep = np.zeros(n, bool)
    keep[0] = keep[far] = True
    stack = [(0, far), (far, n)]
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        A, B = p[a], p[b % n]
        d = B - A
        L = float(np.hypot(*d)) or 1e-9
        seg = p[a + 1:b]
        dist = np.abs(d[1] * (seg[:, 0] - A[0]) - d[0] * (seg[:, 1] - A[1])) / L
        i = int(np.argmax(dist))
        if dist[i] > tol:
            m = a + 1 + i
            keep[m] = True
            stack += [(a, m), (m, b)]
    return p[keep]


def smooth_pixel_outline(ring, tol=1.5, iters=3):
    """The staircase of a pixel outline as a smooth curve: its steps dropped (within tol pixels), the corners cut
    (Chaikin), the straight runs thinned. On QuPath outlines of a colon section 99% of the curve stays within 1.7
    pixels of the original (0.4 um on a Xenium H&E) and 0.07% of the cells change side; only narrow tips are rounded
    more."""
    p = np.asarray(ring, float)
    if len(p) > 1 and np.allclose(p[0], p[-1]):
        p = p[:-1]
    q = _dp_closed(p, tol)
    if len(q) < 4:
        return [[float(x), float(y)] for x, y in p]
    for _ in range(iters):
        r = np.roll(q, -1, axis=0)
        q = np.stack([0.75 * q + 0.25 * r, 0.25 * q + 0.75 * r], axis=1).reshape(-1, 2)
    q = _dp_closed(q, 0.2)
    return [[round(float(x), 3), round(float(y), 3)] for x, y in q]


# ---------------- outlines made by the agent: valid as the viewer's tools save them ----------------
def _clean_ring(r):
    """A ring without its closing point and repeated consecutive points."""
    p = np.asarray(r, float).reshape(-1, 2)
    if len(p) > 1 and np.allclose(p[0], p[-1]):
        p = p[:-1]
    if len(p) > 1:
        p = p[np.r_[True, np.hypot(*np.diff(p, axis=0).T) > 1e-6]]
    if len(p) > 1 and np.hypot(*(p[0] - p[-1])) <= 1e-6:
        p = p[:-1]
    return p


def ring_area(p):
    x, y = np.asarray(p, float).T
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def _edges_meet(rings):
    """Whether two edges of the rings cross or touch (apart from neighbors in one ring), or a point repeats (rounded
    to 1/16 pixel, as the viewer checks): such a shape is not valid for QuPath. Edges are bucketed in a grid, so the
    work grows with the edges, not with their square."""
    k = np.round(np.concatenate(rings) * 16).astype(np.int64)
    if len(np.unique(k, axis=0)) < len(k):
        return True
    A = np.concatenate(rings)
    B = np.concatenate([np.roll(r, -1, 0) for r in rings])
    ring = np.concatenate([np.full(len(r), i) for i, r in enumerate(rings)])
    pos = np.concatenate([np.arange(len(r)) for r in rings])
    size = np.array([len(r) for r in rings])[ring]
    lo, hi = np.minimum(A, B), np.maximum(A, B)
    cell = max(4 * float(np.median(np.hypot(*(B - A).T))), 1e-3)
    c0, c1 = np.floor(lo / cell).astype(np.int64), np.floor(hi / cell).astype(np.int64)
    nx, ny = c1[:, 0] - c0[:, 0] + 1, c1[:, 1] - c0[:, 1] + 1
    cnt = nx * ny
    e = np.repeat(np.arange(len(A)), cnt)
    off = np.arange(int(cnt.sum())) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    cx, cy = c0[e, 0] + off % nx[e], c0[e, 1] + off // nx[e]
    key = (cx - cx.min()) * (int(cy.max() - cy.min()) + 1) + (cy - cy.min())
    o = np.argsort(key, kind="stable")
    key, e = key[o], e[o]
    starts = np.r_[0, np.nonzero(np.diff(key))[0] + 1]
    lens = np.diff(np.r_[starts, len(key)])
    pi, pj = [], []
    for m in np.unique(lens[lens > 1]):  # all pairs within each cell, cells of the same size together
        idx = e[starts[lens == m][:, None] + np.arange(m)]
        a, b = np.triu_indices(m, 1)
        pi.append(idx[:, a].ravel()); pj.append(idx[:, b].ravel())
    if not pi:
        return False
    i, j = np.concatenate(pi), np.concatenate(pj)
    ij = np.unique(np.c_[np.minimum(i, j), np.maximum(i, j)], axis=0)
    i, j = ij[:, 0], ij[:, 1]
    same = ring[i] == ring[j]
    d = np.abs(pos[i] - pos[j])
    nb = same & ((d == 1) | (d == size[i] - 1))  # neighbors share a point by construction
    i, j = i[~nb], j[~nb]
    box = (lo[i] <= hi[j]).all(1) & (lo[j] <= hi[i]).all(1)
    i, j = i[box], j[box]
    if not len(i):
        return False
    p, r, q, s = A[i], B[i] - A[i], A[j], B[j] - A[j]
    cr = lambda u, v: u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0]
    d1, d2 = cr(r, q - p), cr(r, q + s - p)
    d3, d4 = cr(s, p - q), cr(s, p + r - q)
    proper = (np.sign(d1) * np.sign(d2) <= 0) & (np.sign(d3) * np.sign(d4) <= 0)
    collinear = (d1 == 0) & (d2 == 0)
    if collinear.any():  # on one line: they meet only where their spans overlap
        t = np.abs(r).argmax(1)
        a0, a1 = np.minimum(p[np.arange(len(t)), t], (p + r)[np.arange(len(t)), t]), np.maximum(p[np.arange(len(t)), t], (p + r)[np.arange(len(t)), t])
        b0, b1 = np.minimum(q[np.arange(len(t)), t], (q + s)[np.arange(len(t)), t]), np.maximum(q[np.arange(len(t)), t], (q + s)[np.arange(len(t)), t])
        proper = np.where(collinear, (a0 <= b1) & (b0 <= a1), proper)
    return bool(proper.any())


def _inside(pt, ring):
    return bool(points_in(np.array([pt[0]]), np.array([pt[1]]), np.asarray(ring, float))[0])


def _valid(pieces):
    """Pieces that the viewer's tools would save as they are: simple rings, holes inside their outline and apart,
    no piece inside another (other than in its hole)."""
    rings = [r for pc in pieces for r in [pc["outer"], *pc["holes"]]]
    if _edges_meet(rings):
        return False
    for pc in pieces:
        if not all(_inside(h[0], pc["outer"]) for h in pc["holes"]):
            return False
        if any(_inside(a[0], b) for a in pc["holes"] for b in pc["holes"] if a is not b):
            return False
    for a in pieces:
        for b in pieces:
            if a is not b and _inside(a["outer"][0], b["outer"]) and not any(_inside(a["outer"][0], h) for h in b["holes"]):
                return False
    return True


def _retrace(pieces, max_side=4000):
    """The area of the pieces (each its outline minus its holes, all together) drawn into a mask and traced again:
    pieces that do not overlap, rings that neither cross nor touch."""
    import contourpy
    from scipy import ndimage as ndi
    allp = np.concatenate([pc["outer"] for pc in pieces])
    x0, y0 = allp.min(0) - 2
    x1, y1 = allp.max(0) + 2
    step = max(1.0, max(x1 - x0, y1 - y0) / max_side)  # image pixels per mask pixel
    w, h = int(np.ceil((x1 - x0) / step)) + 1, int(np.ceil((y1 - y0) / step)) + 1
    gx, gy = np.meshgrid(x0 + (np.arange(w) + 0.5) * step, y0 + (np.arange(h) + 0.5) * step)
    gx, gy = gx.ravel(), gy.ravel()
    m = np.zeros(w * h, bool)
    for pc in pieces:
        inn = points_in(gx, gy, pc["outer"])
        for hl in pc["holes"]:
            inn &= ~points_in(gx, gy, hl)
        m |= inn
    return _trace_mask(m.reshape(h, w), step, x0, y0)


def _trace_mask(m, step, x0=0.0, y0=0.0):
    """The pieces of a mask (pixel i covering x0 + [i, i + 1] * step image pixels) traced between its pixels: outlines
    with their holes, apart, rings neither crossing nor touching; crumbs of 4 square pixels or less dropped."""
    import contourpy
    from scipy import ndimage as ndi

    def trace(mask):  # the longest boundary of a mask, between its pixels
        z = ndi.gaussian_filter(np.pad(mask, 3).astype(np.float32), 0.7)
        lines = contourpy.contour_generator(z=z).lines(0.5)
        if not lines:  # (a hole of a pixel or so: gone with the blur)
            return np.zeros((0, 2))
        return max((np.asarray(ln, float) for ln in lines), key=len)[:-1] - 3
    out = []
    comp, n = ndi.label(m)
    for k, sl in enumerate(ndi.find_objects(comp), 1):
        c = comp[sl] == k
        if c.sum() * step * step <= 4:
            continue
        filled = ndi.binary_fill_holes(c)
        oy, ox = sl[0].start, sl[1].start
        to_img = lambda p: np.c_[x0 + (p[:, 0] + ox + 0.5) * step, y0 + (p[:, 1] + oy + 0.5) * step]
        o = trace(filled)
        if len(o) < 3:
            continue
        hl, nh = ndi.label(filled & ~c)
        holes = [to_img(h) for h in (trace(hl[s2] == j) + [s2[1].start, s2[0].start] if s2 is not None else np.zeros((0, 2))
                                     for j, s2 in enumerate(ndi.find_objects(hl), 1)) if len(h) >= 3]
        out.append({"outer": to_img(o), "holes": [x for x in holes if abs(ring_area(x)) > 4]})
    return out


def normalize_pieces(pieces):
    """Outlines made by the agent (from the H&E, a compartment or tiles), made as the viewer's tools save a shape:
    rings without repeated points, crumbs (4 square pixels or less) dropped, holes inside their outline; when rings
    cross or touch, or pieces overlap, the area is traced again (as one shape: pieces apart, rings strictly simple,
    as QuPath needs them). pieces: [{"outer", "holes"}] or bare rings, in the image's full-resolution pixels."""
    clean = []
    for pc in pieces:
        pc = pc if isinstance(pc, dict) else {"outer": pc, "holes": []}
        outer = _clean_ring(pc["outer"])
        crumb = lambda q: np.prod(q.max(0) - q.min(0)) <= 4  # (by its box: a ring crossing itself has no area of its own)
        if len(outer) < 3 or crumb(outer):
            continue
        holes = [h for h in (_clean_ring(x) for x in pc.get("holes") or []) if len(h) >= 3 and not crumb(h)]
        clean.append({"outer": outer, "holes": holes})
    if not clean:
        raise ValueError("no outline left (pieces of 4 square pixels or less)")
    if not _valid(clean):
        clean = _retrace(clean)
    r = lambda p: [[round(float(x), 2), round(float(y), 2)] for x, y in p]
    return [{"outer": r(pc["outer"]), "holes": [r(h) for h in pc["holes"]]} for pc in clean]
