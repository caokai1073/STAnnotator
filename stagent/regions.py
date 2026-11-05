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
                out.append({"name": " ".join(str(name).split())[:MAX_NAME], "color": color,
                            "outer": [[float(x), float(y)] for x, y, *_ in rings[0]],
                            "holes": [[[float(x), float(y)] for x, y, *_ in r] for r in rings[1:] if len(r) >= 3]})
    if not out:
        raise ValueError("the file has no polygons (QuPath: export the annotations as GeoJSON)")
    return out, skipped


def cells_in_polygon(px, py, polygon):
    from matplotlib.path import Path
    poly = np.asarray(polygon, float)
    if len(poly) < 3:
        raise vocab.VocabError("a region needs at least 3 points")
    x0, y0 = poly.min(0)
    x1, y1 = poly.max(0)
    cand = np.where((px >= x0) & (px <= x1) & (py >= y0) & (py <= y1))[0]
    inside = Path(poly).contains_points(np.c_[px[cand], py[cand]])
    return cand[inside]
