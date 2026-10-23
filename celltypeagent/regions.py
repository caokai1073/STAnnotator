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
    """Return an array of region names, one per cell ("unknown" where no region applies)."""
    n = len(cells)
    if state.get("source") == "pathologist":
        raw = cells.pathologist_region.astype(object)  # SCP2927 names map to the canonical ones, others stay as they are
        reg = raw.map(lambda v: vocab.PATHOLOGIST_TO_REGION.get(v, v)).to_numpy(dtype=object)
        reg[~np.array([isinstance(v, str) and v not in ("", "nan", "Unlabelled", "unknown") for v in raw])] = "unknown"
    else:  # "none"; sessions saved with the removed "inferred" source also start empty
        reg = np.full(n, "unknown", dtype=object)
    ids = cells.cell.to_numpy()
    for op in state.get("ops", []):
        if op["op"] in ("polygon", "cells_override"):
            reg[np.isin(ids, op["cells"])] = op["name"]
    return reg


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
