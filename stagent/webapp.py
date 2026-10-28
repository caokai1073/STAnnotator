"""Main web app: pick a sample, annotate it, view the result on the H&E, and talk to the agent.

Run:  python -m stagent.webapp [--port 8765]
Then open http://localhost:8765 (forward the port from the server first).

Everything runs on this server. The chat box calls the locally installed Claude Code CLI in print mode with
the logged-in account (no API key); its only tools are this app's MCP tools, served at /mcp by this same
process, so the agent and the viewer share one session object. The conflict review page is at /review.
"""
import argparse, asyncio, colorsys, contextlib, io, json, math, os, platform, re, shutil, sys, threading, time, uuid
from collections import OrderedDict
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import FileResponse, Response, StreamingResponse
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles
from . import config, vocab, mcp_server, review_app, render
from . import markers, markerqc, plausibility, claudelink, claudereview, importer, advice, niches, embedding
from . import communication, projects
from .session import Session, list_samples, list_sessions, panel_genes

STATIC = Path(__file__).parent / "static"
TILE, MARGIN_PX, CELL_RADIUS_UM = 512, 400, 3.45  # dot radius: about a cell's radius
CHAT_MCP = config.CHAT_MCP
JOBS: dict = {}
_views: dict = {}


def _json(obj, status=200):
    return review_app._json(obj, status)


class LRU(OrderedDict):
    def __init__(self, n):
        super().__init__(); self.n = n

    def get_or(self, key, fn):
        if key in self:
            self.move_to_end(key)
            return self[key]
        v = fn()
        self[key] = v
        if len(self) > self.n:
            self.popitem(last=False)
        return v


_tiles, _overlays, _colors = LRU(3000), LRU(3000), LRU(40)


# ---------------- H&E pyramid and sample geometry ----------------



class SampleView:
    """Virtual image = the sample's bounding box on its H&E slide (the slide holds two sections)."""

    def __init__(self, sample_id):
        meta = json.load(open(config.SAMPLES / sample_id / "meta.json"))
        self.sample_id = sample_id
        self.um_per_px = float(meta.get("he_um_per_px", config.HE_UM_PER_PX))
        self.slide = render.slide(meta.get("he_path"))  # H&E, fluorescence, or a blank background
        self.image_kind = meta.get("image_kind", "H&E")
        cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
        px, py = cells.px.to_numpy(), cells.py.to_numpy()
        self.px, self.py = px, py
        lo = -math.inf if isinstance(self.slide, render.BlankSlide) else 0  # a blank background has no edge at 0
        self.ox, self.oy = int(max(lo, math.floor(px.min()) - MARGIN_PX)), int(max(lo, math.floor(py.min()) - MARGIN_PX))
        x1, y1 = min(self.slide.W0, int(px.max()) + MARGIN_PX), min(self.slide.H0, int(py.max()) + MARGIN_PX)
        self.W, self.H = x1 - self.ox, y1 - self.oy
        self.max_level = int(math.ceil(math.log2(max(self.W, self.H))))
        self.order = np.argsort(px)
        self.sx = px[self.order]
        self.tree = cKDTree(np.c_[px, py])

    def tile_box(self, L, tx, ty):
        s = 2.0 ** (L - self.max_level)
        lw, lh = int(math.ceil(self.W * s)), int(math.ceil(self.H * s))
        w, h = min(TILE, lw - tx * TILE), min(TILE, lh - ty * TILE)
        if w <= 0 or h <= 0:
            raise ValueError("tile out of range")
        x0, y0 = self.ox + tx * TILE / s, self.oy + ty * TILE / s
        return x0, y0, x0 + w / s, y0 + h / s, w, h, s

    def cells_in(self, x0, y0, x1, y1, pad):
        i0, i1 = np.searchsorted(self.sx, [x0 - pad, x1 + pad])
        idx = self.order[i0:i1]
        return idx[(self.py[idx] >= y0 - pad) & (self.py[idx] <= y1 + pad)]


def view(sample_id) -> SampleView:
    """The sample's viewer geometry, rebuilt (with its image and dot tiles) when the cell positions change on disk."""
    rev = (config.SAMPLES / sample_id / "cells.parquet").stat().st_mtime
    if sample_id not in _views or _views[sample_id].rev != rev:
        if sample_id in _views:
            for cache in (_tiles, _overlays):
                for k in [k for k in cache if k[0] == sample_id or (len(k) > 8 and k[8] == sample_id)]:
                    cache.pop(k)
        _views[sample_id] = SampleView(sample_id)
        _views[sample_id].rev = rev
    return _views[sample_id]


# ---------------- colors ----------------
LINEAGE_HUES = {"epi": (0.02, 0.15), "fib": (0.20, 0.36), "stroma": (0.42, 0.50), "lymph": (0.56, 0.72), "myeloid": (0.78, 0.95)}
LINEAGE_RGB = {"epi": (235, 130, 20), "fib": (60, 165, 60), "stroma": (40, 160, 170), "lymph": (50, 100, 220), "myeloid": (210, 60, 170)}
MODES = ("type", "lineage", "region", "changed", "conflicts", "plausibility", "confidence", "cluster")
PLAUS_RGB = {plausibility.LIKELY: (230, 30, 30), plausibility.CHECK: (255, 170, 0)}
GRAY = (160, 160, 160)


def lineage_hues(s: Session):
    """Hue range per lineage: fixed for the reference taxonomy, evenly spread for a user's marker tree."""
    if s.scheme["kind"] == "reference":
        return dict(LINEAGE_HUES)
    lins = [l for l in s.lineages if l != markers.UNASSIGNED]
    step = 1.0 / max(len(lins), 1)
    return {l: (i * step + 0.1 * step, i * step + 0.8 * step) for i, l in enumerate(lins)}


def lineage_rgb(s: Session):
    if s.scheme["kind"] == "reference":
        return dict(LINEAGE_RGB)
    return {l: tuple(int(255 * c) for c in colorsys.hls_to_rgb((h0 + h1) / 2, 0.48, 0.75)) for l, (h0, h1) in lineage_hues(s).items()}


def lineage_label(s: Session, lin):
    return vocab.LINEAGE_GROUPS.get(lin, lin) if s.scheme["kind"] == "reference" else lin


def type_palette(s: Session):
    pal = {t: GRAY for t in s.types}
    for lin, (h0, h1) in lineage_hues(s).items():
        members = [t for t in s.types if s.lineage_of_type.get(t) == lin]
        for i, t in enumerate(members):
            h = h0 + (h1 - h0) * i / max(len(members) - 1, 1)
            l = (0.42, 0.58, 0.70)[i % 3]
            pal[t] = tuple(int(255 * c) for c in colorsys.hls_to_rgb(h, l, 0.80))
    return pal


def region_rgb(s: Session, version):
    names = s.region_names(version) + ["unknown"]
    hexes = render.region_colors(names, s.region_colors())
    return {n: (GRAY if n == "unknown" else tuple(int(hexes[n][i:i + 2], 16) for i in (1, 3, 5))) for n in names}


def cluster_ids(s: Session, version):
    if s.scheme["kind"] != "markers":
        return None, []
    names, inv = np.unique(s.marker_result(s.state(version))["cluster_of"].astype(str), return_inverse=True)
    return inv, list(names)


def cluster_rgb(i):
    return tuple(int(255 * c) for c in colorsys.hls_to_rgb((i * 0.618034) % 1.0, (0.45, 0.62)[i % 2], 0.75))


def colors_for(s: Session, version, mode, only):
    def make():
        r = s.result(version)
        n = len(s.cells)
        out = np.zeros((n, 4), np.uint8)
        if mode == "type":
            pal = type_palette(s)
            lut = np.array([pal[t] for t in s.types], np.uint8)
            out[:, :3], out[:, 3] = lut[r["pred"]], 230
        elif mode == "lineage":
            rgbs = lineage_rgb(s)
            lut = np.array([rgbs.get(s.lineage_of_type.get(t), GRAY) for t in s.types], np.uint8)
            out[:, :3], out[:, 3] = lut[r["pred"]], 230
        elif mode == "cluster":
            inv, names = cluster_ids(s, version)
            if inv is not None:
                out[:, :3] = np.array([cluster_rgb(i) for i in range(len(names))], np.uint8)[inv]
                out[:, 3] = 225
        elif mode == "region":
            rgbs = region_rgb(s, version)
            out[:, :3] = np.array([rgbs.get(x, GRAY) for x in r["region"]], np.uint8)
            out[:, 3] = 200
            out[r["region"] == "unknown", 3] = 70
        elif mode.startswith("lr:"):
            b, table = lr_bins(s, mode), lr_table(mode)
            out[:, :3] = np.array([x[2] for x in table], np.uint8)[b]
            out[:, 3] = np.array([x[3] for x in table], np.uint8)[b]
        elif mode.startswith("gene:"):
            b = gene_bins(s, mode)
            out[:, :3] = np.array([x[2] for x in GENE_BINS], np.uint8)[b]
            out[:, 3] = np.array([x[3] for x in GENE_BINS], np.uint8)[b]
        elif mode == "leiden":
            lab = embedding.leiden(s.sample_id)
            if lab is not None:
                out[:, :3] = np.array([cluster_rgb(i) for i in range(int(lab.max()) + 1)], np.uint8)[lab]
                out[:, 3] = 225
        elif mode == "niche":
            meta, lab = niches.current(s)
            if lab is not None:
                pal = niche_rgb(s)
                out[:, :3] = np.array(pal, np.uint8)[lab]
                out[:, 3] = 225
        elif mode == "changed":  # compared with the first annotation of the session
            out[:] = (170, 170, 170, 60)
            ch = r["pred"] != s.result(s.versions()[0])["pred"]
            out[ch] = (255, 0, 190, 245)
            out[r["locked_rows"]] = (0, 200, 255, 255)
        elif mode == "conflicts":
            out[:] = (170, 170, 170, 60)
            cf = r["conflicts"]
            if len(cf):
                rows = cf.row.to_numpy(dtype=int)
                out[rows] = np.where((cf.status == "kept").to_numpy()[:, None], [[230, 30, 30, 250]], [[255, 150, 0, 250]])
        elif mode == "plausibility":
            out[:] = (170, 170, 170, 50)
            stt = plausibility.check(s, version)["status"]
            for k, c in PLAUS_RGB.items():
                out[stt == k] = c + (250,)
        elif mode == "confidence":
            c = np.clip((r["conf"] - 0.3) / 0.7, 0, 1)
            out[:, 0] = (255 * np.clip(2 - 2 * c, 0, 1)).astype(np.uint8)
            out[:, 1] = (220 * np.clip(2 * c, 0, 1)).astype(np.uint8)
            out[:, 3] = 220
        if only:  # show only the cells whose legend entry is selected
            out[~np.isin(cell_keys(s, version, mode), list(only)), 3] = 0
        return out
    return _colors.get_or((s.id, version, mode_key(s, mode), tuple(only)), make)


def mode_key(s, mode):
    """Cache key of a color mode; niche colors change with every niche run."""
    if mode == "leiden":
        return f"leiden:{(embedding.current(s.sample_id) or {}).get('run_id')}"
    return f"niche:{(niches.current(s)[0] or {}).get('run_id')}" if mode == "niche" else mode


GENE_BINS = (("0", "not detected", (200, 200, 200), 35), ("1", "1 transcript", (254, 217, 118), 235),
             ("2-3", "2–3 transcripts", (253, 141, 60), 240), ("4-7", "4–7 transcripts", (227, 26, 28), 245),
             ("8+", "8 or more", (128, 0, 38), 250))


def gene_bins(s, mode):
    """Transcript-count bin (index into GENE_BINS) of every cell for color mode "gene:<name>"."""
    g = mode.split(":", 1)[1]
    i = {x.upper(): k for k, x in enumerate(s.panel)}.get(g.upper())
    if i is None:
        raise ValueError(f"{g} is not in this sample's gene panel")
    c = np.asarray(s.counts[:, i].todense()).ravel()
    return np.digitize(c, [1, 2, 4, 8])


_lr_cache = LRU(6)


def lr_table(mode):
    """The bins of a color mode "lr:<ligand>~<receptor>~<both|receiver|sender>": the cells' roles (both), or the
    percentiles of one side."""
    return communication.ROLE_BINS if mode.endswith("~both") else communication.BINS


def lr_bins(s, mode):
    """Bin (index into lr_table(mode)) of every cell."""
    lig, rec, side = (mode.split(":", 1)[1].split("~") + ["both"])[:3]

    def make():
        a = communication.activity(s, lig, rec, side)
        return communication.role_bins(*a) if side == "both" else communication.bins(a)
    return _lr_cache.get_or((s.sample_id, s._cells_mtime, mode), make)


def niche_rgb(s):
    meta, _ = niches.current(s)
    return [tuple(int(h[i:i + 2], 16) for i in (1, 3, 5)) for h in
            (render.DISTINCT * 3)[:meta["n_niches"]]] if meta else []


def cell_keys(s: Session, version, mode):
    """Legend key of every cell for a color mode (the keys legend_for returns)."""
    r = s.result(version)
    n = len(s.cells)
    if mode == "type":
        return np.array(s.types, dtype=object)[r["pred"]]
    if mode == "lineage":
        return np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)[r["pred"]]
    if mode == "region":
        return r["region"].astype(object)
    if mode == "changed":
        k = np.full(n, "same", dtype=object)
        k[r["pred"] != s.result(s.versions()[0])["pred"]] = "changed"
        k[r["locked_rows"]] = "hand"
        return k
    if mode == "conflicts":
        k = np.full(n, "none", dtype=object)
        cf = r["conflicts"]
        if len(cf):
            k[cf.row.to_numpy(dtype=int)] = np.where((cf.status == "kept").to_numpy(), "kept", "overridden")
        return k
    if mode == "plausibility":
        k = plausibility.check(s, version)["status"].copy()
        k[k == ""] = "ok"
        return k
    if mode.startswith("gene:"):
        return np.array([x[0] for x in GENE_BINS], dtype=object)[gene_bins(s, mode)]
    if mode.startswith("lr:"):
        return np.array([x[0] for x in lr_table(mode)], dtype=object)[lr_bins(s, mode)]
    if mode == "leiden":
        lab = embedding.leiden(s.sample_id)
        return np.array([f"L{i + 1}" for i in lab], dtype=object) if lab is not None else np.full(n, "", dtype=object)
    if mode == "niche":
        meta, lab = niches.current(s)
        return np.array([f"N{i + 1}" for i in range(meta["n_niches"])], dtype=object)[lab] if meta else np.full(n, "", dtype=object)
    if mode == "confidence":
        return np.where(r["conf"] < 0.5, "low", np.where(r["conf"] < 0.8, "mid", "high")).astype(object)
    if mode == "cluster":
        inv, names = cluster_ids(s, version)
        return np.array(names, dtype=object)[inv] if inv is not None else np.full(n, "", dtype=object)
    raise ValueError(f"unknown mode {mode}")


def legend_for(s: Session, version, mode):
    r = s.result(version)
    lin_of = s.lineage_of_type
    counts = np.bincount(r["pred"], minlength=len(s.types))
    order = list(lineage_hues(s)) + [markers.UNASSIGNED]
    if mode == "type":
        pal = type_palette(s)
        items = [{"key": t, "label": t, "color": pal[t], "count": int(counts[i]), "group": lineage_label(s, lin_of.get(t, t))}
                 for i, t in enumerate(s.types)]
        return sorted(items, key=lambda x: (order.index(lin_of.get(x["key"], markers.UNASSIGNED)) if lin_of.get(x["key"]) in order else 99, -x["count"]))
    if mode == "lineage":
        return [{"key": k, "label": lineage_label(s, k), "color": c,
                 "count": int(sum(counts[i] for i, t in enumerate(s.types) if lin_of.get(t) == k))} for k, c in lineage_rgb(s).items()]
    if mode == "cluster":
        inv, names = cluster_ids(s, version)
        if inv is None:
            return []
        rec = {x["id"]: x for x in s.clusters(version)}
        n = np.bincount(inv, minlength=len(names))
        return sorted([{"key": c, "label": f"{c} → {rec.get(c, {}).get('assigned') or 'unassigned'}", "color": cluster_rgb(i),
                        "count": int(n[i])} for i, c in enumerate(names)], key=lambda x: -x["count"])
    if mode == "region":
        vc = dict(zip(*np.unique(r["region"].astype(str), return_counts=True)))
        return [{"key": k, "label": k.replace("_", " ") if k != "unknown" else "no region", "color": c, "count": int(vc.get(k, 0))}
                for k, c in region_rgb(s, version).items()]
    if mode == "leiden":
        lab = embedding.leiden(s.sample_id)
        if lab is None:
            return []
        n = np.bincount(lab)
        return sorted([{"key": f"L{i + 1}", "label": f"L{i + 1}", "color": cluster_rgb(i), "count": int(n[i])}
                       for i in range(len(n))], key=lambda x: -x["count"])
    if mode.startswith("lr:"):
        table = lr_table(mode)
        n = np.bincount(lr_bins(s, mode), minlength=len(table))
        return [{"key": k, "label": lab, "color": c, "count": int(n[i])} for i, (k, lab, c, _) in enumerate(table)][::-1]
    if mode.startswith("gene:"):
        n = np.bincount(gene_bins(s, mode), minlength=len(GENE_BINS))
        return [{"key": k, "label": lab, "color": c, "count": int(n[i])} for i, (k, lab, c, _) in enumerate(GENE_BINS)][::-1]
    if mode == "niche":
        meta, lab = niches.current(s)
        if meta is None:
            return []
        n = np.bincount(lab, minlength=meta["n_niches"])
        names, pal = niches.niche_names(s), niche_rgb(s)
        return [{"key": f"N{i + 1}", "label": names[i] if names[i] == f"N{i + 1}" else f"N{i + 1} · {names[i]}",
                 "color": pal[i], "count": int(n[i])} for i in range(meta["n_niches"])]
    if mode == "changed":
        ch = r["pred"] != s.result(s.versions()[0])["pred"]
        ch[r["locked_rows"]] = False
        ch = int(ch.sum())
        return [{"key": "hand", "label": "labeled by you", "color": (0, 200, 255), "count": int(len(r["locked_rows"]))},
                {"key": "changed", "label": "changed since first annotation", "color": (255, 0, 190), "count": ch},
                {"key": "same", "label": "unchanged", "color": (170, 170, 170), "count": len(s.cells) - ch - int(len(r["locked_rows"]))}]
    if mode == "conflicts":
        cf = r["conflicts"]
        kept = int((cf.status == "kept").sum()) if len(cf) else 0
        return [{"key": "kept", "label": "conflict, expression kept", "color": (230, 30, 30), "count": kept},
                {"key": "overridden", "label": "conflict, prior won", "color": (255, 150, 0), "count": len(cf) - kept},
                {"key": "none", "label": "no conflict", "color": (170, 170, 170), "count": len(s.cells) - len(cf)}]
    if mode == "plausibility":
        stt = plausibility.check(s, version)["status"]
        n = {k: int((stt == k).sum()) for k in PLAUS_RGB}
        names = {plausibility.LIKELY: "out of place, model unsure", plausibility.CHECK: "out of place, model confident"}
        return [{"key": k, "label": names[k], "color": c, "count": n[k]} for k, c in PLAUS_RGB.items()] + \
               [{"key": "ok", "label": "fits its surroundings or not checked", "color": (170, 170, 170),
                 "count": len(s.cells) - sum(n.values())}]
    return [{"key": "low", "label": "low confidence", "color": (255, 0, 0), "count": int((r["conf"] < 0.5).sum())},
            {"key": "mid", "label": "medium", "color": (255, 220, 0), "count": int(((r["conf"] >= 0.5) & (r["conf"] < 0.8)).sum())},
            {"key": "high", "label": "high", "color": (0, 220, 0), "count": int((r["conf"] >= 0.8).sum())}]


# ---------------- tiles ----------------
def he_tile(req: Request):
    sid, L = req.path_params["sample"], int(req.path_params["level"])
    tx, ty = (int(v) for v in req.path_params["xy"].split(".")[0].split("_"))
    v = view(sid)

    def make():
        x0, y0, x1, y1, w, h, _ = v.tile_box(L, tx, ty)
        buf = io.BytesIO()
        v.slide.read(x0, y0, x1, y1, w, h).save(buf, "JPEG", quality=85)
        return buf.getvalue()
    try:
        data = _tiles.get_or((sid, L, tx, ty), make)
    except ValueError:
        return Response(status_code=404)
    return Response(data, media_type="image/jpeg", headers={"Cache-Control": "max-age=86400"})


def overlay_tile(req: Request):
    p = req.path_params
    s = review_app._session(p["sid"])
    version, mode, L = int(p["version"]), p["mode"], int(p["level"])
    tx, ty = (int(v) for v in p["xy"].split(".")[0].split("_"))
    only = tuple(sorted(filter(None, req.query_params.get("only", "").split("|"))))
    size = round(float(np.clip(float(req.query_params.get("size", 1.0)), 0.2, 3.0)), 2)  # dot size factor
    v = view(s.sample_id)

    def make():
        x0, y0, x1, y1, w, h, sc = v.tile_box(L, tx, ty)
        col = colors_for(s, version, mode, only)
        # take cells a little beyond the tile so dots on its border are drawn whole
        rad = CELL_RADIUS_UM / v.um_per_px  # in level-0 image pixels
        idx = v.cells_in(x0, y0, x1, y1, pad=max(rad * size, 2.0 * size / sc) + 2)
        idx = idx[col[idx, 3] > 0]
        u, w_ = (v.px[idx] - x0) * sc, (v.py[idx] - y0) * sc
        r = rad * sc * size
        if only or mode in ("conflicts", "changed", "plausibility"):  # few cells shown: keep them findable when zoomed out
            r = max(r, 2.0 * size)
        if r < 1.6:  # zoomed out: one pixel (or 2x2) per cell
            arr = np.zeros((h, w, 4), np.uint8)
            iu, iv = np.floor(u).astype(int), np.floor(w_).astype(int)
            for du in ((0,) if r < 1.1 else (0, 1)):
                for dv in ((0,) if r < 1.1 else (0, 1)):
                    m = (iu + du >= 0) & (iu + du < w) & (iv + dv >= 0) & (iv + dv < h)
                    arr[iv[m] + dv, iu[m] + du] = col[idx[m]]
            img = Image.fromarray(arr, "RGBA")
        else:
            img = Image.new("RGBA", (w, h), (0, 0, 0, 0))
            d = ImageDraw.Draw(img)
            edge = r >= 3.5
            for x, y, c in zip(u, w_, col[idx]):
                d.ellipse([x - r, y - r, x + r, y + r], fill=tuple(int(k) for k in c),
                          outline=(255, 255, 255, 160) if edge else None)
        buf = io.BytesIO()
        img.save(buf, "PNG", compress_level=1)
        return buf.getvalue()
    try:
        data = _overlays.get_or((s.id, version, mode_key(s, mode), only, size, L, tx, ty, s.sample_id, v.rev), make)
    except ValueError:
        return Response(status_code=404)
    return Response(data, media_type="image/png", headers={"Cache-Control": "max-age=86400"})


# ---------------- cells for the page's WebGL layer ----------------
def api_cells_bin(req):
    """Cell centers in viewer image coordinates, float32 x, y pairs (drawn by the page with WebGL)."""
    v = view(req.path_params["sample"])
    data = np.c_[v.px - v.ox, v.py - v.oy].astype(np.float32).tobytes()
    return Response(data, media_type="application/octet-stream",
                    headers={"X-Cells": str(len(v.px)), "X-Radius-Px": f"{CELL_RADIUS_UM / v.um_per_px:.4f}", **NO_CACHE})


def api_colors_bin(req):
    """RGBA per cell (uint8) for a version, color mode and legend filter, in the order of cells.bin."""
    s = review_app._session(req.path_params["sid"])
    q = req.query_params
    only = tuple(sorted(filter(None, q.get("only", "").split("|"))))
    try:
        col = colors_for(s, int(q["version"]), q.get("mode", "type"), only)
    except (KeyError, ValueError) as e:
        return _json({"error": str(e)}, 400)
    return Response(np.ascontiguousarray(col).tobytes(), media_type="application/octet-stream")


# ---------------- API ----------------
NO_CACHE = {"Cache-Control": "no-cache, no-store, must-revalidate"}  # a refresh must load the current page


async def index(req):
    return FileResponse(STATIC / "app.html", headers=NO_CACHE)


async def review_page(req):
    return FileResponse(STATIC / "review.html", headers=NO_CACHE)


def desktop_mode():
    return os.environ.get("STAGENT_DESKTOP") == "1"


def api_config(req):
    """The built-in chat box (Claude Code CLI) is only for the server viewer; the desktop app connects to the user's
    own Claude Code instead."""
    return _json({"chat": not desktop_mode() and shutil.which("claude") is not None, "desktop": desktop_mode(),
                  "home": str(config.HOME), "platform": sys.platform, "build": config.build(),
                  "remote": os.environ.get("STAGENT_REMOTE") == "1", "host": platform.node()})


FS_KINDS = {"h5ad": (".h5ad", ".rds", ".rda", ".rdata"),
            "xenium": (".parquet", ".csv.gz", ".h5", ".xenium", ".mtx.gz", ".tsv.gz", ".ome.tif"),
            "xcells": (".parquet", ".csv.gz"), "xcounts": (".h5", ".mtx.gz", ".tsv.gz"), "xexp": (".xenium",),
            "ximage": (".ome.tif", ".tif", ".tiff"), "image": (".tif", ".tiff", ".btf"), "csv": (".csv", ".tsv", ".txt"),
            "markers": (".csv", ".tsv", ".txt")}


def api_fs_list(req):
    """Folders and files on the computer that runs this server, to choose inputs in the app when it works on a
    remote server (the app's own file dialogs show the user's computer)."""
    if not desktop_mode():
        return _json({"error": "only in the app"}, 403)
    q = req.query_params
    p = Path(q.get("path") or "~").expanduser()
    if p.is_file():
        p = p.parent
    if not p.is_dir():
        return _json({"error": f"no folder {p}"}, 404)
    kind, dirs, files = q.get("kind", ""), [], []
    exts = FS_KINDS.get(kind)
    try:
        entries = sorted(os.scandir(p), key=lambda e: e.name.lower())
    except OSError as e:
        return _json({"error": f"cannot read {p}: {e.strerror}"}, 403)
    for e in entries[:5000]:
        if e.name.startswith("."):
            continue
        try:
            if e.is_dir():
                dirs.append(e.name)
            elif kind != "folder" and (exts is None or e.name.lower().endswith(exts)):
                files.append({"name": e.name, "size": e.stat().st_size})
        except OSError:
            continue
    return _json({"path": str(p), "parent": str(p.parent) if p.parent != p else None, "dirs": dirs, "files": files,
                  "home": str(Path.home())})


# ---------------- viewer <-> agent, Claude connection ----------------
async def api_viewer_state(req: Request):
    """The page reports what it shows; the agent's use_viewer_session reads it."""
    b = await req.json()
    claudelink.write_viewer_state(sample=b.get("sample"), session=b.get("session"), version=b.get("version"))
    return _json({"ok": True})


def api_viewer_request(req):
    """A session the agent asked the viewer to show (open_viewer)."""
    return _json(claudelink.viewer_request())


def api_claude_status(req):
    return _json(claudelink.status())


REVIEW_REQUEST = ("Check the review candidates of version {v} (I clicked Agent Review: Annotation): look at every "
                  "candidate area on the H&E with show_area and judge every marker check, then record each verdict "
                  "with record_review, keeping only real annotation problems.")
SCOPE_REQUEST = {
    "annotation": ("review the cell types of version {v}: gather sample_context, type_confusion, type_structure, "
                   "review_markers, qc_overview and check_registration; advice in the categories taxonomy, context, qc"),
    "organization": ("review the spatial organization of version {v}: gather show_tissue_overview, preview_compartment, "
                     "spatial_profile, list_regions and composition(rows='region'); advice in the categories "
                     "compartments (outlines checked on the H&E) and organization"),
    "embedding": ("review the PCA + UMAP and its Leiden clusters: gather embedding_facts and show_embedding (by cell "
                  "type, leiden and, with several slices, slice); advice in the categories clusters, integration, "
                  "parameters"),
    "niches": ("review the niches: gather niche_summary, show_niches, group_markers(by='niche'), "
               "composition(rows='niche', cols='region') and niche_context; advice in the categories niches and "
               "parameters"),
    "compare": ("review the comparisons I ran in the Compare tab: recent_analyses(section='compare'), then rerun "
                "them in full (de_genes, group_markers, composition); advice in the categories findings, design, next"),
    "signaling": ("review the cell-cell communication I ran in the Signaling tab: recent_analyses(section='signaling'), "
                  "then cell_communication / compare_communication and show_communication for the striking pairs; "
                  "advice in the categories interactions, design, next"),
}
SCOPE_LEAD = "I clicked Agent Review: {label}. Review only this part: {what}. "
SCOPE_CHECK = ("Draft your advice with draft_advice(scope='{scope}'), then check every draft as a skeptical reviewer "
               "(rerun the tool that would confirm or refute it, look at the figure, weigh other explanations) and "
               "record check_advice: keep only what holds up, with the evidence; drop the rest with the reason.")
REVIEW_END = " Do not change the annotation or the analyses. Then sum up in three lines."
REPORT_REQUEST = ("Write the annotation report of version {v} (I clicked Report): gather annotation_facts, "
                  "session_summary, sample_context, list_advice, check_plausibility and review_markers, and save it "
                  "with write_report. Do not change the annotation.")


def _ask_claude(s, text):
    """Send a request to the user's Claude Code (desktop app) or hand it to the server viewer's chat box."""
    if not desktop_mode():
        return {"message": text}
    return {"claude_code": claudelink.open_claude_code(s.sample_id, s.id, text, send_if_open=True)}


def review_scopes(s):
    """The parts of the analysis the agent can review, whether each has something to review, and its status."""
    from . import analysis, joint
    emb = embedding.current(s.sample_id)
    pid = (emb or {}).get("project")
    n_comp, n_sig = len(analysis.recent(s, "compare")), len(analysis.recent(s, "signaling"))
    avail = {"annotation": None, "organization": None,
             "embedding": None if emb else "run PCA + UMAP first (Embed tab)",
             "niches": None if niches.current(s)[0] else "run the niches first (Niches tab)",
             "compare": None if n_comp else "run a comparison first (Compare tab)",
             "signaling": None if n_sig else "run the cell-cell communication first (Signaling tab)"}
    what = {"annotation": "cell types and markers: the computed review candidates, types to merge, split or add, "
                          "sample context, technical problems",
            "organization": "compartments on the H&E (proposed as regions), gradients and neighbors",
            "embedding": "Leiden clusters against the cell types" + (", slices mixing (joint embedding)" if pid else "")
                         + ", parameters",
            "niches": "what each niche is, names, k and number of niches",
            "compare": f"the {n_comp} comparison(s) you ran: findings, caveats, next analyses",
            "signaling": f"the {n_sig} communication analysis(es) you ran: interactions, caveats, next analyses"}
    return [{"scope": k, "label": v["label"], "categories": v["categories"], "what": what[k], "unavailable": avail[k],
             "status": advice.status(s, k)} for k, v in advice.SCOPES.items()]


async def api_claude_review(req: Request):
    """Agent Review of one part of the analysis (scope). Annotation: the candidates the agent has not judged yet and,
    when the labels changed since its last advice, advice on the cell types. Other scopes: always asked."""
    b = await req.json()
    s = review_app._session(b["session"])
    scope = b.get("scope", "annotation")
    if scope not in advice.SCOPES:
        return _json({"error": f"unknown scope {scope!r}"}, 400)
    v = s.state()["version"]
    label = advice.SCOPES[scope]["label"]
    parts = []
    if scope == "annotation":
        st = claudereview.status(s)
        want_advice = b.get("advice", True) and not advice.status(s, scope)["current"]
        if not st["unchecked"] and not want_advice:
            return _json({"needed": False, **st, "advice": advice.status(s, scope)})
        if st["unchecked"]:
            claudereview.request(s)
            parts.append(REVIEW_REQUEST.format(v=v))
        if want_advice:
            advice.request(s, scope)
            parts.append(f"Then {SCOPE_REQUEST[scope].format(v=v)}. {SCOPE_CHECK.format(scope=scope)}" if parts else
                         SCOPE_LEAD.format(label=label, what=SCOPE_REQUEST[scope].format(v=v)) + SCOPE_CHECK.format(scope=scope))
    else:
        advice.request(s, scope)
        parts.append(SCOPE_LEAD.format(label=label, what=SCOPE_REQUEST[scope].format(v=v)) + SCOPE_CHECK.format(scope=scope))
    text = " ".join(parts) + REVIEW_END
    return _json({"needed": True, **_ask_claude(s, text), **claudereview.status(s), "advice": advice.status(s, scope)})


NICHE_K_REQUEST = ("Choose the neighborhood size k for the niche analysis of version {v} (I clicked 'Claude picks "
                   "k'): look at niche_context (cell density, the distance k neighbors cover and the lineages per "
                   "region) and the tissue, choose k per region (or radius_um) so each neighborhood spans the tissue "
                   "scale that matters there (with no regions drawn: radius_um adapts to density, or one k; you may "
                   "offer region outlines with propose_region for a later run), run_niches with it ({n} niches, "
                   "identity {identity}), then look at "
                   "niche_summary, name the niches with name_niches and explain your choice of k in two lines.")


async def api_claude_niches(req: Request):
    b = await req.json()
    s = review_app._session(b["session"])
    return _json(_ask_claude(s, NICHE_K_REQUEST.format(v=s.state()["version"], n=int(b.get("n_niches", 8)),
                                                       identity=b.get("identity", "nmf"))))


async def api_niches_run(req: Request):
    """Run the niche analysis in the background; progress under /api/job/niche:<session>."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    job = f"niche:{s.id}"
    JOBS[job] = {"stage": "starting", "progress": 0.0, "done": False, "error": None}

    def work():
        try:
            niches.run(s, k=int(b.get("k", 15)), n_niches=int(b.get("n_niches", 8)), identity=b.get("identity", "nmf"),
                       n_programs=int(b.get("n_programs", 20)),
                       progress=lambda f, m: JOBS[job].update(progress=round(f, 3), stage=m))
            JOBS[job].update(progress=1.0, stage="done", done=True)
        except Exception as e:
            JOBS[job].update(error=f"{type(e).__name__}: {e}", done=True)
    threading.Thread(target=work, daemon=True).start()
    return _json({"job": job})


async def api_embedding_run(req: Request):
    """PCA, UMAP and Leiden of a sample in the background; progress under /api/job/embed:<sample>."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    job = f"embed:{s.sample_id}"
    if job in JOBS and not JOBS[job]["done"]:
        return _json({"error": "already running"}, 409)
    JOBS[job] = {"stage": "starting", "progress": 0.0, "done": False, "error": None}

    def work():
        try:
            embedding.run(s, n_pcs=int(b.get("n_pcs", 30)), n_neighbors=int(b.get("n_neighbors", 15)),
                          min_dist=float(b.get("min_dist", 0.3)), resolution=float(b.get("resolution", 1.0)),
                          progress=lambda f, m: JOBS[job].update(progress=round(f, 3), stage=m))
            JOBS[job].update(progress=1.0, stage="done", done=True)
        except Exception as e:
            JOBS[job].update(error=f"{type(e).__name__}: {e}", done=True)
    threading.Thread(target=work, daemon=True).start()
    return _json({"job": job})


def api_embedding(req):
    return _json({"embedding": embedding.current(req.path_params["sample"])})


def api_embedding_bin(req):
    """Embedding coordinates (float32 x, y pairs in cell order)."""
    try:
        Y = embedding.coords(req.path_params["sample"], req.query_params.get("kind", "umap"))
    except (ValueError, OSError) as e:
        return _json({"error": str(e)}, 404)
    return Response(np.ascontiguousarray(Y, np.float32).tobytes(), media_type="application/octet-stream", headers=NO_CACHE)


def api_niches(req):
    s = review_app._session(req.path_params["sid"])
    return _json({"summary": niches.summary(s)})


async def api_claude_report(req: Request):
    b = await req.json()
    s = review_app._session(b["session"])
    return _json(_ask_claude(s, REPORT_REQUEST.format(v=s.state()["version"])))


# ---------------- high-level advice and report ----------------
def api_advice(req):
    """Claude's high-level advice (outlines in viewer image coordinates), its status, and the report's date."""
    s = review_app._session(req.path_params["sid"])
    v = view(s.sample_id)
    out = []
    for x in advice.items(s):
        if x["status"] == "draft":  # the agent has not checked it yet
            continue
        a = dict(x["action"] or {})
        if "polygons_px" in a:
            a["outlines"] = [_viewer_poly(v, p) for p in a.pop("polygons_px")]
        a.pop("correction_px", None)
        out.append(x | {"action": a or None})
    rep = advice.report(s)
    return _json({"items": out, "categories": advice.CATEGORIES, "status": advice.status(s), "scopes": review_scopes(s), "signature": advice.signature(s),
                  "report": rep and {k: rep[k] for k in ("version", "time") if k in rep},
                  "registration_fixed": bool(s.meta.get("registration_fixes"))})


def _reload_sample(sample_id):
    """Drop cached sessions and views of a sample whose cell positions changed (view() and _session() also notice
    it by themselves, from the file's modification time)."""
    view(sample_id)
    for sid in [k for k, x in review_app._sessions.items() if x.sample_id == sample_id]:
        review_app._sessions.pop(sid, None)
    if mcp_server._active.get("session") is not None and mcp_server._active["session"].sample_id == sample_id:
        mcp_server._active["session"] = Session(mcp_server._active["session"].id)


async def api_advice_act(req: Request):
    """The user takes a piece of advice: step region (add the outline as a region), prior_preview / prior (preview
    or add its prior), registration / undo_registration, claude (send its request to Claude Code)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        x = advice.get(s, b.get("id", ""))
        a = x["action"] or {}
        step = b.get("step")
        if step == "region":
            name = str(b.get("name") or a["name"]).strip()
            ver = s.add_region_pieces(name, a["polygons_px"], author="Claude advice")
            if not a.get("prior"):
                advice.set_status(s, x["id"], "done", f"region {name}")
            return _json({"new_version": ver, "message": s.state(ver)["message"], "prior": bool(a.get("prior"))})
        if step in ("prior_preview", "prior"):
            spec = dict(a["prior"], author="Claude advice (accepted in the viewer)")
            if step == "prior_preview":
                p = s.preview_constraint(spec)
                return _json({k: p[k] for k in ("interpretation", "n_label_changes", "top_transitions")})
            ver = s.add_constraint(spec)
            advice.set_status(s, x["id"], "done", "prior added")
            return _json({"new_version": ver, "message": s.state(ver)["message"]})
        if step == "registration":
            r = importer.apply_registration_fix(s.sample_id, a["correction_px"])
            advice.set_status(s, x["id"], "done", "registration corrected")
            _reload_sample(s.sample_id)
            return _json(r | {"reload": True})
        if step == "undo_registration":
            r = importer.undo_registration_fix(s.sample_id)
            advice.set_status(s, x["id"], "open")
            _reload_sample(s.sample_id)
            return _json(r | {"reload": True})
        if step == "claude":
            return _json(_ask_claude(s, a["request"]))
        if step == "niche_names":
            names = niches.set_names(s, a["names"])
            advice.set_status(s, x["id"], "done", "niches named")
            return _json({"names": names, "message": "Niches named"})
        if step in ("show", "rerun"):  # done in the page; noted so the agent sees it was used
            advice.set_status(s, x["id"], "done", "shown on the map" if step == "show" else f"settings put in the {a['tab']} tab")
            return _json({"ok": True})
        raise ValueError(f"unknown step {step!r}")
    except (ValueError, KeyError) as e:
        return _json({"error": str(e)}, 400)


async def api_advice_status(req: Request):
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        advice.set_status(s, b.get("id", ""), b.get("status", ""))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"ok": True})


def api_report(req):
    s = review_app._session(req.path_params["sid"])
    return _json({"report": advice.report(s)})


def api_report_md(req):
    s = review_app._session(req.path_params["sid"])
    rep = advice.report(s)
    if not rep:
        return _json({"error": "no report yet"}, 404)
    return FileResponse(rep["path"], media_type="text/markdown", filename=f"{s.sample_id}_annotation_report.md")


async def api_report_save(req: Request):
    """Write the report to a path the user chose (desktop app: Save dialog)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    rep = advice.report(s)
    try:
        if not rep:
            raise ValueError("no report yet")
        shutil.copy(rep["path"], Path(b["path"]).expanduser())
    except (OSError, KeyError, ValueError) as e:
        return _json({"error": str(e)}, 400)
    return _json({"saved": b["path"]})


async def api_claude_code(req: Request):
    b = await req.json()
    return _json(claudelink.open_claude_code(b.get("sample"), b.get("session"), (b.get("request") or "").strip() or None))


# ---------------- import ----------------
async def api_import_inspect(req: Request):
    b = await req.json()
    try:
        if b.get("kind") == "xenium":
            return _json(importer.inspect_xenium_many(b.get("entries") or [b["folder"]]))
        paths = [p for p in b.get("h5ad", []) if p]
        if not paths:
            raise ValueError("choose one or more h5ad or RDS files")
        return _json(importer.inspect_h5ad(paths))
    except (ValueError, OSError, KeyError) as e:
        return _json({"error": str(e)}, 400)


async def api_import(req: Request):
    """Import a sample in the background; progress under /api/job/import:<name>."""
    b = await req.json()
    name = str(b.get("name", "")).strip()
    job = f"import:{name}"
    JOBS[job] = {"stage": "starting", "progress": 0.0, "done": False, "error": None, "log": []}

    def say(msg):
        JOBS[job]["log"].append(msg)
        JOBS[job].update(stage=msg, progress=min(0.95, JOBS[job]["progress"] + 0.15))

    sections = [s for s in (b.get("sections") or []) if s]
    folders = [f for f in (b.get("folders") or []) if f]  # Xenium: folders and / or files

    def work():
        try:
            n_xenium = len(importer.xenium_sources(folders)) if b.get("kind") == "xenium" and folders else 0
            if len(sections) > 1 or n_xenium > 1:  # several slices: a project, one sample per slice
                importer.import_project(name, b.get("kind", "h5ad"), b.get("h5ad") if b.get("kind") == "h5ad" else folders,
                                        sections=sections, section_column=b.get("section_column") or None,
                                        he=b.get("he") or None, registration=b.get("registration") or None,
                                        he_kind=b.get("he_kind") or "H&E", region_column=b.get("region_column") or None,
                                        description=b.get("description", ""), overwrite=bool(b.get("overwrite")),
                                        progress=say)
            elif b.get("kind") == "xenium":
                importer.import_xenium(name, importer.xenium_sources(folders)[0] if folders else b["folder"], he=b.get("he") or None, registration=b.get("registration") or None,
                                       description=b.get("description", ""), overwrite=bool(b.get("overwrite")), progress=say,
                                       he_kind=b.get("he_kind") or "H&E", registration_section=b.get("registration_section") or None)
            else:
                importer.import_sample(name, b["h5ad"], b.get("he") or None, b.get("registration") or None, b.get("section_column") or None,
                                       b.get("section") or (sections[0] if sections else None),
                                       region_column=b.get("region_column") or None,
                                       description=b.get("description", ""), overwrite=bool(b.get("overwrite")), progress=say,
                                       he_kind=b.get("he_kind") or "H&E", registration_section=b.get("registration_section") or None)
            JOBS[job].update(progress=1.0, stage="done", done=True)
        except Exception as e:  # shown in the dialog
            JOBS[job].update(error=f"{type(e).__name__}: {e}", done=True)
    threading.Thread(target=work, daemon=True).start()
    return _json({"job": job})


def api_samples(req):
    t = list_samples()
    ss = list_sessions()
    rows = []
    for r in t.to_dict("records"):
        meta = json.load(open(config.SAMPLES / r["sample_id"] / "meta.json"))
        mine = ss[ss.sample_id == r["sample_id"]] if len(ss) else ss
        rows.append(r | {"bundle_dir": meta.get("bundle_dir"), "source_dir": meta.get("source_dir"),
                         "he_path": meta.get("he_path"), "image_kind": meta.get("image_kind"),
                         "image_pending": meta.get("image_pending"), "registration_by": meta.get("registration_by"),
                         "has_registration": bool(meta.get("registration_affine")),
                         "can_undo_placement": bool(meta.get("placement_previous")),
                         "project": meta.get("project"), "name": meta.get("name"), "metadata": meta.get("metadata") or {},
                         "sessions": mine.sort_values("created_at", ascending=False)[["session_id", "created_at", "n_versions", "method"]]
                         .to_dict("records") if len(mine) else []})
    return _json(rows)


def _project_job(pid, fn):
    """Run a project-wide analysis in the background; progress under /api/job/project:<id>."""
    job = f"project:{pid}"
    if job in JOBS and not JOBS[job]["done"]:
        return _json({"error": "an analysis of this project is still running"}, 409)
    JOBS[job] = {"stage": "starting", "progress": 0.0, "done": False, "error": None}

    def work():
        try:
            fn(lambda f, m: JOBS[job].update(progress=round(f, 3), stage=m))
            JOBS[job].update(progress=1.0, stage="done", done=True)
        except Exception as e:
            JOBS[job].update(error=f"{type(e).__name__}: {e}", done=True)
    threading.Thread(target=work, daemon=True).start()
    return _json({"job": job})


async def api_project_annotate(req: Request):
    from . import joint
    pid, b = req.path_params["pid"], await req.json()
    return _project_job(pid, lambda say: joint.annotate(pid, b.get("markers_text", ""), b.get("resolution"),
                                                        bool(b.get("carry")), say))


async def api_project_embed(req: Request):
    from . import joint
    pid, b = req.path_params["pid"], await req.json()
    return _project_job(pid, lambda say: joint.embed(pid, int(b.get("n_pcs", 30)), int(b.get("n_neighbors", 15)),
                                                     float(b.get("min_dist", 0.3)), float(b.get("resolution", 1.0)),
                                                     bool(b.get("harmony", True)), say))


async def api_project_niches(req: Request):
    from . import joint
    pid, b = req.path_params["pid"], await req.json()
    return _project_job(pid, lambda say: joint.niches(pid, int(b.get("k", 15)), int(b.get("n_niches", 8)),
                                                      b.get("identity", "nmf"), int(b.get("n_programs", 20)), say))


# ---------------- differential expression and composition ----------------
ANALYSES = LRU(8)   # recent results, for their CSV


_cell_tables = LRU(2)  # the last cell tables (dot plots of one result ask for them again and again)


def _analysis_cells(b):
    """The cell table of an analysis request: {session, version, scope ("slice" or "project"), selection}; the
    selection ({sample_id: rows}) is also kept for the agent's tools."""
    s = review_app._session(b["session"])
    sel = b.get("selection") or None
    claudelink.write_viewer_selection(sel)
    key = (s.id, b.get("version"), b.get("scope"), json.dumps(sel, sort_keys=True) if sel else "",
           (niches.current(s)[0] or {}).get("run_id"), (embedding.current(s.sample_id) or {}).get("run_id"), s.versions()[-1])
    return _cell_tables.get_or(key, lambda: _make_cells(s, b, sel))


def _make_cells(s, b, sel):
    from . import analysis
    if b.get("scope") == "project":
        from . import projects
        pid = projects.project_of(s.sample_id)
        if not pid:
            raise ValueError("this sample is not part of a project")
        return analysis.cells(pid=pid, selection=sel, get_session=review_app._session)
    v = b.get("version")
    return analysis.cells(s, None if v is None else int(v), selection=(sel or {}).get(s.sample_id) if sel else None)


async def api_analysis_options(req: Request):
    b = await req.json()
    try:
        c = await asyncio.to_thread(_analysis_cells, b)
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"groupings": c.options(), "missing": c.missing, "n_cells": len(c.labels)})


async def api_analysis_run(req: Request):
    """kind "de" (by, group, reference, within), "markers" (by, within) or "composition" (rows, cols, within)."""
    from . import analysis
    b = await req.json()

    def work():
        c = _analysis_cells(b)
        within = {k: v for k, v in (b.get("within") or {}).items() if v not in (None, "", [])}
        if b["kind"] == "de":
            r = analysis.de(c, b["by"], b["group"], b.get("reference") or "rest", within, top=60, full=True)
            table = r.pop("table")
        elif b["kind"] == "markers":
            r = analysis.markers(c, b["by"], top=15, within=within)
            table = pd.DataFrame([{"group": g["group"], "n_cells": g["n_cells"], **x} for g in r["groups"] for x in g["genes"]])
        elif b["kind"] == "composition":
            r = analysis.crosstab(c, b["rows"], b["cols"], within)
            table = pd.DataFrame([{b["rows"]: x["row"], "n_cells": x["n_cells"], **x["shares"]} for x in r["table"]]).fillna(0)
        else:
            raise ValueError("kind must be de, markers or composition")
        key = uuid.uuid4().hex[:12]
        ANALYSES[key] = table
        while len(ANALYSES) > ANALYSES.n:
            ANALYSES.popitem(last=False)
        analysis.log_run(review_app._session(b["session"]), b, r)
        return {**r, "id": key}
    try:
        return _json(await asyncio.to_thread(work))
    except (ValueError, KeyError) as e:
        return _json({"error": str(e).strip("'")}, 400)


async def api_communication_run(req: Request):
    """Cell-cell communication in the background (kind "run" or "compare"); progress and result under
    /api/job/ccc:<session>."""
    b = await req.json()
    job = f"ccc:{b['session']}"
    if job in JOBS and not JOBS[job]["done"]:
        return _json({"error": "a communication analysis of this annotation is still running"}, 409)
    JOBS[job] = {"stage": "starting", "progress": 0.0, "done": False, "error": None}

    def work():
        try:
            c = _analysis_cells(b)
            say = lambda f, m: JOBS[job].update(progress=round(f, 3), stage=m)
            within = {k: v for k, v in (b.get("within") or {}).items() if v not in (None, "", [])}
            kw = dict(by=b.get("by") or "cell type", within=within, resource=b.get("resource") or "all",
                      contact_um=float(b.get("contact_um") or communication.CONTACT_UM),
                      secreted_um=float(b.get("secreted_um") or communication.SECRETED_UM),
                      min_pct=float(b.get("min_pct", 10)) / 100, nonprotein=bool(b.get("nonprotein")))
            if b.get("kind") == "compare":
                r = communication.compare(c, b["split"], top=80, **kw, progress=say)
            else:
                r = communication.run(c, top=80, **kw, progress=say)
            table = r.pop("table")
            key = uuid.uuid4().hex[:12]
            ANALYSES[key] = table
            while len(ANALYSES) > ANALYSES.n:
                ANALYSES.popitem(last=False)
            from . import analysis
            analysis.log_run(review_app._session(b["session"]), b, r)
            JOBS[job].update(progress=1.0, stage="done", done=True, result={**r, "id": key, "kind": b.get("kind", "run")})
        except Exception as e:
            JOBS[job].update(error=f"{type(e).__name__}: {e}" if not isinstance(e, (ValueError, KeyError)) else str(e).strip("'"),
                             done=True)
    threading.Thread(target=work, daemon=True).start()
    return _json({"job": job})


ANALYSIS_ASK = {
    "de": ("Interpret the DE genes I found in the app: de_genes(by={by}, group={group}, reference={reference}, "
           "within={within}, scope={scope}). Say in a few lines what sets the group apart, which genes point to which "
           "cell states, programs or neighbors (possibly spillover), how consistent it is across slices, and what to "
           "check next."),
    "markers": ("Interpret the top genes per group I found in the app: group_markers(by={by}, within={within}, "
                "scope={scope}). Say in a few lines what each group seems to be and which groups look alike."),
    "composition": ("Interpret the composition I looked at in the app: composition(rows={rows}, cols={cols}, "
                    "within={within}, scope={scope}). Say in a few lines which differences stand out (enrichment), "
                    "whether they hold in every slice, and what they may mean for the tissue."),
    "run": ("Interpret the cell-cell communication I found in the app: cell_communication(by={by}, within={within}, "
            "scope={scope}, resource={resource}). Say in a few lines which signaling axes and pathways stand out between "
            "which groups, which fit the tissue, which may be artifacts (spillover, broadly expressed genes), and what to "
            "check next (e.g. show_communication on the H&E)."),
    "compare": ("Interpret the comparison of cell-cell communication I made in the app: compare_communication("
                "split={split}, by={by}, within={within}, scope={scope}, resource={resource}). Say in a few lines which "
                "interactions differ between the groups, how much to trust it (replicate slices), and what it may mean."),
}
ANALYSIS_ASK_END = " Separate what the data show from what general knowledge suggests. Do not change the annotation."


async def api_claude_analysis(req: Request):
    b = await req.json()
    s = review_app._session(b["session"])
    q = b["request"]
    args = {k: json.dumps(q.get(k), ensure_ascii=False) for k in ("by", "group", "reference", "rows", "cols", "scope",
                                                                  "split", "resource")}
    args["within"] = json.dumps(q.get("within") or None, ensure_ascii=False)
    return _json(_ask_claude(s, ANALYSIS_ASK[q["kind"]].format(**args) + ANALYSIS_ASK_END))


async def api_analysis_dots(req: Request):
    """Dot plot data of genes over the groups of a grouping: {session, version, scope, selection, by, genes, within,
    first}."""
    from . import analysis
    b = await req.json()

    def work():
        c = _analysis_cells(b)
        within = {k: v for k, v in (b.get("within") or {}).items() if v not in (None, "", [])}
        return analysis.dots(c, b["by"], list(b["genes"])[:200], within, tuple(b.get("first") or ()))
    try:
        return _json(await asyncio.to_thread(work))
    except (ValueError, KeyError) as e:
        return _json({"error": str(e).strip("'")}, 400)


def api_analysis_csv(req):
    t = ANALYSES.get(req.path_params["key"])
    if t is None:
        return _json({"error": "this result is gone: run it again"}, 404)
    return Response(t.to_csv(index=False), media_type="text/csv",
                    headers={"Content-Disposition": f'attachment; filename="{req.query_params.get("name", "analysis.csv")}"'})


def api_communication_bubble(req):
    t = ANALYSES.get(req.path_params["key"])
    if t is None or "significant" not in t:
        return _json({"error": "this result is gone: run it again"}, 404)
    n = max(1, min(60, int(req.query_params.get("n", 20))))
    return _json(communication.bubble(t, n))


async def api_save_text(req: Request):
    """Desktop app: write a text file the page made to the path chosen in the Save dialog."""
    if not desktop_mode():
        return _json({"error": "only in the desktop app"}, 403)
    b = await req.json()
    try:
        Path(b["path"]).expanduser().write_text(b.get("text", ""))
    except (OSError, KeyError) as e:
        return _json({"error": str(e)}, 400)
    return _json({"saved": b["path"]})


def api_projects(req):
    from . import projects
    return _json(projects.list_projects())


def api_he_info(req):
    v = view(req.path_params["sample"])
    return _json({"width": v.W, "height": v.H, "tile_size": TILE, "origin": [v.ox, v.oy], "um_per_px": v.um_per_px,
                  "image_kind": v.image_kind, "rev": int(v.rev)})


async def api_annotate(req: Request):
    b = await req.json()
    sample, mode = str(b["sample_id"]), b.get("mode", "new")
    if mode == "continue" and b.get("session_id"):
        sid = b["session_id"]
    else:
        try:
            sid = Session.create(sample, b.get("region_source", "none"), method=b.get("method", "markers"),
                                 markers_text=b.get("markers_text"), resolution=b.get("resolution"),
                                 carry_from=b.get("carry_from"), reference=b.get("reference"),
                                 existing=b.get("existing")).id
        except (ValueError, OSError, KeyError) as e:
            return _json({"error": str(e)}, 400)
    s = review_app._session(sid)
    JOBS[sid] = {"stage": "starting", "progress": 0.0, "done": False, "error": None}

    def work():
        def cb(f, msg):
            JOBS[sid].update(progress=round(0.9 * f, 3), stage=msg)
        try:
            s.base(progress=cb)
            JOBS[sid].update(progress=0.93, stage="assigning regions, applying priors, finding conflicts")
            s.result()
            view(s.sample_id)
            JOBS[sid].update(progress=1.0, stage="done", done=True)
        except Exception as e:  # report to the page instead of dying silently
            JOBS[sid].update(error=f"{type(e).__name__}: {e}", done=True)
    threading.Thread(target=work, daemon=True).start()
    return _json({"session_id": sid})


def api_job(req):
    return _json(JOBS.get(req.path_params["sid"], {"stage": "idle", "progress": 1.0, "done": True, "error": None}))


def api_legend(req):
    s = review_app._session(req.path_params["sid"])
    v = int(req.query_params.get("version") or s.versions()[-1])
    mode = req.query_params.get("mode", "type")
    return _json({"version": v, "mode": mode, "items": legend_for(s, v, mode)})


def api_cell_at(req):
    s = review_app._session(req.path_params["sid"])
    q = req.query_params
    v = view(s.sample_id)
    x, y = float(q["x"]) + v.ox, float(q["y"]) + v.oy
    d, row = v.tree.query([x, y])
    if d > max(float(q.get("tol", 12)), 12):
        return _json({"cell": None})
    info = s.cell_info(int(row), int(q["version"]) if q.get("version") else None)
    if q.get("gene"):  # the gene the map is colored by
        i = {x.upper(): k for k, x in enumerate(s.panel)}.get(q["gene"].upper())
        if i is not None:
            info["gene"] = {"name": s.panel[i], "count": int(s.counts[int(row), i])}
    return _json(info | {"x_img": float(v.px[row] - v.ox), "y_img": float(v.py[row] - v.oy)})


def api_regions(req):
    s = review_app._session(req.path_params["sid"])
    v = view(s.sample_id)
    ver = int(req.query_params["version"]) if req.query_params.get("version") else None
    info = s.regions_info(ver)
    colors = region_rgb(s, ver)
    for p in info["polygons"]:
        p["polygon"] = [[x - v.ox, y - v.oy] for x, y in p["polygon"]]  # viewer image coordinates
        p["color"] = colors.get(p["name"], GRAY)
    for r in info["regions"]:
        r["color"] = colors.get(r["name"], GRAY)
    return _json(info)


def _viewer_poly(v, poly):
    return [[x - v.ox, y - v.oy] for x, y in poly]


def api_plausibility(req):
    """Candidate areas and Claude's area suggestions for one version (outlines in viewer image coordinates)."""
    s = review_app._session(req.path_params["sid"])
    v = view(s.sample_id)
    ver = int(req.query_params["version"]) if req.query_params.get("version") else None
    out = plausibility.summary(s, ver, top=int(req.query_params.get("top", 60)))
    full = {a["id"]: a for a in plausibility.scan(s, ver)["areas"]}
    area_items, check_items = claudereview.items(s, ver)
    verdict = {a["id"]: x for a, x in area_items}
    for a in out["candidate_areas"]:
        a["outline"] = _viewer_poly(v, full[a["id"]]["outline_px"])
        x = verdict.get(a["id"])
        a["claude"] = None if x is None else {"keep": x["keep"], "reason": x["reason"]}
    try:
        flags, unassigned = markerqc.open_checks(s, ver)
        judged = {k["key"]: x for k, x in check_items}
        mk = {"flags": [f | {"key": markerqc.check_key(f), "claude": judged.get(markerqc.check_key(f))} for f in flags],
              "unassigned": unassigned and (unassigned | {"claude": judged.get("unassigned")}),
              "n_solved": len(markerqc.solved(s))}
    except ValueError:  # annotations made without a marker list
        mk = None
    review = claudereview.status(s, ver)
    solved = [{"kind": "area", "index": i, "label": f"{x['lineage']} out of place · tiles {x['tiles']}",
               "by": x.get("by", "user"), "reason": x.get("reason", ""), "time": x["time"]}
              for i, x in enumerate(plausibility.solved_areas(s))]
    solved += [{"kind": "marker", "key": x["key"], "label": "Unassigned cells" if x["key"] == "unassigned"
                else f"{x['node']} · {x['marker']}", "by": x.get("by", "user"), "reason": x.get("reason", ""),
                "time": x["time"]} for x in markerqc.solved(s)] if mk is not None else []
    return _json(out | {"suggestions": _suggestions(s, v), "markers": mk, "proposal": markerqc.proposal(s),
                        "n_areas_solved": len(out["solved_by_user"]), "solved": solved, "review": review,
                        "review_count": len(claudereview.load(s)["areas"]) + len(claudereview.load(s)["markers"])})


def _suggestions(s, v):
    keep = ("id", "title", "reason", "action", "region_name", "tiles", "status", "lineage")
    return [{k: x.get(k) for k in keep} | {"outline": _viewer_poly(v, x["outline_px"])} for x in plausibility.suggestions(s)]


def api_area_suggestions(req):
    """Claude's area suggestions and marker list proposal (polled by the viewer)."""
    s = review_app._session(req.path_params["sid"])
    return _json({"areas": _suggestions(s, view(s.sample_id)), "proposal": markerqc.proposal(s),
                  "solved": [len(plausibility.solved_areas(s)), len(markerqc.solved(s))],
                  "review": len(claudereview.load(s)["areas"]) + len(claudereview.load(s)["markers"]),
                  "advice": advice.signature(s),
                  "report": (advice.report(s) or {}).get("time"),
                  "niche": [(niches.current(s)[0] or {}).get("run_id"), niches.niche_names(s)]})


async def api_marker_check_solve(req: Request):
    """The user marks a marker check as solved: it leaves the Review panel and the agent sees it as solved."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        rec = markerqc.solve(s, b.get("key", ""), b.get("version"))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json(rec)


async def api_area_solve(req: Request):
    """The user marks a candidate area as solved: it leaves the Review panel and the map; the agent sees it."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        rec = plausibility.solve_area(s, b.get("area_id", ""), b.get("version"))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json(rec)


async def api_reopen(req: Request):
    """Undo a 'solved' mark on an area (index) or a marker check (key)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        if b.get("kind") == "area":
            plausibility.reopen_area(s, b.get("index", -1))
        else:
            markerqc.reopen(s, b.get("key", ""))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"ok": True})


async def api_proposal_status(req: Request):
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        markerqc.set_proposal_status(s, b.get("status", ""))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"ok": True})


async def api_area_adopt(req: Request):
    """Save a suggestion's or a candidate area's outline as region `name` (one new version)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        if b.get("suggestion_id"):
            sg = next((x for x in plausibility.suggestions(s) if x["id"] == b["suggestion_id"]), None)
            if sg is None:
                raise ValueError(f"no suggestion {b['suggestion_id']}")
            poly = sg["outline_px"]
        else:
            poly = plausibility.area_info(s, b.get("area_id", ""), b.get("version"))["outline_px"]
        ver = s.add_region(b.get("name", ""), poly, author=b.get("author") or "viewer")
        if b.get("suggestion_id"):
            plausibility.set_suggestion_status(s, b["suggestion_id"], "done", region_name=s.state(ver)["regions"]["ops"][-1]["name"])
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver, "message": s.state(ver)["message"]})


async def api_area_dismiss(req: Request):
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        plausibility.set_suggestion_status(s, b.get("suggestion_id", ""), "dismissed")
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"ok": True})


async def api_region_add(req: Request):
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    v = view(s.sample_id)
    try:
        ver = s.add_region(b.get("name", ""), [[x + v.ox, y + v.oy] for x, y in b.get("polygon", [])],
                           author=b.get("author") or "viewer")
    except ValueError as e:  # VocabError is a ValueError
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver, "message": s.state(ver)["message"]})


async def api_region_edit(req: Request):
    """Save the edited outline of one drawn polygon (viewer image coordinates)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    v = view(s.sample_id)
    try:
        ver = s.edit_region(req.path_params["rid"], [[x + v.ox, y + v.oy] for x, y in b.get("polygon", [])],
                            author=b.get("author") or "viewer")
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver, "message": s.state(ver)["message"]})


def api_genes(req):
    return _json({"genes": panel_genes(req.path_params["sample"])})


async def api_constraint_delete(req: Request):
    """Remove one prior (a new version; Undo brings it back)."""
    s = review_app._session(req.path_params["sid"])
    try:
        ver = s.remove_constraint(req.path_params["cid"], message=f"remove {req.path_params['cid']} (in the app)")
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver, "message": s.state(ver)["message"]})


def api_diagnostics(req):
    """About / diagnostics: what runs where, whether Claude Code is found and connected, the data folder."""
    import importlib.metadata as md, platform
    def ver(pkg):
        try:
            return md.version(pkg)
        except Exception:
            return None
    n_samples = sum(1 for d in config.SAMPLES.iterdir() if (d / "meta.json").exists()) if config.SAMPLES.exists() else 0
    n_sessions = sum(1 for d in config.SESSIONS.iterdir() if (d / "session.json").exists()) if config.SESSIONS.exists() else 0
    size = sum(f.stat().st_size for f in config.HOME.rglob("*") if f.is_file()) if config.HOME.exists() else 0
    return _json({"build": config.build(), "app": "desktop app" if desktop_mode() else "web viewer",
                  "standalone": bool(config.FROZEN), "system": f"{platform.system()} {platform.release()} ({platform.machine()})",
                  "python": platform.python_version(), "data_folder": str(config.HOME),
                  "n_samples": n_samples, "n_annotations": n_sessions, "data_size_mb": round(size / 2 ** 20),
                  "claude_code": claudelink.claude_code_path(), "claude_workspace": str(claudelink.workspace()),
                  "claude_connection": claudelink.mcp_seen(),
                  "libraries": {k: ver(k) for k in ("numpy", "scipy", "scikit-learn", "igraph", "pywebview", "mcp")}})


async def api_region_color(req: Request):
    """The user picks a region's color (or null for the automatic one); maps colored by region are redrawn."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        s.set_region_color(str(b.get("name", "")), b.get("color"))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    for cache in (_overlays, _colors):
        for k in [k for k in cache if k[0] == s.id and k[2] == "region"]:
            cache.pop(k)
    return _json({"ok": True})


async def api_region_rename(req: Request):
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        ver = s.rename_region(b.get("old", ""), b.get("new", ""))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver})


async def api_region_delete(req: Request):
    s = review_app._session(req.path_params["sid"])
    try:
        ver = s.delete_region(req.path_params["rid"])  # priors that use it keep applying to its cells
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver, "message": s.state(ver)["message"]})


# Claude Code keeps chat transcripts per project directory; only those listed in a session's chat metadata are touched
CLAUDE_PROJECT = Path.home() / ".claude" / "projects" / str(config.ROOT).replace("/", "-")


def _sample_sessions(sample):
    out = []
    for d in sorted(config.SESSIONS.iterdir()) if config.SESSIONS.exists() else []:
        if (d / "session.json").exists() and json.load(open(d / "session.json")).get("sample_id") == sample:
            out.append(d)
    return out


async def api_sample_reset(req: Request):
    """Remove every annotation and chat of one sample: its sessions and their Claude Code transcripts are moved to
    work/trash/<time>_<sample>/ (recoverable by hand, invisible to the app and the agent)."""
    sample = req.path_params["sample"]
    try:
        with_embedding = bool((await req.json()).get("embedding"))
    except ValueError:
        with_embedding = False
    sess = _sample_sessions(sample)
    busy = [d.name for d in sess if d.name in JOBS and not JOBS[d.name]["done"]]
    if busy:
        return _json({"error": "an annotation of this sample is still running; wait until it finishes"}, 409)
    emb = config.SAMPLES / sample / "embedding"
    if not sess and not (with_embedding and emb.exists()):
        return _json({"deleted_sessions": 0, "deleted_conversations": 0, "deleted_embedding": False, "trash": None})
    trash = config.TRASH / f"{time.strftime('%Y%m%d-%H%M%S')}_{sample}"
    trash.mkdir(parents=True)
    n_chats = _trash_sessions(sess, trash)
    if with_embedding and emb.exists():  # PCA / UMAP / Leiden of the sample (expression only)
        shutil.move(str(emb), str(trash / "embedding"))
    return _json({"deleted_sessions": len(sess), "deleted_conversations": n_chats,
                  "deleted_embedding": with_embedding and (trash / "embedding").exists(), "trash": str(trash)})


def _trash_sessions(sess, trash):
    """Move sessions and their Claude Code transcripts into trash (None: delete them for good); forget them in
    memory. Returns the number of transcripts moved or deleted."""
    n_chats = 0
    for d in sess:
        ids = set()
        for f in [d / "chat_meta.json", *(d / "chat_archive").glob("*_chat_meta.json")]:
            if f.exists():
                ids.add(json.load(open(f)).get("claude_session"))
        for cid in filter(None, ids):
            for t in (CLAUDE_PROJECT / f"{cid}.jsonl", CLAUDE_PROJECT / cid):
                if t.exists():
                    n_chats += t.suffix == ".jsonl"
                    if trash is None:
                        shutil.rmtree(t) if t.is_dir() else t.unlink()
                        continue
                    (trash / "claude_transcripts").mkdir(exist_ok=True)
                    shutil.move(str(t), str(trash / "claude_transcripts" / t.name))
        shutil.rmtree(d) if trash is None else shutil.move(str(d), str(trash / d.name))
        review_app._sessions.pop(d.name, None)
        JOBS.pop(d.name, None)
        if mcp_server._active.get("session") is not None and mcp_server._active["session"].id == d.name:
            mcp_server._active.pop("session")
        for cache in (_overlays, _colors):
            for k in [k for k in cache if k[0] == d.name]:
                cache.pop(k)
    return n_chats


async def api_sample_image(req: Request):
    """Add or replace a sample's image; it is shown once there is a registration (file or alignment)."""
    sample = req.path_params["sample"]
    b = await req.json()
    try:
        r = importer.set_image(sample, b["path"], b.get("kind") or "H&E", progress=lambda *a: None)
    except (ValueError, OSError, KeyError) as e:
        return _json({"error": str(e)}, 400)
    _reload_sample(sample)
    return _json(r)


async def api_sample_registration(req: Request):
    """Add a registration file; with an image the cells move onto it at once."""
    sample = req.path_params["sample"]
    b = await req.json()
    try:
        r = importer.set_registration(sample, b["path"], b.get("section") or None, progress=lambda *a: None)
    except (ValueError, OSError, KeyError) as e:
        return _json({"error": str(e)}, 400)
    _reload_sample(sample)
    return _json(r)


async def api_sample_rename(req: Request):
    b = await req.json()
    try:
        r = projects.rename(req.path_params["sample"], b.get("name"), b.get("project_name"))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    _cell_tables.clear()  # slice names are labels of the analyses
    return _json(r)


async def api_sample_placement_undo(req: Request):
    sample = req.path_params["sample"]
    try:
        r = importer.undo_placement(sample)
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    _reload_sample(sample)
    return _json(r)


ALIGN_REQUEST = ("Align the cells of sample {sample} to its image (I clicked 'Align with Claude'; there is no "
                 "registration file): alignment_status, search_alignment, compare the best candidates on their figures "
                 "(show_alignment), refine_alignment the right one, check the close-ups and numbers, then "
                 "apply_alignment. Do not change the annotation. Then tell me in three lines what you chose and why.")


async def api_claude_align(req: Request):
    b = await req.json()
    sample = b["sample"]
    text = ALIGN_REQUEST.format(sample=sample)
    if not desktop_mode():
        return _json({"message": text})
    return _json({"claude_code": claudelink.open_claude_code(sample, b.get("session"), text, send_if_open=True)})


def api_samples_clear(req):
    """Desktop app: delete every sample with all its annotations and chats for good, to import from scratch (and the
    *_all_samples folders an earlier version of this button left in the trash); marker lists are kept."""
    if not desktop_mode():
        return _json({"error": "only in the desktop app"}, 403)
    if any(not j["done"] for j in JOBS.values()):
        return _json({"error": "an import, annotation or niche run is still running; wait until it finishes"}, 409)
    samples = [d for d in sorted(config.SAMPLES.iterdir()) if d.is_dir()] if config.SAMPLES.exists() else []
    sess = [d for d in sorted(config.SESSIONS.iterdir()) if (d / "session.json").exists()] if config.SESSIONS.exists() else []
    old = sorted(config.TRASH.glob("*_all_samples")) if config.TRASH.exists() else []
    if not samples and not sess and not old:
        return _json({"deleted_samples": 0, "deleted_sessions": 0})
    n_chats = _trash_sessions(sess, None)
    for d in samples + old:
        shutil.rmtree(d)
    for cache in (_views, _tiles, _overlays, _colors):
        cache.clear()
    review_app._sessions.clear()
    mcp_server._active.pop("session", None)
    return _json({"deleted_samples": len(samples), "deleted_sessions": len(sess), "deleted_conversations": n_chats})


async def api_undo(req):
    s = review_app._session(req.path_params["sid"])
    try:
        v = s.undo()
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": v, "message": s.state(v)["message"]})


async def api_export(req):
    return _json({"folder": review_app._session(req.path_params["sid"]).export()})


def api_export_csv(req):
    """Download the per-cell table of a version (default: latest); a copy is written to the session's exports."""
    s = review_app._session(req.path_params["sid"])
    v = int(req.query_params["version"]) if req.query_params.get("version") else s.versions()[-1]
    folder = Path(s.export(v))
    name = f"{s.sample_id}_v{v}_cells.csv"
    return FileResponse(folder / "cells.csv", media_type="text/csv", filename=name)


def api_registration_csv(req):
    """The registration in use (file, Claude's alignment or a correction) as an a11..ty CSV."""
    from . import notebook
    sample = req.path_params["sample"]
    return Response(notebook.registration_csv(sample), media_type="text/csv",
                    headers={"Content-Disposition": f'attachment; filename="{sample}_registration.csv"'})


def api_export_notebook(req):
    """Download the analysis notebook of a version (default: latest)."""
    from . import notebook
    s = review_app._session(req.path_params["sid"])
    v = int(req.query_params["version"]) if req.query_params.get("version") else s.versions()[-1]
    return FileResponse(notebook.build(s, v), media_type="application/x-ipynb+json",
                        filename=f"{s.sample_id}_v{v}_analysis.ipynb")


def reply_language(text):
    """The language to answer in, decided here instead of by the model: with an empty history the model sometimes
    picked a third language (German) for an English message."""
    if re.search(r"[\u4e00-\u9fff]", text or ""):
        return "Chinese (Simplified)"
    latin = len(re.findall(r"[A-Za-z]", text or ""))
    return "English" if latin >= 0.5 * max(1, len(re.sub(r"\s", "", text or ""))) else None


def language_line(text):
    lang = reply_language(text)
    return f" Reply in {lang}." if lang else " Reply in the language of the user's message."


# ---------------- marker lists ----------------
LIBRARY = config.USER_MARKERS  # marker lists the user imported or saved (there is no built-in default)


def api_panel(req):
    return _json({"genes": panel_genes(req.query_params.get("sample") or None)})


async def api_markers_check(req: Request):
    b = await req.json()
    p = markers.parse(b.get("text", ""), panel_genes(b.get("sample") or None))
    return _json({k: p[k] for k in ("ok", "errors", "warnings", "n_types", "depth", "genes", "canonical")})


def api_markers_organs(req):
    from . import organs
    return _json({"organs": organs.organs()})


async def api_markers_organ_draft(req: Request):
    """The built-in list of an organ pruned to the sample's panel: {organ, sample}."""
    from . import organs
    b = await req.json()
    try:
        return _json(organs.draft(b["organ"], panel_genes(b.get("sample") or None)))
    except (ValueError, KeyError) as e:
        return _json({"error": str(e).strip("'")}, 400)


ORGAN_REQUEST = ("Draft a marker list for {label} for sample {sample} (I chose this organ in Markers…). Start from "
                 "organ_marker_draft(organ='{organ}', sample_id='{sample}'): the built-in list pruned to the panel, with "
                 "the types it had to leave out. {extra}Then add the cell types specific to this tissue that the panel "
                 "can resolve and the draft lacks, and remove types the panel cannot separate (get_panel_genes). Use the "
                 "built-in markers where they fit; for other types search the literature online for their markers when "
                 "you are not sure. Validate with check_marker_list, save the list with save_marker_list(name="
                 "'{name}', sample_id='{sample}'), and tell me in three lines what you changed; I load it in Markers… "
                 "(load a list…).")


async def api_claude_markers(req: Request):
    """Ask the agent to draft the marker list of an organ for a sample, from the built-in list."""
    from . import organs
    b = await req.json()
    sample, organ = b.get("sample") or "", (b.get("organ") or "").strip()
    if not sample or not organ:
        return _json({"error": "choose a sample and an organ"}, 400)
    built = {o["organ"]: o["label"] for o in organs.organs()}
    if organ in built and organ != "other":
        label, key, extra = built[organ], organ, ""
    else:
        label, key = organ, "other"
        extra = (f"There is no built-in list for {organ}: that draft holds only immune and stromal cells; add the "
                 f"organ's own (parenchymal) cell types and their markers from the literature (search online). ")
    text = ORGAN_REQUEST.format(label=label, organ=key, sample=sample, extra=extra,
                                name=f"{organs.safe_name(organ)}_{sample}"[:60])
    if not desktop_mode():  # the server viewer drafts it itself (Markers… > Draft), without web search
        return _json({"message": f"Organ: {label}. Start from organ_marker_draft(organ='{key}', sample_id='{sample}') "
                                 "and complete it for this tissue: add the cell types specific to it that the panel can "
                                 "resolve, and remove those it cannot separate."})
    return _json({"claude_code": claudelink.open_claude_code(sample, None, text, send_if_open=True)})


def _library_files():
    return sorted([f for f in LIBRARY.glob("*") if f.suffix in (".txt", ".csv")]) if LIBRARY.exists() else []


def api_markers_library(req):
    return _json([{"name": f.stem, "where": f.suffix[1:]} for f in _library_files()])


def _as_text(raw):
    """Marker list in the editor's indented text format (CSV files are converted)."""
    return markers.csv_to_text(raw) if markers.is_csv(raw) else raw


def api_markers_get(req):
    name = req.path_params["name"]
    f = next((f for f in _library_files() if f.stem == name), None)
    if f is None:
        return _json({"error": "not found"}, 404)
    return _json({"name": name, "text": _as_text(f.read_text())})


def _safe_name(name):
    return "".join(ch for ch in name if ch.isalnum() or ch in "-_ .").strip().replace(" ", "_")


async def api_markers_save(req: Request):
    b = await req.json()
    name = _safe_name(b.get("name", ""))
    if not name:
        return _json({"error": "give the list a name"}, 400)
    LIBRARY.mkdir(parents=True, exist_ok=True)
    (LIBRARY / f"{name}.csv").write_text(markers.text_to_csv(b.get("text", "")))
    return _json({"name": name})


async def api_markers_import(req: Request):
    """Import a marker list file (CSV or text) from a path or its content; a copy goes to the library."""
    b = await req.json()
    try:
        if b.get("path"):
            src = Path(b["path"]).expanduser()
            raw, name = src.read_text(), src.stem
        else:
            raw, name = b.get("text", ""), Path(b.get("name") or "imported").stem
        text = _as_text(raw)
    except (OSError, UnicodeDecodeError, ValueError) as e:
        return _json({"error": str(e)}, 400)
    name = _safe_name(name) or "imported"
    LIBRARY.mkdir(parents=True, exist_ok=True)
    (LIBRARY / f"{name}.csv").write_text(markers.text_to_csv(text) if markers.parse(text, None)["tree"] else raw)
    return _json({"name": name, "text": text})


async def api_reference_inspect(req: Request):
    """A single-cell reference h5ad: cells, label columns, gene overlap with the sample's panel."""
    from . import scref
    b = await req.json()
    try:
        return _json(scref.inspect(Path(b["path"]).expanduser(), panel_genes(b.get("sample"))))
    except (OSError, KeyError, ValueError) as e:
        return _json({"error": f"cannot read the reference: {e}"}, 400)


async def api_existing_inspect(req: Request):
    """An annotation CSV (by path, or uploaded text kept in <data folder>/annotations/): columns and matches."""
    from . import existing
    from .session import sample_cell_ids
    b = await req.json()
    try:
        if b.get("text") is not None:
            folder = config.HOME / "annotations"
            folder.mkdir(parents=True, exist_ok=True)
            path = folder / (_safe_name(Path(b.get("name") or "annotation.csv").name) or "annotation.csv")
            path.write_text(b["text"])
        else:
            path = Path(b["path"]).expanduser()
        return _json(existing.inspect(path, sample_cell_ids(b["sample"])))
    except (OSError, KeyError, ValueError, UnicodeDecodeError) as e:
        return _json({"error": f"cannot read the annotation: {e}"}, 400)


async def api_markers_to_csv(req: Request):
    b = await req.json()
    return Response(markers.text_to_csv(b.get("text", "")), media_type="text/csv")


async def api_markers_save_file(req: Request):
    """Write the marker list as CSV to a path the user chose (desktop app: Save dialog)."""
    b = await req.json()
    try:
        Path(b["path"]).expanduser().write_text(markers.text_to_csv(b.get("text", "")))
    except (OSError, KeyError) as e:
        return _json({"error": str(e)}, 400)
    return _json({"saved": b["path"]})


async def api_export_to(req: Request):
    """Write the per-cell CSV (kind "csv") or the analysis notebook (kind "notebook") of a version to a path the user
    chose (desktop app: Save dialog)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    v = int(b.get("version") or s.versions()[-1])
    try:
        from . import notebook
        if b.get("kind") == "registration":
            Path(b["path"]).expanduser().write_text(notebook.registration_csv(s.sample_id))
        else:
            src = notebook.build(s, v) if b.get("kind") == "notebook" else Path(s.export(v)) / "cells.csv"
            shutil.copy(src, Path(b["path"]).expanduser())
    except (OSError, KeyError) as e:
        return _json({"error": str(e)}, 400)
    return _json({"saved": b["path"]})


DRAFT_PROMPT = """You draft a hierarchical marker gene list for cell typing of an imaging-based spatial
transcriptomics sample (fixed gene panel). Call get_panel_genes first and use ONLY genes from the panel.
Format, one type per line, children indented by two spaces under their parent (any depth):
  Name: GENE+ GENE- A|B|C+
GENE+ = expressed, GENE- = not expressed (use these to separate siblings), A|B|C+ = any one of these genes
(useful for lineage level). Markers on one line are averaged, so give 2-4 specific markers per type.
Only include types the panel can resolve. Names must be unique. Cell types are the most specific entries: a type
with subtypes is only a group, and every one of its cells will get one of its subtypes, so the subtypes must cover
all of its cells (add a catch-all subtype, e.g. "Effector/memory CD4 T cells", where needed). Validate with check_marker_list and fix
every error before answering. End your answer with the final list in a fenced block that starts with
```markers and ends with ```. Before the block, write at most three short sentences on the main choices and
any types the panel cannot resolve, in the language named in the request (type names stay as in the list)."""


async def api_markers_draft(req: Request):
    if desktop_mode():
        return _json({"error": "ask Claude in Claude Code to draft the list"}, 400)
    b = await req.json()
    instr = (b.get("instructions") or "").strip() or "Draft a complete list for this tissue."
    current = (b.get("current") or "").strip()
    prompt = f"Tissue / context: {b.get('context') or 'see the request'}.{language_line(b.get('instructions'))}\nRequest: {instr}"
    if current:
        prompt += f"\nStart from this existing list and revise it:\n```\n{current}\n```"
    args = ["claude", "-p", prompt, "--output-format", "stream-json", "--verbose", "--tools", "",
            "--mcp-config", str(CHAT_MCP), "--strict-mcp-config",
            "--allowedTools", "mcp__stagent__get_panel_genes mcp__stagent__check_marker_list mcp__stagent__organ_marker_draft",
            "--append-system-prompt", DRAFT_PROMPT]

    async def stream():
        proc = await asyncio.create_subprocess_exec(*args, cwd=str(config.ROOT), stdout=asyncio.subprocess.PIPE,
                                                    stderr=asyncio.subprocess.PIPE, limit=2 ** 24)
        texts = []
        async for line in proc.stdout:
            try:
                e = json.loads(line)
            except json.JSONDecodeError:
                continue
            if e.get("type") == "assistant":
                for c in e["message"].get("content", []):
                    if c.get("type") == "text" and c["text"].strip():
                        texts.append(c["text"])
                    elif c.get("type") == "tool_use":
                        yield json.dumps({"kind": "tool", "name": c["name"].split("__")[-1]}) + "\n"
        await proc.wait()
        full = "\n".join(texts)
        block = full.split("```markers", 1)
        if len(block) == 2:
            yield json.dumps({"kind": "text", "text": block[0].strip()}) + "\n"
            yield json.dumps({"kind": "markers", "text": block[1].split("```", 1)[0].strip("\n")}) + "\n"
        else:
            err = (await proc.stderr.read()).decode(errors="replace")[-600:]
            yield json.dumps({"kind": "error", "text": full[-800:] or err or "no marker list returned"}) + "\n"
        yield json.dumps({"kind": "done"}) + "\n"
    return StreamingResponse(stream(), media_type="application/x-ndjson")


# ---------------- chat bridge to Claude Code ----------------
SYSTEM_PROMPT = """You are the annotation assistant inside the stagent web viewer.
The user looks at their sample on the H&E with cells colored by type and types into a chat box.
The active session is {sid} (sample {sample}). It is already loaded: never call start_session,
resume_session or list_sessions. Use only the stagent tools.

Reply in the language the context line names (it is detected from the user's message), briefly (2-6
sentences, no tables unless asked); keep cell type and region names as they are. Do not announce what you are
about to do or narrate tool calls; act, then report the outcome. The map refreshes by itself
after every change; do not describe figures you have not seen.

Two kinds of input:
1. A general prior ("epithelial cells are only in the mucosa", "Tfh cells sit next to B cells"): map it to a
   constraint (kind region / depth / neighborhood; names from get_vocabulary), call preview_constraint, tell the
   user in one or two sentences how you understood it and how many labels would change, and ask for a go-ahead.
   Only after the user agrees, call add_constraint with the same arguments and text = the user's words.
2. A correction of cell labels in a region ("the smooth muscle in mucosa_2 is muscularis mucosae", "the cells
   in region X are goblet cells"): apply it directly with label_cells(region=..., cell_type=..., only_types=...
   when only some cells of the region are meant); the instruction itself is the go-ahead. Regions are drawn and
   named by the user with "Draw region" on the map. If the user asks you to draw one, compute the outline
   (preview_compartment, or tiles) and offer it with propose_region: the user adds it from the Review panel. After label_cells, report how many cells you
   labeled and how many other cells changed because the model was updated. If the cells differ mainly by location
   (e.g. muscularis mucosae vs muscularis propria smooth muscle), suggest a region or depth prior as well.
Region priors only act on cells inside drawn regions; cells in no region are never penalized. When the user adds
an "only in X" prior, mention which other regions exist so they know where it takes effect.
If the user asks to undo, use revert to the previous version or remove_correction.
Marker step: list_clusters shows the clusters with their marker scores, assignment and Leiden resolution, and
the resolution chosen at each level (it is picked automatically). label_cluster labels a whole cluster when the
user names one ("root/c27 is neurons"). When a cluster is mixed (unassigned with close scores, or the user says it
contains two kinds of cells), when the user wants a type split further, or when the user keeps correcting cells
of one cluster, suggest a local re-clustering: call preview_recluster (target = that cluster id or type;
resolution "auto" or higher for finer; marker_edits to sharpen existing subtypes' markers), report how many marker
labels would change and the new clusters, and on a go-ahead call recluster. Report changes inside and outside the
re-clustered cells separately. Adding or removing cell types needs a revised marker list: the user edits it in
"Markers..." and clicks Annotate with "keep regions, priors and corrections" checked.
Reviewing the annotation (the user asks whether it makes sense): use
your knowledge of tissue biology, marker genes and histology. Check (1) markers with review_markers: do each type's
top genes fit its label, are listed markers weak or unspecific here with a better panel gene available, do
Unassigned cells or unassigned clusters express a type missing from the list, is the composition plausible for the
tissue; (2) space with show_tissue_overview, check_plausibility and show_area: identify the compartments on the H&E
and areas where the annotation contradicts them (vascular smooth muscle, muscularis mucosae, deep crypts or
fistulas in IBD are real). Reply with a short prioritized list of suggestions, each with its evidence and action.
Make them actionable in the viewer: propose_marker_list for a revised marker list (the user opens it from the
Review panel in "Markers..." and annotates again keeping regions and priors), suggest_areas for areas (drawn on the
map; the user adopts them as regions), preview_recluster for a mixed cluster. Once a suggested region exists, apply
its action after a go-ahead: usually an exclude prior on it (preview_constraint, then add_constraint), or
label_cells(region=..., only_types=...) when the user names the target type. Change nothing without agreement.
When asked to check the review candidates (after Review annotation): look at every candidate area on the H&E with
show_area and judge every marker check, then record each verdict with record_review (keep only real problems; cells at
the edge of a sheet, next to mislabeled neighbors or where few cells were segmented are false alarms). Do not change
the annotation while checking.
When asked to review one part of the analysis (annotation, organization, embedding, niches, compare, signaling):
gather its evidence, draft a few items with draft_advice(scope=...), check every draft as a skeptical reviewer and
record check_advice (keep only what holds up, with the evidence). When asked for a report: annotation_facts plus
the review results, saved with write_report.
If a request cannot be expressed with the tools, say so plainly."""

BLOCKED = ["start_session", "resume_session", "list_sessions", "list_samples", "set_region_source"]


async def api_chat(req: Request):
    sid = req.path_params["sid"]
    s = review_app._session(sid)
    b = await req.json()
    msg = (b.get("message") or "").strip()
    if not msg:
        return _json({"error": "empty message"}, 400)
    mcp_server._active["session"] = s
    ctx = [f"[viewer context] session {sid}, current version v{s.versions()[-1]}."]
    reg = s.regions_info()
    if reg["regions"]:
        ctx.append("Regions drawn by the user: " + ", ".join(f"{r['name']} ({r['n_cells']} cells)" for r in reg["regions"])
                   + f"; {reg['n_unknown']} cells are in no region.")
    else:
        ctx.append("The user has not drawn any region yet.")
    prompt = " ".join(ctx) + language_line(msg) + "\n\nUser: " + msg
    meta_p = s.dir / "chat_meta.json"
    meta = json.load(open(meta_p)) if meta_p.exists() else {}
    args = ["claude", "-p", prompt, "--output-format", "stream-json", "--verbose", "--tools", "",
            "--mcp-config", str(CHAT_MCP), "--strict-mcp-config", "--allowedTools", "mcp__stagent",
            "--disallowedTools", " ".join(f"mcp__stagent__{t}" for t in BLOCKED),
            "--append-system-prompt", SYSTEM_PROMPT.format(sid=sid, sample=s.sample_id)]
    if meta.get("claude_session"):
        args += ["--resume", meta["claude_session"]]
    else:
        meta["claude_session"] = str(uuid.uuid4())
        args += ["--session-id", meta["claude_session"]]
        json.dump(meta, open(meta_p, "w"))
    s._append("chat.jsonl", {"role": "user", "text": msg, "version": s.versions()[-1]})

    async def stream():
        proc = await asyncio.create_subprocess_exec(*args, cwd=str(config.ROOT), stdout=asyncio.subprocess.PIPE,
                                                    stderr=asyncio.subprocess.PIPE, limit=2 ** 24)
        texts = []
        async for line in proc.stdout:
            try:
                e = json.loads(line)
            except json.JSONDecodeError:
                continue
            if e.get("type") == "assistant":
                for c in e["message"].get("content", []):
                    if c.get("type") == "text" and c["text"].strip():
                        texts.append(c["text"])
                        yield json.dumps({"kind": "text", "text": c["text"]}) + "\n"
                    elif c.get("type") == "tool_use":
                        yield json.dumps({"kind": "tool", "name": c["name"].split("__")[-1]}) + "\n"
            elif e.get("type") == "result" and e.get("is_error"):
                yield json.dumps({"kind": "error", "text": str(e.get("result"))}) + "\n"
        err = (await proc.stderr.read()).decode(errors="replace")
        await proc.wait()
        if proc.returncode and not texts:
            yield json.dumps({"kind": "error", "text": err[-800:] or f"claude exited with {proc.returncode}"}) + "\n"
        s._append("chat.jsonl", {"role": "assistant", "text": "\n\n".join(texts), "version": s.versions()[-1]})
        yield json.dumps({"kind": "done", "version": s.versions()[-1]}) + "\n"
    return StreamingResponse(stream(), media_type="application/x-ndjson")


def api_chat_log(req):
    if not (config.SESSIONS / req.path_params["sid"] / "session.json").exists():
        return _json([])  # an annotation that was just removed
    p = review_app._session(req.path_params["sid"]).dir / "chat.jsonl"
    return _json([json.loads(l) for l in open(p)] if p.exists() else [])


def build_app(port):
    claudelink.write_viewer_state(port=port, pid=os.getpid(), mode="desktop" if desktop_mode() else "web")
    CHAT_MCP.write_text(json.dumps({"mcpServers": {"stagent": {"type": "http", "url": f"http://127.0.0.1:{port}/mcp"}}}))
    mcp_http = mcp_server.mcp.streamable_http_app()

    @contextlib.asynccontextmanager
    async def lifespan(app):
        async with mcp_server.mcp.session_manager.run():
            yield

    routes = [
        Route("/", index), Route("/review", review_page),
        Mount("/static", StaticFiles(directory=STATIC)),
        Route("/api/samples", api_samples),
        Route("/api/projects", api_projects),
        Route("/api/save_text", api_save_text, methods=["POST"]),
        Route("/api/project/{pid}/annotate", api_project_annotate, methods=["POST"]),
        Route("/api/project/{pid}/embed", api_project_embed, methods=["POST"]),
        Route("/api/project/{pid}/niches", api_project_niches, methods=["POST"]),
        Route("/api/analysis/options", api_analysis_options, methods=["POST"]),
        Route("/api/analysis/run", api_analysis_run, methods=["POST"]),
        Route("/api/analysis/dots", api_analysis_dots, methods=["POST"]),
        Route("/api/analysis/{key}.csv", api_analysis_csv),
        Route("/api/communication/{key}/bubble", api_communication_bubble),
        Route("/api/claude/analysis", api_claude_analysis, methods=["POST"]),
        Route("/api/communication/run", api_communication_run, methods=["POST"]),
        Route("/api/config", api_config),
        Route("/api/viewer/state", api_viewer_state, methods=["POST"]),
        Route("/api/viewer/request", api_viewer_request),
        Route("/api/claude/status", api_claude_status),
        Route("/api/claude/code", api_claude_code, methods=["POST"]),
        Route("/api/claude/review", api_claude_review, methods=["POST"]),
        Route("/api/claude/report", api_claude_report, methods=["POST"]),
        Route("/api/claude/niches", api_claude_niches, methods=["POST"]),
        Route("/api/s/{sid}/niches", api_niches),
        Route("/api/s/{sid}/embedding/run", api_embedding_run, methods=["POST"]),
        Route("/api/sample/{sample}/embedding", api_embedding),
        Route("/api/sample/{sample}/embedding.bin", api_embedding_bin),
        Route("/api/s/{sid}/niches/run", api_niches_run, methods=["POST"]),
        Route("/api/s/{sid}/advice", api_advice),
        Route("/api/s/{sid}/advice/act", api_advice_act, methods=["POST"]),
        Route("/api/s/{sid}/advice/status", api_advice_status, methods=["POST"]),
        Route("/api/s/{sid}/report", api_report),
        Route("/api/s/{sid}/report.md", api_report_md),
        Route("/api/s/{sid}/report/save", api_report_save, methods=["POST"]),
        Route("/api/import/inspect", api_import_inspect, methods=["POST"]),
        Route("/api/fs/list", api_fs_list),
        Route("/api/import", api_import, methods=["POST"]),
        Route("/api/sample/{sample}/he_info", api_he_info),
        Route("/api/sample/{sample}/he/{level:int}/{xy}", he_tile),
        Route("/api/annotate", api_annotate, methods=["POST"]),
        Route("/api/job/{sid}", api_job),
        Route("/api/s/{sid}/legend", api_legend),
        Route("/api/s/{sid}/overlay/{version:int}/{mode}/{level:int}/{xy}", overlay_tile),
        Route("/api/sample/{sample}/cells.bin", api_cells_bin),
        Route("/api/s/{sid}/colors.bin", api_colors_bin),
        Route("/api/s/{sid}/cell_at", api_cell_at),
        Route("/api/s/{sid}/regions", api_regions),
        Route("/api/s/{sid}/plausibility", api_plausibility),
        Route("/api/s/{sid}/plausibility/suggestions", api_area_suggestions),
        Route("/api/s/{sid}/plausibility/adopt", api_area_adopt, methods=["POST"]),
        Route("/api/s/{sid}/plausibility/dismiss", api_area_dismiss, methods=["POST"]),
        Route("/api/s/{sid}/plausibility/proposal", api_proposal_status, methods=["POST"]),
        Route("/api/s/{sid}/markers/solve", api_marker_check_solve, methods=["POST"]),
        Route("/api/s/{sid}/areas/solve", api_area_solve, methods=["POST"]),
        Route("/api/s/{sid}/review/reopen", api_reopen, methods=["POST"]),
        Route("/api/s/{sid}/regions/add", api_region_add, methods=["POST"]),
        Route("/api/s/{sid}/regions/rename", api_region_rename, methods=["POST"]),
        Route("/api/s/{sid}/regions/color", api_region_color, methods=["POST"]),
        Route("/api/s/{sid}/constraint/{cid}/delete", api_constraint_delete, methods=["POST"]),
        Route("/api/sample/{sample}/genes", api_genes),
        Route("/api/diagnostics", api_diagnostics),
        Route("/api/s/{sid}/regions/{rid}/delete", api_region_delete, methods=["POST"]),
        Route("/api/s/{sid}/regions/{rid}/edit", api_region_edit, methods=["POST"]),
        Route("/api/s/{sid}/undo", api_undo, methods=["POST"]),
        Route("/api/sample/{sample}/reset", api_sample_reset, methods=["POST"]),
        Route("/api/samples/clear", api_samples_clear, methods=["POST"]),
        Route("/api/sample/{sample}/image", api_sample_image, methods=["POST"]),
        Route("/api/sample/{sample}/registration", api_sample_registration, methods=["POST"]),
        Route("/api/sample/{sample}/rename", api_sample_rename, methods=["POST"]),
        Route("/api/sample/{sample}/placement/undo", api_sample_placement_undo, methods=["POST"]),
        Route("/api/claude/align", api_claude_align, methods=["POST"]),
        Route("/api/s/{sid}/export", api_export, methods=["POST"]),
        Route("/api/s/{sid}/export.csv", api_export_csv),
        Route("/api/s/{sid}/export.ipynb", api_export_notebook),
        Route("/api/sample/{sample}/registration.csv", api_registration_csv),
        Route("/api/s/{sid}/chat", api_chat, methods=["POST"]),
        Route("/api/s/{sid}/chat_log", api_chat_log),
        Route("/api/panel", api_panel),
        Route("/api/markers/check", api_markers_check, methods=["POST"]),
        Route("/api/markers/library", api_markers_library),
        Route("/api/markers/library/{name}", api_markers_get),
        Route("/api/markers/save", api_markers_save, methods=["POST"]),
        Route("/api/markers/import", api_markers_import, methods=["POST"]),
        Route("/api/markers/save_file", api_markers_save_file, methods=["POST"]),
        Route("/api/markers/to_csv", api_markers_to_csv, methods=["POST"]),
        Route("/api/reference/inspect", api_reference_inspect, methods=["POST"]),
        Route("/api/existing/inspect", api_existing_inspect, methods=["POST"]),
        Route("/api/s/{sid}/export_to", api_export_to, methods=["POST"]),
        Route("/api/markers/draft", api_markers_draft, methods=["POST"]),
        Route("/api/markers/organs", api_markers_organs),
        Route("/api/markers/organ_draft", api_markers_organ_draft, methods=["POST"]),
        Route("/api/claude/markers", api_claude_markers, methods=["POST"]),
    ] + [r for r in review_app.app.routes if r.path != "/"] + list(mcp_http.routes)
    return Starlette(routes=routes, lifespan=lifespan)


def main():
    import uvicorn
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args()
    uvicorn.run(build_app(a.port), host=a.host, port=a.port, log_level="warning")


if __name__ == "__main__":
    main()
