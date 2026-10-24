"""Main web app: pick a sample, annotate it, view the result on the H&E, and talk to the agent.

Run:  python -m celltypeagent.webapp [--port 8765]
Then open http://localhost:8765 (forward the port from the server first).

Everything runs on this server. The chat box calls the locally installed Claude Code CLI in print mode with
the logged-in account (no API key); its only tools are this app's MCP tools, served at /mcp by this same
process, so the agent and the viewer share one session object. The conflict review page is at /review.
"""
import argparse, asyncio, colorsys, contextlib, io, json, math, re, shutil, threading, time, uuid
from collections import OrderedDict
from pathlib import Path
import numpy as np, pandas as pd, tifffile, zarr
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import FileResponse, Response, StreamingResponse
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles
from . import config, vocab, mcp_server, review_app, render
from . import markers, plausibility
from .session import Session, list_samples, list_sessions, panel_genes

STATIC = Path(__file__).parent / "static"
TILE, MARGIN_PX, CELL_RADIUS_PX = 512, 400, 5.0
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
class Slide:
    def __init__(self, path):
        self.path = path
        with tifffile.TiffFile(path) as t:
            self.levels = []
            for lv in t.series[0].levels:
                chw = lv.shape[0] == 3
                h, w = (lv.shape[1], lv.shape[2]) if chw else (lv.shape[0], lv.shape[1])
                self.levels.append((w, h, chw))
        self.W0, self.H0 = self.levels[0][:2]
        self._z, self._lock = {}, threading.Lock()

    def _zarr(self, k):
        if k not in self._z:
            self._z[k] = zarr.open(tifffile.imread(self.path, aszarr=True, level=k), mode="r")
        return self._z[k]

    def read(self, x0, y0, x1, y1, out_w, out_h):
        """Level-0 pixel box -> RGB array of (out_h, out_w), read from the coarsest adequate pyramid level."""
        need = (x1 - x0) / out_w
        k = max(i for i, (w, _, _) in enumerate(self.levels) if i == 0 or self.W0 / w <= need * 1.001)
        w, h, chw = self.levels[k]
        ds = self.W0 / w
        a0, b0 = max(0, int(x0 / ds)), max(0, int(y0 / ds))
        a1, b1 = min(w, int(math.ceil(x1 / ds))), min(h, int(math.ceil(y1 / ds)))
        canvas = Image.new("RGB", (out_w, out_h), (245, 245, 245))
        if a1 <= a0 or b1 <= b0:
            return canvas
        with self._lock:  # tifffile stores share a file handle
            z = self._zarr(k)
            arr = np.asarray(z[:, b0:b1, a0:a1]).transpose(1, 2, 0) if chw else np.asarray(z[b0:b1, a0:a1])
        img = Image.fromarray(np.ascontiguousarray(arr[..., :3]))
        tw, th = max(1, round((a1 - a0) * ds / need)), max(1, round((b1 - b0) * ds / need))
        canvas.paste(img.resize((tw, th), Image.BILINEAR), (round((a0 * ds - x0) / need), round((b0 * ds - y0) / need)))
        return canvas


class SampleView:
    """Virtual image = the sample's bounding box on its H&E slide (the slide holds two sections)."""

    def __init__(self, sample_id):
        meta = json.load(open(config.SAMPLES / sample_id / "meta.json"))
        self.sample_id = sample_id
        self.um_per_px = float(meta.get("he_um_per_px", config.HE_UM_PER_PX))
        self.slide = Slide(meta["he_path"])
        cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
        px, py = cells.px.to_numpy(), cells.py.to_numpy()
        self.px, self.py = px, py
        self.ox, self.oy = max(0, int(px.min()) - MARGIN_PX), max(0, int(py.min()) - MARGIN_PX)
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
    if sample_id not in _views:
        _views[sample_id] = SampleView(sample_id)
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
    hexes = render.region_colors(names)
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
    return _colors.get_or((s.id, version, mode, tuple(only)), make)


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
        return [{"key": k, "label": k, "color": c, "count": n[k]} for k, c in PLAUS_RGB.items()] + \
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
        idx = v.cells_in(x0, y0, x1, y1, pad=max(CELL_RADIUS_PX * size, 2.0 * size / sc) + 2)
        idx = idx[col[idx, 3] > 0]
        u, w_ = (v.px[idx] - x0) * sc, (v.py[idx] - y0) * sc
        r = CELL_RADIUS_PX * sc * size
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
        data = _overlays.get_or((s.id, version, mode, only, size, L, tx, ty), make)
    except ValueError:
        return Response(status_code=404)
    return Response(data, media_type="image/png", headers={"Cache-Control": "max-age=86400"})


# ---------------- API ----------------
NO_CACHE = {"Cache-Control": "no-cache, no-store, must-revalidate"}  # a refresh must load the current page


async def index(req):
    return FileResponse(STATIC / "app.html", headers=NO_CACHE)


async def review_page(req):
    return FileResponse(STATIC / "review.html", headers=NO_CACHE)


def api_config(req):
    """The viewer's chat box needs the Claude Code CLI; without it the user chats in Claude Desktop instead."""
    return _json({"chat": shutil.which("claude") is not None, "home": str(config.HOME)})


def api_samples(req):
    t = list_samples()
    ss = list_sessions()
    rows = []
    for r in t.to_dict("records"):
        meta = json.load(open(config.SAMPLES / r["sample_id"] / "meta.json"))
        mine = ss[ss.sample_id == r["sample_id"]] if len(ss) else ss
        rows.append(r | {"bundle_dir": meta.get("bundle_dir"), "source_dir": meta.get("source_dir"),
                         "he_path": meta.get("he_path"),
                         "sessions": mine.sort_values("created_at", ascending=False)[["session_id", "created_at", "n_versions", "method"]]
                         .to_dict("records") if len(mine) else []})
    return _json(rows)


def api_he_info(req):
    v = view(req.path_params["sample"])
    return _json({"width": v.W, "height": v.H, "tile_size": TILE, "origin": [v.ox, v.oy], "um_per_px": v.um_per_px})


async def api_annotate(req: Request):
    b = await req.json()
    sample, mode = str(b["sample_id"]), b.get("mode", "new")
    if mode == "continue" and b.get("session_id"):
        sid = b["session_id"]
    else:
        try:
            sid = Session.create(sample, b.get("region_source", "none"), method=b.get("method", "markers"),
                                 markers_text=b.get("markers_text"), resolution=b.get("resolution"),
                                 carry_from=b.get("carry_from")).id
        except ValueError as e:
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


def api_plausibility(req):
    """Groups of isolated cells of cohesive types for one version (image coordinates for the viewer)."""
    s = review_app._session(req.path_params["sid"])
    v = view(s.sample_id)
    ver = int(req.query_params["version"]) if req.query_params.get("version") else None
    top = int(req.query_params.get("top", 40))
    res = plausibility.check(s, ver)
    out = plausibility.summary(s, ver, top)
    for g, full in zip(out["groups"], res["groups"]):
        x0, y0, x1, y1 = full["bbox_px"]
        g["bbox"] = [x0 - v.ox, y0 - v.oy, x1 - v.ox, y1 - v.oy]
    return _json(out)


async def api_plausibility_label(req: Request):
    """Accept suggestions: label the cells of groups computed on `version` (refused if the session moved on)."""
    s = review_app._session(req.path_params["sid"])
    b = await req.json()
    try:
        ver = plausibility.label_groups(s, b.get("group_ids", []), b.get("cell_type") or None, version=b.get("version"),
                                        author=b.get("author") or "viewer")
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver, "message": s.state(ver)["message"]})


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
    body = await req.json() if await req.body() else {}
    try:
        op, users = s.region_users(req.path_params["rid"])
        if users and not body.get("remove_priors"):  # the page asks before removing the priors too
            return _json({"error": f"region {op['name']} is used by priors", "region": op["name"],
                          "priors": [{"id": c["id"], "text": c.get("text", "")} for c in users]}, 409)
        ver = s.delete_region(req.path_params["rid"], remove_priors=bool(body.get("remove_priors")))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"new_version": ver})


# Claude Code keeps chat transcripts per project directory; only those listed in a session's chat metadata are touched
CLAUDE_PROJECT = Path.home() / ".claude" / "projects" / str(config.ROOT).replace("/", "-")


def _sample_sessions(sample):
    out = []
    for d in sorted(config.SESSIONS.iterdir()) if config.SESSIONS.exists() else []:
        if (d / "session.json").exists() and json.load(open(d / "session.json")).get("sample_id") == sample:
            out.append(d)
    return out


def api_sample_reset(req):
    """Remove every annotation and chat of one sample: its sessions and their Claude Code transcripts are moved to
    work/trash/<time>_<sample>/ (recoverable by hand, invisible to the app and the agent)."""
    sample = req.path_params["sample"]
    sess = _sample_sessions(sample)
    busy = [d.name for d in sess if d.name in JOBS and not JOBS[d.name]["done"]]
    if busy:
        return _json({"error": "an annotation of this sample is still running; wait until it finishes"}, 409)
    if not sess:
        return _json({"deleted_sessions": 0, "deleted_conversations": 0, "trash": None})
    trash = config.TRASH / f"{time.strftime('%Y%m%d-%H%M%S')}_{sample}"
    trash.mkdir(parents=True)
    n_chats = 0
    for d in sess:
        ids = set()
        for f in [d / "chat_meta.json", *(d / "chat_archive").glob("*_chat_meta.json")]:
            if f.exists():
                ids.add(json.load(open(f)).get("claude_session"))
        for cid in filter(None, ids):
            for t in (CLAUDE_PROJECT / f"{cid}.jsonl", CLAUDE_PROJECT / cid):
                if t.exists():
                    (trash / "claude_transcripts").mkdir(exist_ok=True)
                    shutil.move(str(t), str(trash / "claude_transcripts" / t.name))
                    n_chats += t.suffix == ".jsonl"
        shutil.move(str(d), str(trash / d.name))
        review_app._sessions.pop(d.name, None)
        JOBS.pop(d.name, None)
        if mcp_server._active.get("session") is not None and mcp_server._active["session"].id == d.name:
            mcp_server._active.pop("session")
        for cache in (_overlays, _colors):
            for k in [k for k in cache if k[0] == d.name]:
                cache.pop(k)
    return _json({"deleted_sessions": len(sess), "deleted_conversations": n_chats, "trash": str(trash)})


async def api_undo(req):
    s = review_app._session(req.path_params["sid"])
    cur = s.state()
    if cur["parent"] is None:
        return _json({"error": "nothing to undo"}, 400)
    v = s.revert(cur["parent"], f"undo v{cur['version']} ({cur['message'][:60]})")
    return _json({"new_version": v})


async def api_export(req):
    return _json({"folder": review_app._session(req.path_params["sid"]).export()})


def api_export_csv(req):
    """Download the per-cell table of a version (default: latest); a copy is written to the session's exports."""
    s = review_app._session(req.path_params["sid"])
    v = int(req.query_params["version"]) if req.query_params.get("version") else s.versions()[-1]
    folder = Path(s.export(v))
    name = f"{s.sample_id}_v{v}_cells.csv"
    return FileResponse(folder / "cells.csv", media_type="text/csv", filename=name)


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
LIBRARY_DIRS = [config.PACKAGE / "data", config.USER_MARKERS]  # shipped lists, then lists saved by the user


def api_panel(req):
    return _json({"genes": panel_genes(req.query_params.get("sample") or None)})


async def api_markers_check(req: Request):
    b = await req.json()
    p = markers.parse(b.get("text", ""), panel_genes(b.get("sample") or None))
    return _json({k: p[k] for k in ("ok", "errors", "warnings", "n_types", "depth", "genes", "canonical")})


def api_markers_library(req):
    rows = []
    for d in LIBRARY_DIRS:
        for f in sorted(d.glob("*.txt")) if d.exists() else []:
            rows.append({"name": f.stem, "where": "project" if d == LIBRARY_DIRS[0] else "saved"})
    return _json(rows)


def api_markers_get(req):
    name = req.path_params["name"]
    for d in LIBRARY_DIRS[::-1]:
        f = d / f"{name}.txt"
        if f.exists() and f.resolve().parent == d.resolve():
            return _json({"name": name, "text": f.read_text()})
    return _json({"error": "not found"}, 404)


async def api_markers_save(req: Request):
    b = await req.json()
    name = "".join(ch for ch in b.get("name", "") if ch.isalnum() or ch in "-_ .").strip().replace(" ", "_")
    if not name:
        return _json({"error": "give the list a name"}, 400)
    LIBRARY_DIRS[1].mkdir(parents=True, exist_ok=True)
    (LIBRARY_DIRS[1] / f"{name}.txt").write_text(b.get("text", ""))
    return _json({"name": name})


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
    b = await req.json()
    instr = (b.get("instructions") or "").strip() or "Draft a complete list for this tissue."
    current = (b.get("current") or "").strip()
    prompt = f"Tissue / context: {b.get('context', 'human intestine (IBD)')}.{language_line(b.get('instructions'))}\nRequest: {instr}"
    if current:
        prompt += f"\nStart from this existing list and revise it:\n```\n{current}\n```"
    args = ["claude", "-p", prompt, "--output-format", "stream-json", "--verbose", "--tools", "",
            "--mcp-config", str(CHAT_MCP), "--strict-mcp-config",
            "--allowedTools", "mcp__celltypeagent__get_panel_genes mcp__celltypeagent__check_marker_list",
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
SYSTEM_PROMPT = """You are the annotation assistant inside the celltypeagent web viewer.
The user looks at their sample on the H&E with cells colored by type and types into a chat box.
The active session is {sid} (sample {sample}). It is already loaded: never call start_session,
resume_session or list_sessions. Use only the celltypeagent tools.

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
   named by the user with "Draw region" on the map; you cannot draw or move them. If the user points at an area
   that is not a region yet, ask them to draw and name it first. After label_cells, report how many cells you
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
Plausibility check (the user asks whether the spatial distribution makes sense, or presses "Ask Claude" in the
plausibility panel): call check_plausibility, look at the top groups with show_plausibility_groups (up to 4 per
call), and judge them with biology: epithelium in submucosa or muscle is suspicious in normal tissue but can be real
in IBD (deep crypts, fistulas, ulcer edges); vascular smooth muscle around vessels and the thin muscularis mucosae
are real; B cells and naive T cells outside follicles are normal. "likely mislabel" groups were wrong about half
the time on test data, confident isolated cells much less often. Reply with three short lists by group id: likely
mislabels (with the suggested type), unusual but possibly real (why), unsure (for the user to look at). Relabel only
after the user agrees, with label_plausibility_groups.
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
            "--mcp-config", str(CHAT_MCP), "--strict-mcp-config", "--allowedTools", "mcp__celltypeagent",
            "--disallowedTools", " ".join(f"mcp__celltypeagent__{t}" for t in BLOCKED),
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
    p = review_app._session(req.path_params["sid"]).dir / "chat.jsonl"
    return _json([json.loads(l) for l in open(p)] if p.exists() else [])


def build_app(port):
    CHAT_MCP.write_text(json.dumps({"mcpServers": {"celltypeagent": {"type": "http", "url": f"http://127.0.0.1:{port}/mcp"}}}))
    mcp_http = mcp_server.mcp.streamable_http_app()

    @contextlib.asynccontextmanager
    async def lifespan(app):
        async with mcp_server.mcp.session_manager.run():
            yield

    routes = [
        Route("/", index), Route("/review", review_page),
        Mount("/static", StaticFiles(directory=STATIC)),
        Route("/api/samples", api_samples),
        Route("/api/config", api_config),
        Route("/api/sample/{sample}/he_info", api_he_info),
        Route("/api/sample/{sample}/he/{level:int}/{xy}", he_tile),
        Route("/api/annotate", api_annotate, methods=["POST"]),
        Route("/api/job/{sid}", api_job),
        Route("/api/s/{sid}/legend", api_legend),
        Route("/api/s/{sid}/overlay/{version:int}/{mode}/{level:int}/{xy}", overlay_tile),
        Route("/api/s/{sid}/cell_at", api_cell_at),
        Route("/api/s/{sid}/regions", api_regions),
        Route("/api/s/{sid}/plausibility", api_plausibility),
        Route("/api/s/{sid}/plausibility/label", api_plausibility_label, methods=["POST"]),
        Route("/api/s/{sid}/regions/add", api_region_add, methods=["POST"]),
        Route("/api/s/{sid}/regions/rename", api_region_rename, methods=["POST"]),
        Route("/api/s/{sid}/regions/{rid}/delete", api_region_delete, methods=["POST"]),
        Route("/api/s/{sid}/undo", api_undo, methods=["POST"]),
        Route("/api/sample/{sample}/reset", api_sample_reset, methods=["POST"]),
        Route("/api/s/{sid}/export", api_export, methods=["POST"]),
        Route("/api/s/{sid}/export.csv", api_export_csv),
        Route("/api/s/{sid}/chat", api_chat, methods=["POST"]),
        Route("/api/s/{sid}/chat_log", api_chat_log),
        Route("/api/panel", api_panel),
        Route("/api/markers/check", api_markers_check, methods=["POST"]),
        Route("/api/markers/library", api_markers_library),
        Route("/api/markers/library/{name}", api_markers_get),
        Route("/api/markers/save", api_markers_save, methods=["POST"]),
        Route("/api/markers/draft", api_markers_draft, methods=["POST"]),
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
