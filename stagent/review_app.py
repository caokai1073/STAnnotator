"""Review web app for pathologists: browse conflict groups, inspect H&E, record verdicts.

Run:  python -m stagent.review_app [--host 127.0.0.1] [--port 8765]
Data never leaves the server. Remote reviewers connect through an SSH tunnel:
      ssh -L 8765:localhost:8765 <server>   then open http://localhost:8765
Verdicts go to the session's verdicts.jsonl through Session.record_verdicts, the same API the agent uses.
"""
import argparse, io, json
from collections import OrderedDict
from pathlib import Path
import numpy as np, pandas as pd
from PIL import Image, ImageDraw
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import FileResponse, Response
from starlette.routing import Route
from . import config, vocab, render
from .session import Session, list_sessions

STATIC = Path(__file__).parent / "static"
VIEWS_UM = {"close": 60.0, "context": 150.0, "wide": 400.0}
DISPLAY_PX = 640
_sessions: OrderedDict = OrderedDict()  # the annotations in memory, most recently used last
# an annotation in memory holds its counts, expression model and results (a few hundred MB for 350k cells): keep only
# the last few; the others are read again from disk when they are opened
MAX_OPEN = int(__import__("os").environ.get("STAGENT_OPEN_ANNOTATIONS", 3))
_crops: OrderedDict = OrderedDict()
CROP_CACHE = 400


def _nan_to_none(o):
    """Browsers reject NaN in JSON; metadata such as missing inflammation grades come through as NaN."""
    if isinstance(o, float) and o != o:
        return None
    if isinstance(o, dict):
        return {k: _nan_to_none(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_nan_to_none(v) for v in o]
    return o


def _json(obj, status=200):
    body = json.dumps(_nan_to_none(json.loads(json.dumps(obj, default=_default))))
    return Response(body, status_code=status, media_type="application/json")


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return None if np.isnan(o) else float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _session(sid) -> Session:
    if sid not in _sessions or _sessions[sid].stale():
        _sessions[sid] = Session(sid)
    _sessions.move_to_end(sid)
    if len(_sessions) > MAX_OPEN:
        _forget_old()
    return _sessions[sid]


def _forget_old():
    """Drop the least recently used annotations beyond MAX_OPEN, except one that is computing right now."""
    dropped = False
    for old in list(_sessions)[:-1]:
        if len(_sessions) <= MAX_OPEN:
            break
        s = _sessions[old]
        if s._lock.acquire(blocking=False):  # not in the middle of a computation
            s._lock.release()
            del _sessions[old]
            dropped = True
    if dropped:
        config.release_memory()


def _version(req):
    v = req.query_params.get("version")
    return int(v) if v else None


# ---------------- crops ----------------
def crop_jpeg(s: Session, cell: int, view: str, overlay: bool) -> bytes:
    key = (s.id, s.versions()[-1], cell, view, overlay)
    if key in _crops:
        _crops.move_to_end(key)
        return _crops[key]
    half_um = VIEWS_UM[view]
    upp = s.he_um_per_px
    half_px = half_um / upp
    row = s.row_of(cell)
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    from .importer import image_of
    im, x0, y0 = render.he_crop(image_of(s.sample_id), px[row], py[row], half_px)
    img = Image.fromarray(im).convert("RGB")
    scale = DISPLAY_PX / max(img.size)
    img = img.resize((round(img.width * scale), round(img.height * scale)), Image.LANCZOS)
    d = ImageDraw.Draw(img)
    if overlay:
        m = (px >= x0) & (px < x0 + im.shape[1]) & (py >= y0) & (py < y0 + im.shape[0])
        pal, lin = s.lineage_colors(), s.lineage_labels()
        r = min(4.0, max(2.0, 1.1 * scale))  # small dots so nuclei stay visible
        edge = (255, 255, 255) if r >= 3 else None  # a white rim swamps the color on small dots
        for ln, x, y in zip(lin[m], (px[m] - x0) * scale, (py[m] - y0) * scale):
            d.ellipse([x - r, y - r, x + r, y + r], fill=pal.get(ln, (0, 0, 0)), outline=edge)
    # target marker: four ticks around the cell (dark outline for contrast), leaving the nucleus visible
    cx, cy = (px[row] - x0) * scale, (py[row] - y0) * scale
    g, L = max(10.0, 8 / upp * scale), 18
    for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
        seg = [cx + dx * g, cy + dy * g, cx + dx * (g + L), cy + dy * (g + L)]
        d.line(seg, fill=(0, 0, 0), width=7)
        d.line(seg, fill=(0, 230, 255), width=4)
    # 50 um scale bar on a white box
    bar = 50 / upp * scale
    d.rectangle([6, img.height - 44, 18 + bar, img.height - 8], fill=(255, 255, 255))
    d.rectangle([12, img.height - 22, 12 + bar, img.height - 15], fill=(0, 0, 0))
    d.text((12, img.height - 40), "50 um", fill=(0, 0, 0))
    buf = io.BytesIO()
    img.save(buf, "JPEG", quality=88)
    _crops[key] = buf.getvalue()
    if len(_crops) > CROP_CACHE:
        _crops.popitem(last=False)
    return _crops[key]


# ---------------- API ----------------
async def index(req):
    return FileResponse(STATIC / "review.html")


async def api_sessions(req):
    t = list_sessions()
    rows = []
    for r in t.to_dict("records"):
        meta = json.load(open(config.SAMPLES / r["sample_id"] / "meta.json"))
        rows.append(r | {"disease": meta.get("disease"), "location": meta.get("location")})
    return _json(rows[::-1])


async def api_info(req):
    s = _session(req.path_params["sid"])
    return _json({"summary": s.summary(_version(req)), "history": s.history(), "verdicts": vocab.VERDICTS,
                  "lineage_colors": [{"name": n, "color": c} for n, c in s.lineage_colors().items()]})


async def api_groups(req):
    s = _session(req.path_params["sid"])
    by = tuple(req.query_params.get("by", "constraint,region").split(","))
    reviewer = req.query_params.get("reviewer", "")
    v = _version(req)
    t = s.conflict_summary(by, v)
    cf = s.result(v)["conflicts"]
    vd = s.verdicts(latest=True)
    mine = set(vd[vd.reviewer == reviewer].cell) if len(vd) else set()
    anyone = set(vd.cell) if len(vd) else set()
    rows = []
    for r in t.to_dict("records"):
        cells = set(s.review_list(r["group_id"], v).cell)
        rows.append(r | {"reviewed_by_me": len(cells & mine), "reviewed_by_anyone": len(cells & anyone)})
    return _json({"n_conflicts": int(len(cf)), "groups": rows})


async def api_queue(req):
    s = _session(req.path_params["sid"])
    gid = req.query_params["group"]
    reviewer = req.query_params.get("reviewer", "")
    q = s.review_list(gid, _version(req))
    labels = s.labels(_version(req))
    vd = s.verdicts(latest=True)
    out = []
    for r in q.to_dict("records"):
        mine = vd[(vd.cell == r["cell"]) & (vd.reviewer == reviewer)] if len(vd) else vd
        others = vd[(vd.cell == r["cell"]) & (vd.reviewer != reviewer)] if len(vd) else vd
        out.append({k: r[k] for k in ["cell", "constraint", "region", "category", "base_label", "base_conf",
                                      "status", "d_allowed_um", "nbr_frac"]}
                   | {"label": labels[r["row"]],
                      "my_verdict": mine.verdict.iloc[-1] if len(mine) else None,
                      "my_note": mine.note.iloc[-1] if len(mine) else "",
                      "others": others[["reviewer", "verdict"]].to_dict("records") if len(others) else []})
    return _json(out)


async def api_crop(req):
    s = _session(req.path_params["sid"])
    view = req.query_params.get("view", "context")
    if view not in VIEWS_UM:
        return _json({"error": f"view must be one of {list(VIEWS_UM)}"}, 400)
    data = crop_jpeg(s, int(req.path_params["cell"]), view, req.query_params.get("overlay") == "1")
    return Response(data, media_type="image/jpeg", headers={"Cache-Control": "max-age=3600"})


async def api_verdict(req: Request):
    s = _session(req.path_params["sid"])
    b = await req.json()
    if not b.get("reviewer"):
        return _json({"error": "reviewer name is required"}, 400)
    try:
        n = s.record_verdicts(b["cells"], b["verdict"], b.get("note", ""), b["reviewer"], b.get("version"), b.get("group"))
    except vocab.VocabError as e:
        return _json({"error": str(e)}, 400)
    return _json({"recorded": n})


async def api_vocab(req):
    return _json(_session(req.path_params["sid"]).vocabulary())


def _clean_spec(b):
    keys = ("text", "kind", "types", "regions", "max_depth_um", "neighbor_types", "min_frac", "mode", "weight")
    return {k: b[k] for k in keys if b.get(k) not in (None, "", [])}


async def api_constraint_preview(req: Request):
    s = _session(req.path_params["sid"])
    b = await req.json()
    try:
        return _json(s.preview_constraint(_clean_spec(b)))
    except (vocab.VocabError, ValueError) as e:
        return _json({"error": str(e)}, 400)


async def api_constraint_add(req: Request):
    s = _session(req.path_params["sid"])
    b = await req.json()
    if not b.get("reviewer"):
        return _json({"error": "reviewer name is required"}, 400)
    if not (b.get("text") or "").strip():
        return _json({"error": "please describe the rule in your own words"}, 400)
    spec = _clean_spec(b) | {"author": b["reviewer"]}
    try:
        v = s.add_constraint(spec, f"add rule by {b['reviewer']} (review page): {spec['text']}",
                             expected_parent=b.get("based_on_version"))
    except vocab.VocabError as e:
        return _json({"error": str(e)}, 400)
    except ValueError as e:
        return _json({"error": str(e)}, 409 if "moved on" in str(e) else 400)
    c = s.state(v)["constraints"][-1]
    return _json({"new_version": v, "constraint_id": c["id"], "diff": s.diff(v - 1, v)})


async def api_suggestion_add(req: Request):
    s = _session(req.path_params["sid"])
    b = await req.json()
    if not b.get("reviewer"):
        return _json({"error": "reviewer name is required"}, 400)
    try:
        sid = s.add_suggestion(b.get("text", ""), b["reviewer"], b.get("cell"), b.get("group"), b.get("version"))
    except ValueError as e:
        return _json({"error": str(e)}, 400)
    return _json({"id": sid})


async def api_suggestions(req):
    return _json(_session(req.path_params["sid"]).suggestions(req.query_params.get("status") or None))


async def api_verdicts_csv(req):
    s = _session(req.path_params["sid"])
    v = s.verdicts(latest=True)
    if len(v):
        info = s.cells.set_index("cell").loc[v.cell, ["pathologist_region"]].reset_index(drop=True)
        v = pd.concat([v.reset_index(drop=True), info], axis=1)
    return Response(v.to_csv(index=False), media_type="text/csv",
                    headers={"Content-Disposition": f"attachment; filename=verdicts_{s.id}.csv"})


app = Starlette(routes=[
    Route("/", index),
    Route("/api/sessions", api_sessions),
    Route("/api/s/{sid}/info", api_info),
    Route("/api/s/{sid}/groups", api_groups),
    Route("/api/s/{sid}/queue", api_queue),
    Route("/api/s/{sid}/crop/{cell:int}", api_crop),
    Route("/api/s/{sid}/verdict", api_verdict, methods=["POST"]),
    Route("/api/s/{sid}/verdicts.csv", api_verdicts_csv),
    Route("/api/s/{sid}/vocab", api_vocab),
    Route("/api/s/{sid}/constraint/preview", api_constraint_preview, methods=["POST"]),
    Route("/api/s/{sid}/constraint", api_constraint_add, methods=["POST"]),
    Route("/api/s/{sid}/suggestion", api_suggestion_add, methods=["POST"]),
    Route("/api/s/{sid}/suggestions", api_suggestions),
])


def main():
    import uvicorn
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args()
    uvicorn.run(app, host=a.host, port=a.port, log_level="warning")


if __name__ == "__main__":
    main()
