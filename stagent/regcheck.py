"""The registration checked over the whole tissue, and its local correction (Agent Review: Spatial organization).

windows(): the cells are cut into square windows (WINDOW_UM) over all the tissue; in each, the nuclei on the image
(hematoxylin of the H&E, or DAPI) are cross-correlated with the cell centroids, giving the offset of the nuclei from
the cells (um; insights.registration does this for a sample of windows). A robust affine fit of the offsets is the
part one registration fixes (scale, rotation, shift); what is left are local shifts (tissue stretched, folded or torn
between the run and the staining). Neighboring windows off their nuclei form areas (L1, L2, ...). Not trusted: a
window whose best offset lies on the edge of the search (densely packed nuclei repeat, so a far shift can fit as
well), and a weak one that disagrees with its neighbors; an offset that fits barely better than none counts as none.
A coarse pass (COARSE_WINDOW_UM windows, image read at COARSE_STEP_UM, search COARSE_SEARCH_UM) finds the areas
shifted further than the fine search reaches (there the fine windows find no peak, so they are not trusted); its shift
counts only where its peak is unique (crypts repeat, so a shift by one crypt can fit almost as well) and most fine
windows in it found no offset they trust.
summary(): the areas, and the offsets by structure (focus groups of types: crypts are where the epithelium is,
vessels where the endothelium is; their nuclei have a clear pattern, so a shift shows first there) and by region.
correction(): first the large shifts of the coarse pass, as a smooth field ("pre"); then, with the cells moved by
it (the fine windows they moved in measured again), affine + a smooth field of the local shifts (Gaussian-weighted,
fading out away from the windows), as importer.apply_registration_fix applies it. Nothing here moves a cell.

The offsets of a placement are kept in the sample's folder (regcheck.json) until the cells move or the image changes.
"""
import json, os, re, time
import numpy as np, pandas as pd
from scipy import ndimage
from . import config, render
from .osutil import dump_json

FORMAT = 3
WINDOW_UM = 200.0   # side of a window
SEARCH_UM = 25      # largest offset looked for
MIN_CELLS = 40      # windows with fewer cells are not measured
MIN_CORR = 0.2      # below: the offset is not trusted
GAIN = 0.85         # an offset counts when the correlation at zero is at most this share of the one at the offset
OFF_UM = 5.0        # a window's cells are off their nuclei above this (a window's mean offset is ~1 um when fine)
FIELD_FLOOR_UM = 2.0  # local shifts below this are left alone by the correction
FIELD_SIGMA = 0.5     # the local field's Gaussian, in windows
COARSE_WINDOW_UM = 600.0  # the coarse pass: side of a window,
COARSE_STEP_UM = 2        # the image read at this many um per pixel,
COARSE_SEARCH_UM = 100    # largest offset looked for,
COARSE_MIN_UM = 20.0      # shorter shifts are left to the fine pass,
UNIQUE = 0.9              # and a peak counts when the best other peak (10 um away or more) reaches at most this share of it
CHECK_ROUNDS = 2          # the correction checked: the cells moved by it, their windows measured again. A window whose
STAYED_SHARE = 0.6        # offset stays (still off, at least this share of it, pointing the same way) is a false reading
STAYED_COS = 0.7          # (a repeating structure matched one period off): left alone, its cells stay, and the field is
                          # made again without it. A reading that jumps elsewhere proves nothing: the prediction stands
DEFAULT_FOCUS = {"crypts / epithelium": r"epith|crypt|colonocyte|enterocyte|goblet|stem|transit|paneth",
                 "vessels": r"endoth|vascul|vessel|pericyte",
                 "smooth muscle": r"smooth|muscle|myocyte"}


def xcorr(sl, kind, cx, cy, half, W, R, qx, qy):
    """Cross-correlation (2R+1 square, in steps of the window side / W) of the nuclei in the image window around
    (cx, cy) (half: its half side in image pixels, read at W x W pixels) with the cell centroids qx, qy (image pixels)
    in it; the peak at (R + dy, R + dx) means the nuclei sit dx, dy steps from the cells (1 um when W is the side in
    um)."""
    from .insights import _nuclei
    img = np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, W, W)).astype(np.float32)
    nuc = ndimage.gaussian_filter(_nuclei(img, kind), 2)
    pts = np.zeros((W, W))
    u = ((qx - (cx - half)) / (2 * half) * W).astype(int).clip(0, W - 1)
    v = ((qy - (cy - half)) / (2 * half) * W).astype(int).clip(0, W - 1)
    np.add.at(pts, (v, u), 1)
    pts = ndimage.gaussian_filter(pts, 3)
    nuc, pts = (nuc - nuc.mean()) / (nuc.std() + 1e-9), (pts - pts.mean()) / (pts.std() + 1e-9)
    cc = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(nuc) * np.conj(np.fft.fft2(pts))))) / W ** 2
    c0 = W // 2
    return cc[c0 - R:c0 + R + 1, c0 - R:c0 + R + 1]


def _file(sample_id):
    return config.SAMPLES / sample_id / "regcheck.json"


def _key(sample_id, image):
    st = (config.SAMPLES / sample_id / "cells.parquet").stat()
    return f"{FORMAT}:{st.st_mtime_ns}:{st.st_size}:{image}"


def windows(sample_id, progress=None, wait_s=300):
    """Offsets per window of the sample's current placement (computed once per placement; another process computing
    them is waited for up to wait_s)."""
    from .importer import image_of, _meta
    say = progress or (lambda f, m: None)
    m = _meta(sample_id)
    image, kind = image_of(sample_id), m.get("image_kind", "H&E")
    from . import visium
    if visium.is_bins(m):
        raise ValueError(visium.NO_REGISTRATION_CHECK)
    if not image or kind == "none":
        raise ValueError("no image to check against (none, or unloaded: Load H&E in the Data tab)")
    key, f = _key(sample_id, image), _file(sample_id)
    busy = f.with_suffix(".computing")
    end = time.time() + wait_s
    while busy.exists() and time.time() - busy.stat().st_mtime < 900 and time.time() < end:
        time.sleep(1)
    if f.exists():
        try:
            old = json.load(open(f))
            if old.get("key") == key:
                return old
        except (OSError, ValueError):
            pass
    busy.write_text(str(os.getpid()), encoding="utf-8")
    try:
        cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
        upp = float(m.get("he_um_per_px") or config.HE_UM_PER_PX)
        sl = render.slide(image)
        px, py = cells.px.to_numpy(), cells.py.to_numpy()
        x0, y0 = float(px.min()), float(py.min())
        tile = WINDOW_UM / upp
        nx = int((px.max() - x0) // tile) + 1
        res = {"key": key, "image": image, "kind": kind, "um_per_px": upp, "window_um": WINDOW_UM, "search_um": SEARCH_UM,
               "origin_px": [x0, y0], "tile_px": tile, "nx": nx, "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "columns": ["i", "j", "cx", "cy", "dx_um", "dy_um", "corr", "n_cells", "corr_at_zero"]}
        res["windows"] = _scan(sl, kind, px, py, res, lambda f: say(0.8 * f, "registration: fine windows"))
        ct = COARSE_WINDOW_UM / upp
        res["coarse"] = {"window_um": COARSE_WINDOW_UM, "search_um": COARSE_SEARCH_UM, "tile_px": ct,
                         "origin_px": [x0, y0], "nx": int((px.max() - x0) // ct) + 1,
                         "columns": res["columns"] + ["second"]}
        res["coarse"]["windows"] = _scan(sl, kind, px, py, res["coarse"] | {"um_per_px": upp},
                                         lambda f: say(0.8 + 0.2 * f, "registration: large shifts"), coarse=True)
        dump_json(res, f)
        return res
    finally:
        busy.unlink(missing_ok=True)


def _scan(sl, kind, px, py, grid, say=None, keys=None, coarse=False):
    """The offset rows of the windows of grid ({"origin_px", "tile_px", "nx", "um_per_px"}) holding the cells at px,
    py (only the window ids in keys, when given); coarse: the coarse pass, its rows ending with the share of the
    best peak that the best other peak (at least 10 um away) reaches."""
    (x0, y0), tile, nx = grid["origin_px"], grid["tile_px"], grid["nx"]
    tj, ti = ((px - x0) // tile).astype(int), ((py - y0) // tile).astype(int)
    inside = (tj >= 0) & (tj < nx) & (ti >= 0)
    groups = pd.Series(np.where(inside)[0]).groupby((ti * nx + tj)[inside]).indices
    step = COARSE_STEP_UM if coarse else 1
    W, R = int((COARSE_WINDOW_UM if coarse else WINDOW_UM) / step), int((COARSE_SEARCH_UM if coarse else SEARCH_UM) / step)
    least = MIN_CELLS * (COARSE_WINDOW_UM / WINDOW_UM) ** 2 / 4 if coarse else MIN_CELLS
    half = tile / 2
    out, todo = [], [(k, rows) for k, rows in groups.items() if len(rows) >= least and (keys is None or k in keys)]
    for n, (k, rows) in enumerate(todo):
        i, j = divmod(int(k), nx)
        cx, cy = x0 + (j + 0.5) * tile, y0 + (i + 0.5) * tile
        if cx - half < 0 or cy - half < 0 or cx + half > sl.W0 or cy + half > sl.H0:
            continue
        win = xcorr(sl, kind, cx, cy, half, W, R, px[rows], py[rows])
        dy, dx = np.unravel_index(np.argmax(win), win.shape)
        row = [i, j, round(cx, 1), round(cy, 1), int(dx - R) * step, int(dy - R) * step, round(float(win.max()), 3),
               len(rows), round(float(win[R, R]), 3)]
        if coarse:
            yy, xx = np.mgrid[-R:R + 1, -R:R + 1]
            far = np.hypot(xx - (dx - R), yy - (dy - R)) * step >= 10
            peaks = far & (win == ndimage.maximum_filter(win, size=5))   # other peaks, not the slope of the best one
            row.append(round(float(win[peaks].max() / max(win.max(), 1e-9)) if peaks.any() else 0.0, 3))
        out.append(row)
        if say and n % 50 == 0:
            say(n / len(todo))
    return out


def _robust_affine(P, D, w=None):
    """offset = (P - mean) @ B + b fitted by least squares, refitted twice without the windows far off the fit."""
    keep = np.ones(len(P), bool)
    mu = P.mean(0)
    for _ in range(3):
        X = np.c_[P[keep] - mu, np.ones(keep.sum())]
        coef = np.linalg.lstsq(X, D[keep], rcond=None)[0]
        r = np.hypot(*(D - np.c_[P - mu, np.ones(len(P))] @ coef).T)
        cut = max(3.0, 3 * 1.4826 * np.median(np.abs(r[keep] - np.median(r[keep]))) + np.median(r[keep]))
        keep = r <= cut
        if keep.sum() < 8:
            break
    return mu, coef


def _analyse(w):
    """The windows as arrays, the affine part, the local shifts and which windows are noise."""
    a = np.asarray(w["windows"], float).reshape(-1, 9)
    edge = np.abs(a[:, 4:6]).max(1) >= w["search_um"]   # no peak within the search: repeating nuclei, not a shift
    ok = (a[:, 6] >= MIN_CORR) & ~edge
    P = a[:, 2:4] * w["um_per_px"] / 1000      # mm
    D = np.where((a[:, 8] <= GAIN * a[:, 6])[:, None], a[:, 4:6], 0.0)  # um; a shift that fits barely better is none
    if ok.sum() >= 8:
        mu, coef = _robust_affine(P[ok], D[ok])
    else:
        mu, coef = P.mean(0), np.zeros((3, 2))
    aff = np.c_[P - mu, np.ones(len(P))] @ coef
    if np.hypot(*aff.T).max() < OFF_UM:  # nothing a registration would fix: all of it is local
        coef, aff = np.zeros((3, 2)), np.zeros_like(aff)
    res = D - aff
    # noise: a weak window disagreeing with the median of its trusted neighbors
    pos = {(int(i), int(j)): n for n, (i, j) in enumerate(a[:, :2])}
    noise = ~ok
    for n, (i, j) in enumerate(a[:, :2].astype(int)):
        if not ok[n]:
            continue
        nb = [pos[(i + di, j + dj)] for di in (-1, 0, 1) for dj in (-1, 0, 1)
              if (di or dj) and (i + di, j + dj) in pos and ok[pos[(i + di, j + dj)]]]
        if nb and np.hypot(*(res[n] - np.median(res[nb], 0))) > OFF_UM and a[n, 6] < 0.3:
            noise[n] = True
    return a, mu, coef, aff, res, noise, D


def _areas(a, D, noise):
    """Connected groups (8-neighbors) of trusted windows whose cells are off their nuclei."""
    off = (np.hypot(*D.T) >= OFF_UM) & ~noise
    grid = {}
    for n in np.where(off)[0]:
        grid[(int(a[n, 0]), int(a[n, 1]))] = n
    lab, seen = [], set()
    for start in grid:
        if start in seen:
            continue
        stack, comp = [start], []
        seen.add(start)
        while stack:
            i, j = stack.pop()
            comp.append(grid[(i, j)])
            for di in (-1, 0, 1):
                for dj in (-1, 0, 1):
                    q = (i + di, j + dj)
                    if q in grid and q not in seen:
                        seen.add(q)
                        stack.append(q)
        lab.append(sorted(comp))
    return sorted(lab, key=lambda c: -len(c))


def _focus_masks(s, focus, version=None):
    from .insights import _lineage, _members
    _, types, lin = _lineage(s, version)
    if focus:
        return {name: _members(s, list(v) if isinstance(v, (list, tuple)) else [v], types, lin) for name, v in focus.items()}
    names = sorted(set(s.types) | set(s.lineage_of_type.values()))
    out = {}
    for name, pat in DEFAULT_FOCUS.items():
        hit = [n for n in names if re.search(pat, str(n), re.I)]
        if hit:
            out[name] = _members(s, hit, types, lin)
    return out


def _window_of_cells(s, w):
    tile = w["tile_px"]
    x0, y0 = w["origin_px"]
    return (((s.cells.py.to_numpy() - y0) // tile).astype(int) * w["nx"]
            + ((s.cells.px.to_numpy() - x0) // tile).astype(int))


def summary(s, focus=None, version=None, progress=None):
    """What windows() found, for the agent and the review: the affine part, the areas off their nuclei (with the
    structures and regions in them), offsets by structure and by region, and what the correction would do."""
    try:
        w = windows(s.sample_id, progress)
    except ValueError as e:  # no image, or Visium bins
        return {"error": str(e)}
    if not w["windows"]:
        return {"error": "no window with enough cells on the image"}
    a, mu, coef, aff, res, noise, D = _analyse(w)
    tot, raw = np.hypot(*D.T), np.hypot(a[:, 4], a[:, 5])   # offsets that count; as measured
    wid = (a[:, 0] * w["nx"] + a[:, 1]).astype(int)
    cw = _window_of_cells(s, w)
    masks = _focus_masks(s, focus, version)
    reg = np.asarray(s.result(version)["region"], dtype=object)
    share = {}
    for name, mk in masks.items():
        num = pd.Series(mk.astype(float)).groupby(cw).mean()
        share[name] = num.reindex(wid).fillna(0).to_numpy()
    rs = pd.DataFrame({"w": cw, "r": reg}).groupby("w").r.agg(lambda x: x.value_counts().index[0])
    top_region = rs.reindex(wid).fillna("unknown").to_numpy()
    tag = np.full(len(a), None, dtype=object)
    if share:
        names = list(share)
        S = np.c_[[share[n] for n in names]].T
        best = S.argmax(1)
        tag = np.where(S.max(1) >= 0.3, np.array(names, dtype=object)[best], None)

    def stats(sel):
        sel = sel & ~noise
        if not sel.any():
            return None
        return {"windows": int(sel.sum()), "median_offset_um": round(float(np.median(raw[sel])), 1),
                "max_offset_um": round(float(tot[sel].max()), 1),
                "share_off": round(float((tot[sel] >= OFF_UM).mean()), 3)}
    b, big, _ = _large(w)
    fix = correction(s.sample_id, w)
    fw = {tuple(x) for x in fix["false_windows"]}   # offsets that stayed with the cells moved: not a shift
    false = np.array([(int(i), int(j)) in fw for i, j in a[:, :2]], bool)
    in_big = np.zeros(len(a), bool)   # fine windows inside a large shift: reported with it
    if big.any():
        c = w["coarse"]
        cid = ((a[:, 3] - c["origin_px"][1]) // c["tile_px"]).astype(int) * c["nx"] + ((a[:, 2] - c["origin_px"][0]) // c["tile_px"]).astype(int)
        in_big = np.isin(cid, (b[big, 0] * c["nx"] + b[big, 1]).astype(int))
    areas = []
    for k, comp in enumerate(_areas(a, D, noise | in_big | false)):
        c = np.array(comp)
        tags = pd.Series([t for t in tag[c] if t]).value_counts()
        regs = pd.Series(top_region[c]).value_counts()
        areas.append({"id": f"L{k + 1}", "windows": len(c), "box_px": [round(float(f(a[c, k2]))) for f, k2 in ((np.min, 2), (np.min, 3), (np.max, 2), (np.max, 3))],
                      "median_offset_um": round(float(np.median(tot[c])), 1), "max_offset_um": round(float(tot[c].max()), 1),
                      "direction_um": [round(float(x), 1) for x in np.median(D[c], 0)],
                      "local_after_affine_um": round(float(np.median(np.hypot(*res[c].T))), 1),
                      "structures": {k2: int(v) for k2, v in tags.items()}, "regions": {str(k2): int(v) for k2, v in regs.items()},
                      "correlation": round(float(np.median(a[c, 6])), 2),
                      "center_px": [round(float(a[c, 2].mean())), round(float(a[c, 3].mean()))]})
    if big.any():  # the large shifts: groups of neighboring coarse windows
        Db = np.where(big[:, None], b[:, 4:6], 0.0)
        for comp in _areas(b, Db, ~big):
            c = np.array(comp)
            off = np.hypot(*Db[c].T)
            mine = np.isin(cid, (b[c, 0] * w["coarse"]["nx"] + b[c, 1]).astype(int))
            tags = pd.Series([t for t in tag[mine] if t]).value_counts()
            regs = pd.Series(top_region[mine]).value_counts()
            areas.append({"id": "", "windows": len(c), "window_um": w["coarse"]["window_um"], "half_px": w["coarse"]["tile_px"] / 2,
                          "box_px": [round(float(f(b[c, k2]))) for f, k2 in ((np.min, 2), (np.min, 3), (np.max, 2), (np.max, 3))],
                          "median_offset_um": round(float(np.median(off)), 1), "max_offset_um": round(float(off.max()), 1),
                          "direction_um": [round(float(x), 1) for x in np.median(Db[c], 0)],
                          "local_after_affine_um": round(float(np.median(off)), 1),
                          "structures": {k2: int(v) for k2, v in tags.items()}, "regions": {str(k2): int(v) for k2, v in regs.items()},
                          "correlation": round(float(np.median(b[c, 6])), 2), "large_shift": True,
                          "center_px": [round(float(b[c, 2].mean())), round(float(b[c, 3].mean()))]})
    areas.sort(key=lambda x: -x["windows"] * (x.get("window_um", w["window_um"]) / w["window_um"]) ** 2)
    for k, x in enumerate(areas):
        x["id"] = f"L{k + 1}"
    A = coef[:2].T / 1000
    affine = {"scale_x": round(1 + float(A[0, 0]), 5), "scale_y": round(1 + float(A[1, 1]), 5),
              "rotation_deg": round(float(np.degrees((A[1, 0] - A[0, 1]) / 2)), 3),
              "max_offset_um": round(float(np.hypot(*aff.T).max()), 1)}
    after, ta = np.hypot(*fix["after_um"].T), fix["trusted_after"]
    trusted = ~noise
    n_off = int(((tot >= OFF_UM) & trusted & ~false).sum()) + int(big.sum())
    loc = [x for x in areas if x["local_after_affine_um"] >= OFF_UM]
    verdict = ("fine" if not n_off else
               "global" if affine["max_offset_um"] >= OFF_UM and not loc else
               "local" if affine["max_offset_um"] < OFF_UM else "global and local")
    return {"verdict": verdict, "windows": int(len(a)), "trusted_windows": int(trusted.sum()), "window_um": w["window_um"],
            "windows_off": n_off, "median_offset_um": round(float(np.median(raw[trusted])), 1) if trusted.any() else None,
            "affine_part": affine, "areas": areas[:15], "n_areas": len(areas),
            "noise_windows": int((noise & (a[:, 6] >= MIN_CORR)).sum()),
            "beyond_search_share": round(float((np.abs(a[:, 4:6]).max(1) >= w["search_um"]).mean()), 3),
            "by_structure": {k: stats(tag == k) for k in masks},
            "by_region": {str(r): stats(top_region == r) for r in pd.unique(top_region)},
            "correction": {"moves_cells_um_median": fix["moved_um_median"], "moves_cells_um_max": fix["moved_um_max"],
                           "checked_windows": fix["checked_windows"],
                           "false_readings": {"windows": int(false.sum()),
                                              "centers_px": [[round(float(x)), round(float(y))] for x, y in a[false, 2:4]],
                                              "offsets_um": [[float(x), float(y)] for x, y in a[false, 4:6]]},
                           "windows_off_after": int(((after >= OFF_UM) & ta).sum()),
                           "median_offset_after_um": round(float(np.median(np.hypot(*fix["after_raw_um"][ta].T))), 1)
                           if ta.any() else None},
            "note": "offset = where the nuclei sit relative to the cell centroids (+dx right, +dy down, um); a window's mean "
                    f"offset is about 1 um when the cells sit on their nuclei; areas: neighboring windows off by >= {OFF_UM:g} um; "
                    f"large_shift areas: {COARSE_WINDOW_UM:g} um windows off by {COARSE_MIN_UM:g}-{COARSE_SEARCH_UM} um, beyond the "
                    f"{SEARCH_UM} um the fine windows look; "
                    "beyond_search_share above ~0.3: the cells are far off everywhere, redo the registration instead; correction: "
                    "measured again with the cells moved by it (windows_off_after, median_offset_after_um); "
                    "false_readings: windows whose offset stayed (as large, the same way) when their cells were moved onto the nuclei "
                    "(a repeating structure matched one period off), left out of the areas and left alone"}


def _spread(C, r, sig):
    """Shifts at the points C whose Gaussian field (_field) comes close to r there (a Gaussian average alone smooths
    them away)."""
    S = np.zeros_like(r, dtype=float)
    for _ in range(4):
        S += r - _field(C, C, S, sig)
    return S


def _large(w):
    """The coarse windows (rows), and which of them are shifted further than the fine pass reaches (a unique peak,
    clearly better than none, where most fine windows found no offset they trust) or are anchors (their cells stay
    for the fine pass: its windows there trusted, or a trusted short shift)."""
    c = w.get("coarse") or {}
    b = np.asarray(c.get("windows") or [], float).reshape(-1, 10)
    if not len(b):
        return b, np.zeros(0, bool), np.zeros(0, bool)
    a, *_, noise, D = _analyse(w)
    (x0, y0), ct = c["origin_px"], c["tile_px"]
    cid = ((a[:, 3] - y0) // ct).astype(int) * c["nx"] + ((a[:, 2] - x0) // ct).astype(int)
    share = pd.Series((~noise).astype(float)).groupby(cid).mean()
    fine_ok = share.reindex((b[:, 0] * c["nx"] + b[:, 1]).astype(int)).fillna(0).to_numpy() >= 0.5
    ok = (b[:, 6] >= MIN_CORR) & (np.abs(b[:, 4:6]).max(1) < c["search_um"]) & (b[:, 9] <= UNIQUE)
    big = ok & (np.hypot(b[:, 4], b[:, 5]) >= COARSE_MIN_UM) & (b[:, 8] <= GAIN * b[:, 6]) & ~fine_ok
    anchor = ~big & (ok | fine_ok)
    return b, big, anchor


def _pre(w):
    """The field of the large shifts (image pixels), or None when there are none."""
    b, big, anchor = _large(w)
    if not big.any():
        return None
    use = big | anchor
    C, sig = b[use, 2:4], FIELD_SIGMA * w["coarse"]["tile_px"]
    r = np.where(big[use, None], b[use, 4:6], 0.0) / w["um_per_px"]
    return {"points_px": C.round(1).tolist(), "shift_px": _spread(C, r, sig).round(3).tolist(), "sigma_px": sig}


_moved = {}  # (placement key, pre) -> the fine windows with the cells moved by pre, measured again


def _rescan(w, p, q):
    """The fine windows of w measured again with the cells at q instead of p, where they moved by more than 1 um:
    the windows they left or entered (rows keyed (i, j); a window left with too few cells has no row)."""
    (x0, y0), tile, nx = w["origin_px"], w["tile_px"], w["nx"]
    wid = lambda z: ((z[:, 1] - y0) // tile).astype(int) * nx + ((z[:, 0] - x0) // tile).astype(int)
    hit = np.hypot(*(q - p).T) * w["um_per_px"] > 1
    if not hit.any():
        return {}, set()
    keys = set(wid(p[hit]).tolist()) | set(wid(q[hit]).tolist())
    new = _scan(render.slide(w["image"]), w["kind"], q[:, 0], q[:, 1], w, keys=keys)
    return {(int(r[0]), int(r[1])): r for r in new}, keys


def _remeasure(sample_id, w, pre):
    """w with its fine windows measured again where the large shifts move the cells (by more than 1 um)."""
    k = (w["key"], json.dumps(pre))
    if k not in _moved:
        cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
        p = cells[["px", "py"]].to_numpy()
        q = p + _field(p, np.asarray(pre["points_px"]), np.asarray(pre["shift_px"]), pre["sigma_px"])
        new, keys = _rescan(w, p, q)
        nx = w["nx"]
        rows = {(r[0], r[1]): r for r in w["windows"] if r[0] * nx + r[1] not in keys}
        rows.update(new)
        _moved.clear()
        _moved[k] = dict(w, windows=list(rows.values()))
    return _moved[k]


_checked = {}  # (placement key, pre) -> the correction, checked by measuring again


def correction(sample_id, w=None):
    """The correction of the current placement: "pre": the large shifts (or none), then {"origin", "A", "b"} (affine
    part, image pixels, as importer.apply_registration_fix takes it) and "local": the local shifts at the trusted
    windows (image pixels), spread with a Gaussian of one window and fading out away from them. Also how far it moves
    the cells.

    Checked before it is offered: the cells are moved by it and the windows they move in measured again. A real shift
    shrinks; an offset that stays (as large, the same way) was a false reading, mostly a repeating structure (villi,
    crypts) matched one period off: that window keeps its cells (a zero shift in the field) and the field is made
    again without it (CHECK_ROUNDS times at most). A reading that jumps to another peak proves nothing either way:
    there the prediction stands. "false_windows": those windows [i, j]; "checked_windows": how many were checked;
    "after_um" / "after_raw_um": the offsets with the cells moved (that count / as measured; measured where the
    measure is clear, else predicted; windows the correction does not move keep theirs); "trusted_after": the windows those count at (not noise, not false)."""
    w = w or windows(sample_id)
    pre = _pre(w)
    k = (w["key"], json.dumps(pre))
    if k in _checked:
        return _checked[k]
    w2 = _remeasure(sample_id, w, pre) if pre else w
    a, mu, coef, aff, res, noise, D = _analyse(w2)
    upp = w["um_per_px"]
    r = res.copy()
    n = np.hypot(*r.T)
    r *= (np.maximum(n - FIELD_FLOOR_UM, 0) / np.maximum(n, 1e-9))[:, None]   # leave small shifts alone
    use = ~noise
    C, sig = a[use, 2:4], FIELD_SIGMA * w["tile_px"]
    origin = (mu * 1000 / upp).tolist()                      # px
    A = coef[:2].T / 1000                                    # um offset per um position
    cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
    p = cells[["px", "py"]].to_numpy()
    ij = [(int(i), int(j)) for i, j in a[:, :2]]
    false, checked = np.zeros(len(a), bool), np.zeros(len(a), bool)
    for rnd in range(CHECK_ROUNDS + 1):
        rr = np.where(false[:, None], 0.0, r)                # a false reading: its cells stay where they are
        S = _spread(C, rr[use] / upp, sig)
        fix = {"origin": origin, "A": A.tolist(), "b": (coef[2] / upp).tolist(),
               "local": {"points_px": C.round(1).tolist(), "shift_px": S.round(3).tolist(), "sigma_px": sig}}
        at = displacement(a[:, 2:4], fix) * upp              # at the windows as measured (after the large shifts)
        if pre:
            fix = {"pre": pre} | fix
        mv = displacement(p, fix)
        new, _ = _rescan(w, p, p + mv)
        pred_raw, pred = a[:, 4:6] - at, D - at              # as the field predicts
        meas_raw, meas, got = pred_raw.copy(), pred.copy(), np.zeros(len(a), bool)
        for m, key in enumerate(ij):
            row = new.get(key)
            if row is not None and row[6] >= MIN_CORR:
                got[m] = True
                meas_raw[m] = row[4:6]
                meas[m] = (0.0, 0.0) if row[8] > GAIN * row[6] else row[4:6]
        nm, nd = np.hypot(*meas.T), np.hypot(*D.T)
        cos = (meas * D).sum(1) / np.maximum(nm * nd, 1e-9)
        stayed = got & (nm >= OFF_UM) & (nm >= STAYED_SHARE * nd) & (cos >= STAYED_COS)
        sure = got & ((nm < OFF_UM) | stayed)                # on its nuclei now, or the offset stayed: else it jumped
        after, after_raw = np.where(sure[:, None], meas, pred), np.where(sure[:, None], meas_raw, pred_raw)
        size = np.hypot(*rr.T)
        judge = use & ~false & got & (size > 0) & (np.hypot(*at.T) >= 0.5 * size)   # moved most of the way
        checked |= judge
        bad = judge & stayed
        if not bad.any() or rnd == CHECK_ROUNDS:
            break
        false |= bad
    out = fix | {"moved_um_median": round(float(np.median(np.hypot(*mv.T)) * upp), 1),
                 "moved_um_max": round(float(np.hypot(*mv.T).max() * upp), 1),
                 "after_um": after, "after_raw_um": after_raw, "trusted_after": ~noise & ~false,
                 "false_windows": [list(ij[m]) for m in np.where(false)[0]], "checked_windows": int(checked.sum())}
    _checked.clear()
    _checked[k] = out
    return out


FIX_KEYS = ("pre", "origin", "A", "b", "local")  # what a correction keeps (and importer.apply_registration_fix takes)


def displacement(p, fix, chunk=20000):
    """How far the correction moves points p (image pixels): the large shifts (pre), then, where they put the
    points, A (p - origin) + b + the local field."""
    p = np.asarray(p, float)
    pre = fix.get("pre")
    d0 = _field(p, np.asarray(pre["points_px"]), np.asarray(pre["shift_px"]), float(pre["sigma_px"]), chunk) if pre \
        else np.zeros_like(p)
    q = p + d0
    out = d0 + (q - np.asarray(fix["origin"])) @ np.asarray(fix["A"]).T + np.asarray(fix["b"])
    loc = fix.get("local")
    if loc and loc["points_px"]:
        out += _field(q, np.asarray(loc["points_px"]), np.asarray(loc["shift_px"]), float(loc["sigma_px"]), chunk)
    return out


def _field(p, C, S, sig, chunk=20000):
    """Gaussian-weighted average at p of the shifts S at the points C, fading out where no point is near."""
    out = np.zeros((len(p), 2))
    for k in range(0, len(p), chunk):
        wt = np.exp(-((p[k:k + chunk, None, :] - C[None]) ** 2).sum(-1) / (2 * sig ** 2))
        out[k:k + chunk] = (wt @ S) / (wt.sum(1) + 0.1)[:, None]
    return out


def figure(s, info, version=None):
    """The section with the offset per window (colored where >= 3 um) and the areas off their nuclei (left; large
    shifts dashed); close-ups of up to four areas with the cells now (cyan) and after the correction (yellow) on the
    nuclei."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    w = windows(s.sample_id)
    a = np.asarray(w["windows"], float).reshape(-1, 9)
    sl = render.slide(w["image"])
    upp, tile = w["um_per_px"], w["tile_px"]
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    bx0, by0, bx1, by1 = px.min() - tile, py.min() - tile, px.max() + tile, py.max() + tile
    sc = 1000 / max(bx1 - bx0, by1 - by0)
    tw, th = max(int((bx1 - bx0) * sc), 1), max(int((by1 - by0) * sc), 1)
    areas = info.get("areas", [])[:4]
    fig = plt.figure(figsize=(15, 7))
    ax = fig.add_axes([0.01, 0.03, 0.47, 0.88])
    ax.imshow(np.asarray(sl.read(bx0, by0, bx1, by1, tw, th)), extent=(bx0, bx1, by1, by0), alpha=0.6)
    off = np.hypot(a[:, 4], a[:, 5])
    good = a[:, 6] >= MIN_CORR
    x0, y0 = w["origin_px"]
    ny = int(a[:, 0].max()) + 1
    G = np.full((ny, w["nx"]), np.nan)
    G[a[good, 0].astype(int), a[good, 1].astype(int)] = off[good]
    G[G < 3] = np.nan  # windows on their nuclei stay clear
    im = ax.imshow(G, extent=(x0, x0 + w["nx"] * tile, y0 + ny * tile, y0), cmap="plasma", vmin=3, vmax=15,
                   interpolation="nearest", alpha=0.85)
    cb = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.01)
    cb.set_label("offset (um), windows >= 3 um", fontsize=8)
    for x in info.get("areas", []):
        x0, y0, x1, y1 = x["box_px"]  # window centers
        h = x.get("half_px", tile / 2)
        ax.add_patch(Rectangle((x0 - h, y0 - h), x1 - x0 + 2 * h, y1 - y0 + 2 * h, fill=False, ec="red", lw=1.2,
                               ls="--" if x.get("large_shift") else "-"))
        ax.text(x0 - h, y0 - h, x["id"], color="red", fontsize=8, va="bottom")
    fr = ((info.get("correction") or {}).get("false_readings") or {}).get("centers_px") or []
    for fx_, fy_ in fr:  # offsets that stayed with the cells moved: not a shift, left alone
        ax.add_patch(Rectangle((fx_ - tile / 2, fy_ - tile / 2), tile, tile, fill=False, ec="0.35", lw=1.0, ls=":"))
    ax.set_xlim(bx0, bx1); ax.set_ylim(by1, by0); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"offset of the nuclei from the cells per {w['window_um']:g} um window; "
                 f"{info['windows_off']} of {info['trusted_windows']} windows off; verdict: {info['verdict']}"
                 + (f"\n{len(fr)} false reading(s), dotted: the offset stayed with the cells moved" if fr else ""), fontsize=9)
    if areas:
        fix = correction(s.sample_id, w)
        for k, x in enumerate(areas):
            cx, cy = x["center_px"]
            if not x.get("large_shift"):
                d2 = (a[:, 2] - cx) ** 2 + (a[:, 3] - cy) ** 2
                cx, cy = a[np.argmin(d2 + 1e12 * (off < OFF_UM)), 2:4]  # the area's window nearest its center
            inw = (np.abs(px - cx) < tile / 2) & (np.abs(py - cy) < tile / 2)
            if inw.any():
                cx, cy = np.median(px[inw]), np.median(py[inw])     # on its cells
            half = max(50.0, 1.5 * x["median_offset_um"]) / upp
            m = (np.abs(px - cx) < half) & (np.abs(py - cy) < half)
            moved = np.c_[px[m], py[m]] + displacement(np.c_[px[m], py[m]], fix)
            a2 = fig.add_axes([0.53 + 0.235 * (k % 2), 0.5 - 0.47 * (k // 2), 0.22, 0.42])
            a2.imshow(np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, 400, 400)),
                      extent=(cx - half, cx + half, cy + half, cy - half))
            a2.scatter(px[m], py[m], s=10, facecolors="none", edgecolors="#00e6ff", linewidths=0.9)
            a2.scatter(moved[:, 0], moved[:, 1], s=5, c="#ffd400", linewidths=0)
            a2.set_xlim(cx - half, cx + half); a2.set_ylim(cy + half, cy - half); a2.set_xticks([]); a2.set_yticks([])
            st = ", ".join(x["structures"]) or "no focus structure"
            a2.set_title(f"{x['id']} ({st}): {x['median_offset_um']} um off; {2 * half * upp:.0f} um\ncyan now, yellow corrected", fontsize=7)
    out = s.dir / "figures" / f"registration_map_v{s.state(version)['version']}.jpg"
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=90, pil_kwargs={"quality": 85}); plt.close(fig)
    return str(out)
