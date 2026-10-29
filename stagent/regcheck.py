"""The registration checked over the whole tissue, and its local correction (Agent Review: Spatial organization).

windows(): the cells are cut into square windows (WINDOW_UM) over all the tissue; in each, the nuclei on the image
(hematoxylin of the H&E, or DAPI) are cross-correlated with the cell centroids, giving the offset of the nuclei from
the cells (um; insights.registration does this for a sample of windows). A robust affine fit of the offsets is the
part one registration fixes (scale, rotation, shift); what is left are local shifts (tissue stretched, folded or torn
between the run and the staining). Neighboring windows off their nuclei form areas (L1, L2, ...). Not trusted: a
window whose best offset lies on the edge of the search (densely packed nuclei repeat, so a far shift can fit as
well), and a weak one that disagrees with its neighbors; an offset that fits barely better than none counts as none.
summary(): the areas, and the offsets by structure (focus groups of types: crypts are where the epithelium is,
vessels where the endothelium is; their nuclei have a clear pattern, so a shift shows first there) and by region.
correction(): affine + a smooth field of the local shifts (Gaussian-weighted, fading out away from the windows), as
importer.apply_registration_fix applies it. Nothing here moves a cell.

The offsets of a placement are kept in the sample's folder (regcheck.json) until the cells move or the image changes.
"""
import json, os, re, time
import numpy as np, pandas as pd
from scipy import ndimage
from . import config, render

FORMAT = 2
WINDOW_UM = 200.0   # side of a window
SEARCH_UM = 25      # largest offset looked for
MIN_CELLS = 40      # windows with fewer cells are not measured
MIN_CORR = 0.2      # below: the offset is not trusted
GAIN = 0.85         # an offset counts when the correlation at zero is at most this share of the one at the offset
OFF_UM = 5.0        # a window's cells are off their nuclei above this (a window's mean offset is ~1 um when fine)
FIELD_FLOOR_UM = 2.0  # local shifts below this are left alone by the correction
FIELD_SIGMA = 0.5     # the local field's Gaussian, in windows
DEFAULT_FOCUS = {"crypts / epithelium": r"epith|crypt|colonocyte|enterocyte|goblet|stem|transit|paneth",
                 "vessels": r"endoth|vascul|vessel|pericyte",
                 "smooth muscle": r"smooth|muscle|myocyte"}


def xcorr(sl, kind, cx, cy, half, W, R, qx, qy):
    """Cross-correlation (2R+1 square, 1 um steps) of the nuclei in the image window around (cx, cy) (half: its half
    side in image pixels, read at W x W pixels) with the cell centroids qx, qy (image pixels) in it; the peak at
    (R + dy, R + dx) means the nuclei sit dx, dy um from the cells."""
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
    busy.write_text(str(os.getpid()))
    try:
        cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
        upp = float(m.get("he_um_per_px") or config.HE_UM_PER_PX)
        sl = render.slide(image)
        tile = WINDOW_UM / upp
        px, py = cells.px.to_numpy(), cells.py.to_numpy()
        x0, y0 = float(px.min()), float(py.min())
        tj, ti = ((px - x0) // tile).astype(int), ((py - y0) // tile).astype(int)
        nx = int(tj.max()) + 1
        groups = pd.Series(np.arange(len(px))).groupby(ti * nx + tj).indices
        W, R, half = int(WINDOW_UM), SEARCH_UM, tile / 2
        out, todo = [], [(k, rows) for k, rows in groups.items() if len(rows) >= MIN_CELLS]
        for n, (k, rows) in enumerate(todo):
            i, j = divmod(int(k), nx)
            cx, cy = x0 + (j + 0.5) * tile, y0 + (i + 0.5) * tile
            if cx - half < 0 or cy - half < 0 or cx + half > sl.W0 or cy + half > sl.H0:
                continue
            win = xcorr(sl, kind, cx, cy, half, W, R, px[rows], py[rows])
            dy, dx = np.unravel_index(np.argmax(win), win.shape)
            out.append([i, j, round(cx, 1), round(cy, 1), int(dx - R), int(dy - R), round(float(win.max()), 3), len(rows),
                        round(float(win[R, R]), 3)])
            if n % 50 == 0:
                say(n / len(todo), f"registration: window {n} of {len(todo)}")
        res = {"key": key, "image": image, "kind": kind, "um_per_px": upp, "window_um": WINDOW_UM, "search_um": R,
               "origin_px": [x0, y0], "tile_px": tile, "nx": nx, "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "columns": ["i", "j", "cx", "cy", "dx_um", "dy_um", "corr", "n_cells", "corr_at_zero"], "windows": out}
        json.dump(res, open(f, "w"))
        return res
    finally:
        busy.unlink(missing_ok=True)


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
    w = windows(s.sample_id, progress)
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
    areas = []
    for k, comp in enumerate(_areas(a, D, noise)):
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
    A = coef[:2].T / 1000
    affine = {"scale_x": round(1 + float(A[0, 0]), 5), "scale_y": round(1 + float(A[1, 1]), 5),
              "rotation_deg": round(float(np.degrees((A[1, 0] - A[0, 1]) / 2)), 3),
              "max_offset_um": round(float(np.hypot(*aff.T).max()), 1)}
    fix = correction(s.sample_id, w)
    after = np.hypot(*(D - fix["at_windows_um"]).T)
    trusted = ~noise
    n_off = int(((tot >= OFF_UM) & trusted).sum())
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
                           "windows_off_after": int(((after >= OFF_UM) & trusted).sum()),
                           "median_offset_after_um": round(float(np.median(np.hypot(*(a[trusted, 4:6] - fix["at_windows_um"][trusted]).T))), 1)
                           if trusted.any() else None},
            "note": "offset = where the nuclei sit relative to the cell centroids (+dx right, +dy down, um); a window's mean "
                    f"offset is about 1 um when the cells sit on their nuclei; areas: neighboring windows off by >= {OFF_UM:g} um; "
                    "beyond_search_share above ~0.3: the cells are far off everywhere, redo the registration instead"}


def correction(sample_id, w=None):
    """The correction of the current placement: {"origin", "A", "b"} (affine part, image pixels, as
    importer.apply_registration_fix takes it) and "local": the local shifts at the trusted windows (image pixels),
    spread with a Gaussian of one window and fading out away from them."""
    w = w or windows(sample_id)
    a, mu, coef, aff, res, noise, _ = _analyse(w)
    upp = w["um_per_px"]
    r = res.copy()
    n = np.hypot(*r.T)
    r *= (np.maximum(n - FIELD_FLOOR_UM, 0) / np.maximum(n, 1e-9))[:, None]   # leave small shifts alone
    use = ~noise
    C, sig = a[use, 2:4], FIELD_SIGMA * w["tile_px"]
    S = np.zeros((use.sum(), 2))
    for _ in range(4):  # the field at the windows comes close to their shifts (a Gaussian average alone smooths them away)
        S += r[use] / upp - _field(C, C, S, sig)
    origin = (mu * 1000 / upp).tolist()                      # px
    A = coef[:2].T / 1000                                    # um offset per um position
    fix = {"origin": origin, "A": A.tolist(), "b": (coef[2] / upp).tolist(),
           "local": {"points_px": C.round(1).tolist(), "shift_px": S.round(3).tolist(), "sigma_px": sig}}
    cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["px", "py"])
    p = cells[["px", "py"]].to_numpy()
    mv = displacement(p, fix)
    at = displacement(a[:, 2:4], fix) * upp
    return fix | {"moved_um_median": round(float(np.median(np.hypot(*mv.T)) * upp), 1),
                  "moved_um_max": round(float(np.hypot(*mv.T).max() * upp), 1), "at_windows_um": at}


def displacement(p, fix, chunk=20000):
    """How far the correction moves points p (image pixels): A (p - origin) + b + the local field."""
    p = np.asarray(p, float)
    out = (p - np.asarray(fix["origin"])) @ np.asarray(fix["A"]).T + np.asarray(fix["b"])
    loc = fix.get("local")
    if loc and loc["points_px"]:
        out += _field(p, np.asarray(loc["points_px"]), np.asarray(loc["shift_px"]), float(loc["sigma_px"]), chunk)
    return out


def _field(p, C, S, sig, chunk=20000):
    """Gaussian-weighted average at p of the shifts S at the points C, fading out where no point is near."""
    out = np.zeros((len(p), 2))
    for k in range(0, len(p), chunk):
        wt = np.exp(-((p[k:k + chunk, None, :] - C[None]) ** 2).sum(-1) / (2 * sig ** 2))
        out[k:k + chunk] = (wt @ S) / (wt.sum(1) + 0.1)[:, None]
    return out


def figure(s, info, version=None):
    """The section with the offset per window (colored where >= 3 um) and the areas off their nuclei (left); close-ups of up to four areas with the cells now (cyan) and after the correction (yellow) on the nuclei."""
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
        ax.add_patch(Rectangle((x0 - tile / 2, y0 - tile / 2), x1 - x0 + tile, y1 - y0 + tile, fill=False, ec="red", lw=1.2))
        ax.text(x0 - tile / 2, y0 - tile / 2, x["id"], color="red", fontsize=8, va="bottom")
    ax.set_xlim(bx0, bx1); ax.set_ylim(by1, by0); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"offset of the nuclei from the cells per {w['window_um']:g} um window; "
                 f"{info['windows_off']} of {info['trusted_windows']} windows off; verdict: {info['verdict']}", fontsize=9)
    if areas:
        fix = correction(s.sample_id, w)
        for k, x in enumerate(areas):
            cx, cy = x["center_px"]
            d2 = (a[:, 2] - cx) ** 2 + (a[:, 3] - cy) ** 2
            cx, cy = a[np.argmin(d2 + 1e12 * (off < OFF_UM)), 2:4]  # the area's window nearest its center
            inw = (np.abs(px - cx) < tile / 2) & (np.abs(py - cy) < tile / 2)
            cx, cy = np.median(px[inw]), np.median(py[inw])         # on its cells
            half = 50 / upp
            m = (np.abs(px - cx) < half) & (np.abs(py - cy) < half)
            moved = np.c_[px[m], py[m]] + displacement(np.c_[px[m], py[m]], fix)
            a2 = fig.add_axes([0.53 + 0.235 * (k % 2), 0.5 - 0.47 * (k // 2), 0.22, 0.42])
            a2.imshow(np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, 400, 400)),
                      extent=(cx - half, cx + half, cy + half, cy - half))
            a2.scatter(px[m], py[m], s=10, facecolors="none", edgecolors="#00e6ff", linewidths=0.9)
            a2.scatter(moved[:, 0], moved[:, 1], s=5, c="#ffd400", linewidths=0)
            a2.set_xlim(cx - half, cx + half); a2.set_ylim(cy + half, cy - half); a2.set_xticks([]); a2.set_yticks([])
            st = ", ".join(x["structures"]) or "no focus structure"
            a2.set_title(f"{x['id']} ({st}): {x['median_offset_um']} um off; 100 um\ncyan now, yellow corrected", fontsize=7)
    out = s.dir / "figures" / f"registration_map_v{s.state(version)['version']}.jpg"
    out.parent.mkdir(exist_ok=True)
    plt.savefig(out, dpi=90, pil_kwargs={"quality": 85}); plt.close(fig)
    return str(out)
