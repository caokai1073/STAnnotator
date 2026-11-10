"""Aligning a sample's cells to its image without a registration file (the app's "Align with Claude").

1. search(): the tissue on a thumbnail of the image is compared with the footprint of the cells (where cells are,
   from their um coordinates) for 8 orientations (4 rotations, mirrored or not) and a range of image pixel sizes
   (around the size stored in an OME-TIFF, else 0.1-2 um/px), by FFT cross-correlation. The score is the agreement
   of tissue and footprint over the footprint's box (+1 both or neither, -1 one only), so a footprint that is too
   small or too large scores low. The best few are kept as candidates (affines cell um -> image level-0 px).
2. Claude looks at the candidates (show()), picks one (several sections on one slide, symmetric tissue).
3. refine(): moves the candidate to micrometer precision: nuclei on the image vs cell centroids in windows over the
   tissue (as insights.registration), an affine fit of the window offsets, three rounds with a shrinking search.
4. deviation(): the refined candidate checked over the whole tissue, as the registration check does it (regcheck: the
   nuclei vs the cells in 200 um windows everywhere, and 600 um windows for shifts up to 100 um): how far the cells
   sit from their nuclei, and whether areas sit off them (tissue stretched, folded or torn between the run and the
   staining), which one alignment cannot fix: then the local registration (Review: Spatial organization) is advised.
   The review window of apply_alignment shows its histogram.
5. apply(): puts the cells on the image (importer.place_cells; undo with importer.undo_placement); the deviation's
   windows become the registration check of the new placement (regcheck.json), so the local registration starts
   from them.
"""
import json, time
import numpy as np
from scipy import ndimage
from . import config, importer, render, insights, regcheck

LOCAL_SHARE = 0.05    # the local registration is advised when this share of the tissue's windows sits off its nuclei in
                      # areas (2+ neighboring windows off by >= regcheck.OFF_UM), or where a large shift is found
READ_SHARE = 0.15     # the 600 um windows with a clear reading (a unique offset, or none): below this share the alignment
                      # itself is wrong (measured: about 0.02 for wrong candidates, 0.5-0.9 for right ones, torn tissue too)
THUMB = 1000          # longest side of the image thumbnail for the search (pixels)
SCALES = np.exp(np.linspace(np.log(0.1), np.log(2.0), 40))  # image um per pixel tried when the file does not say


def _orient(k):
    """2x2 matrix of orientation k (0-7): rotations by 0/90/180/270 degrees, k >= 4 mirrored first."""
    R = np.array([[1, 0], [0, 1]], float)
    for _ in range(k % 4):
        R = np.array([[0, -1], [1, 0]]) @ R
    return R @ (np.diag([-1.0, 1.0]) if k >= 4 else np.eye(2))


ORIENT_NAMES = ["as is", "rotated 90", "rotated 180", "rotated 270", "mirrored", "mirrored, rotated 90",
                "mirrored, rotated 180", "mirrored, rotated 270"]


def image_um_per_px(path):
    """Pixel size stored in an OME-TIFF, or None."""
    try:
        import tifffile, re
        with tifffile.TiffFile(path) as t:
            ome = t.ome_metadata or ""
        m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', ome)
        unit = re.search(r'PhysicalSizeXUnit="([^"]+)"', ome)
        if m and (not unit or unit.group(1) in ("µm", "um", "micron")):
            return float(m.group(1))
    except Exception:
        pass
    return None


def _image(s):
    m = importer._meta(s.sample_id)
    path = m.get("image_pending") or m.get("he_path")
    kind = (m.get("image_pending_kind") or "H&E") if m.get("image_pending") else m.get("image_kind", "H&E")
    if not path or kind == "none":
        raise ValueError("this sample has no image yet: import one first (Import > Image)")
    return path, kind


def _tissue(img, kind):
    """Soft tissue mask (0..1) of a thumbnail."""
    a = np.asarray(img).astype(np.float32)
    if kind == "H&E":
        gray = a.mean(2)
        sat = a.max(2) - a.min(2)
        t = np.clip((235 - gray) / 40, 0, 1) * 0.6 + np.clip((sat - 12) / 30, 0, 1) * 0.6
    else:
        v = ndimage.gaussian_filter(a[..., 0], 2)
        thr = np.percentile(v, 50) + 0.5 * (np.percentile(v, 99) - np.percentile(v, 50))
        t = np.clip(v / max(thr, 1e-6), 0, 1)
    t = ndimage.gaussian_filter(np.clip(t, 0, 1), 1.5)
    return np.clip(t, 0, 1)


def _footprint(u, shape):
    """Soft footprint (0..1) of points u (thumbnail pixels, >= 0) on a grid of `shape`."""
    g = np.zeros(shape, np.float32)
    iu = np.clip(u.astype(int), 0, [shape[1] - 1, shape[0] - 1])
    np.add.at(g, (iu[:, 1], iu[:, 0]), 1)
    g = ndimage.gaussian_filter(g, 1.5)
    return np.clip(g / max(np.percentile(g[g > 0], 30) if (g > 0).any() else 1, 1e-6), 0, 1)


def search(s, scales=None, top=6, progress=None):
    """Candidate placements of the cells on the image, best first (stored for show / refine / apply)."""
    say = progress or (lambda f, m: None)
    path, kind = _image(s)
    sl = render.slide(path)
    W0, H0 = sl.W0, sl.H0
    f = max(W0, H0) / THUMB  # level-0 pixels per thumbnail pixel
    w, h = int(round(W0 / f)), int(round(H0 / f))
    I = _tissue(sl.read(0, 0, W0, H0, w, h), kind)
    J = 2 * I - 1
    known = image_um_per_px(path)
    if scales is None:
        scales = known * np.exp(np.linspace(-0.08, 0.08, 9)) if known else SCALES
    xy = s.cells[["x", "y"]].to_numpy()
    xy = xy[np.random.default_rng(0).choice(len(xy), min(len(xy), 200000), replace=False)]
    found = []
    say(0.05, f"searching {len(scales)} pixel sizes x 8 orientations")
    for k in range(8):
        O = _orient(k)
        for sc in scales:
            u = (xy @ O.T) / (sc * f)
            umin = u.min(0)
            u = u - umin
            th, tw = int(u[:, 1].max()) + 2, int(u[:, 0].max()) + 2
            if th > 1.6 * h or tw > 1.6 * w or th < 8 or tw < 8:
                continue
            K = 2 * _footprint(u, (th, tw)) - 1
            P = (h + th, w + tw)
            Jp = -np.ones(P, np.float32); Jp[:h, :w] = J  # outside the image is background
            C = np.fft.irfft2(np.fft.rfft2(Jp) * np.conj(np.fft.rfft2(K, P)), P)
            C = C / K.size
            C = np.roll(np.roll(C, 0, 0), 0, 1)
            i, j = np.unravel_index(np.argmax(C), C.shape)
            dy, dx = (i if i <= h else i - P[0]), (j if j <= w else j - P[1])
            found.append({"orientation": k, "um_per_px": float(sc), "score": float(C[i, j]),
                          "A": np.c_[O / sc, f * (np.array([dx, dy]) - umin)].tolist()})
        say(0.1 + 0.8 * (k + 1) / 8, f"orientation {k + 1} of 8")
    found.sort(key=lambda x: -x["score"])
    out = []
    for c in found:  # one candidate per orientation and place
        A = np.asarray(c["A"])
        if any(c["orientation"] == o["orientation"] and np.hypot(*(A[:, 2] - np.asarray(o["A"])[:, 2])) < 0.1 * max(W0, H0)
               for o in out):
            continue
        out.append(c)
        if len(out) >= top:
            break
    for i, c in enumerate(out):
        c.update(id=f"K{i + 1}", orientation_name=ORIENT_NAMES[c["orientation"]], score=round(c["score"], 3),
                 um_per_px=round(c["um_per_px"], 4))
    st = {"image": path, "kind": kind, "image_um_per_px_in_file": known, "candidates": out, "time": time.time()}
    json.dump(st, open(_file(s), "w"), indent=1)
    say(1.0, "done")
    return st


def _file(s):
    d = config.SAMPLES / s.sample_id / "align"
    d.mkdir(exist_ok=True)
    return d / "candidates.json"


def _cands(s):
    if not _file(s).exists():
        raise ValueError("no candidates yet: run the search first")
    return json.load(open(_file(s)))


def _get(s, cid):
    st = _cands(s)
    hit = [c for c in st["candidates"] if c["id"] == str(cid).upper()]
    if not hit:
        raise ValueError(f"no candidate {cid}; candidates: {[c['id'] for c in st['candidates']]}")
    return st, hit[0]


def _window_offsets(sl, kind, pxy, um_per_px, n_windows=30, window_um=300.0, max_shift_um=25.0):
    """Offsets (um) between nuclei on the image and cells, in windows of dense tissue: [(cx, cy, dx, dy, corr)]."""
    half = window_um / um_per_px / 2
    W = int(window_um)
    R = int(max_shift_um)
    tile = window_um / um_per_px
    gx, gy = (pxy[:, 0] // tile).astype(int), (pxy[:, 1] // tile).astype(int)
    key = gx * 100000 + gy
    vals, cnt = np.unique(key, return_counts=True)
    dense = vals[cnt >= 0.4 * cnt.max()]
    pick = np.random.default_rng(0).choice(dense, min(n_windows, len(dense)), replace=False)
    out = []
    for k in pick:
        cx, cy = (k // 100000 + 0.5) * tile, (k % 100000 + 0.5) * tile
        if cx - half < 0 or cy - half < 0 or cx + half > sl.W0 or cy + half > sl.H0:
            continue
        img = np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, W, W)).astype(np.float32)
        nuc = ndimage.gaussian_filter(insights._nuclei(img, kind), 2)
        m = (pxy[:, 0] > cx - half) & (pxy[:, 0] < cx + half) & (pxy[:, 1] > cy - half) & (pxy[:, 1] < cy + half)
        if m.sum() < 80:
            continue
        pts = np.zeros((W, W))
        u = ((pxy[m, 0] - (cx - half)) / (2 * half) * W).astype(int).clip(0, W - 1)
        v = ((pxy[m, 1] - (cy - half)) / (2 * half) * W).astype(int).clip(0, W - 1)
        np.add.at(pts, (v, u), 1)
        pts = ndimage.gaussian_filter(pts, 3)
        nuc, pts = (nuc - nuc.mean()) / (nuc.std() + 1e-9), (pts - pts.mean()) / (pts.std() + 1e-9)
        cc = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(nuc) * np.conj(np.fft.fft2(pts))))) / W ** 2
        c0 = W // 2
        win = cc[c0 - R:c0 + R + 1, c0 - R:c0 + R + 1]
        dy, dx = np.unravel_index(np.argmax(win), win.shape)
        out.append((cx, cy, float(dx - R), float(dy - R), float(win.max())))
    return out


def refine(s, cid, rounds=(60.0, 25.0, 12.0)):
    """Bring a candidate to micrometer precision (see the module doc); stores the refined affine on the candidate."""
    st, c = _get(s, cid)
    sl = render.slide(st["image"])
    xy = s.cells[["x", "y"]].to_numpy()
    A = np.asarray(c["A"], float)
    log = []
    for R in rounds:
        upp = float(1 / np.sqrt(abs(np.linalg.det(A[:, :2]))))
        pxy = xy @ A[:, :2].T + A[:, 2]
        win = _window_offsets(sl, st["kind"], pxy, upp, window_um=max(300.0, 5 * R), max_shift_um=R)
        if len(win) < 6:
            raise ValueError("too few windows of dense tissue overlap the image with this candidate: try another")
        P = np.array([[w[0], w[1]] for w in win])
        D = np.array([[w[2], w[3]] for w in win]) / upp           # offsets in image pixels
        X = np.c_[P - P.mean(0), np.ones(len(P))]
        coef = np.linalg.lstsq(X, D, rcond=None)[0]               # offset = (p - mean) @ B + b
        B, b = coef[:2].T, coef[2]
        resid = float(np.median(np.hypot(*((D - X @ coef) * upp).T)))
        M = np.eye(2) + B                                         # new px = M (px - mean) + mean + b
        A = np.c_[M @ A[:, :2], M @ A[:, 2] - M @ P.mean(0) + P.mean(0) + b]
        log.append({"search_um": R, "windows": len(win), "median_shift_um": round(float(np.median(np.hypot(*(D * upp).T))), 1),
                    "fit_residual_um": round(resid, 1), "correlation": round(float(np.median([w[4] for w in win])), 3)})
    c["A_refined"] = A.tolist()
    c["refine_log"] = log
    json.dump(st, open(_file(s), "w"), indent=1)
    d = deviation(s, cid)
    return {"id": c["id"], "rounds": log, "um_per_px": round(float(1 / np.sqrt(abs(np.linalg.det(A[:, :2])))), 4),
            "deviation": {k: v for k, v in d.items() if k not in ("A", "offsets_um")}, "deviation_advice": advice_lines(d),
            "note": "a good alignment ends with median shifts of a few um and correlations clearly above 0.2; deviation: "
                    "the whole tissue checked (as the review window of apply_alignment shows it): 'far off' is not the "
                    "right alignment, 'local shifts' needs the local registration after it is applied"}


def show(s, cid, refined=True):
    """Figure of a candidate: the whole image with the cells' footprint outlined (left) and two close-ups of dense
    tissue with the cell centroids on the image (right)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    st, c = _get(s, cid)
    A = np.asarray(c.get("A_refined") if refined and c.get("A_refined") else c["A"], float)
    sl = render.slide(st["image"])
    W0, H0 = sl.W0, sl.H0
    f = max(W0, H0) / THUMB
    w, h = int(round(W0 / f)), int(round(H0 / f))
    xy = s.cells[["x", "y"]].to_numpy()
    pxy = xy @ A[:, :2].T + A[:, 2]
    fig = plt.figure(figsize=(14, 6.5))
    ax = fig.add_axes([0.01, 0.02, 0.48, 0.9])
    ax.imshow(np.asarray(sl.read(0, 0, W0, H0, w, h)), extent=(0, W0, H0, 0))
    sub = np.random.default_rng(0).choice(len(pxy), min(len(pxy), 40000), replace=False)
    ax.scatter(pxy[sub, 0], pxy[sub, 1], s=0.3, c="#00c8ff", alpha=0.35, linewidths=0, rasterized=True)
    ax.set_xlim(0, W0); ax.set_ylim(H0, 0); ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(f"{c['id']}: {c['orientation_name']}, {1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f} um/px, "
                 f"score {c['score']}" + (" (refined)" if refined and c.get("A_refined") else ""), fontsize=9)
    upp = float(1 / np.sqrt(abs(np.linalg.det(A[:, :2]))))
    tile = 200 / upp
    key = (pxy[:, 0] // tile).astype(int) * 100000 + (pxy[:, 1] // tile).astype(int)
    vals, cnt = np.unique(key, return_counts=True)
    for i, k in enumerate(vals[np.argsort(-cnt)][[0, len(vals) // 6]]):
        cx, cy = (k // 100000 + 0.5) * tile, (k % 100000 + 0.5) * tile
        half = 60 / upp
        a2 = fig.add_axes([0.51 + 0.245 * i, 0.08, 0.235, 0.8])
        a2.imshow(np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, 480, 480)), extent=(cx - half, cx + half, cy + half, cy - half))
        m = (np.abs(pxy[:, 0] - cx) < half) & (np.abs(pxy[:, 1] - cy) < half)
        a2.scatter(pxy[m, 0], pxy[m, 1], s=14, facecolors="none", edgecolors="#00e6ff", linewidths=1.0)
        a2.set_xlim(cx - half, cx + half); a2.set_ylim(cy + half, cy - half); a2.set_xticks([]); a2.set_yticks([])
        a2.set_title("close-up 120 um: cell centroids (cyan) should sit on nuclei", fontsize=8)
    out = config.SAMPLES / s.sample_id / "align" / f"{c['id']}{'_refined' if refined and c.get('A_refined') else ''}.jpg"
    plt.savefig(out, dpi=90, pil_kwargs={"quality": 85}); plt.close(fig)
    return str(out)


def _affine(c):
    return np.asarray(c.get("A_refined") or c["A"], float)


def _wfile(s, cid):
    return config.SAMPLES / s.sample_id / "align" / f"{cid}_windows.json"


def deviation(s, cid, progress=None):
    """How far the cells sit from their nuclei with candidate cid (refined if it was), over the whole tissue: the
    offsets of regcheck's windows (fine and coarse) with the cells where the candidate puts them, a verdict ("fine",
    "local shifts": advise the local registration after the alignment, "far off": not the right alignment) and the
    numbers behind it. Kept on the candidate (and its windows next to it) until the candidate changes."""
    st, c = _get(s, cid)
    A = _affine(c)
    if (c.get("deviation") or {}).get("A") == A.tolist() and _wfile(s, cid).exists():
        return c["deviation"]
    sl = render.slide(st["image"])
    pxy = s.cells[["x", "y"]].to_numpy() @ A[:, :2].T + A[:, 2]
    px, py = pxy[:, 0], pxy[:, 1]
    upp = float(1 / np.sqrt(abs(np.linalg.det(A[:, :2]))))
    x0, y0 = float(px.min()), float(py.min())
    tile, ct = regcheck.WINDOW_UM / upp, regcheck.COARSE_WINDOW_UM / upp
    cols = ["i", "j", "cx", "cy", "dx_um", "dy_um", "corr", "n_cells", "corr_at_zero"]
    w = {"image": st["image"], "kind": st["kind"], "um_per_px": upp, "window_um": regcheck.WINDOW_UM,
         "search_um": regcheck.SEARCH_UM, "origin_px": [x0, y0], "tile_px": tile, "nx": int((px.max() - x0) // tile) + 1,
         "time": time.strftime("%Y-%m-%dT%H:%M:%S"), "columns": cols}
    say = progress or (lambda f: None)
    w["windows"] = regcheck._scan(sl, st["kind"], px, py, w, lambda f: say(0.8 * f))
    w["coarse"] = {"window_um": regcheck.COARSE_WINDOW_UM, "search_um": regcheck.COARSE_SEARCH_UM, "tile_px": ct,
                   "origin_px": [x0, y0], "nx": int((px.max() - x0) // ct) + 1, "columns": cols + ["second"]}
    w["coarse"]["windows"] = regcheck._scan(sl, st["kind"], px, py, w["coarse"] | {"um_per_px": upp},
                                            lambda f: say(0.8 + 0.2 * f), coarse=True)
    if not w["windows"]:
        raise ValueError("no window of dense tissue overlaps the image with this candidate: try another")
    a, mu, coef, aff, res, noise, D = regcheck._analyse(w)
    tot, raw = np.hypot(*D.T), np.hypot(a[:, 4], a[:, 5])
    trusted = ~noise
    areas = [x for x in regcheck._areas(a, D, noise) if len(x) >= 2]
    in_areas = int(sum(len(x) for x in areas))
    b, big, anchor = regcheck._large(w)
    read = float((big | anchor).mean()) if len(b) else 0.0
    beyond = float((np.abs(a[:, 4:6]).max(1) >= w["search_um"]).mean())
    share_areas = in_areas / max(int(trusted.sum()), 1)
    verdict = ("far off" if read < READ_SHARE else
               "local shifts" if share_areas >= LOCAL_SHARE or big.any() or beyond > 0.3 else "fine")
    t = raw[trusted]
    out = {"verdict": verdict, "windows": int(len(a)), "trusted_windows": int(trusted.sum()), "window_um": regcheck.WINDOW_UM,
           "median_um": round(float(np.median(t)), 1) if len(t) else None,
           "p90_um": round(float(np.percentile(t, 90)), 1) if len(t) else None,
           "share_off": round(float(((tot >= regcheck.OFF_UM) & trusted).sum() / max(int(trusted.sum()), 1)), 3),
           "off_um": regcheck.OFF_UM, "areas_off": len(areas), "share_in_areas": round(share_areas, 3),
           "largest_area_windows": max((len(x) for x in areas), default=0),
           "large_shift_windows": int(big.sum()), "large_shift_um": round(float(np.hypot(b[big, 4], b[big, 5]).max()), 0) if big.any() else None,
           "beyond_search_share": round(beyond, 3), "clear_reading_share": round(read, 3), "affine_part_max_um": round(float(np.hypot(*aff.T).max()), 1),
           "offsets_um": [round(float(x), 1) for x in t], "A": A.tolist()}
    json.dump(w, open(_wfile(s, cid), "w"))
    st, c = _get(s, cid)   # (as stored now)
    c["deviation"] = out
    json.dump(st, open(_file(s), "w"), indent=1)
    return out


def advice_lines(d):
    """The deviation in words, for the review window and the agent."""
    lines = [f"Deviation of the cells from the nuclei on the image ({d['trusted_windows']} windows of "
             f"{d['window_um']:g} um over the tissue): median {d['median_um']} um, 90% within {d['p90_um']} um; "
             f"{d['share_off']:.0%} of the windows off by {d['off_um']:g} um or more."]
    if d["verdict"] == "far off":
        lines.append(f"Far off: only {d['clear_reading_share']:.0%} of the {regcheck.COARSE_WINDOW_UM:g} um windows find "
                     f"the nuclei within {regcheck.COARSE_SEARCH_UM} um of the cells. This does not look like the right "
                     "alignment: try another candidate rather than accepting it.")
    elif d["verdict"] == "local shifts":
        what = []
        if d["areas_off"]:
            what.append(f"{d['areas_off']} area(s) off their nuclei ({d['share_in_areas']:.0%} of the tissue's windows)")
        if d["beyond_search_share"] > 0.3:
            what.append(f"{d['beyond_search_share']:.0%} of the windows off by more than {regcheck.SEARCH_UM} um")
        if d["large_shift_windows"]:
            what.append(f"{d['large_shift_windows']} window(s) of {regcheck.COARSE_WINDOW_UM:g} um shifted up to "
                        f"{d['large_shift_um']:g} um")
        lines.append("Shifts too large for one alignment: " + "; ".join(what) + " (tissue stretched, folded or torn "
                     "between the run and the staining). Accept, then do the local registration: Steps > Check "
                     "registration (Review: Spatial organization); the agent measures the shifts and offers the local "
                     "correction.")
    else:
        lines.append("No area sits off its nuclei: no local registration needed.")
    return lines


def deviation_figure(s, cid):
    """Histogram of the windows' deviation (offset of the nuclei from the cells, um) and the windows on the tissue
    colored by it (large shifts dashed)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    d = deviation(s, cid)
    w = json.load(open(_wfile(s, cid)))
    a, mu, coef, aff, res, noise, D = regcheck._analyse(w)
    b, big, _ = regcheck._large(w)
    raw = np.hypot(a[:, 4], a[:, 5])
    off = regcheck.OFF_UM
    fig = plt.figure(figsize=(13, 4.6))
    ax = fig.add_axes([0.06, 0.13, 0.42, 0.72])
    t = raw[~noise]
    bins = np.arange(0, regcheck.SEARCH_UM + 1, 1.0)
    ax.hist(t[t < off], bins=bins, color="#4c9a6a", label=f"on their nuclei (< {off:g} um)")
    ax.hist(t[t >= off], bins=bins, color="#d9822b", label=f"off (>= {off:g} um)")
    ax.axvline(off, color="#555", lw=1, ls="--")
    ax.set_xlim(0, regcheck.SEARCH_UM)
    ax.set_xlabel(f"deviation: nuclei on the image vs cell centroids, per {w['window_um']:g} um window (um)")
    ax.set_ylabel("windows")
    ax.legend(frameon=False, fontsize=9, loc="best")
    extra = (f"; + {d['large_shift_windows']} large shift(s) of {regcheck.COARSE_WINDOW_UM:g} um windows up to "
             f"{d['large_shift_um']:g} um") if d["large_shift_windows"] else ""
    ax.set_title(f"median {d['median_um']} um, {d['share_off']:.0%} of windows >= {off:g} um{extra}", fontsize=9.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax2 = fig.add_axes([0.52, 0.04, 0.42, 0.84])
    sl = render.slide(w["image"])
    half = w["tile_px"] / 2
    X0, Y0 = a[:, 2].min() - 3 * half, a[:, 3].min() - 3 * half
    X1, Y1 = a[:, 2].max() + 3 * half, a[:, 3].max() + 3 * half
    X0, Y0, X1, Y1 = max(X0, 0), max(Y0, 0), min(X1, sl.W0), min(Y1, sl.H0)
    f = max(X1 - X0, Y1 - Y0) / 900
    ax2.imshow(np.asarray(sl.read(X0, Y0, X1, Y1, int((X1 - X0) / f), int((Y1 - Y0) / f))), extent=(X0, X1, Y1, Y0), alpha=0.55)
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("dev", [(0, "#4c9a6a"), (off / 15, "#f2d15c"), (0.6, "#d9822b"), (1, "#a01010")])
    for n in range(len(a)):
        if noise[n]:
            ax2.add_patch(Rectangle((a[n, 2] - half, a[n, 3] - half), 2 * half, 2 * half, fc="none", ec="#999", lw=0.3))
        else:
            v = min(raw[n] / 15.0, 1.0)  # (the offsets beyond the fine search: 15+)
            ax2.add_patch(Rectangle((a[n, 2] - half, a[n, 3] - half), 2 * half, 2 * half, fc=cmap(v), ec="none", alpha=0.75))
    hb = w["coarse"]["tile_px"] / 2
    for n in np.where(big)[0]:
        ax2.add_patch(Rectangle((b[n, 2] - hb, b[n, 3] - hb), 2 * hb, 2 * hb, fc="none", ec="#b00000", lw=1.4, ls="--"))
    ax2.set_xlim(X0, X1); ax2.set_ylim(Y1, Y0); ax2.set_xticks([]); ax2.set_yticks([]); ax2.set_aspect("equal")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 15))
    cb = fig.colorbar(sm, ax=ax2, fraction=0.035, pad=0.01)
    cb.set_label("deviation (um)", fontsize=8)
    cb.set_ticks([0, off, 10, 15]); cb.set_ticklabels(["0", f"{off:g}", "10", "15+"])
    ax2.set_title("windows on the tissue (gray outline: no clear reading; dashed: large shift)", fontsize=9)
    head = {"fine": "Cells on their nuclei: no local registration needed",
            "local shifts": "Some areas sit off their nuclei: do the local registration after accepting",
            "far off": "Far off: probably not the right alignment"}[d["verdict"]]
    fig.suptitle(f"{c_id(cid)}: {head}", fontsize=11, fontweight="bold",
                 color={"fine": "#2d6a45", "local shifts": "#9a5200", "far off": "#a01010"}[d["verdict"]])
    out = config.SAMPLES / s.sample_id / "align" / f"{c_id(cid)}_deviation.jpg"
    plt.savefig(out, dpi=90, pil_kwargs={"quality": 88}); plt.close(fig)
    return str(out)


def c_id(cid):
    return str(cid).upper()


def apply(s, cid):
    """Put the cells on the image with candidate cid (refined if it was); its deviation check, when made for this
    affine, becomes the registration check of the placement (and the sample notes its verdict for the steps)."""
    st, c = _get(s, cid)
    A = _affine(c)
    r = importer.place_cells(s.sample_id, A.tolist(), st["image"], st["kind"], by="Claude", progress=lambda *a: None)
    d = c.get("deviation")
    if d and d.get("A") == A.tolist() and _wfile(s, cid).exists():
        image = importer.image_of(s.sample_id)
        key = regcheck._key(s.sample_id, image)
        w = json.load(open(_wfile(s, cid)))
        json.dump(dict(w, key=key, image=image), open(regcheck._file(s.sample_id), "w"))
        m = importer._meta(s.sample_id)
        m["alignment_check"] = {"key": key, **{k: v for k, v in d.items() if k not in ("A", "offsets_um")},
                                "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
        json.dump(m, open(config.SAMPLES / s.sample_id / "meta.json", "w"), indent=1)
        r["deviation"] = {"verdict": d["verdict"], "advice": advice_lines(d)[1:]}
    return r


class Sample:
    """What the alignment needs of a sample (it may have no annotation yet): its id and its cells."""
    def __init__(self, sample_id):
        import pandas as pd
        self.sample_id = sample_id
        self.cells = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["x", "y", "px", "py"])


def status(sample_id):
    """What the sample has: the image shown, an image waiting for alignment, the registration and who made it."""
    m = importer._meta(sample_id)
    c = Sample(sample_id).cells
    return {"sample": sample_id, "image_shown": m.get("he_path") if m.get("image_kind") != "none" else None,
            "image_kind": m.get("image_kind"), "image_waiting_for_alignment": m.get("image_pending"),
            "registration": {"file": bool(m.get("registration_affine")), "placed_by": m.get("registration_by")},
            "can_undo": bool(m.get("placement_previous")),
            "cells_extent_um": [round(float(c.x.max() - c.x.min())), round(float(c.y.max() - c.y.min()))],
            "image_um_per_px_in_file": image_um_per_px(m.get("image_pending") or m.get("he_path") or "")}
