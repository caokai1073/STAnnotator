"""Aligning a sample's cells to its image without a registration file (the app's "Align with Claude").

1. search(): the tissue on a thumbnail of the image is compared with the footprint of the cells (where cells are,
   from their um coordinates) for 8 orientations (4 rotations, mirrored or not) and a range of image pixel sizes
   (around the size stored in an OME-TIFF, else 0.1-2 um/px), by FFT cross-correlation. The score is the agreement
   of tissue and footprint over the footprint's box (+1 both or neither, -1 one only), so a footprint that is too
   small or too large scores low. The best few are kept as candidates (affines cell um -> image level-0 px).
2. Claude looks at the candidates (show()), picks one (several sections on one slide, symmetric tissue).
3. refine(): moves the candidate to micrometer precision: nuclei on the image vs cell centroids in windows over the
   tissue (as insights.registration), an affine fit of the window offsets, three rounds with a shrinking search.
4. apply(): puts the cells on the image (importer.place_cells; undo with importer.undo_placement).
"""
import json, time
import numpy as np
from scipy import ndimage
from . import config, importer, render, insights

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
    return {"id": c["id"], "rounds": log, "um_per_px": round(float(1 / np.sqrt(abs(np.linalg.det(A[:, :2])))), 4),
            "note": "a good alignment ends with median shifts of a few um and correlations clearly above 0.2"}


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


def apply(s, cid):
    """Put the cells on the image with candidate cid (refined if it was)."""
    st, c = _get(s, cid)
    A = c.get("A_refined") or c["A"]
    return importer.place_cells(s.sample_id, A, st["image"], st["kind"], by="Claude", progress=lambda *a: None)


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
