"""Compartments of the tissue read from the H&E itself, for outlining regions without cells (a sample of an H&E
only) or before an annotation. The agent reads the image classes together with the organ's priors
(data/priors/tissue_priors.md: how each compartment looks on the H&E) and proposes outlines; nothing is applied
without the user's Accept.

segment(): a thumbnail of the whole image at about STEP_UM per pixel. Tissue = pixels darker than the background
  (optical density above TISSUE_OD). Color deconvolution into hematoxylin and eosin (stain vectors of Ruifrok and
  Johnston 2001); tissue = stained pixels (colored, not the gray of the glass) without dust. Features of every tissue pixel: hematoxylin and eosin smoothed at two scales, the texture of the
  hematoxylin (local standard deviation: nuclei-rich epithelium and lymphoid tissue are grainy, smooth muscle is not);
  smoothed over the tissue only, so the glass does not pull the edges' values. k-means of the tissue pixels into n classes,
  numbered by their hematoxylin (1 = the darkest blue), then smoothed (the class most common in a small window).
  Per class: area, mean hematoxylin, eosin and texture, the distance to the background (tissue edge or a lumen),
  the share next to it.
outline(): the pixels of chosen classes (optionally only within a box or a polygon), thin seams opened, closed over
  gaps of CLOSE_UM, holes under MAX_HOLE_UM2 filled, pieces under min_area dropped, traced into polygons (contourpy) and simplified, in the image's
  full-resolution pixels.
view(): the image, or a box of it, with a grid labelled in full-resolution pixels, so an outline can be given point
  by point.
Results are kept in the sample's folder (he_segment/) for outline() and the review window.
"""
import json, re, time
import numpy as np
from . import config

STEP_UM = 8.0          # thumbnail resolution
MAX_SIDE = 3000        # thumbnail pixels on its longer side, at most
TISSUE_OD = 0.05       # summed optical density above this (and colored: MIN_SAT) is tissue
MIN_SAT = 0.05         # saturation of stained tissue, at least (the gray of glass is near 0; pale submucosa ~0.1)
MIN_SPECK_UM2 = 20000.0  # tissue pieces smaller than this are dust
MIN_HOLE_UM2 = 3000.0    # holes in the tissue smaller than this are filled; larger ones (lumens) stay out
CLOSE_UM = 40.0        # outline(): gaps narrower than this are closed
MIN_AREA_UM2 = 20000.0 # outline(): smaller pieces are dropped
MAX_HOLE_UM2 = 50000.0 # outline(): holes smaller than this are filled
FIT_PIXELS = 200000    # k-means is fitted on a random sample of the tissue pixels
SEED = 0
# hematoxylin, eosin and residual stain vectors (optical density, RGB), Ruifrok & Johnston 2001
STAINS = np.array([[0.65, 0.70, 0.29], [0.07, 0.99, 0.11], [0.27, 0.57, 0.78]])
PALETTE = ["#1f4e9c", "#c0392b", "#27ae60", "#8e44ad", "#e67e22", "#16a085", "#7f8c8d", "#d4ac0d", "#2c3e50", "#e84393"]


def _dir(sample_id):
    d = config.SAMPLES / sample_id / "he_segment"
    d.mkdir(exist_ok=True)
    return d


def _image(sample_id):
    from .importer import image_of, _meta
    path = image_of(sample_id)
    m = _meta(sample_id)
    if not path or m.get("image_kind") in (None, "none"):
        raise ValueError("this sample shows no image (none, still waiting for its alignment, or unloaded)")
    return path, float(m.get("he_um_per_px") or 0.5), m.get("image_kind", "H&E"), m.get("um_per_px_known", True)


def thumbnail(sample_id, step_um=STEP_UM, box=None):
    """(RGB array, level-0 pixels per thumbnail pixel, box) of the image or of box = (x0, y0, x1, y1) in level-0 px."""
    from .render import slide
    path, upp, _, _ = _image(sample_id)
    s = slide(str(path))
    x0, y0, x1, y1 = box or (0, 0, s.W0, s.H0)
    scale = max(step_um / upp, (x1 - x0) / MAX_SIDE, (y1 - y0) / MAX_SIDE, 1.0)
    w, h = max(1, int((x1 - x0) / scale)), max(1, int((y1 - y0) / scale))
    return np.asarray(s.read(x0, y0, x1, y1, w, h).convert("RGB")), scale, (x0, y0, x1, y1)


def _deconvolve(rgb):
    od = -np.log10(np.clip(rgb.astype(np.float32), 1, 255) / 255.0)
    M = STAINS / np.linalg.norm(STAINS, axis=1, keepdims=True)
    c = od.reshape(-1, 3) @ np.linalg.inv(M)
    return od.sum(2), np.clip(c[:, 0], 0, None).reshape(od.shape[:2]), np.clip(c[:, 1], 0, None).reshape(od.shape[:2])


def _tissue(rgb, odsum, step):
    """Stained tissue: colored (saturation above MIN_SAT) and not white; specks under MIN_SPECK_UM2
    dropped, holes under MIN_HOLE_UM2 filled (lumens and glass stay out)."""
    from scipy import ndimage as ndi
    a = rgb.astype(np.float32)
    mx, mn = a.max(2), a.min(2)
    sat = ndi.gaussian_filter((mx - mn) / np.maximum(mx, 1), 1)
    t = sat > MIN_SAT
    t &= ndi.gaussian_filter(odsum, 1) > TISSUE_OD
    t = ndi.binary_opening(t, iterations=1)
    comp, n = ndi.label(t)
    if n:
        sizes = ndi.sum(t, comp, range(1, n + 1)) * step ** 2
        t = np.isin(comp, 1 + np.nonzero(sizes >= MIN_SPECK_UM2)[0])
    holes = ndi.binary_fill_holes(t) & ~t
    comp, n = ndi.label(holes)
    if n:
        sizes = ndi.sum(holes, comp, range(1, n + 1)) * step ** 2
        t |= np.isin(comp, 1 + np.nonzero(sizes < MIN_HOLE_UM2)[0])
    return t


def segment(sample_id, n_classes=6, progress=None):
    """Classes of the tissue on the image (see the module). Returns the stats; the label image is kept."""
    from scipy import ndimage as ndi
    from sklearn.cluster import KMeans
    say = progress or (lambda m: None)
    n_classes = int(n_classes)
    if not 2 <= n_classes <= len(PALETTE):
        raise ValueError(f"n_classes must be 2 to {len(PALETTE)}")
    t0 = time.time()
    rgb, scale, _ = thumbnail(sample_id)
    path, upp, kind, known = _image(sample_id)
    if kind != "H&E":
        raise ValueError("the image is not an H&E (DAPI / fluorescence): its colors cannot be read as stains")
    step = scale * upp  # um per thumbnail pixel
    say(f"thumbnail {rgb.shape[1]} x {rgb.shape[0]} px at {step:.1f} um per pixel")
    odsum, H, E = _deconvolve(rgb)
    px = lambda um: max(1.0, um / step)
    tissue = _tissue(rgb, odsum, step)
    if tissue.sum() < 100:
        raise ValueError("almost no tissue found on the image")
    tw = tissue.astype(np.float32)
    # smoothing over the tissue only (normalized convolution): the glass around it does not pull the edges' values
    sm = lambda a, um: ndi.gaussian_filter(a * tw, px(um)) / np.maximum(ndi.gaussian_filter(tw, px(um)), 1e-3)
    feats = [sm(H, 10), sm(H, 40), sm(E, 10), sm(E, 40)]
    hm = sm(H, 8)
    tex = np.sqrt(np.maximum(sm(H ** 2, 8) - hm ** 2, 0))
    feats.append(tex)
    dist = ndi.distance_transform_edt(tissue) * step  # a statistic of the classes, not a feature (it would split
    F = np.stack([f[tissue] for f in feats], 1).astype(np.float32)  # every edge off as a class of its own)
    mu, sd = F.mean(0), F.std(0) + 1e-6
    Fz = (F - mu) / sd
    rng = np.random.default_rng(SEED)
    fit = Fz[rng.choice(len(Fz), min(FIT_PIXELS, len(Fz)), replace=False)]
    say(f"k-means of {len(Fz):,} tissue pixels into {n_classes} classes")
    km = KMeans(n_classes, n_init=4, random_state=SEED).fit(fit)
    lab0 = km.predict(Fz)
    order = np.argsort(-km.cluster_centers_[:, 0])  # class 1 = the most hematoxylin
    rank = np.empty(n_classes, int)
    rank[order] = np.arange(n_classes)
    lab = np.full(tissue.shape, -1, np.int16)
    lab[tissue] = rank[lab0]
    # smoothing: the class most common in a window of about 24 um
    votes = np.stack([ndi.uniform_filter((lab == k).astype(np.float32), int(px(24)) | 1) for k in range(n_classes)])
    lab = np.where(tissue, votes.argmax(0), -1).astype(np.int16)
    edge = tissue & ~ndi.binary_erosion(tissue, iterations=max(1, int(px(30))))
    classes = []
    for k in range(n_classes):
        m = lab == k
        if not m.any():
            continue
        ys, xs = np.nonzero(m)
        classes.append({"class": k + 1, "color": PALETTE[k], "area_mm2": round(float(m.sum() * step ** 2 / 1e6), 3),
                        "share_of_tissue": round(float(m.sum() / tissue.sum()), 3),
                        "hematoxylin": round(float(H[m].mean()), 3), "eosin": round(float(E[m].mean()), 3),
                        "texture": round(float(tex[m].mean()), 3),
                        "median_um_from_edge": round(float(np.median(dist[m])), 0),
                        "share_at_edge": round(float((m & edge).sum() / m.sum()), 3),
                        "center_px": [round(float(xs.mean() * scale)), round(float(ys.mean() * scale))]})
    d = _dir(sample_id)
    np.save(d / "labels.npy", lab)
    meta = {"time": time.strftime("%Y-%m-%dT%H:%M:%S"), "n_classes": n_classes, "scale_px": scale, "um_per_px": step,
            "image": str(path), "shape": list(lab.shape), "classes": classes, "pixel_size_known": known}
    json.dump(meta, open(d / "meta.json", "w"), indent=1)
    say(f"done in {time.time() - t0:.1f} s")
    return meta


def load(sample_id):
    d = config.SAMPLES / sample_id / "he_segment"
    if not (d / "meta.json").exists():
        raise ValueError("no classes of the H&E yet: run segment_he first")
    return json.load(open(d / "meta.json")), np.load(d / "labels.npy")


def _simplify(pts, tol):
    """Douglas-Peucker of a closed ring (N x 2, the last point = the first): split at the point farthest from the
    first, each half simplified."""
    if len(pts) < 5:
        return pts
    if np.allclose(pts[0], pts[-1]):
        far = int(np.argmax(((pts - pts[0]) ** 2).sum(1)))
        a, b = _simplify_open(pts[:far + 1], tol), _simplify_open(pts[far:], tol)
        return np.vstack([a, b[1:]])
    return _simplify_open(pts, tol)


def _simplify_open(pts, tol):
    if len(pts) < 3:
        return pts
    keep = np.zeros(len(pts), bool)
    keep[0] = keep[-1] = True
    stack = [(0, len(pts) - 1)]
    while stack:
        a, b = stack.pop()
        if b <= a + 1:
            continue
        A, B = pts[a], pts[b]
        dx, dy = B - A
        L = np.hypot(dx, dy) or 1e-9
        seg = pts[a + 1:b]
        dist = np.abs(dy * (seg[:, 0] - A[0]) - dx * (seg[:, 1] - A[1])) / L
        i = int(np.argmax(dist))
        if dist[i] > tol:
            keep[a + 1 + i] = True
            stack += [(a, a + 1 + i), (a + 1 + i, b)]
    return pts[keep]


def _ring_area(p):
    x, y = p[:, 0], p[:, 1]
    return abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1))) / 2


def outline(sample_id, classes=None, box_px=None, within_px=None, min_area_um2=MIN_AREA_UM2, close_um=CLOSE_UM,
            max_hole_um2=MAX_HOLE_UM2):
    """Outlines of the pixels of `classes` (numbers from segment), or of the tissue when classes is None, only within
    box_px = [x0, y0, x1, y1] and / or the polygon within_px: [{"outer": ring, "holes": [ring, ...]}] in level-0
    pixels (as regions.parse_geojson gives them; a hole is other tissue or a lumen inside the outline)."""
    import contourpy
    from scipy import ndimage as ndi
    from matplotlib.path import Path
    meta, lab = load(sample_id)
    scale, step = meta["scale_px"], meta["um_per_px"]
    m = lab >= 0 if not classes else np.isin(lab, [int(c) - 1 for c in classes])
    if box_px is not None:
        x0, y0, x1, y1 = [float(v) / scale for v in box_px]
        keep = np.zeros_like(m)
        keep[max(0, int(y0)):max(0, int(np.ceil(y1))), max(0, int(x0)):max(0, int(np.ceil(x1)))] = True
        m &= keep
    if within_px is not None:
        h, w = m.shape
        yy, xx = np.mgrid[0:h, 0:w]
        inside = Path(np.asarray(within_px, float) / scale).contains_points(np.c_[xx.ravel() + 0.5, yy.ravel() + 0.5])
        m &= inside.reshape(h, w)
    m = ndi.binary_opening(m, iterations=2)  # seams up to about 4 thumbnail pixels go (blended pixels along edges)
    r = max(1, int(round(close_um / step / 2)))
    m = ndi.binary_closing(np.pad(m, r + 1), iterations=r)[r + 1:-r - 1, r + 1:-r - 1]
    holes = ndi.binary_fill_holes(m) & ~m  # small holes filled; large ones (other tissue, lumens) stay out
    comp, n = ndi.label(holes)
    if n:
        sizes = ndi.sum(holes, comp, range(1, n + 1)) * step ** 2
        m |= np.isin(comp, 1 + np.nonzero(sizes < max_hole_um2)[0])
    comp, n = ndi.label(m)
    sizes = ndi.sum(m, comp, range(1, n + 1)) * step ** 2
    m = np.isin(comp, 1 + np.nonzero(sizes >= min_area_um2)[0])
    if not m.any():
        raise ValueError("nothing of that size: choose other classes, a larger box, or a smaller min_area_um2")
    def ring(mask, oy, ox):  # the longest closed boundary of a mask, in level-0 pixels
        lines = contourpy.contour_generator(z=np.pad(mask, 1).astype(np.float32)).lines(0.5)
        p = max((np.asarray(ln, float) for ln in lines), key=len) - 1 + 0.5 + [ox, oy]  # pixel centers
        p = _simplify(p, 0.8)
        return [[round(float(x * scale), 1), round(float(y * scale), 1)] for x, y in p[:-1]]
    pieces = []
    comp, n = ndi.label(m)
    for i, sl in enumerate(ndi.find_objects(comp), 1):
        sl = (slice(max(0, sl[0].start - 1), sl[0].stop + 1), slice(max(0, sl[1].start - 1), sl[1].stop + 1))
        c = comp[sl] == i
        filled = ndi.binary_fill_holes(c)
        if filled.sum() * step ** 2 < min_area_um2:
            continue
        outer = ring(filled, sl[0].start, sl[1].start)
        if len(outer) < 3:
            continue
        hl, nh = ndi.label(filled & ~c)
        holes = [ring(hl == j, sl[0].start, sl[1].start) for j in range(1, nh + 1)]
        pieces.append({"outer": outer, "holes": [h for h in holes if len(h) >= 3]})
    return pieces


def tissue_priors(organ=None):
    """The organ's section of the priors document (with the shared notes of its group), or the list of organs."""
    text = (config.PACKAGE / "data" / "priors" / "tissue_priors.md").read_text(encoding="utf-8")
    parts = re.split(r"(?m)^## ", text)
    head, secs = parts[0], {p.split("\n", 1)[0].strip(): "## " + p for p in parts[1:]}
    organs = [k for k in secs if not k.lower().startswith(("how to", "using", "contents")) and "shared notes" not in k]
    if not organ:
        return {"organs": organs, "how_to_read": secs.get("How to read a rule", ""),
                "care": secs.get("Using priors with care", "")}
    q = organ.lower().strip()
    hit = next((k for k in organs if k.lower() == q), None) or next((k for k in organs if q in k.lower()), None) \
        or next((k for k in organs if any(w in k.lower() for w in q.split())), None)
    if hit is None:
        return {"error": f"no priors for {organ!r}", "organs": organs}
    shared = [k for k in secs if "shared notes" in k and any(w.lower() in k.lower() for w in re.split(r"[ ,()]+", hit) if len(w) > 3)]
    if hit in ("Colon", "Small intestine", "Stomach"):
        shared = [k for k in secs if k.startswith("Gastrointestinal")]
    return {"organ": hit, "text": "\n\n".join([secs[k] for k in shared] + [secs[hit]])}


# ---------------- figures ----------------
def _grid(ax, scale, box, n=8):
    """Ticks in level-0 pixels over a thumbnail of box."""
    x0, y0, x1, y1 = box
    step = 10 ** np.floor(np.log10(max(x1 - x0, y1 - y0) / n))
    step *= 5 if max(x1 - x0, y1 - y0) / step > 4 * n else 2 if max(x1 - x0, y1 - y0) / step > 2 * n else 1
    xt = np.arange(np.ceil(x0 / step) * step, x1, step)
    yt = np.arange(np.ceil(y0 / step) * step, y1, step)
    ax.set_xticks((xt - x0) / scale); ax.set_xticklabels([f"{int(v)}" for v in xt], fontsize=6, rotation=90)
    ax.set_yticks((yt - y0) / scale); ax.set_yticklabels([f"{int(v)}" for v in yt], fontsize=6)
    ax.grid(True, color="k", alpha=0.25, lw=0.4)


def _out(sample_id, s, name):
    d = (s.dir / "figures") if s is not None else _dir(sample_id)
    d.mkdir(exist_ok=True)
    return d / f"{name}_{time.strftime('%H%M%S')}.jpg"


def fig_classes(sample_id, s=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    meta, lab = load(sample_id)
    rgb, scale, box = thumbnail(sample_id)
    lab = lab[:rgb.shape[0], :rgb.shape[1]]
    asp = rgb.shape[0] / rgb.shape[1]
    fig, axs = plt.subplots(1, 2, figsize=(14, max(4, min(14, 7 * asp))))
    axs[0].imshow(rgb); axs[0].set_title("H&E (ticks: full-resolution pixels)", fontsize=9); _grid(axs[0], scale, box)
    show = np.ma.masked_less(lab, 0)
    axs[1].imshow(np.full_like(rgb, 255)); axs[1].imshow(show, cmap=ListedColormap(PALETTE[:meta["n_classes"]]),
                                                      vmin=0, vmax=meta["n_classes"] - 1, interpolation="nearest")
    _grid(axs[1], scale, box)
    for c in meta["classes"]:
        axs[1].plot([], [], "s", color=c["color"], label=f"{c['class']}: H {c['hematoxylin']:.2f} E {c['eosin']:.2f} "
                    f"tex {c['texture']:.2f} · {c['area_mm2']:.1f} mm² · {c['median_um_from_edge']:.0f} µm from edge")
    axs[1].legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    axs[1].set_title("classes of the tissue (1 = most hematoxylin)", fontsize=9)
    out = _out(sample_id, s, "he_classes")
    fig.savefig(out, dpi=110, bbox_inches="tight"); plt.close(fig)
    return str(out)


def fig_outline(sample_id, pieces, name, s=None, box=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rgb, scale, b = thumbnail(sample_id, box=box)
    fig, ax = plt.subplots(figsize=(9, max(4, min(12, 9 * rgb.shape[0] / rgb.shape[1]))))
    ax.imshow(rgb); _grid(ax, scale, b)
    for pc in pieces:
        for j, p in enumerate([pc["outer"]] + list(pc.get("holes") or [])):
            q = (np.asarray(p + p[:1], float) - [b[0], b[1]]) / scale
            ax.plot(q[:, 0], q[:, 1], color="#00e6ff", lw=1.6, ls="-" if j == 0 else "--")
            ax.plot(q[:, 0], q[:, 1], color="k", lw=0.5)
    nh = sum(len(pc.get("holes") or []) for pc in pieces)
    ax.set_title(f"proposed {name}: {len(pieces)} outline(s)" + (f", {nh} hole(s) dashed" if nh else ""), fontsize=9)
    out = _out(sample_id, s, f"he_outline_{re.sub(r'[^A-Za-z0-9]+', '_', name)}")
    fig.savefig(out, dpi=110, bbox_inches="tight"); plt.close(fig)
    return str(out)


def fig_view(sample_id, box_px=None, s=None, regions=None):
    """The image, or box_px of it, with the pixel grid and the regions drawn so far (outlines)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rgb, scale, b = thumbnail(sample_id, step_um=1.0 if box_px else STEP_UM, box=tuple(box_px) if box_px else None)
    fig, ax = plt.subplots(figsize=(10, max(4, min(13, 10 * rgb.shape[0] / rgb.shape[1]))))
    ax.imshow(rgb); _grid(ax, scale, b, n=10)
    for k, (nm, polys) in enumerate((regions or {}).items()):
        for p in polys:
            q = (np.asarray(p + p[:1], float) - [b[0], b[1]]) / scale
            ax.plot(q[:, 0], q[:, 1], color=PALETTE[k % len(PALETTE)], lw=1.4, label=nm if p is polys[0] else None)
    if regions:
        ax.legend(fontsize=7, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    ax.set_title("ticks: full-resolution pixels" + (f" · box {[int(v) for v in b]}" if box_px else ""), fontsize=9)
    out = _out(sample_id, s, "he_view")
    fig.savefig(out, dpi=110, bbox_inches="tight"); plt.close(fig)
    return str(out)
