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
  gaps of CLOSE_UM, holes under MAX_HOLE_UM2 filled, pieces under min_area dropped, traced smoothly into polygons
  (the mask blurred and its half level traced between pixels, corners rounded), in the image's
  full-resolution pixels.
view(): the image, or a box of it, with a grid labelled in full-resolution pixels, so an outline can be given point
  by point.
Results are kept in the sample's folder (he_segment/) for outline() and the review window.
"""
import json, re, time
import numpy as np
from . import config
from .osutil import dump_json

STEP_UM = 6.0          # thumbnail resolution
MAX_SIDE = 4500        # thumbnail pixels on its longer side, at most
TISSUE_OD = 0.05       # summed optical density above this (and colored: MIN_SAT) is tissue
MIN_SAT = 0.05         # saturation of stained tissue, at least (the gray of glass is near 0; pale submucosa ~0.1)
MIN_SPECK_UM2 = 20000.0  # tissue pieces smaller than this are dust
MIN_HOLE_UM2 = 3000.0    # holes in the tissue smaller than this are filled; larger ones (lumens) stay out
CLOSE_UM = 40.0        # outline(): gaps narrower than this are closed
MIN_AREA_UM2 = 20000.0 # outline(): smaller pieces are dropped
MAX_HOLE_UM2 = 50000.0 # outline(): holes smaller than this are filled
SMOOTH_PX = 1.2        # outline(): the mask is blurred by this many thumbnail pixels before it is traced
SIMPLIFY_PX = 0.15     # outline(): points closer than this (thumbnail pixels) to the line are dropped
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


def segment(sample_id, n_classes=6, progress=None, organ=None):
    """Classes of the tissue on the image (see the module). Returns the stats; the label image is kept. organ: the
    tissue whose priors guide the choice of classes (kept with them, and remembered for the sample)."""
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
            "image": str(path), "shape": list(lab.shape), "classes": classes, "pixel_size_known": known, "organ": organ}
    dump_json(meta, d / "meta.json", indent=1)
    if organ:
        _remember_organ(sample_id, organ)
    say(f"done in {time.time() - t0:.1f} s")
    return meta


def segmented(sample_id):
    return (config.SAMPLES / sample_id / "he_segment" / "meta.json").exists()


def sample_organ(sample_id):
    """The organ given for the sample's tissue earlier (segment_he), or None."""
    try:
        organ = json.load(open(config.SAMPLES / sample_id / "meta.json", encoding="utf-8")).get("organ")
    except (OSError, ValueError):
        return None
    return str(organ or "").strip() or None


# words that name an organ in a sample's name or description, by built-in marker list
ORGAN_WORDS = {"colon": ["colon", "colorectal", "rectum", "rectal", "large intestine", "sigmoid", "cecum", "caecum"],
               "small_intestine": ["small intestine", "small bowel", "ileum", "ileal", "jejunum", "jejunal", "duodenum",
                                   "duodenal"],
               "stomach": ["stomach", "gastric", "antrum", "fundus"], "liver": ["liver", "hepatic", "hepatocellular"],
               "pancreas": ["pancreas", "pancreatic"], "kidney": ["kidney", "renal", "nephron"],
               "lung": ["lung", "pulmonary", "alveolar", "bronch"], "heart": ["heart", "cardiac", "myocard"],
               "breast": ["breast", "mammary"], "skin": ["skin", "cutaneous", "epidermis", "dermis"],
               "brain": ["brain", "cerebr", "hippocamp", "neocortex"],   # not "cortex": renal / adrenal cortex
               "lymph_node": ["lymph node", "lymph-node", "tonsil"]}   # matched at the start of a word ("adrenal"
# is not "renal"); "ventric" is left out: the brain has ventricles too


GENERIC_TYPES = {"T cells", "Endothelial cells", "Stromal", "Proliferating cells", "Endocrine", "Exocrine", "Vascular EC",
                 "VSMCs"}   # names other lists may use too: no evidence of an organ


def _distinct_types():
    """Per built-in organ list, the cell types no other organ's list has (its organ's own evidence)."""
    from . import organs
    names = {o["organ"]: {n["name"] for n in organs.tree(o["organ"])} for o in organs.organs() if o["organ"] != "other"}
    return {o: v - set().union(*[w for k, w in names.items() if k != o]) - GENERIC_TYPES for o, v in names.items()}


def guess_organ(s):
    """The organ of the sample's tissue for its priors: {"organ" (a section of tissue_priors, or None), "from",
    "evidence"}. In order: the organ given earlier (segment_he), the built-in marker list the annotation's cell types
    come from, words in the sample's name or description."""
    evidence = []
    given = sample_organ(s.sample_id)
    if given and given != "other":
        evidence.append(f"given for this sample: {given}")
    names = set(getattr(s, "node_names", None) or []) | set(s.types)
    hits = {o: sorted(d & names) for o, d in _distinct_types().items() if d & names}
    by_list = None
    if hits:
        best = sorted(hits.items(), key=lambda kv: -len(kv[1]))
        if len(best) == 1 or len(best[0][1]) > len(best[1][1]):
            by_list = best[0][0]
            evidence.append(f"the marker list has cell types only the built-in {by_list} list has: {', '.join(best[0][1][:6])}")
        else:
            evidence.append("the marker list matches several organs: " + ", ".join(f"{k} ({', '.join(v[:3])})" for k, v in best[:3]))
    meta = s.meta
    text = " ".join(str(x) for x in [s.sample_id, meta.get("name"), meta.get("description"),
                                     *(meta.get("metadata") or {}).values()] if x).lower().replace("_", " ")
    by_text = [o for o, words in ORGAN_WORDS.items() if any(re.search(r"\b" + re.escape(w), text) for w in words)]
    if by_text:
        evidence.append("the sample's name or description: " + ", ".join(by_text))
    pick, src = ((given, "given") if given and given != "other" else (by_list, "marker list") if by_list else
                 (by_text[0], "name or description") if len(by_text) == 1 else (None, None))
    hit = tissue_priors(pick).get("organ") if pick else None
    if hit is None and src == "given":
        hit = given   # an organ the priors do not cover, chosen for the sample: kept (the user may add its priors)
    return {"organ": hit, "from": src if hit else None, "covered": bool(hit) and not default_priors(hit).get("error"),
            "evidence": evidence or ["nothing names the organ"]}


def _remember_organ(sample_id, organ):
    import os
    f = config.SAMPLES / sample_id / "meta.json"
    m = json.load(open(f))
    if m.get("organ") == organ:
        return
    m["organ"] = organ
    dump_json(m, f, indent=1)   # atomic, with a temporary name of its own


NEXT_HEADINGS = ("Priors", "References", "Depth rule", "Verification note", "The depth rule", "How to read the depth")


def compartments(organ):
    """The organ's "Compartments to draw" from the priors (how each compartment looks on the H&E), to match the classes
    of segment() to; organ "other": the organ-independent ones. {"organ", "compartments_to_draw"} or {"error", ...}."""
    q = str(organ or "").strip()
    if not q:
        return {"error": "which organ is this tissue? Name it", "organs": default_priors()["organs"],
                "hint": "use one of the organs, or organ='other' for a tissue none of them covers"}
    r = tissue_priors("General cross-organ priors" if q.lower() in ("other", "none", "unknown") else q)
    if r.get("error"):
        return r | {"hint": "use one of the organs, or organ='other' for a tissue none of them covers"}
    text = r["text"]   # plain text: "=== Organ ===" sections, a "Compartments to draw" line with its "- " items
    blocks = []
    parts = re.split(r"(?m)^=== (.+?) ===[ \t]*$", text)
    for title, body in zip(parts[1::2], parts[2::2]):
        lines = body.splitlines()
        i = next((k for k, l in enumerate(lines) if l.lower().startswith("compartments to draw")), None)
        if i is None:
            continue
        items = [lines[i][len("compartments to draw"):].lstrip(":").strip()]
        for line in lines[i + 1:]:
            if line.startswith(NEXT_HEADINGS):
                break   # the section's next part (its priors, references, depth notes)
            items.append(line)
        blocks.append((f"=== {title} ===\n" if blocks or "shared" in title.lower() else "") + "\n".join(items).strip())
    return {"organ": r["organ"], "compartments_to_draw": "\n\n".join(blocks) or text[:4000],
            "note": "match each class (hematoxylin, eosin, texture, place) to these descriptions; one compartment is "
                    "often 2-3 classes; a class that fits none is left out"}


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


def _chaikin(p, rounds):
    """Corner cutting of a closed ring (the last point = the first): every edge becomes its 1/4 and 3/4 points."""
    q = p[:-1] if np.allclose(p[0], p[-1]) else p
    for _ in range(rounds):
        nxt = np.roll(q, -1, axis=0)
        q = np.stack([0.75 * q + 0.25 * nxt, 0.25 * q + 0.75 * nxt], 1).reshape(-1, 2)
    return np.vstack([q, q[:1]])


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
    def ring(mask, oy, ox):  # the longest closed boundary of a mask, smooth, in level-0 pixels
        # the mask blurred, its 0.5 level traced between pixels (no pixel stairs), corners rounded (Chaikin), then
        # only points closer than SIMPLIFY_PX to the line dropped
        z = ndi.gaussian_filter(np.pad(mask, 4).astype(np.float32), SMOOTH_PX)
        lines = contourpy.contour_generator(z=z).lines(0.5)
        p = max((np.asarray(ln, float) for ln in lines), key=len) - 4 + 0.5 + [ox, oy]  # pixel centers
        p = _simplify(_chaikin(p, 2), SIMPLIFY_PX)
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


def _user_dir():
    return config.HOME / "priors"   # the user's changes to the priors: one file per organ (removed and added priors)


def to_txt(md):
    """The priors as plain text: "=== Title ===" headings, no ** or ` marks (rules, lists and references stay)."""
    out = []
    for line in md.splitlines():
        if line.strip() == "---":
            out.append("")
            continue
        m = re.match(r"^#{1,3} (.+)$", line)
        if m:
            line = f"=== {m.group(1).strip()} ==="
        out.append(line.replace("**", "").replace("`", "").rstrip())
    return re.sub(r"\n{3,}", "\n\n", "\n".join(out)).strip() + "\n"


WINDOWS_RESERVED = {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)), *(f"lpt{i}" for i in range(1, 10))}


def _slug(organ):
    """The organ's file name: its letters and digits in any script (a Chinese name keeps its own file), safe on
    Windows; a name without any: a short hash of it."""
    import hashlib
    name = str(organ).strip()
    slug = re.sub(r"[\W_]+", "_", name.lower()).strip("_")[:80]
    if not slug:
        return "organ_" + hashlib.sha1(name.encode("utf-8")).hexdigest()[:8]
    return slug + "_" if slug in WINDOWS_RESERVED else slug


def _user_file(organ):
    return _user_dir() / f"{_slug(organ)}.json"


def _changes(organ):
    f = _user_file(organ)
    try:
        c = json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}
    except ValueError:
        c = {}
    if not isinstance(c, dict) or (c.get("organ") and _slug(c["organ"]) != _slug(organ)):
        c = {}   # another organ's file (an older version gave every non-ASCII name the same one)
    return {"organ": c.get("organ") or organ, "removed": list(c.get("removed") or []), "added": list(c.get("added") or []),
            "next_id": int(c.get("next_id") or 0)}


def user_organs():
    """Organs the user added priors for that the app's document does not cover."""
    out = []
    for f in sorted(_user_dir().glob("*.json")) if _user_dir().exists() else []:
        try:
            c = json.loads(f.read_text(encoding="utf-8"))
        except ValueError:
            continue
        if not isinstance(c, dict) or not c.get("organ") or f.stem != _slug(c["organ"]):
            continue   # not the organ's own file
        if default_priors(c["organ"]).get("error") and c.get("added"):
            out.append(c["organ"])
    return out


PRIOR_START = re.compile(r"(?m)^- \*\*([A-Z]+-\d+[a-z]?)\. ")


def _prior_blocks(md):
    """The priors of an organ's section: [{"id", "title", "rules", "start", "end"}] (positions in md)."""
    sec0 = md.rfind("\n## ") + 1 if md.count("\n## ") else 0   # the organ's own section (after its group's notes)
    ref = md.find("**References**", sec0)
    stop = ref if ref > 0 else len(md)
    ms = [m for m in PRIOR_START.finditer(md, sec0) if m.start() < stop]
    out = []
    for k, m in enumerate(ms):
        end = ms[k + 1].start() if k + 1 < len(ms) else stop
        block = md[m.start():end]
        title = re.sub(r"\s+", " ", block[m.end() - m.start():].split("**", 1)[0]).strip()
        rp = re.search(r"(?ms)^  - Rule.*?(?=^  - |\Z)", block)   # the Rule item with its sub-items
        rules = [re.sub(r"\s+", " ", x).strip() for x in re.findall(r"`(kind=[^`]+)`", rp.group(0) if rp else "")]
        short = re.split(r"(?<=[a-z0-9)\]])\. (?=[A-Z])", title, maxsplit=1)[0].rstrip(".") + "."
        out.append({"id": m.group(1), "title": title, "short": short, "rules": rules, "start": m.start(), "end": end})
    return out


def organ_priors(organ):
    """What the Prior window shows: the organ's priors (id, statement, rules), those the user removed and added."""
    organ = str(organ or "").strip()
    if not organ:
        raise ValueError("name the organ")
    d = default_priors(organ)
    name = d.get("organ") or next((o for o in user_organs() if _slug(o) == _slug(organ)), str(organ).strip())
    ch = _changes(name)
    blocks = [] if d.get("error") else _prior_blocks(d["text"])
    return {"organ": name, "covered": not d.get("error"),
            "priors": [{k: b[k] for k in ("id", "title", "short", "rules")} for b in blocks],
            "removed": [x for x in ch["removed"] if any(b["id"] == x for b in blocks)], "added": ch["added"],
            "edited": bool(ch["removed"] or ch["added"])}


def tissue_priors(organ=None):
    """The organ's priors as plain text: the app's document (with its group's shared notes) less the priors the user
    removed, plus those the user added (Prior in the app; "edited": true); without organ, the list of organs."""
    organ = str(organ or "").strip()
    if not organ:
        r = default_priors()
        return {"organs": r["organs"] + user_organs(), "how_to_read": to_txt(r["how_to_read"]),
                "care": to_txt(r["care"])}
    d = default_priors(organ)
    v = organ_priors(organ)
    if d.get("error") and not v["added"]:
        return d
    md = "" if d.get("error") else d["text"]
    for b in sorted([b for b in _prior_blocks(md) if b["id"] in v["removed"]], key=lambda b: -b["start"]):
        md = md[:b["start"]] + md[b["end"]:]
    text = to_txt(md) if md else f"=== {v['organ']} ===\n"
    if v["removed"]:
        text += f"\nRemoved by the user (do not apply): {', '.join(v['removed'])}\n"
    if v["added"]:
        text += "\nPriors added by the user\n" + "".join(
            f"- {a['id']}. {a['text']}\n" + (f"  - Rule: {a['rule']}\n" if a.get("rule") else
                                              f"  - Rule: none ({a['no_rule']})\n" if a.get("no_rule") else
                                              "  - Rule: not written yet (set_prior_rule)\n") for a in v["added"])
    out = {"organ": v["organ"], "text": text, "edited": v["edited"]}
    if v["edited"]:
        out["note"] = "the user removed or added priors in the app: follow this version, not the app's default"
    return out


def set_priors(organ, removed=None, added=None):
    """Keep the user's changes to an organ's priors: the ids removed, and the added priors [{"text"}, or as returned
    with their "id" and rule]: a new one gets the next id U-n and no rule yet (the agent writes it, set_prior_rule).
    Nothing removed or added: back to the app's priors."""
    name = str(organ or "").strip()
    if not name or len(name) > 80:
        raise ValueError("name the organ (at most 80 characters)")
    d = default_priors(name)
    name = d.get("organ") or name
    ids = {b["id"] for b in _prior_blocks(d["text"])} if not d.get("error") else set()
    removed = sorted({str(x) for x in (removed or []) if str(x) in ids})
    ch = _changes(name)
    old = {a["id"]: a for a in ch["added"]}
    # the ids U-n are never given twice (a deleted prior's id may still be in a rule request to the agent)
    n = max([int(i[2:]) for i in old if re.fullmatch(r"U-\d+", i)] + [ch["next_id"] - 1, 0])
    clean = []
    for a in added or []:
        if isinstance(a, str):
            a = {"text": a}
        elif not isinstance(a, dict):
            continue
        text = " ".join(str(a.get("text") or "").split())
        if not text:
            continue
        if len(text) > 600:
            raise ValueError("a prior is longer than 600 characters")
        prev = old.get(a.get("id"))
        if prev and prev["text"] == text:
            clean.append(prev)   # unchanged: its rule stays
        else:
            n += 1
            clean.append({"id": f"U-{n}", "text": text})
    _write_changes(name, removed, clean, n + 1)
    return organ_priors(name)


def _write_changes(name, removed, added, next_id):
    f = _user_file(name)
    if not removed and not added:
        f.unlink(missing_ok=True)
    else:
        f.parent.mkdir(parents=True, exist_ok=True)
        dump_json({"organ": name, "removed": removed, "added": added, "next_id": next_id}, f, indent=1)


RULE_KINDS, RULE_MODES, RULE_STRENGTHS = ("region", "depth", "neighborhood"), ("only", "exclude", "enrich"), ("strong", "moderate", "weak")


def set_prior_rule(organ, prior_id, types=None, mode=None, strength="moderate", kind="region", regions=None,
                   max_depth_um=None, neighbor_types=None, min_frac=None, no_rule=None):
    """The rule of a prior the user added (in the constraint language), or no_rule: why it cannot be one."""
    v = organ_priors(organ)
    ch = _changes(v["organ"])
    a = next((x for x in ch["added"] if x["id"] == prior_id), None)
    if a is None:
        raise ValueError(f"no added prior {prior_id!r} for {v['organ']}: " + ", ".join(x["id"] for x in ch["added"]))
    a.pop("rule", None); a.pop("no_rule", None)
    if no_rule:
        a["no_rule"] = " ".join(str(no_rule).split())[:300]
    else:
        types = [str(x).strip() for x in (types or []) if str(x).strip()]
        if not types:
            raise ValueError("give the cell types (or no_rule with the reason)")
        if kind not in RULE_KINDS or mode not in RULE_MODES or strength not in RULE_STRENGTHS:
            raise ValueError(f"kind one of {RULE_KINDS}, mode one of {RULE_MODES}, strength one of {RULE_STRENGTHS}")
        parts = [f"kind={kind}", f"types=[{', '.join(types)}]"]
        if kind == "region":
            regions = [str(x).strip() for x in (regions or []) if str(x).strip()]
            if not regions:
                raise ValueError("a region rule needs regions (the names the user draws)")
            parts.append(f"regions=[{', '.join(regions)}]")
        elif kind == "depth":
            if max_depth_um is None:
                raise ValueError("a depth rule needs max_depth_um")
            parts.append(f"max_depth_um={float(max_depth_um):g}")
        else:
            if not neighbor_types or min_frac is None:
                raise ValueError("a neighborhood rule needs neighbor_types and min_frac")
            parts += [f"neighbor_types=[{', '.join(neighbor_types)}]", f"min_frac={float(min_frac):g}"]
        a["rule"] = ", ".join(parts + [f"mode={mode}", f"strength={strength}"])
    _write_changes(v["organ"], ch["removed"], ch["added"], ch["next_id"])
    return next(x for x in ch["added"] if x["id"] == prior_id)


# other names of the priors' sections, besides their own (kept short: a name not known is an organ not covered)
PRIOR_ALIASES = {"Colon": ["large intestine", "large bowel", "colorectal", "colorectum", "rectum"],
                 "Small intestine": ["small bowel", "ileum", "jejunum", "duodenum"], "Lung": ["lungs"],
                 "Lymph node and tonsil": ["lymph nodes", "tonsils"]}


def _key(name):
    """A name compared without case, spaces, "_", "-" or "`" (the built-in lists' names: small_intestine, lymph_node)."""
    return re.sub(r"[\s_`-]+", "", str(name).lower())


def _section_keys(k):
    """The names that mean section k: its title, with and without the part in brackets, each organ of an "X and Y"
    title ("Brain (cortex)": "brain"; "Lymph node and tonsil": "lymph node", "tonsil"), its aliases."""
    base = re.sub(r"\s*\(.*\)", "", k).strip()
    names = [k, base, k.replace("(", " ").replace(")", " "), *re.split(r",\s*|\s+and\s+", base),
             *PRIOR_ALIASES.get(k, [])]
    return {_key(n) for n in names if n.strip()}


def default_priors(organ=None):
    """The organ's section of the app's priors document (with the shared notes of its group), or the list of organs."""
    text = (config.PACKAGE / "data" / "priors" / "tissue_priors.md").read_text(encoding="utf-8")
    parts = re.split(r"(?m)^## ", text)
    head, secs = parts[0], {p.split("\n", 1)[0].strip(): "## " + p for p in parts[1:]}
    organs = [k for k in secs if not k.lower().startswith(("how to", "using", "contents")) and "shared notes" not in k]
    organ = str(organ or "").strip()
    if not organ:
        return {"organs": organs, "how_to_read": secs.get("How to read a rule", ""),
                "care": secs.get("Using priors with care", "")}
    q = _key(organ)   # a whole name only: "Large intestine" is not "Small intestine", "Ear" not "Heart"
    hit = next((k for k in organs if q in _section_keys(k)), None)
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
    ax.set_xlim(-0.5, rgb.shape[1] - 0.5); ax.set_ylim(rgb.shape[0] - 0.5, -0.5)  # a box: the outline beyond is cut
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
