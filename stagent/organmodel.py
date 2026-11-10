"""Organ models: the compartments of a tissue on the H&E (mucosa, muscularis mucosae, submucosa, ...) learned from
outlines a pathologist drew (GeoJSON, e.g. QuPath annotations in the image's full-resolution pixels), for the agent
to outline them on new slides of the organ (propose_he_region model_class=...; the user reviews every outline).

Features (FEATURE_VERSION): the image read at UM um per pixel in tiles; for every pixel hematoxylin and eosin (color
deconvolution, Ruifrok & Johnston 2001) smoothed at 2-32 um, the share of pale pixels (lumens, goblet mucin, fat) at
4-32 um, the texture of the hematoxylin (nuclei) and the coherence of the structure tensor of the eosin at 6-24 um
(oriented fibers: smooth muscle). A thin band such as the muscularis mucosae differs from the mucosa by these at a
few um, which a coarse thumbnail blurs away.

Training: the pixels inside each class's outlines (up to PER_CLASS_SLIDE per class and slide), a gradient-boosted
tree classifier. Its accuracy is measured on data it did not see: each slide held out in turn (two slides or more),
else one half of the annotated area. A model is one file per organ in the data folder's organ_models/ (all organs
together, kept across app updates): meta.json, the training pixels (train.npz: so a model is fitted again when the
app's scikit-learn changes, and new slides can be added to it) and the fitted classifier.

Prediction: class probabilities on every second pixel, then on every pixel where the two most likely classes are
close (UNSURE) (the same classes as everywhere, in a fraction of the time), smoothed over SMOOTH_UM, the most likely
class (none where the tissue is not stained). Features of tiles are computed in WORKERS threads. outline(): one class's pixels, small gaps closed, small holes filled, small pieces dropped, traced.
"""
import io, json, os, re, time, zipfile
import numpy as np
from scipy import ndimage as ndi
from . import config, render

FEATURE_VERSION = 1
UM = 2.0                # um per pixel of the features
TILE = 1024             # pixels of a tile (at UM)
PAD = 64                # context read around a tile (pixels at UM)
STAINED_OD = 0.06       # summed optical density of tissue (glass below)
PER_CLASS_SLIDE = 25000 # training pixels per class and slide, at most
PER_CLASS = 100000      # training pixels per class in a model, at most
SMOOTH_UM = 4.0         # the class probabilities are smoothed over this before the most likely class is taken
MIN_AREA_UM2 = 5000.0   # outline(): smaller pieces are dropped
MAX_HOLE_UM2 = 5000.0   # outline(): smaller holes are filled
CLOSE_UM = 6.0          # outline(): gaps narrower than this are closed
SUFFIX = ".stmodel"
UNSURE = 0.35           # predict(): pixels whose two most likely classes differ by less than this are predicted one by one
WORKERS = max(1, min(6, (os.cpu_count() or 2) - 1))  # threads computing the features of tiles


def models_dir():
    d = config.HOME / "organ_models"
    d.mkdir(parents=True, exist_ok=True)
    return d


def slug(organ):
    return re.sub(r"[^a-z0-9]+", "_", str(organ).lower()).strip("_") or "organ"


# ---------------- features ----------------
def _deconvolve(rgb):
    from .hetissue import _deconvolve as d
    return d(rgb)


def features(rgb, um=UM):
    """Per-pixel features (h, w, n) of an RGB tile read at um um per pixel, and its summed optical density."""
    od, H, E = _deconvolve(rgb.astype(np.float32))
    px = lambda u: max(0.5, u / um)
    g = lambda a, u: ndi.gaussian_filter(a, px(u))
    F = []
    for u in (2, 4, 8, 16, 32):
        F += [g(H, u), g(E, u)]
    white = (od < 0.25).astype(np.float32)
    for u in (4, 8, 16, 32):
        F.append(g(white, u))
    for u in (4, 8):
        m = g(H, u)
        F.append(np.sqrt(np.maximum(g(H * H, u) - m * m, 0)))
    I = g(E, 1.5)
    Ix, Iy = ndi.sobel(I, 1), ndi.sobel(I, 0)
    for u in (6, 12, 24):
        a, b, c = g(Ix * Ix, u), g(Ix * Iy, u), g(Iy * Iy, u)
        F.append(np.sqrt((a - c) ** 2 + 4 * b * b) / (a + c + 1e-6))
    return np.stack(F, -1).astype(np.float32), od


def image_um_per_px(path, given=None):
    if given:
        return float(given)
    from .align import image_um_per_px as f
    u = f(str(path))
    if not u:
        raise ValueError(f"the pixel size of {os.path.basename(str(path))} is not in the file: give it (um per pixel)")
    return float(u)


def _tiles(path, upp, want=None, say=None):
    """(ty, tx, h, w, rgb) tiles of the image read at UM um per pixel (pixels of that grid); want(ty, tx, h, w) picks."""
    sl = render.slide(str(path))
    s = UM / upp  # image pixels per grid pixel
    W, H = int(sl.W0 / s) + 1, int(sl.H0 / s) + 1
    todo = [(ty, tx, min(TILE, H - ty), min(TILE, W - tx)) for ty in range(0, H, TILE) for tx in range(0, W, TILE)]
    todo = [t for t in todo if want is None or want(*t)]
    for n, (ty, tx, h, w) in enumerate(todo):
        X0, Y0 = (tx - PAD) * s, (ty - PAD) * s
        X1, Y1 = (tx + w + PAD) * s, (ty + h + PAD) * s
        rgb = np.full((h + 2 * PAD, w + 2 * PAD, 3), 255, np.uint8)
        a0, b0, a1, b1 = max(X0, 0), max(Y0, 0), min(X1, sl.W0), min(Y1, sl.H0)
        if a1 > a0 and b1 > b0:
            r = np.asarray(sl.read(a0, b0, a1, b1, max(1, int(round((a1 - a0) / s))), max(1, int(round((b1 - b0) / s)))))[..., :3]
            oy, ox = int(round((b0 - Y0) / s)), int(round((a0 - X0) / s))
            r = r[: h + 2 * PAD - oy, : w + 2 * PAD - ox]
            rgb[oy:oy + r.shape[0], ox:ox + r.shape[1]] = r
        if say:
            say(n / max(len(todo), 1))
        yield ty, tx, h, w, rgb
    return


def _feature_tiles(path, upp, want=None, say=None):
    """(ty, tx, h, w, F, od) of the tiles (the tile itself, its context cut off): the image is read in order, the
    features computed in WORKERS threads (the filters release the GIL), a few tiles ahead."""
    from concurrent.futures import ThreadPoolExecutor
    from collections import deque
    from . import jobs

    def work(ty, tx, h, w, rgb):
        F, od = features(rgb)
        return ty, tx, h, w, F[PAD:PAD + h, PAD:PAD + w], od[PAD:PAD + h, PAD:PAD + w]
    with ThreadPoolExecutor(WORKERS) as ex:
        q = deque()
        for t in _tiles(path, upp, want, say):
            q.append(ex.submit(work, *t))
            if len(q) > 2 * WORKERS:
                jobs.checkpoint()
                yield q.popleft().result()
        while q:
            jobs.checkpoint()
            yield q.popleft().result()


def _grid(path, upp):
    sl = render.slide(str(path))
    s = UM / upp
    return s, int(sl.H0 / s) + 1, int(sl.W0 / s) + 1


def _raster(pieces, ty, tx, h, w, s):
    """Which pixels of a tile (grid at s image pixels per pixel) lie in the pieces (outline minus holes)."""
    from .regions import points_in
    m = np.zeros(h * w, bool)
    x0, y0, x1, y1 = tx * s, ty * s, (tx + w) * s, (ty + h) * s
    gx = gy = None
    for pc in pieces:
        o = np.asarray(pc["outer"], float)
        if o[:, 0].max() < x0 or o[:, 0].min() > x1 or o[:, 1].max() < y0 or o[:, 1].min() > y1:
            continue
        if gx is None:
            gx, gy = np.meshgrid((tx + np.arange(w) + 0.5) * s, (ty + np.arange(h) + 0.5) * s)
            gx, gy = gx.ravel(), gy.ravel()
        inn = points_in(gx, gy, o)
        for hl in pc.get("holes") or []:
            inn &= ~points_in(gx, gy, np.asarray(hl, float))
        m |= inn
    return m.reshape(h, w)


# ---------------- training data ----------------
def inspect_pair(image, geojson, um_per_px=None):
    """What a training pair holds: the GeoJSON's classes (outlines, area in mm2) and whether they lie on the image."""
    from .regions import parse_geojson
    sl = render.slide(str(image))
    try:
        upp = image_um_per_px(image, um_per_px)
    except ValueError:
        upp = None
    pieces, skipped = parse_geojson(open(geojson, encoding="utf-8").read())
    cls = {}
    allp = []
    for p in pieces:
        o = np.asarray(p["outer"], float)
        allp.append(o)
        a = abs(0.5 * float(np.dot(o[:, 0], np.roll(o[:, 1], -1)) - np.dot(o[:, 1], np.roll(o[:, 0], -1))))
        c = cls.setdefault(p["name"], {"outlines": 0, "area_px": 0.0})
        c["outlines"] += 1
        c["area_px"] += a
    P = np.concatenate(allp)
    inside = float(((P[:, 0] >= 0) & (P[:, 0] <= sl.W0) & (P[:, 1] >= 0) & (P[:, 1] <= sl.H0)).mean())
    return {"image": str(image), "geojson": str(geojson), "image_px": [sl.W0, sl.H0], "um_per_px": upp,
            "classes": [{"name": k, "outlines": v["outlines"],
                         "area_mm2": round(v["area_px"] * upp * upp / 1e6, 3) if upp else None} for k, v in cls.items()],
            "points_on_image": round(inside, 3), "skipped": skipped}


def training_pixels(image, geojson, classes, um_per_px=None, slide_id=0, progress=None, seed=0):
    """Features of pixels inside the outlines of the chosen classes ({GeoJSON class: model class}): X (float16), y
    (index into the model's classes, given by order of first appearance in `classes` values), the slide id and the
    pixels' place (grid rows, cols: for the held-out half)."""
    from .regions import parse_geojson
    upp = image_um_per_px(image, um_per_px)
    names = list(dict.fromkeys(classes.values()))
    pieces, _ = parse_geojson(open(geojson, encoding="utf-8").read())
    by = {}
    for p in pieces:
        if p["name"] in classes:
            by.setdefault(names.index(classes[p["name"]]), []).append(p)
    if not by:
        raise ValueError(f"{os.path.basename(str(geojson))}: none of the chosen classes is in it")
    s, H, W = _grid(image, upp)
    boxes = {k: np.array([[np.min(np.asarray(p["outer"])[:, 0]), np.min(np.asarray(p["outer"])[:, 1]),
                           np.max(np.asarray(p["outer"])[:, 0]), np.max(np.asarray(p["outer"])[:, 1])] for p in v])
             for k, v in by.items()}

    def want(ty, tx, h, w):
        x0, y0, x1, y1 = tx * s, ty * s, (tx + w) * s, (ty + h) * s
        return any(((b[:, 0] <= x1) & (b[:, 2] >= x0) & (b[:, 1] <= y1) & (b[:, 3] >= y0)).any() for b in boxes.values())
    rng = np.random.default_rng(seed)
    parts = {k: [] for k in by}
    for ty, tx, h, w, F, od in _feature_tiles(image, upp, want, progress):
        for k, pcs in by.items():
            m = _raster(pcs, ty, tx, h, w, s) & (od > STAINED_OD)
            idx = np.flatnonzero(m)
            if len(idx):
                idx = rng.choice(idx, min(len(idx), 4000), replace=False)  # (a tile's share; capped per slide below)
                r, c = np.divmod(idx, w)
                parts[k].append((F.reshape(-1, F.shape[-1])[idx].astype(np.float16), r + ty, c + tx))
    X, y, rc = [], [], []
    for k, ps in parts.items():
        if not ps:
            continue
        Fk = np.concatenate([p[0] for p in ps])
        rck = np.c_[np.concatenate([p[1] for p in ps]), np.concatenate([p[2] for p in ps])]
        keep = rng.choice(len(Fk), min(len(Fk), PER_CLASS_SLIDE), replace=False)
        X.append(Fk[keep]); y.append(np.full(len(keep), k, np.int16)); rc.append(rck[keep])
    if not X:
        raise ValueError(f"{os.path.basename(str(image))}: no stained pixels inside the outlines (is the GeoJSON in this "
                         "image's full-resolution pixels?)")
    X, y, rc = np.concatenate(X), np.concatenate(y), np.concatenate(rc)
    return {"X": X, "y": y, "slide": np.full(len(y), slide_id, np.int16), "rc": rc.astype(np.int32), "names": names}


# ---------------- the classifier ----------------
def _fit(X, y, seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    # 60 rounds of 15-leaf trees: as accurate on held-out slides as 150 rounds of 31 leaves, predicts 2.5 times faster
    return HistGradientBoostingClassifier(max_iter=60, max_leaf_nodes=15, learning_rate=0.1, early_stopping=False,
                                          random_state=seed).fit(X.astype(np.float32), y)


def _balance(y, cap, seed=0):
    rng = np.random.default_rng(seed)
    return np.concatenate([rng.choice(np.flatnonzero(y == k), min(int((y == k).sum()), cap), replace=False)
                           for k in np.unique(y)])


def held_out(X, y, slide, rc, names, say=None):
    """Accuracy on pixels the classifier did not see: each slide held out in turn (2+ slides), else the two halves
    of the annotated area along its longer side. Per class: recall (its pixels found) and precision (pixels called
    it that are it), and the classes it is taken for."""
    slides = np.unique(slide)
    if len(slides) >= 2:
        folds = [(slide != k, slide == k) for k in slides]
        how = f"each of the {len(slides)} slides held out in turn"
    else:
        ax = 0 if np.ptp(rc[:, 0]) >= np.ptp(rc[:, 1]) else 1
        mid = np.median(rc[:, ax])
        folds = [(rc[:, ax] < mid, rc[:, ax] >= mid), (rc[:, ax] >= mid, rc[:, ax] < mid)]
        how = "one half of the annotated area held out, then the other (one slide only: add slides for a fairer test)"
    pred, true = [], []
    for n, (tr, te) in enumerate(folds):
        if say:
            say(n / len(folds))
        if len(np.unique(y[tr])) < 2 or not te.any():
            continue
        keep = np.flatnonzero(tr)[_balance(y[tr], 40000)]
        clf = _fit(X[keep], y[keep])
        pred.append(clf.classes_[clf.predict_proba(X[te].astype(np.float32)).argmax(1)]); true.append(y[te])
    if not pred:
        return {"how": how, "classes": []}
    p, t = np.concatenate(pred), np.concatenate(true)
    out = []
    for k, nm in enumerate(names):
        if not (t == k).any():
            continue
        conf = {names[j]: round(float(np.mean(p[t == k] == j)), 3) for j in range(len(names)) if j != k and np.mean(p[t == k] == j) >= 0.03}
        out.append({"class": nm, "recall": round(float(np.mean(p[t == k] == k)), 3),
                    "precision": round(float(np.mean(t[p == k] == k)), 3) if (p == k).any() else None,
                    "taken_for": dict(sorted(conf.items(), key=lambda x: -x[1]))})
    return {"how": how, "pixel_accuracy": round(float(np.mean(p == t)), 3), "classes": out}


# ---------------- the model files ----------------
def _path(organ):
    return models_dir() / f"{slug(organ)}{SUFFIX}"


def _write(path, meta, data, clf):
    import pickle, sklearn
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    buf = io.BytesIO()
    np.savez_compressed(buf, **data)
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("meta.json", json.dumps(meta | {"sklearn": sklearn.__version__}, indent=1))
        z.writestr("train.npz", buf.getvalue())
        z.writestr("classifier.pkl", pickle.dumps(clf))
    os.replace(tmp, path)


def _read(path, part):
    with zipfile.ZipFile(path) as z:
        if part == "meta":
            return json.loads(z.read("meta.json"))
        if part == "data":
            d = np.load(io.BytesIO(z.read("train.npz")))
            return {k: d[k] for k in d.files}
        return z.read("classifier.pkl")


def _check(path):
    """A model file the app can use: its parts, and features of this app's version."""
    try:
        with zipfile.ZipFile(path) as z:
            names = set(z.namelist())
            meta = json.loads(z.read("meta.json"))
    except (OSError, zipfile.BadZipFile, KeyError, ValueError) as e:
        raise ValueError(f"not an organ model file ({e})")
    if not {"meta.json", "train.npz"} <= names or not meta.get("organ") or not meta.get("classes"):
        raise ValueError("not an organ model file (meta.json and train.npz needed)")
    if meta.get("feature_version") != FEATURE_VERSION:
        raise ValueError(f"made with features of version {meta.get('feature_version')}; this app uses {FEATURE_VERSION}")
    return meta


def list_models():
    out = []
    for f in sorted(models_dir().glob(f"*{SUFFIX}")):
        try:
            m = _check(f)
        except ValueError:
            continue
        out.append({k: m.get(k) for k in ("organ", "classes", "slides", "n_pixels", "trained", "held_out", "app_version")}
                   | {"file": str(f), "size_mb": round(f.stat().st_size / 1e6, 1)})
    return out


def model_of(organ):
    """The meta of the organ's model, or None."""
    f = _path(organ)
    if not f.exists():
        return None
    try:
        return _check(f)
    except ValueError:
        return None


_loaded = {}  # path -> (mtime, classifier, meta)


def load(organ):
    """The organ's classifier (fitted again from its training pixels when this app's scikit-learn cannot read it)."""
    import pickle, sklearn
    f = _path(organ)
    if not f.exists():
        raise ValueError(f"no model for {organ}: train one (Outline with Agent > Train a model) or load a model file")
    mt = f.stat().st_mtime_ns
    if f in _loaded and _loaded[f][0] == mt:
        return _loaded[f][1], _loaded[f][2]
    meta = _check(f)
    clf = None
    if meta.get("sklearn") == sklearn.__version__:
        try:
            clf = pickle.loads(_read(f, "clf"))
        except Exception:
            clf = None
    if clf is None:  # another scikit-learn: fitted again from the stored pixels (the same data, the same seed)
        d = _read(f, "data")
        keep = _balance(d["y"], PER_CLASS)
        clf = _fit(d["X"][keep], d["y"][keep])
        _write(f, meta, d, clf)
        mt = f.stat().st_mtime_ns
    _loaded[f] = (mt, clf, meta)
    return clf, meta


def train(organ, pairs, classes, add=False, progress=None, app_version=None):
    """Train (or extend, add=True) the organ's model from pairs [{"image", "geojson", "um_per_px"?}] and the classes
    to learn {GeoJSON class: model class} (several GeoJSON classes may become one). Returns the meta with the held-out
    accuracy. progress(f, message)."""
    say = progress or (lambda f, m: None)
    if not pairs:
        raise ValueError("add at least one H&E with its GeoJSON")
    if not classes:
        raise ValueError("choose the classes to learn")
    f = _path(organ)
    old = _read(f, "data") if add and f.exists() else None
    old_meta = _check(f) if old is not None else None
    names = list(old_meta["classes"]) if old_meta else []
    for c in classes.values():
        if c not in names:
            names.append(c)
    first = int(old["slide"].max()) + 1 if old is not None else 0
    parts = []
    for i, p in enumerate(pairs):
        lo, hi = 0.8 * i / len(pairs), 0.8 * (i + 1) / len(pairs)
        d = training_pixels(p["image"], p["geojson"], classes, p.get("um_per_px"), slide_id=first + i,
                            progress=lambda x, i=i, lo=lo, hi=hi: say(lo + (hi - lo) * x, f"slide {i + 1} of {len(pairs)}: features"))
        d["y"] = np.array([names.index(d["names"][k]) for k in d["y"]], np.int16)
        parts.append(d)
    X = np.concatenate(([old["X"]] if old is not None else []) + [d["X"] for d in parts])
    y = np.concatenate(([old["y"]] if old is not None else []) + [d["y"] for d in parts])
    sl = np.concatenate(([old["slide"]] if old is not None else []) + [d["slide"] for d in parts])
    rc = np.concatenate(([old["rc"]] if old is not None else []) + [d["rc"] for d in parts])
    if len(np.unique(y)) < 2:
        raise ValueError("a model needs at least two classes")
    say(0.82, "accuracy on held-out data")
    ho = held_out(X, y, sl, rc, names, lambda x: say(0.82 + 0.13 * x, "accuracy on held-out data"))
    say(0.96, "fitting the model")
    keep = _balance(y, PER_CLASS)
    clf = _fit(X[keep], y[keep])
    slides = (old_meta["slides"] if old_meta else []) + [
        {"image": os.path.basename(str(p["image"])), "geojson": os.path.basename(str(p["geojson"]))} for p in pairs]
    meta = {"organ": organ, "classes": names, "feature_version": FEATURE_VERSION, "um": UM, "slides": slides,
            "class_map": {**(old_meta or {}).get("class_map", {}), **classes}, "n_pixels": int(len(y)),
            "pixels_per_class": {names[k]: int((y == k).sum()) for k in np.unique(y)},
            "trained": time.strftime("%Y-%m-%dT%H:%M:%S"), "held_out": ho, "app_version": app_version}
    _write(f, meta, {"X": X, "y": y, "slide": sl, "rc": rc}, clf)
    _loaded.pop(f, None)
    say(1.0, "done")
    return meta


def import_model(path):
    """Put a model file (e.g. from a colleague) among the organ models; it replaces the organ's model."""
    meta = _check(path)
    dst = _path(meta["organ"])
    tmp = dst.with_name(f".{dst.name}.{os.getpid()}.tmp")
    with open(path, "rb") as a, open(tmp, "wb") as b:
        b.write(a.read())
    os.replace(tmp, dst)
    _loaded.pop(dst, None)
    return meta


def delete_model(organ):
    f = _path(organ)
    if not f.exists():
        raise ValueError(f"no model for {organ}")
    trash = config.TRASH / "organ_models"
    trash.mkdir(parents=True, exist_ok=True)
    os.replace(f, trash / f"{f.stem}-{time.strftime('%Y%m%d-%H%M%S')}{SUFFIX}")
    _loaded.pop(f, None)


# ---------------- prediction on a sample's image ----------------
def _cache(sample_id, organ):
    meta = _check(_path(organ))
    d = config.SAMPLES / sample_id / "organ_model"
    d.mkdir(exist_ok=True)
    key = f"{slug(organ)}-{meta['trained'].replace(':', '')}"
    return d / f"{key}.npz"


def predict(sample_id, organ, progress=None):
    """The most likely class of every stained pixel of the sample's H&E (grid at UM; -1: glass), the classes and the
    grid's image pixels per pixel. Kept in the sample's folder until the model changes."""
    from .hetissue import _image
    path, upp, kind, known = _image(sample_id)
    if kind != "H&E":
        raise ValueError("the image is not an H&E")
    f = _cache(sample_id, organ)
    busy = f.with_suffix(".computing")  # the app may have started it (Outline with Agent) for the agent's call
    end = time.time() + 900
    while busy.exists() and time.time() - busy.stat().st_mtime < 900 and time.time() < end and not f.exists():
        time.sleep(1)
    if f.exists():
        d = np.load(f)
        return d["lab"], json.loads(str(d["classes"])), float(d["s"])
    busy.write_text(str(os.getpid()))
    try:
        return _predict(sample_id, organ, f, path, upp, progress)
    finally:
        busy.unlink(missing_ok=True)


def _predict(sample_id, organ, f, path, upp, progress):
    clf, meta = load(organ)
    s, H, W = _grid(path, upp)
    lab = np.full((H, W), -1, np.int8)
    say = progress or (lambda x, m: None)
    from .hetissue import load as seg_load, segmented
    tissue = None
    if segmented(sample_id):  # tiles without tissue on the thumbnail are skipped
        m, l = seg_load(sample_id)
        tissue, ts = ndi.binary_dilation(l >= 0, iterations=3), m["scale_px"]
    def want(ty, tx, h, w):
        if tissue is None:
            return True
        f_ = ts / s
        return tissue[int(ty / f_):int(np.ceil((ty + h) / f_)) + 1, int(tx / f_):int(np.ceil((tx + w) / f_)) + 1].any()
    sig = SMOOTH_UM / UM
    for ty, tx, h, w, F, od in _feature_tiles(path, upp, want, lambda x: say(x, "the model reads the H&E")):
        stained = od > STAINED_OD
        if not stained.any():
            continue
        lab[ty:ty + h, tx:tx + w] = _classify(clf, F, stained, sig)
        busy = f.with_suffix(".computing")
        if busy.exists():
            os.utime(busy)  # (still at it)
    tmp = f.with_name(f".{f.stem}.{os.getpid()}.npz")
    np.savez_compressed(tmp, lab=lab, classes=json.dumps(meta["classes"]), s=s)
    os.replace(tmp, f)
    for old in f.parent.glob(f"{slug(organ)}-*.npz"):
        if old != f:
            old.unlink(missing_ok=True)
    return lab, meta["classes"], s


def _classify(clf, F, stained, sig):
    """The most likely class of each stained pixel of a tile (-1 elsewhere): probabilities on every second pixel,
    interpolated, then on each pixel where the two most likely classes are close; smoothed over sig pixels."""
    h, w, k = F.shape[0], F.shape[1], len(clf.classes_)
    c = np.zeros(((h + 1) // 2, (w + 1) // 2, k), np.float32)
    sc = stained[::2, ::2]
    if sc.any():
        c[sc] = clf.predict_proba(F[::2, ::2][sc])
    P = ndi.zoom(c, (h / c.shape[0], w / c.shape[1], 1), order=1, grid_mode=True, mode="nearest")[:h, :w]
    top = np.sort(P, -1)
    unsure = ndi.binary_dilation((top[..., -1] - top[..., -2] < UNSURE) & stained, iterations=2) & stained
    if unsure.any():
        P[unsure] = clf.predict_proba(F[unsure])
    P = ndi.gaussian_filter(P, (sig, sig, 0))
    return np.where(stained, clf.classes_[P.argmax(-1)], -1).astype(np.int8)


def outline(sample_id, organ, class_name, box_px=None, within_px=None, min_area_um2=MIN_AREA_UM2, progress=None):
    """Outlines of one class of the organ's model on the sample's H&E: [{"outer", "holes"}] in level-0 pixels."""
    from .regions import _trace_mask
    lab, names, s = predict(sample_id, organ, progress)
    hit = [i for i, n in enumerate(names) if n.lower() == str(class_name).lower()]
    if not hit:
        raise ValueError(f"the {organ} model has no class {class_name!r}: {names}")
    m = lab == hit[0]
    if box_px is not None:
        x0, y0, x1, y1 = [float(v) / s for v in box_px]
        keep = np.zeros_like(m)
        keep[max(0, int(y0)):max(0, int(np.ceil(y1))), max(0, int(x0)):max(0, int(np.ceil(x1)))] = True
        m &= keep
    if within_px is not None:
        from .regions import points_in
        ys, xs = np.nonzero(m)
        inn = points_in((xs + 0.5) * s, (ys + 0.5) * s, np.asarray(within_px, float))
        m[ys[~inn], xs[~inn]] = False
    a = UM * UM
    m = ndi.binary_opening(m, iterations=1)
    r = max(1, int(round(CLOSE_UM / UM / 2)))
    m = ndi.binary_closing(np.pad(m, r + 1), iterations=r)[r + 1:-r - 1, r + 1:-r - 1]
    holes = ndi.binary_fill_holes(m) & ~m
    hl, n = ndi.label(holes)
    if n:
        m |= np.isin(hl, 1 + np.nonzero(ndi.sum(holes, hl, range(1, n + 1)) * a < MAX_HOLE_UM2)[0])
    cl, n = ndi.label(m)
    if n:
        m = np.isin(cl, 1 + np.nonzero(ndi.sum(m, cl, range(1, n + 1)) * a >= float(min_area_um2))[0])
    if not m.any():
        raise ValueError(f"no {class_name} of that size found by the {organ} model here")
    return _trace_mask(m, s)
