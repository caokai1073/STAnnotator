"""Figures: region maps (PNG), H&E conflict crops and tissue overviews (JPEG, smaller for MCP clients). Each returns
its path. Slide reads any box of a pyramidal H&E at a given output size (also used by the viewer's tiles)."""
import contextlib, logging, math, re, threading
from pathlib import Path
import numpy as np, tifffile, zarr
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import config

FIXED_LINEAGE_RGB = {"epi": (235, 130, 20), "fib": (60, 165, 60), "stroma": (40, 160, 170), "lymph": (50, 100, 220),
                     "myeloid": (210, 60, 170), "Unassigned": (150, 150, 150), "Low quality": (190, 190, 190), "Not annotated": (150, 150, 150)}


# distinct colors that stay visible on H&E (no pale tints); assigned in marker-tree order
DISTINCT = ["#f58231", "#4363d8", "#3cb44b", "#f032e6", "#42d4f4", "#9a6324", "#911eb4", "#e6194b", "#469990",
            "#000075", "#bfef45", "#800000", "#808000", "#ffe119", "#a9a9a9"]


def lineage_palette(lineages):
    """RGB per lineage name: fixed for the reference lineages, a distinct qualitative list for a user's marker tree."""
    out, k = {}, 0
    for n in lineages:
        if n in FIXED_LINEAGE_RGB:
            out[n] = FIXED_LINEAGE_RGB[n]
        else:
            h = DISTINCT[k % len(DISTINCT)]
            out[n] = tuple(int(h[i:i + 2], 16) for i in (1, 3, 5))
            k += 1
    return out


REGION_COLORS = {"mucosa": "#e6862a", "submucosa": "#57a35a", "muscularis_propria": "#3d76b8",
                 "lymphoid_aggregate": "#8e5bb5", "fibromuscular_hyperplasia": "#c9b63c", "granuloma": "#b07aa1",
                 "unknown": "#cfcfcf"}
class IncompleteImage(ValueError):
    """The image file ends before the data it lists (a copy or sync still going, or a file cut short)."""


class _Problems(logging.Handler):
    def __init__(self):
        super().__init__(logging.WARNING)
        self.found = []

    def emit(self, record):
        self.found.append(record.getMessage())


@contextlib.contextmanager
def _tiff_problems():
    """What tifffile logs while reading a file's structure: it logs (does not raise) a tag whose value lies past the
    end of the file, then goes on without it (e.g. the pyramid's levels when they are stored last)."""
    h, log = _Problems(), logging.getLogger("tifffile")
    log.addHandler(h)
    try:
        yield h.found
    finally:
        log.removeHandler(h)


def _check_complete(t, series, problems):
    """Raise IncompleteImage when the file is shorter than its structure says: a tag's value or a level's pixels past
    the end. A pyramid stores its smaller levels last, so a file still being copied opens at full resolution but its
    zoomed-out views have nothing to read."""
    size = t.filehandle.size
    cut = [m for m in problems if "offset" in m.lower() or "truncat" in m.lower() or "corrupt" in m.lower()]
    if cut:
        raise IncompleteImage(f"the file ({size:,} bytes) is shorter than its own structure says")
    for k, lv in enumerate(series.levels):
        for pg in list(lv.pages)[:3]:  # (the pages read: the first, or the 3 colors of a planar RGB)
            offs, counts = getattr(pg, "dataoffsets", None), getattr(pg, "databytecounts", None)
            if offs is None or counts is None or not len(offs):
                continue
            end = max(int(o) + int(c) for o, c in zip(offs, counts))
            if end > size:
                raise IncompleteImage(f"the file ({size:,} bytes) ends before the pixels of pyramid level {k} "
                                      f"(up to {end:,} bytes)")


class Slide:
    """A pyramidal TIFF / OME-TIFF: RGB (H&E) or a fluorescence channel (e.g. DAPI, shown in gray with its display
    range taken from the coarsest level). channel: which channel of a file of several (default the first)."""

    def __init__(self, path, channel=None):
        self.path = path
        with _tiff_problems() as problems, tifffile.TiffFile(path) as t:
            s = t.series[0]
            self.levels = []
            for lv in s.levels:
                ax, sh = lv.axes, lv.shape
                self.levels.append((sh[ax.index("X")], sh[ax.index("Y")], ax))
            rgb = "S" in s.axes and s.shape[s.axes.index("S")] in (3, 4)
            rgb = rgb or ("C" in s.axes and s.shape[s.axes.index("C")] == 3 and s.dtype == np.uint8)
            self.rgb, self.dtype = rgb, s.dtype
            n_ch = s.shape[s.axes.index("C")] if "C" in s.axes and not rgb else 1
            self.channel = int(channel or 0)
            if not 0 <= self.channel < n_ch:
                raise ValueError(f"{Path(path).name} has {n_ch} channel(s): no channel {self.channel + 1}")
            _check_complete(t, s, problems)
        self.W0, self.H0 = self.levels[0][:2]
        self._z, self._stores, self._lock = {}, [], threading.Lock()
        self.fill = (245, 245, 245) if self.rgb else (0, 0, 0)
        self._range = None

    def _zarr(self, k):
        if k not in self._z:
            store = tifffile.imread(self.path, aszarr=True, level=k)
            self._stores.append(store)
            self._z[k] = zarr.open(store, mode="r")
        return self._z[k]

    def close(self):
        """Close the file (it opens again when read): Windows cannot move or delete a file that is open."""
        with self._lock:
            for st in self._stores:
                try:
                    st.close()
                except Exception:
                    pass
            self._stores, self._z = [], {}

    def _box(self, k, a0, b0, a1, b1, disp=None):
        """Pixels of level k in a box as (h, w, 3) uint8. disp: (low, high, gamma) of a fluorescence image (its raw
        values shown from black at low to white at high, with that gamma; the user's), else its own range."""
        ax = self.levels[k][2]
        idx = tuple(slice(b0, b1) if c == "Y" else slice(a0, a1) if c == "X" else
                    slice(None) if (c in "SC" and self.rgb) else self.channel if c == "C" else 0 for c in ax)
        with self._lock:  # tifffile stores share a file handle
            arr = np.asarray(self._zarr(k)[idx])
        kept = [c for c in ax if c in "YX" or (c in "SC" and self.rgb)]
        arr = np.moveaxis(arr, [kept.index(c) for c in ("Y", "X")], [0, 1])
        if self.rgb:
            return np.ascontiguousarray(arr[..., :3].astype(np.uint8))
        lo, hi = self.display_range()
        gamma = 1.0
        if disp:
            lo, hi, gamma = disp
        g = np.clip((arr.astype(np.float32) - lo) / max(hi - lo, 1e-6), 0, 1)
        if gamma != 1.0:
            g = g ** np.float32(gamma)
        return np.repeat((255 * g).astype(np.uint8)[..., None], 3, axis=2)

    def stats(self):
        """For the image's display settings: RGB or not; a fluorescence image's own range (display_range), its
        values' extent on the coarsest level and the largest value its type holds."""
        if self.rgb:
            return {"rgb": True}
        k = len(self.levels) - 1
        ax = self.levels[k][2]
        idx = tuple(slice(None) if c in "YX" else self.channel if c == "C" else 0 for c in ax)
        with self._lock:
            a = np.asarray(self._zarr(k)[idx])
        lo, hi = self.display_range()
        top = float(np.iinfo(self.dtype).max) if np.issubdtype(self.dtype, np.integer) else float(a.max())
        return {"rgb": False, "auto": [round(lo, 3), round(hi, 3)], "min": float(a.min()), "max": float(a.max()),
                "type_max": top}

    def display_range(self):
        if self._range is None:
            k = len(self.levels) - 1
            w, h, _ = self.levels[k]
            ax = self.levels[k][2]
            idx = tuple(slice(None) if c in "YX" else self.channel if c == "C" else 0 for c in ax)
            with self._lock:
                a = np.asarray(self._zarr(k)[idx]).astype(np.float32).ravel()
            a = a[a > 0] if (a > 0).any() else a
            self._range = (float(np.percentile(a, 1)), float(np.percentile(a, 99.7)))
        return self._range

    def read(self, x0, y0, x1, y1, out_w, out_h, disp=None):
        """Level-0 pixel box -> RGB image of (out_w, out_h), read from the coarsest adequate pyramid level (disp: see
        _box)."""
        need = (x1 - x0) / out_w
        k = max(i for i, (w, _, _) in enumerate(self.levels) if i == 0 or self.W0 / w <= need * 1.001)
        w, h, _ = self.levels[k]
        ds = self.W0 / w
        a0, b0 = max(0, int(x0 / ds)), max(0, int(y0 / ds))
        a1, b1 = min(w, int(math.ceil(x1 / ds))), min(h, int(math.ceil(y1 / ds)))
        canvas = Image.new("RGB", (out_w, out_h), self.fill)
        if a1 <= a0 or b1 <= b0:
            return canvas
        img = Image.fromarray(self._box(k, a0, b0, a1, b1, disp))
        tw, th = max(1, round((a1 - a0) * ds / need)), max(1, round((b1 - b0) * ds / need))
        canvas.paste(img.resize((tw, th), Image.BILINEAR), (round((a0 * ds - x0) / need), round((b0 * ds - y0) / need)))
        return canvas


def read_warped(sl, T, x0, y0, x1, y1, out_w, out_h, disp=None):
    """A box of another image's pixels shown in the frame's place: T (2x3) maps frame level-0 pixels to the image's
    level-0 pixels (any rotation, mirroring or scale between them); the box (in frame pixels) -> RGB (out_w, out_h)."""
    T = np.asarray(T, float)
    L, t = T[:, :2], T[:, 2]
    fx, fy = (x1 - x0) / out_w, (y1 - y0) / out_h  # frame pixels per output pixel
    corners = np.array([[x0, y0], [x1, y0], [x0, y1], [x1, y1]], float) @ L.T + t
    bx0, by0 = np.floor(corners.min(0))
    bx1, by1 = np.ceil(corners.max(0))
    bx0, by0, bx1, by1 = max(bx0, 0), max(by0, 0), min(bx1, sl.W0), min(by1, sl.H0)
    if bx1 <= bx0 or by1 <= by0:
        return Image.new("RGB", (out_w, out_h), sl.fill)
    s = max(np.sqrt(abs(np.linalg.det(L))) * np.sqrt(fx * fy), 1e-9)  # image pixels per output pixel
    rw, rh = max(1, int(math.ceil((bx1 - bx0) / s))), max(1, int(math.ceil((by1 - by0) / s)))
    region = sl.read(bx0, by0, bx1, by1, rw, rh, disp)
    sx, sy = (bx1 - bx0) / rw, (by1 - by0) / rh  # image pixels per region pixel
    # output (u, v) -> frame (x0 + u fx, y0 + v fy) -> image (L @ frame + t) -> region ((ix - bx0) / sx, (iy - by0) / sy)
    a, b_, c = L[0, 0] * fx / sx, L[0, 1] * fy / sx, (L[0, 0] * x0 + L[0, 1] * y0 + t[0] - bx0) / sx
    d, e, f = L[1, 0] * fx / sy, L[1, 1] * fy / sy, (L[1, 0] * x0 + L[1, 1] * y0 + t[1] - by0) / sy
    return region.transform((out_w, out_h), Image.AFFINE, (a, b_, c, d, e, f), resample=Image.BILINEAR,
                            fillcolor=tuple(sl.fill))


class BlankSlide:
    """No image: a plain background of unlimited size (the cells are drawn on it)."""
    W0 = H0 = 10 ** 9
    fill = (250, 250, 250)

    def read(self, x0, y0, x1, y1, out_w, out_h, disp=None):
        return Image.new("RGB", (out_w, out_h), self.fill)

    def stats(self):
        return {"rgb": True}


_slides = {}


def _key_path(k):
    return k[0] if isinstance(k, tuple) else k


def unload(path):
    """Forget an open image (its file handles and stores), every channel of it."""
    for k in [k for k in _slides if _key_path(k) == path]:
        _slides.pop(k, None)


def forget(folder):
    """Close and drop the images open from inside folder (before it is moved or deleted)."""
    folder = str(Path(folder).resolve())
    for k in [k for k in _slides if _key_path(k) and str(Path(_key_path(k)).resolve()).startswith(folder)]:
        s = _slides.pop(k)
        if hasattr(s, "close"):
            s.close()


def slide(path, channel=None):
    """The open image of path (one channel of a file of several: channel, from 0)."""
    if not path:
        return BlankSlide()
    key = (path, int(channel)) if channel else path
    if key not in _slides:
        _slides[key] = Slide(path, channel)
    return _slides[key]


_CHANNEL_NAMES = {}


def channel_names(path):
    """The channels of a fluorescence image of several (its OME names, else "channel 1", ...), or None: one channel,
    or RGB (an H&E)."""
    try:
        st = Path(path).stat()
    except OSError:
        return None
    key = (str(path), st.st_mtime_ns, st.st_size)
    if key not in _CHANNEL_NAMES:
        names = None
        try:
            with tifffile.TiffFile(path) as t:
                sr = t.series[0]
                rgb = "S" in sr.axes and sr.shape[sr.axes.index("S")] in (3, 4)
                rgb = rgb or ("C" in sr.axes and sr.shape[sr.axes.index("C")] == 3 and sr.dtype == np.uint8)
                n = sr.shape[sr.axes.index("C")] if "C" in sr.axes and not rgb else 1
                if n > 1:
                    found = re.findall(r'<Channel[^>]*?Name="([^"]+)"', t.ome_metadata or "")
                    names = found if len(found) == n else [f"channel {k + 1}" for k in range(n)]
        except Exception:
            names = None
        _CHANNEL_NAMES[key] = names
    return _CHANNEL_NAMES[key]


def _grid_axes(ax, g, box, labels=True, color="k", alpha=0.3):
    """Tile grid lines and spreadsheet-style labels (rows A.. on the left, columns 1.. on top) inside box."""
    x0, y0, x1, y1 = box
    tp = g["tile_px"]
    js = [j for j in range(g["nx"] + 1) if x0 - 1 <= g["x0"] + j * tp <= x1 + 1]
    is_ = [i for i in range(g["ny"] + 1) if y0 - 1 <= g["y0"] + i * tp <= y1 + 1]
    for j in js:
        ax.axvline(g["x0"] + j * tp, color=color, lw=0.4, alpha=alpha)
    for i in is_:
        ax.axhline(g["y0"] + i * tp, color=color, lw=0.4, alpha=alpha)
    if labels:
        cj = [j for j in range(g["nx"]) if x0 <= g["x0"] + (j + 0.5) * tp <= x1]
        ci = [i for i in range(g["ny"]) if y0 <= g["y0"] + (i + 0.5) * tp <= y1]
        ax.set_xticks([g["x0"] + (j + 0.5) * tp for j in cj], [g["cols"][j] for j in cj], fontsize=6)
        ax.set_yticks([g["y0"] + (i + 0.5) * tp for i in ci], [g["rows"][i] for i in ci], fontsize=6)
        ax.xaxis.tick_top()
        ax.tick_params(length=0, pad=1)
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)


def _lineage_dots(ax, px, py, lin, palette, highlight, s=0.4, max_n=250000, seed=0):
    idx = np.arange(len(px))
    if len(idx) > max_n:
        idx = np.random.default_rng(seed).choice(idx, max_n, replace=False)
    l = lin[idx]
    if highlight:
        m = ~np.isin(l, list(highlight))
        ax.scatter(px[idx][m], py[idx][m], s=s, c=[(0.78, 0.78, 0.78)], linewidths=0, rasterized=True)
    for k, col in palette.items():
        if highlight and k not in highlight:
            continue
        m = l == k
        ax.scatter(px[idx][m], py[idx][m], s=s * (2.5 if highlight else 1), c=[tuple(v / 255 for v in col)],
                   linewidths=0, rasterized=True)


def _legend(fig, palette, highlight, fontsize=7):
    keys = [k for k in palette if not highlight or k in highlight]
    hs = [plt.Line2D([], [], marker="o", ls="", color=tuple(v / 255 for v in palette[k]), label=k) for k in keys]
    if highlight:
        hs.append(plt.Line2D([], [], marker="o", ls="", color=(0.78, 0.78, 0.78), label="other cells"))
    fig.legend(handles=hs, loc="lower center", ncol=min(len(hs), 7), fontsize=fontsize, frameon=False,
               markerscale=fontsize / 7)
    return len(hs)


OVERVIEW_FONT = 13  # pt on a figure about 20 in wide: readable when the figure is shown a page wide


def tissue_overview(he_path, px, py, lin, palette, g, areas, out, highlight=None, title="", image_kind="H&E",
                    area_label="candidate areas", colored_by="lineage"):
    """Whole sample: the image (H&E or DAPI) with the tile grid (left), cells colored by lineage with candidate areas
    outlined (right). Without an image only the right panel is drawn."""
    box = (g["x0"], g["y0"], g["x0"] + g["nx"] * g["tile_px"], g["y0"] + g["ny"] * g["tile_px"])
    w, h = box[2] - box[0], box[3] - box[1]
    W = 1000
    Hh = max(1, round(W * h / w))
    n = 2 if he_path else 1
    F = OVERVIEW_FONT
    n_keys = len([k for k in palette if not highlight or k in highlight]) + bool(highlight)
    leg_in = -(-n_keys // 7) * F * 1.7 / 72 + 0.25  # legend rows below the panels
    fig_h = Hh / 100 + 1.6 + leg_in
    fig, axs = plt.subplots(1, n, figsize=(n * W / 100 + 0.6, fig_h), dpi=100, squeeze=False)
    ax = list(axs[0])
    if he_path:
        im = slide(he_path).read(*box, W, Hh)
        ax[0].imshow(np.asarray(im), extent=(box[0], box[2], box[3], box[1]))
        _grid_axes(ax[0], g, box)
        ax[0].tick_params(labelsize=F * 0.7)
        ax[0].set_title(f"{image_kind} ({g['tile_um']} um tiles)", fontsize=F, pad=16)
    ax = [None] * (2 - n) + ax
    ax[1].set_facecolor("white")
    _lineage_dots(ax[1], px, py, lin, palette, highlight)
    _grid_axes(ax[1], g, box)
    for a in areas:
        poly = np.array(a["outline_px"] + a["outline_px"][:1])
        for k in ([0, 1] if he_path else [1]):
            ax[k].plot(poly[:, 0], poly[:, 1], color="red", lw=1.1, ls="--")
        ax[1].text(poly[:, 0].min(), poly[:, 1].min(), a["id"], color="red", fontsize=F * 0.85, fontweight="bold",
                   va="bottom", ha="left")
    ax[1].tick_params(labelsize=F * 0.7)
    ax[1].set_title(("" if he_path else f"no image; {g['tile_um']} um tiles; ") + f"cells colored by {colored_by}" + (" (highlighted: " + ", ".join(sorted(highlight)) + ")" if highlight else "")
                    + (f"; red dashed = {area_label}" if areas else ""), fontsize=F, pad=16)
    _legend(fig, palette, highlight, fontsize=F)
    fig.suptitle(title, fontsize=F * 1.25, fontweight="bold")
    plt.tight_layout(rect=(0, leg_in / fig_h, 1, 0.97)); plt.savefig(out, dpi=100, pil_kwargs={"quality": 82}); plt.close()
    return str(out)


def area_view(he_path, px, py, lin, palette, out_of_place, box, outline_px, out, g, highlight=None, um_per_px=None,
              title="", image_kind="H&E"):
    """One area: the image (left) and cells colored by lineage on the faded image (right); out-of-place cells ringed.
    Without an image only the right panel is drawn."""
    x0, y0, x1, y1 = box
    W = 900
    Hh = max(1, round(W * (y1 - y0) / (x1 - x0)))
    if Hh > 900:
        W, Hh = max(1, round(900 * (x1 - x0) / (y1 - y0))), 900
    im = np.asarray(slide(he_path).read(x0, y0, x1, y1, W, Hh))
    m = (px >= x0) & (px <= x1) & (py >= y0) & (py <= y1)
    n = 2 if he_path else 1
    fig, axs = plt.subplots(1, n, figsize=(n * W / 100 + 0.6, Hh / 100 + 1.3), dpi=100, squeeze=False)
    ax = [None] * (2 - n) + list(axs[0])
    if he_path:
        ax[0].imshow(im, extent=(x0, x1, y1, y0))
        ax[1].imshow(im, extent=(x0, x1, y1, y0), alpha=0.45)
    sel = np.where(m)[0]
    _lineage_dots(ax[1], px[sel], py[sel], lin[sel], palette, highlight, s=5 if len(sel) < 20000 else 2)
    o = sel[out_of_place[sel]]
    ax[1].scatter(px[o], py[o], s=40, facecolors="none", edgecolors="black", linewidths=0.8)
    poly = np.array(outline_px + outline_px[:1])
    for a in ax:
        if a is None:
            continue
        _grid_axes(a, g, box, alpha=0.25)
        a.plot(poly[:, 0], poly[:, 1], color="red", lw=1.2, ls="--")
    width_um = (x1 - x0) * (um_per_px or config.HE_UM_PER_PX)
    if he_path:
        ax[0].set_title(f"{image_kind}, {width_um:.0f} um wide", fontsize=8, pad=14)
    ax[1].set_title(("" if he_path else f"no image, {width_um:.0f} um wide; ")
                    + "cells colored by lineage; black rings = out-of-place cells", fontsize=8, pad=14)
    _legend(fig, palette, highlight)
    fig.suptitle(title, fontsize=9)
    plt.tight_layout(rect=(0, 0.05, 1, 0.98)); plt.savefig(out, dpi=100, pil_kwargs={"quality": 82}); plt.close()
    return str(out)


def he_crop(path, cx, cy, half_px):
    """Full-resolution crop around (cx, cy), clipped to the image: (array, x0, y0)."""
    sl = slide(path)
    x0, x1 = int(max(cx - half_px, 0)), int(min(cx + half_px, sl.W0))
    y0, y1 = int(max(cy - half_px, 0)), int(min(cy + half_px, sl.H0))
    return np.asarray(sl.read(x0, y0, x1, y1, max(1, x1 - x0), max(1, y1 - y0))), x0, y0


def region_colors(names, chosen=None):
    """Fixed colors for the canonical names, a generated palette for user-drawn region names; chosen = colors the
    user picked ({name: "#rrggbb"}) win."""
    import colorsys
    out, k = {}, 0
    chosen = chosen or {}
    for n in names:
        if n in REGION_COLORS:
            out[n] = REGION_COLORS[n]
        else:  # every user-drawn name takes its turn, so picking one color leaves the others as they were
            h = (0.07 + k * 0.618034) % 1.0
            if h > 0.93 or h < 0.065:  # no red: the viewer outlines the selected region in bright red
                h = (h + 0.12) % 1.0
            r, g, b = colorsys.hls_to_rgb(h, 0.5, 0.65)
            out[n] = "#%02x%02x%02x" % (int(255 * r), int(255 * g), int(255 * b))
            k += 1
        out[n] = _not_red(chosen.get(n, out[n]))
    return out


def _not_red(hexcol):
    """A color kept clear of red (the viewer outlines the selected region in bright red): a reddish one, e.g. QuPath's
    "Tumor" class or a color picked before, is turned toward orange or magenta, keeping its lightness."""
    import colorsys
    try:
        r, g, b = (int(hexcol[i:i + 2], 16) / 255 for i in (1, 3, 5))
    except (ValueError, TypeError, IndexError):
        return hexcol
    h, l, sat = colorsys.rgb_to_hls(r, g, b)
    if sat < 0.35 or not (h > 0.93 or h < 0.065):
        return hexcol
    h = 0.08 if h < 0.5 else 0.88  # toward orange, or toward magenta
    r, g, b = colorsys.hls_to_rgb(h, l, sat)
    return "#%02x%02x%02x" % (int(255 * r), int(255 * g), int(255 * b))


def region_map(cells, region, out, title="", highlight=None, highlight_label="conflicts", chosen=None):
    fig, ax = plt.subplots(figsize=(9, 8))
    names = [n for n in dict.fromkeys(region.tolist()) if n != "unknown"] + ["unknown"]
    for name, col in region_colors(names, chosen).items():
        m = region == name
        if m.any():
            ax.scatter(cells.x.to_numpy()[m], cells.y.to_numpy()[m], s=0.05, c=col, label=f"{name} ({m.sum()})", rasterized=True)
    if highlight is not None and len(highlight):
        ax.scatter(cells.x.to_numpy()[highlight], cells.y.to_numpy()[highlight], s=4, c="k",
                   label=f"{highlight_label} ({len(highlight)})")
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.set_title(title, fontsize=9)
    ax.legend(markerscale=12, fontsize=7, loc="lower left")
    plt.tight_layout(); plt.savefig(out, dpi=80); plt.close()
    return str(out)


def finding_map(x, y, rows, of_type, out, title="", max_n=150000, seed=0, width_in=4.6, dpi=110):
    """A small map of the whole sample: every cell grey, the other cells of the finding's type blue, the finding's
    cells (rows) red (the Refine window's findings against the priors). Drawn without pyplot: the window asks for
    several at once, in parallel threads, and pyplot's shared state is not thread-safe (it drew blank images)."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    w, h = float(np.ptp(x)) or 1.0, float(np.ptp(y)) or 1.0
    fw = width_in
    fig = Figure(figsize=(fw, min(max(fw * h / w, 1.6 * fw / 4.6), 5.0 * fw / 4.6) + 0.3))
    FigureCanvasAgg(fig)
    ax = fig.add_subplot(111)
    n = len(x)
    keep = np.random.default_rng(seed).choice(n, max_n, replace=False) if n > max_n else slice(None)
    k = (fw / 4.6) ** 2 * (110.0 / dpi) ** 2  # (a larger map: every dot as large on it)
    b0 = max(0.1, min(1.5, 30000.0 / n))  # (dot size by the number of cells: the tissue reads as a grey shape)
    base = b0 * k
    ax.scatter(x[keep], y[keep], s=base, c="#b8b8b8", rasterized=True, linewidths=0)
    rest = np.setdiff1d(of_type, rows) if of_type is not None and rows is not None else np.zeros(0, int)
    if len(rest):
        ax.scatter(x[rest], y[rest], s=max(1.0, 3 * b0) * k, c="#4a90c8", rasterized=True, linewidths=0)
    if rows is not None and len(rows):
        ax.scatter(x[rows], y[rows], s=(max(2.0, 4 * b0) if len(rows) > 500 else 8) * k, c="#d62728", rasterized=True, linewidths=0)
    ax.set_aspect("equal"); ax.invert_yaxis(); ax.axis("off")
    ax.set_title(title, fontsize=8 * max(1.0, fw / 6))
    fig.tight_layout(pad=0.3)
    fig.savefig(out, dpi=dpi)
    return str(out)


def conflict_sheet(cells, rows, he_path, labels, out, half_um=120.0, title="", lineage=None, palette=None, um_per_px=None):
    """Two panels per conflict cell: raw H&E (top) and the cells colored by their current lineage (bottom)."""
    half = half_um / (um_per_px or config.HE_UM_PER_PX)
    n = len(rows)
    fig, ax = plt.subplots(2, n, figsize=(3 * n, 6.6), squeeze=False)
    px, py = cells.px.to_numpy(), cells.py.to_numpy()
    for j, r in enumerate(rows):
        im, x0, y0 = he_crop(he_path, px[r["row"]], py[r["row"]], half)
        m = (px >= x0) & (px < x0 + im.shape[1]) & (py >= y0) & (py < y0 + im.shape[0])
        for a in ax[:, j]:
            a.imshow(im); a.set_xticks([]); a.set_yticks([])
        lin = lineage[m]
        for k, col in palette.items():
            mm = lin == k
            ax[1, j].scatter(px[m][mm] - x0, py[m][mm] - y0, s=4, c=[tuple(v / 255 for v in col)])
        for a in ax[:, j]:
            a.scatter([px[r["row"]] - x0], [py[r["row"]] - y0], s=250, facecolors="none", edgecolors="cyan", linewidths=1.5)
        ax[0, j].set_title(f"#{j} cell {r['cell']}\n{r['region']} {r['category']}", fontsize=7)
        ax[1, j].set_title(f"expr: {r['base_label']} ({r['base_conf']:.2f})\nnow: {labels[r['row']]}", fontsize=7)
    fig.suptitle(title + "  [cyan = cell; dots = cells colored by current lineage]", fontsize=8)
    plt.tight_layout(); plt.savefig(out, dpi=70, pil_kwargs={"quality": 80}); plt.close()
    return str(out)


def embedding_plot(Y, labels, out, title="", max_n=200000, seed=0):
    """A UMAP / PCA scatter colored by a grouping (the largest groups in DISTINCT colors, the rest gray)."""
    idx = np.arange(len(Y))
    if len(idx) > max_n:
        idx = np.random.default_rng(seed).choice(idx, max_n, replace=False)
    lab = np.asarray(labels, dtype=object)[idx].astype(str)
    keys = pd_value_counts(lab)
    pal = {k: DISTINCT[i % len(DISTINCT)] for i, k in enumerate(keys[:len(DISTINCT) * 2])}
    fig, ax = plt.subplots(figsize=(9, 7.5), dpi=100)
    for k in keys:
        m = lab == k
        ax.scatter(Y[idx][m, 0], Y[idx][m, 1], s=0.6, c=pal.get(k, "#cccccc"), linewidths=0, rasterized=True)
    hs = [plt.Line2D([], [], marker="o", ls="", color=pal[k], label=str(k)) for k in keys if k in pal]
    ax.legend(handles=hs, loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7, frameon=False, markerscale=1.5)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    plt.tight_layout(); plt.savefig(out, dpi=100, pil_kwargs={"quality": 85}); plt.close()
    return str(out)


def pd_value_counts(a):
    """Distinct values of a, most frequent first."""
    v, n = np.unique(np.asarray(a, dtype=str), return_counts=True)
    return [x for _, x in sorted(zip(-n, v))]
