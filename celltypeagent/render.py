"""Figures: region maps (PNG) and H&E conflict crops (JPEG, smaller for MCP clients). Each returns its path."""
import numpy as np, tifffile, zarr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from . import config

FIXED_LINEAGE_RGB = {"epi": (235, 130, 20), "fib": (60, 165, 60), "stroma": (40, 160, 170), "lymph": (50, 100, 220),
                     "myeloid": (210, 60, 170), "Unassigned": (150, 150, 150)}


def lineage_palette(lineages):
    """RGB per lineage name: fixed for the reference lineages, golden-ratio hues for a user's marker tree."""
    import colorsys
    out, k = {}, 0
    for n in lineages:
        if n in FIXED_LINEAGE_RGB:
            out[n] = FIXED_LINEAGE_RGB[n]
        else:
            out[n] = tuple(int(255 * c) for c in colorsys.hls_to_rgb((0.05 + k * 0.618034) % 1.0, 0.5, 0.75))
            k += 1
    return out
REGION_COLORS = {"mucosa": "#e6862a", "submucosa": "#57a35a", "muscularis_propria": "#3d76b8",
                 "lymphoid_aggregate": "#8e5bb5", "fibromuscular_hyperplasia": "#c9b63c", "granuloma": "#c4473f",
                 "unknown": "#cfcfcf"}
_he = {}


def _he_array(path):
    if path not in _he:
        _he[path] = zarr.open(tifffile.imread(path, aszarr=True, level=0), mode="r")
    return _he[path]


def he_crop(path, cx, cy, half_px):
    z = _he_array(path)
    chw = z.shape[0] == 3
    H, W = (z.shape[1], z.shape[2]) if chw else (z.shape[0], z.shape[1])
    x0, x1 = int(max(cx - half_px, 0)), int(min(cx + half_px, W))
    y0, y1 = int(max(cy - half_px, 0)), int(min(cy + half_px, H))
    im = np.asarray(z[:, y0:y1, x0:x1]).transpose(1, 2, 0) if chw else np.asarray(z[y0:y1, x0:x1])
    return im, x0, y0


def region_colors(names):
    """Fixed colors for the canonical names, a generated palette for user-drawn region names."""
    import colorsys
    out, k = {}, 0
    for n in names:
        if n in REGION_COLORS:
            out[n] = REGION_COLORS[n]
        else:
            r, g, b = colorsys.hls_to_rgb((0.07 + k * 0.618034) % 1.0, 0.5, 0.65)
            out[n] = "#%02x%02x%02x" % (int(255 * r), int(255 * g), int(255 * b))
            k += 1
    return out


def region_map(cells, region, out, title="", highlight=None, highlight_label="conflicts"):
    fig, ax = plt.subplots(figsize=(9, 8))
    names = [n for n in dict.fromkeys(region.tolist()) if n != "unknown"] + ["unknown"]
    for name, col in region_colors(names).items():
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
