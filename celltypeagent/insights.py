"""Evidence for Claude's high-level review of an annotation (the six kinds of advice in advice.py):

  compartments    compartment(): outlines of where a lineage forms tissue (mucosa around the epithelium, the
                  muscularis where smooth muscle is dense, follicles where B cells cluster), from the cells
  taxonomy        type_confusion(): pairs of types the expression model cannot tell apart (merge candidates);
                  type_structure(): sub-clusters of one type with their genes and where they sit (split candidates)
  organization    spatial_profile(): per type, distance to the epithelium and to the muscle, regions, and the types
                  it neighbors more (or less) than expected
  context         sample_context(): what the sample is and what it is made of
  QC              qc_overview(): transcripts, cell density, Unassigned and mixed-lineage cells per tile, and marker
                  spillover next to the epithelium; registration(): offset between the nuclei on the image and the
                  cell centroids
  report          annotation_facts(): how the annotation was made, for the report
All of it reads the annotation; nothing is changed.
"""
import re
import numpy as np, pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree
from . import plausibility, markers, model, render

GRID_UM = 40.0  # raster of compartment outlines


def _lineage(s, version=None):
    r = s.result(version)
    lut = np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)
    return r, np.array(s.types, dtype=object)[r["pred"]], lut[r["pred"]]


def _members(s, names, types, lineages):
    """Cells whose type, group or lineage is one of names."""
    want = set()
    for n in names:
        want |= set(s.type_groups.get(n, [n]))
    return np.isin(types, list(want)) | np.isin(lineages, list(names))


def _xy(s):
    return s.cells[["x", "y"]].to_numpy()


def _core_distance(s, mask, min_share=0.5):
    """Distance (um) of every cell to the nearest 'core' cell of mask (a member whose 15 nearest cells are mostly
    members), i.e. to where that lineage forms tissue; inf if there is none."""
    share = mask[s.knn()].mean(1)
    core = mask & (share >= min_share)
    if not core.any():
        return np.full(len(mask), np.inf)
    return cKDTree(_xy(s)[core]).query(_xy(s))[0]


# ---------------- compartments ----------------
def compartment(s, lineages, max_dist_um=30.0, min_share=0.5, min_cells=200, version=None):
    """Where the given lineages / types form tissue: cells within max_dist_um of their core, closed and hole-filled
    on a GRID_UM raster; one polygon per piece holding >= min_cells cells (H&E level-0 pixels)."""
    r, types, lin = _lineage(s, version)
    mask = _members(s, lineages, types, lin)
    if not mask.any():
        raise ValueError(f"no cells of {lineages}")
    inside = _core_distance(s, mask, min_share) <= max_dist_um
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    step = GRID_UM / s.he_um_per_px
    x0, y0 = px.min() - 2 * step, py.min() - 2 * step
    ix, iy = ((px - x0) / step).astype(int), ((py - y0) / step).astype(int)
    grid = np.zeros((iy.max() + 3, ix.max() + 3), bool)
    grid[iy[inside], ix[inside]] = True
    grid = ndimage.binary_fill_holes(ndimage.binary_closing(grid, iterations=2))
    lab, n = ndimage.label(grid)
    cell_lab = lab[iy, ix]
    pieces, polys = [], []
    for k in range(1, n + 1):
        cells = cell_lab == k
        if cells.sum() < min_cells:
            continue
        pix = {(int(i), int(j)) for i, j in zip(*np.where(lab == k))}
        poly = plausibility._simplify(plausibility.outline({"x0": x0, "y0": y0, "tile_px": step}, pix), step)
        polys.append(poly)
        pieces.append(cells)
    if not polys:
        raise ValueError("no piece of that compartment holds enough cells; lower min_cells or raise max_dist_um")
    covered = np.any(pieces, axis=0)
    comp = pd.Series(lin[covered]).value_counts(normalize=True).head(6).round(3)
    area_mm2 = float(sum((lab == k).sum() for k in range(1, n + 1) if (cell_lab == k).sum() >= min_cells) * (GRID_UM / 1000) ** 2)
    return {"polygons_px": polys, "n_pieces": len(polys), "n_cells": int(covered.sum()),
            "share_of_members_inside": round(float(mask[covered].sum() / mask.sum()), 3),
            "composition": comp.to_dict(), "area_mm2": round(area_mm2, 2), "cells": covered}


def fig_compartment(s, comp, name, version=None):
    """The proposed outlines on the overview (H&E left, lineages right), for checking against the image."""
    sc = plausibility.scan(s, version)
    _, _, lin = _lineage(s, version)
    areas = [{"id": f"{name} {i + 1}" if len(comp["polygons_px"]) > 1 else name, "outline_px": p}
             for i, p in enumerate(comp["polygons_px"])]
    out = s.dir / "figures" / f"compartment_{re.sub(r'[^A-Za-z0-9]+', '_', name)}_v{sc['version']}.jpg"
    return render.tissue_overview(s.meta["he_path"], s.cells.px.to_numpy(), s.cells.py.to_numpy(), lin,
                                  s.lineage_colors(), sc["grid"], areas, out, title=f"proposed {name}",
                                  image_kind=s.meta.get("image_kind", "H&E"), area_label="proposed outline")


# ---------------- taxonomy ----------------
def type_confusion(s, top=15, version=None):
    """Type pairs the expression model mixes up: for cells of type A, the mean probability of type B, both ways."""
    r = s.result(version)
    P = np.exp(s.base(r["state"])["logp"].astype(np.float32))
    T = len(s.types)
    one = np.zeros((T, len(r["pred"])), np.float32)
    one[r["pred"], np.arange(len(r["pred"]))] = 1
    n = one.sum(1)
    M = (one @ P) / np.maximum(n, 1)[:, None]                       # M[a, b] = mean p(b) over cells labelled a
    second = np.argsort(-P, axis=1)[:, 1]
    rows = []
    for a in range(T):
        for b in range(a + 1, T):
            if n[a] < 30 or n[b] < 30:
                continue
            score = M[a, b] + M[b, a]
            rows.append({"types": [s.types[a], s.types[b]], "mutual_confusion": round(float(score), 3),
                         "a_as_b": round(float(M[a, b]), 3), "b_as_a": round(float(M[b, a]), 3),
                         "share_a_second_b": round(float(np.mean(second[r["pred"] == a] == b)), 3),
                         "same_lineage": s.lineage_of_type.get(s.types[a]) == s.lineage_of_type.get(s.types[b]),
                         "n": [int(n[a]), int(n[b])]})
    rows.sort(key=lambda x: -x["mutual_confusion"])
    return {"pairs": rows[:top], "note": "mutual_confusion = mean p(B) in A cells + mean p(A) in B cells; above "
                                         "~0.3 the panel hardly separates them (merge, or sharper markers)"}


def type_structure(s, cell_type, resolution=0.5, max_cells=30000, version=None):
    """Sub-clusters of one type: size, genes that set each apart, and where it lies (distance to the epithelium,
    regions, how strongly its cells neighbor each other)."""
    r, types, lin = _lineage(s, version)
    rows = np.where(types == cell_type)[0]
    if len(rows) < 100:
        raise ValueError(f"{cell_type}: {len(rows)} cells, too few to sub-cluster")
    if len(rows) > max_cells:
        rows = np.sort(np.random.default_rng(0).choice(rows, max_cells, replace=False))
    Xn = model.lognorm(s.counts[rows])
    g, nbr = markers._graph(Xn)
    sub = markers._leiden(g, nbr, resolution, 0)
    det = (s.counts[rows] > 0).astype(np.float32).tocsr()
    depth = s.cells.depth_um.to_numpy()[rows]
    reg = r["region"][rows]
    xy = _xy(s)[rows]
    nn = cKDTree(xy).query(xy, k=min(11, len(rows)))[1][:, 1:]
    out = []
    for c in np.unique(sub):
        m = sub == c
        if m.sum() < 0.05 * len(rows):
            continue
        a = np.asarray(det[m].mean(0)).ravel()
        b = np.asarray(det[~m].mean(0)).ravel()
        sc = (a - b) * np.log2((a + 0.02) / (b + 0.02))
        genes = [{"gene": s.panel[i], "in": round(float(a[i]), 2), "rest": round(float(b[i]), 2)}
                 for i in np.argsort(-sc)[:6] if a[i] > b[i]]
        coherence = float((sub[nn[m]] == c).mean() / max(m.mean(), 1e-6))
        out.append({"subcluster": f"{cell_type} /{int(c)}", "n_cells": int(m.sum()), "share": round(float(m.mean()), 3),
                    "top_genes": genes,
                    "depth_um_median": None if np.isnan(depth[m]).all() else round(float(np.nanmedian(depth[m])), 1),
                    "regions": pd.Series(reg[m]).value_counts(normalize=True).head(3).round(2).to_dict(),
                    "spatial_clustering": round(coherence, 2)})
    return {"cell_type": cell_type, "n_cells": int((types == cell_type).sum()), "subclusters": out,
            "note": "spatial_clustering = how much more often a sub-cluster's cells neighbor each other (among this "
                    "type's cells) than chance (1 = mixed, >1.5 = sits apart: a candidate subtype)"}


# ---------------- spatial organization ----------------
def spatial_profile(s, version=None, min_cells=50):
    """Per type: distance to the epithelium and to the smooth muscle (quartiles, um), regions, and the types among
    its 15 nearest cells more / less often than expected."""
    r, types, lin = _lineage(s, version)
    depth = s.cells.depth_um.to_numpy()
    smc = _members(s, [l for l in s.lineages if "muscle" in l.lower()] or ["Smooth muscle cells"], types, lin)
    dmus = _core_distance(s, smc) if smc.any() else np.full(len(types), np.nan)
    code = r["pred"]
    T = len(s.types)
    nn = s.knn()
    counts = np.zeros((T, T))
    np.add.at(counts, (np.repeat(code, nn.shape[1]), code[nn].ravel()), 1)
    freq = np.bincount(code, minlength=T) / len(code)
    enrich = np.log2((counts / np.maximum(counts.sum(1, keepdims=True), 1) + 1e-4) / (freq[None, :] + 1e-4))
    q = lambda v: None if np.all(~np.isfinite(v)) else [round(float(x), 1) for x in np.nanpercentile(v[np.isfinite(v)], [25, 50, 75])]
    out = []
    for t in range(T):
        m = code == t
        if m.sum() < min_cells:
            continue
        order = np.argsort(-enrich[t])
        out.append({"cell_type": s.types[t], "n": int(m.sum()), "lineage": s.lineage_of_type.get(s.types[t]),
                    "dist_epithelium_um_q": q(depth[m]), "dist_muscle_um_q": q(dmus[m]),
                    "regions": pd.Series(r["region"][m]).value_counts(normalize=True).head(3).round(2).to_dict(),
                    "neighbors_more": [{"type": s.types[o], "log2": round(float(enrich[t, o]), 2)} for o in order[:4]
                                       if o != t and enrich[t, o] > 0.5][:3],
                    "neighbors_less": [{"type": s.types[o], "log2": round(float(enrich[t, o]), 2)}
                                       for o in order[::-1][:6] if freq[o] > 0.01 and enrich[t, o] < -1][:3]})
    return {"types": out, "note": "distances are from the cell to where epithelium / smooth muscle form tissue; "
                                  "log2 = neighbor share vs the type's share of all cells"}


# ---------------- context ----------------
def sample_context(s, version=None):
    r, types, lin = _lineage(s, version)
    meta = {k: v for k, v in s.meta.items() if k not in ("he_path", "bundle_dir", "source", "has_pathologist_regions")}
    return {"sample": meta, "method": s.scheme["kind"], "n_cells": len(types),
            "lineages": pd.Series(lin).value_counts(normalize=True).round(3).to_dict(),
            "top_types": pd.Series(types).value_counts(normalize=True).head(15).round(3).to_dict(),
            "regions_drawn": s.region_names(version),
            "source": s.meta.get("source", {}).get("kind")}


# ---------------- QC ----------------
def qc_overview(s, version=None):
    """Per tile: median transcripts per cell, cells per mm2, Unassigned and mixed-lineage shares; hotspots; and
    epithelial-marker spillover into nearby non-epithelial cells."""
    r, types, lin = _lineage(s, version)
    g = plausibility.grid(s)
    ti, tj = plausibility.tile_index(s, g)
    key = ti * g["nx"] + tj
    P = np.exp(s.base(r["state"])["logp"].astype(np.float32))
    lut = np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)
    names = sorted(set(lut))
    L = np.stack([P[:, lut == n].sum(1) for n in names], 1)
    top2 = np.sort(L, 1)[:, -2:]
    mixed = (top2[:, 0] >= 0.3) & (top2[:, 1] >= 0.3)
    df = pd.DataFrame({"tile": key, "counts": s.cells.total_counts.to_numpy(), "unassigned": types == markers.UNASSIGNED,
                       "low_conf": r["conf"] < 0.5, "mixed": mixed})
    agg = df.groupby("tile").agg(n=("counts", "size"), counts=("counts", "median"), unassigned=("unassigned", "mean"),
                                 low_conf=("low_conf", "mean"), mixed=("mixed", "mean"))
    agg = agg[agg.n >= 50]
    tile_mm2 = (g["tile_um"] / 1000) ** 2
    agg["density"] = agg.n / tile_mm2
    name = lambda k: plausibility.tile_name(g, int(k) // g["nx"], int(k) % g["nx"])
    med = agg.median()
    have = set(agg.index)  # tiles at the tissue edge are only partly covered: no density call there
    interior = np.array([all(k + d in have for d in (-1, 1, -g["nx"], g["nx"])) for k in agg.index])
    hot = {"low_transcripts": [name(k) for k in agg.index[agg.counts < 0.4 * med.counts]][:12],
           "low_density": [name(k) for k in agg.index[(agg.density < 0.35 * med.density) & interior]][:12],
           "many_unassigned": [name(k) for k in agg.index[agg.unassigned > max(3 * med.unassigned, 0.15)]][:12],
           "many_mixed": [name(k) for k in agg.index[agg.mixed > max(3 * med.mixed, 0.1)]][:12]}
    spill = None
    gi = {x.upper(): i for i, x in enumerate(s.panel)}
    idx = [gi[x] for x in ("EPCAM", "KRT8", "KRT19", "KRT20") if x in gi]
    epi = _members(s, [l for l in s.lineages if "epith" in l.lower()], types, lin)
    if idx and epi.any():
        d = s.cells.depth_um.to_numpy()
        hit = np.asarray((s.counts[:, idx] > 0).sum(1)).ravel() > 0
        bins = {"<10 um": d < 10, "10-30 um": (d >= 10) & (d < 30), ">60 um": d > 60}
        spill = {k: round(float(hit[~epi & m].mean()), 3) for k, m in bins.items() if (~epi & m).sum() > 100}
    return {"tile_um": g["tile_um"], "sample_median": {"transcripts_per_cell": float(med.counts),
                                                       "cells_per_mm2": round(float(med.density)),
                                                       "unassigned": round(float(med.unassigned), 3),
                                                       "mixed_lineage": round(float(med.mixed), 3)},
            "hotspot_tiles": hot,
            "epithelial_marker_in_non_epithelial_cells_by_distance": spill,
            "note": "mixed = the model gives two lineages >= 0.3 each (doublets, spillover or transitions); a share "
                    "of non-epithelial cells with epithelial markers that falls with distance indicates spillover"}


def fig_qc(s, version=None):
    """Four tile maps: transcripts per cell, cell density, Unassigned share, mixed-lineage share."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    r, types, lin = _lineage(s, version)
    g = plausibility.grid(s)
    ti, tj = plausibility.tile_index(s, g)
    P = np.exp(s.base(r["state"])["logp"].astype(np.float32))
    lut = np.array([s.lineage_of_type.get(t, t) for t in s.types], dtype=object)
    L = np.stack([P[:, lut == n].sum(1) for n in sorted(set(lut))], 1)
    top2 = np.sort(L, 1)[:, -2:]
    vals = {"median transcripts per cell": s.cells.total_counts.to_numpy().astype(float),
            "cells per tile": np.ones(len(types)), "Unassigned share": (types == markers.UNASSIGNED).astype(float),
            "mixed-lineage share": ((top2[:, 0] >= 0.3) & (top2[:, 1] >= 0.3)).astype(float)}
    fig, axs = plt.subplots(2, 2, figsize=(12, 12 * g["ny"] / max(g["nx"], 1) * 0.55 + 1.5))
    n = np.zeros((g["ny"], g["nx"]))
    np.add.at(n, (ti, tj), 1)
    for ax, (title, v) in zip(axs.ravel(), vals.items()):
        acc = np.zeros((g["ny"], g["nx"]))
        if title.startswith("median"):
            df = pd.DataFrame({"i": ti, "j": tj, "v": v}).groupby(["i", "j"]).v.median()
            for (i, j), x in df.items():
                acc[i, j] = x
        elif title == "cells per tile":
            acc = n.copy()
        else:
            np.add.at(acc, (ti, tj), v)
            acc = acc / np.maximum(n, 1)
        acc[n < 50] = np.nan
        im = ax.imshow(acc, cmap="viridis", interpolation="nearest")
        ax.set_title(title, fontsize=9)
        ax.set_xticks(range(g["nx"]), g["cols"], fontsize=5); ax.set_yticks(range(g["ny"]), g["rows"], fontsize=5)
        ax.xaxis.tick_top()
        plt.colorbar(im, ax=ax, fraction=0.03)
    out = s.dir / "figures" / f"qc_v{r['state']['version']}.jpg"
    plt.tight_layout(); plt.savefig(out, dpi=90, pil_kwargs={"quality": 82}); plt.close()
    return str(out)


def _nuclei(img, kind):
    """Nuclear signal of an image tile: hematoxylin optical density (Ruifrok-Johnston) of an H&E, or the DAPI."""
    if kind != "H&E":
        return img[..., 0]
    od = -np.log((img + 1) / 256.0)
    M = np.array([[0.65, 0.70, 0.29], [0.07, 0.99, 0.11], [0.27, 0.57, 0.78]])
    M = M / np.linalg.norm(M, axis=1, keepdims=True)
    return (od.reshape(-1, 3) @ np.linalg.inv(M))[:, 0].reshape(img.shape[:2])


def registration(s, n_windows=40, window_um=400.0, max_shift_um=25.0):
    """Offset between the nuclei on the image (hematoxylin of the H&E, or DAPI) and the cell centroids, from the
    cross-correlation of the two in windows spread over the tissue: the mean over windows (one window alone is noisy,
    crypts repeat), and an affine fit of the per-window offsets against position, which shows a scale or rotation
    error (offsets that grow across the section) that a single shift misses."""
    kind = s.meta.get("image_kind", "H&E")
    if not s.meta.get("he_path") or kind == "none":
        return {"error": "no image"}
    sl = render.slide(s.meta["he_path"])
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    g = plausibility.grid(s)
    ti, tj = plausibility.tile_index(s, g)
    cnt = pd.Series(ti * g["nx"] + tj).value_counts()
    tiles = cnt.index[cnt.to_numpy() >= 0.3 * cnt.iloc[0]].to_numpy()  # tissue-dense tiles all over the section
    tiles = np.random.default_rng(0).choice(tiles, min(n_windows, len(tiles)), replace=False)
    half = window_um / s.he_um_per_px / 2
    W = int(window_um)            # 1 um per pixel
    R = int(max_shift_um)
    c0 = W // 2
    total = np.zeros((2 * R + 1, 2 * R + 1))
    shifts = []
    for k in tiles:
        i, j = int(k) // g["nx"], int(k) % g["nx"]
        cx, cy = g["x0"] + (j + 0.5) * g["tile_px"], g["y0"] + (i + 0.5) * g["tile_px"]
        img = np.asarray(sl.read(cx - half, cy - half, cx + half, cy + half, W, W)).astype(np.float32)
        nuc = ndimage.gaussian_filter(_nuclei(img, kind), 2)
        m = (px > cx - half) & (px < cx + half) & (py > cy - half) & (py < cy + half)
        if m.sum() < 100:
            continue
        pts = np.zeros((W, W))
        u = ((px[m] - (cx - half)) / (2 * half) * W).astype(int).clip(0, W - 1)
        v = ((py[m] - (cy - half)) / (2 * half) * W).astype(int).clip(0, W - 1)
        np.add.at(pts, (v, u), 1)
        pts = ndimage.gaussian_filter(pts, 3)
        nuc, pts = (nuc - nuc.mean()) / (nuc.std() + 1e-9), (pts - pts.mean()) / (pts.std() + 1e-9)
        cc = np.fft.fftshift(np.real(np.fft.ifft2(np.fft.fft2(nuc) * np.conj(np.fft.fft2(pts))))) / W ** 2
        win = cc[c0 - R:c0 + R + 1, c0 - R:c0 + R + 1]
        total += win
        dy, dx = np.unravel_index(np.argmax(win), win.shape)
        shifts.append({"tile": plausibility.tile_name(g, i, j), "dx_um": int(dx - R), "dy_um": int(dy - R),
                       "n_cells": int(m.sum()), "correlation": round(float(win.max()), 3),
                       "_mm": ((cx - g["x0"]) * s.he_um_per_px / 1000, (cy - g["y0"]) * s.he_um_per_px / 1000)})
    if not shifts:
        return {"error": "no window with enough cells"}
    total /= len(shifts)
    dy, dx = np.unravel_index(np.argmax(total), total.shape)
    dx, dy = int(dx - R), int(dy - R)
    edge = abs(dx) == R or abs(dy) == R
    off = np.array([[x["dx_um"], x["dy_um"]] for x in shifts], float)
    pos = np.array([x.pop("_mm") for x in shifts])
    spread = float(np.median(np.hypot(*(off - [dx, dy]).T)))
    fit = None
    if len(shifts) >= 8:
        X = np.c_[pos, np.ones(len(pos))]
        coef = np.linalg.lstsq(X, off, rcond=None)[0]   # offset = pos_mm @ A + b
        resid = float(np.median(np.hypot(*(off - X @ coef).T)))
        A = coef[:2].T / 1000                            # um of offset per um of position
        w, h = (np.array([g["nx"], g["ny"]]) * g["tile_um"]) / 1000
        corners = np.array([[0, 0], [w, 0], [0, h], [w, h]])
        at = np.c_[corners, np.ones(4)] @ coef
        fit = {"scale_x": round(1 + float(A[0, 0]), 5), "scale_y": round(1 + float(A[1, 1]), 5),
               "rotation_deg": round(float(np.degrees((A[1, 0] - A[0, 1]) / 2)), 3),
               "median_residual_um": round(resid, 1),
               "offset_at_corners_um": {n: [round(float(a), 1) for a in v]
                                        for n, v in zip(("top left", "top right", "bottom left", "bottom right"), at)},
               "max_offset_um": round(float(np.hypot(*at.T).max()), 1),
               # moving every cell by A @ (p - origin) + b (H&E level-0 pixels) puts it on its nucleus
               "correction_px": {"origin": [g["x0"], g["y0"]], "A": A.tolist(),
                                 "b": (coef[2] / s.he_um_per_px).tolist()}}
    return {"image": kind, "offset_um": {"dx": dx, "dy": dy, "length": round(float(np.hypot(dx, dy)), 1)},
            "correlation_at_offset": round(float(total.max()), 3), "correlation_at_zero": round(float(total[R, R]), 3),
            "offset_at_search_limit": edge, "window_spread_um": round(spread, 1), "affine_fit": fit, "windows": shifts,
            "note": "offset = where the image's nuclei sit relative to the cell centroids (image minus cells, um; +dx = "
                    "right, +dy = down). Below ~5 um is normal (nucleus vs cell centroid). If the affine fit's residual is much "
                    "smaller than window_spread_um, the offset changes across the section (scale_x/scale_y != 1 or a "
                    "rotation); a max offset above ~8 um suggests redoing the registration"}


# ---------------- report ----------------
def annotation_facts(s, version=None):
    """How this annotation was made: method and its inputs, priors, corrections, regions, history."""
    st = s.state(version)
    sch = {k: v for k, v in s.scheme.items() if k in ("kind", "resolution", "reference", "existing", "n_shared_genes")}
    if s.scheme["kind"] == "markers":
        sch["marker_list"] = s.scheme.get("text")
    from . import constraints as cmod
    return {"sample": {k: v for k, v in s.meta.items() if k in ("sample_id", "description", "n_cells", "n_genes",
                                                                   "image_kind", "he_um_per_px")},
            "session": s.id, "version": st["version"], "scheme": sch, "n_types": len(s.types),
            "priors": [cmod.describe(c) + f'  ("{c.get("text", "")}")' for c in st["constraints"]],
            "corrections": [{"label": c["label"], "n_cells": len(c["cells"]), "text": c.get("text", "")}
                            for c in st.get("corrections", [])],
            "regions": s.region_names(version), "history": s.history()[-15:]}
