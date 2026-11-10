"""Moving a torn piece of tissue onto its place on the image (the H&E tore or folded after the run).

A smooth correction (regcheck) cannot follow a piece that tore off, turned or folded over: the user outlines the
piece's cells as a region, and fit() searches the rigid move (rotation, shift, and mirrored: a piece that folded over
shows its back) that puts them on the image, within MAX_SHIFT_UM of where they are. The score is the normalized
cross-correlation, over the piece's box, of the cells' density with the nuclei (hematoxylin of the H&E, or DAPI), and
of the cells' footprint with the tissue (TISSUE: its weight; a folded piece is often out of focus, its nuclei faint).
The other cells nearby are taken off both first: they already explain the nuclei and the tissue they sit on, so a
piece is not pulled onto tissue whose cells are there. Every angle is tried at COARSE_UM per pixel (in steps that
move the piece's edge by about a pixel), then the best few at FINE_UM, scored mostly by the nuclei (FINE_TISSUE). The fit counts when it is clearly better
than the cells where they are, reaches MIN_CORR, and is unique (the best other placement, one that puts the cells
APART_UM away, scores at most UNIQUE of it).

Fits are kept in the sample's folder (pieces.json) with the placement they were computed on; the agent looks at the
figure and proposes the move (mcp_server.move_piece), which importer.move_piece applies on the user's Accept.
"""
import json, time
import numpy as np, pandas as pd
from scipy import ndimage
from . import config, render
from .osutil import dump_json

MAX_SHIFT_UM = 300.0
COARSE_UM, FINE_UM = 4.0, 2.0  # image read at these um per pixel
SIGMA_UM = 6.0                 # the cells' density and the nuclei are smoothed by this much
MAX_RADIUS_UM = 1500.0         # a larger outline is not one piece
MIN_CELLS = 30
UNIQUE, LIKELY = 0.97, 0.9     # the best other placement scores at most UNIQUE of the fit (crypts repeat: often
                               # above LIKELY when the fit is right)
MIN_CORR = 0.45                # a fit below this does not match the tissue
APART_UM = 30.0                # another placement: one that puts the cells this far from the best (median)
TISSUE, FINE_TISSUE = 0.5, 0.2  # weight of the tissue term (the nuclei get the rest), in the search and the refinement
MIN_GAIN = 0.1                 # the fit must beat the cells where they are by this much (correlation)


def _file(sample_id):
    return config.SAMPLES / sample_id / "pieces.json"


def load(sample_id):
    try:
        return json.load(open(_file(sample_id)))
    except (OSError, ValueError):
        return {"pieces": []}


def _rot(deg, mirrored):
    a = np.deg2rad(deg)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    return R @ np.diag([-1.0, 1.0]) if mirrored else R


def _cover(D):
    """Soft footprint (0..1) of a cell density: 1 where cells are as dense as in sparse tissue."""
    pos = D[D > 1e-3]
    return np.clip(D / max(np.percentile(pos, 30) if len(pos) else 1, 1e-6), 0, 1)


class _Image:
    """Around the piece at g um per pixel, zero mean: the nuclei with the other cells' density taken off (I), and the
    tissue no other cell covers (U)."""

    def __init__(self, sl, kind, center, half_um, upp, g, others):
        from .insights import _nuclei
        from .align import _tissue
        self.g, self.upp, self.c = g, upp, np.asarray(center, float)
        self.n = n = int(np.ceil(2 * half_um / g))
        h = n * g / upp / 2
        img = np.asarray(sl.read(center[0] - h, center[1] - h, center[0] + h, center[1] + h, n, n)).astype(np.float32)
        nuc = ndimage.gaussian_filter(np.clip(_nuclei(img, kind), 0, None), SIGMA_UM / g / 2)
        nuc = nuc / (nuc.std() + 1e-9)
        U = _tissue(img, kind)
        if len(others):
            d = self.density(others)
            # take off what the other cells explain: their density scaled to the nuclei where they are
            k = float((nuc * d).sum() / max((d * d).sum(), 1e-9))
            nuc = nuc - np.clip(k, 0, None) * d
            U = U * (1 - _cover(d))
        self.I = nuc - nuc.mean()
        self.U = (U - U.mean()) / (U.std() + 1e-9)

    def grid(self, p):
        """Image pixels (level 0) -> this grid's pixels (float)."""
        return (p - self.c) * self.upp / self.g + self.n / 2

    def density(self, p):
        n, u = self.n, self.grid(p)
        D = np.zeros((n, n), np.float32)
        ok = (u[:, 0] >= 0) & (u[:, 0] < n) & (u[:, 1] >= 0) & (u[:, 1] < n)
        np.add.at(D, (u[ok, 1].astype(int), u[ok, 0].astype(int)), 1)
        return ndimage.gaussian_filter(D, SIGMA_UM / self.g)


def _boxsum(A, m):
    """Sum of A over every m x m box that fits whole; out[i, j]: the box with its corner at (i, j)."""
    S = np.pad(A, ((1, 0), (1, 0))).cumsum(0).cumsum(1)
    return S[m:, m:] - S[:-m, m:] - S[m:, :-m] + S[:-m, :-m]


def _ncc(I, T, S1, S2):
    """Normalized cross-correlation of template T (m x m, zero mean) at every place in I (n x n) where it fits whole
    (S1, S2: _boxsum of I and I^2); out[i, j]: T's corner at (i, j)."""
    n, m = I.shape[0], T.shape[0]
    num = np.fft.irfft2(np.fft.rfft2(I) * np.conj(np.fft.rfft2(T, (n, n))), (n, n))[:n - m + 1, :n - m + 1]
    var = np.maximum(S2 - S1 ** 2 / (m * m), 0.25 * (I * I).mean() * m * m)  # a blank box does not score by chance
    return num / (np.sqrt(var) * np.sqrt((T * T).sum()) + 1e-9)


def _template(q, g, upp, r_px):
    """Density and footprint of the piece's cells q (level-0 pixels, centered on 0), on an m x m grid centered on
    0, zero mean."""
    m = int(np.ceil(2 * r_px * upp / g)) + 2 * int(np.ceil(2 * SIGMA_UM / g)) + 1
    u = q * upp / g + m / 2
    D = np.zeros((m, m), np.float32)
    np.add.at(D, (u[:, 1].astype(int).clip(0, m - 1), u[:, 0].astype(int).clip(0, m - 1)), 1)
    D = ndimage.gaussian_filter(D, SIGMA_UM / g)
    K = _cover(D)
    return D - D.mean(), K - K.mean(), m


def _search(im, P, c, angles, mirrors, max_shift_um, r_px, tissue=TISSUE):
    """[(score, angle, mirrored, shift_px)] best place per orientation."""
    out, box = [], {}
    for mi in mirrors:
        for a in angles:
            T, K, m = _template((P - c) @ _rot(a, mi).T, im.g, im.upp, r_px)
            if m not in box:
                box[m] = [(_boxsum(X, m), _boxsum(X * X, m)) for X in (im.I, im.U)]
            C = (1 - tissue) * _ncc(im.I, T, *box[m][0]) + tissue * _ncc(im.U, K, *box[m][1])
            k = C.shape[0]
            # T's center at im pixel (j + m/2, i + m/2); the piece's center c sits at im pixel n/2
            off = (np.arange(k) + m / 2 - im.n / 2) * im.g  # um
            far = np.hypot(off[None, :], off[:, None]) > max_shift_um
            C = np.where(far, -np.inf, C)
            i, j = np.unravel_index(np.argmax(C), C.shape)
            out.append((float(C[i, j]), float(a), bool(mi), np.array([off[j], off[i]]) / im.upp))
    return out


def fit(s, region, max_shift_um=MAX_SHIFT_UM, mirror=True, progress=None):
    """Fit the rigid move of the cells of `region` (in the open annotation) onto the image; stored as P1, P2, ..."""
    from .importer import image_of
    from . import regcheck
    say = progress or (lambda f, m: None)
    image = image_of(s.sample_id)
    kind = s.meta.get("image_kind", "H&E")
    if not image or kind == "none":
        raise ValueError("no image to fit the piece on (none, or unloaded: Load H&E in the Data tab)")
    ids = s._cells_for(region, None, None)
    if len(ids) < MIN_CELLS:
        raise ValueError(f"region {region!r} holds {len(ids)} cells: outline a piece of {MIN_CELLS} cells or more")
    upp = s.he_um_per_px
    allp = s.cells[["px", "py"]].to_numpy()
    inp = np.isin(s.cells.cell.to_numpy(), ids)
    P = allp[inp]
    c = P.mean(0)
    r_px = float(np.percentile(np.hypot(*(P - c).T), 99)) + SIGMA_UM / upp
    if r_px * upp > MAX_RADIUS_UM:
        raise ValueError(f"region {region!r} spans {2 * r_px * upp:.0f} um: outline one torn piece at a time "
                         f"(up to {2 * MAX_RADIUS_UM:.0f} um)")
    sl = render.slide(image)
    half = r_px * upp + max_shift_um + 2 * SIGMA_UM
    near = np.abs(allp - c).max(1) < half / upp
    others = allp[near & ~inp]
    say(0.05, f"region {region}: {len(P):,} cells")
    im = _Image(sl, kind, c, half, upp, COARSE_UM, others)
    step = np.rad2deg(1.0 * COARSE_UM / max(r_px * upp, 10))
    angles = np.arange(-180, 180, max(step, 0.5))
    mirrors = [False, True] if mirror else [False]
    found = []
    for k, mi in enumerate(mirrors):
        found += _search(im, P, c, angles, [mi], max_shift_um, r_px)
        say(0.1 + 0.6 * (k + 1) / len(mirrors), "trying every angle" + (", mirrored" if mi else ""))
    found.sort(key=lambda x: -x[0])
    best = found[0]

    sub = P[np.random.default_rng(0).choice(len(P), min(len(P), 500), replace=False)]

    def apart(a, b):  # the two placements put the cells APART_UM apart (median)
        qa, qb = ((sub - c) @ _rot(x[1], x[2]).T + x[3] for x in (a, b))
        return np.median(np.hypot(*(qa - qb).T)) * upp > APART_UM
    second = next((f for f in found[1:] if apart(f, best)), None)
    # fine: around the best few orientations, at FINE_UM, within a few coarse pixels
    fine_im = _Image(sl, kind, c, r_px * upp + max_shift_um + 2 * SIGMA_UM, upp, FINE_UM, others)
    cands = [best] + [f for f in found[1:6] if not apart(f, best)][:2]
    fine = []
    for sc, a, mi, t in cands:
        fa = a + np.arange(-1, 1.01, 0.25) * max(step, 0.5)
        for x in _search(fine_im, P, c, fa, [mi], max_shift_um, r_px, FINE_TISSUE):
            if np.hypot(*(x[3] - t)) * upp <= 3 * COARSE_UM:
                fine.append(x)
    say(0.9, "refining")
    sc, a, mi, t = max(fine, key=lambda x: x[0]) if fine else best
    now = _search(fine_im, P, c, [0.0], [False], 3 * FINE_UM, r_px, FINE_TISSUE)[0][0]  # the cells where they are
    M = _rot(a, mi)
    Q = (P - c) @ M.T + c + t
    mv = np.hypot(*(Q - P).T) * upp
    pieces = [p for p in load(s.sample_id)["pieces"] if p.get("key") == regcheck._key(s.sample_id, image)]
    pid = f"P{len(pieces) + 1}"
    rec = {"id": pid, "region": region, "session_id": s.id, "n_cells": int(len(P)), "key": regcheck._key(s.sample_id, image),
           "angle_deg": round(a, 2), "mirrored": mi, "shift_um": [round(float(v) * upp, 1) for v in t],
           "moved_um_median": round(float(np.median(mv)), 1), "moved_um_max": round(float(mv.max()), 1),
           "corr_now": round(now, 3), "corr": round(sc, 3),
           "second": round(second[0] / max(best[0], 1e-9), 3) if second else 0.0,
           "second_at": {"angle_deg": round(second[1], 1), "mirrored": second[2],
                         "shift_um": [round(float(v) * upp, 1) for v in second[3]]} if second else None,
           "center_px": c.tolist(), "M": M.tolist(), "t_px": t.tolist(), "radius_um": round(r_px * upp, 1),
           "max_shift_um": max_shift_um, "time": time.strftime("%Y-%m-%dT%H:%M:%S"), "cells": [int(x) for x in ids]}
    rec["verdict"] = ("no better than where the cells are" if sc - now < MIN_GAIN else
                      "weak: the cells do not match the tissue anywhere near" if sc < MIN_CORR else
                      "not unique: another placement fits as well" if rec["second"] > UNIQUE else
                      "fits" + (f"; another placement scores {rec['second']:.0%} of it (repeating structures): check the "
                                "figure" if rec["second"] > LIKELY else ""))
    if max(abs(v) for v in rec["shift_um"]) > 0.95 * max_shift_um:
        rec["verdict"] += "; at the edge of the search (try a larger max_shift_um)"
    dump_json({"pieces": pieces + [rec]}, _file(s.sample_id))
    say(1.0, "done")
    return rec


def get(sample_id, pid):
    from .importer import image_of
    from . import regcheck
    p = next((x for x in load(sample_id)["pieces"] if x["id"] == pid.upper()), None)
    if p is None:
        raise ValueError(f"no fitted piece {pid}: fit it first (fit_piece)")
    if p["key"] != regcheck._key(sample_id, image_of(sample_id)):
        raise ValueError(f"the cells moved since piece {pid} was fitted: fit it again")
    return p


def public(p):
    """A fit without its cell ids (for the agent)."""
    return {k: v for k, v in p.items() if k not in ("cells", "key", "M", "t_px")}


def figure(s, p):
    """Left: the H&E around the piece and where it goes (the piece's cells now cyan, moved yellow, other cells gray).
    Right: the H&E where the piece goes, with the moved cells (yellow should sit on the nuclei of the torn piece)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .importer import image_of
    sl = render.slide(image_of(s.sample_id))
    upp = s.he_um_per_px
    allp = s.cells[["px", "py"]].to_numpy()
    inp = np.isin(s.cells.cell.to_numpy(), p["cells"])
    P = allp[inp]
    c = np.asarray(p["center_px"])
    Q = (P - c) @ np.asarray(p["M"]).T + c + np.asarray(p["t_px"])
    fig, ax = plt.subplots(1, 2, figsize=(15, 7.6))
    both = np.vstack([P, Q])
    lo, hi = both.min(0), both.max(0)
    mid, half = (lo + hi) / 2, max(hi - lo) / 2 * 1.15 + 20 / upp
    for a, (ctr, h, title) in zip(ax, [(mid, half, "now (cyan) and moved (yellow)"),
                                        ((Q.min(0) + Q.max(0)) / 2, max(Q.max(0) - Q.min(0)) / 2 * 1.1 + 15 / upp,
                                         "moved (yellow): on the nuclei of the piece?")]):
        a.imshow(np.asarray(sl.read(ctr[0] - h, ctr[1] - h, ctr[0] + h, ctr[1] + h, 900, 900)),
                 extent=(ctr[0] - h, ctr[0] + h, ctr[1] + h, ctr[1] - h))
        m = (np.abs(allp - ctr).max(1) < h) & ~inp
        a.scatter(*allp[m].T, s=1.5, c="#9a9a9a", alpha=0.5, linewidths=0, rasterized=True)
        a.set_xlim(ctr[0] - h, ctr[0] + h); a.set_ylim(ctr[1] + h, ctr[1] - h); a.set_xticks([]); a.set_yticks([])
        a.set_title(f"{title}; {2 * h * upp:.0f} um wide", fontsize=9)
    ax[0].scatter(*P.T, s=3, c="#00c8ff", linewidths=0, alpha=0.8, rasterized=True)
    ax[0].scatter(*Q.T, s=3, c="#ffd400", linewidths=0, alpha=0.8, rasterized=True)
    ax[1].scatter(*Q.T, s=5, c="#ffd400", linewidths=0, alpha=0.9, rasterized=True)
    fig.suptitle(f"{p['id']}: region {p['region']}, {p['n_cells']:,} cells; rotated {p['angle_deg']:.1f} deg"
                 + (", mirrored (folded over)" if p["mirrored"] else "")
                 + f", shifted {np.hypot(*p['shift_um']):.0f} um\ncorrelation with the image {p['corr_now']:.2f} -> "
                 f"{p['corr']:.2f}, best other placement {p['second']:.0%} of it: {p['verdict']}", fontsize=10)
    out = config.SAMPLES / s.sample_id / "figures" / f"piece_{p['id']}.jpg"
    out.parent.mkdir(exist_ok=True)
    plt.tight_layout(rect=(0, 0, 1, 0.93))
    plt.savefig(out, dpi=90, pil_kwargs={"quality": 85}); plt.close(fig)
    return str(out)
