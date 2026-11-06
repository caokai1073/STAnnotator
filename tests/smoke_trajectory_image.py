"""Smoke test of v0.2.45's two additions, on synthetic data (no dataset needed):
1. trajectory.py: a spatial axis, diffusion pseudotime (expression only and mixed with space) and the genes along it,
   on cells along a gradient whose genes rise, fall or peak in the middle;
2. an H&E-only sample (importer.import_image): regions drawn on an image without cells, written as GeoJSON and read
   back, and handed to a sample of cells on the same image;
3. outlines read from the H&E itself (hetissue.py) and the tissue priors document.
Run: python tests/smoke_trajectory_image.py  (in a temporary data folder, removed at the end)."""
import os, tempfile, time
from pathlib import Path

HOME = Path(tempfile.mkdtemp(prefix="stagent_smoke_"))
os.environ["STAGENT_HOME"] = str(HOME)
import numpy as np, scipy.sparse as sp  # noqa: E402
from stagent import trajectory as tr, importer, regions  # noqa: E402

t0 = time.time()
rng = np.random.default_rng(1)

# ---- 1. trajectory on a synthetic gradient: 3,000 cells along x (0..3000 um); s = x / 3000 is the truth
n = 3000
x = np.sort(rng.uniform(0, 3000, n))
xy = np.c_[x, rng.uniform(0, 200, n)]
s_true = x / 3000
rate = np.c_[20 * s_true, 20 * (1 - s_true), 20 * np.exp(-((s_true - 0.5) / 0.12) ** 2), np.full((n, 5), 3.0)]
counts = sp.csr_matrix(rng.poisson(rate + 0.2).astype(np.float32))
genes = ["RISE", "FALL", "MID"] + [f"FLAT{i}" for i in range(5)]
from stagent import model  # noqa: E402
Xn = model.lognorm(counts)
one = np.zeros(n, dtype=object)
start, end = s_true < 0.05, s_true > 0.95

t, raw = tr.spatial_axis(xy, one, start, end)
r_sp = np.corrcoef(t, s_true)[0, 1]
d, _ = tr.spatial_axis(xy, one, start)
assert r_sp > 0.95, r_sp
assert raw is not None and np.nanmax(d) <= 1.0
trends, tab = tr.gene_trends(Xn, t, genes)
print(f"spatial axis: r with the truth {r_sp:.3f}; rising {trends['rising'][0]['gene']}, falling "
      f"{trends['falling'][0]['gene']}, transient {[g['gene'] for g in trends['transient'][:2]]}")
assert trends["rising"][0]["gene"] == "RISE" and trends["falling"][0]["gene"] == "FALL"
assert "MID" in [g["gene"] for g in trends["transient"]]

from stagent import markers  # noqa: E402
Z = markers.pca(Xn, 5)[:, ::-1].copy()
for alpha in (0.0, 0.5):
    pt, info = tr.pseudotime(Z, xy, one, start, alpha=alpha)
    rho = np.corrcoef(np.argsort(np.argsort(pt)), np.argsort(np.argsort(s_true)))[0, 1]
    print(f"pseudotime (space {alpha}): rank correlation with the truth {rho:.3f}; unreached {info['n_unreached']}")
    assert rho > (0.8 if alpha == 0 else 0.9), rho  # noisy counts of 8 genes; the tissue's order helps
cache = {}
pt, _ = tr.pseudotime(Z, xy, one, start, max_cells=1000, cache=cache)  # a sample of 1,000, the others placed
rho = np.corrcoef(np.argsort(np.argsort(pt)), np.argsort(np.argsort(s_true)))[0, 1]
print(f"pseudotime on a sample of 1,000: rank correlation {rho:.3f}; unplaced {int(np.isnan(pt).sum())}")
assert rho > 0.75 and not np.isnan(pt).any(), rho
n_keys = len(cache)
pt2, _ = tr.pseudotime(Z, xy, one, end, max_cells=1000, cache=cache)  # another root: the graph is reused
assert len(cache) == n_keys and np.corrcoef(pt, pt2)[0, 1] < -0.5, (len(cache), n_keys)
comp = tr.composition(t, np.where(s_true < 0.5, "A", "B"))
assert comp["shares"][0][0] > 0.9 and comp["shares"][-1][1] > 0.9
assert tr.legend({"unit": "um", "max_um": 1000})[0] == "0–100 µm"

# ---- 2. H&E only: an image without cells, regions on it, GeoJSON out and in, handed to cells on the same image
from PIL import Image  # noqa: E402
img = HOME / "he.png"
Image.fromarray(rng.integers(150, 255, (600, 800, 3), dtype=np.uint8)).save(img)
m = importer.import_image("he_only", img, um_per_px=0.5, progress=lambda *_: None)
assert m["image_only"] and m["n_cells"] == 0 and abs(m["he_um_per_px"] - 0.5) < 1e-9
m2 = importer.import_image("he_unknown", img, progress=lambda *_: None)
assert m2["um_per_px_known"] is False
from stagent.session import Session  # noqa: E402
s = Session.create("he_only", "none", method="none")
s.add_region("mucosa", [[10, 10], [400, 10], [400, 300], [10, 300]])
s.add_region("submucosa", [[400, 10], [790, 10], [790, 590], [400, 590]])
info = s.regions_info()
assert [r["name"] for r in info["regions"]] == ["mucosa", "submucosa"], info["regions"]
text, gi = regions.to_geojson(s.state()["regions"], {"mucosa": "#ff8800"}, (800, 600))
pieces, skipped = regions.parse_geojson(text)
assert [p["name"] for p in pieces] == ["mucosa", "submucosa"] and pieces[0]["color"] == "#ff8800" and not skipped
try:
    Session.create("he_only", "none", method="markers", markers_text="A: X+")
    raise AssertionError("annotating an image only should fail")
except ValueError:
    pass

# a sample of cells on the same image (pixels = um / 0.5): the regions of he_only come along
cells_px = rng.uniform([20, 20], [780, 580], (500, 2))
ids = np.array([f"c{i}" for i in range(500)], dtype=object)
X = sp.csr_matrix(rng.poisson(2, (500, 4)).astype(np.float32))
A = np.array([[2.0, 0, 0], [0, 2.0, 0]])
tif = m["he_path"]  # the PNG as the pyramidal OME-TIFF the import made of it
importer._write_bundle("cells", ids, cells_px * 0.5, X, ["G1", "G2", "G3", "G4"], tif, A, None, "", {"kind": "test"},
                       lambda *_: None, "H&E")
gj, _ = importer.image_sample_regions("he_only", tif)
per = importer.attach_regions("cells", text=gj, source="sample he_only")
inside_m = ((cells_px[:, 0] < 400) & (cells_px[:, 1] < 300) & (cells_px[:, 1] > 10)).sum()
print(f"H&E only: regions {[p['name'] for p in pieces]}; handed to 500 cells: {per} (mucosa expected {inside_m})")
assert per["mucosa"] == inside_m
ph = Session.create("cells", "none", method="none")
assert {r["name"] for r in ph.regions_info()["regions"]} == {"mucosa", "submucosa"}
# ---- 3. outlines from the H&E itself: a synthetic "H&E" with a purple band (mucosa) over a pink one (muscle)
from stagent import hetissue  # noqa: E402
he = np.full((1200, 1600, 3), 235, np.uint8)
he[200:500, 200:1400] = (120, 60, 150)      # hematoxylin-rich
he[500:1000, 200:1400] = (230, 90, 140)     # eosin-rich
he += rng.integers(0, 12, he.shape, dtype=np.uint8)
Image.fromarray(he).save(HOME / "he2.png")
importer.import_image("he2", HOME / "he2.png", um_per_px=2.0, progress=lambda *_: None)
seg = hetissue.segment("he2", n_classes=2)
blue = min(seg["classes"], key=lambda c: -c["hematoxylin"])["class"]
polys = [pc["outer"] for pc in hetissue.outline("he2", [blue])]
xs = np.array([x for p in polys for x, _ in p]); ys = np.array([y for p in polys for _, y in p])
print(f"H&E classes: {[(c['class'], c['hematoxylin'], c['eosin']) for c in seg['classes']]}; outline box "
      f"{xs.min():.0f}-{xs.max():.0f} x {ys.min():.0f}-{ys.max():.0f}")
assert len(polys) == 1 and abs(ys.min() - 200) < 30 and abs(ys.max() - 500) < 30, (len(polys), ys.min(), ys.max())
assert "Compartments to draw" in hetissue.tissue_priors("colon")["text"]

import shutil  # noqa: E402
shutil.rmtree(HOME, ignore_errors=True)
print(f"ok ({time.time() - t0:.1f} s)")
