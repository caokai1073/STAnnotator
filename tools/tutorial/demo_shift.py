"""A demo slice for the Local registration chapter: a window of 147328 with one patch of cells shifted on purpose.

Both demo sections sit on their H&E (their only offsets are the cell-centroid one of the villus epithelium at the
luminal edge, not a drift), so a local registration has nothing to fix there. This takes a window of the 147328 crop
away from that edge (crypts, lamina propria, submucosa), moves the cells in one patch by a smooth bump (as a section
stretched between the run and the staining would be), and writes the registration the agent found for 147328, so the
slice imports already placed: right everywhere but in the patch. The tutorial says the patch was moved on purpose.

    .venv/bin/python tools/tutorial/demo_shift.py work/tut/data work/tut/home/samples/147328/meta.json

Writes <data>/scp2927_147328_shifted.h5ad, _shifted_HE.ome.tif and registration_shifted.csv.
"""
import json, sys
import numpy as np, anndata as ad, tifffile, zarr

D, META = sys.argv[1].rstrip("/") + "/", sys.argv[2]
X0, X1, Y0, Y1 = 5600, 7600, 8300, 10200  # the window, um (cell coordinates): crypts down to the submucosa
CENTER, SIGMA, SHIFT = (6600, 9500), 260.0, (14.0, -12.0)  # the bump: center and width (um), peak shift (um)
M = 8  # H&E margin around the cells (px): nuclei on the image beyond the cells would pull the edge windows off

A = np.array(json.load(open(META))["placement_affine"])  # um -> px of the H&E crop: the agent's alignment of 147328
a = ad.read_h5ad(D + "scp2927_147328_crop.h5ad")
xy = a.obsm["spatial"]
b = a[(xy[:, 0] >= X0) & (xy[:, 0] < X1) & (xy[:, 1] >= Y0) & (xy[:, 1] < Y1)].copy()
xy = np.asarray(b.obsm["spatial"], float)
w = np.exp(-((xy[:, 0] - CENTER[0]) ** 2 + (xy[:, 1] - CENTER[1]) ** 2) / (2 * SIGMA ** 2))
moved = xy + w[:, None] * np.array(SHIFT)
px = xy @ A[:, :2].T + A[:, 2]  # the true positions on the image, for the crop
b.obsm["spatial"] = moved.astype(np.float32)
b.obs["x_centroid"], b.obs["y_centroid"] = moved[:, 0], moved[:, 1]
print(f"{b.n_obs} cells; {int((w * np.hypot(*SHIFT) > 5).sum())} moved more than 5 um, at most {np.hypot(*SHIFT):.1f} um;"
      f" bump center at image px {np.round(np.array(CENTER) @ A[:, :2].T + A[:, 2])}")

cx0, cy0 = int(max(0, px[:, 0].min() - M)), int(max(0, px[:, 1].min() - M))
cx1, cy1 = int(px[:, 0].max() + M), int(px[:, 1].max() + M)
store = tifffile.imread(D + "scp2927_147328_HE_crop.ome.tif", aszarr=True, level=0)
crop = np.asarray(zarr.open(store, mode="r")[cy0:cy1, cx0:cx1, :]); store.close()
print("H&E", crop.shape, "offset", cx0, cy0)
with tifffile.TiffWriter(D + "scp2927_147328_shifted_HE.ome.tif", ome=True) as tw:
    meta = {"axes": "YXS", "PhysicalSizeX": 0.69, "PhysicalSizeXUnit": "µm", "PhysicalSizeY": 0.69, "PhysicalSizeYUnit": "µm"}
    opts = dict(photometric="rgb", tile=(512, 512), compression="jpeg", compressionargs={"level": 88})
    tw.write(crop, subifds=3, metadata=meta, **opts)
    lv = crop
    for _ in range(3):
        lv = lv[::2, ::2]
        tw.write(np.ascontiguousarray(lv), subfiletype=1, **opts)
b.write_h5ad(D + "scp2927_147328_shifted.h5ad", compression="gzip")
with open(D + "registration_shifted.csv", "w") as f:
    f.write("a11,a12,a21,a22,tx,ty\n")
    f.write(f"{A[0,0]},{A[0,1]},{A[1,0]},{A[1,1]},{A[0,2] - cx0},{A[1,2] - cy0}\n")
