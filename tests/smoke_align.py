"""Import the cells alone, add the H&E, align it without a registration file, compare with the registration file.
Runs in a temporary data folder. Needs the SCP2927 example files (paths below)."""
import os, shutil, sys, tempfile, time
import numpy as np
H5AD = os.environ.get("SMOKE_H5AD", "/data2/celltypeagent/share/validation/scp2927/scp2927_2sections.h5ad")
HE = os.environ.get("SMOKE_HE", "/data/kai/endothelial/data/SCP2927/HE/human_slideid-0050423-108731-147328.ome.tif")
REG = os.path.join(os.path.dirname(__file__), "..", "examples", "scp2927", "registration.csv")
home = tempfile.mkdtemp(prefix="celltypeagent_align_")
os.environ["CELLTYPEAGENT_HOME"] = home
from celltypeagent import importer, align
try:
    t = time.time()
    importer.import_sample("t", [H5AD], section_column="donor_id", section="147328", progress=lambda *a: None)
    print("cells only:", importer._meta("t")["image_kind"], f"{time.time() - t:.0f}s")
    print("image:", importer.set_image("t", HE))
    s = align.Sample("t")
    st = align.search(s)
    print("candidates:", [(c["id"], c["orientation_name"], c["um_per_px"], c["score"]) for c in st["candidates"][:3]])
    print("refine:", align.refine(s, "K1")["rounds"][-1])
    print("apply:", align.apply(s, "K1"))
    A = importer.placement("t")
    R = importer.load_registration(REG, "147328")
    xy = align.Sample("t").cells[["x", "y"]].to_numpy()
    d = np.hypot(*((xy @ A[:, :2].T + A[:, 2]) - (xy @ R[:, :2].T + R[:, 2])).T) * 0.69
    print(f"vs the registration file: median {np.median(d):.1f} um (the file itself is ~13 um off the nuclei)")
    print("undo:", importer.undo_placement("t"), "| registration file:", importer.set_registration("t", REG, "147328"))
finally:
    shutil.rmtree(home)
