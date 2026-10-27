"""R files without R: two small objects made in R from 1,500 cells of the example MERFISH slice (tests/data:
a Seurat v5 object with an FOV of centroids and the counts split into two layers, and an xz-compressed
SpatialExperiment) must give the same counts, cells and coordinates as the h5ad, and import as a sample."""
import os, shutil, tempfile
import numpy as np
HERE = os.path.dirname(__file__)
home = tempfile.mkdtemp(prefix="stagent_rds_")
os.environ["STAGENT_HOME"] = home
from stagent import importer, rds
try:
    ref = importer.read_h5ad(os.path.join(HERE, "..", "examples", "cadinu_merfish", "cadinu_merfish_D0_m3a_1_slice_2.h5ad"))
    pos = {c: i for i, c in enumerate(ref["cell_id"])}
    for f in ("seurat_v5_small.rds", "spe_small.rds"):
        t = rds.table(os.path.join(HERE, "data", f))
        rows = [pos[c] for c in t["cell_id"]]
        assert t["genes"] == ref["genes"], f
        assert abs(t["X"] - ref["X"][rows]).max() == 0, f
        assert np.allclose(t["xy"], ref["xy"][rows]), f
        print(f, t["kind"], t["X"].shape, "obs:", list(t["obs"].columns))
    m = importer.import_sample("small", [os.path.join(HERE, "data", "seurat_v5_small.rds")], section_column="half",
                               section="left", progress=lambda m: None)
    print("imported the left half:", m["n_cells"], "cells")
finally:
    shutil.rmtree(home)
