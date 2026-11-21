"""10x Atera output read as a Xenium output: a synthetic bundle laid out as Atera Onboard Analysis 1.0 writes it
(cell_feature_matrix.zarr.zip: AnnData zarr v2 in a zip, obs centroid_column / centroid_row in um, the barcode as
index, var feature_name / feature_type with control features; csc_cell_feature_matrix.zarr.zip beside it;
experiment.spatial with the pixel size). Found as a folder, field by field (the app's dialog: the matrix in Counts
or in Cells) and as loose files; imported with only the Gene Expression features and the cells' own positions.
ATERA_DIR=<a real Atera output folder> also imports that one (e.g. 10x's tiny_atera_dataset_outs)."""
import json, os, shutil, tempfile, zipfile
os.environ.setdefault("STAGENT_SETTINGS_DIR", tempfile.mkdtemp(prefix="stagent_cfg_"))
home = tempfile.mkdtemp(prefix="stagent_atera_")
os.environ["STAGENT_HOME"] = home
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd, scipy.sparse as sp

rng = np.random.default_rng(0)
root = Path(home) / "atera_out"
root.mkdir()
n, genes = 300, [f"G{i}" for i in range(20)]
ctrl = [f"NegControlProbe_{i}" for i in range(3)] + [f"UnassignedCodeword_{i}" for i in range(2)]
X = sp.csr_matrix(rng.poisson(1.0, (n, len(genes) + len(ctrl))).astype(np.int32))
bc = [f"{''.join(rng.choice(list('abcdefghijklmnop'), 8))}-1" for _ in range(n)]
obs = pd.DataFrame({"barcode": bc, "cell_id": rng.integers(1e9, 9e9, n), "filtered": True,
                    "centroid_column": rng.uniform(0, 500, n).astype(np.float32),
                    "centroid_row": rng.uniform(0, 300, n).astype(np.float32),
                    "azimuth_broad": pd.Categorical(rng.choice(["T", "B"], n))}, index=pd.Index(bc, name="barcode"))
var = pd.DataFrame({"feature_name": genes + ctrl,
                    "feature_type": ["Gene Expression"] * len(genes) + ["Negative Control Probe"] * 3 + ["Unassigned Codeword"] * 2},
                   index=pd.Index([f"ENSG{i:05d}" for i in range(len(genes) + len(ctrl))], name="feature_id"))
ad.settings.zarr_write_format = 2  # (as Atera writes it)
a = ad.AnnData(X=X, obs=obs.drop(columns=["barcode"]), var=var)
a.obs["barcode"] = bc
tmp = Path(home) / "z"
a.write_zarr(tmp)


def zipdir(d, out):
    with zipfile.ZipFile(out, "w", zipfile.ZIP_STORED) as z:
        for p in sorted(d.rglob("*")):
            if p.is_file():
                z.write(p, p.relative_to(d).as_posix())


zipdir(tmp, root / "cell_feature_matrix.zarr.zip")
zipdir(tmp, root / "csc_cell_feature_matrix.zarr.zip")
shutil.rmtree(tmp)
(root / "experiment.spatial").write_text(json.dumps({"pixel_size": 0.21192053, "panel_wta": {"panel_name": "WTA"},
                                                     "chemistry_version": "Atera v1"}))

from stagent import importer

# a folder
s = importer.inspect_xenium_many([str(root)])["slices"]
assert len(s) == 1 and s[0]["n_cells"] == n and s[0]["platform"] == "Atera", s
assert abs(s[0]["pixel_size_um"] - 0.21192053) < 1e-9 and s[0]["panel_name"] == "WTA", s
# field by field: the matrix in Counts (Cells empty), or in Cells
for fields in ({"cells": [], "counts": [str(root / "cell_feature_matrix.zarr.zip")], "experiment": [str(root / "experiment.spatial")]},
               {"cells": [str(root / "cell_feature_matrix.zarr.zip")], "counts": []}):
    s = importer.inspect_xenium_many([{"fields": fields}])["slices"]
    assert len(s) == 1 and s[0]["n_cells"] == n, s
# loose files: csc_ copy not taken for a second slice
s = importer.xenium_sources([str(root / "cell_feature_matrix.zarr.zip"), str(root / "csc_cell_feature_matrix.zarr.zip"),
                             str(root / "experiment.spatial")])
assert len(s) == 1 and "experiment" in s[0]["files"], s
# an unzipped bundle
unz = Path(home) / "atera_unzipped"
unz.mkdir()
with zipfile.ZipFile(root / "cell_feature_matrix.zarr.zip") as z:
    z.extractall(unz / "cell_feature_matrix.zarr")
assert importer.inspect_xenium_many([str(unz)])["slices"][0]["n_cells"] == n

meta = importer.import_xenium("atera1", str(root), progress=lambda *_: None)
assert meta["n_cells"] == n and meta["n_genes"] == len(genes), meta
assert meta["source"]["kind"] == "xenium" and meta["source"]["platform"] == "Atera"
assert meta["image_kind"] == "none"  # (no morphology image in this bundle)
c = pd.read_parquet(Path(home) / "samples" / "atera1" / "cells.parquet")
assert list(c.cell_id) == bc
assert np.allclose(c.x, obs.centroid_column.to_numpy(float), atol=1e-4) and np.allclose(c.y, obs.centroid_row.to_numpy(float), atol=1e-4)
assert (c.total_counts.to_numpy() == np.asarray(X[:, :len(genes)].sum(1)).ravel()).all()
assert open(Path(home) / "samples" / "atera1" / "genes.txt").read().split() == genes
# the notebook's replay reads the same files
src = meta["source"]
bc2, g2, X2 = importer._xenium_matrix({"key": src["key"], "folder": src["folder"], "files": src["files"]})
assert list(g2) == genes and X2.shape == (n, len(genes))
from stagent import protein
assert protein.read(src) is None  # no protein features: no stains

if os.environ.get("ATERA_DIR"):  # a real output
    m = importer.import_xenium("atera_real", os.environ["ATERA_DIR"], progress=print)
    print("real:", m["n_cells"], "cells,", m["n_genes"], "genes, image", m["image_kind"], m["he_um_per_px"])
shutil.rmtree(home, ignore_errors=True)
print("smoke_atera ok")
