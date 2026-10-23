"""Lab server only: build SCP2927 per-sample bundles from the analysis pipeline outputs.

Usage: python -m celltypeagent.prepare
Each bundle (work/samples/<sample_id>/) holds:
  cells.parquet   one row per cell: global index, coordinates (um and H&E level-0 px), lineage,
                  reference label (dataset annotation, shown for comparison only), pathologist region,
                  depth / marker neighborhoods, H&E domain, total counts
  base_logp.npy   float16 (n_cells, n_types) calibrated base log-probabilities. For SCP2927 these are the
                  leave-one-slide-out predictions of scripts/02_key_test.py (model never saw this slide).
  tiles.parquet   H&E tiles of this sample with their domain (for domain example sheets)
  counts.npz      raw counts (CSR, rows in cells.parquet order) for re-prediction after model updates
  meta.json       sample metadata, type list, H&E path, source folder
Also builds the reference pool for the expression model (celltypeagent.model.build_reference).
"""
import glob, json
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from . import config, model

# Lab-server inputs for SCP2927 (Pokatayev et al., Nature 2026); other users import data with `celltypeagent import`
SOURCE_DIR = Path("/data/kai/endothelial/data/SCP2927")
HE_DIR = SOURCE_DIR / "HE"
REGISTRATION = Path("/data/kai/endothelial/hne_veins2/registration/transforms_final_lvl0.csv")


def main():
    W = config.WORK
    obs = pd.read_parquet(W / "obs_all.parquet")
    types = list(np.load(W / "types.npy", allow_pickle=True))
    lp = np.load(W / "logprobs_full.npy", mmap_mode="r")
    X = sp.load_npz(W / "counts.npz").tocsr()
    dep = pd.read_parquet(W / "depth.parquet")
    dom = pd.read_parquet(W / "domains/regions_inferred.parquet")
    tiles = pd.read_parquet(W / "domains/tile_domains.parquet")
    aff = pd.read_csv(REGISTRATION).set_index("donor")
    lineage_of_type = obs.groupby("celltype").file.first().to_dict()
    domain_names = json.load(open(W / "domains/domain_names_auto.json"))
    config.SAMPLES.mkdir(parents=True, exist_ok=True)
    json.dump({"types": types, "lineage_of_type": lineage_of_type, "domain_names_auto": domain_names},
              open(config.SAMPLES / "dataset.json", "w"), indent=1)

    for donor, g in obs.groupby("donor_id"):
        idx = g.index.values
        a = aff.loc[donor]
        x, y = g.x_centroid.values, g.y_centroid.values
        cells = pd.DataFrame({
            "cell": idx, "x": x, "y": y,
            "px": a.a11 * x + a.a12 * y + a.tx, "py": a.a21 * x + a.a22 * y + a.ty,
            "lineage": g.file.values, "reference_label": g.celltype.values,
            "pathologist_region": g.histopathology.values, "total_counts": g.total_counts.values,
            "depth_um": dep.depth_um.values[idx], "epi_nbr": dep.epi_nbr.values[idx],
            "smc_nbr": dep.smc_nbr.values[idx], "lym_nbr": dep.lym_nbr.values[idx],
            "domain": dom.domain.values[idx],
        })
        sid = str(donor)
        d = config.SAMPLES / sid
        d.mkdir(exist_ok=True)
        cells.to_parquet(d / "cells.parquet")
        np.save(d / "base_logp.npy", np.asarray(lp[idx]))
        sp.save_npz(d / "counts.npz", X[idx])
        tiles[tiles.donor_id == donor][["gx", "gy", "px", "py", "n_cells", "domain"]].to_parquet(d / "tiles.parquet")
        slide = int(g.slide_id.iloc[0])
        he = glob.glob(str(HE_DIR / f"human_slideid-{slide:07d}-*.tif"))
        first = g.iloc[0]
        meta = {"sample_id": sid, "slide_id": slide, "disease": first.disease_diagnosis, "location": first.location,
                "inflammation": first.inflammation, "ulceration": first.ulceration, "fibrosis": first.fibrosis,
                "n_cells": int(len(g)), "he_path": he[0] if he else None,
                "source_dir": str(SOURCE_DIR), "bundle_dir": str(d),
                "base_model": "logistic regression on 480 genes, temperature-scaled, trained on the other 7 slides",
                "has_pathologist_regions": bool((g.histopathology != "Unlabelled").any())}
        json.dump(meta, open(d / "meta.json", "w"), indent=1)
        print("prepared", sid, len(g), flush=True)
    model.build_reference()


if __name__ == "__main__":
    main()
