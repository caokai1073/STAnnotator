"""An existing annotation loaded from a CSV (one row per cell): the starting point of an annotation session.

The CSV needs a cell id column (the sample's cell ids, e.g. the h5ad obs names or Xenium cell_id) and a cell type
column; an optional lineage column groups the types (lineages for colors and priors). Cells of the sample that are
not in the file are Unassigned. The loaded labels stay the labels: the expression model trained on them only adds
a confidence and lets priors and corrections act (see Session._train_base). The app's own Export CSV (cell_id,
cell_type, lineage, ...) can be loaded back.
"""
import numpy as np, pandas as pd
from . import markers

ID_NAMES = ("cell_id", "cell", "barcode", "cellid", "cell_barcode", "obs_names", "index", "")
LABEL_NAMES = ("cell_type", "celltype", "cell_types", "annotation", "label", "labels", "cluster_name", "type")
LINEAGE_NAMES = ("lineage", "compartment", "group", "major_type", "level_1", "class")
MAX_TYPES = 500


def read(path):
    """The CSV (or TSV) as strings."""
    sep = "\t" if str(path).lower().endswith((".tsv", ".txt")) else ","
    return pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)


def _guess(cols, names):
    low = {c.lower(): c for c in cols}
    return next((low[n] for n in names if n in low), None)


def inspect(path, cell_ids):
    """For the dialog: columns, the guessed id / label / lineage columns and how many rows match the sample."""
    df = read(path)
    cols = list(df.columns)
    id_col = _guess(cols, ID_NAMES) or cols[0]
    label_col = _guess(cols, LABEL_NAMES) or next((c for c in cols if c != id_col), None)
    lineage_col = _guess(cols, LINEAGE_NAMES)
    ids = set(map(str, cell_ids))
    n_match = int(df[id_col].astype(str).isin(ids).sum())
    return {"path": str(path), "n_rows": int(len(df)), "columns": cols, "id_column": id_col, "label_column": label_col,
            "lineage_column": lineage_col, "n_matched": n_match, "n_cells": len(ids),
            "n_types": int(df[label_col].nunique()) if label_col else 0,
            "examples": {c: df[c].astype(str).head(3).tolist() for c in cols}}


def load(path, cell_ids, id_column, label_column, lineage_column=None):
    """Labels per cell of the sample (type indices; Unassigned where the file has no label), the types and their
    lineages."""
    df = read(path)
    for c in (id_column, label_column, lineage_column):
        if c and c not in df.columns:
            raise ValueError(f"the file has no column {c!r}")
    df = df[df[label_column].astype(str).str.strip() != ""].drop_duplicates(id_column, keep="last")
    lab = pd.Series(df[label_column].astype(str).str.strip().to_numpy(), index=df[id_column].astype(str).to_numpy())
    per_cell = lab.reindex(pd.Index(np.asarray(cell_ids).astype(str)))
    n_matched = int(per_cell.notna().sum())
    if n_matched == 0:
        raise ValueError(f"no cell id in column {id_column!r} matches this sample's cells")
    types = sorted(per_cell.dropna().unique().tolist())
    if len(types) > MAX_TYPES:
        raise ValueError(f"{len(types)} different labels in {label_column!r}: is it the cell type column?")
    if per_cell.isna().any() and markers.UNASSIGNED not in types:
        types = types + [markers.UNASSIGNED]
    code = {t: i for i, t in enumerate(types)}
    labels = np.array([code[t] for t in per_cell.fillna(markers.UNASSIGNED)], dtype=np.int64)
    lineage = {t: t for t in types}
    if lineage_column:
        lin = pd.Series(df[lineage_column].astype(str).to_numpy(), index=df[id_column].astype(str).to_numpy())
        tab = pd.crosstab(lab.to_numpy(), lin.reindex(lab.index).to_numpy())
        lineage.update({t: str(tab.loc[t].idxmax()) for t in tab.index if t in code})
    return {"labels": labels, "types": types, "lineage_of_type": lineage, "n_matched": n_matched,
            "n_rows": int(len(df))}


def mix(logp, labels, weight=0.9):
    """The loaded label keeps `weight` of the probability; the expression model shares the rest. So the loaded
    labels stay the labels, the model adds a confidence, and strong priors or corrections can still move a cell."""
    P = (1 - weight) * np.exp(logp.astype(np.float32))
    P[np.arange(len(labels)), labels] += weight
    return np.log(np.maximum(P, 1e-12)).astype(np.float16)
