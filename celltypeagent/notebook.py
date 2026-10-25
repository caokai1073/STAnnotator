"""Export an annotation as a self-contained Jupyter notebook: it needs only the data files and common Python packages
(numpy, pandas, scipy, scikit-learn, igraph, pynndescent, umap-learn, h5py, matplotlib), not this software.

The standard steps are written out as plain code: reading the cells, placing them on the image (the registration
affine), PCA / neighbors / Leiden / UMAP, NMF niches, figures and the cell table. The annotation algorithm of the app
(hierarchical marker clustering, the expression model, priors) is carried verbatim as source code in one folded cell
and run by a readable driver, so the replay uses the same code and gives the same labels. Every step's results at
export time are filled in (tables, figures), and running the notebook checks its labels, Leiden clusters and niches
against them. UMAP runs on all cores, so its coordinates differ slightly on a rerun.
"""
import base64, inspect, io, json, time, zlib
import numpy as np, pandas as pd
from . import config, embedding, niches, render

# app modules carried into the notebook, in import order (relative imports resolve inside the notebook)
CARRIED = ["vocab", "markers", "model", "constraints", "regions", "existing", "importer", "scref"]
CONFIG_STUB = '''"""Stand-in for the app's settings: only what the carried modules use."""
from pathlib import Path
HOME = WORK = SAMPLES = SESSIONS = Path(".")
def device():
    """As in the app: exact neighbors on a GPU when torch is installed, else pynndescent on the CPU. Run on the same
    kind of machine as the app (GPU or not) to get the same neighbors, hence the same clusters and labels."""
    try:
        import torch
    except ImportError:
        return "cpu"
    if torch.cuda.is_available():
        return "cuda"
    if getattr(torch.backends, "mps", None) is not None and torch.backends.mps.is_available():
        return "mps"
    return "cpu"
'''


def _png(fig):
    import matplotlib.pyplot as plt
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")}


def _code(src, outputs=(), hidden=False):
    md = {"jupyter": {"source_hidden": True}, "collapsed": True} if hidden else {}
    return {"cell_type": "code", "metadata": md, "execution_count": None, "outputs": list(outputs),
            "source": src.strip("\n")}


def _text(t):
    return {"output_type": "stream", "name": "stdout", "text": t if t.endswith("\n") else t + "\n"}


def _table(df):
    return {"output_type": "execute_result", "execution_count": None, "metadata": {},
            "data": {"text/plain": df.to_string(), "text/html": df.to_html()}}


def _image(b64):
    return {"output_type": "display_data", "metadata": {}, "data": {"image/png": b64, "text/plain": "<Figure>"}}


def _pack(a, dtype):
    """An integer array as a short string (zlib + base64)."""
    return base64.b64encode(zlib.compress(np.asarray(a, dtype).tobytes(), 9)).decode()


def _scatter(x, y, keys, title, invert=True, max_n=200000):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    idx = rng.choice(len(x), min(len(x), max_n), replace=False)
    names = pd.Series(keys).value_counts().index.tolist()
    pal = render.DISTINCT * 4
    fig, ax = plt.subplots(figsize=(8, 7))
    for i, n in enumerate(names[:40]):
        m = keys[idx] == n
        ax.scatter(x[idx][m], y[idx][m], s=0.3, c=pal[i % len(pal)], label=n, linewidths=0, rasterized=True)
    if invert:
        ax.invert_yaxis()
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10)
    ax.legend(markerscale=12, fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    return _png(fig)


def _sources():
    import importlib
    out = {"config": CONFIG_STUB}
    for name in CARRIED:
        out[name] = inspect.getsource(importlib.import_module(f"celltypeagent.{name}"))
    return out


PLOT = '''
def plot(x, y, keys, title, invert=True):
    fig, ax = plt.subplots(figsize=(8, 7))
    for n in pd.Series(keys).value_counts().index[:40]:
        m = keys == n
        ax.scatter(x[m], y[m], s=0.3, label=n, linewidths=0, rasterized=True)
    if invert:
        ax.invert_yaxis()
    ax.set_aspect("equal"); ax.set_title(title)
    ax.legend(markerscale=12, fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    plt.show()
'''


def build(s, version=None):
    """Write the notebook of a version into the session's exports and return its path."""
    from . import importer
    st = s.state(version)
    v = st["version"]
    meta, sch = s.meta, s.scheme
    labels = s.labels(v)
    n = len(s.cells)
    nb = []
    count = iter(range(1, 20))
    sec = lambda: next(count)  # section numbers in order

    dev = config.device()
    dev_text = {"cuda": "an NVIDIA GPU (exact neighbors)", "mps": "an Apple GPU (exact neighbors)"}.get(
        dev, "the CPU (approximate neighbors for large samples, pynndescent)")
    nb.append(_md(f"""
# Spatial Annotator analysis: sample `{s.sample_id}`, annotation `{s.id}`, version {v}

Exported {time.strftime('%Y-%m-%d %H:%M')} (Spatial Annotator build `{config.build()}`). This notebook needs only the
data files and common Python packages: `numpy pandas scipy scikit-learn python-igraph pynndescent umap-learn h5py
matplotlib`. The outputs below are the results at export time, so it can be read without running it; running it
from the top repeats the analysis and checks that the labels, Leiden clusters and niches come out the same.
The annotation algorithm of the app is carried as source code in the folded cell of section 1.

The neighbor graphs were computed with **{dev_text}**. Run the notebook on the same kind of machine to get
identical results: with a GPU, install PyTorch (`torch`) so that neighbors are exact as in the app; without one,
large samples use approximate neighbors (pynndescent), which change a share of the clusters and labels.
"""))
    # ---- 1. code carried from the app
    srcs = _sources()
    loader = (f"SOURCES = json.loads({json.dumps(json.dumps(srcs))})\n" + '''
import sys, types
sa = types.ModuleType("sa"); sa.__path__ = []; sys.modules["sa"] = sa   # a package made from the sources above
for name in ["config"] + ''' + json.dumps(CARRIED) + ''':
    mod = types.ModuleType(f"sa.{name}"); mod.__package__ = "sa"; sys.modules[f"sa.{name}"] = mod
    exec(compile(SOURCES[name], f"sa/{name}.py", "exec"), mod.__dict__)
    setattr(sa, name, mod)
print("algorithm modules:", ", ".join(SOURCES))''')
    nb.append(_md(f"## {sec()}. Packages and the app's annotation algorithm\nThe folded cell holds the source of the app's "
                  "modules used here (marker clustering, expression model, priors, regions, data readers), run as a "
                  "small package `sa`; click it to read the code."))
    nb.append(_code('''
import json, base64, zlib
import numpy as np, pandas as pd, scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
unpack = lambda s, dtype: np.frombuffer(zlib.decompress(base64.b64decode(s)), dtype)
''' + PLOT))
    nb.append(_code(loader, [_text("algorithm modules: " + ", ".join(srcs))], hidden=True))

    # ---- 2. data
    src = meta.get("source") or {}
    nb.append(_md(f"## {sec()}. The cells\nRead from the original files (edit the paths if they moved). The distance of "
                  "every cell to the epithelium (`depth_um`, used by depth priors) is computed as at import."))
    if src.get("kind") == "xenium":
        read = f'''
XENIUM = r"{src.get('folder')}"
from pathlib import Path
c = sa.importer._xenium_cells(Path(XENIUM))
bc, genes, X = sa.importer._xenium_matrix(Path(XENIUM))
pos = pd.Series(np.arange(len(bc)), index=bc)
c = c[c.cell_id.isin(pos.index)]
counts = X[pos[c.cell_id].to_numpy()].tocsr()
xy = c[["x_centroid", "y_centroid"]].to_numpy(float)
cell_ids = c.cell_id.to_numpy()'''
    elif src.get("kind") == "h5ad":
        read = f'''
H5AD = {json.dumps(src.get('h5ad'))}
SECTION_COLUMN, SECTION = {json.dumps(src.get('section_column'))}, {json.dumps(src.get('section'))}
parts = [sa.importer.read_h5ad(p, SECTION_COLUMN, SECTION) for p in H5AD]
parts = [q for q in parts if len(q["cell_id"])]
genes = parts[0]["genes"]
for q in parts[1:]:   # the same gene order in every file
    order = {{g: i for i, g in enumerate(q["genes"])}}
    q["X"] = q["X"][:, [order[g] for g in genes]]
counts = sp.vstack([q["X"] for q in parts]).tocsr()
xy = np.vstack([q["xy"] for q in parts])
cell_ids = np.concatenate([q["cell_id"] for q in parts])'''
    else:
        read = "raise RuntimeError('this sample was imported before its source files were recorded')"
    nb.append(_code(read + '''
depth_um = sa.importer.epithelial_depth(counts, genes, xy)
cells = pd.DataFrame({"cell": np.arange(len(cell_ids)), "cell_id": cell_ids, "x": xy[:, 0], "y": xy[:, 1], "depth_um": depth_um})
print(f"{len(cells):,} cells, {len(genes)} genes")''', [_text(f"{n:,} cells, {len(s.panel)} genes")]))

    # ---- 3. registration
    A = importer.placement(s.sample_id)
    by = meta.get("registration_by") or ("file" if src.get("registration") else "none")
    nb.append(_md(f"## {sec()}. Cells on the image\nThe registration: an affine from cell coordinates (um) to the image's "
                  f"full-resolution pixels, as used by the app (made by: {by}). The same numbers are in "
                  "**Export > Registration (CSV)**."))
    nb.append(_code(f'''
IMAGE = {json.dumps(meta.get("he_path"))}
A = np.array({json.dumps(np.round(A, 10).tolist())})   # px = A[0,0]*x + A[0,1]*y + A[0,2];  py = A[1,0]*x + A[1,1]*y + A[1,2]
pxy = xy @ A[:, :2].T + A[:, 2]
cells["px"], cells["py"] = pxy[:, 0], pxy[:, 1]
print(f"{{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f}} um per image pixel")''',
                    [_text(f"{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f} um per image pixel")]))

    # ---- 4. annotation
    rows = {"regions": dict(st["regions"], ops=[dict({k: x for k, x in o.items() if k != "cells"},
                                                     cells_packed=_pack(o.get("cells", []), np.int32))
                                                for o in st["regions"].get("ops", [])]),
            "constraints": st.get("constraints", []),
            "corrections": [dict({k: x for k, x in c.items() if k != "cells"}, cells_packed=_pack(c["cells"], np.int32))
                            for c in st.get("corrections", [])],
            "marker_ops": st.get("marker_step", {}).get("ops", [])}
    scheme = {k: sch.get(k) for k in ("kind", "text", "tree", "types", "resolution", "top", "lineage_of_type",
                                      "reference", "reference_counts", "existing", "leaf_types")}
    types = list(s.types)
    lineage = {t: s.lineage_of_type.get(t, t) for t in types}
    how = {"markers": "marker genes (a hierarchical marker list: Leiden clusters per level scored by their markers, "
                      "then an expression model trained on those labels)",
           "scref": "a single-cell reference (a CellTypist-like classifier on the genes shared with the panel)",
           "existing": "an existing annotation loaded from a CSV (its labels stay; the expression model adds "
                       "confidences)"}.get(sch["kind"], sch["kind"])
    priors = [c.get("text") or c["id"] for c in st["constraints"]]
    nb.append(_md(f"""## {sec()}. Annotation
Method: {how}. On top: {len(rows['regions']['ops'])} drawn region polygons, {len(priors)} priors
({'; '.join(priors) or 'none'}), {len(rows['corrections'])} hand corrections and {len(rows['marker_ops'])}
re-clusterings, as in version {v}. The regions and corrected cells are stored as cell lists."""))
    nb.append(_code(f'''
SCHEME = json.loads(r\'\'\'{json.dumps(scheme)}\'\'\')
STATE = json.loads(r\'\'\'{json.dumps(rows)}\'\'\')
TYPES = SCHEME["types"]; T = len(TYPES)
LINEAGE = json.loads(r\'\'\'{json.dumps(lineage)}\'\'\')
for o in STATE["regions"]["ops"]:
    o["cells"] = unpack(o.pop("cells_packed"), np.int32).tolist()
corr = {{}}
for c in STATE["corrections"]:          # hand-labeled cells (a later correction of the same cell wins)
    for cell in unpack(c["cells_packed"], np.int32):
        corr[int(cell)] = TYPES.index(c["label"])
lock_rows = np.array(list(corr), int); lock_labels = np.array(list(corr.values()), int)'''))
    if sch["kind"] == "markers":
        base = '''
import copy
def marker_tree(ops):        # the marker list with the edits of re-clusterings applied
    tree = copy.deepcopy(SCHEME["tree"])
    for op in ops:
        for t, mk in (op.get("markers") or {}).items():
            sa.markers._find(tree, t)["markers"] = sa.markers.parse(f"{t}: {mk}", genes)["tree"][0]["markers"]
    return tree
res = sa.markers.annotate(counts, genes, SCHEME["tree"], resolution=SCHEME["resolution"])   # marker step
for j, op in enumerate(STATE["marker_ops"]):
    res = sa.markers.rerun(res, counts, genes, marker_tree(STATE["marker_ops"][:j + 1]), op["target"], op["resolution"], prefix=f"{op['id']}:")
index = {t: i for i, t in enumerate(TYPES)}
marker_labels = np.array([index.get(x, -1) for x in res["labels"]])   # -1: left at a lineage / parent level
model = sa.model.train_pseudo(counts, marker_labels, T)                   # expression model
if len(lock_rows):
    rX, ry, _ = sa.model.pseudo_pool(counts, marker_labels, T)
    model = sa.model.finetune(model, counts[lock_rows], lock_labels, rX, ry)
logp = sa.model.predict(model, counts).astype(np.float32)
groups = sa.markers.descendants(SCHEME["tree"])  # cells left at a group get a type of that group only
for node in set(res["labels"].tolist()) - set(TYPES):
    rows_ = np.where(res["labels"] == node)[0]
    logp[np.ix_(rows_, np.where(~np.isin(TYPES, groups.get(node, [])))[0])] += np.float32(np.log(1e-4))
logp -= logp.max(1, keepdims=True)
logp = (logp - np.log(np.exp(logp).sum(1, keepdims=True))).astype(np.float16)'''
    elif sch["kind"] == "existing":
        e = sch["existing"]
        base = f'''
EXISTING = {json.dumps({k: e[k] for k in ("path", "id_column", "label_column", "lineage_column")})}
loaded = sa.existing.load(EXISTING["path"], cell_ids, EXISTING["id_column"], EXISTING["label_column"], EXISTING["lineage_column"])
train = loaded["labels"].copy()
if "Unassigned" in TYPES:
    train[train == TYPES.index("Unassigned")] = -1        # cells without a label do not train the model
model = sa.model.train_pseudo(counts, train, T)
if len(lock_rows):
    rX, ry, _ = sa.model.pseudo_pool(counts, train, T)
    model = sa.model.finetune(model, counts[lock_rows], lock_labels, rX, ry)
logp = sa.existing.mix(sa.model.predict(model, counts), loaded["labels"])   # the loaded label keeps 0.9'''
    elif sch["kind"] == "scref":
        r = sch["reference"]
        base = f'''
REFERENCE = {json.dumps(r)}
tr = sa.scref.load_training(REFERENCE["path"], REFERENCE["label_column"], genes)
model = sa.model._fit(sa.model.lognorm(tr["X"]), tr["y"], T, [SCHEME["reference_counts"][t] for t in TYPES])
mcounts = counts[:, tr["panel_idx"]]
if len(lock_rows):
    model = sa.model.finetune(model, mcounts[lock_rows], lock_labels, sa.model.lognorm(tr["X"]), tr["y"])
logp = sa.model.predict(model, mcounts)
if REFERENCE.get("majority_voting"):
    g, nbr = sa.markers._graph(sa.model.lognorm(counts))
    logp = sa.scref.majority_vote(logp, sa.markers._leiden(g, nbr, sa.scref.over_resolution(counts.shape[0]), 0))'''
    else:
        base = "raise RuntimeError('this annotation method cannot be replayed outside the lab')"
    tnames = sorted(set(labels))
    tcodes = np.searchsorted(tnames, labels)
    comp = pd.Series(labels).value_counts().rename_axis("cell type").to_frame("cells").head(30)
    nb.append(_code(base + f'''
base_pred, _ = sa.constraints.posterior(logp, [], None, TYPES)          # expression only
region = sa.regions.assign(cells, STATE["regions"])                    # drawn regions
nn = cKDTree(xy).query(xy, k=sa.constraints.NEIGHBORS_K + 1)[1][:, 1:]
neighbor_frac = lambda names: np.isin(base_pred, [TYPES.index(t) for t in names])[nn].mean(1)
ctx = {{"region": region, "depth_um": depth_um, "cell": cells.cell.to_numpy(), "neighbor_frac": neighbor_frac}}
pred, conf, _, _ = sa.constraints.posterior(logp, STATE["constraints"], ctx, TYPES, top=3)   # priors
pred[lock_rows], conf[lock_rows] = lock_labels, 1.0                     # hand corrections
labels = np.array(TYPES, dtype=object)[pred]
ORIGINAL = np.array({json.dumps(tnames)}, dtype=object)[unpack("{_pack(tcodes, np.uint16)}", np.uint16)]
same = labels == ORIGINAL
print("labels identical to the original" if same.all() else f"{{same.mean():.2%}} of the labels as in the original "
      f"({{(~same).sum():,}} cells differ; an annotation made with an earlier app version can differ)")
pd.Series(labels).value_counts().rename_axis("cell type").to_frame("cells").head(30)''',
                    [_text("labels identical to the original"), _table(comp)]))

    # ---- 5. spatial map
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    lin = np.array([lineage[t] for t in labels], dtype=object)
    nb.append(_md(f"## {sec()}. Spatial map\nCells at their place on the image, colored by lineage."))
    nb.append(_code('''
lineages = np.array([LINEAGE[t] for t in labels], dtype=object)
plot(cells.px.to_numpy(), cells.py.to_numpy(), lineages, "lineages")''',
                    [_image(_scatter(px, py, lin, f"sample {s.sample_id}, lineages"))]))

    # ---- 6. embedding
    em = embedding.current(s.sample_id)
    if em:
        Y = embedding.coords(s.sample_id, "umap")
        lab = embedding.leiden(s.sample_id)
        nb.append(_md(f"## {sec()}. PCA, neighbors, Leiden and UMAP\nFrom expression only: log-normalized counts, PCA of the "
                      f"standardized genes ({em['n_pcs']} components), {em['n_neighbors']} nearest neighbors, Leiden "
                      f"(resolution {em['resolution']}: {em['n_clusters']} clusters) and UMAP (min_dist {em['min_dist']}) "
                      "on the same neighbors, started from PC1/PC2."))
        nb.append(_code(f'''
import umap
N_PCS, K, RESOLUTION, MIN_DIST = {em['n_pcs']}, {em['n_neighbors']}, {em['resolution']}, {em['min_dist']}
Z = sa.markers.pca(sa.model.lognorm(counts), N_PCS)[:, ::-1].copy()     # PC1 first
dist, nbr = sa.markers.knn(Z, K)
leiden = sa.markers._leiden(sa.markers.knn_graph(dist, nbr), nbr, RESOLUTION, 0)
print(leiden.max() + 1, "Leiden clusters; identical to the original:",
      bool((leiden == unpack("{_pack(lab, np.int32)}", np.int32)).all()))
n = len(Z)
init = Z[:, :2].astype(np.float64); init = 10 * (init - init.mean(0)) / np.abs(init - init.mean(0)).max()
umap_xy = umap.UMAP(n_neighbors=K + 1, min_dist=MIN_DIST, init=init, n_jobs=-1, low_memory=True,
                    precomputed_knn=(np.c_[np.arange(n), nbr], np.c_[np.zeros(n, np.float32), dist].astype(np.float32), None)
                    ).fit_transform(Z)
plot(umap_xy[:, 0], umap_xy[:, 1], labels, "UMAP, cell types", invert=False)''',
                        [_text(f"{em['n_clusters']} Leiden clusters; identical to the original: True"),
                         _image(_scatter(Y[:, 0], Y[:, 1], np.asarray(labels, dtype=object), "UMAP, cell types", invert=False))]))

    # ---- 7. niches
    nmeta, nlab = niches.current(s)
    if nmeta:
        summ = niches.summary(s, v)
        tab = pd.DataFrame([{"niche": x["name"], "main cell types": ", ".join(
            f"{k} {round(100 * p)}%" for k, p in list(x["cell_types"].items())[:3])} for x in summ["niches"]]).set_index("niche")
        names = niches.niche_names(s)
        if nmeta["identity"] == "nmf":
            ident = f'''
from sklearn.decomposition import NMF
Xn = sa.model.lognorm(counts).astype(np.float32).tocsr()
fit = np.sort(np.random.default_rng(0).choice(Xn.shape[0], min(60000, Xn.shape[0]), replace=False))
nmf = NMF(n_components={nmeta['n_programs']}, init="nndsvda", max_iter=300, random_state=0, tol=1e-4).fit(Xn[fit])
W = np.vstack([nmf.transform(Xn[i:i + 50000]) for i in range(0, Xn.shape[0], 50000)]).astype(np.float32)
identity = W / np.maximum(W.sum(1, keepdims=True), 1e-9)    # each cell's composition of gene programs'''
        else:
            ident = '''
identity = np.zeros((len(pred), T), np.float32); identity[np.arange(len(pred)), pred] = 1   # cell types'''
        if nmeta.get("radius_um"):
            kk = f'''
kk = np.clip(cKDTree(xy).query_ball_point(xy, r={nmeta['radius_um']}, return_length=True) - 1, 5, 80)   # by density'''
        else:
            kk = f'''
kk = np.full(len(xy), {nmeta['k']})
for name, k in {json.dumps(nmeta.get('k_by_region') or {})}.items():   # k per drawn region
    kk[region == name] = k'''
        nb.append(_md(f"## {sec()}. Niches\nEach cell's identity (" + ("its composition of NMF gene programs" if nmeta["identity"] == "nmf"
                      else "its cell type") + ") is averaged over the cell and its nearest cells, and the averages are "
                      f"clustered (k-means) into {nmeta['n_niches']} niches, numbered by size."))
        nb.append(_code(ident + kk + f'''
from sklearn.cluster import KMeans
nn_ = cKDTree(xy).query(xy, k=int(kk.max()) + 1)[1]
avg = np.zeros_like(identity)
for k in np.unique(kk):
    r = np.where(kk == k)[0]
    for i in range(0, len(r), 20000):
        avg[r[i:i + 20000]] = identity[nn_[r[i:i + 20000], :k + 1]].mean(1)
km = KMeans(n_clusters={nmeta['n_niches']}, n_init=4, random_state=0).fit(avg)
order = np.argsort(-np.bincount(km.labels_, minlength={nmeta['n_niches']}))
niche = np.argsort(order)[km.labels_]
NAMES = {json.dumps(names)}
print("niches identical to the original:", bool((niche == unpack("{_pack(nlab, np.int16)}", np.int16)).all()))
niche_name = np.array(NAMES, dtype=object)[niche]
pd.DataFrame({{"niche": niche_name, "cell type": labels}}).groupby("niche")["cell type"].agg(
    lambda x: ", ".join(f"{{k}} {{round(100 * v)}}%" for k, v in x.value_counts(normalize=True).head(3).items())).to_frame("main cell types")''',
                        [_text("niches identical to the original: True"), _table(tab),
                         _image(_scatter(px, py, np.array(names, dtype=object)[nlab], "niches"))]))

    # ---- 8. table
    table = s.export_table(v)
    cols = ["cell_id", "region", "cell_type", "lineage", "confidence"] + [c for c in ("niche", "leiden", "umap_1", "umap_2")
                                                                           if c in table] + ["x_um", "y_um"]
    nb.append(_md(f"## {sec()}. Cell table\nOne row per cell, written to `cells.csv`."))
    extra = []
    if nmeta:
        extra.append('table["niche"] = niche_name')
    if em:
        extra += ['table["leiden"] = ["L" + str(i + 1) for i in leiden]',
                  'table["umap_1"], table["umap_2"] = umap_xy[:, 0].round(4), umap_xy[:, 1].round(4)']
    nb.append(_code('''
table = pd.DataFrame({"cell_id": cell_ids, "region": np.where(region == "unknown", "", region), "cell_type": labels,
                      "lineage": [LINEAGE[t] for t in labels], "confidence": conf.round(4)})
''' + "\n".join(extra) + '''
table["x_um"], table["y_um"] = xy[:, 0].round(2), xy[:, 1].round(2)
table.to_csv("cells.csv", index=False)
table.head()''', [_table(table[cols].head())]))

    doc = {"nbformat": 4, "nbformat_minor": 5, "cells": nb,
           "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                        "language_info": {"name": "python"}}}
    for i, c in enumerate(doc["cells"]):
        c["id"] = f"c{i:02d}"
    d = s.dir / "exports" / f"v{v:04d}"
    d.mkdir(parents=True, exist_ok=True)
    out = d / f"{s.sample_id}_v{v}_analysis.ipynb"
    json.dump(doc, open(out, "w"), indent=1)
    return str(out)


def registration_csv(sample_id):
    """The affine that places the sample's cells on its image, as an a11..ty CSV (importable as a registration)."""
    from . import importer
    m = importer._meta(sample_id)
    A = importer.placement(sample_id)
    return pd.DataFrame([{"sample": sample_id, "a11": A[0, 0], "a12": A[0, 1], "a21": A[1, 0], "a22": A[1, 1],
                          "tx": A[0, 2], "ty": A[1, 2], "um_per_px": float(1 / np.sqrt(abs(np.linalg.det(A[:, :2])))),
                          "made_by": m.get("registration_by") or "file", "image": m.get("he_path")}]).to_csv(index=False)
