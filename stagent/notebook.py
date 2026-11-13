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
# visium and stereoseq before importer, which imports them; rds (R files, read without R) for samples and references
# imported from an .rds / .rda: importer converts it to an h5ad again, in converted/ next to the notebook
CARRIED = ["osutil", "csrutil", "harmony", "timing", "vocab", "markers", "model", "constraints", "regions", "existing", "visium", "stereoseq", "rds", "importer", "scref",
           "trajectory"]
# printed after a comparison with the original that a replay on one slice cannot match (computed on all slices)
_JOINT_NOTE = (', "(computed on all slices of the project together: this slice alone gives other ones; '
               'Export > Save project replays them exactly)"')
CONFIG_STUB = '''"""Stand-in for the app's settings: only what the carried modules use."""
from pathlib import Path
HOME = WORK = SAMPLES = SESSIONS = Path(".")
RECORD_TIMINGS = False  # the app keeps how long its steps take; the notebook does not

def device():
    """As in the app: the GPU that exact neighbor searches run on when torch is installed (it gives the same
    neighbors as the CPU, only faster)."""
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


def _py(x):
    """A value as Python source (plain types: no numpy scalars; None, True and False as Python writes them)."""
    return repr(json.loads(json.dumps(x)))


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
    """Source of the carried modules: from the package, or (standalone app: only compiled code) from the copies the
    app bundle keeps in algorithm_src/."""
    import importlib
    out = {"config": CONFIG_STUB}
    for name in CARRIED:
        try:
            out[name] = inspect.getsource(importlib.import_module(f"stagent.{name}"))
        except (OSError, TypeError):
            f = config.PACKAGE / "algorithm_src" / f"{name}.py"
            if not f.exists():
                raise RuntimeError(f"the source of {name}.py is missing from this app build")
            out[name] = f.read_text(encoding="utf-8")  # not the Windows code page (GBK on Chinese Windows)
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

    from . import markers
    dev = config.device()
    knn = s.knn_mode()
    if knn == "auto":  # annotations made before the choice of neighbor search
        dev_text = {"cuda": "an NVIDIA GPU (exact neighbors)", "mps": "an Apple GPU (exact neighbors)"}.get(
            dev, "the CPU (approximate neighbors for large samples, pynndescent)")
        knn_text = f"""The neighbor graphs were computed with **{dev_text}**. Run the notebook on the same kind of machine to get
identical results: with a GPU, install PyTorch (`torch`) so that neighbors are exact as in the app; without one,
large samples use approximate neighbors (pynndescent), which change a share of the clusters and labels."""
    elif knn == "exact":
        knn_text = """The marker step used **exact neighbors**, which are the same on every computer, with or without a GPU: the
notebook gives the same clusters and labels anywhere. Installing PyTorch (`torch`) only makes the search faster."""
    else:
        used = s.marker_result(st).get("knn") or markers.knn_backend("approx")[1]
        knn_text = f"""The marker step used **approximate neighbors** ({used.split("(", 1)[-1].rstrip(")")}). With those
versions the notebook gives the same clusters and labels; other versions can change a share of them."""
    nb.append(_md(f"""
# STAnnotator analysis: sample `{s.sample_id}`, annotation `{s.id}`, version {v}

Exported {time.strftime('%Y-%m-%d %H:%M')} (STAnnotator build `{config.build()}`). This notebook needs only the
data files and common Python packages: `numpy pandas scipy scikit-learn python-igraph pynndescent umap-learn h5py
matplotlib{" pillow" if (meta.get("source") or {}).get("segmentation") else ""}`, installed by the first code cell. The outputs below are the results at export time, so it can be read without running it; running it
from the top repeats the analysis and checks that the labels, Leiden clusters and niches come out the same.
The annotation algorithm of the app is carried as source code in the folded cell of section 1.

{knn_text}
"""))
    nb += _package_cells(sec, ["pillow"] if (meta.get("source") or {}).get("segmentation") else [])
    nb += _data_cells(s, sec)

    # ---- 4. annotation
    rows = {"regions": dict(st["regions"], ops=[dict({k: x for k, x in o.items() if k != "cells"},
                                                     cells_packed=_pack(o.get("cells", []), np.int32))
                                                for o in st["regions"].get("ops", [])]),
            "constraints": st.get("constraints", []),
            "corrections": [dict({k: x for k, x in c.items() if k != "cells"}, cells_packed=_pack(c["cells"], np.int32))
                            for c in st.get("corrections", [])],
            "marker_ops": st.get("marker_step", {}).get("ops", [])}
    scheme = {k: sch.get(k) for k in ("kind", "text", "tree", "types", "resolution", "top", "lineage_of_type",
                                      "reference", "reference_counts", "existing", "leaf_types", "neighbors",
                                      "spillover")}
    scheme["knn"] = s.knn_mode()
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
SCHEME = {_py(scheme)}
STATE = {_py(rows)}
TYPES = SCHEME["types"]; T = len(TYPES)
LINEAGE = {_py(lineage)}
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
SPILL = {"xy": xy} if SCHEME.get("spillover") else None   # transcripts spilled over from neighbors left out (option)
res = sa.markers.annotate(counts, genes, SCHEME["tree"], resolution=SCHEME["resolution"], knn=SCHEME["knn"], spill=SPILL)   # marker step
for j, op in enumerate(STATE["marker_ops"]):
    res = sa.markers.rerun(res, counts, genes, marker_tree(STATE["marker_ops"][:j + 1]), op["target"], op["resolution"], prefix=f"{op['id']}:", knn=SCHEME["knn"], spill=SPILL)
marker_labels = sa.markers.training_labels(res, TYPES)   # -1: left at a lineage / parent level
model = sa.model.train_pseudo(counts, marker_labels, T, keep=sa.markers.gene_columns(marker_tree(STATE["marker_ops"]), genes))   # expression model
if len(lock_rows):
    rX, ry, _ = sa.model.pseudo_pool(counts, marker_labels, T)
    model = sa.model.finetune(model, counts[lock_rows], lock_labels, rX, ry)
logp = sa.model.predict(model, counts)
# cells left at a group get a type of that group; types the marker step found no cell of are not given; the cells it
# could not type get one from their marker genes and their neighbors' types
logp = sa.markers.restrict(logp, res, marker_tree(STATE["marker_ops"]), TYPES, lock_labels, counts=counts, panel=genes,
                           xy=xy, neighbors=bool(SCHEME.get("neighbors")), spill=bool(SCHEME.get("spillover")),
                           legacy=SCHEME.get("neighbors") is None)   # (an annotation made before 0.2.57: its rules then)'''
    elif sch["kind"] == "existing":
        e = sch["existing"]
        base = f'''
EXISTING = {_py({k: e[k] for k in ("path", "id_column", "label_column", "lineage_column")})}
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
REFERENCE = {_py(r)}
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
ORIGINAL = np.array({_py(tnames)}, dtype=object)[unpack("{_pack(tcodes, np.uint16)}", np.uint16)]
same = labels == ORIGINAL
print("labels identical to the original" if same.all() else f"{{same.mean():.2%}} of the labels as in the original "
      f"({{(~same).sum():,}} cells differ; {"it was annotated on all slices of the project together, and this replays it on this slice alone: Export > Save project replays it exactly" if sch.get("joint") else "an annotation made with an earlier app version can differ"})")
pd.Series(labels).value_counts().rename_axis("cell type").to_frame("cells").head(30)''',
                    [_text("labels identical to the original" if not sch.get("joint") else
                           "(the annotation was made on all slices of the project together; a replay on this slice "
                           "alone can label some cells differently)"), _table(comp)]))
    if sch.get("joint"):
        nb.insert(len(nb) - 1, _md(f"*This annotation was made on all {sch['joint'].get('slices', 'the')} slices of project "
                                   f"`{sch['joint'].get('project')}` together (one marker step and one model). The cell below "
                                   "replays it on this slice alone, so its clusters and some labels differ from the app's; "
                                   "**Export > Save project…** replays it exactly.*"))

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
FIT_MAX = {em.get('pca_fit_cells') or None}     # the cells the components were fitted on (None: every cell)
Z = sa.markers.pca(sa.model.lognorm(counts), N_PCS, fit_max=FIT_MAX)[:, ::-1].copy()     # PC1 first
dist, nbr = sa.markers.knn(Z, K)
leiden = sa.markers._leiden(sa.markers.knn_graph(dist, nbr), nbr, RESOLUTION, 0)
print(leiden.max() + 1, "Leiden clusters; identical to the original:",
      bool((leiden == unpack("{_pack(lab, np.int32)}", np.int32)).all()){_JOINT_NOTE if em.get("project") else ""})
n = len(Z)
init = Z[:, :2].astype(np.float64); init = 10 * (init - init.mean(0)) / np.abs(init - init.mean(0)).max()
umap_xy = umap.UMAP(n_neighbors=K + 1, min_dist=MIN_DIST, init=init, n_jobs=-1, low_memory=True,
                    **(dict(compatibility_layout=True) if "compatibility_layout" in umap.UMAP.__init__.__code__.co_varnames else dict()),  # umap 0.6: the 0.5 layout
                    precomputed_knn=(np.c_[np.arange(n), nbr], np.c_[np.zeros(n, np.float32), dist].astype(np.float32), None)
                    ).fit_transform(Z)
plot(umap_xy[:, 0], umap_xy[:, 1], labels, "UMAP, cell types", invert=False)''',
                        [_text(f"{em['n_clusters']} Leiden clusters; identical to the original: True" if not em.get("project")
                               else "(computed on all slices of the project together"
                               + (", corrected with Harmony" if em.get("harmony") else "")
                               + "; a replay on this slice alone gives other clusters)"),
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
for name, k in {_py(nmeta.get('k_by_region') or {})}.items():   # k per drawn region
    kk[region == name] = k'''
        other_v = nmeta.get("version") not in (None, v) and (nmeta["identity"] != "nmf" or nmeta.get("k_by_region"))
        caveat = (f" *They were computed on all slices of project `{nmeta['project']}` together, so a replay on this "
                  "slice alone differs.*" if nmeta.get("project") else
                  f" *They were computed on version {nmeta.get('version')}; this notebook replays version {v}'s "
                  "labels and regions, so they can differ.*" if other_v else "")
        nb.append(_md(f"## {sec()}. Niches\nEach cell's identity (" + ("its composition of NMF gene programs" if nmeta["identity"] == "nmf"
                      else "its cell type") + ") is averaged over the cell and its nearest cells, and the averages are "
                      f"clustered (k-means) into {nmeta['n_niches']} niches, numbered by size." + caveat))
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
NAMES = {_py(names)}
print("niches identical to the original:", bool((niche == unpack("{_pack(nlab, np.int16)}", np.int16)).all()){_JOINT_NOTE if nmeta.get("project") else ""})
niche_name = np.array(NAMES, dtype=object)[niche]
pd.DataFrame({{"niche": niche_name, "cell type": labels}}).groupby("niche")["cell type"].agg(
    lambda x: ", ".join(f"{{k}} {{round(100 * v)}}%" for k, v in x.value_counts(normalize=True).head(3).items())).to_frame("main cell types")''',
                        [_text("niches identical to the original: True" if not caveat else "(see the note above: a replay "
                               "can differ)"), _table(tab),
                         _image(_scatter(px, py, np.array(names, dtype=object)[nlab], "niches"))]))

    # ---- trajectory: the annotation's last run (Trajectory tab or the agent), replayed from its cells
    _trajectory_section(s, nb, sec)

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
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
    return str(out)


TRAJ_CODE = '''
TRAJ = {traj}
rows_of = lambda s: unpack(s, np.int32)
MASK = {{k: np.isin(np.arange(len(cells)), rows_of(v)) for k, v in {masks}.items()}}
one = np.zeros(len(cells), dtype=object)   # one slice
Xn = sa.model.lognorm(counts)
if TRAJ["method"] == "spatial":
    t, raw = sa.trajectory.spatial_axis(xy, one, MASK["start"], MASK.get("end"))
    t[~MASK["within"]] = np.nan
    r = np.where(MASK["within"], raw, np.nan)
    if MASK.get("end") is None and np.isfinite(r).any():   # a distance: scaled over the cells ordered, as in the app
        t = np.clip(r / max(float(np.nanpercentile(r, 99)), 1e-9), 0, 1)
else:
    rows = np.where(MASK["within"])[0]
    Z = sa.markers.pca(Xn[rows], sa.trajectory.N_PCS)[:, ::-1].copy()
    tt, info = sa.trajectory.pseudotime(Z, xy[rows], one[rows], MASK["root"][rows], alpha=TRAJ["alpha"])
    t = np.full(len(cells), np.nan); t[rows] = tt
T_APP = unpack({t_app}, np.float32)   # the order at export time (-1: not ordered)
both = (T_APP >= 0) & np.isfinite(t)
print(f"{{np.isfinite(t).sum():,}} cells ordered; correlation with the app's order: {{np.corrcoef(t[both], T_APP[both])[0, 1]:.4f}}")
trends, gene_table = sa.trajectory.gene_trends(Xn, t, list(genes))
for k in ("rising", "transient", "falling"):
    print(k + ":", ", ".join(x["gene"] for x in trends[k][:10]))
H = trends["heatmap"]
fig, ax = plt.subplots(figsize=(6, 0.18 * len(H["genes"]) + 1))
ax.imshow(np.clip(np.array(H["z"]), -2, 2), cmap="RdBu_r", aspect="auto", vmin=-2, vmax=2)
ax.set_yticks(range(len(H["genes"]))); ax.set_yticklabels(H["genes"], fontsize=7); ax.set_xlabel("bins along the trajectory")
plt.show()'''


def _trajectory_section(s, nb, sec):
    """The last trajectory of each method of the annotation, replayed: its cells (the masks of the run, packed rows),
    the same functions (trajectory.py, carried), the genes along it, checked against the order at export time."""
    from . import trajectory
    for method in trajectory.all_last(s):
        _trajectory_one(s, nb, sec, method)


def _trajectory_one(s, nb, sec, method):
    from . import trajectory
    meta, t, raw = trajectory.last(s, method)
    if meta is None:
        return
    z = np.load(trajectory._files(s, method)[1])
    masks = {k[5:]: np.where(z[k])[0] for k in z.files if k.startswith("mask_")}
    if "within" not in masks or meta.get("pooled"):
        nb.append(_md(f"## {sec()}. Trajectory ({trajectory.METHOD_NAMES[method]})\nThe last trajectory ({meta['label']}) "
                      + ("pooled the slices of the project: run it again on this slice to replay it here."
                         if meta.get("pooled") else "was made before the app kept its cells: run it again to replay it here.")))
        return
    if meta["method"] == "spatial":
        how = ("each cell's place between the start and end groups in the tissue: t = d_start / (d_start + d_end)"
               if "end" in masks else "the distance to the start group (um)")
    else:
        how = "diffusion pseudotime from the root group on the PCA neighbor graph" + (
            f", with tissue neighbors joined in (weight of space {meta.get('alpha')})" if meta["method"] == "mixed" else "")
    nb.append(_md(f"""## {sec()}. Trajectory ({trajectory.METHOD_NAMES[method]})
The annotation's last trajectory of this method: **{meta['label']}** ({how}). The cells it orders and its groups are given as their
rows (as at the time of the run); the functions are the app's (`trajectory.py`). The genes along it: the correlation of
each gene with the order (rising, falling) and its share of variance between {trajectory.NBINS} bins (transient).
The pseudotime uses approximate neighbors for large samples, so it can differ slightly from the app's."""))
    code = TRAJ_CODE.format(traj=_py({"method": meta["method"], "alpha": meta.get("alpha") or 0.0}),
                            masks=_py({k: _pack(v, np.int32) for k, v in masks.items()}),
                            t_app=_py(_pack(np.nan_to_num(t, nan=-1).astype(np.float32), np.float32)))
    out = (f"{int(np.isfinite(t).sum()):,} cells ordered; correlation with the app's order: 1.0000 (at export)\n"
           + "\n".join(f"{k}: " + ", ".join(meta.get("top", {}).get(k, [])) for k in ("rising", "transient", "falling")))
    nb.append(_code(code, [_text(out)]))


PINNED = ("pynndescent", "numba", "igraph")   # their versions change the clusters (pinned in the app too)


def _install_cell(extra=()):
    """A cell that installs the packages into the notebook's Python: the app's versions of those that change results."""
    from importlib import metadata
    def spec(name):
        try:
            return f"{name}=={metadata.version(name)}" if name in PINNED else name
        except metadata.PackageNotFoundError:
            return name
    pkgs = ["numpy", "pandas", "scipy", "scikit-learn", "igraph", "pynndescent", "numba", "umap-learn", "h5py",
            "matplotlib"] + list(extra)
    return _code("# installs the packages this notebook needs into its Python (once; nothing to do when they are there);\n"
                 "# pynndescent, numba and igraph at the app's versions, so the clusters and labels come out the same\n"
                 "%pip install -q " + " ".join(spec(x) for x in pkgs))


def _package_cells(sec, extra=()):
    """The cells that install and import the packages and run the app's carried modules as the small package `sa`."""
    out = []
    # ---- 1. code carried from the app
    srcs = _sources()
    loader = (f"SOURCES = json.loads({json.dumps(json.dumps(srcs))})\n" + '''
import sys, types
sa = types.ModuleType("sa"); sa.__path__ = []; sys.modules["sa"] = sa   # a package made from the sources above
for name in ["config"] + ''' + _py(CARRIED) + ''':
    mod = types.ModuleType(f"sa.{name}"); mod.__package__ = "sa"; sys.modules[f"sa.{name}"] = mod
    exec(compile(SOURCES[name], f"sa/{name}.py", "exec"), mod.__dict__)
    setattr(sa, name, mod)
print("algorithm modules:", ", ".join(SOURCES))''')
    out.append(_md(f"## {sec()}. Packages and the app's annotation algorithm\nThe first cell installs the Python packages "
                  "(with pip, into the notebook's Python). The folded cell holds the source of the app's "
                  "modules used here (marker clustering, expression model, priors, regions, data readers), run as a "
                  "small package `sa`; click it to read the code."))
    out.append(_install_cell(extra))
    out.append(_code('''
import os
# the approximate neighbor search runs on 8 threads, as in the app (its result depends on the thread count): numba must
# allow 8 threads even with fewer cores. Set before numba is imported.
os.environ["NUMBA_NUM_THREADS"] = str(max(int(os.environ.get("NUMBA_NUM_THREADS") or 0), os.cpu_count() or 1, 8))
# on a Mac with MLX the exact neighbor search runs on the GPU: in full float32, as in the app (not TF32)
os.environ["MLX_ENABLE_TF32"] = "0"
import json, base64, zlib
import numpy as np, pandas as pd, scipy.sparse as sp
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
unpack = lambda s, dtype: np.frombuffer(zlib.decompress(base64.b64decode(s)), dtype)
''' + PLOT))
    out.append(_code(loader, [_text("algorithm modules: " + ", ".join(srcs))], hidden=True))

    return out


def _data_cells(s, sec):
    """The notebook cells that read a slice's cells from its original files (with its filter) and place them on its
    image (registration and corrections), leaving cells, counts, genes, xy, cell_ids and depth_um defined."""
    from . import importer
    meta, n, out = s.meta, len(s.cells), []
    # ---- 2. data
    src = meta.get("source") or {}
    unit = {"visium_hd": "cells built from the 2 um bins by the cell outlines" if src.get("unit") == "cells" else
            f"bins of {src.get('bin_um')} um", "visium": "spots", "stereoseq": "cells (cell bins)" if src.get("unit") == "cells"
            else f"bins of {src.get('bin_um')} um"}.get(src.get("kind"))
    out.append(_md(f"## {sec()}. The cells\nRead from the original files (edit the paths if they moved)" +
                  (f", as at import: {unit}" if unit else "") + ". The distance of "
                  "every cell to the epithelium (`depth_um`, used by depth priors) is computed as at import."))
    if src.get("kind") == "xenium":
        read = f'''
XENIUM = {_py({"key": src.get("key") or "", "folder": src.get("folder"), "files": src["files"]}) if src.get("files") else _py(src.get("folder"))}
c = sa.importer._xenium_cells(XENIUM)
bc, genes, X = sa.importer._xenium_matrix(XENIUM)
pos = pd.Series(np.arange(len(bc)), index=bc)
c = c[c.cell_id.isin(pos.index)]
counts = X[pos[c.cell_id].to_numpy()].tocsr()
xy = c[["x_centroid", "y_centroid"]].to_numpy(float)
cell_ids = c.cell_id.to_numpy()'''
    elif src.get("kind") == "h5ad":
        read = f'''
H5AD = {_py(src.get('h5ad'))}
SECTION_COLUMN, SECTION = {_py(src.get('section_column'))}, {_py(src.get('section'))}
X_COLUMN, Y_COLUMN, COUNTS_LAYER = {_py(src.get('x_column'))}, {_py(src.get('y_column'))}, {_py(src.get('counts_layer'))}   # None: found as at import
parts = [sa.importer.read_h5ad(p, SECTION_COLUMN, SECTION, X_COLUMN, Y_COLUMN, COUNTS_LAYER) for p in H5AD]
parts = [q for q in parts if len(q["cell_id"])]
genes = parts[0]["genes"]
for q in parts[1:]:   # the same gene order in every file
    order = {{g: i for i, g in enumerate(q["genes"])}}
    q["X"] = q["X"][:, [order[g] for g in genes]]
counts = sp.vstack([q["X"] for q in parts]).tocsr()
xy = np.vstack([q["xy"] for q in parts])
cell_ids = np.concatenate([q["cell_id"] for q in parts])'''
    elif src.get("kind") in ("cosmx", "merfish"):
        read = f'''
KIND, FILES = {_py(src["kind"])}, {_py(src.get("files"))}
cell_ids, xy, genes, counts = sa.importer.flat_cells_counts(KIND, FILES)   # metadata order, cells with counts'''
    elif src.get("kind") in ("visium_hd", "visium"):
        read = f'''
FOLDER = {_py(src.get("folder"))}            # the square_XXXum folder (Visium HD) or outs/ (Visium)
BIN_UM = {_py(src.get("source_bin_um") if src.get("kind") == "visium_hd" else None)}   # the bins read (None: Visium spots)
SEGMENTATION = {_py(src.get("segmentation"))}   # cell outlines (Visium HD): cells from the bins, split by area
if SEGMENTATION:
    cell_ids, pxy, genes, counts, _ = sa.visium.read_cells(FOLDER, BIN_UM, SEGMENTATION)
else:
    cell_ids, pxy, genes, counts, _, _ = sa.visium.read(FOLDER, BIN_UM)
counts = counts.tocsr()
xy = pxy * sa.visium.scale(FOLDER, BIN_UM)[1]   # full-resolution pixels -> um'''
    elif src.get("kind") == "stereoseq":
        read = f'''
FILE, CELLS = {_py(src.get("file"))}, {_py(src.get("unit") == "cells")}   # a cell-bin GEF or a GEM; cells, or 8 um bins
if FILE.endswith(".gef"):
    cell_ids, dnb, genes, counts = sa.stereoseq.read_cellbin_gef(FILE)
else:
    cell_ids, dnb, genes, counts = sa.stereoseq.read_gem(FILE, CELLS)
ok = np.asarray(counts.sum(1)).ravel() > 0
cell_ids, dnb, counts = cell_ids[ok], dnb[ok], counts[ok].tocsr()
xy = dnb * sa.stereoseq.DNB_UM   # DNBs -> um'''
    else:
        read = "raise RuntimeError('this sample was imported before its source files were recorded')"
    read += '''
ok = np.isfinite(xy).all(1)   # cells without a position (NaN) are left out, as at import
counts, xy, cell_ids = sp.csr_matrix(counts)[np.where(ok)[0]], xy[ok], np.asarray(cell_ids)[ok]'''
    flt = meta.get("filter")
    if flt:  # the sample's QC filter (qcfilter.py), applied to the cells as imported
        out.append(_md(f"The sample was filtered in the app ({', '.join(f'{k} {x:g}' for k, x in flt['params'].items())}): "
                      f"{flt['n_cells_removed']:,} cells and {len(flt.get('genes_removed') or [])} genes were removed; "
                      "the same filter is applied below."))
        read += f'''
FILTER, GENES_REMOVED = {_py(flt["params"])}, {_py(flt.get("genes_removed") or [])}   # the app's QC filter
counts = counts.tocsr(); counts.eliminate_zeros()
total, n_genes = np.asarray(counts.sum(1)).ravel(), np.diff(counts.indptr)
keep = np.ones(len(total), bool)
if "min_counts" in FILTER: keep &= total >= FILTER["min_counts"]
if "max_counts" in FILTER: keep &= total <= FILTER["max_counts"]
if "min_genes" in FILTER: keep &= n_genes >= FILTER["min_genes"]
mito = np.array([str(g).lower().startswith("mt-") for g in genes])
if "max_mito_pct" in FILTER and mito.any():
    keep &= 100 * np.asarray(counts[:, mito].sum(1)).ravel() / np.maximum(total, 1) <= FILTER["max_mito_pct"]
gkeep = ~np.isin(np.asarray(genes, dtype=object), GENES_REMOVED)   # genes detected in too few cells (of all slices)
counts, genes = counts[keep][:, gkeep].tocsr(), [g for g, k in zip(genes, gkeep) if k]
xy, cell_ids = xy[keep], np.asarray(cell_ids)[keep]'''
    out.append(_code(read + '''
depth_um = sa.importer.epithelial_depth(counts, genes, xy)
cells = pd.DataFrame({"cell": np.arange(len(cell_ids)), "cell_id": cell_ids, "x": xy[:, 0], "y": xy[:, 1], "depth_um": depth_um})
print(f"{len(cells):,} cells, {len(genes)} genes")''', [_text(f"{n:,} cells, {len(s.panel)} genes")]))

    # ---- 3. registration
    A = importer.placement(s.sample_id)
    by = meta.get("registration_by") or ("file" if src.get("registration") else "none")
    fixes = meta.get("registration_fixes") or []
    if not fixes:
        out.append(_md(f"## {sec()}. Cells on the image\nThe registration: an affine from cell coordinates (um) to the "
                      f"image's full-resolution pixels, as used by the app (made by: {by}). The same numbers are in "
                      "**Export > Registration (CSV)**."))
        out.append(_code(f'''
IMAGE = {_py(meta.get("he_path"))}
A = np.array({_py(np.round(A, 10).tolist())})   # px = A[0,0]*x + A[0,1]*y + A[0,2];  py = A[1,0]*x + A[1,1]*y + A[1,2]
pxy = xy @ A[:, :2].T + A[:, 2]
cells["px"], cells["py"] = pxy[:, 0], pxy[:, 1]
print(f"{{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f}} um per image pixel")''',
                        [_text(f"{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f} um per image pixel")]))
    else:  # the registration as imported, then its corrections in order (local shifts and moved pieces are no affine)
        A0 = np.asarray(meta.get("placement_affine_import") or A, float)
        ids = s.cells.cell_id.to_numpy()  # a piece's cells as rows of the notebook's table (the same order)
        def positions_fix():  # a cell positions file: the rows it places and where (the file is in the sample)
            pos = importer._positions_of(s.sample_id)
            if pos is None:
                return {"positions": None}
            at = pd.Index(pos.cell_id.astype(str)).get_indexer(pd.Index(ids).astype(str))
            rows = np.flatnonzero(at >= 0)
            xy_ = pos[["px", "py"]].to_numpy(float)[at[rows]]
            return {"positions": {"rows_packed": _pack(rows, np.int32), "px_packed": _pack(xy_[:, 0], np.float64),
                                  "py_packed": _pack(xy_[:, 1], np.float64)}}
        packed = [{"piece": dict({k: x for k, x in f["piece"].items() if k not in ("cells", "cell_ids")},
                                 cells_packed=_pack(importer.piece_rows(f["piece"], ids), np.int32))} if "piece" in f
                  else positions_fix() if "positions" in f
                  else {"correction_px": f["correction_px"]} for f in fixes]
        off = _replay_error(s.cells, meta)
        out.append(_md(f"## {sec()}. Cells on the image\nThe registration as imported: an affine from cell coordinates "
                      "(um) to the image's full-resolution pixels; then the corrections made in the app, in order "
                      f"({len(fixes)}: made by {by}). A correction moves the cells by an affine plus local shifts that "
                      "fade with distance (stretched or torn tissue); a piece move turns and shifts the cells of one torn "
                      "piece. The resulting positions are in **Export > Registration > Cell positions (CSV)**; the "
                      "affine alone (**Save registration**) leaves the local parts out." +
                      ("" if off < 0.01 else f" Replayed here, they differ from the app's positions by up to {off:.1f} "
                       "image pixels; the exported cell positions are the app's.")))
        out.append(_code(f'''
IMAGE = {_py(meta.get("he_path"))}
A_IMPORT = np.array({_py(np.round(A0, 10).tolist())})   # px = A[0,0]*x + A[0,1]*y + A[0,2];  py = A[1,0]*x + A[1,1]*y + A[1,2]
FIXES = {_py(packed)}

def field(p, C, S, sig, chunk=20000):     # Gaussian-weighted average at p of the shifts S at points C, fading out far away
    out = np.zeros((len(p), 2))
    for k in range(0, len(p), chunk):
        wt = np.exp(-((p[k:k + chunk, None, :] - C[None]) ** 2).sum(-1) / (2 * sig ** 2))
        out[k:k + chunk] = (wt @ S) / (wt.sum(1) + 0.1)[:, None]
    return out

def displacement(p, fix):                  # large shifts first, then A (p - origin) + b and the local shifts
    pre, loc = fix.get("pre"), fix.get("local")
    d0 = field(p, np.asarray(pre["points_px"]), np.asarray(pre["shift_px"]), float(pre["sigma_px"])) if pre else np.zeros_like(p)
    q = p + d0
    out = d0 + (q - np.asarray(fix["origin"])) @ np.asarray(fix["A"]).T + np.asarray(fix["b"])
    if loc and loc["points_px"]:
        out += field(q, np.asarray(loc["points_px"]), np.asarray(loc["shift_px"]), float(loc["sigma_px"]))
    return out

pxy = xy @ A_IMPORT[:, :2].T + A_IMPORT[:, 2]
for f in FIXES:
    if "correction_px" in f:
        pxy = pxy + displacement(pxy, f["correction_px"])
    elif "positions" in f:                  # a cell positions file (the registration): its cells where it put them
        if f["positions"]:
            rows = unpack(f["positions"]["rows_packed"], np.int32)
            pxy[rows, 0], pxy[rows, 1] = unpack(f["positions"]["px_packed"], np.float64), unpack(f["positions"]["py_packed"], np.float64)
    else:                                   # a torn piece: q = M (p - center) + center + t
        pc = f["piece"]; sel = np.isin(np.arange(len(pxy)), unpack(pc["cells_packed"], np.int32))
        M, c, t = (np.asarray(pc[k], float) for k in ("M", "center_px", "t_px"))
        pxy[sel] = (pxy[sel] - c) @ M.T + c + t
cells["px"], cells["py"] = pxy[:, 0], pxy[:, 1]
print(f"{{len(FIXES)}} correction(s) applied")''', [_text(f"{len(fixes)} correction(s) applied")]))

    return out


# ---------------- a project: all slices together ----------------
PROJECT_SLICE = '''
SLICES[{sid}] = dict(cells=cells, counts=counts, genes=list(genes), xy=xy, cell_ids=np.asarray(cell_ids), depth_um=depth_um)
print({sid}, f"{{len(cells):,}} cells")'''

PROJECT_ANNOTATE = '''
ORDER = {order}                          # the slices pooled, in the order the app pooled them
SCHEME = {scheme}
TYPES = SCHEME["types"]; T = len(TYPES)
LINEAGE = {lineage}
genes = SLICES[ORDER[0]]["genes"]
assert all(SLICES[x]["genes"] == genes for x in ORDER), "the slices must share one gene panel"
OFF = np.cumsum([0] + [SLICES[x]["counts"].shape[0] for x in ORDER])
pooled = sp.vstack([SLICES[x]["counts"] for x in ORDER]).tocsr()
SPILL = ({{"xy": np.vstack([SLICES[x]["xy"] for x in ORDER]), "group": np.repeat(np.arange(len(ORDER)), np.diff(OFF))}}
         if SCHEME.get("spillover") else None)   # transcripts spilled over from neighbors (of the same slice) left out
res = sa.markers.annotate(pooled, genes, SCHEME["tree"], resolution=SCHEME["resolution"], knn=SCHEME["knn"], spill=SPILL)  # one marker step
index = {{t: i for i, t in enumerate(TYPES)}}
pooled_labels = sa.markers.training_labels(res, TYPES)
joint_model = sa.model.train_pseudo(pooled, pooled_labels, T, keep=sa.markers.gene_columns(SCHEME["tree"], genes))  # one model

def slice_res(res, o, n):                # the pooled marker step's part of one slice (rows made local)
    members = {{}}
    for k, rows in res["members"].items():
        r = np.asarray(rows); r = r[(r >= o) & (r < o + n)] - o
        if len(r):
            members[k] = r
    return {{"labels": res["labels"][o:o + n], "cluster_of": res["cluster_of"][o:o + n], "members": members,
            "records": [dict(r, n_cells=int(len(members[r["id"]]))) for r in res["records"] if r["id"] in members],
            "nodes": res["nodes"], "knn": res.get("knn")}}

import copy
def marker_tree(ops):                    # the marker list with the edits of a slice's re-clusterings applied
    tree = copy.deepcopy(SCHEME["tree"])
    for op in ops:
        for t, mk in (op.get("markers") or {{}}).items():
            sa.markers._find(tree, t)["markers"] = sa.markers.parse(f"{{t}}: {{mk}}", genes)["tree"][0]["markers"]
    return tree
print(f"{{pooled.shape[0]:,}} cells of {{len(ORDER)}} slices; marker step and expression model done")'''

PROJECT_SLICE_LABELS = '''
sid = {sid}
STATE = {state}
counts, xy, cells, depth_um = (SLICES[sid][k] for k in ("counts", "xy", "cells", "depth_um"))
i = ORDER.index(sid)
for o in STATE["regions"]["ops"]:
    o["cells"] = unpack(o.pop("cells_packed"), np.int32).tolist()
corr = {{}}
for c in STATE["corrections"]:           # hand-labeled cells (a later correction of the same cell wins)
    for cell in unpack(c["cells_packed"], np.int32):
        corr[int(cell)] = TYPES.index(c["label"])
lock_rows = np.array(list(corr), int); lock_labels = np.array(list(corr.values()), int)
r = slice_res(res, OFF[i], OFF[i + 1] - OFF[i])
for j, op in enumerate(STATE["marker_ops"]):   # this slice's own re-clusterings
    r = sa.markers.rerun(r, counts, genes, marker_tree(STATE["marker_ops"][:j + 1]), op["target"], op["resolution"], prefix=f"{{op['id']}}:", knn=SCHEME["knn"],
                         spill={{"xy": xy}} if SCHEME.get("spillover") else None)
marker_labels = sa.markers.training_labels(r, TYPES)
model = joint_model if not STATE["marker_ops"] else sa.model.train_pseudo(
    counts, marker_labels, T, keep=sa.markers.gene_columns(marker_tree(STATE["marker_ops"]), genes))
if len(lock_rows):
    rX, ry, _ = sa.model.pseudo_pool(counts, marker_labels, T)
    model = sa.model.finetune(model, counts[lock_rows], lock_labels, rX, ry)
logp = sa.markers.restrict(sa.model.predict(model, counts), r, marker_tree(STATE["marker_ops"]), TYPES, lock_labels,
                           counts=counts, panel=genes, xy=xy, neighbors=bool(SCHEME.get("neighbors")),
                           spill=bool(SCHEME.get("spillover")), legacy=SCHEME.get("neighbors") is None)
base_pred, _ = sa.constraints.posterior(logp, [], None, TYPES)          # expression only
region = sa.regions.assign(cells, STATE["regions"])                    # drawn regions
nn = cKDTree(xy).query(xy, k=sa.constraints.NEIGHBORS_K + 1)[1][:, 1:]
neighbor_frac = lambda names: np.isin(base_pred, [TYPES.index(t) for t in names])[nn].mean(1)
ctx = {{"region": region, "depth_um": depth_um, "cell": cells.cell.to_numpy(), "neighbor_frac": neighbor_frac}}
pred, conf, _, _ = sa.constraints.posterior(logp, STATE["constraints"], ctx, TYPES, top=3)   # priors
pred[lock_rows], conf[lock_rows] = lock_labels, 1.0                     # hand corrections
labels = np.array(TYPES, dtype=object)[pred]
ORIGINAL = np.array({tnames}, dtype=object)[unpack("{tcodes}", np.uint16)]
same = labels == ORIGINAL
SLICES[sid].update(labels=labels, pred=pred, conf=conf, region=region)
print(sid, "labels identical to the app's" if same.all() else f"{{same.mean():.2%}} of the labels as in the app "
      f"({{(~same).sum():,}} cells differ)")'''

PROJECT_EMBED = '''
EORDER = {order}
N_PCS, K, RESOLUTION, MIN_DIST, HARMONY = {n_pcs}, {k}, {resolution}, {min_dist}, {harmony}
FIT_MAX = {fit_max}                      # the cells the components were fitted on (None: every cell)
EOFF = np.cumsum([0] + [SLICES[x]["counts"].shape[0] for x in EORDER])
Z = sa.markers.pca(sa.model.lognorm(sp.vstack([SLICES[x]["counts"] for x in EORDER]).tocsr()), N_PCS,
                   fit_max=FIT_MAX)[:, ::-1].copy()
if HARMONY:                              # slice-to-slice differences removed from the PCA, as the app did
    batch = np.repeat(EORDER, np.diff(EOFF))
    if {numpy_harmony}:                  # the app ran on Windows: its NumPy Harmony
        Z = sa.harmony.harmonize(Z, batch, random_state=0)
    else:
        import harmonypy
        ho = harmonypy.run_harmony(Z, pd.DataFrame({{"slice": batch}}), ["slice"], verbose=False, random_state=0)
        Z = np.ascontiguousarray(np.asarray(ho.Z_corr).T if np.asarray(ho.Z_corr).shape[0] == Z.shape[1] else ho.Z_corr, dtype=np.float32)
dist, nbr = sa.markers.knn(Z, K)
leiden = sa.markers._leiden(sa.markers.knn_graph(dist, nbr), nbr, RESOLUTION, 0)
print(leiden.max() + 1, "Leiden clusters of all slices; identical to the app's:",
      bool((leiden == unpack("{leiden}", np.int32)).all()))
import umap
n = len(Z)
init = Z[:, :2].astype(np.float64); init = 10 * (init - init.mean(0)) / max(float(np.abs(init - init.mean(0)).max()), 1e-9)
umap_xy = umap.UMAP(n_neighbors=K + 1, min_dist=MIN_DIST, init=init, n_jobs=-1, low_memory=True,
                    **(dict(compatibility_layout=True) if "compatibility_layout" in umap.UMAP.__init__.__code__.co_varnames else dict()),  # umap 0.6: the 0.5 layout
                    precomputed_knn=(np.c_[np.arange(n), nbr], np.c_[np.zeros(n, np.float32), dist].astype(np.float32), None)
                    ).fit_transform(Z)
for j, x in enumerate(EORDER):
    SLICES[x].update(leiden=leiden[EOFF[j]:EOFF[j + 1]], umap=umap_xy[EOFF[j]:EOFF[j + 1]])
plot(umap_xy[:, 0], umap_xy[:, 1], np.repeat(EORDER, np.diff(EOFF)).astype(object), "UMAP of all slices, by slice", invert=False)'''

PROJECT_NICHES = '''
from sklearn.cluster import KMeans
K_NICHE, N_NICHES, IDENTITY, N_PROGRAMS = {k}, {n_niches}, {identity}, {n_programs}
if IDENTITY == "nmf":                    # NMF programs fitted on the pooled cells
    from sklearn.decomposition import NMF
    Xn = sa.model.lognorm(pooled).astype(np.float32).tocsr()
    fit = np.sort(np.random.default_rng(0).choice(Xn.shape[0], min(60000, Xn.shape[0]), replace=False))
    nmf = NMF(n_components=N_PROGRAMS, init="nndsvda", max_iter=300, random_state=0, tol=1e-4).fit(Xn[fit])
    W = np.vstack([nmf.transform(Xn[i:i + 50000]) for i in range(0, Xn.shape[0], 50000)]).astype(np.float32)
    identity = W / np.maximum(W.sum(1, keepdims=True), 1e-9)
else:                                    # the cell types of each slice
    pred_all = np.concatenate([SLICES[x]["pred"] for x in ORDER])
    identity = np.zeros((len(pred_all), T), np.float32); identity[np.arange(len(pred_all)), pred_all] = 1
avg = np.zeros_like(identity)
for i, x in enumerate(ORDER):            # averaged over each cell's neighbors in its own slice
    nn_ = cKDTree(SLICES[x]["xy"]).query(SLICES[x]["xy"], k=K_NICHE + 1)[1]
    part = identity[OFF[i]:OFF[i + 1]]
    avg[OFF[i]:OFF[i + 1]] = np.vstack([part[nn_[j:j + 20000]].mean(1) for j in range(0, len(nn_), 20000)])
km = KMeans(n_clusters=N_NICHES, n_init=4, random_state=0).fit(avg)
order_ = np.argsort(-np.bincount(km.labels_, minlength=N_NICHES))
niche = np.argsort(order_)[km.labels_]
NAMES = {names}
print("niches of all slices identical to the app's:", bool((niche == unpack("{labels}", np.int16)).all()))
for i, x in enumerate(ORDER):
    SLICES[x]["niche"] = np.array(NAMES, dtype=object)[niche[OFF[i]:OFF[i + 1]]]'''


def joint_run_sessions(pid, run=None):
    """{sample_id: session id} of one joint annotation of the project (Annotate with All slices together), found as
    the app's notebook dialog lists them: run (its time, scheme["joint"]["run"]), or the latest run that covers every
    slice of the project (else the latest run); {} when the project has none."""
    from . import projects
    from .session import list_sessions
    p, ss, runs = projects.load(pid), list_sessions(), {}
    for sid, smp in (zip(ss.session_id, ss.sample_id) if len(ss) else ()):
        if smp not in p["samples"]:
            continue
        try:
            with open(config.SESSIONS / sid / "session.json", encoding="utf-8") as fh:
                j = (json.load(fh).get("scheme") or {}).get("joint") or {}
        except (OSError, ValueError):
            continue
        if j.get("project") == pid and j.get("run"):
            runs.setdefault(j["run"], {})[smp] = sid
    if run is not None:
        return runs.get(run, {})
    full = [r for r, sess in runs.items() if set(sess) >= set(p["samples"])]
    return runs[max(full or runs)] if runs else {}


def build_project(pid, sessions=None):
    """A notebook of a whole project annotated jointly (Annotate with All slices together): every slice read from its
    own files and placed on its image, the joint marker step and expression model on the pooled cells, then each
    slice's own regions, priors, corrections and re-clusterings, checked against the app slice by slice; the joint PCA /
    UMAP (Harmony) and niches when they were made for the project. Returns its path."""
    import os
    from . import joint, projects
    from .session import Session
    p = projects.load(pid)
    # the slices of one joint run only (a slice annotated again on its own, or added to the project after the run,
    # is left out rather than taken with an annotation of its own); no sessions given: the latest run of all slices
    sessions = sessions or joint_run_sessions(pid)
    latest = joint.latest_sessions(pid, sessions)
    S = {x: Session(latest[x]) for x in p["samples"] if x in latest and x in sessions}
    runs = {(s.scheme.get("joint") or {}).get("run") for s in S.values()}
    if len(S) < 2 or None in runs or len(runs) != 1:
        raise ValueError("the slices' annotations are not one joint annotation: annotate with 'All slices together' "
                         "(or export each slice's notebook)")
    af = joint._dir(pid) / "annotation.json"
    ann = json.load(open(af)) if af.exists() else {}
    order = list(ann["sessions"]) if ann.get("run") in runs and set(ann.get("sessions", {})) == set(S) else list(S)
    S = {x: S[x] for x in order}
    s0 = S[order[0]]
    sch = s0.scheme
    nb, count = [], iter(range(1, 50))
    sec = lambda: next(count)
    ver = {x: s.state()["version"] for x, s in S.items()}
    nb.append(_md(f"""
# STAnnotator analysis: project `{p.get('name', pid)}`, {len(S)} slices annotated together

Exported {time.strftime('%Y-%m-%d %H:%M')} (STAnnotator build `{config.build()}`). Slices, with the annotation and version
of each: {", ".join(f"`{x}` ({S[x].id}, v{ver[x]})" for x in order)}. The notebook reads every slice from its own files,
places it on its image, runs the app's joint annotation on all cells together (one marker step, one expression model),
then each slice's own regions, priors, hand corrections and re-clusterings, and checks the labels slice by slice;
then the PCA / UMAP and the niches of all slices when they were made for the project. It needs the data files and
common Python packages (`numpy pandas scipy scikit-learn python-igraph pynndescent umap-learn h5py matplotlib`, and
`harmonypy` when the embedding used Harmony), installed by the first code cell; the app's own algorithm is carried in the folded cell of section 1.
The outputs below are the results at export time."""))
    extra = (["pillow"] if any((S[x].meta.get("source") or {}).get("segmentation") for x in order) else []) + (
        ["harmonypy"] if os.name != "nt" and any((embedding.current(x) or {}).get("harmony") for x in order) else [])
    nb += _package_cells(sec, extra)
    nb.append(_code("SLICES = {}   # every slice's cells, counts and results, by sample"))
    for x in order:
        k = sec()
        sub = iter(range(1, 20))
        nb.append(_md(f"## {k}. Slice `{x}`"))
        nb += _data_cells(S[x], lambda: f"{k}.{next(sub)}")
        nb.append(_code(PROJECT_SLICE.format(sid=_py(x)), [_text(f"{x} {len(S[x].cells):,} cells")]))

    # ---- the joint annotation
    scheme = {k: sch.get(k) for k in ("kind", "text", "tree", "types", "resolution", "top", "lineage_of_type", "leaf_types",
                                      "neighbors", "spillover")}
    scheme["knn"] = s0.knn_mode()
    lineage = {t: s0.lineage_of_type.get(t, t) for t in s0.types}
    n_all = sum(len(S[x].cells) for x in order)
    nb.append(_md(f"## {sec()}. The joint annotation\nThe marker list's hierarchy clustered on all {n_all:,} cells of the "
                  "slices together (Leiden per level, scored by the markers), one expression model trained on those labels."))
    nb.append(_code(PROJECT_ANNOTATE.format(order=_py(order), scheme=_py(scheme), lineage=_py(lineage)),
                    [_text(f"{n_all:,} cells of {len(order)} slices; marker step and expression model done")]))
    for x in order:
        s, st = S[x], S[x].state()
        rows = {"regions": dict(st["regions"], ops=[dict({kk: vv for kk, vv in o.items() if kk != "cells"},
                                                         cells_packed=_pack(o.get("cells", []), np.int32))
                                                    for o in st["regions"].get("ops", [])]),
                "constraints": st.get("constraints", []),
                "corrections": [dict({kk: vv for kk, vv in c.items() if kk != "cells"}, cells_packed=_pack(c["cells"], np.int32))
                                for c in st.get("corrections", [])],
                "marker_ops": st.get("marker_step", {}).get("ops", [])}
        labels = s.labels(st["version"])
        tnames = sorted(set(labels))
        priors = [c.get("text") or c["id"] for c in st["constraints"]]
        nb.append(_md(f"### Slice `{x}`: its regions and priors\n{len(rows['regions']['ops'])} region polygons, "
                      f"{len(priors)} priors ({'; '.join(priors) or 'none'}), {len(rows['corrections'])} hand corrections "
                      f"and {len(rows['marker_ops'])} re-clusterings of its own, as in version {st['version']}."))
        nb.append(_code(PROJECT_SLICE_LABELS.format(sid=_py(x), state=_py(rows), tnames=_py(tnames),
                                                    tcodes=_pack(np.searchsorted(tnames, labels), np.uint16)),
                        [_text(f"{x} labels identical to the app's")]))

    # ---- the joint embedding
    ems = [embedding.current(x) or {} for x in order]
    if all(e.get("project") == pid for e in ems) and len({e.get("run_id") for e in ems}) == 1:
        em = ems[0]
        eorder = [x for x in em.get("slices") or order if x in S]
        lab = np.concatenate([embedding.leiden(x) for x in eorder])
        nb.append(_md(f"## {sec()}. PCA, neighbors, Leiden and UMAP of all slices\nPCA of the pooled log-normalized "
                      f"expression ({em['n_pcs']} components)" + (", Harmony across slices" if em.get("harmony") else "")
                      + f", {em['n_neighbors']} nearest neighbors, Leiden (resolution {em['resolution']}: "
                      f"{em['n_clusters']} clusters) and UMAP on the same neighbors."))
        nb.append(_code(PROJECT_EMBED.format(order=_py(eorder), n_pcs=em["n_pcs"], k=em["n_neighbors"],
                                             fit_max=em.get("pca_fit_cells") or None,
                                             resolution=em["resolution"], min_dist=em["min_dist"],
                                             harmony=bool(em.get("harmony")), numpy_harmony=os.name == "nt",
                                             leiden=_pack(lab, np.int32)),
                        [_text(f"{em['n_clusters']} Leiden clusters of all slices; identical to the app's: True")]))
    # ---- the joint niches
    nms = [niches.current(S[x])[0] or {} for x in order]
    if all(m.get("project") == pid for m in nms) and len({m.get("run_id") for m in nms}) == 1:
        nm = nms[0]
        lab = np.concatenate([niches.current(S[x])[1] for x in order])
        nb.append(_md(f"## {sec()}. Niches of all slices\nEach cell's identity (" + ("NMF gene programs fitted on the "
                      "pooled cells" if nm["identity"] == "nmf" else "its cell type") + f") averaged over the cell and its "
                      f"{nm['k']} nearest cells in its own slice, then one k-means over all slices into {nm['n_niches']} niches."))
        nb.append(_code(PROJECT_NICHES.format(k=int(nm["k"]), n_niches=int(nm["n_niches"]), identity=_py(nm["identity"]),
                                              n_programs=int(nm.get("n_programs") or 20), names=_py(niches.niche_names(S[order[0]])),
                                              labels=_pack(lab, np.int16)),
                        [_text("niches of all slices identical to the app's: True")]))
    # ---- the cell table
    nb.append(_md(f"## {sec()}. Cell table\nOne row per cell of every slice, written to `cells.csv`."))
    nb.append(_code('''
parts = []
for x in ORDER:
    d = SLICES[x]
    t = pd.DataFrame({"slice": x, "cell_id": d["cell_ids"], "region": np.where(d["region"] == "unknown", "", d["region"]),
                      "cell_type": d["labels"], "lineage": [LINEAGE[t] for t in d["labels"]], "confidence": d["conf"].round(4)})
    for k in ("niche", "leiden"):
        if k in d:
            t[k] = d[k] if k == "niche" else ["L" + str(i + 1) for i in d[k]]
    if "umap" in d:
        t["umap_1"], t["umap_2"] = d["umap"][:, 0].round(4), d["umap"][:, 1].round(4)
    t["x_um"], t["y_um"] = d["xy"][:, 0].round(2), d["xy"][:, 1].round(2)
    parts.append(t)
table = pd.concat(parts, ignore_index=True)
table.to_csv("cells.csv", index=False)
table.groupby("slice").cell_type.value_counts().unstack(0).fillna(0).astype(int).head(30)'''))

    doc = {"nbformat": 4, "nbformat_minor": 5, "cells": nb,
           "metadata": {"kernelspec": {"name": "python3", "display_name": "Python 3", "language": "python"},
                        "language_info": {"name": "python"}}}
    for i, c in enumerate(doc["cells"]):
        c["id"] = f"c{i:02d}"
    d = joint._dir(pid) / "exports"
    d.mkdir(parents=True, exist_ok=True)
    out = d / f"{pid}_all_slices_analysis.ipynb"
    with open(out, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=1)
    return str(out)


def _replay_error(cells, meta):
    """How far (image pixels, at most) replaying the corrections from the import affine lands from the stored positions."""
    from .importer import replay_placement
    _, p = replay_placement(meta, cells[["x", "y"]].to_numpy(float), cells.cell_id.to_numpy())
    return float(np.abs(p - cells[["px", "py"]].to_numpy(float)).max()) if len(p) else 0.0


def positions_csv(sample_id):
    """Every cell's place on the image as the app uses it now, local corrections and moved pieces included: cell_id,
    x_um, y_um (as imported) and px, py (the image's full-resolution pixels)."""
    c = pd.read_parquet(config.SAMPLES / sample_id / "cells.parquet", columns=["cell_id", "x", "y", "px", "py"])
    return pd.DataFrame({"cell_id": c.cell_id, "x_um": c.x.round(3), "y_um": c.y.round(3),
                         "px": c.px.round(3), "py": c.py.round(3)}).to_csv(index=False)


def registration_csv(sample_id):
    """The affine that places the sample's cells on its image, as an a11..ty CSV (importable as a registration)."""
    from . import importer
    m = importer._meta(sample_id)
    A = importer.placement(sample_id)
    return pd.DataFrame([{"sample": sample_id, "a11": A[0, 0], "a12": A[0, 1], "a21": A[1, 0], "a22": A[1, 1],
                          "tx": A[0, 2], "ty": A[1, 2], "um_per_px": float(1 / np.sqrt(abs(np.linalg.det(A[:, :2])))),
                          "made_by": m.get("registration_by") or "file", "image": m.get("he_path")}]).to_csv(index=False)
