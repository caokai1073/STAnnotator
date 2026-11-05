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
# visium and stereoseq before importer, which imports them
CARRIED = ["timing", "vocab", "markers", "model", "constraints", "regions", "existing", "visium", "stereoseq", "importer", "scref"]
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
matplotlib{" pillow" if (meta.get("source") or {}).get("segmentation") else ""}`. The outputs below are the results at export time, so it can be read without running it; running it
from the top repeats the analysis and checks that the labels, Leiden clusters and niches come out the same.
The annotation algorithm of the app is carried as source code in the folded cell of section 1.

{knn_text}
"""))
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
    unit = {"visium_hd": "cells built from the 2 um bins by the cell outlines" if src.get("unit") == "cells" else
            f"bins of {src.get('bin_um')} um", "visium": "spots", "stereoseq": "cells (cell bins)" if src.get("unit") == "cells"
            else f"bins of {src.get('bin_um')} um"}.get(src.get("kind"))
    nb.append(_md(f"## {sec()}. The cells\nRead from the original files (edit the paths if they moved)" +
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
parts = [sa.importer.read_h5ad(p, SECTION_COLUMN, SECTION) for p in H5AD]
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
ids, genes, X = sa.importer._flat_counts(KIND, FILES["counts"])
cell_ids, xy = sa.importer._flat_cells(KIND, FILES["cells"])
pos = pd.Series(np.arange(len(ids)), index=ids)
have = pd.Index(cell_ids).isin(pos.index)   # cells of the metadata with counts
cell_ids, xy = cell_ids[have], xy[have]
counts = X[pos[cell_ids].to_numpy()].tocsr()'''
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
    flt = meta.get("filter")
    if flt:  # the sample's QC filter (qcfilter.py), applied to the cells as imported
        nb.append(_md(f"The sample was filtered in the app ({', '.join(f'{k} {x:g}' for k, x in flt['params'].items())}): "
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
    nb.append(_code(read + '''
depth_um = sa.importer.epithelial_depth(counts, genes, xy)
cells = pd.DataFrame({"cell": np.arange(len(cell_ids)), "cell_id": cell_ids, "x": xy[:, 0], "y": xy[:, 1], "depth_um": depth_um})
print(f"{len(cells):,} cells, {len(genes)} genes")''', [_text(f"{n:,} cells, {len(s.panel)} genes")]))

    # ---- 3. registration
    A = importer.placement(s.sample_id)
    by = meta.get("registration_by") or ("file" if src.get("registration") else "none")
    fixes = meta.get("registration_fixes") or []
    if not fixes:
        nb.append(_md(f"## {sec()}. Cells on the image\nThe registration: an affine from cell coordinates (um) to the "
                      f"image's full-resolution pixels, as used by the app (made by: {by}). The same numbers are in "
                      "**Export > Registration (CSV)**."))
        nb.append(_code(f'''
IMAGE = {_py(meta.get("he_path"))}
A = np.array({_py(np.round(A, 10).tolist())})   # px = A[0,0]*x + A[0,1]*y + A[0,2];  py = A[1,0]*x + A[1,1]*y + A[1,2]
pxy = xy @ A[:, :2].T + A[:, 2]
cells["px"], cells["py"] = pxy[:, 0], pxy[:, 1]
print(f"{{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f}} um per image pixel")''',
                        [_text(f"{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f} um per image pixel")]))
    else:  # the registration as imported, then its corrections in order (local shifts and moved pieces are no affine)
        A0 = np.asarray(meta.get("placement_affine_import") or A, float)
        ids = s.cells.cell_id.to_numpy()  # a piece's cells as rows of the notebook's table (the same order)
        packed = [{"piece": dict({k: x for k, x in f["piece"].items() if k not in ("cells", "cell_ids")},
                                 cells_packed=_pack(importer.piece_rows(f["piece"], ids), np.int32))} if "piece" in f
                  else {"correction_px": f["correction_px"]} for f in fixes]
        off = _replay_error(s.cells, meta)
        nb.append(_md(f"## {sec()}. Cells on the image\nThe registration as imported: an affine from cell coordinates "
                      "(um) to the image's full-resolution pixels; then the corrections made in the app, in order "
                      f"({len(fixes)}: made by {by}). A correction moves the cells by an affine plus local shifts that "
                      "fade with distance (stretched or torn tissue); a piece move turns and shifts the cells of one torn "
                      "piece. The resulting positions are in **Export > Registration > Cell positions (CSV)**; the "
                      "affine alone (**Save registration**) leaves the local parts out." +
                      ("" if off < 0.01 else f" Replayed here, they differ from the app's positions by up to {off:.1f} "
                       "image pixels; the exported cell positions are the app's.")))
        nb.append(_code(f'''
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
    else:                                   # a torn piece: q = M (p - center) + center + t
        pc = f["piece"]; sel = np.isin(np.arange(len(pxy)), unpack(pc["cells_packed"], np.int32))
        M, c, t = (np.asarray(pc[k], float) for k in ("M", "center_px", "t_px"))
        pxy[sel] = (pxy[sel] - c) @ M.T + c + t
cells["px"], cells["py"] = pxy[:, 0], pxy[:, 1]
print(f"{{len(FIXES)}} correction(s) applied")''', [_text(f"{len(fixes)} correction(s) applied")]))

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
res = sa.markers.annotate(counts, genes, SCHEME["tree"], resolution=SCHEME["resolution"], knn=SCHEME["knn"])   # marker step
for j, op in enumerate(STATE["marker_ops"]):
    res = sa.markers.rerun(res, counts, genes, marker_tree(STATE["marker_ops"][:j + 1]), op["target"], op["resolution"], prefix=f"{op['id']}:", knn=SCHEME["knn"])
index = {t: i for i, t in enumerate(TYPES)}
marker_labels = np.array([index.get(x, -1) for x in res["labels"]])   # -1: left at a lineage / parent level
model = sa.model.train_pseudo(counts, marker_labels, T, keep=sa.markers.gene_columns(SCHEME["tree"], genes))   # expression model
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
for name, k in {_py(nmeta.get('k_by_region') or {})}.items():   # k per drawn region
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
NAMES = {_py(names)}
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
