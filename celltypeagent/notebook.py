"""Export an annotation as a Jupyter notebook: the code of every step with the parameters this annotation used, and
its results already filled in (tables, figures), so the analysis can be read without running it and repeated by
running it (with the celltypeagent package and the sample's data folder).

Steps: data (and how to re-import it from the original files), placement on the image, annotation (method and its
inputs, then the regions, priors, corrections and re-clusterings replayed with Session.apply_state; a fingerprint of
the labels checks the replay), PCA / UMAP / Leiden, niches, the spatial map and the cell table. UMAP runs on all
cores, so its coordinates are not bit-for-bit identical on a rerun; everything else is.
"""
import base64, io, json, time, zlib
import numpy as np, pandas as pd
from . import config, embedding, niches, render


def _png(fig):
    import matplotlib.pyplot as plt
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=90, bbox_inches="tight")
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def _md(text):
    return {"cell_type": "markdown", "metadata": {}, "source": text.strip("\n")}


def _code(src, outputs=()):
    return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": list(outputs),
            "source": src.strip("\n")}


def _text(t):
    return {"output_type": "stream", "name": "stdout", "text": t if t.endswith("\n") else t + "\n"}


def _table(df):
    return {"output_type": "execute_result", "execution_count": None, "metadata": {},
            "data": {"text/plain": df.to_string(), "text/html": df.to_html()}}


def _image(b64):
    return {"output_type": "display_data", "metadata": {}, "data": {"image/png": b64, "text/plain": "<Figure>"}}


def _scatter(x, y, keys, title, invert=True, max_n=200000):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    idx = rng.choice(len(x), min(len(x), max_n), replace=False)
    names = pd.Series(keys).value_counts().index.tolist()
    pal = (render.DISTINCT * 4)
    fig, ax = plt.subplots(figsize=(8, 7))
    for i, n in enumerate(names[:40]):
        m = keys[idx] == n
        ax.scatter(x[idx][m], y[idx][m], s=0.3, c=pal[i % len(pal)], label=n, linewidths=0, rasterized=True)
    if invert:
        ax.invert_yaxis()
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=10)
    ax.legend(markerscale=12, fontsize=6, loc="center left", bbox_to_anchor=(1, 0.5), frameon=False)
    return _png(fig)


def build(s, version=None):
    """Write the notebook of a version into the session's exports and return its path."""
    from . import importer
    st = s.state(version)
    v = st["version"]
    meta = s.meta
    sch = s.scheme
    labels = s.labels(v)
    nb = []
    count = iter(range(1, 20))
    sec = lambda: next(count)  # section numbers in order
    nb.append(_md(f"""
# Spatial Annotator analysis: sample `{s.sample_id}`, annotation `{s.id}`, version {v}

Exported {time.strftime('%Y-%m-%d %H:%M')} from Spatial Annotator build `{config.build()}`. Every step below is the
code of what was done, with the parameters used; the outputs are the results at export time, so the notebook can be
read without running it. To repeat the analysis, run it from the top with the `celltypeagent` package installed
(the app's Python environment, e.g. `~/Applications/STCellTypeAgent/.venv`) and `DATA` set to a data folder that
holds the sample; the annotation is replayed into a new annotation of the sample and checked against the original.
"""))
    nb.append(_code(f"""
import os, json
os.environ["CELLTYPEAGENT_HOME"] = r"{config.HOME}"   # the data folder with the sample (edit if it moved)
import numpy as np, pandas as pd
from celltypeagent import config, importer, embedding, niches
from celltypeagent.session import Session
SAMPLE = "{s.sample_id}"
"""))
    # ---- data
    src = meta.get("source") or {}
    if src.get("kind") == "xenium":
        reimport = f'importer.import_xenium(SAMPLE, r"{src.get("folder")}")'
    elif src.get("kind") == "h5ad":
        reimport = (f'importer.import_sample(SAMPLE, {json.dumps(src.get("h5ad"))}, section_column={json.dumps(src.get("section_column"))}, '
                    f'section={json.dumps(src.get("section"))})')
    else:
        reimport = "# the sample was imported before its source was recorded"
    nb.append(_md(f"## {sec()}. Data\nThe sample's cells, their expression and positions. To rebuild the sample from the "
                  "original files in an empty data folder, uncomment the import."))
    nb.append(_code(f"""
# {reimport}
cells = pd.read_parquet(config.SAMPLES / SAMPLE / "cells.parquet")
genes = (config.SAMPLES / SAMPLE / "genes.txt").read_text().split()
print(f"{{len(cells):,}} cells, {{len(genes)}} genes")
""", [_text(f"{len(s.cells):,} cells, {len(s.panel)} genes")]))
    # ---- placement
    A = importer.placement(s.sample_id)
    by = meta.get("registration_by") or ("file" if src.get("registration") else "none")
    nb.append(_md(f"## {sec()}. Cells on the image\nThe affine that maps cell coordinates (um) to the image's full-resolution "
                  f"pixels (placed by: {by}" + (f"; registration file `{src.get('registration')}`" if src.get("registration") else "")
                  + "). `importer.place_cells` puts the cells there."))
    nb.append(_code(f"""
IMAGE = {json.dumps(meta.get("he_path"))}
A = np.array({json.dumps(np.round(A, 8).tolist())})   # cell um -> image px
# importer.place_cells(SAMPLE, A, IMAGE, {json.dumps(meta.get("image_kind", "H&E"))})
print(f"{{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f}} um per image pixel")
""", [_text(f"{1 / np.sqrt(abs(np.linalg.det(A[:, :2]))):.4f} um per image pixel")]))
    # ---- annotation
    if sch["kind"] == "markers":
        create = (f'MARKERS = """{sch["text"]}"""\n'
                  f's = Session.create(SAMPLE, method="markers", markers_text=MARKERS, resolution={json.dumps(sch.get("resolution"))})')
        how = "marker genes (hierarchical marker list, Leiden clusters per level, marker scores, expression model)"
    elif sch["kind"] == "scref":
        create = f's = Session.create(SAMPLE, method="scref", reference={json.dumps(sch["reference"])})'
        how = "a single-cell reference (CellTypist-like classifier on the shared genes)"
    elif sch["kind"] == "existing":
        e = {k: sch["existing"][k] for k in ("path", "id_column", "label_column", "lineage_column")}
        create = f's = Session.create(SAMPLE, method="existing", existing={json.dumps(e)})'
        how = "an existing annotation loaded from a CSV"
    else:
        create = f'# method {sch["kind"]!r} cannot be replayed outside the lab\ns = Session({json.dumps(s.id)})'
        how = sch["kind"]
    replay = {k: st.get(k) for k in ("regions", "constraints", "corrections", "marker_step")}
    replay["regions"] = dict(replay["regions"], ops=[{k: x for k, x in o.items() if k != "cells"}
                                                    for o in replay["regions"].get("ops", [])])
    comp = pd.Series(labels).value_counts().rename_axis("cell type").to_frame("cells").head(30)
    priors = [c.get("text") or c["id"] for c in st["constraints"]]
    nb.append(_md(f"## {sec()}. Annotation\nMethod: {how}. Then the regions drawn ({len(replay['regions'].get('ops', []))} "
                  f"polygons), the priors ({len(priors)}: {'; '.join(priors) or 'none'}), the hand corrections "
                  f"({len(st.get('corrections', []))}) and the re-clusterings ({len(st.get('marker_step', {}).get('ops', []))}) "
                  "of this version are replayed; the fingerprint of the labels must match the original."))
    types = sorted(set(labels))
    codes = np.searchsorted(types, labels).astype(np.uint16)
    packed = base64.b64encode(zlib.compress(codes.tobytes(), 9)).decode()
    made = s.state(1)["created_at"][:16].replace("T", " ")
    nb.append(_code(f"""
{create}
STATE = json.loads(r'''{json.dumps(replay)}''')
v = s.apply_state(STATE, "replayed from the notebook")
labels = s.labels(v)
import base64, zlib   # the original labels (made {made}), to check the replay
ORIGINAL = np.array({json.dumps(types)}, dtype=object)[np.frombuffer(zlib.decompress(base64.b64decode("{packed}")), np.uint16)]
same = labels == ORIGINAL
print("labels identical to the original" if same.all() else
      f"{{same.mean():.1%}} of the labels as in the original ({{(~same).sum():,}} cells differ: the original was made "
      "with an earlier version of the software)")
pd.Series(labels).value_counts().rename_axis("cell type").to_frame("cells").head(30)
""", [_text("labels identical to the original"), _table(comp)]))
    # ---- spatial map
    px, py = s.cells.px.to_numpy(), s.cells.py.to_numpy()
    lin = s.lineage_labels(v)
    nb.append(_md(f"## {sec()}. Spatial map\nCells at their place on the image, colored by lineage."))
    nb.append(_code("""
import matplotlib.pyplot as plt
lin = s.lineage_labels(v)
fig, ax = plt.subplots(figsize=(8, 7))
for n in pd.Series(lin).value_counts().index[:40]:
    m = lin == n
    ax.scatter(s.cells.px[m], s.cells.py[m], s=0.3, label=n, linewidths=0, rasterized=True)
ax.invert_yaxis(); ax.set_aspect("equal"); ax.legend(markerscale=12, fontsize=6, bbox_to_anchor=(1, 0.5), loc="center left")
""", [_image(_scatter(px, py, lin, f"sample {s.sample_id}, lineages"))]))
    # ---- embedding
    em = embedding.current(s.sample_id)
    if em:
        Y = embedding.coords(s.sample_id, "umap")
        lab = embedding.leiden(s.sample_id)
        nb.append(_md(f"## {sec()}. PCA, UMAP and Leiden\nFrom expression only: sparse PCA ({em['n_pcs']} components), "
                      f"{em['n_neighbors']} nearest neighbors, UMAP (min_dist {em['min_dist']}), Leiden (resolution "
                      f"{em['resolution']}: {em['n_clusters']} clusters)."))
        nb.append(_code(f"""
em = embedding.run(s, n_pcs={em['n_pcs']}, n_neighbors={em['n_neighbors']}, min_dist={em['min_dist']}, resolution={em['resolution']})
umap_xy, leiden = embedding.coords(SAMPLE, "umap"), embedding.leiden(SAMPLE)
print(em["n_clusters"], "Leiden clusters")
fig, ax = plt.subplots(figsize=(8, 7))
for n in pd.Series(labels).value_counts().index[:40]:
    m = labels == n
    ax.scatter(umap_xy[m, 0], umap_xy[m, 1], s=0.3, label=n, linewidths=0, rasterized=True)
ax.legend(markerscale=12, fontsize=6, bbox_to_anchor=(1, 0.5), loc="center left")
""", [_text(f"{em['n_clusters']} Leiden clusters"),
      _image(_scatter(Y[:, 0], Y[:, 1], np.asarray(labels), "UMAP, cell types", invert=False)),
      _image(_scatter(Y[:, 0], Y[:, 1], np.char.add("L", (lab + 1).astype(str)), "UMAP, Leiden clusters", invert=False))]))
    # ---- niches
    nmeta, nlab = niches.current(s)
    if nmeta:
        summ = niches.summary(s, v)
        tab = pd.DataFrame([{"niche": x["niche"], "name": x["name"], "cells": x["n_cells"],
                             "main cell types": ", ".join(f"{k} {round(100 * p)}%" for k, p in list(x["cell_types"].items())[:3])}
                            for x in summ["niches"]]).set_index("niche")
        names = {k: n for k, n in nmeta.get("names", {}).items()}
        kw = f"k={nmeta['k']}, n_niches={nmeta['n_niches']}, identity={json.dumps(nmeta['identity'])}"
        if nmeta.get("n_programs"):
            kw += f", n_programs={nmeta['n_programs']}"
        if nmeta.get("k_by_region"):
            kw += f", k_by_region={json.dumps(nmeta['k_by_region'])}"
        if nmeta.get("radius_um"):
            kw += f", radius_um={nmeta['radius_um']}"
        nb.append(_md(f"## {sec()}. Niches\nEach cell's identity (" + ("NMF gene programs" if nmeta["identity"] == "nmf" else "cell type")
                      + ") averaged over its neighbors, clustered with k-means."))
        name_line = f"niches.set_names(s, {json.dumps(names)})\n" if names else ""
        nb.append(_code(f"""
niches.run(s, {kw})
{name_line}pd.DataFrame([{{"niche": x["niche"], "name": x["name"], "cells": x["n_cells"],
               "main cell types": ", ".join(f"{{k}} {{round(100 * p)}}%" for k, p in list(x["cell_types"].items())[:3])}}
              for x in niches.summary(s, v)["niches"]]).set_index("niche")
""", [_table(tab), _image(_scatter(px, py, np.array(niches.niche_names(s), dtype=object)[nlab], "niches"))]))
    # ---- table
    table = s.export_table(v)
    nb.append(_md(f"## {sec()}. Cell table\nOne row per cell, as **Export > Cells table (CSV)** writes it."))
    nb.append(_code("""
table = s.export_table(v)
table.head()
""", [_table(table.head())]))
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
