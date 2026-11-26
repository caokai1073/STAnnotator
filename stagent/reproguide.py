"""Export > Reproduction guide for agents: one zip that an outside AI agent (or a person) reads to repeat every step of
an annotation's analysis and redraw its figures without this app's interface.

  REPRODUCE.md    the steps in order, written for an agent: software, input files, cell filter, placement on the
                  image, annotation (its list, settings, regions, priors, corrections, re-clusterings), PCA / UMAP /
                  Leiden, niches, every analysis kept with the annotation (DE genes, top genes, composition, DESeq2,
                  cell-cell communication, co-localization, programs, spatially variable genes, trajectories), the
                  groupings, and how the figures are drawn. For each step: its settings, where the algorithm is in the
                  source (a link at this version) and the results to check a replay against.
  manifest.json   the same as data
  data/           cells.csv (one row per cell, every label), the marker list, regions.geojson, registration.csv,
                  cell_positions.csv, corrections.csv, state.json, results/<kind>.json + .csv, analyses_log.json,
                  colors.json, trajectory_<method>.csv, types_not_found.csv, CHANGELOG.md
  analysis.ipynb  (optional) the notebook: the exact replay of the data, annotation, embedding and niches

The counts and images are not copied (too large): the guide names their files with size and checksum.
"""
import hashlib, importlib, inspect, io, json, os, platform, shutil, sys, time, zipfile
from pathlib import Path
import numpy as np, pandas as pd
from . import __version__, config

REPO = "https://github.com/caokai1073/STAnnotator"


def _ref():
    """The git ref of this code: the commit of the build (its line numbers match), else the version's tag."""
    b = (config.build() or "").split()
    return b[0] if b and len(b[0]) >= 7 and all(c in "0123456789abcdef" for c in b[0]) else f"v{__version__}"
FULL_HASH_MAX = 256 << 20   # files up to this size get a full SHA-256; larger ones of their first and last 8 MiB
PART = 8 << 20
PACKAGES = ("numpy", "pandas", "scipy", "scikit-learn", "igraph", "leidenalg", "pynndescent", "numba", "umap-learn",
            "h5py", "anndata", "matplotlib", "torch", "mlx", "pydeseq2", "statsmodels")
RESULT_TITLES = {"de": "DE genes (Wilcoxon / AUC)", "deseq2": "DE genes (DESeq2 on pseudobulk)", "markers": "Top genes of every group",
                 "composition": "Composition", "spatial": "Spatially variable genes (Moran's I)", "program": "Gene program scores",
                 "run": "Cell-cell communication", "compare": "Cell-cell communication compared across a split",
                 "coloc": "Co-localization (neighborhood enrichment)", "trajectory_spatial": "Trajectory (spatial axis)",
                 "trajectory_expression": "Trajectory (expression pseudotime)"}
RESULT_CODE = {"de": ["analysis.de", "analysis.de_all"], "deseq2": ["analysis.deseq2"], "markers": ["analysis.markers"],
               "composition": ["analysis.crosstab", "analysis.proportions"], "spatial": ["svg.spatial_genes"],
               "program": ["programs.run", "programs.score"], "run": ["communication.run"], "compare": ["communication.compare"],
               "coloc": ["coloc.run"], "trajectory_spatial": ["trajectory.run", "trajectory.spatial_axis", "trajectory.gene_trends"],
               "trajectory_expression": ["trajectory.run", "trajectory.pseudotime", "trajectory.gene_trends"]}
RESULT_MODULE = {"de": "analysis", "deseq2": "analysis", "markers": "analysis", "composition": "analysis", "spatial": "svg",
                 "program": "programs", "run": "communication", "compare": "communication", "coloc": "coloc",
                 "trajectory_spatial": "trajectory", "trajectory_expression": "trajectory"}


# ---------------- helpers ----------------
def _src(qualname):
    """Where a function is in the source of this version: {"function", "file", "line", "url"}."""
    mod, *attrs = qualname.split(".")
    out = {"function": f"stagent.{qualname}", "file": f"stagent/{mod}.py"}
    try:
        obj = importlib.import_module(f"stagent.{mod}")
        for a in attrs:
            obj = getattr(obj, a)
        out["line"] = inspect.getsourcelines(obj)[1]
        out["url"] = f"{REPO}/blob/{_ref()}/stagent/{mod}.py#L{out['line']}"
    except Exception:
        out["url"] = f"{REPO}/blob/{_ref()}/stagent/{mod}.py"
    return out


def _code_lines(names):
    return "\n".join(f"  - [`{x['function']}`]({x['url']}) ({x['file']}" + (f", line {x['line']})" if x.get("line") else ")")
                     for x in (_src(n) for n in names))


def _module_doc(mod):
    try:
        return inspect.cleandoc(importlib.import_module(f"stagent.{mod}").__doc__ or "")
    except Exception:
        return ""


def _file_info(p):
    """Name, size, time and checksum of an input file (a folder: only its name)."""
    p = Path(p)
    out = {"path": str(p)}
    try:
        if not p.exists():
            return {**out, "exists": False}
        if p.is_dir():
            return {**out, "folder": True}
        size = p.stat().st_size
        h = hashlib.sha256()
        with open(p, "rb") as f:
            if size <= FULL_HASH_MAX:
                for chunk in iter(lambda: f.read(1 << 20), b""):
                    h.update(chunk)
                how = "sha256"
            else:
                h.update(f.read(PART))
                f.seek(size - PART)
                h.update(f.read(PART))
                how = "sha256 of the first and last 8 MiB"
        return {**out, "bytes": size, "modified": time.strftime("%Y-%m-%d %H:%M", time.localtime(p.stat().st_mtime)),
                "checksum": h.hexdigest(), "checksum_kind": how}
    except OSError as e:
        return {**out, "error": str(e)}


def _paths_in(x, out):
    """Absolute file paths anywhere in a (nested) metadata value."""
    if isinstance(x, dict):
        for v in x.values():
            _paths_in(v, out)
    elif isinstance(x, (list, tuple)):
        for v in x:
            _paths_in(v, out)
    elif isinstance(x, str) and len(x) < 1000 and os.path.isabs(x) and x not in out:
        out.append(x)
    return out


def _versions():
    from importlib import metadata
    out = {"python": platform.python_version()}
    for p in PACKAGES:
        try:
            out[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            pass
    return out


def _hex(rgb):
    return "#%02x%02x%02x" % tuple(int(c) for c in rgb[:3])


def _jdefault(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, Path):
        return str(o)
    return str(o)


def _dump(x):
    return json.dumps(x, indent=1, default=_jdefault, ensure_ascii=False)


def _short(x, n=1500):
    t = json.dumps(x, default=_jdefault, ensure_ascii=False)
    return t if len(t) <= n else t[:n] + " …"


def _table_md(df, max_rows=60):
    df = df.head(max_rows)
    cols = [str(c) for c in df.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        lines.append("| " + " | ".join("" if pd.isna(v) else str(v) for v in r.to_numpy()) + " |")
    return "\n".join(lines)


# ---------------- the guide ----------------
def build(s, version=None, with_notebook=True, progress=None):
    """Write the guide of a version (default: latest) as a zip into the session's exports; return its path."""
    say = progress or (lambda f, m: None)
    from . import notebook
    st = s.state(version)
    v = st["version"]
    work = s.dir / "exports" / f"guide_v{v:04d}"
    if work.exists():
        shutil.rmtree(work)
    data = work / "data"
    (data / "results").mkdir(parents=True)
    M = {"format": 1, "app": {"name": "STAnnotator", "version": __version__, "build": config.build(),
                              "source": f"{REPO}/tree/{_ref()}"},
         "exported": time.strftime("%Y-%m-%d %H:%M"), "sample_id": s.sample_id, "annotation_id": s.id, "version": v,
         "software": _versions(), "machine": f"{platform.system()} {platform.machine()}, {os.cpu_count()} cores",
         "files": {}, "steps": {}, "checks": [], "skipped": []}
    md = []

    def skip(what, e):
        M["skipped"].append({"part": what, "reason": f"{type(e).__name__}: {e}"})

    # ---- data files of the annotation
    say(0.05, "the cell table")
    exp = Path(s.export(v))
    for f in ("cells.csv", "state.json", "conflicts.csv", "types_not_found.csv", "CHANGELOG.md"):
        if (exp / f).exists():
            shutil.copy(exp / f, data / f)
    cells = pd.read_csv(data / "cells.csv", low_memory=False)
    ids = s.cell_ids()
    sch = s.scheme
    if sch.get("text"):
        (data / "marker_list.txt").write_text(sch["text"] + "\n", encoding="utf-8")
        try:
            from . import markers
            (data / "marker_list.csv").write_text(markers.text_to_csv(sch["text"]), encoding="utf-8")
        except Exception as e:
            skip("marker_list.csv", e)
    try:
        (data / "registration.csv").write_text(notebook.registration_csv(s.sample_id), encoding="utf-8")
    except Exception as e:
        skip("registration.csv", e)
    say(0.15, "cell positions")
    try:
        (data / "cell_positions.csv").write_text(notebook.positions_csv(s.sample_id), encoding="utf-8")
    except Exception as e:
        skip("cell_positions.csv", e)
    try:
        from . import webapp
        (data / "regions.geojson").write_text(webapp._regions_geojson(s, v)[0], encoding="utf-8")
    except Exception as e:
        skip("regions.geojson", e)
    corr_rows = []
    for c in st.get("corrections", []):
        for i in c.get("cells", []):
            corr_rows.append({"cell_id": ids[int(i)], "label": c.get("label"), "correction": c.get("id"),
                              "note": c.get("text") or c.get("reason") or ""})
    if corr_rows:
        pd.DataFrame(corr_rows).to_csv(data / "corrections.csv", index=False)
    colors = _colors(s, v, skip)
    (data / "colors.json").write_text(_dump(colors), encoding="utf-8")

    # ---- kept analysis results
    say(0.3, "the analyses' results")
    results = []
    rd = s.dir / "results"
    for f in sorted(rd.glob("*.json")) if rd.exists() else []:
        try:
            rec = json.load(open(f))
        except ValueError as e:
            skip(f"results/{f.name}", e)
            continue
        kind = f.stem
        shutil.copy(f, data / "results" / f.name)
        r = rec.get("result") or {}
        tables = []
        for tid, suffix in ((r.get("id"), ""), ((r.get("test") or {}).get("id"), "_test")):
            if tid and (rd / f"{tid}.parquet").exists():
                try:
                    pd.read_parquet(rd / f"{tid}.parquet").to_csv(data / "results" / f"{kind}{suffix}.csv", index=False)
                    tables.append(f"data/results/{kind}{suffix}.csv")
                except Exception as e:
                    skip(f"results/{kind}{suffix}.csv", e)
        results.append({"kind": kind, "time": rec.get("time"), "request": rec.get("request") or {},
                        "selection_cells": sum(len(x) for x in (rec.get("selection") or {}).values()) if isinstance(rec.get("selection"), dict) else None,
                        "result": r, "tables": tables})
    if (s.dir / "analyses.json").exists():
        shutil.copy(s.dir / "analyses.json", data / "analyses_log.json")

    # ---- trajectories: every cell's place along the last run of each method
    try:
        from . import trajectory
        for m in trajectory.all_last(s):
            meta, t, _ = trajectory.last(s, m)
            if t is not None:
                tt = np.asarray(t, float)
                pd.DataFrame({"cell_id": ids, "t": np.round(tt, 5), "bin": trajectory.bin_of(tt)}).to_csv(
                    data / f"trajectory_{m}.csv", index=False)
    except Exception as e:
        skip("trajectory_<method>.csv", e)

    # ---- the guide
    say(0.45, "writing the guide")
    _write_guide(s, st, v, cells, results, colors, M, md, skip)
    (work / "REPRODUCE.md").write_text("\n".join(md).rstrip() + "\n", encoding="utf-8")
    if with_notebook:
        say(0.55, "the notebook (figures included)")
        try:
            shutil.copy(notebook.build(s, v), work / "analysis.ipynb")
            M["files"]["analysis.ipynb"] = "the notebook: the exact replay of the data, annotation, embedding and niches"
        except Exception as e:
            skip("analysis.ipynb", e)
    for f in sorted(data.rglob("*")):
        if f.is_file():
            M["files"].setdefault(str(f.relative_to(work)), f.stat().st_size)
    (work / "manifest.json").write_text(_dump(M), encoding="utf-8")
    say(0.9, "packing the zip")
    from . import projects
    out = s.dir / "exports" / f"{projects.file_stem(s.sample_id)}_v{v}_reproduce.zip"
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for f in sorted(work.rglob("*")):
            if f.is_file():
                z.write(f, f"{out.stem}/{f.relative_to(work)}")
    shutil.rmtree(work, ignore_errors=True)
    say(1.0, "done")
    return out


def _colors(s, v, skip):
    """The map's colors: cell types, lineages, niches, Leiden clusters, regions and the gene / program color scale."""
    out = {}
    try:
        from . import webapp, embedding, niches
        out["cell_type"] = {t: _hex(c) for t, c in webapp.type_palette(s).items()}
        out["lineage"] = {l: _hex(c) for l, c in webapp.lineage_rgb(s).items()}
        names = niches.niche_names(s)
        if names:
            out["niche"] = {n: _hex(c) for n, c in zip(names, webapp.niche_rgb(s))}
        lab = embedding.leiden(s.sample_id)
        if lab is not None:
            out["leiden"] = {f"L{i + 1}": _hex(webapp.cluster_rgb(i)) for i in range(int(lab.max()) + 1)}
        out["region"] = {n: _hex(c) for n, c in webapp.region_rgb(s, v).items()}
        out["gene_scale"] = {"value": webapp.GENE_UNIT, "ramp_low_to_high": [_hex(c) for c in webapp.GENE_RAMP],
                             "not_detected_rgba": list(webapp.GENE_NONE),
                             "top": "the 99th percentile of the value among the cells where the gene is detected; in a "
                                    "project the largest of its slices' tops, so all slices share one scale",
                             "alpha": "215 + 40 * min(value / top, 1) of 255 for detected cells",
                             "order": "cells where the gene is not detected are drawn first (under the others)"}
        out["program_scale"] = "as gene_scale, on the program score above background (scores at or below 0: not detected)"
        out["unassigned"] = _hex(webapp.GRAY)
    except Exception as e:
        skip("colors.json", e)
    return out


def _write_guide(s, st, v, cells, results, colors, M, md, skip):
    meta, sch = s.meta, s.scheme
    n = len(cells)
    sec = iter(range(1, 40))
    add = md.append
    add(f"# Reproduction guide: sample `{s.sample_id}`, annotation `{s.id}`, version {v}")
    add("")
    add(f"Exported {M['exported']} by STAnnotator {__version__} (build `{M['app']['build']}`). This file is written for an "
        "AI agent (or a person) who has to repeat this analysis and redraw its figures **without the STAnnotator "
        "interface**. Read it to the end before starting.")
    add("")
    add("## How to use this guide")
    add("""
- **Goal**: starting from the input files listed below, produce the same cells, cell types, clusters, niches and
  analysis results, and figures drawn the same way. Every step gives its settings, the code that does it (a link to
  the public source at this version) and the results to check against (`manifest.json` → `checks`).
- **Route A, exact** (recommended): run `analysis.ipynb` from the top, if it is in this folder. It carries the app's
  own algorithm code and replays the data, placement, annotation, PCA / UMAP / Leiden and niches, and checks the
  labels against the originals. Or install the same version of the software and call the functions named below:
  `pip install "git+""" + REPO + f"""@{_ref()}"` (the package is `stagent`, version {__version__}; the source is
  public under the PolyForm Noncommercial license).
- **Route B, re-implementation**: follow the steps and the algorithm descriptions (the appendix quotes the source's own
  descriptions); expect small differences where a step is stochastic or depends on library versions (see
  *Determinism*). Compare with the checks; report any step that differs and by how much.
- **Files here**: `data/cells.csv` holds every cell's final labels (cell type, lineage, levels, confidence, region,
  niche, Leiden cluster, UMAP, coordinates in µm); the other files of `data/` are named in the steps. Paths of the
  original input files are those of the computer this was exported from; ask the user where they are now and check
  them by size and checksum.
- Cell ids: `cell_id` in `data/cells.csv` is the id of the original data. Rows of `data/cells.csv` are in the
  app's cell order (the order after the cell filter), which the notebook and the state files use.
""".strip("\n"))
    add("")

    # ---- software
    add(f"## {next(sec)}. Software")
    vers = M["software"]
    add(f"STAnnotator {__version__}, build `{M['app']['build']}`, source: {M['app']['source']} (the links below "
        "point into this commit, so their line numbers match). Exported on "
        f"{M['machine']}. Packages of the app at export (those that change results are marked):")
    add("")
    pinned = {"pynndescent", "numba", "igraph", "leidenalg", "numpy", "scipy", "umap-learn", "scikit-learn"}
    for k, x in vers.items():
        add(f"- {k} {x}" + (" (can change results: use this version)" if k in pinned else
                            " (only speed: exact neighbors on the GPU give the same result)" if k in ("torch", "mlx") else ""))
    add("")

    # ---- input data
    add(f"## {next(sec)}. Input data")
    src = meta.get("source") or {}
    add(f"Sample `{s.sample_id}`" + (f" ({meta.get('description')})" if meta.get("description") else "") +
        f": imported from **{src.get('kind') or 'files'}** on {meta.get('created_at') or '?'}; "
        f"{meta.get('n_cells', n):,} cells and {meta.get('n_genes', '?')} genes in the app.")
    paths = _paths_in({"source": src, "he_path": meta.get("he_path"), "images": meta.get("images")}, [])
    infos = []
    for p in paths[:40]:
        infos.append(_file_info(p))
    M["inputs"] = {"source": src, "image": {"he_path": meta.get("he_path"), "um_per_px": meta.get("he_um_per_px"),
                                             "kind": meta.get("image_kind"), "images": meta.get("images"),
                                             "hidden": meta.get("images_hidden")},
                   "files": infos, "organ": meta.get("organ"), "project": meta.get("project"),
                   "project_section": meta.get("project_section"), "metadata": meta.get("metadata")}
    add("")
    add("| file | bytes | checksum |")
    add("|---|---|---|")
    for x in infos:
        add(f"| `{x['path']}` | {x.get('bytes', 'folder' if x.get('folder') else 'missing' if x.get('exists') is False else '')} | "
            + (f"{x['checksum_kind']} `{x['checksum']}`" if x.get("checksum") else "") + " |")
    add("")
    add("Reading the data: the importer of its kind, " + ("Xenium output (cells and cell_feature_matrix.h5)" if src.get("kind") == "xenium"
        else "an h5ad / R file (raw counts, coordinates in µm)" if src.get("kind") in (None, "h5ad", "rds") else f"{src.get('kind')} output") + ":")
    add(_code_lines(["importer.import_sample", "importer.import_xenium", "importer.read_h5ad", "importer.flat_cells_counts",
                     "importer.import_visium", "importer.import_stereoseq"]))
    add("")

    # ---- filter
    flt = meta.get("filter") or {}
    add(f"## {next(sec)}. Cells and genes kept (filter)")
    if flt:
        add(f"Thresholds `{_short(flt.get('params'))}`: {flt.get('n_cells_imported', '?')} cells imported, "
            f"{flt.get('n_cells_removed', 0)} removed, {flt.get('n_genes_imported', '?')} genes imported, "
            f"{len(flt.get('genes_removed') or [])} genes removed" + (f" ({', '.join(flt['genes_removed'][:40])})" if flt.get("genes_removed") else "") +
            (f". Filtered together with the slices {', '.join(flt['slices'])} (one gene panel)" if len(flt.get("slices") or []) > 1 else "") + ".")
    else:
        add("No filter was applied: every imported cell and gene is used.")
    add("")
    add(_code_lines(["qcfilter.apply", "qcfilter.params"]))
    add("")
    M["steps"]["filter"] = flt
    M["checks"].append({"step": "filter", "what": "cells after the filter", "expected": n})
    add(f"**Check**: {n:,} cells remain (the rows of `data/cells.csv`, same ids).")
    add("")

    # ---- placement
    add(f"## {next(sec)}. Placement of the cells on the image")
    reg = {k: meta.get(k) for k in ("placement_affine", "registration_affine", "registration_by", "he_um_per_px", "image_kind")}
    M["steps"]["placement"] = reg
    if meta.get("he_path"):
        add(f"Image `{meta.get('he_path')}` ({meta.get('image_kind') or 'H&E'}, {meta.get('he_um_per_px')} µm per pixel). "
            f"Registration by: {meta.get('registration_by') or 'none'}. `data/registration.csv` holds the affine from the cells' "
            "µm coordinates to the image's full-resolution pixels; `data/cell_positions.csv` every cell's final place on "
            "the image (local corrections and moved pieces included). To overlay cells on the image use "
            "`data/cell_positions.csv`; the affine alone misses the local corrections when there are any.")
    else:
        add("No image: the cells are drawn at their µm coordinates.")
    add("")
    add(_code_lines(["importer.load_registration", "importer.positions_affine", "importer.apply_registration_fix",
                     "importer.replay_placement", "align.search", "align.refine"]))
    add("")

    # ---- annotation
    add(f"## {next(sec)}. Annotation")
    kind = sch.get("kind")
    how = {"markers": "the marker list (hierarchical marker clustering, then an expression model trained on those labels)",
           "scref": "a single-cell reference (a classifier on the genes shared with the panel)",
           "existing": "an existing annotation loaded from a CSV", "none": "no annotation"}.get(kind, kind)
    ann = {k: sch.get(k) for k in ("kind", "resolution", "knn", "neighbors", "spillover", "untyped_bf", "leaf_types",
                                   "reference", "existing")}
    ann["knn"] = s.knn_mode()
    M["steps"]["annotation"] = {**ann, "types": list(s.types), "regions": {k: x for k, x in st["regions"].items() if k != "ops"},
                                "n_region_ops": len(st["regions"].get("ops", [])),
                                "constraints": st.get("constraints", []),
                                "marker_ops": st.get("marker_step", {}).get("ops", []),
                                "n_corrections": len(st.get("corrections", []))}
    add(f"Method: {how}. Settings: `{_short(ann)}`.")
    add("")
    if kind == "markers":
        add("The marker list is `data/marker_list.txt` (indented text; the CSV form is `data/marker_list.csv`). "
            "Types with subtypes are groups; the cell types are the leaves. A gene with `!` is decisive; `-GENE` is a "
            "negative marker. The resolution `auto` means the app chose each level's Leiden resolution (the choice is "
            "in `data/state.json` and the notebook); a number is a fixed resolution.")
        add("")
        add("```")
        add((sch.get("text") or "").rstrip())
        add("```")
        add("")
        try:
            lv = s.marker_levels(v)
            rows = [{"level": x.get("node"), "cells": x.get("n_cells"), "resolution": x.get("resolution"),
                     "share assigned": x.get("assigned_frac")} for x in lv]
            if rows:
                M["steps"]["annotation"]["levels"] = rows
                add("Leiden resolution used at every level of the marker step (with `auto` the app tried several and "
                    "kept this one; to replay it, use these values):")
                add("")
                add(_table_md(pd.DataFrame(rows), 200))
                add("")
        except Exception as e:
            skip("marker levels", e)
    pri = st.get("constraints", [])
    try:
        from . import constraints as cmod
        add(f"Priors ({len(pri)}): " + ("; ".join(f"`{cmod.describe(c)}` (\"{c.get('text', '')}\")" for c in pri) or "none") + ".")
    except Exception:
        add(f"Priors ({len(pri)}): see `data/state.json` → constraints.")
    add(f"Regions: {len(st['regions'].get('ops', []))} drawn outlines (`data/regions.geojson`, image pixels; the cells' "
        "region is the `region` column of `data/cells.csv`).")
    add(f"Hand corrections: {len(st.get('corrections', []))} (cells and labels in `data/corrections.csv`; a later "
        "correction of the same cell wins; corrected cells keep their label).")
    ops = st.get("marker_step", {}).get("ops", [])
    add(f"Re-clusterings: {len(ops)}" + (": " + "; ".join(f"{o.get('id')} {o.get('target')} at {o.get('resolution')}" for o in ops) if ops else "") +
        " (`data/state.json` → marker_step).")
    add("History of the versions: `data/CHANGELOG.md`.")
    add("")
    add("Code:")
    add(_code_lines(["markers.parse", "markers.annotate", "markers.marker_assign", "markers.spill_detected",
                     "markers.training_labels", "model.lognorm", "model.train", "model.predict", "markers.restrict",
                     "markers.neighbor_rounds", "constraints.posterior", "regions.assign", "session.Session.result"]
                    if kind == "markers" else ["scref.inspect", "model.train", "model.predict", "constraints.posterior",
                                               "regions.assign"]))
    add("")
    counts = cells["cell_type"].value_counts()
    M["checks"].append({"step": "annotation", "what": "cells per cell type", "expected": counts.to_dict(),
                        "tolerance": "identical with route A; route B: each type within a few percent"})
    add("**Check**: cells per cell type (`data/cells.csv` → `cell_type`):")
    add("")
    add(_table_md(counts.rename_axis("cell type").reset_index(name="cells"), 200))
    add("")
    try:
        nf = pd.read_csv(Path(s.dir / "exports" / f"v{v:04d}" / "types_not_found.csv"))
        if len(nf):
            add("Types of the list without cells (`data/types_not_found.csv`): " +
                "; ".join(f"{a}: {b}" for a, b in zip(nf.cell_type, nf.reason)) + ".")
            add("")
    except Exception:
        pass
    fp = s.dir / "base" / "fingerprints.json"
    if fp.exists():
        try:
            M["steps"]["annotation"]["fingerprints"] = json.load(open(fp))
            add("Checksums of the marker step's intermediate results (PCA, neighbors, root clusters, labels) are in "
                "`manifest.json` → steps.annotation.fingerprints: where a replay first differs shows which step did.")
            add("")
        except ValueError:
            pass

    # ---- embedding
    from . import embedding, niches
    em = embedding.current(s.sample_id)
    add(f"## {next(sec)}. PCA, UMAP and Leiden clusters")
    if em:
        M["steps"]["embedding"] = em
        add(f"Settings `{_short(em)}`. On the log-normalized counts (ln(1 + counts per 100)) of all genes: PCA, "
            "nearest neighbors, Leiden, UMAP. Independent of the annotation. Columns `leiden`, `umap_1`, `umap_2` of "
            "`data/cells.csv`" + (" (computed on all slices of the project together)" if em.get("project") else "") + ".")
        add("")
        add(_code_lines(["embedding.run", "markers.pca", "markers.knn", "markers.knn_graph"]))
        add("")
        if "leiden" in cells:
            lc = cells["leiden"].value_counts()
            M["checks"].append({"step": "embedding", "what": "Leiden clusters and their sizes", "expected": lc.to_dict(),
                                "tolerance": "identical with route A and the same library versions; UMAP coordinates differ slightly"})
            add(f"**Check**: {len(lc)} Leiden clusters; sizes: " + ", ".join(f"{k} {x}" for k, x in lc.items()) + ".")
            add("")
    else:
        add("Not computed.")
        add("")

    # ---- niches
    add(f"## {next(sec)}. Niches")
    nmeta, nlab = niches.current(s)
    if nmeta is not None:
        M["steps"]["niches"] = nmeta
        add(f"Settings `{_short({k: nmeta.get(k) for k in ('k', 'k_by_region', 'radius_um', 'identity', 'n_programs', 'n_niches', 'version')})}`. "
            "Each cell's identity (" + ("its composition of NMF gene programs" if nmeta.get("identity") == "nmf" else "its cell type") +
            ") averaged over its k nearest neighbors in space (the cell included), then k-means (n_init 4, "
            "random_state 0) into the niches, numbered by size (N1 the largest). Names given by the user: "
            f"{_short(nmeta.get('names') or {})}. Column `niche` of `data/cells.csv`.")
        add("")
        add(_code_lines(["niches.programs", "niches.run"]))
        add("")
        if "niche" in cells:
            nc = cells["niche"].value_counts()
            M["checks"].append({"step": "niches", "what": "cells per niche", "expected": nc.to_dict()})
            add("**Check**: cells per niche: " + ", ".join(f"{k} {x}" for k, x in nc.items()) + ".")
            add("")
    else:
        add("Not computed.")
        add("")

    # ---- groupings
    add(f"## {next(sec)}. Groupings (slices and regions)")
    try:
        from . import projects
        pid = projects.project_of(s.sample_id)
        g = projects.groupings(pid) if pid else None
        rg = projects.region_groups(s.sample_id)
        M["steps"]["groupings"] = {"project": pid, "slices": g, "region_groups": rg}
        if g:
            add(f"Project `{g['name']}` ({pid}), slices: " + ", ".join(f"`{x['sample_id']}`" for x in g["slices"]) + ".")
            for k, col in g["columns"].items():
                add(f"- slice grouping **{k}**: " + ", ".join(f"{a} = {b}" for a, b in col.items()))
        else:
            add("Not part of a project.")
        for k, col in (rg or {}).items():
            add(f"- region grouping **{k}**: " + ", ".join(f"{a} → {b}" for a, b in col.items()) + " (other regions: none)")
        add("")
        add("A grouping is a column of the cell table used by the analyses below (`by`, `within`, `rows`, `cols`, "
            "`split`, `replicate`): a slice grouping gives every cell of a slice the slice's value; a region grouping "
            "gives a cell the group of its region.")
    except Exception as e:
        skip("groupings", e)
        add("See `manifest.json`.")
    add("")

    # ---- analyses
    add(f"## {next(sec)}. Analyses kept with this annotation")
    add("The last result of each kind of analysis the user ran (in the app or through its agent), with the request "
        "that made it. `data/results/<kind>.json` holds the request and the full result; `data/results/<kind>.csv` "
        "its full table. `data/analyses_log.json` lists the last 30 runs (any kind, with headlines).")
    add("")
    add("The cells of an analysis: `scope` slice = this annotation's cells; project = the cells of every slice named "
        "in `sessions` (sample → annotation id), with their niches / Leiden runs in `runs`. `within` filters the cells "
        "({grouping: values}); a `selection` (cells drawn on the map) is counted in `selection_cells` and its cells "
        "are in the JSON. Groupings available: cell type, lineage, region, niche, leiden, selection, slice and the "
        "groupings above. Settings missing from a request took their defaults; the result JSON records the values "
        "used (e.g. distances, minimum shares, permutations).")
    add("")
    M["steps"]["analyses"] = []
    from . import analysis
    for r in results:
        kind = r["kind"]
        title = RESULT_TITLES.get(kind, kind)
        req = {k: x for k, x in r["request"].items() if k not in ("selection",)}
        add(f"### {title} (`{kind}`), run {r.get('time') or '?'}")
        add(f"Request: `{_short(req, 2500)}`" + (f"; selection of {r['selection_cells']} cells" if r.get("selection_cells") else "") + ".")
        res = r["result"]
        if res.get("note"):
            add(f"Note of the result: {res['note']}")
        hk = "trajectory" if kind.startswith("trajectory") else kind
        try:
            head = analysis.headline(hk, res)
        except Exception:
            head = {}
        if head:
            add(f"**Check** (headline): `{_short(head, 2500)}`")
            M["checks"].append({"step": f"analysis {kind}", "what": "headline of the result", "expected": head,
                                "tolerance": "the same top genes / groups; p-values may differ slightly with library versions"})
        add("Files: `data/results/" + kind + ".json`" + "".join(f", `{t}`" for t in r["tables"]) + ".")
        add("Code:")
        add(_code_lines(RESULT_CODE.get(kind, [])))
        add("")
        M["steps"]["analyses"].append({"kind": kind, "title": title, "time": r.get("time"), "request": req,
                                       "files": [f"data/results/{kind}.json", *r["tables"]],
                                       "code": [_src(x) for x in RESULT_CODE.get(kind, [])]})
    if not results:
        add("None kept.")
        add("")

    # ---- figures
    add(f"## {next(sec)}. Figures: how the app draws them")
    add("""
- **Coordinates**: cells at their µm coordinates (`x_um`, `y_um`), y pointing down as in the image (invert the y
  axis). On the image: `data/cell_positions.csv` (image pixels). Equal aspect ratio; a scale bar in µm.
- **Colors** (`data/colors.json`): cell types (`cell_type`), lineages, niches, Leiden clusters, regions, each as hex.
  Unassigned / Low quality / not in a region: gray. Cells are drawn as filled dots, alpha about 0.9, small (about
  the cell size at full view; 1–4 px in a whole-slide figure).
- **Gene expression**: the value is ln(1 + transcripts per 100) of the cell (`model.lognorm`); color from the ramp
  in `colors.json` → gene_scale, between 0 and the top (99th percentile of detected cells; in a project one shared
  top: the largest of the slices' tops); cells where the gene is not detected are drawn first in faint gray, the
  others on top in order of value; a continuous color bar with the unit.
- **Program scores**: like a gene, on the score above background (`programs.score`).
- **UMAP**: `umap_1`, `umap_2` colored like the map.
- **DE genes**: volcano-like or dot plots from `data/results/de.csv` (log2 fold change, adjusted p, shares of cells
  expressing); top genes per group as dot plots (mean and share per group) from `data/results/markers.csv`.
- **Composition**: stacked bars of shares (`data/results/composition*.csv`).
- **Communication**: interactions table and the network of groups (`data/results/run.csv`).
- **Trajectories**: cells colored by `t` (0 to 1) of `data/trajectory_<method>.csv` in 10 bins; gene trends along t.
- The notebook draws the map, UMAP and niche figures with these rules (section "Packages", function `plot`).
""".strip("\n"))
    add("")
    add(_code_lines(["webapp.type_palette", "webapp.lineage_rgb", "webapp.gene_rgba", "webapp.slice_tops",
                     "webapp.project_top", "render.region_colors"]))
    add("")

    # ---- determinism
    add(f"## {next(sec)}. Determinism")
    knn = s.knn_mode()
    add(f"- Neighbor search of the marker step: **{knn}**. "
        + ("Exact neighbors are the same on every computer, so the labels replay exactly. " if knn == "exact" else
           "Approximate neighbors (pynndescent) depend on its version, numba and the thread count (the app uses 8 "
           "threads: set NUMBA_NUM_THREADS to at least 8 before numba is imported). "))
    add("- Leiden runs with seed 0 (igraph's Leiden, the app's igraph version). UMAP runs on all cores, so its "
        "coordinates differ slightly on a rerun (the clusters do not depend on UMAP).")
    add("- On a Mac with MLX the app turns off TF32 (MLX_ENABLE_TF32=0) so exact neighbors are computed in full float32.")
    add("- Niches: k-means with random_state 0 on float64; NMF with init nndsvda, random_state 0.")
    add("")

    # ---- appendix: the source's own descriptions
    add("## Appendix: the algorithms, as the source describes them")
    mods = ["qcfilter", "markers", "model", "constraints", "regions", "embedding", "niches", "analysis"]
    for r in results:
        m = RESULT_MODULE.get(r["kind"])
        if m and m not in mods:
            mods.append(m)
    for m in mods:
        d = _module_doc(m)
        if d:
            add(f"### `stagent/{m}.py` ({REPO}/blob/{_ref()}/stagent/{m}.py)")
            add("```text")
            add(d)
            add("```")
            add("")
    if M["skipped"]:
        add("## Parts that could not be written")
        for x in M["skipped"]:
            add(f"- {x['part']}: {x['reason']}")
        add("")
