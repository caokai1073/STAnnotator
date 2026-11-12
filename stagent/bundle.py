"""An annotation as one file (.stannot, a zip): save it, and load it again here or in another data folder.

The file holds the annotation's folder (every version with its regions, priors and corrections, the marker step and
expression model, niches with their names, review verdicts and advice, area suggestions, the report, the last result
of each analysis: DE genes, top genes, composition, cell-cell communication) and the sample's PCA / UMAP / Leiden
clusters. Chats are left out (they are conversations with Claude, not the annotation), and so are exports and
figures (made again on demand).

It also holds the sample's state the annotation was made on (sample/): the cell and gene filter (thresholds and the
genes left out) and the registration (the placement, its corrections and piece moves, every cell's place on the
image, the image's size), so that the same data imported again needs nothing else.

Loading needs the same data: the sample with the same cell ids and genes in the same order. A sample imported again
from the same files (before any filter) also matches: plan() says what loading would change there (filter it as
saved, which renumbers its cells, so its own annotations and PCA / UMAP go to the trash first; put its cells where
they were on the same image), the caller asks the user, apply_sample() does it, then load(). The annotation comes
back as a new annotation of that sample; a different PCA / UMAP of the sample is moved to the trash before the saved
one replaces it.
"""
import hashlib, io, json, shutil, time, uuid, zipfile
from pathlib import PurePosixPath
from . import config
from .osutil import dump_json

FORMAT = 3  # 2: the sample's filter and registration (sample/); 3: with the cell positions file of a registration
#   by cell positions (sample/registration_positions.parquet); a file is marked 3 only when it holds that one, so the
#   others still load in an app of format 2; an older app asks to be updated
SUFFIX = ".stannot"
LEAVE_OUT = {"chat.jsonl", "chat_meta.json"}
LEAVE_OUT_DIRS = {"chat_archive", "exports", "figures"}
TYPES = {".json", ".jsonl", ".npy", ".npz", ".parquet", ".md", ".txt", ".csv"}  # what a bundle may hold


def cells_key(sample_id):
    """A fingerprint of the sample's cells (ids in order) and of its genes (in order)."""
    from .session import sample_cell_ids, panel_genes
    ids, genes = sample_cell_ids(sample_id), panel_genes(sample_id)
    return {"n_cells": int(len(ids)), "sha1": hashlib.sha1("\n".join(ids).encode()).hexdigest(),
            "n_genes": len(genes), "genes_sha1": hashlib.sha1("\n".join(genes).encode()).hexdigest()}


REG_KEYS = ("placement_affine", "placement_affine_import", "registration_affine", "registration_by", "registration_fixes")


def _imported_dir(sample_id):
    d = config.SAMPLES / sample_id
    return d / "original" if (d / "original" / "cells.parquet").exists() else d


def imported_key(sample_id):
    """The fingerprint (as cells_key) of the sample's cells and genes as imported, before any filter."""
    import pandas as pd
    src = _imported_dir(sample_id)
    c = pd.read_parquet(src / "cells.parquet")
    ids = c.cell_id.astype(str).to_numpy() if "cell_id" in c else c.index.astype(str).to_numpy()
    genes = config.read_genes(src / "genes.txt")
    return {"n_cells": int(len(ids)), "sha1": hashlib.sha1("\n".join(ids).encode()).hexdigest(),
            "n_genes": len(genes), "genes_sha1": hashlib.sha1("\n".join(genes).encode()).hexdigest()}


def _meta(sample_id):
    return json.load(open(config.SAMPLES / sample_id / "meta.json"))


def _image(meta):
    """The sample's image: the one its cells are on, else the one waiting for a registration."""
    return meta.get("he_path") or meta.get("image_pending")


def _sample_state(sample_id):
    from .importer import _he_size
    m = _meta(sample_id)
    f = m.get("filter")
    st = {"imported": imported_key(sample_id),
          "filter": {"params": f.get("params"), "genes_removed": f.get("genes_removed") or []} if f else None,
          "registration": None}
    if m.get("he_path") and not m.get("image_only") and m.get("placement_affine"):
        try:
            size = list(_he_size(m["he_path"]))
        except Exception:
            size = None
        st["registration"] = {**{k: m.get(k) for k in REG_KEYS}, "image_kind": m.get("image_kind"),
                              "image_name": PurePosixPath(m["he_path"].replace("\\", "/")).name, "image_size": size}
        from .importer import POSITIONS_FILE
        if m.get("registration_positions") and (config.SAMPLES / sample_id / POSITIONS_FILE).exists():
            st["registration"]["registration_positions"] = m["registration_positions"]
    return st


def _files(d, prefix):
    for p in sorted(d.rglob("*")):
        rel = p.relative_to(d)
        if p.is_file() and p.name not in LEAVE_OUT and not (set(rel.parts[:-1]) & LEAVE_OUT_DIRS) \
                and p.suffix in TYPES and not p.name.startswith("."):
            yield p, f"{prefix}/{rel.as_posix()}"


def save(session_id, out):
    """Write the annotation session_id to out (a path or a binary file object); returns its description."""
    from .session import Session
    s = Session(session_id)
    emb = config.SAMPLES / s.sample_id / "embedding"
    info = {"format": FORMAT, "app": "STAnnotator", "build": config.build(), "saved_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "session_id": s.id, "sample_id": s.sample_id, "sample_name": s.meta.get("name"),
            "cells": cells_key(s.sample_id), "filter": (s.meta.get("filter") or {}).get("params"),
            "method": s.scheme["kind"], "n_versions": len(s.versions()),
            "has_niches": (s.dir / "niches" / "run.json").exists(), "has_embedding": (emb / "umap.npy").exists(),
            "results": sorted(p.stem for p in (s.dir / "results").glob("*.json")) if (s.dir / "results").exists() else [],
            "sample_state": _sample_state(s.sample_id)}
    if not (info["sample_state"]["registration"] or {}).get("registration_positions"):
        info["format"] = 2  # (nothing of format 3 in it: an app of format 2 loads it)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=4) as z:
        z.writestr("bundle.json", json.dumps(info, indent=1))
        if info["sample_state"]["registration"]:  # every cell's place on the image (corrections are no affine)
            import pandas as pd
            c = pd.read_parquet(config.SAMPLES / s.sample_id / "cells.parquet")
            cols = [x for x in ("cell_id", "px", "py", "px_import", "py_import") if x in c]
            buf = io.BytesIO(); c[cols].to_parquet(buf); z.writestr("sample/positions.parquet", buf.getvalue())
            pf = config.SAMPLES / s.sample_id / "pieces.json"
            if pf.exists():
                z.write(pf, "sample/pieces.json")
            if info["sample_state"]["registration"].get("registration_positions"):  # the whole file: a filter, Restore
                from .importer import POSITIONS_FILE  # all and the notebook replay it for cells left out above
                z.write(config.SAMPLES / s.sample_id / POSITIONS_FILE, "sample/registration_positions.parquet")
        for p, name in _files(s.dir, "annotation"):
            z.write(p, name)
        if info["has_embedding"]:
            for p, name in _files(emb, "embedding"):
                z.write(p, name)
    return info


def _check(z):
    if "bundle.json" not in z.namelist():
        raise ValueError("this is not a saved annotation (no bundle.json)")
    info = json.loads(z.read("bundle.json"))
    if info.get("format", 0) > FORMAT:
        raise ValueError("this annotation was saved by a newer STAnnotator: update the app to load it")
    for n in z.namelist():
        p = PurePosixPath(n)
        if n == "bundle.json" or n.endswith("/"):
            continue
        if p.parts[0] == "sample" and len(p.parts) == 2 and p.name in ("positions.parquet", "pieces.json",
                                                                         "registration_positions.parquet"):
            continue
        if p.is_absolute() or ".." in p.parts or p.parts[0] not in ("annotation", "embedding") or p.suffix not in TYPES:
            raise ValueError(f"unexpected file in the annotation file: {n}")
    if "annotation/session.json" not in z.namelist():
        raise ValueError("the annotation file has no session.json")
    return info


def inspect(src):
    with zipfile.ZipFile(src) as z:
        return _check(z)


def _same(a, b):
    return bool(a) and bool(b) and a.get("sha1") == b.get("sha1") and a.get("genes_sha1") == b.get("genes_sha1")


def plan(src, sample_id=None):
    """What loading src would do, before anything changes: the target sample, and whether its filter and its
    registration are made as saved. {"target", "filter": None | "apply", "registration": None | "apply" | "same" |
    "no image" | "other image", "n_annotations" (to the trash, when filtered), "embedding" (idem), "notes": [...]}."""
    from .session import list_samples
    from .importer import _he_size
    with zipfile.ZipFile(src) as z:
        info = _check(z)
        has_pos = "sample/positions.parquet" in z.namelist()
    want, st = info["cells"], info.get("sample_state") or {}
    tried, target, refilter = [], None, False
    order = [x for x in (sample_id, info["sample_id"]) if x] + [x for x in list_samples().get("sample_id", [])]
    for cand in dict.fromkeys(order):
        if not (config.SAMPLES / cand / "meta.json").exists():
            continue
        if _meta(cand).get("n_cells") == want["n_cells"]:
            k = cells_key(cand)
            if _same(k, want):
                target = cand
                break
            if k["sha1"] == want["sha1"]:
                tried.append(f"{cand} has the same cells but another gene panel ({k['n_genes']} genes, saved with {want['n_genes']})")
        if st.get("imported") and _same(imported_key(cand), st["imported"]):  # the same data, filtered otherwise
            target, refilter = cand, True
            break
    if target is None:
        raise ValueError(f"no sample here has the cells of this annotation: it was saved on sample {info['sample_id']} "
                         f"with {want['n_cells']:,} cells" + (f" (filtered: {info['filter']})" if info.get("filter") else "")
                         + ". Import that sample from the same files first."
                         + (f" ({'; '.join(tried)})" if tried else ""))
    out = {"target": target, "filter": "apply" if refilter else None, "registration": None, "notes": [],
           "saved_filter": st.get("filter"), "n_annotations": 0, "embedding": False}
    if refilter:
        from .session import list_sessions
        t = list_sessions()
        out["n_annotations"] = int(((t.sample_id == target) & (t.method != "none")).sum()) if len(t) else 0
        out["n_drafts"] = int(((t.sample_id == target) & (t.method == "none")).sum()) if len(t) else 0
        out["embedding"] = (config.SAMPLES / target / "embedding").exists()
    reg = st.get("registration")
    if reg and has_pos:
        m = _meta(target)
        img = _image(m)
        if not img:
            out["registration"] = "no image"
            out["notes"].append(f"the registration was not applied: {target} has no image; add the same image "
                                f"({reg.get('image_name')}) and load again")
        else:
            try:
                size = list(_he_size(img))
            except Exception:
                size = None
            if reg.get("image_size") and size and size != list(reg["image_size"]):
                out["registration"] = "other image"
                out["notes"].append(f"the registration was not applied: the saved one is on an image of {reg['image_size'][0]} x "
                                    f"{reg['image_size'][1]} px ({reg.get('image_name')}), {target}'s is {size[0]} x {size[1]} px")
            elif (not refilter and m.get("he_path") and not m.get("image_pending")
                  and m.get("placement_affine") == reg.get("placement_affine")
                  and (m.get("registration_fixes") or []) == (reg.get("registration_fixes") or [])):
                out["registration"] = "same"
            else:
                out["registration"] = "apply"
                out["registration_by"] = reg.get("registration_by")
                out["image"] = img
    return out


def apply_sample(src, p, progress=lambda *a: None):
    """Make the target of plan p as saved: filter (p["filter"] == "apply"; its annotations and PCA / UMAP must be in
    the trash already), then the registration (p["registration"] == "apply"): the cells placed on the image with the
    saved affine, then every cell put where it was (local corrections and piece moves are no affine)."""
    import numpy as np, pandas as pd
    from . import importer, qcfilter
    target = p["target"]
    with zipfile.ZipFile(src) as z:
        info = _check(z)
        st = info.get("sample_state") or {}
        if p.get("filter") == "apply":
            f = st.get("filter")
            qcfilter.apply([target], (f or {}).get("params") or {}, progress,
                           genes_removed=(f or {}).get("genes_removed") if f else None)
        if p.get("registration") == "apply":
            reg = st["registration"]
            base = reg.get("placement_affine_import") if reg.get("registration_fixes") else reg.get("placement_affine")
            by = (reg.get("registration_by") or "file").split(" + ")[0]
            importer.place_cells(target, base or reg["placement_affine"], p["image"], reg.get("image_kind") or "H&E", by=by,
                                 progress=lambda *a: None)
            m = _meta(target)
            m.update({k: reg.get(k) for k in REG_KEYS if reg.get(k) is not None})
            if reg.get("registration_fixes"):  # every cell where it was: local shifts and piece moves
                d = config.SAMPLES / target
                cells = pd.read_parquet(d / "cells.parquet")
                pos = pd.read_parquet(io.BytesIO(z.read("sample/positions.parquet"))).set_index("cell_id")
                ids = cells.cell_id.astype(str)
                pos.index = pos.index.astype(str)
                hit = ids.isin(pos.index).to_numpy()
                for col in ("px", "py", "px_import", "py_import"):
                    if col in pos:
                        v = np.array(cells[col], float) if col in cells else np.full(len(cells), np.nan)  # (a copy: pandas' may be read-only)
                        v[hit] = pos.loc[ids[hit], col].to_numpy(float)
                        cells[col] = v
                cells.to_parquet(d / "cells.parquet")
                if "sample/pieces.json" in z.namelist():
                    (d / "pieces.json").write_bytes(z.read("sample/pieces.json"))
            if reg.get("registration_positions") and "sample/registration_positions.parquet" in z.namelist():
                f = config.SAMPLES / target / importer.POSITIONS_FILE  # (the cells a filter left out placed by it later)
                f.parent.mkdir(parents=True, exist_ok=True)
                f.write_bytes(z.read("sample/registration_positions.parquet"))
                m["registration_positions"] = reg["registration_positions"]
            dump_json(m, config.SAMPLES / target / "meta.json", indent=1)


def load(src, sample_id=None):
    """Load a saved annotation (a path or a binary file object) as a new annotation of the sample with the same
    cells (sample_id first, then the saved sample, then any other). Returns what was loaded."""
    from .session import list_samples
    with zipfile.ZipFile(src) as z:
        info = _check(z)
        want = info["cells"]
        tried, target = [], None  # tried: samples with the same cells but another panel
        order = [x for x in (sample_id, info["sample_id"]) if x] + [x for x in list_samples().get("sample_id", [])]
        for cand in dict.fromkeys(order):
            if not (config.SAMPLES / cand / "meta.json").exists():
                continue
            if json.load(open(config.SAMPLES / cand / "meta.json")).get("n_cells") != want["n_cells"]:
                continue
            k = cells_key(cand)
            if k["sha1"] == want["sha1"]:
                if want.get("genes_sha1") and k["genes_sha1"] != want["genes_sha1"]:
                    tried.append(f"{cand} has the same cells but another gene panel ({k['n_genes']} genes, saved with {want['n_genes']})")
                    continue
                target = cand
                break
        if target is None:
            raise ValueError(f"no sample here has the cells of this annotation: it was saved on sample {info['sample_id']} "
                             f"with {want['n_cells']:,} cells" + (f" (filtered: {info['filter']})" if info.get("filter") else "")
                             + ". Import that sample from the same files first."
                             + (f" ({'; '.join(tried)})" if tried else ""))
        sid = f"{target}_{time.strftime('%Y%m%d-%H%M%S')}_{uuid.uuid4().hex[:4]}"
        d = config.SESSIONS / sid
        d.mkdir(parents=True)
        try:
            for n in z.namelist():
                if n.startswith("annotation/") and not n.endswith("/"):
                    f = d / n[len("annotation/"):]
                    f.parent.mkdir(parents=True, exist_ok=True)
                    f.write_bytes(z.read(n))
            sj = json.loads((d / "session.json").read_text(encoding="utf-8"))
            sj.update(session_id=sid, sample_id=target, loaded_from={"session_id": info["session_id"], "sample_id": info["sample_id"],
                                                                      "saved_at": info["saved_at"]})
            (d / "session.json").write_text(json.dumps(sj, indent=1), encoding="utf-8")
        except Exception:
            shutil.rmtree(d, ignore_errors=True)
            raise
        replaced = False
        emb_files = [n for n in z.namelist() if n.startswith("embedding/") and not n.endswith("/")]
        if emb_files:
            emb = config.SAMPLES / target / "embedding"
            saved_run = json.loads(z.read("embedding/meta.json")).get("run_id") if "embedding/meta.json" in emb_files else None
            cur_run = json.load(open(emb / "meta.json")).get("run_id") if (emb / "meta.json").exists() else None
            if not emb.exists() or cur_run != saved_run:
                if emb.exists():  # the sample's own PCA / UMAP, kept in the trash
                    trash = config.TRASH / f"{time.strftime('%Y%m%d-%H%M%S')}_{target}_embedding"
                    trash.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(emb), str(trash))
                    replaced = True
                emb.mkdir(parents=True)
                for n in emb_files:
                    (emb / n[len("embedding/"):]).write_bytes(z.read(n))
    return {"session_id": sid, "sample_id": target, "saved": info, "replaced_embedding": replaced,
            "embedding_loaded": bool(emb_files)}


def save_bytes(session_id):
    buf = io.BytesIO()
    save(session_id, buf)
    return buf.getvalue()
