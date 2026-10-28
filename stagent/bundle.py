"""An annotation as one file (.stannot, a zip): save it, and load it again here or in another data folder.

The file holds the annotation's folder (every version with its regions, priors and corrections, the marker step and
expression model, niches with their names, review verdicts and advice, area suggestions, the report, the last result
of each analysis: DE genes, top genes, composition, cell-cell communication) and the sample's PCA / UMAP / Leiden
clusters. Chats are left out (they are conversations with Claude, not the annotation), and so are exports and
figures (made again on demand).

Loading needs the same data: the sample with the same cell ids and genes in the same order (a filtered sample matches
only the same filter). The annotation comes back as a new annotation of that sample; a different PCA / UMAP of the sample is
moved to the trash before the saved one replaces it.
"""
import hashlib, io, json, shutil, time, uuid, zipfile
from pathlib import PurePosixPath
from . import config

FORMAT = 1
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
            "results": sorted(p.stem for p in (s.dir / "results").glob("*.json")) if (s.dir / "results").exists() else []}
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED, compresslevel=4) as z:
        z.writestr("bundle.json", json.dumps(info, indent=1))
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
        if p.is_absolute() or ".." in p.parts or p.parts[0] not in ("annotation", "embedding") or p.suffix not in TYPES:
            raise ValueError(f"unexpected file in the annotation file: {n}")
    if "annotation/session.json" not in z.namelist():
        raise ValueError("the annotation file has no session.json")
    return info


def inspect(src):
    with zipfile.ZipFile(src) as z:
        return _check(z)


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
                             + ". Import that sample (with the same cell filter) first."
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
            sj = json.loads((d / "session.json").read_text())
            sj.update(session_id=sid, sample_id=target, loaded_from={"session_id": info["session_id"], "sample_id": info["sample_id"],
                                                                      "saved_at": info["saved_at"]})
            (d / "session.json").write_text(json.dumps(sj, indent=1))
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
