"""Many samples without the app: `stagent batch manifest.csv --out DIR` imports, registers and annotates every sample
of a manifest in a server's terminal and writes CSV files, deciding by itself where the app would ask (auto, the
default). Each sample runs in a process of its own (--jobs at once; --shard I/N for a cluster's array jobs) and
keeps a status file, so a stopped run goes on where it stopped and a sample that failed does not stop the others.

The manifest (CSV, one row per sample; only `sample` and `xenium` or `h5ad` are needed):
  sample        the sample's id (letters, digits, - _ .)
  xenium        a Xenium output folder              | or h5ad: an h5ad file (raw counts, x / y in um), with
                                                    |   section_column / section when it holds several sections,
                                                    |   x_column / y_column when the coordinates are not x / y
  he            the H&E image (optional; a pyramidal TIFF / OME-TIFF is read in place)
  registration  the H&E's alignment file (optional): Xenium Explorer's, an affine CSV, or cell positions saved by the
                app. Without it the cells are aligned to the H&E by the app's search (auto) or left for the app.
  organ         a built-in marker list (colon, small_intestine, lung, ...), pruned to the sample's panel
  markers       or a marker list of your own (CSV or text, as the app's Markers... dialog reads it)
  (--organ / --markers give one for every row that names none.)

Decisions in auto mode (--no-auto: no registration is chosen or corrected; those samples are flagged for the app):
  - registration without a file: align.search's best CANDIDATES candidates are refined and checked over the whole
    tissue (align.deviation); the one whose cells sit closest to their nuclei wins (clear readings first, then the
    median offset), "far off" ones never; two right-looking candidates in different orientations: flagged.
  - local shifts (tissue stretched or torn between the run and the staining, found by the registration check):
    the local correction (regcheck.correction) is applied when its own check moves the cells closer to their nuclei,
    else left and flagged.
  - annotation: the marker list as given; no agent (no Refine rounds).

DIR/data is the batch's data folder (open flagged samples in the app: STAGENT_HOME=DIR/data stagent remote-serve),
DIR/results/<sample>/ holds annotation.csv (one row per cell: cell type, lineage, levels, confidence, x_um / y_um),
cell_positions.csv (cell_id, x_um, y_um and px, py on the H&E), registration.csv (the affine, how it was made, the
check), qc.json; DIR/summary.csv one row per sample (status, flags, numbers); DIR/logs/<sample>.log.
"""
import csv, hashlib, json, os, re, shutil, sys, time, traceback
from pathlib import Path

CANDIDATES = 3          # alignment candidates refined and checked (each about a minute on a large section)
COLUMNS = ("sample", "xenium", "h5ad", "section_column", "section", "x_column", "y_column", "he", "registration",
           "organ", "markers")
SUMMARY = ("sample", "status", "flags", "error", "n_cells", "n_genes", "median_transcripts", "registration",
           "registration_verdict", "median_offset_um", "p90_offset_um", "share_off", "local_correction",
           "candidate", "n_types_with_cells", "types_without_cells", "unassigned_share", "seconds", "session")
STEPS = ("import", "register", "annotate", "export")
NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


def read_manifest(path, organ=None, markers=None):
    """The manifest's rows, checked before anything runs (every problem at once)."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        rows = [{(k or "").strip().lower(): (v or "").strip() for k, v in r.items()} for r in csv.DictReader(f)]
    errs, seen = [], set()
    if not rows:
        raise ValueError(f"{path}: no rows")
    unknown = sorted(set(rows[0]) - set(COLUMNS))
    if unknown:
        errs.append(f"unknown columns {unknown} (known: {', '.join(COLUMNS)})")
    for i, r in enumerate(rows, 2):
        s = r.get("sample", "")
        if not NAME.match(s):
            errs.append(f"row {i}: sample {s!r} (letters, digits, - _ . only)")
        elif s in seen:
            errs.append(f"row {i}: sample {s!r} twice")
        seen.add(s)
        if bool(r.get("xenium")) == bool(r.get("h5ad")):
            errs.append(f"row {i}: give xenium or h5ad (one of them)")
        for k in ("xenium", "h5ad", "he", "registration", "markers"):
            if r.get(k) and not Path(r[k]).expanduser().exists():
                errs.append(f"row {i}: {k} {r[k]} not found")
        r["organ"] = r.get("organ") or ("" if r.get("markers") else (organ or ""))
        r["markers"] = r.get("markers") or ("" if r.get("organ") else (markers or ""))
        if not r["organ"] and not r["markers"]:
            errs.append(f"row {i}: no marker list (organ or markers column, or --organ / --markers)")
    if errs:
        raise ValueError("the manifest has problems:\n  " + "\n  ".join(errs))
    return rows


def _hash(*parts):
    return hashlib.sha256(json.dumps(parts, sort_keys=True, default=str).encode()).hexdigest()[:16]


def _input_key(r):
    """What the import reads: the row's paths with their files' sizes and times (a replaced file runs it again)."""
    def stamp(p):
        if not p:
            return None
        p = Path(p).expanduser()
        try:
            st = p.stat()
            return [str(p.resolve()), st.st_size, int(st.st_mtime)]
        except OSError:
            return [str(p), None, None]
    return _hash({k: r.get(k) for k in ("xenium", "h5ad", "section_column", "section", "x_column", "y_column")},
                 [stamp(r.get(k)) for k in ("xenium", "h5ad", "he", "registration")])


def _write_atomic(path, text):
    path = Path(path)
    import socket
    tmp = path.with_name(f".{path.name}.{socket.gethostname()}.{os.getpid()}.tmp")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


# ---------------- one sample (in a worker process; STAGENT_HOME is the batch's data folder) ----------------
def _one(r, out, opts):
    out = Path(out)
    sid = r["sample"]
    stf = out / "status" / f"{sid}.json"
    try:
        st = json.load(open(stf))
    except (OSError, ValueError):
        st = {}
    st = {"sample": sid, "steps": st.get("steps", {}), "flags": [], "error": None}
    log = open(out / "logs" / f"{sid}.log", "a", encoding="utf-8", buffering=1)
    sys.stdout = sys.stderr = log
    print(f"==== {time.strftime('%Y-%m-%d %H:%M:%S')} sample {sid}", flush=True)
    t0 = time.time()
    try:
        keys = _keys(r, opts)
        stale = {}
        for step in STEPS:  # (a step runs again when what it reads changed, or a step it depends on runs again)
            stale[step] = (opts["redo"] or st["steps"].get(step, {}).get("key") != keys[step]
                           or any(stale[d] for d in DEPENDS[step]))
        exported = not stale["export"] and (out / "results" / sid / "annotation.csv").exists()
        if not _sample_there(sid) and not exported:  # (its data removed: --drop-data, or by hand) from the start
            stale = dict.fromkeys(STEPS, True)
        for step in STEPS:
            done = st["steps"].get(step, {})
            if not stale[step]:
                print(f"{step}: done before, kept", flush=True)
                st["flags"] += done.get("flags", [])
                continue
            st["steps"].pop(step, None)
            ts = time.time()
            print(f"{step} ...", flush=True)
            res = globals()[f"_step_{step}"](r, out, opts, st)
            st["steps"][step] = {"key": keys[step], "seconds": round(time.time() - ts), **res}
            st["flags"] += res.get("flags", [])
            _write_atomic(stf, json.dumps(st, indent=1, default=str))
        st["status"] = "flagged" if st["flags"] else "ok"
    except Exception as e:  # (one sample's failure never stops the batch)
        traceback.print_exc()
        st["status"], st["error"] = "failed", f"{type(e).__name__}: {e}"
    st["seconds"] = round(time.time() - t0)
    _write_atomic(stf, json.dumps(st, indent=1, default=str))
    print(f"==== {st['status']} in {st['seconds']} s {('; '.join(st['flags']))}", flush=True)
    log.close()
    return {"sample": sid, "status": st["status"], "flags": st["flags"], "error": st["error"], "seconds": st["seconds"]}


DEPENDS = {"import": (), "register": ("import",), "annotate": ("import",), "export": ("register", "annotate")}


def _keys(r, opts):
    """What each step reads: import the files; register those and the options of the decisions; annotate those and the
    marker list's content; export both."""
    k = {"import": _input_key(r)}
    k["register"] = _hash(k["import"], opts["auto"], opts["local_fix"], opts["candidates"])
    k["annotate"] = _hash(k["import"], _marker_text_key(r), opts["knn"] or "default")
    k["export"] = _hash(k["register"], k["annotate"])
    return k


def _sample_there(sid):
    from . import config
    return (config.SAMPLES / sid / "meta.json").exists()


def _marker_text_key(r):
    if r.get("markers"):
        return _hash(Path(r["markers"]).expanduser().read_text(encoding="utf-8"))
    return r.get("organ")


def _step_import(r, out, opts, st):
    from . import importer
    say = lambda *a: print(*a, flush=True)
    he = str(Path(r["he"]).expanduser()) if r.get("he") else None
    reg = str(Path(r["registration"]).expanduser()) if r.get("registration") else None
    if r.get("xenium"):
        importer.import_xenium(r["sample"], [str(Path(r["xenium"]).expanduser())], he, reg, overwrite=True, progress=say)
    else:
        importer.import_sample(r["sample"], str(Path(r["h5ad"]).expanduser()), he=he, registration=reg,
                               section_column=r.get("section_column") or None, section=r.get("section") or None,
                               x_column=r.get("x_column") or None, y_column=r.get("y_column") or None,
                               overwrite=True, progress=say)
    m = importer._meta(r["sample"])
    return {"n_cells": m.get("n_cells"), "n_genes": m.get("n_genes")}


def _check(sid):
    """The registration check of the cells' current place on the H&E (align.judge on regcheck's windows)."""
    from . import align, regcheck
    d = align.judge(regcheck.windows(sid))
    d.pop("offsets_um", None)
    return d


def _pick(cands):
    """The best checked candidate: never a "far off" one; clear readings first (the share of the coarse windows that
    find their nuclei), then the median offset. Ambiguous: another candidate in a different orientation that is not
    "far off" either (the check's own line between right and wrong placements, align.READ_SHARE)."""
    ok = [c for c in cands if c["check"]["verdict"] != "far off"]
    if not ok:
        return None, False
    ok.sort(key=lambda c: (-c["check"]["clear_reading_share"],
                           1e9 if c["check"]["median_um"] is None else c["check"]["median_um"]))
    best = ok[0]
    return best, any(c["orientation"] != best["orientation"] for c in ok[1:])


def _step_register(r, out, opts, st):
    from . import align, importer, regcheck
    sid = r["sample"]
    m = importer._meta(sid)
    if not r.get("he"):
        return {"registration": "none (no H&E)"}
    try:  # (a correction of an earlier try is undone: this one starts from the cells as imported / placed)
        importer.undo_registration_fix(sid)
    except ValueError:
        pass
    flags, res = [], {}
    if r.get("registration"):
        res["registration"] = "file"
    elif not opts["auto"]:
        a = align.search(align.Sample(sid))
        return {"registration": "waiting (no file; auto off)", "candidates": len(a["candidates"]),
                "flags": ["registration: no file, pick an alignment in the app"]}
    else:  # the app's search, decided here
        a = align.search(align.Sample(sid))
        cands = []
        for c in a["candidates"][:opts["candidates"]]:
            try:
                align.refine(align.Sample(sid), c["id"])
                d = align.deviation(align.Sample(sid), c["id"])
                cands.append({"id": c["id"], "orientation": c["orientation_name"], "score": c["score"],
                              "check": {k: v for k, v in d.items() if k not in ("A", "offsets_um")}})
                print(f"  candidate {c['id']} ({c['orientation_name']}): {d['verdict']}, median {d['median_um']} um, "
                      f"clear {d['clear_reading_share']:.0%}", flush=True)
            except ValueError as e:
                print(f"  candidate {c['id']}: {e}", flush=True)
        best, ambiguous = _pick(cands)
        res["candidates"] = cands
        if best is None:
            return {"registration": "failed (no candidate fits)", "candidates": cands,
                    "flags": ["registration: no alignment found, align in the app"]}
        align.apply(align.Sample(sid), best["id"])
        res |= {"registration": "auto", "candidate": f"{best['id']} {best['orientation']}"}
        if ambiguous:
            flags.append("registration: another orientation fits almost as well, check it")
    try:
        chk = _check(sid)
    except (ValueError, OSError) as e:  # (e.g. too few cells on the image for a window): flagged, the annotation goes on
        print(f"  registration check: {e}", flush=True)
        return res | {"verdict": None, "flags": flags + [f"registration: the check could not run ({e})"]}
    res["check"] = chk
    if chk["verdict"] == "local shifts" and opts["local_fix"] and opts["auto"]:
        try:
            res["local_correction"] = _local_fix(sid)
            if res["local_correction"].get("applied"):
                chk = res["check_after"] = _check(sid)
        except (ValueError, OSError) as e:
            print(f"  local correction: {e}", flush=True)
            res["local_correction"] = {"applied": False, "why": str(e)}
    if chk["verdict"] == "far off":
        flags.append("registration: the cells are far off their nuclei")
    elif chk["verdict"] == "local shifts":
        flags.append("registration: local shifts left")
    return res | {"verdict": chk["verdict"], "flags": flags}


def _local_fix(sid):
    """The local correction, applied only when its own check (the windows measured again with the cells moved) puts
    the cells closer to their nuclei than now (median offset and share of windows off; the windows of large shifts,
    which the fine windows cannot read, count as off before, as the app's check counts them). Applied as the app
    applies it: the placement before kept for its Undo, the sample noting the correction."""
    import numpy as np
    from . import importer, regcheck, config
    from .osutil import dump_json
    w = regcheck.windows(sid)
    fix = regcheck.correction(sid, w)
    a, mu, coef, aff, res, noise, D = regcheck._analyse(w)
    b, big, anchor = regcheck._large(w)
    tb = np.r_[np.hypot(*D[~noise].T), np.hypot(b[big, 4], b[big, 5]) if big.any() else []]
    ta = np.hypot(*np.asarray(fix["after_um"])[np.asarray(fix["trusted_after"])].T)
    if not len(tb) or not len(ta):
        return {"applied": False, "why": "no trusted windows"}
    mb, ma = float(np.median(tb)), float(np.median(ta))
    ob, oa = float((tb >= regcheck.OFF_UM).mean()), float((ta >= regcheck.OFF_UM).mean())
    out = {"median_before_um": round(mb, 1), "median_after_um": round(ma, 1), "off_before": round(ob, 3),
           "off_after": round(oa, 3), "moved_median_um": fix["moved_um_median"], "moved_max_um": fix["moved_um_max"]}
    if not (ma < mb and oa <= ob):
        return out | {"applied": False, "why": "it would not bring the cells closer to their nuclei"}
    sd = config.SAMPLES / sid
    keep = sd / "align" / "before_local"  # (the app's Undo of the local registration)
    keep.mkdir(parents=True, exist_ok=True)
    shutil.copy2(sd / "cells.parquet", keep / "cells.parquet")
    shutil.copy2(sd / "meta.json", keep / "meta.json")
    if regcheck._file(sid).exists():
        shutil.copy2(regcheck._file(sid), keep / "regcheck.json")
    importer.apply_registration_fix(sid, {k: fix[k] for k in regcheck.FIX_KEYS if k in fix})
    m = json.load(open(sd / "meta.json", encoding="utf-8"))
    m["registration_local"] = {"key": regcheck._key(sid, importer.image_of(sid)), "time": time.strftime("%Y-%m-%d %H:%M"),
                               "moved_um_max": fix["moved_um_max"], "moved_um_median": fix["moved_um_median"],
                               "by": "batch"}
    dump_json(m, sd / "meta.json", indent=1)
    return out | {"applied": True}


def _marker_text(r):
    from . import organs
    from .session import panel_genes
    if r.get("markers"):
        return Path(r["markers"]).expanduser().read_text(encoding="utf-8")
    return organs.draft(r["organ"], panel_genes(r["sample"]))["csv"]


def _step_annotate(r, out, opts, st):
    if os.environ.get("STAGENT_BATCH_TEST_CRASH") == r["sample"]:  # (tests: a process that dies, as out of memory)
        os._exit(9)
    import numpy as np
    from .session import Session
    s = Session.create(r["sample"], "none", method="markers", markers_text=_marker_text(r), knn=opts["knn"])
    lb = np.asarray(s.labels()).astype(str)
    from . import markers
    leaves = [t for t in s.types if t not in markers.NO_TYPE and not (markers._find(s.scheme["tree"], t) or {}).get("children")]
    have = set(lb.tolist())
    return {"session": s.id, "n_types_with_cells": len(have - set(markers.NO_TYPE)),
            "types_without_cells": [t for t in leaves if t not in have],
            "unassigned_share": round(float(np.isin(lb, markers.NO_TYPE).mean()), 4)}


def _step_export(r, out, opts, st):
    import numpy as np
    from . import notebook, importer
    from .session import Session
    sid = r["sample"]
    d = out / "results" / sid
    d.mkdir(parents=True, exist_ok=True)
    s = Session(st["steps"]["annotate"]["session"])
    tab = s.export_table()
    tmp = d / ".annotation.csv.tmp"
    tab.to_csv(tmp, index=False)
    os.replace(tmp, d / "annotation.csv")
    _write_atomic(d / "cell_positions.csv", notebook.positions_csv(sid))
    reg0 = st["steps"].get("register", {})
    if r.get("he") and reg0.get("registration") in ("file", "auto"):  # (placed on the H&E: not waiting, not failed)
        reg = st["steps"].get("register", {})
        import pandas as pd, io
        t = pd.read_csv(io.StringIO(notebook.registration_csv(sid)))
        chk = reg.get("check_after") or reg.get("check") or {}
        for k in ("verdict", "median_um", "p90_um", "share_off"):
            t[f"check_{k}"] = chk.get(k)
        t["local_correction"] = (reg.get("local_correction") or {}).get("applied", False)
        _write_atomic(d / "registration.csv", t.to_csv(index=False))
    tc = np.asarray(s.counts.sum(1)).ravel()
    qc = {"sample": sid, "n_cells": int(len(tab)), "n_genes": len(s.panel), "median_transcripts": float(np.median(tc)),
          "steps": {k: {x: y for x, y in v.items() if x not in ("key",)} for k, v in st["steps"].items()},
          "flags": st["flags"]}
    _write_atomic(d / "qc.json", json.dumps(qc, indent=1, default=str))
    if opts["drop_data"]:  # (results written: the sample's data folder and annotations go, to save disk)
        from . import config
        for p in (config.SAMPLES / sid, config.SESSIONS / st["steps"]["annotate"]["session"]):
            shutil.rmtree(p, ignore_errors=True)
    return {"median_transcripts": qc["median_transcripts"]}


# ---------------- the batch (the parent process) ----------------
def summary(out):
    """DIR/summary.csv from every sample's status file."""
    out = Path(out)
    rows = []
    for f in sorted((out / "status").glob("*.json")):
        try:
            st = json.load(open(f))
        except (OSError, ValueError):
            continue
        s = st.get("steps", {})
        reg, ann, imp, exp = (s.get(k, {}) for k in ("register", "annotate", "import", "export"))
        chk = reg.get("check_after") or reg.get("check") or {}
        rows.append({"sample": st["sample"], "status": st.get("status", "running"), "flags": "; ".join(st.get("flags") or []),
                     "error": st.get("error") or "", "n_cells": imp.get("n_cells"), "n_genes": imp.get("n_genes"),
                     "median_transcripts": exp.get("median_transcripts"), "registration": reg.get("registration"),
                     "registration_verdict": chk.get("verdict"), "median_offset_um": chk.get("median_um"),
                     "p90_offset_um": chk.get("p90_um"), "share_off": chk.get("share_off"),
                     "local_correction": (reg.get("local_correction") or {}).get("applied"),
                     "candidate": reg.get("candidate"), "n_types_with_cells": ann.get("n_types_with_cells"),
                     "types_without_cells": "; ".join(ann.get("types_without_cells") or []),
                     "unassigned_share": ann.get("unassigned_share"), "seconds": st.get("seconds"),
                     "session": ann.get("session")})
    import io
    buf = io.StringIO()
    w = csv.DictWriter(buf, SUMMARY)
    w.writeheader()
    w.writerows(rows)
    _write_atomic(out / "summary.csv", buf.getvalue())
    return rows


def run(manifest, out, jobs=1, shard=None, auto=True, redo=False, organ=None, markers=None, drop_data=False,
        local_fix=True, candidates=CANDIDATES, only=None, knn=None):
    rows = read_manifest(manifest, organ, markers)
    if shard:  # (by the row in the whole manifest: a sample's shard never depends on --only)
        i, n = shard
        rows = [r for k, r in enumerate(rows) if k % n == i]
    if only:
        known = {r["sample"] for r in read_manifest(manifest, organ, markers)}
        missing = sorted(set(only) - known)
        if missing:
            raise ValueError(f"--only: not in the manifest: {', '.join(missing)}")
        rows = [r for r in rows if r["sample"] in set(only)]
    out = Path(out).expanduser().resolve()
    for sub in ("data", "results", "status", "logs"):
        (out / sub).mkdir(parents=True, exist_ok=True)
    jobs = max(1, int(jobs))
    threads = max(1, (os.cpu_count() or 1) // jobs)
    env = {"STAGENT_HOME": str(out / "data"), "STAGENT_NO_QC_BACKFILL": "1"}
    for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[v] = os.environ.get(v) or str(threads)
    if jobs > 1:  # several samples at once: their threads share the cores (numba's count is set by stagent itself,
        env["NUMBA_NUM_THREADS"] = str(max(threads, 8))  # at least the 8 the approximate search needs), and no Leiden
        env["STAGENT_LEIDEN_PARALLEL"] = "0"              # worker processes inside each
    os.environ.update(env)  # (inherited by the worker processes, before they import the app)
    opts = {"auto": auto, "redo": redo, "drop_data": drop_data, "local_fix": local_fix, "candidates": int(candidates),
            "knn": knn}
    print(f"{len(rows)} samples, {jobs} at a time ({threads} threads each); data folder {out / 'data'}; "
          f"auto {'on' if auto else 'off'}", flush=True)
    import multiprocessing
    ctx = multiprocessing.get_context("spawn")
    t0, n_done, counts = time.time(), 0, {}
    todo, running = list(rows), {}
    # one process per sample (its memory given back at its end); a process that dies (out of memory, a crash) fails
    # only its own sample: the others go on
    while todo or running:
        while todo and len(running) < jobs:
            r = todo.pop(0)
            pr = ctx.Process(target=_one, args=(r, str(out), opts), name=f"batch-{r['sample']}")
            pr.start()
            running[r["sample"]] = pr
        time.sleep(1)
        for sid, pr in list(running.items()):
            if pr.is_alive():
                continue
            pr.join()
            del running[sid]
            n_done += 1
            x = _result(out, sid, pr.exitcode)
            counts[x["status"]] = counts.get(x["status"], 0) + 1
            el = time.time() - t0
            print(f"[{n_done}/{len(rows)}] {sid}: {x['status']}"
                  + (f" ({'; '.join(x['flags'])})" if x["flags"] else "") + (f" {x['error']}" if x["error"] else "")
                  + f"  [{el / 60:.0f} min, about {el / n_done * (len(rows) - n_done) / 60:.0f} min left]", flush=True)
            summary(out)
    summary(out)
    print(f"done: {counts}; {out / 'summary.csv'}", flush=True)
    return counts


def _result(out, sid, exitcode):
    """A finished sample's outcome from its status file; a process that died before writing it (killed for memory, a
    crash): failed, with the steps it had finished kept for the next run."""
    f = out / "status" / f"{sid}.json"
    try:
        st = json.load(open(f))
    except (OSError, ValueError):
        st = {"sample": sid, "steps": {}}
    if exitcode != 0 or st.get("status") not in ("ok", "flagged", "failed"):
        why = (f"the process died (exit code {exitcode}: out of memory?)" if exitcode
               else "the process ended without a result")
        st |= {"status": "failed", "error": why, "flags": st.get("flags") or []}
        _write_atomic(f, json.dumps(st, indent=1, default=str))
    return {"status": st["status"], "flags": st.get("flags") or [], "error": st.get("error")}
