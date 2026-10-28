"""The data folder's caches, kept under a size limit (default 5 GB; settings.json "cache_limit_gb", or $STAGENT_CACHE_GB).

Caches are files the app can make again when they are needed: figures the agent looked at, notebook exports, saved
annotation files waiting for their download, R files converted to h5ad, the alignment figures, and the expression
model of earlier versions of an annotation (sessions/*/base: each set of corrections or re-clusterings has its own).
When they add up to more than the limit the oldest go first. Never touched: the cells and counts of a sample, the
versions of an annotation, its niches and results, anything the latest version of an annotation uses (its expression
model, marker step and the first model), the trash, and files changed in the last 10 minutes (perhaps still in use).
"""
import json, os, threading, time
from . import config

DEFAULT_GB = 5.0
RECENT_S = 600
KEEP_BASE = {"init.npy", "init_model.npz", "markers_init.npz", "markers_init.json", "reference_cells.npz", "overclusters.npy"}


def limit_gb():
    if os.environ.get("STAGENT_CACHE_GB"):
        return float(os.environ["STAGENT_CACHE_GB"])
    try:
        return float(json.load(open(config.HOME / "settings.json")).get("cache_limit_gb", DEFAULT_GB))
    except (OSError, ValueError, TypeError):
        return DEFAULT_GB


def set_limit_gb(gb):
    gb = float(gb)
    if not 0.5 <= gb <= 10000:
        raise ValueError("the cache limit is between 0.5 and 10000 GB")
    f = config.HOME / "settings.json"
    try:
        cur = json.load(open(f)) if f.exists() else {}
    except ValueError:
        cur = {}
    cur["cache_limit_gb"] = gb
    json.dump(cur, open(f, "w"), indent=1)
    return gb


def _latest_state(d):
    vs = sorted((d / "versions").glob("v*.json"))
    return json.load(open(vs[-1])) if vs else None


def _session_keep(d):
    """Base files the latest version of the annotation in folder d uses."""
    from .session import Session, _tag
    st = _latest_state(d)
    if st is None:
        return None  # unreadable: keep all of it
    keep = set(KEEP_BASE)
    keep.add(f"{Session._base_key(None, st)}.npy")
    ops = st.get("marker_step", {}).get("ops", [])
    for i in range(len(ops) + 1):
        k = Session._mkey(ops[:i])
        keep |= {f"markers_{k}.npz", f"markers_{k}.json", f"model_{k}.npz"}
    return keep


def files():
    """[(path, kind, size, mtime)] of every cache file."""
    out = []

    def add(p, kind):
        try:
            st = p.stat()
        except OSError:
            return
        out.append((p, kind, st.st_size, st.st_mtime))
    for d in sorted(config.SESSIONS.iterdir()) if config.SESSIONS.exists() else []:
        if not (d / "session.json").exists():
            continue
        for sub, kind in (("figures", "figures"), ("exports", "exports")):
            for p in (d / sub).rglob("*") if (d / sub).exists() else []:
                if p.is_file():
                    add(p, kind)
        try:
            keep = _session_keep(d)
        except Exception:
            keep = None
        if keep is not None and (d / "base").exists():
            for p in (d / "base").iterdir():
                if p.is_file() and p.name not in keep:
                    add(p, "models of earlier versions")
    for root, kind, pattern in ((config.HOME / "cache" / "bundles", "saved files", "*"), (config.HOME / "converted", "converted R files", "*"),
                                (config.SAMPLES, "alignment figures", "*/align/*.jpg")):
        for p in root.glob(pattern) if root.exists() else []:
            if p.is_file():
                add(p, kind)
    return out


def usage():
    fs = files()
    by = {}
    for _, kind, size, _ in fs:
        by[kind] = by.get(kind, 0) + size
    trash = sum(f.stat().st_size for f in config.TRASH.rglob("*") if f.is_file()) if config.TRASH.exists() else 0
    return {"bytes": sum(by.values()), "limit_gb": limit_gb(), "by_kind": by, "n_files": len(fs), "trash_bytes": trash}


def prune(limit=None, everything=False):
    """Delete the oldest cache files until they hold at most the limit (everything: all of them, but never files
    changed in the last 10 minutes). Returns (files deleted, bytes freed)."""
    cap = (limit_gb() if limit is None else float(limit)) * 1e9
    fs = sorted(files(), key=lambda x: x[3])
    total = sum(x[2] for x in fs)
    n = freed = 0
    now = time.time()
    for p, _, size, mtime in fs:
        if not everything and total - freed <= cap:
            break
        if now - mtime < RECENT_S:
            continue
        try:
            p.unlink()
            n += 1
            freed += size
        except OSError:
            pass
    return n, freed


def start(every_s=600):
    """Prune now and every every_s seconds in the background (the app's lifetime)."""
    def loop():
        while True:
            try:
                prune()
            except Exception as e:
                print(f"cache: {e}")
            time.sleep(every_s)
    threading.Thread(target=loop, daemon=True, name="cache-prune").start()
