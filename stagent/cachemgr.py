"""The data folder's caches and trash, each kept under a size limit (settings.json "cache_limit_gb", default 5 GB, or
$STAGENT_CACHE_GB; "trash_limit_gb", default 10 GB, or $STAGENT_TRASH_GB).

Caches are files the app can make again when they are needed: figures the agent looked at, notebook exports, saved
annotation files waiting for their download, R files converted to h5ad, the alignment figures, and the expression
model of earlier versions of an annotation (sessions/*/base: each set of corrections or re-clusterings has its own).
When they add up to more than the limit the oldest go first. Never touched: the cells and counts of a sample, the
versions of an annotation, its niches and results, anything the latest version of an annotation uses (its expression
model, marker step and the first model), the trash, and files changed in the last 10 minutes (perhaps still in use).

The trash (trash/ of the data folder) holds what Reset, Remove sample and the cell filter moved away, one folder per
removal named by its time; past its limit the oldest removals are deleted for good, whole, never one of the last
10 minutes.
"""
import json, os, threading, time
from . import config
from .osutil import dump_json

DEFAULT_GB = 5.0
DEFAULT_TRASH_GB = 10.0
RECENT_S = 600
KEEP_BASE = {"init.npy", "init_model.npz", "markers_init.npz", "markers_init.json", "reference_cells.npz", "overclusters.npy"}


def _setting(key, env, default):
    if os.environ.get(env):
        return float(os.environ[env])
    try:
        return float(json.load(open(config.HOME / "settings.json")).get(key, default))
    except (OSError, ValueError, TypeError):
        return default


def limit_gb():
    return _setting("cache_limit_gb", "STAGENT_CACHE_GB", DEFAULT_GB)


def trash_limit_gb():
    return _setting("trash_limit_gb", "STAGENT_TRASH_GB", DEFAULT_TRASH_GB)


def _store(key, gb, what):
    gb = float(gb)
    if not 0.5 <= gb <= 10000:
        raise ValueError(f"the {what} limit is between 0.5 and 10000 GB")
    f = config.HOME / "settings.json"
    try:
        cur = json.load(open(f)) if f.exists() else {}
    except ValueError:
        cur = {}
    cur[key] = gb
    dump_json(cur, f, indent=1)
    return gb


def set_limit_gb(gb):
    return _store("cache_limit_gb", gb, "cache")


def set_trash_limit_gb(gb):
    return _store("trash_limit_gb", gb, "trash")


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
    tr = trash_entries()
    return {"bytes": sum(by.values()), "limit_gb": limit_gb(), "by_kind": by, "n_files": len(fs),
            "trash_bytes": sum(x[1] for x in tr), "trash_entries": len(tr), "trash_limit_gb": trash_limit_gb(),
            "trash_path": str(config.TRASH)}


def _size(p):
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file() and not f.is_symlink())


def removed_at(p):
    """When an entry went into the trash: the time its name starts with (every removal is named
    <YYYYmmdd-HHMMSS>_...), else None. The files inside keep their own, older times when they are moved."""
    import re
    m = re.match(r"(\d{8}-\d{6})_", p.name)
    if not m:
        return None
    try:
        return time.mktime(time.strptime(m.group(1), "%Y%m%d-%H%M%S"))
    except ValueError:
        return None


def trash_entries():
    """[(path, bytes, time)] of the removals in the trash, oldest first (time: when it was removed, from its name;
    for an entry not named so, the newest file inside)."""
    out = []
    for p in config.TRASH.iterdir() if config.TRASH.exists() else []:
        try:
            t = removed_at(p)
            if t is None:
                files = [p] if p.is_file() else [f for f in p.rglob("*") if f.is_file()]
                t = max((f.stat().st_mtime for f in files), default=p.stat().st_mtime)
            out.append((p, _size(p), t))
        except OSError:
            continue
    return sorted(out, key=lambda x: x[2])


def prune_trash(limit=None, everything=False):
    """Delete the oldest removals until the trash holds at most the limit (everything: all of it, but never one of
    the last 10 minutes). Returns (removals deleted, bytes freed)."""
    import shutil
    cap = (trash_limit_gb() if limit is None else float(limit)) * 1e9
    tr = trash_entries()
    total = sum(x[1] for x in tr)
    n = freed = 0
    now = time.time()
    for p, size, t in tr:
        if not everything and total - freed <= cap:
            break
        if now - t < RECENT_S:
            continue
        try:
            shutil.rmtree(p) if p.is_dir() else p.unlink()
            n += 1
            freed += size
        except OSError:
            pass
    return n, freed


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
    """Prune the caches and the trash now and every every_s seconds in the background (the app's lifetime)."""
    def loop():
        while True:
            for fn in (prune, prune_trash):
                try:
                    fn()
                except Exception as e:
                    print(f"{fn.__name__}: {e}")
            time.sleep(every_s)
    threading.Thread(target=loop, daemon=True, name="cache-prune").start()
