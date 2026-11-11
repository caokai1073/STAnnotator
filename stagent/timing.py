"""How long a marker annotation takes on this computer.

Every step of a run is timed and kept as a rate (seconds per unit of work, e.g. per million cells, or per 10^9 cell
pairs for the exact neighbor search) in timings.json in the data folder, one set per computer. An estimate before a
run multiplies the work the sample implies (cells, genes, the levels of the marker list, the neighbor search) by
those rates; a step never timed here uses the rate measured on the reference workstation (20 cores, RTX PRO 4500),
scaled by a short speed test of this computer (under a second). The first run here corrects the estimate.
"""
import json, platform, threading, time
from contextlib import contextmanager
import numpy as np
from . import config
from .osutil import dump_json

# seconds per unit on the reference workstation; units in the comments
REF = {"prepare": 0.65,         # per 10^6 cells: normalize, detection matrix
       "pca": 11.0,             # per 10^9 cells x genes
       "knn_exact_gpu": 0.055,  # per 10^9 cell pairs
       "knn_exact_torch": 1.93,
       "knn_exact_numpy": 0.21,  # (8 threads)
       "knn_exact_mlx": 0.1,    # Apple GPU (MLX; no reference workstation has one): 2.5 x the "mlx" speed test, which
                                # searches 20,000 points (0.4 x 10^9 pairs; small, so a rate on the slow side)
       "knn_rerank": 1.8,       # per 10^6 cells: float64 re-ranking of the exact candidates
       "knn_approx": 9.2,       # per 10^6 cells (pynndescent)
       "approx_compile": 10.0,  # once per app start: numba compiles pynndescent (less once its cache is on disk)
       "graph": 3.1,            # per 10^6 cells
       "leiden": 13.0,          # per 10^6 cells, per resolution tried
       "model": 430.0,          # per 10^12 training cells x genes x cell types
       "finish": 5.7}           # per 10^6 cells: the expression model on every cell, regions, priors, conflicts and
                                # whatever else is not timed separately
# which speed test scales each step: "cores" (matrix product on all cores), "one" (one core), or the neighbor
# search itself on a small set ("torch", "numpy", "gpu", "mlx"; "mlx" has no reference run: its REF_PROBE is nominal)
SCALE = {"knn_exact_torch": "torch", "knn_exact_numpy": "numpy", "knn_approx": "cores", "model": "cores",
         "knn_exact_gpu": "gpu", "knn_exact_mlx": "mlx"}
REF_PROBE = {"cores": 0.135, "one": 0.0164, "torch": 0.23, "numpy": 0.21, "gpu": 0.0113, "mlx": 0.04}  # their seconds there
VERSION = {"knn_exact_numpy": 2, "knn_rerank": 2, "knn_exact_mlx": 2}  # a step whose code changed: its rates timed before are dropped
MIN_SECONDS = 1.0  # shorter steps are not recorded: fixed overheads would dominate their rate
DECAY = 0.8        # earlier records count this much less with every new one; the rate is seconds / units of what
                   # is kept, so large runs weigh more than small ones (whose fixed overheads differ)

_lock = threading.Lock()
_local = threading.local()
_warm = set()  # one-time steps (approx_compile) already paid in this process
_probe = {}


def _file():
    return config.HOME / "timings.json"


def _host():
    return platform.node() or "this computer"


def _load():
    try:
        return json.load(open(_file()))
    except (OSError, ValueError):
        return {}


def _current(key, m):
    return bool(m) and m.get("v", 1) == VERSION.get(key, 1)


def measured():
    """{step: {"rate", "runs"}} timed on this computer (by the step's current code)."""
    return {k: m for k, m in _load().get(_host(), {}).items() if _current(k, m)}


def record(key, units, seconds):
    if seconds < MIN_SECONDS or units <= 0 or not getattr(config, "RECORD_TIMINGS", True):
        return
    with _lock:
        d = _load()
        mine = d.setdefault(_host(), {})
        old = mine.get(key) if _current(key, mine.get(key)) else {}
        sec = DECAY * old.get("seconds", 0.0) + seconds
        unit = DECAY * old.get("units", 0.0) + units
        mine[key] = {"rate": sec / unit, "seconds": sec, "units": unit, "runs": old.get("runs", 0) + 1}
        if key in VERSION:
            mine[key]["v"] = VERSION[key]
        try:
            dump_json(d, _file(), indent=1)
        except OSError:
            pass


@contextmanager
def step(key, units, warmup=None):
    """Time a step of `units` work. warmup: a one-time cost paid inside the first such step of the process (its
    time beyond the usual rate is recorded under that name)."""
    t = time.perf_counter()
    yield
    dt = time.perf_counter() - t
    if getattr(_local, "timed", None) is not None:
        _local.timed += dt
    if warmup and warmup not in _warm:
        _warm.add(warmup)
        rate = measured().get(key, {}).get("rate") or REF[key] * speed_factor(key)
        record(warmup, 1.0, max(dt - rate * units, 0.0))
        return
    record(key, units, dt)


@contextmanager
def run(n_cells):
    """A whole annotation of n_cells: the time not taken by the steps timed inside is recorded as "finish"."""
    _local.timed = 0.0
    t = time.perf_counter()
    try:
        yield
        record("finish", n_cells / 1e6, time.perf_counter() - t - _local.timed)
    finally:
        _local.timed = None


# ---------------- speed test ----------------
def _best(fn, reps=3):
    out = []
    for _ in range(reps):
        t = time.perf_counter()
        fn()
        out.append(time.perf_counter() - t)
    return min(out)


def _points(n):
    return np.random.default_rng(0).standard_normal((n, 30)).astype(np.float32)


def _probe_cores():
    A = np.random.default_rng(0).standard_normal((4096, 512)).astype(np.float32)
    return _best(lambda: A @ A.T)


def _probe_one():
    x = np.random.default_rng(0).random(2_000_000)
    return _best(lambda: np.sort(x.copy()))


def _probe_numpy():
    from .markers import _topk_numpy
    Z = _points(20000)
    return _best(lambda: _topk_numpy(Z, 24), 2)


def _probe_mlx():
    from .markers import _topk_mlx
    Z = _points(20000)
    return _best(lambda: _topk_mlx(Z, 24), 2)  # (the first run also compiles MLX's GPU kernels)


def _probe_knn(dev):
    import torch
    Z = torch.tensor(_points(50000), device=dev)

    def f():
        torch.topk(torch.cdist(Z[:4096], Z), 24, largest=False)
        if dev == "cuda":
            torch.cuda.synchronize()
        elif dev == "mps":
            torch.mps.synchronize()
    f()
    return _best(f)


def speed_factor(key):
    """How much slower this computer is than the reference workstation for a step (1 = as fast)."""
    kind = SCALE.get(key, "one")
    if kind not in _probe:
        try:
            if kind == "gpu":
                from .config import device
                _probe[kind] = _probe_knn(device())
            else:
                _probe[kind] = {"cores": _probe_cores, "one": _probe_one, "numpy": _probe_numpy, "mlx": _probe_mlx,
                                "torch": lambda: _probe_knn("cpu")}[kind]()
        except Exception:
            _probe[kind] = REF_PROBE[kind]
    return max(_probe[kind] / REF_PROBE[kind], 0.05)


def rate(key, mine=None):
    """(seconds per unit, "timed here" or "speed test") for a step."""
    m = (measured() if mine is None else mine).get(key)
    if m:
        return m["rate"], "timed here"
    return REF[key] * speed_factor(key), "speed test"


# ---------------- estimate ----------------
def level_sizes(n, tree, min_cells=60):
    """Cells clustered at each level of the marker step, before it runs: all cells at the top, and each type with
    subtypes is assumed to get an equal share of its parent's cells."""
    out = [n]

    def walk(nodes, m):
        share = m / max(len(nodes), 1)
        for nd in nodes:
            if nd["children"] and share >= min_cells:
                out.append(share)
                walk(nd["children"], share)
    walk(tree or [], n)
    return out


def estimate(n, n_genes, tree, resolution="auto", knn="exact", n_types=None):
    """Seconds a marker annotation of n cells should take here: {"seconds", "steps": [[name, seconds]], "basis",
    "backend"}. n_types: the number of cell types (leaves of the tree)."""
    from . import markers, model
    mine = measured()
    sizes = level_sizes(n, tree, markers.MIN_CELLS)
    n_res = len(markers.RESOLUTIONS) if resolution in (None, "", markers.AUTO) else 1
    key, backend = markers.knn_backend(knn)
    T = n_types or max(len(markers.descendants(tree or [])) // 2, 2)
    basis = set()

    def r(k):
        v, b = rate(k, mine)
        basis.add(b)
        return v
    nb = sum(r(key) * m * m / 1e9 + r("knn_rerank") * m / 1e6 for m in sizes) if knn != "approx" else \
        sum(r("knn_approx") * m / 1e6 for m in sizes) + (0 if "approx_compile" in _warm else r("approx_compile"))
    pool = min(n, model.CAP * T)
    steps = [["reading the cells", r("prepare") * n / 1e6],
             ["PCA", sum(r("pca") * m * n_genes / 1e9 for m in sizes)],
             ["neighbor graph", nb + sum(r("graph") * m / 1e6 for m in sizes)],
             ["Leiden clustering", sum(r("leiden") * n_res * m / 1e6 for m in sizes)],
             ["expression model", r("model") * pool * n_genes * T / 1e12],
             ["expression on every cell, regions, priors, conflicts", r("finish") * n / 1e6]]
    return {"seconds": round(sum(s for _, s in steps), 1), "steps": [[k, round(s, 1)] for k, s in steps],
            "levels": len(sizes), "backend": backend,
            "basis": "from runs on this computer" if basis == {"timed here"} else
            "from a speed test of this computer" if basis == {"speed test"} else
            "partly from runs on this computer"}


def words(seconds):
    """A duration for people: "under a minute", "about 4 min", "about 1 h 20 min"."""
    if seconds < 50:
        return "under a minute"
    m = max(1, round(seconds / 60))
    return f"about {m} min" if m < 60 else f"about {m // 60} h" + (f" {m % 60} min" if m % 60 else "")
