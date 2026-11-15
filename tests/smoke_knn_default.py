"""The default neighbor search (session.knn_mode with knn None) is the same every time on a computer: it does not
follow the runs timed here (timings.json), the free memory or whether pynndescent is compiled in this process; the
estimate endpoint gives the app that same default. Run: python tests/smoke_knn_default.py (a temporary data folder)."""
import os, tempfile
home = tempfile.mkdtemp()
os.environ["STAGENT_HOME"] = home; os.environ["XDG_CONFIG_HOME"] = home + "/cfg"
from stagent import markers, timing
from stagent.session import knn_mode

assert knn_mode(None, 500_000, 480) == markers.KNN == "exact"  # below a million cells: exact
assert knn_mode("approx", 10) == "approx" and knn_mode("exact", 5_000_000) == "exact"
cases = [(1_200_000, 480), (3_000_000, 5000), (8_000_000, 480)]
before = [knn_mode(None, n, g, None, 20, n * 100) for n, g in cases]
# runs timed here that would flip the choice either way, and pynndescent compiled in this process
for key, units, sec in (("knn_exact_gpu", 1e14, 1.5), ("knn_exact_torch", 1e14, 1.5), ("knn_exact_numpy", 1e14, 1.5),
                        ("knn_exact_mlx", 1e14, 1.5), ("knn_approx", 1e6, 500.0)):
    for _ in range(5):
        timing.record(key, units, sec)
assert {"knn_exact_gpu", "knn_approx"} <= set(timing.measured()), timing.measured()
timing._warm.add("approx_compile")
timing.available_memory_gb = lambda: 0.5  # (the free memory now: not used by the default)
after = [knn_mode(None, n, g, None, 20, n * 100) for n, g in cases]
assert before == after, (before, after)
# while the estimate shown does follow those runs
est = timing.estimate(3_000_000, 5000, None, "auto", "exact", 20, 3e8)
assert est["basis"] != "from the reference workstation's rates"
print("default neighbor search:", dict(zip([n for n, _ in cases], before)), "| backend", markers.knn_backend("exact")[1])
print("ok")
