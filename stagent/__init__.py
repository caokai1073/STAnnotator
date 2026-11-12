"""Spatial cell type annotation engine with natural-language-driven constraints.

The LLM only translates user intent into tool calls; every annotation result is computed
deterministically here from (base probabilities, regions, constraints) and versioned.
Entry points: stagent.session.Session, stagent.mcp_server.
"""
__version__ = "0.2.56"  # one step up per release: 0.2.1, 0.2.2, ...

# The approximate neighbor search (pynndescent, markers.knn) always runs on KNN_THREADS threads, because its result
# depends on the number of threads: so that every computer gives the same annotation, numba must allow that many
# threads even with fewer cores (they then share the cores; nothing fails). Set before numba is first imported.
KNN_THREADS = 8  # = markers.KNN_THREADS (kept here too: this runs before anything imports numba)


def _numba_threads():
    import os
    try:
        have = int(os.environ.get("NUMBA_NUM_THREADS") or 0)
    except ValueError:
        have = 0
    os.environ["NUMBA_NUM_THREADS"] = str(max(have, os.cpu_count() or 1, KNN_THREADS))


_numba_threads()

# MLX (the exact neighbor search on the Apple GPU) multiplies float32 matrices in TF32 by default on the GPUs that have
# it (Apple M5): about 10 bits, which drops true neighbors and changes the annotation. Off: full float32, the same
# neighbors as every other computer. MLX reads it once, at its first product (markers._mlx_gpu checks the precision).
__import__("os").environ["MLX_ENABLE_TF32"] = "0"
