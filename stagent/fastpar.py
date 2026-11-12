"""Work spread over threads with the same result as one item after the other, and the loops that let it.

pmap() runs independent items on up to WORKERS threads: each item is computed on its own exactly as in a loop, so
the results do not depend on the threads (only the time does). It pays when the work releases Python's lock: numpy's
sorts and arithmetic, and the numba loops below (nogil).

The loops repeat a numpy / scipy computation in the same order, so they give the same numbers bit for bit (checked
by tests/smoke_fastpar.py):
- lr_edges(): the neighbor-pair products of a ligand-receptor pair and their sums per group pair, as
  communication._Graph.score did with fancy indexing and np.bincount (which holds Python's lock);
- col_means(): the column means of some rows of a CSR matrix, as scipy's X[rows].mean(0) (float32 sums over the rows
  in order), without copying those rows.
"""
import os
import numpy as np

WORKERS = 8
_K = {}


def workers():
    return max(1, min(WORKERS, os.cpu_count() or 1))


def pmap(fn, items, n=None, done=None):
    """[fn(x) for x in items] on up to n threads (default workers()), in the order of items. done(k) is called in the
    caller's thread as items finish (progress). A stop of the run (jobs) reaches the caller within half a second;
    items not started then are dropped."""
    items = list(items)
    n = min(n or workers(), len(items))
    if n <= 1:
        out = []
        for i, x in enumerate(items):
            out.append(fn(x))
            if done:
                done(i + 1)
        return out
    from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
    ex = ThreadPoolExecutor(n, thread_name_prefix="stagent-pmap")
    try:
        fs = [ex.submit(fn, x) for x in items]
        pending, k = set(fs), 0
        while pending:
            fin, pending = wait(pending, timeout=0.5, return_when=FIRST_COMPLETED)  # (a timeout: a stop gets through)
            for f in fin:
                f.result()  # an error: at once
            k += len(fin)
            if done and fin:
                done(k)
        return [f.result() for f in fs]
    finally:
        ex.shutdown(wait=False, cancel_futures=True)


def _compiled():
    if not _K:
        from numba import njit
        cache = bool(globals().get("__file__")) and os.path.isfile(globals()["__file__"])  # (not in a notebook's sa)

        @njit(nogil=True, cache=cache)
        def lr(src, dst, w, pair, L, R, nb, v, wv):
            acc = np.zeros(nb)
            for e in range(len(src)):
                x = L[src[e]] * R[dst[e]]    # float32, as numpy's L[src] * R[dst]
                v[e] = x
                y = w[e] * x                 # float32, as w * v
                wv[e] = y
                acc[pair[e]] += y            # in float64, in order, as np.bincount(pair, wv)
            return acc

        @njit(nogil=True, cache=cache)
        def cols(indptr, indices, data, data2, rows, ncol, s):
            a = np.zeros(ncol, np.float32)
            b = np.zeros(ncol, np.float32)
            c = np.zeros(ncol, np.int64)
            for r in rows:
                for j in range(indptr[r], indptr[r + 1]):
                    k = indices[j]
                    a[k] += data[j] * s      # float32: each value times 1/n, summed over the rows in order
                    b[k] += data2[j] * s
                    c[k] += 1
            return a, b, c
        _K.update(lr=lr, cols=cols)
    return _K


def lr_edges(src, dst, w, pair, L, R, nb):
    """(v, wv, sums): v = L[src] * R[dst], wv = w * v (float32) and the sums of wv per pair value (float64, length
    nb): the same numbers as numpy's products and np.bincount(pair, wv, nb)."""
    v, wv = np.empty(len(src), np.float32), np.empty(len(src), np.float32)
    L, R = np.ascontiguousarray(L, np.float32), np.ascontiguousarray(R, np.float32)
    acc = _compiled()["lr"](src, dst, np.asarray(w, np.float32), pair, L, R, int(nb), v, wv)
    return v, wv, acc


def col_means(X, rows, data2):
    """For the rows of the CSR matrix X (float32, rows ascending): the mean of every column of X, the mean of data2
    (another array of X's nonzeros, e.g. expm1(X.data)) and the share of rows with a nonzero, as scipy's
    X[rows].mean(0), the same of a copy with data2 and X[rows].getnnz(axis=0) / len(rows)."""
    rows = np.ascontiguousarray(rows, np.int64)
    s = np.float32(1.0 / len(rows))  # scipy: (X * (1.0 / n)).sum(0), float32
    a, b, c = _compiled()["cols"](X.indptr, X.indices, X.data, data2, rows, X.shape[1], s)
    return a, b, c / max(len(rows), 1)
