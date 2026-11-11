"""Sparse count matrices built with little memory beyond their own size (a CosMx 6K, Visium HD or Stereo-seq matrix is
gigabytes): row blocks stacked into one matrix, columns picked block by block, rows kept without a copy when all are."""
import os
import numpy as np, scipy.sparse as sp


def stack_rows(data, ind, nnz, n_cols):
    """One CSR of row blocks given as their data, column indices and nonzeros per row. The lists are emptied while
    their blocks are copied in, so memory holds the matrix once, not the blocks and the matrix side by side. The same
    matrix as scipy's vstack of the blocks."""
    per_row = np.concatenate(nnz) if nnz else np.zeros(0, np.int64)
    total = int(per_row.sum())
    it = np.int32 if total < 2 ** 31 else np.int64
    indptr = np.zeros(len(per_row) + 1, it)
    np.cumsum(per_row, out=indptr[1:])
    out_d, out_i, at = np.empty(total, np.float32), np.empty(total, it), 0
    while data:
        d, i = data.pop(0), ind.pop(0)
        out_d[at:at + len(d)], out_i[at:at + len(d)] = d, i
        at += len(d)
        del d, i
    nnz.clear()
    return sp.csr_matrix((out_d, out_i, indptr), shape=(len(per_row), n_cols))


def columns_f32(X, keep, block=50000):
    """X[:, keep] as float32 CSR, made block by block: the matrix and its copy, not two copies besides it."""
    X = sp.csr_matrix(X)
    keep = np.asarray(keep)
    if len(keep) == X.shape[1] and (keep == np.arange(X.shape[1])).all():
        return X if X.dtype == np.float32 else X.astype(np.float32)
    data, ind, nnz = [], [], []
    for a in range(0, X.shape[0], block):
        m = X[a:a + block][:, keep].astype(np.float32).tocsr()
        data.append(m.data); ind.append(m.indices); nnz.append(np.diff(m.indptr))
    return stack_rows(data, ind, nnz, len(keep))


def rows(X, ok):
    """X[ok] (a boolean mask of rows); X itself when every row is kept."""
    ok = np.asarray(ok, bool)
    return X if ok.all() else X[ok]


LIVE = set()  # the folders of the RowBlocks of this process not closed yet


class RowBlocks:
    """Row blocks of a matrix read in chunks (a CSV of counts), kept on disk until the last one: then the matrix is
    read back in one piece, so memory holds it once, not its blocks and itself side by side (a CosMx 6K slide of
    two million cells is some 15 GB). folder: where the two temporary files go (the data folder: room for the
    matrix); they are removed when the matrix is made or the reading fails. A process killed meanwhile leaves its
    rows_* folder: cachemgr.prune_importing removes it later (owner.pid names the process that wrote it)."""

    def __init__(self, n_cols, folder=None):
        import tempfile
        if folder is not None:
            os.makedirs(folder, exist_ok=True)
        self.dir = tempfile.mkdtemp(prefix="rows_", dir=folder)
        LIVE.add(os.path.abspath(self.dir))
        with open(os.path.join(self.dir, "owner.pid"), "w", encoding="utf-8") as f:  # (cachemgr.prune_importing)
            f.write(str(os.getpid()))
        self.n_cols, self.nnz, self.total = n_cols, [], 0
        self.fd = open(os.path.join(self.dir, "data.bin"), "wb")
        self.fi = open(os.path.join(self.dir, "indices.bin"), "wb")

    def add(self, m):
        m = sp.csr_matrix(m)
        self.fd.write(np.ascontiguousarray(m.data, np.float32).data)  # (the array's bytes, not a copy of them)
        self.fi.write(np.ascontiguousarray(m.indices, np.int32).data)
        self.nnz.append(np.diff(m.indptr).astype(np.int64))
        self.total += m.nnz

    def matrix(self):
        try:
            self.fd.close(); self.fi.close()
            per_row = np.concatenate(self.nnz) if self.nnz else np.zeros(0, np.int64)
            it = np.int32 if self.total < 2 ** 31 else np.int64
            indptr = np.zeros(len(per_row) + 1, it)
            np.cumsum(per_row, out=indptr[1:])
            data = np.fromfile(os.path.join(self.dir, "data.bin"), np.float32)
            ind = np.fromfile(os.path.join(self.dir, "indices.bin"), np.int32).astype(it, copy=False)
            return sp.csr_matrix((data, ind, indptr), shape=(len(per_row), self.n_cols))
        finally:
            self.close()

    def close(self):
        import shutil
        for f in (self.fd, self.fi):
            try:
                f.close()
            except OSError:
                pass
        shutil.rmtree(self.dir, ignore_errors=True)
        LIVE.discard(os.path.abspath(self.dir))
