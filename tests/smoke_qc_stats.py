"""The filter dialog's numbers saved beside the counts (qcfilter.STATS): an import writes them, the dialog and a preview
without a minimum of cells per gene then read no matrix, and every result is the same as from the matrix itself;
the numbers follow the data as imported into original/ when a filter keeps it; changed counts are read again.
Runs in a temporary data folder."""
import os, shutil, tempfile, json
import numpy as np, scipy.sparse as sp
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_qcstats_")
os.environ["STAGENT_HOME"] = home; os.environ["XDG_CONFIG_HOME"] = home + "/cfg"
from stagent import config, importer, qcfilter, osutil

quiet = lambda *a: None
try:
    h5 = os.path.join(home, "one.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    importer.import_sample("s1", h5, progress=quiet)
    d = config.SAMPLES / "s1"
    assert (d / qcfilter.STATS).exists(), "the import saved no numbers"
    # the dialog from the saved numbers: no matrix read
    qcfilter._data.clear(); qcfilter._light.clear()
    ov = qcfilter.overview(["s1"])
    assert not qcfilter._data, "the dialog read the matrix"
    pv = qcfilter.preview(["s1"], {"min_counts": 20})
    assert not qcfilter._data, "a preview without min_cells read the matrix"
    # the same from the matrix itself (the saved numbers removed)
    (d / qcfilter.STATS).unlink(); qcfilter._data.clear(); qcfilter._light.clear()
    assert json.dumps(qcfilter.overview(["s1"]), sort_keys=True) == json.dumps(ov, sort_keys=True)
    assert json.dumps(qcfilter.preview(["s1"], {"min_counts": 20}), sort_keys=True) == json.dumps(pv, sort_keys=True)
    assert (d / qcfilter.STATS).exists(), "reading the matrix saved no numbers"
    qcfilter._data.clear(); qcfilter._light.clear()
    p2 = qcfilter.preview(["s1"], {"min_counts": 20, "min_cells": 50})  # (min_cells with cells removed: the matrix)
    qcfilter._data.clear(); qcfilter._light.clear()
    (d / qcfilter.STATS).unlink()
    assert json.dumps(qcfilter.preview(["s1"], {"min_counts": 20, "min_cells": 50}), sort_keys=True) == json.dumps(p2, sort_keys=True)
    # a filter keeps the data as imported in original/, with its numbers
    qcfilter.apply(["s1"], {"min_counts": 20, "min_cells": 50})
    assert (d / "original" / qcfilter.STATS).exists()
    qcfilter._data.clear(); qcfilter._light.clear()
    ov2 = qcfilter.overview(["s1"])
    assert not qcfilter._data and ov2["n_cells"] == ov["n_cells"] and ov2["counts"] == ov["counts"], "original/ numbers"
    # counts changed (another import, a gene reorder): the saved numbers no longer match and are not used
    X = sp.load_npz(d / "original" / "counts.npz").tocsr()
    X = X[: X.shape[0] // 2]
    osutil.save_counts(d / "original" / "counts.npz", X)
    qcfilter._data.clear(); qcfilter._light.clear()
    assert qcfilter.overview(["s1"])["n_cells"] == X.shape[0], "stale numbers used"
    print("ok")
finally:
    shutil.rmtree(home, ignore_errors=True)
