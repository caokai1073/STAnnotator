"""The marker step's rare-type search (markers.RARE) on the bundled mouse colon example:
1. types no cluster of their level went to are found where their cells show their markers (finer clusters, or single
   cells), and no cell of them was found without it;
2. the types still without cells are listed in not_found, with the detection of their markers;
3. the run is the same twice (seeded Leiden, fixed rules); without the search, those types have no cells.
Run: python tests/smoke_rare.py  (no data folder needed)."""
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
from stagent import importer, markers

EX = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish"
r = importer.read_h5ad(str(EX / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), None, None, None, None, None, None)
genes, X = list(r["genes"]), sp.csr_matrix(r["X"])
tree = markers.parse(open(EX / "mouse_colon_merfish_markers.txt").read(), genes)["tree"]

a = markers.annotate(X, genes, tree, knn="exact")
b = markers.annotate(X, genes, tree, knn="exact")
assert np.array_equal(a["labels"], b["labels"]) and a["not_found"] == b["not_found"], "not the same twice"
markers.RARE = False
try:
    c = markers.annotate(X, genes, tree, knn="exact")
finally:
    markers.RARE = True

rare = [x for x in a["records"] if x.get("rare")]
got = {x["assigned"] for x in rare}
assert got, "the search found nothing"
n_with, n_without = pd.Series(a["labels"]).value_counts(), pd.Series(c["labels"]).value_counts()
for t in got:  # each type it found had no cells without it
    assert n_without.get(t, 0) == 0, (t, n_without.get(t))
print("found by the rare-type search:", {t: int(n_with[t]) for t in sorted(got)},
      f"({sum(1 for x in rare if x['rare'] == 'cells')} of {len(rare)} groups as single cells)")

# the evidence of every fine cluster taken: one positive marker raised in it
for nd in a["nodes"]:
    for t, tries in (nd.get("rare") or {}).items():
        for x in tries:
            if x["accepted"] and x.get("positive"):
                assert any(c_ >= 0.25 and c_ >= 3 * p_ or c_ >= 0.4 and c_ >= 1.5 * p_ for c_, p_ in x["positive"].values()), x

leaves = [t for t in markers.descendants(tree) if not markers._find(tree, t)["children"]]
missing = {t for t in leaves if n_with.get(t, 0) == 0}
assert missing == {x["type"] for x in a["not_found"]}, (missing, a["not_found"])
assert all(x["reason"] for x in a["not_found"])
print("not found:", ", ".join(sorted(missing)))
print("ok")
