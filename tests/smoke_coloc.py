"""Spatial co-localization (coloc.run) on synthetic slices with a known layout:
1. two groups in separate halves of a slice sit apart (z < 0), each next to itself (z > 0); a third group spread over
   the whole slice is next to neither more than by chance;
2. the result is symmetric and the same on every run (seeded shuffles), whatever the number of threads;
3. pooled slices: neighbor pairs never cross slices (two slices placed on top of each other add no pairs between
   them), and labels are shuffled within each slice;
4. `within` keeps only the cells it names; an impossible radius or a single group is refused.
Run: python tests/smoke_coloc.py  (no data folder needed)."""
import numpy as np, pandas as pd
from types import SimpleNamespace
from stagent import coloc, fastpar

rng = np.random.default_rng(1)


def slice_cells(n, seed):
    r = np.random.default_rng(seed)
    xy = r.uniform(0, 1000, (n, 2))
    lab = np.where(xy[:, 0] < 500, "left", "right").astype(object)
    mix = r.random(n) < 0.2
    lab[mix] = "spread"
    return pd.DataFrame({"x": xy[:, 0], "y": xy[:, 1]}), lab


def table(slices):
    sess = [SimpleNamespace(cells=c) for c, _ in slices]
    labels = pd.DataFrame({"cell type": np.concatenate([l for _, l in slices]),
                           "slice": np.concatenate([[f"s{i}"] * len(c) for i, (c, _) in enumerate(slices)])})
    return SimpleNamespace(sessions=sess, labels=labels, missing={})


one = table([slice_cells(4000, 2)])
r, t = coloc.run(one, "cell type", radius_um=30, n_perm=100)
g = {x: i for i, x in enumerate(r["groups"])}
Z = np.array(r["z"])
assert Z[g["left"], g["right"]] < -10, Z
assert Z[g["left"], g["left"]] > 5 and Z[g["right"], g["right"]] > 5
assert abs(Z[g["spread"], g["left"]]) < 4 and abs(Z[g["spread"], g["right"]]) < 4, Z
assert np.allclose(Z, Z.T, equal_nan=True)
assert len(t) == 6 and set(t.columns) >= {"group_a", "group_b", "observed", "expected", "z", "fdr"}
print(f"layout found: left/right z {Z[g['left'], g['right']]:.1f}, left/left {Z[g['left'], g['left']]:.1f}, "
      f"spread/left {Z[g['spread'], g['left']]:.1f}")

old = fastpar.WORKERS
fastpar.WORKERS = 1
r1, t1 = coloc.run(one, "cell type", radius_um=30, n_perm=100)
fastpar.WORKERS = old
assert r1["z"] == r["z"] and t1.equals(t), "the result depends on the threads"
print("the same numbers on one thread and on several")

# two slices on top of each other: the pooled count of pairs is the sum of the slices' own
a, b = slice_cells(3000, 3), slice_cells(3000, 4)
ra, _ = coloc.run(table([a]), "cell type", radius_um=25, n_perm=50)
rb, _ = coloc.run(table([b]), "cell type", radius_um=25, n_perm=50)
rp, _ = coloc.run(table([a, b]), "cell type", radius_um=25, n_perm=50)
assert rp["n_pairs"] == ra["n_pairs"] + rb["n_pairs"], (rp["n_pairs"], ra["n_pairs"], rb["n_pairs"])
assert rp["scope"] == "project"
print(f"pooled slices: {rp['n_pairs']:,} pairs = {ra['n_pairs']:,} + {rb['n_pairs']:,} (none across slices)")

rw, _ = coloc.run(table([a, b]), "cell type", within={"slice": "s0"}, radius_um=25, n_perm=50)
assert rw["n_pairs"] == ra["n_pairs"] and rw["n_cells"] == ra["n_cells"]
for bad in (dict(radius_um=0), dict(radius_um=500), dict(n_perm=5)):
    try:
        coloc.run(one, "cell type", **bad)
        raise AssertionError(f"accepted {bad}")
    except ValueError:
        pass
try:
    coloc.run(one, "cell type", within={"cell type": "left"})
    raise AssertionError("one group accepted")
except ValueError:
    pass
print("within keeps its cells; a bad radius, too few shuffles and a single group are refused")
print("ok")
