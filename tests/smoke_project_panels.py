"""Slices join a project when they have the same genes in another order, and are told which genes differ otherwise:
1. a slice with the project's genes in another order (as CosMx exports can list them) is put in the project's order:
   the same counts, gene for gene, also its data as imported (original/, kept by a filter);
2. a few gene names that differ are counted and listed in the error; many are only counted;
3. a slice annotated in its own order is not reordered (an annotation's model uses the gene columns).
Run: python tests/smoke_project_panels.py  (the bundled example, in a temporary data folder removed at the end)."""
import json, os, shutil, tempfile
from pathlib import Path

HOME = Path(tempfile.mkdtemp(prefix="stagent_panels_"))
os.environ["STAGENT_HOME"] = str(HOME)
import numpy as np, scipy.sparse as sp  # noqa: E402
from stagent import config, importer, projects, qcfilter  # noqa: E402

H5 = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish" / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"
quiet = lambda *a: None


def counts(sid, sub=""):
    d = config.SAMPLES / sid / sub
    return sp.load_npz(d / "counts.npz").tocsr(), config.read_genes(d / "genes.txt")


def shuffle(sid, seed):  # the slice's genes in another order (its current data and the data as imported)
    from stagent.osutil import save_counts
    for sub in ("original", ""):
        d = config.SAMPLES / sid / sub
        if (d / "counts.npz").exists():
            X, g = counts(sid, sub)
            p = np.random.default_rng(seed).permutation(len(g))
            save_counts(d / "counts.npz", X[:, p])
            config.write_genes(d / "genes.txt", [g[j] for j in p])


try:
    for s in ("A", "B", "C", "D", "E"):
        importer.import_sample(s, str(H5), progress=quiet)
    qcfilter.apply(["A"], {"min_cells": 400})  # a filter keeps the data as imported in original/
    qcfilter.apply(["B"], {"min_cells": 400})
    XA, gA = counts("A")
    XA0, gA0 = counts("A", "original")
    assert len(gA) < len(gA0), "the filter removed no gene"

    # ---- 1. the same genes in another order
    shuffle("B", 1)
    XB, gB = counts("B")
    assert gB != gA and sorted(gB) == sorted(gA)
    assert counts("B", "original")[1] != gA0
    done = []
    projects.assign({"A": "P", "B": "P"}, done)
    assert done == ["B"], done
    XB, gB = counts("B")
    XB0, gB0 = counts("B", "original")
    assert gB == gA and gB0 == gA0, "not in the project's order"
    assert (XB != XA).nnz == 0 and (XB0 != XA0).nnz == 0, "the counts changed"
    assert XB.has_canonical_format
    print(f"same genes in another order: B put in A's order ({len(gA)} genes, {len(gA0)} as imported), counts equal")

    # ---- 2. differing gene names: few listed, many counted
    g = config.read_genes(config.SAMPLES / "C" / "genes.txt")
    config.write_genes(config.SAMPLES / "C" / "genes.txt", [f"{x}_v2" if i < 3 else x for i, x in enumerate(g)])
    try:
        projects.assign({"D": "R", "C": "R"})
        raise AssertionError("joined with other gene names")
    except ValueError as e:
        msg = str(e)
    print(msg)
    assert "6 gene names differ" in msg and "only in D: " + ", ".join(g[:3]) in msg and f"{g[0]}_v2" in msg
    config.write_genes(config.SAMPLES / "C" / "genes.txt", [f"{x}_v2" if i < 25 else x for i, x in enumerate(g)])
    try:
        projects.assign({"D": "R", "C": "R"})
        raise AssertionError("joined with other gene names")
    except ValueError as e:
        msg = str(e)
    print(msg)
    assert "50 gene names differ" in msg and "(25 only in D, 25 only in C)" in msg and "_v2" not in msg
    assert not (projects.PROJECTS / "R.json").exists(), "a refused project was written"
    config.write_genes(config.SAMPLES / "C" / "genes.txt", g)

    # ---- 3. an annotated slice keeps its order
    shuffle("E", 2)
    (config.SESSIONS / "E_x").mkdir(parents=True)
    json.dump({"session_id": "E_x", "sample_id": "E", "created_at": "2026-11-11T00:00:00", "method": "markers"},
              open(config.SESSIONS / "E_x" / "session.json", "w"))
    projects.assign({"E": "S", "D": "S"})  # E annotated: D (not annotated) takes E's order
    assert counts("D")[1] == counts("E")[1]
    print("an annotated slice keeps its order; the other slice takes it")
    shuffle("C", 3)
    (config.SESSIONS / "C_x").mkdir(parents=True)
    json.dump({"session_id": "C_x", "sample_id": "C", "created_at": "2026-11-11T00:00:00", "method": "markers"},
              open(config.SESSIONS / "C_x" / "session.json", "w"))
    try:
        projects.assign({"C": "S"})
        raise AssertionError("an annotated slice was reordered")
    except ValueError as e:
        assert "annotated in their own order" in str(e), e
    print("two slices annotated in different orders: refused")
    print("ok")
finally:
    shutil.rmtree(HOME, ignore_errors=True)
