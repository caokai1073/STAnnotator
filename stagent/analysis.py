"""Differential expression and composition, for one slice or for all slices of a project.

Every analysis works on a table of cells and the groupings they belong to: cell type, lineage, region (drawn by the
user; "no region" outside them), niche, leiden (clusters of PCA + UMAP), selection (cells selected in the viewer),
slice, and the slices' metadata columns (e.g. condition). Scope "slice" is one annotation; scope "project" is the
latest annotation of every slice of a project, pooled.

de(): genes that differ between a group of cells and the rest (or another group), optionally within a subset (e.g.
  fibroblasts only). Wilcoxon rank-sum test of the log-normalized expression on up to MAX_TEST_CELLS cells per side
  (random, fixed seed); means and detection rates on all cells. Effect size: AUC (the chance that a cell of the
  group expresses the gene more than a cell of the reference; 0.5 = no difference) and log2 fold change of the
  normalized counts. Pooled slices: the direction of the change is checked in every slice that has both groups, or,
  when the groups are made of different slices, the slices of the group are compared with those of the reference.
de_all(): de() of every group of a grouping against the other cells, in one run (group "all"); the test on up to
  ALL_TEST_CELLS cells per side, so that many groups take seconds to a minute or two.
deseq2(): DESeq2 (PyDESeq2) on pseudobulk counts: a group against the rest or another group, the raw counts summed per
  replicate (a slice, or a patient / donor column of the project's metadata) and side; needs several replicates on
  each side, so it is offered for a project's slices only. Paired (~ replicate + side) when the replicates hold both
  sides, else ~ side.
markers(): the top genes of every group of a grouping against the rest of the cells.
crosstab(): how the cells of one grouping split over another (cell types per niche, niches per slice, ...), with the
  enrichment over what independence predicts and, for slice-level rows, every slice's own shares.
proportions(): for slice-level rows (a condition column of the project's slices), whether each value of cols (cell
  type, niche, ...) makes up a different share of the cells between conditions, with the replicates (slices, or
  patients) as the samples, as propeller (Phipson et al. 2022): per replicate the share of each value,
  arcsine-square-root transformed; per value a linear model of the conditions with the residual variances moderated
  across values (limma's empirical Bayes); two conditions: moderated t, more: moderated F; Benjamini-Hochberg FDR.
log_run() / recent(): the analyses the user ran in the app (request and headline), for the agent's review.
embedding_facts(): Leiden clusters against the labels and, in a joint embedding, against the slices.
"""
import json, time
import numpy as np, pandas as pd, scipy.sparse as sp
from . import config, embedding, fastpar, model, niches as nmod, projects

MAX_TEST_CELLS = 20000
ALL_TEST_CELLS = 5000   # per side when every group is compared (de_all): quicker, the AUC hardly moves
ALL_MIN_CELLS = 30      # smaller groups are skipped by de_all
ALL_MAX_ROWS = 1000000  # de_all's full table beyond this many rows (genes x groups): significant genes only
ALL = ("all", "__all__")  # group "all" (unless a group is called so): every group against the other cells
MIN_CELLS = 10       # smallest group compared, and smallest group per slice for the per-slice check
PB_MIN_CELLS = 30    # DESeq2: smallest pseudobulk sample (cells of one side in one replicate)
PB_MIN_COUNT = 10    # DESeq2: genes with fewer counts over all pseudobulk samples are left out
import re as _re
from .osutil import dump_json
DONOR = _re.compile(r"donor|patient|subject|individual|participant|case|mouse|animal|person", _re.I)
PSEUDO = 0.01        # pseudocount (normalized counts per 100) for the fold change
BUILTIN = ["cell type", "lineage", "region", "niche", "leiden", "selection", "slice"]
NOTE = ("p-values treat cells as independent and are tiny for large groups: rank by AUC and fold change, and with "
        "several slices trust changes seen in every slice. Imaging panels pick up transcripts of neighboring cells "
        "(segmentation spillover): genes of the cell types around a group can appear as its DE genes.")


class Cells:
    """Cells of one annotation or of a project's slices: labels (a DataFrame, one column per grouping, None where a
    slice lacks it), the slice of every cell, and the counts (loaded when needed)."""

    def __init__(self, sessions, labels, slices, pid=None, missing=None):
        self.sessions, self.labels, self.slices, self.pid = sessions, labels, slices, pid
        self.missing = missing or {}
        self.genes = sessions[0].panel
        self._X = self._E = None

    @property
    def X(self):
        if self._X is None:
            self._X = model.lognorm(sp.vstack([s.counts for s in self.sessions]).tocsr()).tocsr()
        return self._X

    @property
    def E(self):
        """expm1 of X's nonzeros (the normalized counts), for the means of _stats."""
        if self._E is None:
            self._E = np.expm1(self.X.data)
        return self._E

    def replicates(self):
        """The columns a pseudobulk replicate can be (DESeq2): slice, and the metadata columns with one value per slice;
        and the default: a patient / donor-like column when there is one (several slices of one patient are one
        sample), else slice."""
        if self.pid is None:
            return [], None
        cols = ["slice"]
        for c in self.labels.columns:
            if c in BUILTIN or c == "selection":
                continue
            per = self.labels.groupby(self.labels["slice"])[c].nunique(dropna=False)
            if (per <= 1).all() and self.labels[c].nunique() > 1:
                cols.append(c)
        donor = next((c for c in cols[1:] if DONOR.search(c)), None)
        return cols, donor or "slice"

    def counts_sum(self, groups):
        """Raw counts summed over each group of cells: groups = [row indices into the pooled table] -> groups x genes."""
        off = np.cumsum([0] + [len(s.cells) for s in self.sessions])
        out = np.zeros((len(groups), len(self.genes)), np.float64)
        for i, s in enumerate(self.sessions):
            X = sp.csr_matrix(s.counts)
            for k, rows in enumerate(groups):
                r = rows[(rows >= off[i]) & (rows < off[i + 1])] - off[i]
                if len(r):
                    out[k] += np.asarray(X[r].sum(0)).ravel()
        return out

    def options(self):
        """{grouping: [values by size]}; metadata columns only when they differ between slices."""
        out = {}
        for c in self.labels.columns:
            vc = self.labels[c].dropna().value_counts()
            if len(vc) and (c in BUILTIN or len(vc) > 1):
                out[c] = [str(v) for v in vc.index]
        return out


def _labels(s, version=None, selection=None):
    """Grouping columns of one annotation; selection = rows of the cells selected in the viewer."""
    n = len(s.cells)
    r = s.result(version)
    out = {"cell type": s.labels(version), "lineage": s.lineage_labels(version),
           "region": np.where(r["region"] == "unknown", "no region", r["region"]).astype(object)}
    from . import projects  # the regions' groupings (Regions → Group): a cell's group is its region's; none outside
    for g, mp in projects.region_groups(s.sample_id).items():
        if g not in out:
            v = pd.Series(r["region"], dtype=object).astype(str).map(mp)  # (vectorized: millions of cells)
            out[g] = v.where(v.notna() & (v != ""), None).to_numpy(dtype=object)
    meta, lab = nmod.current(s)
    out["niche"] = np.array(nmod.niche_names(s), dtype=object)[lab] if meta is not None and len(lab) == n else None
    lei = embedding.leiden(s.sample_id)
    out["leiden"] = (np.char.add("L", (lei + 1).astype(str)).astype(object) if lei is not None and len(lei) == n
                     else None)
    if selection is not None:
        sel = np.zeros(n, bool)
        sel[np.asarray(selection, dtype=int)] = True
        out["selection"] = np.where(sel, "selected", "not selected").astype(object)
    else:
        out["selection"] = None
    return out


def cells(session=None, version=None, pid=None, selection=None, get_session=None, sessions=None):
    """The cell table: one annotation (session, version) or every slice of project pid, all annotated (the annotations
    in sessions, {sample_id: session id}, else the latest ones). selection: rows (one annotation) or {sample_id: rows}
    (project)."""
    from .session import Session
    from . import joint
    get = get_session or Session
    if pid is None:
        S, versions, sel = [session], [version], [selection]
        meta = {}
    else:
        p = projects.load(pid)
        latest = joint.latest_sessions(pid, sessions)
        ids = [x for x in p["samples"] if x in latest]
        if not ids:
            raise ValueError("no slice of this project has an annotation yet")
        if len(ids) < len(p["samples"]):  # all slices means all: a slice left out silently would change the result
            from . import config
            import json
            names = [(json.load(open(config.SAMPLES / x / "meta.json", encoding="utf-8")).get("name") or x)
                     for x in p["samples"] if x not in latest]
            raise ValueError(f"all slices together need every slice annotated; not yet: {', '.join(names)}. Annotate "
                             f"{'them' if len(names) > 1 else 'it'} first, or analyze this slice only")
        S = [get(latest[x]) for x in ids]
        versions = [None] * len(S)
        sel = [(selection or {}).get(x) for x in ids] if selection else [None] * len(S)
        meta = {x: p["slices"].get(x, {}).get("metadata", {}) for x in ids}
        if any(s.panel != S[0].panel for s in S[1:]):
            raise ValueError("the slices of this project do not share one gene panel")
    parts, missing = [], {}
    # pooled slices: Leiden clusters and niches mean the same in every slice only when they come from one joint run of
    # the project (slice A's L3 and slice B's L3 of separate runs are unrelated groups)
    unjoint = set()
    if pid is not None and len(S) > 1:
        emb = [embedding.current(s.sample_id) or {} for s in S]
        if len({(e.get("project"), e.get("run_id")) for e in emb}) > 1 or any(e.get("project") != pid for e in emb):
            unjoint.add("leiden")
        nm = [nmod.current(s)[0] or {} for s in S]
        if len({(m.get("project"), m.get("run_id")) for m in nm}) > 1 or any(m.get("project") != pid for m in nm):
            unjoint.add("niche")
    for s, v, se in zip(S, versions, sel):
        lab = _labels(s, v, se)
        for k in unjoint:
            lab[k] = None
        n = len(s.cells)
        if selection is not None and lab["selection"] is None:
            lab["selection"] = np.full(n, "not selected", dtype=object)
        for k, x in list(lab.items()):
            if x is None:
                if k != "selection":
                    missing.setdefault(k, []).append(s.sample_id)
                lab[k] = np.full(n, None, dtype=object)
        section = p["slices"].get(s.sample_id, {}).get("section") if pid else None  # as renamed in the project
        lab["slice"] = np.full(n, str(section or s.meta.get("project_section", s.sample_id) if pid else s.sample_id),
                               dtype=object)  # (as text: a section numbered 1 is chosen as "1" in the page)
        for mk, mv in meta.get(s.sample_id, {}).items():
            lab[mk if mk not in BUILTIN and mk not in lab else f"meta:{mk}"] = np.full(n, str(mv), dtype=object)  # (not over a region grouping)
        parts.append(pd.DataFrame(lab))
    for k in unjoint:
        missing[k] = ["not from one run of all slices together: run it for the project"]
    labels = pd.concat(parts, ignore_index=True)
    labels = labels[[c for c in labels.columns if labels[c].notna().any()]]
    slices = np.concatenate([np.full(len(s.cells), s.sample_id, dtype=object) for s in S])
    return Cells(S, labels, slices, pid, missing)


def _column(c, by):
    if by not in c.labels.columns:
        m = c.missing.get(by) or []
        why = (f" ({m[0]})" if m and m[0].startswith("not from one run") else f" (missing in {', '.join(m)})") if m else ""
        raise ValueError(f"no grouping {by!r}{why}; groupings: {', '.join(c.labels.columns)}")
    return c.labels[by].to_numpy()


def _mask(c, by, values):
    col = _column(c, by)
    values = [values] if isinstance(values, str) else list(values)
    have = set(col[pd.notna(col)])
    bad = [v for v in values if v not in have]
    if bad:
        raise ValueError(f"no {by} {', '.join(map(repr, bad))}; {by} values: {', '.join(map(str, sorted(have, key=str)[:40]))}")
    return np.isin(col, values)


def _within(c, within):
    """Cells inside every filter of within = {grouping: value or [values]}."""
    m = np.ones(len(c.labels), bool)
    for by, v in (within or {}).items():
        m &= _mask(c, by, v)
    return m


def _wilcoxon(X, a, b, rng, max_cells, parallel=True):
    """AUC and two-sided p of every gene, group rows a vs reference rows b (subsampled). The genes go in blocks of 256,
    on several threads when parallel (each block on its own: the same numbers)."""
    import warnings
    from scipy.stats import mannwhitneyu
    sa = np.sort(rng.choice(a, max_cells, replace=False)) if len(a) > max_cells else a
    sb = np.sort(rng.choice(b, max_cells, replace=False)) if len(b) > max_cells else b
    A, B = X[sa].tocsc(), X[sb].tocsc()

    def block(g):
        with warnings.catch_warnings():  # (a gene detected nowhere: no spread, p NaN, set to 1 below)
            warnings.simplefilter("ignore", RuntimeWarning)
            return mannwhitneyu(A[:, g:g + 256].toarray(), B[:, g:g + 256].toarray(), axis=0, method="asymptotic")
    starts = range(0, X.shape[1], 256)
    res = fastpar.pmap(block, starts) if parallel else [block(g) for g in starts]
    auc, p = np.empty(X.shape[1]), np.empty(X.shape[1])
    for g, (u, pv) in zip(starts, res):
        auc[g:g + 256], p[g:g + 256] = u / (len(sa) * len(sb)), pv
    return np.nan_to_num(auc, nan=0.5), np.nan_to_num(p, nan=1.0), len(sa), len(sb)


def _bh(p):
    n = len(p)
    o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[o] = np.minimum(q, 1)
    return out


def _stats(c, rows):
    """Mean log-normalized expression, mean normalized counts and detection rate of every gene over the rows (ascending)
    of c.X: scipy's X[rows].mean(0) and the rest, summed in the same order without copying the rows (fastpar)."""
    return fastpar.col_means(c.X, rows, c.E)


def _compare(c, a_mask, b_mask, max_cells=MAX_TEST_CELLS, parallel=True):
    """DE of the cells in a_mask vs those in b_mask: the full gene table and the per-slice check. parallel: the
    Wilcoxon test's gene blocks on several threads (de_all runs the groups so instead)."""
    a, b = np.where(a_mask)[0], np.where(b_mask)[0]
    if len(a) < MIN_CELLS or len(b) < MIN_CELLS:
        raise ValueError(f"too few cells to compare: {len(a)} in the group, {len(b)} in the reference "
                         f"(at least {MIN_CELLS} each)")
    X = c.X
    rng = np.random.default_rng(0)
    auc, p, na, nb = _wilcoxon(X, a, b, rng, max_cells, parallel)
    ml_a, me_a, pa = _stats(c, a)
    ml_b, me_b, pb = _stats(c, b)
    lfc = np.log2((me_a + PSEUDO) / (me_b + PSEUDO))
    t = pd.DataFrame({"gene": c.genes, "auc": auc, "log2fc": lfc, "mean_group": ml_a, "mean_reference": ml_b,
                      "pct_group": pa, "pct_reference": pb, "p": p, "fdr": _bh(p)})
    out = {"n_group": int(len(a)), "n_reference": int(len(b)), "n_tested": [int(na), int(nb)]}
    ids = list(dict.fromkeys(c.slices))
    if len(ids) > 1:
        sl = {}
        for x in ids:
            m = c.slices == x
            sl[x] = (np.where(a_mask & m)[0], np.where(b_mask & m)[0])
        out["cells_per_slice"] = {x: [int(len(u)), int(len(v))] for x, (u, v) in sl.items()}
        both = [x for x, (u, v) in sl.items() if len(u) >= MIN_CELLS and len(v) >= MIN_CELLS]
        if len(both) >= 2:
            # the fold change within every slice that has enough cells of both groups
            per = np.array([np.log2((_stats(c, sl[x][0])[1] + PSEUDO) / (_stats(c, sl[x][1])[1] + PSEUDO))
                            for x in both])
            sign = np.sign(lfc)
            t["slices_same_direction"] = [f"{int((np.sign(per[:, g]) == sign[g]).sum())}/{len(both)}"
                                          for g in range(len(lfc))]
            out["per_slice_check"] = f"fold change within each of {len(both)} slices with both groups"
        else:
            # groups made of different slices: is every slice of the group beyond every slice of the reference?
            ga = [x for x, (u, _) in sl.items() if len(u) >= MIN_CELLS]
            gb = [x for x, (_, v) in sl.items() if len(v) >= MIN_CELLS]
            if not ga or not gb:  # one side spread thin (under MIN_CELLS in every slice): nothing to check per slice
                out["per_slice_check"] = (f"none: the {'group' if not ga else 'reference'} has fewer than {MIN_CELLS} "
                                          "cells in every slice")
            elif len(ga) + len(gb) > 2:
                ma = np.array([_stats(c, sl[x][0])[0] for x in ga])
                mb = np.array([_stats(c, sl[x][1])[0] for x in gb])
                up = ma.min(0) > mb.max(0)
                down = ma.max(0) < mb.min(0)
                t["slices_separate"] = np.where(lfc > 0, up, down)
                out["per_slice_check"] = (f"every one of {len(ga)} group slice(s) beyond every one of {len(gb)} "
                                          "reference slice(s) (mean expression)")
            elif ga == gb:  # both groups in the same one slice (e.g. two regions of it)
                out["per_slice_check"] = (f"none: both groups are in one slice ({ga[0]}), so there is no other slice "
                                          "to check the difference in")
            else:
                out["per_slice_check"] = ("none: one slice on each side, so a difference between the groups cannot be "
                                          "told from one between the two slices")
    return t, out


def _rows(t, n):
    keep = [c for c in t.columns if c != "p"]
    r = t[keep].head(n).astype({c: float for c in keep if c not in ("gene", "slices_same_direction", "slices_separate")})
    for col in ("auc", "log2fc", "mean_group", "mean_reference", "pct_group", "pct_reference"):
        r[col] = r[col].round(3)
    r["fdr"] = r["fdr"].map(lambda v: float(f"{v:.2e}"))
    return r.to_dict("records")


def de(c, by, group, reference="rest", within=None, top=20, full=False):
    """DE genes of group (a value of grouping by, or a list of values) against reference ("rest" = every other cell
    with a value of by, or a value / list), among the cells in within ({grouping: value(s)}). group "all": every
    group against the other cells (de_all)."""
    col = _column(c, by)
    if isinstance(group, str) and group in ALL and group not in set(col[pd.notna(col)]):
        return de_all(c, by, within, top, full)
    w = _within(c, within)
    a = _mask(c, by, group) & w
    b = (pd.notna(col) & ~_mask(c, by, group) if reference in (None, "", "rest") else _mask(c, by, reference)) & w
    if (a & b).any():
        raise ValueError("the group and the reference overlap")
    t, info = _compare(c, a, b)
    sig = t.fdr < 0.05
    up = t[sig & (t.auc > 0.5) & (t.log2fc > 0)].sort_values("auc", ascending=False)
    down = t[sig & (t.auc < 0.5) & (t.log2fc < 0)].sort_values("auc")
    res = {"scope": "project" if c.pid else "slice", "by": by, "group": group,
           "reference": "rest" if reference in (None, "", "rest") else reference, "within": within or {},
           **info, "n_up": int(len(up)), "n_down": int(len(down)), "up": _rows(up, top), "down": _rows(down, top),
           "note": NOTE}
    if full:
        res["table"] = t.sort_values("auc", ascending=False)
    return res


def de_all(c, by, within=None, top=20, full=False, min_cells=ALL_MIN_CELLS, max_cells=ALL_TEST_CELLS):
    """DE genes of every group of grouping by against the other cells (within the filter), as de() gives them for
    one group: per group the top up and down genes, counts and the per-slice check. full: the gene table of every
    group (column group first)."""
    w = _within(c, within)
    col = _column(c, by)
    vc = pd.Series(col[w & pd.notna(col)]).value_counts()
    if len(vc) < 2:
        raise ValueError(f"{by} has only one value here ({', '.join(map(str, vc.index))}): nothing to compare")
    groups, tables, skipped = [], [], []
    todo = []
    for g, n in vc.items():
        if n < min_cells or (w & pd.notna(col) & (col != g)).sum() < MIN_CELLS:
            skipped.append(str(g))
        else:
            todo.append(g)
    _ = c.X, c.E  # (made once, before the threads)
    done = fastpar.pmap(lambda g: _compare(c, w & (col == g), w & pd.notna(col) & (col != g), max_cells=max_cells,
                                           parallel=False), todo)  # one group per thread: the same tables
    for g, (t, info) in zip(todo, done):
        sig = t.fdr < 0.05
        up = t[sig & (t.auc > 0.5) & (t.log2fc > 0)].sort_values("auc", ascending=False)
        down = t[sig & (t.auc < 0.5) & (t.log2fc < 0)].sort_values("auc")
        groups.append({"group": str(g), **info, "n_up": int(len(up)), "n_down": int(len(down)),
                       "up": _rows(up, top), "down": _rows(down, top)})
        if full:
            tables.append(t.sort_values("auc", ascending=False).assign(group=str(g)))
    if not groups:
        raise ValueError(f"no {by} group has {min_cells} cells or more here")
    res = {"scope": "project" if c.pid else "slice", "by": by, "group": "all", "reference": "rest",
           "within": within or {}, "all": True, "groups": groups, "skipped_small": skipped,
           "max_test_cells": max_cells, "note": NOTE}
    if full:
        tab = pd.concat(tables, ignore_index=True)
        tab = tab[["group"] + [k for k in tab.columns if k != "group"]]
        if len(tab) > ALL_MAX_ROWS:
            tab = tab[tab.fdr < 0.05]
            res["table_note"] = "the table holds the genes with FDR < 0.05 (all genes of every group would be too many rows)"
        res["table"] = tab
    return res


def deseq2(c, by, group, reference="rest", within=None, replicate=None, top=20, full=False, min_cells=PB_MIN_CELLS):
    """DESeq2 on pseudobulk counts: the cells of group (a value of grouping by, or a list) against reference ("rest"
    or a value / list), among the cells in within; the raw counts are summed per replicate (a column of
    Cells.replicates(): a slice, or a patient / donor) and side. Replicate samples with fewer than min_cells cells are
    left out. Needs at least 2 samples on each side. Paired design (~ replicate + side) when at least 2 replicates hold
    both sides (e.g. two regions in the same slices; then only those), else ~ side (e.g. conditions made of different
    slices or patients)."""
    import logging, warnings
    cols, default = c.replicates()
    if not cols:
        raise ValueError("DESeq2 needs several slices (pseudobulk replicates): use it on a project's slices, or use "
                         "DE genes on one slice")
    replicate = replicate or default
    if replicate not in cols:
        raise ValueError(f"the replicate must be one of {', '.join(cols)} (a column with one value per slice)")
    col = _column(c, by)
    if isinstance(group, str) and group in ALL:
        raise ValueError("DESeq2 compares one group with the rest or with another group, not every group at once")
    w = _within(c, within)
    a = _mask(c, by, group) & w
    b = (pd.notna(col) & ~_mask(c, by, group) if reference in (None, "", "rest") else _mask(c, by, reference)) & w
    if (a & b).any():
        raise ValueError("the group and the reference overlap")
    rep = c.labels[replicate].astype(str).to_numpy()
    samples, groups, small = [], [], []
    for side, m in (("group", a), ("reference", b)):
        for r in sorted(set(rep[m])):
            rows = np.where(m & (rep == r))[0]
            (samples if len(rows) >= min_cells else small).append({"replicate": r, "side": side, "n_cells": int(len(rows))})
            if len(rows) >= min_cells:
                groups.append(rows)
    na = sum(x["side"] == "group" for x in samples)
    nb = len(samples) - na
    if na < 2 or nb < 2:
        raise ValueError(f"DESeq2 needs at least 2 replicates ({replicate}) with {min_cells}+ cells on each side; here "
                         f"{na} for the group and {nb} for the reference" + (f" ({len(small)} left out with fewer cells)" if small else ""))
    shared = {x["replicate"] for x in samples if x["side"] == "group"} & {x["replicate"] for x in samples if x["side"] == "reference"}
    paired = len(shared) >= 2
    if paired:  # only the replicates holding both sides
        keep = [i for i, x in enumerate(samples) if x["replicate"] in shared]
        samples, groups = [samples[i] for i in keep], [groups[i] for i in keep]
    counts = c.counts_sum(groups)
    genes = np.asarray(c.genes)
    ok = counts.sum(0) >= PB_MIN_COUNT
    names = [f"{x['replicate']}|{x['side']}" for x in samples]
    cm = pd.DataFrame(np.rint(counts[:, ok]).astype(np.int64), index=names, columns=genes[ok])
    md = pd.DataFrame({"side": [x["side"] for x in samples], "rep": [x["replicate"] for x in samples]}, index=names)
    design = "~rep + side" if paired else "~side"
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    logging.getLogger("pydeseq2").setLevel(logging.ERROR)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        dds = DeseqDataSet(counts=cm, metadata=md, design=design, quiet=True, n_cpus=1)
        dds.deseq2()
        stt = DeseqStats(dds, contrast=["side", "group", "reference"], quiet=True, n_cpus=1)
        stt.summary()
    r = stt.results_df
    norm = pd.DataFrame(dds.layers["normed_counts"], index=names, columns=cm.columns)
    side = md["side"].to_numpy()
    t = pd.DataFrame({"gene": r.index.to_numpy(), "log2fc": r["log2FoldChange"].to_numpy(), "lfc_se": r["lfcSE"].to_numpy(),
                      "base_mean": r["baseMean"].to_numpy(),
                      "mean_group": norm[side == "group"].mean(0).reindex(r.index).to_numpy(),
                      "mean_reference": norm[side == "reference"].mean(0).reindex(r.index).to_numpy(),
                      "stat": r["stat"].to_numpy(), "p": r["pvalue"].to_numpy(), "padj": r["padj"].to_numpy()})
    t["padj"] = t["padj"].fillna(1.0)
    t["p"] = t["p"].fillna(1.0)
    t = t.sort_values(["padj", "stat"], key=lambda v: v if v.name == "padj" else -v.abs())  # then the strongest
    sig = t.padj < 0.05
    up, down = t[sig & (t.log2fc > 0)], t[sig & (t.log2fc < 0)]
    def rows(x):
        out = x.head(top)[["gene", "log2fc", "lfc_se", "base_mean", "mean_group", "mean_reference", "padj"]].copy()
        for k in ("log2fc", "lfc_se"):
            out[k] = out[k].astype(float).round(3)
        for k in ("base_mean", "mean_group", "mean_reference"):
            out[k] = out[k].astype(float).round(1)
        out["padj"] = out["padj"].map(lambda v: float(f"{v:.2e}"))
        return out.to_dict("records")
    res = {"scope": "project", "by": by, "group": group,
           "reference": "rest" if reference in (None, "", "rest") else reference, "within": within or {},
           "replicate": replicate, "design": design.replace("rep", replicate), "paired": paired,
           "samples": samples, "left_out": small, "n_samples": [na if not paired else len(shared), nb if not paired else len(shared)],
           "n_group": int(sum(x["n_cells"] for x in samples if x["side"] == "group")),
           "n_reference": int(sum(x["n_cells"] for x in samples if x["side"] == "reference")),
           "n_genes_tested": int(ok.sum()), "n_up": int(len(up)), "n_down": int(len(down)), "up": rows(up), "down": rows(down),
           # every gene tested, for the volcano plot: [gene, log2 fold change, adjusted p]
           "volcano": [[str(g), round(float(l), 3), float(f"{q:.3g}")] for g, l, q in
                       zip(t.gene, t.log2fc, t.padj) if np.isfinite(l)],
           "note": ("DESeq2 (PyDESeq2) on pseudobulk counts: the raw counts of each side summed per " + replicate +
                    (", paired by " + replicate if paired else "") + ". The replicates carry the statistics: with 2 per side"
                    " few genes are trustworthy, 3 or more are better. Imaging panels pick up transcripts of neighboring"
                    " cells (segmentation spillover), as in DE genes.")}
    if full:
        res["table"] = t
    return res


def markers(c, by, top=8, within=None, min_cells=30):
    """The top up genes (AUC) of every group of grouping by against the other cells (within the filter)."""
    w = _within(c, within)
    col = _column(c, by)
    vc = pd.Series(col[w & pd.notna(col)]).value_counts()
    if len(vc) < 2:
        raise ValueError(f"{by} has only one value here ({', '.join(map(str, vc.index))}): nothing to compare")
    out = []
    todo = [(g, n) for g, n in vc.items() if n >= min_cells and (w & pd.notna(col) & (col != g)).sum() >= MIN_CELLS]
    _ = c.X, c.E  # (made once, before the threads)
    done = fastpar.pmap(lambda gn: _compare(c, w & (col == gn[0]), w & pd.notna(col) & (col != gn[0]), max_cells=5000,
                                            parallel=False)[0], todo)  # 5000 per side: quicker, the AUC hardly moves
    for (g, n), t in zip(todo, done):
        up = t[(t.fdr < 0.05) & (t.log2fc > 0)].sort_values("auc", ascending=False).head(top)
        out.append({"group": str(g), "n_cells": int(n),
                    "genes": [{"gene": r.gene, "auc": round(float(r.auc), 3), "log2fc": round(float(r.log2fc), 2),
                               "pct_group": round(float(r.pct_group), 2), "pct_rest": round(float(r.pct_reference), 2)}
                              for r in up.itertuples()]})
    return {"scope": "project" if c.pid else "slice", "by": by, "within": within or {}, "groups": out,
            "skipped_small": [str(g) for g, n in vc.items() if n < min_cells], "note": NOTE}


PROP_MIN_CELLS = 20   # replicates with fewer cells (in the subset) are left out of the proportion test


def _trigamma_inv(x):
    """The y with trigamma(y) = x (limma's trigammaInverse: Newton's method)."""
    from scipy.special import polygamma
    x = float(x)
    if x > 1e7:
        return 1 / np.sqrt(x)
    if x < 1e-6:
        return 1 / x
    y = 0.5 + 1 / x
    for _ in range(50):
        tri = float(polygamma(1, y))
        dif = tri * (1 - tri / x) / float(polygamma(2, y))
        y += dif
        if -dif / y < 1e-8:
            break
    return y


def _squeeze(s2, df):
    """limma's empirical Bayes: (moderated variances, prior df d0) from residual variances s2 with df each."""
    from scipy.special import digamma, polygamma
    s2 = np.maximum(np.asarray(s2, float), 1e-12)
    if len(s2) < 3 or df <= 0:
        return s2, 0.0
    z = np.log(s2)
    e = z - digamma(df / 2) + np.log(df / 2)
    emean = e.mean()
    evar = e.var(ddof=1) - float(polygamma(1, df / 2))
    if evar > 0:
        d0 = 2 * _trigamma_inv(evar)
        s0 = np.exp(emean + digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0, s0 = np.inf, np.exp(emean)
    if not np.isfinite(d0):
        return np.full(len(s2), s0), 1e6
    return (d0 * s0 + df * s2) / (d0 + df), float(d0)


def proportions(c, condition, cols, within=None, replicate=None):
    """Whether each value of cols makes up a different share of the cells between the conditions (values of the
    slice-level column condition), the replicates (replicate: slice, or a patient / donor column; default
    Cells.replicates()) being the samples; see the module's note (propeller)."""
    from scipy import stats as st
    reps, default = c.replicates()
    if not reps:
        raise ValueError("the proportion test needs a project's slices (the replicates)")
    replicate = replicate or default
    if replicate not in reps:
        raise ValueError(f"the replicate must be one of {', '.join(reps)}")
    if condition in ("slice", replicate):
        raise ValueError(f"compare conditions made of several {replicate}s (a metadata column such as condition), "
                         f"not the {replicate}s themselves")
    w = _within(c, within)
    cond, k = _column(c, condition), _column(c, cols)
    rep = c.labels[replicate].astype(str).to_numpy()
    m = w & pd.notna(cond) & pd.notna(k)
    per = pd.DataFrame({"rep": rep[m], "cond": cond[m].astype(str)}).groupby("rep")["cond"].nunique()
    if (per > 1).any():
        raise ValueError(f"a {replicate} holds several values of {condition} ({', '.join(per[per > 1].index[:4])}): "
                         f"choose slice as the replicate")
    ct = pd.crosstab(pd.Series(rep[m], name="rep"), pd.Series(k[m], name=cols))
    small = ct.index[ct.sum(1) < PROP_MIN_CELLS].tolist()
    ct = ct.drop(index=small)
    rc = pd.Series(cond[m].astype(str), index=None).groupby(rep[m]).first().reindex(ct.index)
    groups = sorted(rc.unique().tolist())
    nper = rc.value_counts().reindex(groups).astype(int)
    if len(groups) < 2:
        raise ValueError(f"only one value of {condition} among the {replicate}s with {PROP_MIN_CELLS}+ cells")
    if (nper < 2).any():
        raise ValueError(f"the test needs 2 or more {replicate}s per value of {condition}; here "
                         + ", ".join(f"{g}: {n}" for g, n in nper.items()))
    N = ct.sum(1).to_numpy(float)[:, None]
    P = ct.to_numpy(float) / N
    Y = np.arcsin(np.sqrt(P))
    g = rc.to_numpy()
    n, G = len(g), len(groups)
    means = np.stack([Y[g == x].mean(0) for x in groups])        # groups x values
    resid = Y - means[[groups.index(x) for x in g]]
    df = n - G
    s2 = (resid ** 2).sum(0) / df
    s2m, d0 = _squeeze(s2, df)
    dft = df + d0
    shares = {x: P[g == x].mean(0) for x in groups}
    out = []
    for j, val in enumerate(ct.columns):
        row = {"value": str(val), "mean_share": {x: round(float(shares[x][j]), 4) for x in groups},
               "n_cells": int(ct.iloc[:, j].sum())}
        if G == 2:
            a, b = groups
            t = (means[1, j] - means[0, j]) / np.sqrt(s2m[j] * (1 / nper[b] + 1 / nper[a]))
            row["t"] = round(float(t), 2)
            row["p"] = float(2 * st.t.sf(abs(t), dft))
            row["log2_ratio"] = round(float(np.log2((shares[b][j] + 1e-4) / (shares[a][j] + 1e-4))), 2)
        else:
            ssb = sum(nper[x] * (means[i, j] - Y[:, j].mean()) ** 2 for i, x in enumerate(groups))
            f = (ssb / (G - 1)) / s2m[j]
            row["F"] = round(float(f), 2)
            row["p"] = float(st.f.sf(f, G - 1, dft))
        out.append(row)
    fdr = _bh(np.array([r["p"] for r in out]))
    for r, q in zip(out, fdr):
        r["fdr"] = float(f"{q:.3g}")
        r["p"] = float(f"{r['p']:.3g}")
    out.sort(key=lambda r: (r["p"], -abs(r.get("t", r.get("F", 0)))))
    return {"condition": condition, "cols": cols, "replicate": replicate, "groups": groups,
            "n_replicates": {x: int(nper[x]) for x in groups}, "left_out": small, "test": "t" if G == 2 else "F",
            "contrast": f"{groups[1]} vs {groups[0]}" if G == 2 else None, "prior_df": round(d0, 1),
            "values": out, "n_significant": int(sum(r["fdr"] < 0.05 for r in out)),
            "note": (f"Each {replicate} is one sample: its share of each {cols} value, arcsine-square-root "
                     "transformed, then a linear model per value with the variances moderated across values (as limma "
                     "/ propeller). log2_ratio: the mean shares' ratio. With 2-3 replicates per condition only large "
                     "changes reach significance; the cell counts are not the statistics, the replicates are.")}


def crosstab(c, rows, cols, within=None, top_cols=None, replicate=None):
    """Cells of grouping rows split over grouping cols: counts, shares within each row, log2 enrichment over
    independence; slice-level rows (slice or a metadata column) of a project also give every slice's shares."""
    w = _within(c, within)
    r, k = _column(c, rows), _column(c, cols)
    m = w & pd.notna(r) & pd.notna(k)
    if not m.any():
        raise ValueError("no cells with both groupings" + (" in that subset" if within else ""))
    ct = pd.crosstab(pd.Series(r[m], name=rows), pd.Series(k[m], name=cols))
    order = ct.sum(0).sort_values(ascending=False).index
    ct = ct[order]
    tot = ct.to_numpy().sum()
    exp = np.outer(ct.sum(1), ct.sum(0)) / tot
    enr = np.log2((ct.to_numpy() + 0.5) / (exp + 0.5))
    out = []
    for i, name in enumerate(ct.index):
        n = int(ct.iloc[i].sum())
        out.append({"row": str(name), "n_cells": n,
                    "shares": {str(col): round(float(v) / n, 4) for col, v in ct.iloc[i].items() if v},
                    "enrichment": {str(col): round(float(enr[i, j]), 2) for j, col in enumerate(ct.columns) if ct.iloc[i, j]}})
    res = {"scope": "project" if c.pid else "slice", "rows": rows, "cols": cols, "within": within or {},
           "n_cells": int(m.sum()), "keys": [str(x) for x in ct.columns], "table": out,
           "note": "enrichment = log2(observed / expected if the two groupings were independent)"}
    per_slice = rows in ("slice",) or rows in BUILTIN or not c.pid or \
        not (c.labels.groupby(c.labels["slice"])[rows].nunique(dropna=False) <= 1).all()
    if not per_slice:  # a slice-level column (one value per slice, e.g. condition): the slices are replicates;
        # a column that varies inside a slice (a region grouping) has no per-slice bars
        reps = []
        for x in dict.fromkeys(c.slices):
            ms = m & (c.slices == x)
            if not ms.any():
                continue
            vc = pd.Series(k[ms]).value_counts(normalize=True)
            reps.append({"slice": str(c.labels.slice.to_numpy()[ms][0]), rows: str(r[ms][0]), "n_cells": int(ms.sum()),
                         "shares": {str(a): round(float(v), 4) for a, v in vc.items()}})
        res["by_slice"] = reps
        try:  # the replicates allow a test of the shares between the conditions
            res["test"] = proportions(c, rows, cols, within, replicate)
        except ValueError as e:
            res["test_note"] = str(e)
    if c.missing:
        res["missing"] = c.missing
    return res


def dots(c, by, genes, within=None, first=(), max_groups=25, compare=None):
    """Dot plot data: for every group of grouping by (within the filter; the groups in first, then the largest, at
    most max_groups, or all of first) the share of its cells detecting each gene and the mean log-normalized expression. compare =
    {"group", "reference"} (values or lists, reference "rest"): only the two sides of that comparison, the rest of
    the cells as one row "all other cells"."""
    col = _column(c, by)
    m = _within(c, within) & pd.notna(col)
    vc = pd.Series(col[m]).value_counts()
    pos = {g: i for i, g in enumerate(c.genes)}
    missing = [g for g in genes if g not in pos]
    if missing:
        raise ValueError(f"not in the panel: {', '.join(missing[:5])}")
    if compare:
        name = lambda v: v if isinstance(v, str) else " + ".join(map(str, v))
        a = m & _mask(c, by, compare["group"])
        ref = compare.get("reference") or "rest"
        b = m & ~a if ref == "rest" else m & _mask(c, by, ref)
        order = [name(compare["group"]), "all other cells" if ref == "rest" else name(ref)]
        sides = [np.where(a)[0], np.where(b)[0]]
        rows = np.concatenate(sides)
        G = sp.csr_matrix((np.ones(len(rows), np.float32), (np.repeat([0, 1], [len(x) for x in sides]), rows)),
                          shape=(2, len(col)))
    else:
        order = [g for g in first if g in vc.index] + [g for g in vc.index if g not in first]
        order = order[:max(max_groups, len(first))]  # every group of a result (first) has its row
        gi = {g: i for i, g in enumerate(order)}
        rows = np.where(m & np.isin(col, order))[0]
        G = sp.csr_matrix((np.ones(len(rows), np.float32), ([gi[x] for x in col[rows]], rows)), shape=(len(order), len(col)))
    X = c.X[:, [pos[g] for g in genes]].tocsc()
    n = np.asarray(G.sum(1)).ravel()
    det = X.copy()
    det.data = (det.data > 0).astype(np.float32)
    pct = np.asarray((G @ det).todense()) / n[:, None]
    mean = np.asarray((G @ X).todense()) / n[:, None]
    return {"by": by, "groups": [str(g) for g in order], "n_cells": [int(x) for x in n], "genes": list(genes),
            "pct": np.round(pct, 3).tolist(), "mean": np.round(mean, 3).tolist(),
            "n_groups_total": len(order) if compare else int(len(vc))}


# ---------------- the analyses the user ran in the app ----------------
def headline(kind, r):
    """A few lines of a result, for the log the agent reads."""
    if kind == "de" and r.get("all"):
        return {"groups": {g["group"]: {"up": [x["gene"] for x in g["up"][:5]], "down": [x["gene"] for x in g["down"][:3]]}
                           for g in r["groups"][:30]}}
    if kind == "de":
        return {"n_group": r["n_group"], "n_reference": r["n_reference"], "n_up": r["n_up"], "n_down": r["n_down"],
                "top_up": [x["gene"] for x in r["up"][:8]], "top_down": [x["gene"] for x in r["down"][:8]],
                "per_slice_check": r.get("per_slice_check")}
    if kind == "markers":
        return {"groups": {g["group"]: [x["gene"] for x in g["genes"][:5]] for g in r["groups"][:20]}}
    if kind == "deseq2":
        return {"replicate": r["replicate"], "design": r["design"], "n_samples": r["n_samples"], "n_up": r["n_up"],
                "n_down": r["n_down"], "top_up": [x["gene"] for x in r["up"][:8]], "top_down": [x["gene"] for x in r["down"][:8]]}
    if kind == "composition":
        return {"rows": [x["row"] for x in r["table"][:20]], "keys": r["keys"][:20], "n_cells": r["n_cells"]}
    if kind == "run":
        return {"n_pairs_panel": r["n_pairs_panel"], "n_significant": r["n_significant"],
                "top": [f'{x["interaction"]}: {x["sender"]} -> {x["receiver"]}' for x in r["interactions"][:8]]}
    if kind == "compare":
        return {"values": r["values"], "n_significant": r["n_significant"],
                "top": [f'{x["interaction"]}: {x["sender"]} -> {x["receiver"]}' for x in r["interactions"][:8]]}
    if kind == "coloc":
        return {"by": r["by"], "radius_um": r["radius_um"],
                "together": [f'{x["group_a"]} + {x["group_b"]} (z {x["z"]})' for x in r["together"][:8]],
                "apart": [f'{x["group_a"]} / {x["group_b"]} (z {x["z"]})' for x in r["apart"][:8]]}
    if kind == "trajectory":
        return {"method": r["method"], "n_on_trajectory": r["n_on_trajectory"],
                "rising": [x["gene"] for x in r["rising"][:8]], "falling": [x["gene"] for x in r["falling"][:8]],
                "transient": [x["gene"] for x in r["transient"][:6]]}
    return {}


def _jdefault(o):
    if isinstance(o, np.generic):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def keep_result(sessions, kind, request, result, table, selection=None, tables=None):
    """The last result of an analysis kind, kept in every given annotation's folder (results/<kind>.json and its full
    table as results/<id>.parquet; tables: more {id: DataFrame}, e.g. a composition's test): the app's Last result
    shows it after a reload, for a project's analysis on each of its slices, and for one the agent ran too."""
    import sys
    tables = {result["id"]: table, **(tables or {})} if table is not None and result.get("id") else dict(tables or {})
    req = {k: v for k, v in request.items() if k not in ("session", "selection")}
    for s in dict.fromkeys(sessions):
        d = s.dir / "results"
        d.mkdir(exist_ok=True)
        f = d / f"{kind}.json"
        try:  # the tables of the result it replaces
            old = json.load(open(f))["result"] if f.exists() else {}
            for oid in (old.get("id"), (old.get("test") or {}).get("id")):
                if oid and oid not in tables and (d / f"{oid}.parquet").exists():
                    (d / f"{oid}.parquet").unlink()
        except (ValueError, KeyError, OSError, AttributeError):
            pass
        for tid, t in tables.items():
            try:
                t.to_parquet(d / f"{tid}.parquet")
            except Exception as e:  # the result stays; only its CSV is missing after a restart
                print(f"could not keep the table of {kind}: {e}", file=sys.stderr)
        dump_json({"time": time.strftime("%Y-%m-%dT%H:%M:%S"), "request": req, "selection": selection,
                   "result": result}, f, default=_jdefault)


def results_signature(s):
    """What changes when a result is kept (the app reloads Last result)."""
    d = s.dir / "results"
    if not d.exists():
        return []
    return sorted([f.stem, int(f.stat().st_mtime)] for f in d.glob("*.json"))


def log_run(s, request, result):
    """Note an analysis the user ran in the app (in the annotation they had open), keeping the last 30."""
    f = s.dir / "analyses.json"
    try:
        log = json.load(open(f)) if f.exists() else []
    except ValueError:
        log = []
    req = {k: v for k, v in request.items() if k not in ("session", "selection", "version")}
    if request.get("selection"):
        req["selection_cells"] = sum(len(v) for v in request["selection"].values())
    kind = request.get("kind")
    log.append({"time": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "section": "signaling" if kind in ("run", "compare") else "trajectory" if kind == "trajectory"
                else "niches" if kind == "coloc" else "compare",
                "request": req, "headline": headline(kind, result)})
    dump_json(log[-30:], f, indent=1)


def recent(s, section=None):
    f = s.dir / "analyses.json"
    log = json.load(open(f)) if f.exists() else []
    return [x for x in log if section is None or x["section"] == section]


# ---------------- embedding ----------------
def embedding_facts(c, meta, top=3):
    """Leiden clusters against the cell types (purity, types split over several clusters) and, with several slices,
    against the slices (clusters made of one slice)."""
    lab, typ = c.labels["leiden"].to_numpy(), c.labels["cell type"].to_numpy()
    ok = pd.notna(lab)
    ct = pd.crosstab(pd.Series(lab[ok], name="leiden"), pd.Series(typ[ok], name="type"))
    clusters = []
    for k in ct.sum(1).sort_values(ascending=False).index:
        row = ct.loc[k]
        n = int(row.sum())
        sh = (row / n).sort_values(ascending=False)
        clusters.append({"cluster": k, "n_cells": n, "purity": round(float(sh.iloc[0]), 2),
                         "types": {str(a): round(float(b), 2) for a, b in sh.head(top).items() if b > 0}})
    types = []
    for t_ in ct.sum(0).sort_values(ascending=False).index:
        col = ct[t_]
        sh = (col / col.sum()).sort_values(ascending=False)
        types.append({"type": str(t_), "n_cells": int(col.sum()),
                      "clusters": {str(a): round(float(b), 2) for a, b in sh.items() if b >= 0.1}})
    out = {"params": meta, "n_cells": int(ok.sum()), "n_clusters": len(clusters),
           "median_purity": round(float(np.median([x["purity"] for x in clusters])), 2) if clusters else None,
           "mixed_clusters": [x["cluster"] for x in clusters if x["purity"] < 0.5],
           "split_types": [x["type"] for x in types if len(x["clusters"]) >= 3],
           "clusters": clusters, "types": types}
    sl = c.labels["slice"].to_numpy()
    if len(set(sl[ok])) > 1:
        overall = pd.Series(sl[ok]).value_counts(normalize=True)
        xs = pd.crosstab(pd.Series(lab[ok], name="leiden"), pd.Series(sl[ok], name="slice"))
        one = []
        for k, row in xs.iterrows():
            sh = row / row.sum()
            top_s = sh.idxmax()
            if sh[top_s] >= 0.8 and overall[top_s] < 0.6 and row.sum() >= 50:
                one.append({"cluster": k, "slice": str(top_s), "share": round(float(sh[top_s]), 2), "n_cells": int(row.sum())})
        out["slices"] = {"overall_shares": {str(a): round(float(b), 3) for a, b in overall.items()},
                         "single_slice_clusters": one}
    return out
