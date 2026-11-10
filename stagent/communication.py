"""Cell-cell communication: which ligand-receptor pairs link groups of cells that sit close enough to signal.

Resource (data/lr_<species>.csv.gz, built by stagent/lr_resource.py): CellChatDB v2, CellPhoneDB v5 (human)
and the OmniPath consensus LIANA uses (CellPhoneDB, CellChatDB, ICELLNET, connectomeDB2020, CellTalkDB), merged; each
pair has its ligand and receptor genes (complexes have several subunits), a pathway and the kind of signaling. A pair
is used when every gene of both sides is in the panel; a complex's expression in a cell is its lowest subunit.
Small-molecule signals (Non-protein Signaling: a synthesizing enzyme stands in for the molecule) are left out unless
asked for.

Scoring, at single-cell resolution (actual neighbors are paired, not whole groups):
1. Neighbors. Cell-Cell Contact: every other cell within contact_um. Secreted, ECM and unannotated signals: the cells
   within secreted_um (the MAX_NB nearest), weighted by a Gaussian of the distance (sigma = secreted_um / 2). Pairs
   never cross slices.
2. Strength of a neighbor pair (sender i, receiver j): L_i * R_j on log-normalized expression. The score of sender
   group A -> receiver group B is its weighted mean over the neighbor pairs from A to B.
3. Significance against random group labels (the CellPhoneDB idea): the expected score is the mean over all neighbor
   pairs, z uses its spread and the number of distinct cells involved (pairs sharing a cell are not independent);
   FDR over every tested (pair, A, B). Tested only when the ligand is detected in min_pct of A and the receptor in
   min_pct of B, and A and B share at least MIN_PAIRS neighbor pairs.
specificity = score / expected score. colocalization = score / (mean L of A * mean R of B): above 1 when the A cells
that express the ligand sit next to the B cells that express the receptor more than the group averages imply.
"""
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import norm
from . import config

CONTACT_UM, SECRETED_UM = 15.0, 50.0
MAX_NB = {"contact": 12, "secreted": 30}
MIN_PAIRS = 20
RESOURCES = ("all", "CellChatDB", "CellPhoneDB", "consensus")
NOTE = ("Co-expression of a ligand and its receptor in neighboring cells suggests, not proves, signaling (protein "
        "levels, processing and activity are not measured). Small panels hold few complete pairs. Transcripts of "
        "neighbors picked up by segmentation (spillover) can make a cell look like it expresses a ligand or receptor.")
_graphs: dict = {}


def species(genes):
    """mouse when most panel genes are written Xxxx, else human."""
    g = [x for x in genes if x[:1].isalpha()]
    return "mouse" if sum(x[1:] != x[1:].upper() for x in g) > len(g) / 2 else "human"


def load_resource(genes, resource="all", nonprotein=False):
    """Pairs whose genes are all in the panel: DataFrame (interaction, ligand, receptor, pathway, category, sources,
    lig_idx, rec_idx) and the number of pairs in the resource. resource = one of RESOURCES or a CSV path with columns
    ligand, receptor (genes of a complex joined by _ or +) and optional pathway, category."""
    if resource in RESOURCES:
        sp = species(genes)
        t = pd.read_csv(config.PACKAGE / "data" / f"lr_{sp}.csv.gz", keep_default_na=False)
        if resource != "all":
            t = t[t.sources.str.split(";").map(lambda s: resource in s)]
    else:
        t = pd.read_csv(resource, keep_default_na=False)
        t.columns = [c.lower() for c in t.columns]
        if not {"ligand", "receptor"} <= set(t.columns):
            raise ValueError("a ligand-receptor CSV needs the columns ligand and receptor")
        for c in ("pathway", "category"):
            if c not in t.columns:
                t[c] = ""
            t[c] = t[c].astype(str).replace({"nan": "", "None": ""})  # a missing value written out as text: no pathway
        t["ligand"] = t.ligand.str.replace("+", "_", regex=False)
        t["receptor"] = t.receptor.str.replace("+", "_", regex=False)
        t["interaction"] = t.ligand.str.replace("_", "+") + " → " + t.receptor.str.replace("_", "+")
        t["sources"] = "custom"
    if not nonprotein:
        t = t[t.category != "Non-protein Signaling"]
    n_all = len(t)
    idx = {g: i for i, g in enumerate(genes)}
    side = lambda s: [idx.get(g) for g in s.split("_")]
    t = t.assign(lig_idx=t.ligand.map(side), rec_idx=t.receptor.map(side))
    ok = t.lig_idx.map(lambda x: None not in x) & t.rec_idx.map(lambda x: None not in x)
    return t[ok].reset_index(drop=True), n_all


def _slice_edges(s, kind, radius):
    """Neighbor pairs of one annotation's cells: (sender rows, receiver rows, weights), cached per sample."""
    key = (s.sample_id, s._cells_mtime, kind, float(radius))
    if key not in _graphs:
        xy = s.cells[["x", "y"]].to_numpy()
        k = MAX_NB[kind]
        d, nb = cKDTree(xy).query(xy, k=k + 1, distance_upper_bound=float(radius))
        d, nb = d[:, 1:], nb[:, 1:]                       # column 0 is the cell itself
        ok = np.isfinite(d)
        dst = np.repeat(np.arange(len(xy), dtype=np.int32), k).reshape(-1, k)[ok]
        src = nb[ok].astype(np.int32)
        w = (np.ones(len(src), np.float32) if kind == "contact" else
             np.exp(-d[ok] ** 2 / (2 * (radius / 2) ** 2)).astype(np.float32))
        while len(_graphs) >= 6:
            _graphs.pop(next(iter(_graphs)))
        _graphs[key] = (src, dst, w)
    return _graphs[key]


def _edges(c, kind, radius):
    """Neighbor pairs of all cells of the table (rows offset slice by slice)."""
    parts, off = [], 0
    for s in c.sessions:
        src, dst, w = _slice_edges(s, kind, radius)
        parts.append((src + off, dst + off, w))
        off += len(s.cells)
    return tuple(np.concatenate(x) for x in zip(*parts))


def _expr(Xc, idx, cache):
    """Expression of one side (lowest subunit of a complex) in every cell."""
    key = tuple(idx)
    if key not in cache:
        cols = [np.asarray(Xc[:, i].todense()).ravel() for i in idx]
        cache[key] = cols[0] if len(cols) == 1 else np.minimum.reduce(cols)
    return cache[key]


class _Graph:
    """The neighbor pairs between groups for one kind of signaling."""

    def __init__(self, c, kind, radius, grp, keep, K):
        src, dst, w = _edges(c, kind, radius)
        m = keep[src] & keep[dst]
        self.src, self.dst, self.w = src[m], dst[m], w[m]
        self.pair = grp[self.src] * K + grp[self.dst]
        self.wsum = np.bincount(self.pair, self.w, K * K)
        self.n = np.bincount(self.pair, minlength=K * K)
        # distinct senders and receivers per group pair: pairs that share a cell are not independent
        us = np.unique(self.pair.astype(np.int64) * len(grp) + self.src) // len(grp)
        ud = np.unique(self.pair.astype(np.int64) * len(grp) + self.dst) // len(grp)
        self.n_eff = np.minimum(self.n, np.minimum(np.bincount(us, minlength=K * K), np.bincount(ud, minlength=K * K)))
        self.W = self.w.sum()

    def score(self, L, R):
        v = L[self.src] * R[self.dst]
        wv = self.w * v
        S = np.bincount(self.pair, wv, len(self.wsum)) / np.maximum(self.wsum, 1e-12)
        m = wv.sum() / self.W
        sd = np.sqrt(max((self.w * (v - m) ** 2).sum() / self.W, 1e-12))
        z = (S - m) / (sd / np.sqrt(np.maximum(self.n_eff, 1)))
        return S, m, z


def run(c, by="cell type", within=None, resource="all", contact_um=CONTACT_UM, secreted_um=SECRETED_UM,
        min_pct=0.1, nonprotein=False, top=50, progress=None):
    """Interactions between the groups of grouping by (cells inside within = {grouping: value(s)})."""
    from . import analysis
    say = progress or (lambda f, m: None)
    col = analysis._column(c, by)
    keep = analysis._within(c, within) & pd.notna(col)
    groups = list(pd.Series(col[keep]).value_counts().index)
    if len(groups) < 1:
        raise ValueError("no cells to analyze")
    K = len(groups)
    gi = {g: i for i, g in enumerate(groups)}
    grp = np.array([gi.get(x, 0) for x in col], np.int64)
    size = np.bincount(grp[keep], minlength=K).astype(float)
    lr, n_all = load_resource(c.genes, resource, nonprotein)
    say(0.05, f"{len(lr)} of {n_all} ligand-receptor pairs have all their genes in the panel")
    if not len(lr):
        raise ValueError(f"none of the {n_all} ligand-receptor pairs has all its genes in this panel")
    graphs = {}
    for kind, radius in (("contact", contact_um), ("secreted", secreted_um)):
        say(0.08 if kind == "contact" else 0.15, f"neighbors within {radius:g} µm ({kind})")
        graphs[kind] = _Graph(c, kind, radius, grp, keep, K)
    Xc = c.X.tocsc()
    cache, rows = {}, []
    kf = keep.astype(float)
    A, B = np.divmod(np.arange(K * K), K)
    for n, r in enumerate(lr.itertuples()):
        if n % 25 == 0:
            say(0.2 + 0.75 * n / len(lr), f"ligand-receptor pair {n + 1} of {len(lr)}")
        L, R = _expr(Xc, r.lig_idx, cache), _expr(Xc, r.rec_idx, cache)
        pl = np.bincount(grp, (L > 0) * kf, K) / np.maximum(size, 1)
        pr = np.bincount(grp, (R > 0) * kf, K) / np.maximum(size, 1)
        ml = np.bincount(grp, L * kf, K) / np.maximum(size, 1)
        mr = np.bincount(grp, R * kf, K) / np.maximum(size, 1)
        g = graphs["contact" if r.category == "Cell-Cell Contact" else "secreted"]
        test = (pl[A] >= min_pct) & (pr[B] >= min_pct) & (g.n >= MIN_PAIRS)
        if not test.any():
            continue
        S, m, z = g.score(L, R)
        for k in np.where(test)[0]:
            a, b = A[k], B[k]
            rows.append((n, a, b, S[k], S[k] / max(m, 1e-12), S[k] / max(ml[a] * mr[b], 1e-12), pl[a], pr[b], int(g.n[k]),
                         float(norm.sf(z[k]))))
    if not rows:
        raise ValueError("no ligand-receptor pair passes the expression filter (ligand and receptor each detected in "
                         f"at least {min_pct:.0%} of a group); lower min_pct or use larger groups")
    t = pd.DataFrame(rows, columns=["i", "a", "b", "score", "specificity", "colocalization", "pct_ligand",
                                    "pct_receptor", "neighbor_pairs", "p"])
    from .analysis import _bh
    t["fdr"] = _bh(t.p.to_numpy())
    info = lr.iloc[t.i.to_numpy()]
    t = pd.concat([info[["interaction", "ligand", "receptor", "pathway", "category"]].reset_index(drop=True),
                   pd.DataFrame({"sender": np.array(groups, dtype=object)[t.a], "receiver": np.array(groups, dtype=object)[t.b]}),
                   t.drop(columns=["i", "a", "b"])], axis=1)
    t["significant"] = (t.fdr < 0.05) & (t.specificity > 1)
    t = t.sort_values(["significant", "score"], ascending=False).reset_index(drop=True)
    sig = t[t.significant]
    net_n = pd.crosstab(sig.sender, sig.receiver).reindex(index=groups, columns=groups, fill_value=0)
    net_s = sig.pivot_table(index="sender", columns="receiver", values="score", aggfunc="sum").reindex(
        index=groups, columns=groups).fillna(0)
    pw = (sig.assign(pathway=sig.pathway.replace("", "other")).groupby("pathway")
          .agg(score=("score", "sum"), n=("score", "size")).sort_values("score", ascending=False))
    pathways = []
    for p_, x in pw.head(25).iterrows():
        tops = sig[sig.pathway.replace("", "other") == p_].groupby(["sender", "receiver"]).score.sum().sort_values(ascending=False)
        pathways.append({"pathway": p_, "score": round(float(x.score), 3), "n": int(x.n),
                         "top": [f"{a} → {b}" for a, b in tops.head(3).index]})
    say(1.0, "done")
    return {"scope": "project" if c.pid else "slice", "by": by, "within": within or {}, "resource": resource,
            "species": species(c.genes), "contact_um": contact_um, "secreted_um": secreted_um, "min_pct": min_pct,
            "n_pairs_resource": int(n_all), "n_pairs_panel": int(len(lr)), "n_tested": int(len(t)),
            "n_significant": int(len(sig)), "groups": groups, "group_sizes": [int(x) for x in size],
            "interactions": _records(sig.head(top)),
            "network": {"count": net_n.to_numpy().tolist(), "score": np.round(net_s.to_numpy(), 3).tolist()},
            "pathways": pathways, "note": NOTE, "table": t}


def bubble(t, n=20, max_pairs=30):
    """A bubble plot of a run's table: the n interactions with the highest significant score (rows) over the
    sender -> receiver pairs where they are strongest (columns: each row's strongest pair, then by summed significant
    score, up to max_pairs), with every tested score in that grid (significant or not)."""
    sig = t[t.significant]
    rows = list(sig.groupby("interaction").score.max().sort_values(ascending=False).head(n).index)
    top = sig[sig.interaction.isin(rows)]
    best = top.sort_values("score", ascending=False).drop_duplicates("interaction")  # every row's strongest pair first
    pairs = list(dict.fromkeys([*zip(best.sender, best.receiver),
                                *top.groupby(["sender", "receiver"]).score.sum().sort_values(ascending=False).index]))
    pairs = sorted(pairs[:max(max_pairs, len(set(zip(best.sender, best.receiver))))], key=lambda p_: -top[
        (top.sender == p_[0]) & (top.receiver == p_[1])].score.sum())
    ri, ci = {r: i for i, r in enumerate(rows)}, {p_: j for j, p_ in enumerate(pairs)}
    x = t[t.interaction.isin(rows)]
    cells = [[ri[r.interaction], ci[(r.sender, r.receiver)], round(float(r.score), 4), float(f"{r.fdr:.2e}"), bool(r.significant)]
             for r in x.itertuples() if (r.sender, r.receiver) in ci]
    pw = sig.drop_duplicates("interaction").set_index("interaction").pathway
    return {"rows": rows, "pathways": [pw.get(r, "") for r in rows], "cols": [list(p_) for p_ in pairs], "cells": cells,
            "n_interactions": int(sig.interaction.nunique()),
            "n_pairs": int(len(sig.groupby(["sender", "receiver"])))}


def _records(t):
    out = []
    for r in t.itertuples():
        out.append({"interaction": r.interaction, "pathway": r.pathway, "category": r.category, "sender": r.sender,
                    "receiver": r.receiver, "score": round(float(r.score), 3), "specificity": round(float(r.specificity), 2),
                    "colocalization": round(float(r.colocalization), 2), "pct_ligand": round(float(r.pct_ligand), 2),
                    "pct_receptor": round(float(r.pct_receptor), 2), "neighbor_pairs": int(r.neighbor_pairs),
                    "fdr": float(f"{r.fdr:.2e}")})
    return out


def compare(c, split, by="cell type", within=None, top=40, progress=None, **kw):
    """The same analysis in every group of split (e.g. condition, slice, region): which interactions are stronger
    where. Interactions significant in at least one group, with the score in each."""
    from . import analysis
    col = analysis._column(c, split)
    w = analysis._within(c, within)
    values = list(pd.Series(col[w & pd.notna(col)]).value_counts().index)
    if len(values) < 2:
        raise ValueError(f"{split} has only one value here: nothing to compare")
    say = progress or (lambda f, m: None)
    res = {}
    for i, v in enumerate(values):
        part = lambda f, m, i=i, v=v: say((i + f) / len(values), f"{split} {v}: {m}")
        try:
            res[v] = run(c, by, {**(within or {}), split: v}, top=0, progress=part, **kw)["table"]
        except ValueError:
            continue
    if len(res) < 2:
        raise ValueError("fewer than two groups have interactions to compare")
    key = ["interaction", "sender", "receiver"]
    sig = pd.concat([t[t.significant][key] for t in res.values()]).drop_duplicates()
    tab = sig.copy()
    for v, t in res.items():
        x = t.set_index(key)
        tab[f"score:{v}"] = [float(x.score.get(k, 0.0)) if k in x.index else 0.0 for k in map(tuple, sig[key].to_numpy())]
        tab[f"significant:{v}"] = [bool(x.significant.get(k, False)) if k in x.index else False for k in map(tuple, sig[key].to_numpy())]
    vs = list(res)
    sc = tab[[f"score:{v}" for v in vs]].to_numpy()
    tab["spread"] = sc.max(1) - sc.min(1)
    if len(vs) == 2:
        tab[f"log2({vs[0]}/{vs[1]})"] = np.log2((sc[:, 0] + 0.01) / (sc[:, 1] + 0.01))
    tab = tab.sort_values("spread", ascending=False).reset_index(drop=True)
    recs = []
    for r in tab.head(top).to_dict("records"):
        recs.append({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    note = NOTE
    if c.pid and split not in ("slice",) and len(set(c.labels.slice)) <= len(vs):
        note = "One slice per group: a difference between groups cannot be told from one between slices. " + note
    return {"scope": "project" if c.pid else "slice", "split": split, "by": by, "within": within or {}, "values": vs,
            "n_significant": {v: int(t.significant.sum()) for v, t in res.items()}, "interactions": recs,
            "note": note, "table": tab}


BINS = (("none", "none", (200, 200, 200), 30), ("low", "lowest half", (254, 217, 118), 220),
        ("mid", "50–80th percentile", (253, 141, 60), 235), ("high", "80–95th percentile", (227, 26, 28), 245),
        ("top", "top 5%", (128, 0, 38), 250))


def bins(a):
    """Bin (index into BINS) of every cell's signal: none, then percentiles of the cells with a signal."""
    b = np.zeros(len(a), np.int64)
    pos = a > 0
    if pos.any():
        b[pos] = 1 + np.searchsorted(np.quantile(a[pos], [0.5, 0.8, 0.95]), a[pos], side="right")
    return b


ROLE_BINS = (("none", "no signal", (200, 200, 200), 30), ("weak", "weak signal", (225, 215, 160), 140),
             ("sender", "sends", (33, 102, 172), 240), ("receiver", "receives", (214, 96, 77), 240),
             ("both", "sends and receives", (118, 42, 131), 250))


def role_bins(send, recv):
    """Role (index into ROLE_BINS) of every cell in one signal: sends or receives strongly (above the median of the
    cells with that signal), both, only weakly, or not at all."""
    strong = lambda a: (a > np.median(a[a > 0])) if (a > 0).any() else np.zeros(len(a), bool)
    s, r = strong(send), strong(recv)
    b = np.where((send > 0) | (recv > 0), 1, 0)
    b[s & ~r], b[r & ~s], b[s & r] = 2, 3, 4
    return b


def activity(s, ligand, receptor, side="receiver", category=None, contact_um=CONTACT_UM, secreted_um=SECRETED_UM):
    """Per cell of one annotation: the signal it receives (its receptor times the ligand of its neighbors, weighted)
    or sends (its ligand times the receptor of its neighbors); side="both": (sends, receives)."""
    from . import model
    idx = {g: i for i, g in enumerate(s.panel)}
    li = [idx.get(g) for g in ligand.split("_")]
    ri = [idx.get(g) for g in receptor.split("_")]
    if None in li or None in ri:
        raise ValueError(f"{ligand} or {receptor} is not in this panel")
    if category is None:
        t, _ = load_resource(s.panel, "all", True)
        hit = t[(t.ligand == ligand) & (t.receptor == receptor)]
        category = hit.category.iloc[0] if len(hit) else ""
    kind, radius = ("contact", contact_um) if category == "Cell-Cell Contact" else ("secreted", secreted_um)
    src, dst, w = _slice_edges(s, kind, radius)
    X = model.lognorm(s.counts).tocsc()
    cache = {}
    L, R = _expr(X, li, cache), _expr(X, ri, cache)
    n = len(s.cells)
    wsum = np.maximum(np.bincount(dst, w, n), 1e-12)
    recv = R * np.bincount(dst, w * L[src], n) / wsum   # mean ligand around each cell, times its receptor
    send = L * np.bincount(dst, w * R[src], n) / wsum   # its ligand times the receptor around it
    return (send, recv) if side == "both" else recv if side == "receiver" else send
