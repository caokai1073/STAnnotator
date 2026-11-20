"""Gene program scores: how strongly each cell expresses a set of genes (a program such as proliferation, the
interferon response or hypoxia) beyond what genes of the same expression level give, as scanpy's score_genes (Satija
et al. 2015): the panel's genes are put into 25 bins by mean log-normalized expression; for every program gene up to
50 genes of its bin (not in the program) are controls; a cell's score = its mean over the program genes minus its mean
over the controls. Above 0: the program is higher than background in that cell. Genes are matched case-insensitively
(human and mouse symbols); a program needs MIN_GENES of its genes in the panel.

The built-in programs are short, common marker sets (not a database); the user's own list works the same way. The
scores of an annotation are kept in its folder (programs/<key>.npy, programs.json) so the map can be colored by them
(color mode "program:<key>").
"""
import json, re, time
import numpy as np, pandas as pd
from .osutil import dump_json, replace_atomic

MIN_GENES = 3
N_BINS, CTRL = 25, 50
BUILTIN = {
    "Proliferation": "MKI67 TOP2A CENPF BIRC5 CCNB1 CCNB2 CDK1 UBE2C PCNA MCM2 MCM5 MCM6 TYMS RRM2 TK1 CDC20 AURKA "
                     "AURKB PLK1 NUSAP1 TPX2 HMGB2 SMC4 KIF11",
    "Interferon type I response": "ISG15 IFI6 IFI44L IFIT1 IFIT2 IFIT3 MX1 MX2 OAS1 OAS2 OAS3 RSAD2 IFI27 XAF1 STAT1 "
                                  "IRF7 HERC5 USP18 BST2 LY6E",
    "Interferon gamma response": "CXCL9 CXCL10 CXCL11 IDO1 GBP1 GBP2 GBP4 GBP5 STAT1 IRF1 CIITA CD274 WARS1 SOCS1 "
                                 "PSMB9 TAP1",
    "Hypoxia": "VEGFA SLC2A1 CA9 PGK1 LDHA BNIP3 ENO1 ALDOA P4HA1 ADM NDRG1 ANKRD37 EGLN3 PDK1 HK2 ERO1A",
    "Inflammation (TNF / NF-kB)": "TNF IL1B IL6 CXCL8 CXCL1 CXCL2 CXCL3 CCL2 CCL20 NFKBIA NFKB1 RELB TNFAIP3 ICAM1 "
                                  "PTGS2 IER3 SOD2 IL1A",
    "Antigen presentation (MHC-II)": "HLA-DRA HLA-DRB1 HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 HLA-DMA HLA-DMB CD74 "
                                     "CIITA H2-AA H2-AB1 H2-EB1 H2-DMA H2-DMB1",
    "Cytotoxicity": "GZMA GZMB GZMH GZMK PRF1 NKG7 GNLY KLRD1 KLRK1 CTSW FASLG IFNG",
    "T cell exhaustion": "PDCD1 HAVCR2 LAG3 TIGIT CTLA4 TOX ENTPD1 CXCL13 LAYN BATF",
    "EMT / mesenchymal": "VIM FN1 SNAI2 TWIST1 ZEB1 ZEB2 CDH2 SPARC TAGLN MMP2 POSTN TGFBI SERPINE1 ITGA5 COL1A1 "
                         "COL3A1",
    "Fibrosis / ECM": "COL1A1 COL1A2 COL3A1 COL5A1 COL6A3 FN1 POSTN SPARC LUM DCN FAP TNC THBS2 CTHRC1",
    "Endothelial activation": "SELE SELP VCAM1 ICAM1 ACKR1 PLVAP ANGPT2 ESM1",
    "Stress / immediate early": "FOS FOSB JUN JUNB EGR1 ATF3 IER2 DUSP1 HSPA1A HSPA1B HSPB1 DNAJB1 KLF6 ZFP36",
}


def key_of(name):
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")[:40] or "program"


def parse_genes(text):
    return [g for g in re.split(r"[\s,;]+", str(text or "").strip()) if g]


def resolve(name=None, genes=None, panel=()):
    """(name, genes asked, genes in the panel (panel spelling), missing) of a built-in program or the user's list."""
    if genes:
        asked = parse_genes(genes) if isinstance(genes, str) else list(genes)
        name = name or "My program"
    else:
        hit = next((k for k in BUILTIN if k.lower() == str(name or "").lower() or key_of(k) == key_of(name or "")), None)
        if hit is None:
            raise ValueError(f"no built-in program {name!r}; built-in: {', '.join(BUILTIN)}; or give the genes")
        name, asked = hit, BUILTIN[hit].split()
    up = {str(g).upper(): g for g in panel}
    used = list(dict.fromkeys(up[g.upper()] for g in asked if g.upper() in up))
    missing = [g for g in asked if g.upper() not in up]
    return name, asked, used, missing


def score(X, genes, program, seed=0):
    """score_genes on X (cells x genes, log-normalized; columns = genes): per cell, the mean over `program` (column
    indices) minus the mean over expression-matched controls."""
    program = np.asarray(sorted(set(int(i) for i in program)), int)
    mean = np.asarray(X.mean(0)).ravel()
    rank = pd.Series(mean).rank(method="min").to_numpy()
    nb = min(N_BINS, max(1, len(mean) // 4))
    bins = pd.cut(rank, nb, labels=False)
    rng = np.random.default_rng(seed)
    ctrl = set()
    pset = set(program.tolist())
    for b in np.unique(bins[program]):
        pool = [g for g in np.flatnonzero(bins == b) if g not in pset]
        if pool:
            ctrl.update(rng.choice(pool, min(CTRL, len(pool)), replace=False).tolist())
    ctrl = np.array(sorted(ctrl), int)
    if not len(ctrl):
        ctrl = np.array([g for g in range(len(mean)) if g not in pset], int)
    m = lambda cols: np.asarray(X[:, cols].mean(1)).ravel()
    return (m(program) - m(ctrl)).astype(np.float32), ctrl


def _dir(s):
    d = s.dir / "programs"
    d.mkdir(exist_ok=True)
    return d


def saved(s):
    """{key: info} of the programs scored on an annotation."""
    f = s.dir / "programs" / "programs.json"
    try:
        return json.load(open(f)) if f.exists() else {}
    except ValueError:
        return {}


def keep(s, key, info, values):
    """The scores of one annotation's cells under programs/<key>.npy."""
    d = _dir(s)
    replace_atomic(d / f"{key}.npy", lambda t: np.save(t, np.asarray(values, np.float32)))
    allp = saved(s)
    allp[key] = info
    dump_json(allp, d / "programs.json", indent=1)


def values(s, key):
    f = s.dir / "programs" / f"{key}.npy"
    if not f.exists():
        raise ValueError(f"the program {key!r} was not scored on this annotation (Compare -> Programs)")
    v = np.load(f)
    if len(v) != len(s.cells):
        raise ValueError(f"the program {key!r} was scored before the cells changed: score it again")
    return v


def run(c, name=None, genes=None, by="cell type", within=None, keep_for=None):
    """Score a program on the cells of `c` (analysis.Cells) and summarize it per value of `by` (mean score, share of
    cells above background); keep_for: the Sessions to keep the scores for (map colors)."""
    from . import analysis
    name, asked, used, missing = resolve(name, genes, c.genes)
    if len(used) < MIN_GENES:
        raise ValueError(f"{name}: only {len(used)} of its {len(asked)} genes are in the panel "
                         f"({', '.join(used) or 'none'}; at least {MIN_GENES})")
    gi = {g: i for i, g in enumerate(c.genes)}
    sc, ctrl = score(c.X, c.genes, [gi[g] for g in used])
    key = key_of(name)
    info = {"name": name, "genes": used, "missing": missing, "n_controls": int(len(ctrl)),
            "run_id": time.strftime("%Y%m%d-%H%M%S") + f"-{int(time.time() * 1000) % 1000:03d}"}
    off = np.cumsum([0] + [len(s.cells) for s in c.sessions])
    for i, s in enumerate(c.sessions):
        if keep_for is None or s in keep_for:
            keep(s, key, info, sc[off[i]:off[i + 1]])
    w = analysis._within(c, within)
    col = analysis._column(c, by)
    m = w & pd.notna(col)
    df = pd.DataFrame({"g": col[m].astype(str), "v": sc[m]})
    tab = df.groupby("g")["v"].agg(["mean", "median", "count"])
    tab["above"] = df.assign(a=df.v > 0).groupby("g")["a"].mean()
    tab = tab.sort_values("mean", ascending=False)
    q = np.quantile(sc[m], [0.5, 0.9, 0.99]) if m.any() else [0, 0, 0]
    return {"scope": "project" if c.pid else "slice", "program": name, "key": key, "mode": f"program:{key}",
            "run_id": info["run_id"],
            "genes": used, "missing": missing, "by": by, "within": within or {},
            "n_cells": int(m.sum()), "quantiles": {"50%": round(float(q[0]), 3), "90%": round(float(q[1]), 3),
                                                   "99%": round(float(q[2]), 3)},
            "groups": [{"group": g, "mean": round(float(r["mean"]), 3), "median": round(float(r["median"]), 3),
                        "pct_above": round(float(r["above"]), 3), "n_cells": int(r["count"])} for g, r in tab.iterrows()],
            "table": tab.reset_index().rename(columns={"g": by, "above": "pct_above", "count": "n_cells"}),
            "note": (f"Score = mean log-normalized expression of the program's {len(used)} genes minus that of "
                     f"expression-matched control genes (as scanpy's score_genes); above 0 = higher than background. "
                     + (f"Not in the panel: {', '.join(missing[:12])}{'…' if len(missing) > 12 else ''}. " if missing else "")
                     + "A short program on a targeted panel rests on few genes: check them with DE genes or the map.")}


BINS = (("low", "at or below background", (200, 200, 200), 40), ("q1", "above background", (254, 217, 118), 225),
        ("q2", "top half above background", (253, 141, 60), 235), ("q3", "top 20%", (227, 26, 28), 245),
        ("q4", "top 5%", (128, 0, 38), 250))


def bins(v):
    """Bin (index into BINS) of every cell: at or below 0, then quantiles of the scores above 0."""
    pos = v[v > 0]
    if not len(pos):
        return np.zeros(len(v), int)
    t = np.quantile(pos, [0.5, 0.8, 0.95])
    return np.where(v <= 0, 0, 1 + np.digitize(v, t)).astype(int)
