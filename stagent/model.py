"""Expression model: multinomial logistic regression on log-normalized counts (numpy; small enough for the CPU).

Two sources of training labels: the reference pool (reference mode) or the sample's own marker-annotation
pseudo-labels (marker mode, no reference needed).

Reference pool (work/reference/): a stratified sample of annotated cells from all slides (at most CAP cells
per slide and type), plus the natural type counts per slide. Training for a sample excludes that sample's
slide (leave-one-slide-out) and corrects the class prior back to the natural distribution, so the stratified
sampling does not inflate rare types. User corrections enter as extra, up-weighted training cells.
"""
import json
import numpy as np, pandas as pd, scipy.sparse as sp
from . import config
from .osutil import dump_json

REF = config.HOME / "reference"  # reference pool (reference mode / benchmarks on the lab server only)
CAP = 5000
EPOCHS, BATCH, LR = 6, 8192, 3e-3
CORRECTION_WEIGHT = 20.0
MIN_TRAIN = 10  # fewest labeled cells the expression model is trained on (10% of them validate its temperature)
# fine-tuning from the initial model when the user corrects cells: changes stay local to what was corrected
# only the classes involved in the corrections may move, and an L2 pull keeps them near the initial model
FT_STEPS, FT_LR, FT_REPLAY, FT_L2 = 100, 5e-4, 6144, 10.0


LOGNORM_BLOCK = 20_000  # cells normalized at once (the repeated per-cell factors are one block long: keep it small)


def lognorm(m):
    """log1p of the counts per 100 per cell. The result SHARES the input's column indices (only the values are new),
    so it costs half of what a copy would: a five-million-cell CosMx 6K matrix is some 13 GB, half of it the indices.
    The values are the same to the bit as diags(100 / total) @ m with log1p over its data. The input is not changed;
    nothing may sort or prune the result's indices in place, as that would reach the matrix it came from."""
    m = sp.csr_matrix(m, dtype=np.float32)
    tot = np.asarray(m.sum(1)).ravel()
    tot[tot == 0] = 1
    if not m.has_canonical_format:  # duplicate entries for one cell and gene, which the product below adds up:
        m = sp.diags(100.0 / tot) @ m      # (a sample's counts never have any -- shared_counts sums them at load)
        m.data = np.log1p(m.data)
        return m
    inv = np.asarray(100.0 / tot, dtype=m.dtype)
    out = np.empty_like(m.data)
    for a in range(0, m.shape[0], LOGNORM_BLOCK):  # a block at a time: no temporary the size of the matrix
        b = min(a + LOGNORM_BLOCK, m.shape[0])
        lo, hi = m.indptr[a], m.indptr[b]
        np.multiply(m.data[lo:hi], np.repeat(inv[a:b], np.diff(m.indptr[a:b + 1])), out=out[lo:hi])
        np.log1p(out[lo:hi], out=out[lo:hi])
    return m._with_data(out, copy=False)


def build_reference():
    """Write the reference pool from work/obs_all.parquet and work/counts.npz."""
    W = config.WORK
    obs = pd.read_parquet(W / "obs_all.parquet", columns=["celltype", "slide_id"])
    types = list(np.load(W / "types.npy", allow_pickle=True))
    X = sp.load_npz(W / "counts.npz").tocsr()
    y = pd.Categorical(obs.celltype, categories=types).codes
    rng = np.random.default_rng(0)
    keep = []
    for _, idx in pd.Series(np.arange(len(obs))).groupby([obs.slide_id.values, y]):
        idx = idx.to_numpy()
        keep.append(rng.choice(idx, min(CAP, len(idx)), replace=False))
    keep = np.sort(np.concatenate(keep))
    REF.mkdir(parents=True, exist_ok=True)
    sp.save_npz(REF / "counts.npz", X[keep])
    np.save(REF / "y.npy", y[keep].astype(np.int16))
    np.save(REF / "slide.npy", obs.slide_id.values[keep])
    natural = pd.crosstab(obs.slide_id, y).reindex(columns=range(len(types)), fill_value=0)
    dump_json({"types": types, "slides": [int(s) for s in natural.index],
               "natural_counts": natural.to_numpy().tolist()}, REF / "meta.json")
    print("reference pool:", len(keep), "cells")


class Reference:
    def __init__(self):
        self.X = lognorm(sp.load_npz(REF / "counts.npz"))
        self.y = np.load(REF / "y.npy").astype(np.int64)
        self.slide = np.load(REF / "slide.npy")
        m = json.load(open(REF / "meta.json"))
        self.types = m["types"]
        self.natural = pd.DataFrame(m["natural_counts"], index=m["slides"])


_ref = None


def reference():
    global _ref
    if _ref is None:
        _ref = Reference()
    return _ref


class _Adam:
    """Adam on a list of numpy arrays, updated in place (weight decay added to the gradient, as in torch.optim.Adam)."""

    def __init__(self, params, lr, weight_decay=0.0, betas=(0.9, 0.999), eps=1e-8):
        self.p, self.lr, self.wd, (self.b1, self.b2), self.eps, self.t = params, lr, weight_decay, betas, eps, 0
        self.m = [np.zeros_like(x) for x in params]
        self.v = [np.zeros_like(x) for x in params]

    def step(self, grads):
        self.t += 1
        c1, c2 = 1 - self.b1 ** self.t, 1 - self.b2 ** self.t
        for x, g, m, v in zip(self.p, grads, self.m, self.v):
            if self.wd:
                g = g + self.wd * x
            m *= self.b1; m += (1 - self.b1) * g
            v *= self.b2; v += (1 - self.b2) * g * g
            x -= self.lr * (m / c1) / (np.sqrt(v / c2) + self.eps)


def _ce(Z, y, w=None):
    """Weighted mean cross-entropy of logits Z and its gradient with respect to Z."""
    Z = Z - Z.max(1, keepdims=True)
    P = np.exp(Z)
    P /= P.sum(1, keepdims=True)
    w = np.ones(len(y), np.float32) if w is None else w
    r = np.arange(len(y))
    loss = float(-(w * np.log(P[r, y] + 1e-12)).sum() / w.sum())
    P[r, y] -= 1
    P *= (w / w.sum())[:, None]
    return loss, P


def _genes(X, keep=None):
    """The genes the model learns from: all of a panel; for whole-transcriptome data (more than
    markers.MAX_PCA_GENES genes, e.g. Visium HD) the markers.N_HVG most variable of the training cells
    (log-normalized) and the genes `keep` (the marker list's), which keeps the training matrix small."""
    from . import markers
    if X.shape[1] <= markers.MAX_PCA_GENES:
        return None
    mu = np.asarray(X.mean(0)).ravel()
    var = np.asarray(X.multiply(X).mean(0)).ravel() - mu ** 2
    return np.array(sorted(set(np.argsort(-var)[:markers.N_HVG].tolist()) | set(keep or [])), np.int64)


def _fit(X, y, T, pi_target, extra_X=None, extra_y=None, progress=None, seed=0, genes=None):
    """Core trainer. X: log-normalized CSR, y: class indices; 10% of (X, y) is held out for temperature scaling.
    pi_target: class distribution to correct the prior back to (training data may be class-capped). genes: the
    columns to learn from (default: _genes); the model keeps them ("genes") and predict() uses only them."""
    from scipy.optimize import minimize_scalar
    if len(y) < MIN_TRAIN:
        raise ValueError(f"only {len(y)} labeled cells to learn the expression model from (at least {MIN_TRAIN}): "
                         "check that the labels match the sample's cells, and that not all of them are Unassigned")
    genes = _genes(X) if genes is None else genes
    if genes is not None:
        X = X[:, genes]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(y))
    n_val = max(min(1000, len(y) // 5), len(y) // 10)
    va, tr = perm[:n_val], perm[n_val:]
    Xtr, ytr = X[tr], y[tr]
    wtr = np.ones(len(ytr), np.float32)
    if extra_X is not None and len(extra_y):
        ex = lognorm(extra_X)
        Xtr = sp.vstack([Xtr, ex if genes is None else ex[:, genes]]).tocsr()
        ytr = np.concatenate([ytr, np.asarray(extra_y, np.int64)])
        wtr = np.concatenate([wtr, np.full(len(extra_y), CORRECTION_WEIGHT, np.float32)])
    pi_train = np.bincount(ytr, weights=wtr, minlength=T) + 1.0
    pi_target = np.asarray(pi_target, float) + 1.0
    prior_shift = np.log(pi_target / pi_target.sum()) - np.log(pi_train / pi_train.sum())

    G = Xtr.toarray().astype(np.float32, copy=False)  # toarray() is float32 already: no second dense copy
    mu = G.mean(0)
    sd = G.std(0, ddof=1) + 1e-3
    G -= mu
    G /= sd
    bound = 1 / np.sqrt(G.shape[1])
    W = rng.uniform(-bound, bound, (T, G.shape[1])).astype(np.float32)
    b = rng.uniform(-bound, bound, T).astype(np.float32)
    opt = _Adam([W, b], LR, weight_decay=1e-6)
    steps = EPOCHS * ((len(ytr) + BATCH - 1) // BATCH)
    k = 0
    for ep in range(EPOCHS):
        idx = rng.permutation(len(ytr))
        for s in range(0, len(ytr), BATCH):
            bi = idx[s:s + BATCH]
            Xb = G[bi]
            _, dZ = _ce(Xb @ W.T + b, ytr[bi], wtr[bi])
            opt.step([dZ.T @ Xb, dZ.sum(0)])
            k += 1
            if progress and k % 20 == 0:
                progress(k / steps, f"training expression model, epoch {ep + 1}/{EPOCHS}")
    del G
    model = {"W": W, "b": b, "mu": mu, "sd": sd, "prior_shift": prior_shift.astype(np.float32), "temp": 1.0}
    lv, yv = _logits(model, X[va]), y[va]
    res = minimize_scalar(lambda lt: _ce(lv / np.exp(lt), yv)[0], bounds=(-4.0, 4.0), method="bounded")
    model["temp"] = float(np.exp(res.x))
    if genes is not None:
        model["genes"] = genes
    return model


def train(exclude_slide, extra_X=None, extra_y=None, progress=None, seed=0):
    """Reference mode: train on the reference pool without `exclude_slide`, plus optional corrected cells."""
    ref = reference()
    m = ref.slide != exclude_slide
    pi_nat = ref.natural.drop(index=exclude_slide, errors="ignore").sum(0).to_numpy()
    return _fit(ref.X[m], ref.y[m], len(ref.types), pi_nat, extra_X, extra_y, progress, seed)


def pseudo_pool(counts, labels, T, seed=0):
    """Marker mode: training pool from the sample's own pseudo-labeled cells (at most CAP per type).
    Labels < 0 (cells the marker step left at a lineage or parent level) are not used for training."""
    rng = np.random.default_rng(seed)
    if not (np.asarray(labels) >= 0).any():
        raise ValueError("no cell has a cell type to learn the expression model from (all Unassigned, or the labels "
                         "match none of the sample's cells)")
    keep = np.concatenate([rng.choice(np.where(labels == t)[0], min(CAP, (labels == t).sum()), replace=False)
                           for t in range(T) if (labels == t).any()])
    keep.sort()
    return lognorm(counts[keep]), np.asarray(labels[keep], np.int64), np.bincount(labels[labels >= 0], minlength=T)


def train_pseudo(counts, labels, T, progress=None, seed=0, keep=None):
    """Marker mode: learn p(type | expression) from marker-annotation pseudo-labels of this sample. keep: columns
    of genes always used (the marker list's), when the data has more genes than a panel (_genes)."""
    from . import timing
    X, y, natural = pseudo_pool(counts, labels, T, seed)
    genes = _genes(X, keep)
    with timing.step("model", len(y) * (X.shape[1] if genes is None else len(genes)) * T / 1e12):
        return _fit(X, y, T, natural, progress=progress, seed=seed, genes=genes)


def finetune(init, extra_X, extra_y, replay_X, replay_y, progress=None, seed=0):
    """Start from `init` and adapt to corrected cells, replaying (replay_X, replay_y) so the rest is preserved.
    Only the weights of involved classes (the corrected labels and what init predicted for those cells) are
    updated, so unrelated cell types cannot drift. Normalization, prior correction and temperature are kept."""
    rng = np.random.default_rng(seed)
    T = len(init["b"])
    mu, sd = init["mu"], init["sd"]
    g = init.get("genes")
    sub = (lambda m: m) if g is None else (lambda m: m[:, g])
    W, b = init["W"].astype(np.float32).copy(), init["b"].astype(np.float32).copy()
    W0, b0 = W.copy(), b.copy()
    opt = _Adam([W, b], FT_LR)
    Xc = (sub(lognorm(extra_X)).toarray() - mu) / sd
    yc = np.asarray(extra_y, np.int64)
    before = (Xc @ W.T + b + init["prior_shift"]).argmax(1)  # the classes init actually predicted
    involved = np.zeros(T, bool)
    involved[yc] = True
    involved[before] = True
    n_rep = min(FT_REPLAY, len(replay_y))
    n_c = len(yc) * CORRECTION_WEIGHT
    for k in range(FT_STEPS):
        rep = rng.choice(len(replay_y), n_rep, replace=False)
        Xr = (sub(replay_X[rep]).toarray() - mu) / sd
        yr = np.asarray(replay_y[rep], np.int64)
        _, dr = _ce(Xr @ W.T + b, yr)
        _, dc = _ce(Xc @ W.T + b, yc)
        a_r, a_c = n_rep / (n_rep + n_c), n_c / (n_rep + n_c)
        gW = a_r * dr.T @ Xr + a_c * dc.T @ Xc + 2 * FT_L2 * (W - W0)
        gb = a_r * dr.sum(0) + a_c * dc.sum(0) + 2 * FT_L2 * (b - b0)
        gW[~involved] = 0
        gb[~involved] = 0
        opt.step([gW.astype(np.float32), gb.astype(np.float32)])
        if progress and k % 10 == 0:
            progress(k / FT_STEPS, f"updating expression model with your corrections ({k}/{FT_STEPS})")
    return dict(init, W=W, b=b)


def reference_replay(exclude_slide):
    ref = reference()
    m = ref.slide != exclude_slide
    return ref.X[m], ref.y[m]


def save(model, path):
    from .osutil import replace_atomic
    replace_atomic(path, lambda t: np.savez(t, **{k: np.asarray(v) for k, v in model.items()}))


def load(path):
    z = np.load(path)
    return {k: (float(z[k]) if z[k].ndim == 0 else z[k]) for k in z.files}


def _logits(model, Xn, chunk=65536):
    out, g = [], model.get("genes")
    for i in range(0, Xn.shape[0], chunk):
        part = Xn[i:i + chunk]
        a = ((part if g is None else part[:, g]).toarray() - model["mu"]) / model["sd"]
        out.append(a @ model["W"].T + model["b"] + model["prior_shift"])
    return np.concatenate(out).astype(np.float32)


def predict(model, counts):
    """Calibrated log-probabilities (float16) for raw counts."""
    z = _logits(model, lognorm(counts)) / model["temp"]
    z -= z.max(1, keepdims=True)
    return (z - np.log(np.exp(z).sum(1, keepdims=True))).astype(np.float16)
