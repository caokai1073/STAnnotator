"""Constraint DSL: validation, normalization and Bayesian reweighting of base probabilities.

A constraint is a dict:
  id      e.g. "C1" (assigned by the session if missing)
  text    the user's own words, kept verbatim for traceability
  author  who added it (optional)
  kind    what the condition is about (default "region"):
            region        cell lies in one of `regions`
            depth         cell lies within `max_depth_um` of the mucosa core (epithelial-marker density;
                          0 inside crypts / surface epithelium)
            neighborhood  >= `min_frac` of the cell's K nearest neighbors are of `neighbor_types`
                          (neighbors typed by the expression-only prediction, so the rule never feeds on itself)
  types   cell types and/or type group names the rule is about
  mode    only    -> types are unlikely where the condition is FALSE (weight < 1 applied there)
          exclude -> types are unlikely where the condition is TRUE  (weight < 1 applied there)
          enrich  -> types are more likely where the condition is TRUE (weight > 1 applied there)
  weight  multiplicative prior factor; strength names map to: strong 0.01, moderate 0.1, weak 0.5
          (enrich: strong 10, moderate 3, weak 1.5)
Region-kind constraints never touch cells in vocab.FREE_REGIONS. When a region a constraint uses is deleted, its
cells at that moment are kept in the constraint (cells_of: {region name: [cell ids]}), so the prior keeps applying
to them.
posterior(t | x, cell) is proportional to p(t | x) * prod_c w_c(t, cell).
"""
import numpy as np
from . import vocab

KINDS = ("region", "depth", "neighborhood")
MODES = ("only", "exclude", "enrich")
STRENGTH = {"only": {"strong": 0.01, "moderate": 0.1, "weak": 0.5},
            "exclude": {"strong": 0.01, "moderate": 0.1, "weak": 0.5},
            "enrich": {"strong": 10.0, "moderate": 3.0, "weak": 1.5}}
NEIGHBORS_K = 15
SPEC_FIELDS = ("id", "text", "author", "kind", "types", "regions", "max_depth_um", "neighbor_types", "min_frac",
               "mode", "weight")


def normalize(c, types, type_groups, region_names=None):
    """Validate a constraint and return a normalized copy with expanded type/region lists.
    region_names: regions that exist in the session (user-drawn names); defaults to the canonical vocabulary."""
    c = {k: v for k, v in dict(c).items() if v is not None and not k.endswith("_resolved")}
    c.setdefault("kind", "region")
    if c.get("kind", "region") not in KINDS:
        raise vocab.VocabError(f"kind must be one of {KINDS}, got {c['kind']!r}")
    if c.get("mode") not in MODES:
        raise vocab.VocabError(f"mode must be one of {MODES}, got {c.get('mode')!r}")
    w = c.get("weight", "strong")
    if isinstance(w, str):
        if w not in STRENGTH[c["mode"]]:
            raise vocab.VocabError(f"weight must be a number or one of {list(STRENGTH[c['mode']])}")
        w = STRENGTH[c["mode"]][w]
    w = float(w)
    if c["mode"] == "enrich" and w <= 1:
        raise vocab.VocabError("enrich needs weight > 1")
    if c["mode"] != "enrich" and not 0 < w < 1:
        raise vocab.VocabError(f"{c['mode']} needs 0 < weight < 1")
    c["weight"] = w
    if not c.get("types"):
        raise vocab.VocabError("types is required")
    c["types_resolved"] = vocab.resolve(c["types"], types, type_groups, "cell type")
    if c.get("kind", "region") == "region":
        if not c.get("regions") and not c.get("cells_of"):
            raise vocab.VocabError("region constraints need regions")
        c["regions_resolved"] = vocab.resolve(c["regions"], region_names if region_names is not None else vocab.REGIONS,
                                              {}, "region") if c.get("regions") else []
        if set(c["regions_resolved"]) & set(vocab.FREE_REGIONS):
            raise vocab.VocabError(f"constraints cannot target {vocab.FREE_REGIONS}")
    elif c.get("kind", "region") == "depth":
        if not c.get("max_depth_um") or float(c["max_depth_um"]) <= 0:
            raise vocab.VocabError("depth constraints need max_depth_um > 0")
        c["max_depth_um"] = float(c["max_depth_um"])
    else:
        if not c.get("neighbor_types"):
            raise vocab.VocabError("neighborhood constraints need neighbor_types")
        c["neighbor_types_resolved"] = vocab.resolve(c["neighbor_types"], types, type_groups, "cell type")
        c["min_frac"] = float(c.get("min_frac", 0.2))
        if not 0 < c["min_frac"] <= 1:
            raise vocab.VocabError("min_frac must be in (0, 1]")
    c.setdefault("text", "")
    return c


def _list(names):
    return ", ".join(names) if len(names) <= 4 else f"{len(names)} types ({', '.join(names[:3])}, ...)"


def condition_text(c):
    if c.get("kind", "region") == "region":
        kept = [f"{len(v):,} cells of deleted region {k}" for k, v in c.get("cells_of", {}).items()]
        return f"in [{', '.join(c['regions_resolved'] + kept)}]"
    if c.get("kind", "region") == "depth":
        return f"within {c['max_depth_um']:g} um of the epithelium"
    return f"with >= {c['min_frac']:.0%} of {NEIGHBORS_K} neighbors being [{_list(c['neighbor_types_resolved'])}]"


def describe(c):
    cond = condition_text(c)
    verb = {"only": f"only found {cond}; elsewhere prior x{c['weight']:g}",
            "exclude": f"not found {cond}; there prior x{c['weight']:g}",
            "enrich": f"enriched {cond}; there prior x{c['weight']:g}"}[c["mode"]]
    return f"{c['id']}: [{_list(c['types_resolved'])}] {verb}"


def condition_mask(c, ctx):
    """Boolean per cell: is the constraint's condition true?"""
    if c.get("kind", "region") == "region":
        m = np.isin(ctx["region"], c["regions_resolved"])
        if c.get("cells_of"):
            m |= _kept_cells(c, ctx)
        return m
    if c.get("kind", "region") == "depth":
        return ctx["depth_um"] <= c["max_depth_um"]
    return ctx["neighbor_frac"](c["neighbor_types_resolved"]) >= c["min_frac"]


def penalized_cells(c, ctx):
    """Boolean per cell: where this constraint multiplies its weight."""
    cond = condition_mask(c, ctx)
    m = ~cond if c["mode"] == "only" else cond
    if c.get("kind", "region") == "region":
        free = np.isin(ctx["region"], vocab.FREE_REGIONS)
        if c.get("cells_of"):
            free &= ~_kept_cells(c, ctx)  # cells of a deleted region are in no region now but still constrained
        m &= ~free
    return m


def _kept_cells(c, ctx):
    return np.isin(ctx["cell"], [x for v in c["cells_of"].values() for x in v])


def posterior(base_logp, constraints, ctx, types, top=0):
    """Return (pred index, confidence) after applying constraints; ctx holds per-cell context.
    With top=k also return the k most likely types per cell and their probabilities."""
    lp = base_logp.astype(np.float32)
    for c in constraints:
        ti = [types.index(t) for t in c["types_resolved"]]
        m = penalized_cells(c, ctx)
        lp[np.ix_(m, ti)] += np.float32(np.log(c["weight"]))
    lp -= lp.max(1, keepdims=True)
    p = np.exp(lp)
    p /= p.sum(1, keepdims=True)
    if not top:
        return p.argmax(1), p.max(1)
    idx = np.argsort(-p, axis=1)[:, :top]
    return idx[:, 0], np.take_along_axis(p, idx[:, :1], 1)[:, 0], idx.astype(np.int16), \
        np.take_along_axis(p, idx, 1).astype(np.float16)
