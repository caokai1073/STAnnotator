"""Cell type and region vocabulary, named groups, and name resolution with suggestions."""
import difflib

# Canonical region names of the pathologist annotation (user-drawn regions can use any name).
REGIONS = ["mucosa", "submucosa", "muscularis_propria", "lymphoid_aggregate",
           "fibromuscular_hyperplasia", "granuloma", "unknown"]
# Regions where no constraint is ever applied.
FREE_REGIONS = ["unknown"]

PATHOLOGIST_TO_REGION = {
    "Normal": "mucosa", "AC": "mucosa", "CC": "mucosa", "Submucosa": "submucosa",
    "Musc.Propria": "muscularis_propria", "LAgg": "lymphoid_aggregate",
    "FMHP": "fibromuscular_hyperplasia", "Granuloma": "granuloma", "Unlabelled": "unknown",
}
PATHOLOGIST_DETAIL = {"Normal": "normal mucosa", "AC": "mucosa with active colitis",
                      "CC": "mucosa with chronic colitis", "FMHP": "fibromuscular hyperplasia"}

LINEAGE_GROUPS = {"epi": "epithelial", "fib": "fibroblast", "lymph": "lymphoid",
                  "myeloid": "myeloid", "stroma": "stromal"}

EXTRA_TYPE_GROUPS = {
    "mp_smc": ["SMC (inner muscularis propria)", "SMC (outer muscularis propria)", "ICC"],
    "smc": ["SMC (inner muscularis propria)", "SMC (outer muscularis propria)", "SMC (muscularis mucosae)"],
}

VERDICTS = {
    "constraint_right": "the spatial prior is right; the expression-based (reference) label is wrong",
    "reference_right": "the expression-based label is right; the prior does not apply here",
    "region_wrong": "the cell's region assignment is wrong",
    "artifact": "technical artifact (tissue edge, debris, segmentation, no tissue on H&E)",
    "new_state": "real cell state not covered by the current cell type taxonomy",
    "unclear": "cannot decide",
}


class VocabError(ValueError):
    pass


def type_groups(types, lineage_of_type):
    groups = {name: sorted(t for t in types if lineage_of_type[t] == lin) for lin, name in LINEAGE_GROUPS.items()}
    groups.update({k: [t for t in v if t in types] for k, v in EXTRA_TYPE_GROUPS.items()})
    return groups


def _suggest(name, options):
    options = list(dict.fromkeys(options))
    m = difflib.get_close_matches(name, options, n=3, cutoff=0.4)
    if not m:
        low = name.lower()
        m = [o for o in options if low in o.lower()][:3]
    return m


def resolve(names, vocabulary, groups, kind):
    """Expand group names and validate items (case-insensitive); raise VocabError with suggestions on unknown names."""
    if isinstance(names, str):
        names = [names]
    vocabulary = list(vocabulary)
    by_lower = {v.lower(): v for v in vocabulary}
    groups_lower = {g.lower(): g for g in groups}
    out, errors = [], []
    for n in names:
        if n in groups or n.lower() in groups_lower:
            out += groups[n if n in groups else groups_lower[n.lower()]]
        elif n in vocabulary or n.lower() in by_lower:
            out.append(n if n in vocabulary else by_lower[n.lower()])
        else:
            errors.append(f"unknown {kind} '{n}'; did you mean: {_suggest(n, vocabulary + list(groups))}")
    if errors:
        raise VocabError("; ".join(errors))
    return sorted(set(out), key=vocabulary.index)
