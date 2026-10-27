"""The agent's review of each part of the analysis, and the annotation report (stored per session).

The user picks what to review (Agent Review): one scope of SCOPES, each with its categories of advice:
  annotation     cell types: types to merge, split or add (taxonomy), the sample's site and disease (context),
                 technical problems such as registration, spillover, segmentation (qc); besides the computed review
                 candidates the agent checks (claudereview.py)
  organization   compartments seen on the image and in the cells, proposed as regions (compartments); gradients and
                 neighbors that are expected or contradict the labels (organization)
  embedding      Leiden clusters against the labels, slices mixing in a joint embedding, its parameters
  niches         what the niches are, their names, k and their number
  compare        the DE genes and compositions the user looked at: findings, design and caveats, next analyses
  signaling      the cell-cell communication the user looked at: interactions, design and caveats, next analyses
The agent reviews its own advice before the user sees it: draft() stores drafts (not shown), check() keeps the
drafts that hold up against the data (with the evidence) and files the others as filtered out (with the reason,
folded away in the panel). Kept items replace the scope's earlier open advice.
Each item can carry one action the user takes in the viewer with one click:
  region         {"kind": "region", "name", "compartment"}: add the compartment outline (from preview_compartment)
                 as a region; an optional "prior" (constraint spec) is offered once the region exists
  prior          {"kind": "prior", "prior": spec}: preview, then add the prior
  marker_list    {"kind": "marker_list"}: open the marker list Claude proposed (propose_marker_list)
  registration   {"kind": "registration"}: correct the cells' place on the image (insights.registration)
  claude         {"kind": "claude", "request": "..."}: ask Claude Code to do it (e.g. a re-clustering), as a request
  niche_names    {"kind": "niche_names", "names": {"N1": "...", ...}}: name the niches
  show           {"kind": "show", "mode": "gene:<gene>" | "lr:<ligand>~<receptor>~receiver|sender" | "niche" | ...}:
                 color the map
  rerun          {"kind": "rerun", "tab": "embed" | "niches" | "compare" | "signal", "params": {...}}: put these
                 settings in the tab (the user clicks its Run button)
Advice is advice: nothing changes until the user clicks. A new set of advice replaces the open items; done and
dismissed ones are kept. The report is Markdown Claude writes (write_report); the viewer shows it and saves it.
"""
import hashlib, json, re, time
import numpy as np
from . import constraints as cmod

SCOPES = {
    "annotation": {"label": "Annotation", "categories": {"taxonomy": "Cell types", "context": "Sample context",
                                                         "qc": "Technical"}},
    "organization": {"label": "Spatial organization", "categories": {"compartments": "Compartments",
                                                                     "organization": "Gradients and neighbors"}},
    "embedding": {"label": "Embedding", "categories": {"clusters": "Clusters and labels", "integration": "Slices",
                                                       "parameters": "Parameters"}},
    "niches": {"label": "Niches", "categories": {"niches": "Niches", "parameters": "Parameters"}},
    "compare": {"label": "Compare", "categories": {"findings": "Findings", "design": "Design and caveats",
                                                   "next": "Next analyses"}},
    "signaling": {"label": "Signaling", "categories": {"interactions": "Interactions", "design": "Design and caveats",
                                                       "next": "Next analyses"}},
}
CATEGORIES = {c: l for sc in SCOPES.values() for c, l in sc["categories"].items()}
ACTIONS = ("region", "prior", "marker_list", "registration", "claude", "niche_names", "show", "rerun")
RERUN_TABS = ("embed", "niches", "compare", "signal")
MAX_ITEMS = 12


def scope_of(item):
    """An item's scope (advice made before scopes existed: by its category)."""
    if item.get("scope"):
        return item["scope"]
    return "organization" if item.get("category") in ("compartments", "organization") else "annotation"


def _file(s):
    return s.dir / "advice.json"


def load(s):
    p = _file(s)
    return json.load(open(p)) if p.exists() else {"items": [], "runs": [], "requests": [], "compartments": {}}


def _save(s, st):
    json.dump(st, open(_file(s), "w"), indent=1)


def labels_key(s, version=None):
    """Fingerprint of a version's labels: advice is redone only when the labels changed."""
    return hashlib.sha1(np.ascontiguousarray(s.result(version)["pred"]).tobytes()).hexdigest()[:16]


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%S")


# ---------------- compartment outlines Claude previewed ----------------
def store_compartment(s, name, comp, spec):
    st = load(s)
    st["compartments"][name] = {"polygons_px": comp["polygons_px"], "n_cells": comp["n_cells"], "spec": spec,
                                "time": _now()}
    _save(s, st)


# ---------------- advice ----------------
def _next_id(used):
    k = 1
    while f"H{k}" in used:
        k += 1
    used.add(f"H{k}")
    return f"H{k}"


def propose_region(s, name, polygons_px, n_cells, text, prior=None, source=None, version=None):
    """A region the user asked Claude for ("draw the mucosa"): added to the advice without touching the other items
    or the advice runs, and shown in the Review panel, where the user adds it as a region with one click."""
    name = re.sub(r"\s+", "_", str(name).strip())
    if not name:
        raise ValueError("the region needs a name")
    if not polygons_px:
        raise ValueError("no outline")
    st = load(s)
    act = {"kind": "region", "name": name, "compartment": source, "polygons_px": polygons_px, "n_cells": int(n_cells)}
    if prior:
        act["prior"], act["prior_text"] = _check_prior(s, prior, name)
    item = {"id": _next_id({x["id"] for x in st["items"]}), "scope": "organization", "category": "compartments", "requested": True,
            "title": f"Region “{name}” (you asked for it)", "text": str(text).strip() or f"Outline of {name}.",
            "action": act, "status": "open", "version": s.state(version)["version"], "time": _now()}
    st["items"].append(item)
    _save(s, st)
    return item

def _check_prior(s, spec, new_region=None):
    """A prior spec as the user will add it; regions may include the region the same item creates."""
    spec = {k: v for k, v in dict(spec).items() if v is not None}
    names = s.region_names() + ([new_region] if new_region else [])
    c = cmod.normalize(dict(spec, id="C?"), s.types, s.type_groups, names)
    rule = cmod.describe(c).split(": ", 1)[1]
    return {k: spec[k] for k in cmod.SPEC_FIELDS if k in spec and k != "id"}, \
        (f'{spec["text"]} ({rule.split("] ", 1)[-1]})' if spec.get("text") else rule)


def _action(s, st, act):
    """An item's action as stored (validated; outlines and corrections filled in)."""
    act = dict(act or {})
    kind = act.get("kind")
    if kind is not None and kind not in ACTIONS:
        raise ValueError(f"action kind must be one of {ACTIONS}, got {kind!r}")
    out = {"kind": kind} if kind else None
    if kind == "region":
        comp = st["compartments"].get(str(act.get("compartment", "")))
        if comp is None:
            raise ValueError(f"no previewed compartment {act.get('compartment')!r}: call preview_compartment first")
        name = re.sub(r"\s+", "_", str(act.get("name") or act["compartment"]).strip())
        out |= {"name": name, "compartment": act["compartment"], "polygons_px": comp["polygons_px"],
                "n_cells": comp["n_cells"]}
        if act.get("prior"):
            out["prior"], out["prior_text"] = _check_prior(s, act["prior"], name)
    elif kind == "prior":
        if not act.get("prior"):
            raise ValueError("a prior action needs 'prior' (a constraint spec, as for preview_constraint)")
        out["prior"], out["prior_text"] = _check_prior(s, act["prior"])
    elif kind == "registration":
        from . import insights
        fit = insights.registration(s).get("affine_fit")
        if not fit:
            raise ValueError("no registration fit for this sample")
        out |= {"correction_px": fit["correction_px"], "max_offset_um": fit["max_offset_um"]}
    elif kind == "claude":
        if not str(act.get("request", "")).strip():
            raise ValueError("a claude action needs 'request' (what Claude Code should do when the user asks)")
        out["request"] = str(act["request"]).strip()
    elif kind == "niche_names":
        from . import niches
        meta, _ = niches.current(s)
        if meta is None:
            raise ValueError("no niches to name")
        names = {str(k): str(v).strip()[:40] for k, v in dict(act.get("names") or {}).items()}
        bad = [k for k in names if k not in {f"N{i + 1}" for i in range(meta["n_niches"])}]
        if not names or bad:
            raise ValueError(f"names must map niches N1..N{meta['n_niches']} to names" + (f" (unknown: {bad})" if bad else ""))
        out["names"] = names
    elif kind == "show":
        mode = str(act.get("mode", "")).strip()
        if not mode:
            raise ValueError("a show action needs 'mode' (e.g. gene:EPCAM, lr:CCL19~CCR7~receiver, niche)")
        out["mode"] = mode
    elif kind == "rerun":
        if act.get("tab") not in RERUN_TABS:
            raise ValueError(f"a rerun action needs 'tab', one of {RERUN_TABS}")
        out |= {"tab": act["tab"], "params": dict(act.get("params") or {})}
    return out


def draft(s, scope, items, version=None):
    """The agent's draft advice on one scope (replaces its earlier drafts there); not shown until checked."""
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {list(SCOPES)}")
    if not items:
        raise ValueError("no drafts")
    if len(items) > MAX_ITEMS:
        raise ValueError(f"at most {MAX_ITEMS} drafts: keep the ones that matter")
    st = load(s)
    st["items"] = [x for x in st["items"] if not (x["status"] == "draft" and scope_of(x) == scope)]
    used = {x["id"] for x in st["items"]}
    v = s.state(version)["version"]
    cats = SCOPES[scope]["categories"]
    new = []
    for it in items:
        cat = str(it.get("category", "")).strip().lower()
        if cat not in cats:
            raise ValueError(f"category of {scope} advice must be one of {list(cats)}, got {cat!r}")
        title, text = str(it.get("title", "")).strip(), str(it.get("text", "")).strip()
        if not title or not text:
            raise ValueError("every piece of advice needs a title and a text")
        new.append({"id": _next_id(used), "scope": scope, "category": cat, "title": title, "text": text,
                    "action": _action(s, st, it.get("action")), "status": "draft", "version": v, "time": _now()})
    st["items"] += new
    _save(s, st)
    return new


def check(s, verdicts, version=None):
    """The agent's check of its drafts: [{"id", "keep", "reason", "evidence"}]. Kept drafts become the scope's open
    advice (replacing the earlier open items there, except regions the user asked for); dropped ones are filed as
    filtered out with the reason. Returns what was done and the drafts still unchecked."""
    st = load(s)
    by_id = {x["id"]: x for x in st["items"] if x["status"] == "draft"}
    scopes = set()
    for vd in verdicts:
        x = by_id.get(str(vd.get("id", "")).strip())
        if x is None:
            raise ValueError(f"no draft {vd.get('id')!r} (drafts: {', '.join(by_id) or 'none'})")
        reason = str(vd.get("reason", "")).strip()
        if not reason:
            raise ValueError(f"{x['id']}: give the reason for keeping or dropping it")
        scopes.add(scope_of(x))
    for sc in scopes:  # a checked scope: its earlier open advice and filtered drafts make way for the new ones
        st["items"] = [x for x in st["items"] if not (scope_of(x) == sc and (x["status"] == "filtered" or
                                                       (x["status"] == "open" and not x.get("requested"))))]
    for vd in verdicts:
        x = by_id[str(vd["id"]).strip()]
        keep = bool(vd.get("keep"))
        x["status"] = "open" if keep else "filtered"
        x["check"] = {"keep": keep, "reason": str(vd["reason"]).strip(), "evidence": str(vd.get("evidence", "")).strip(),
                      "time": _now()}
    left = [x["id"] for x in st["items"] if x["status"] == "draft"]
    v = s.state(version)["version"]
    for sc in scopes:
        if not any(scope_of(x) == sc and x["status"] == "draft" for x in st["items"]):
            st["runs"].append({"scope": sc, "version": v, "labels_key": labels_key(s, version), "time": _now(),
                               "n_items": sum(scope_of(x) == sc and x["status"] == "open" and "check" in x for x in st["items"]),
                               "n_filtered": sum(scope_of(x) == sc and x["status"] == "filtered" for x in st["items"])})
    _save(s, st)
    return {"kept": [vd["id"] for vd in verdicts if vd.get("keep")], "filtered": [vd["id"] for vd in verdicts if not vd.get("keep")],
            "unchecked": left}


def items(s, status=None, scope=None):
    return [x | {"scope": scope_of(x)} for x in load(s)["items"]
            if (status is None or x["status"] == status) and (scope is None or scope_of(x) == scope)]


def get(s, aid):
    hit = [x for x in load(s)["items"] if x["id"] == aid]
    if not hit:
        raise ValueError(f"no advice {aid}")
    return hit[0]


def set_status(s, aid, status, note=""):
    if status not in ("open", "done", "dismissed"):
        raise ValueError("status must be open, done or dismissed")
    st = load(s)
    hit = [x for x in st["items"] if x["id"] == aid]
    if not hit:
        raise ValueError(f"no advice {aid}")
    hit[0]["status"] = status
    if note:
        hit[0]["note"] = note
    hit[0]["status_time"] = _now()
    _save(s, st)
    return hit[0]


def request(s, scope="annotation", version=None):
    """Note that the agent was asked to review a scope of this version."""
    st = load(s)
    st["requests"].append({"scope": scope, "version": s.state(version)["version"], "labels_key": labels_key(s, version),
                           "time": _now()})
    _save(s, st)


def status(s, scope="annotation", version=None):
    """A scope's review: whether it is current (annotation: for these labels), being prepared (asked, not finished),
    its drafts, open and filtered items, and its last run."""
    st = load(s)
    runs = [r for r in st["runs"] if r.get("scope", "annotation") == scope]
    reqs = [r for r in st["requests"] if r.get("scope", "annotation") == scope]
    mine = [x for x in st["items"] if scope_of(x) == scope]
    key = labels_key(s, version)
    last_req, last_run = (reqs[-1]["time"] if reqs else None), (runs[-1]["time"] if runs else None)
    return {"current": any(r["labels_key"] == key for r in runs),
            "preparing": bool(last_req) and (last_run is None or last_run < last_req),
            "n_drafts": sum(x["status"] == "draft" for x in mine), "n_open": sum(x["status"] == "open" for x in mine),
            "n_filtered": sum(x["status"] == "filtered" for x in mine), "last_run": runs[-1] if runs else None,
            "asked_at": last_req, "n_runs": len(runs)}


def signature(s):
    """Changes the viewer polls for."""
    st = load(s)
    return [len(st["runs"]), len(st["items"]), sum(x["status"] == "draft" for x in st["items"]), len(st["requests"])]


# ---------------- report ----------------
def _rfile(s):
    return s.dir / "report.md"


def write_report(s, markdown, version=None):
    text = str(markdown).strip()
    if len(text) < 40:
        raise ValueError("the report is empty")
    meta = {"version": s.state(version)["version"], "time": _now()}
    _rfile(s).write_text(text + "\n")
    json.dump(meta, open(s.dir / "report.json", "w"))
    return meta | {"path": str(_rfile(s))}


def report(s):
    p = _rfile(s)
    if not p.exists():
        return None
    meta = json.load(open(s.dir / "report.json")) if (s.dir / "report.json").exists() else {}
    return meta | {"markdown": p.read_text(), "path": str(p)}
