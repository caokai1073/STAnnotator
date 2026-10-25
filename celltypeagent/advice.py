"""Claude's high-level advice on an annotation, and the annotation report (stored per session).

Beyond single areas and marker checks, Claude looks at the annotation as a whole (evidence from insights.py) and
gives a few pieces of advice, each in one category:
  compartments   tissue compartments seen on the image and in the cells (mucosa, muscularis, lymphoid follicles):
                 proposed as regions with their outlines, often with a prior to use them
  taxonomy       types to merge, split or add (a revised marker list, a re-clustering)
  organization   spatial organization (gradients, neighbors) that is expected or contradicts the labels; priors
  context        what the sample is (site, disease) and what that implies for the annotation
  qc             technical problems: registration offset, low-quality areas, spillover, segmentation
Each item can carry one action the user takes in the viewer with one click:
  region         {"kind": "region", "name", "compartment"}: add the compartment outline (from preview_compartment)
                 as a region; an optional "prior" (constraint spec) is offered once the region exists
  prior          {"kind": "prior", "prior": spec}: preview, then add the prior
  marker_list    {"kind": "marker_list"}: open the marker list Claude proposed (propose_marker_list)
  registration   {"kind": "registration"}: correct the cells' place on the image (insights.registration)
  claude         {"kind": "claude", "request": "..."}: ask Claude Code to do it (e.g. a re-clustering), as a request
Advice is advice: nothing changes until the user clicks. A new set of advice replaces the open items; done and
dismissed ones are kept. The report is Markdown Claude writes (write_report); the viewer shows it and saves it.
"""
import hashlib, json, re, time
import numpy as np
from . import constraints as cmod

CATEGORIES = {"compartments": "Compartments", "taxonomy": "Cell types", "organization": "Spatial organization",
              "context": "Sample context", "qc": "Technical"}
ACTIONS = ("region", "prior", "marker_list", "registration", "claude")
MAX_ITEMS = 12


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
def _check_prior(s, spec, new_region=None):
    """A prior spec as the user will add it; regions may include the region the same item creates."""
    spec = {k: v for k, v in dict(spec).items() if v is not None}
    names = s.region_names() + ([new_region] if new_region else [])
    c = cmod.normalize(dict(spec, id="C?"), s.types, s.type_groups, names)
    rule = cmod.describe(c).split(": ", 1)[1]
    return {k: spec[k] for k in cmod.SPEC_FIELDS if k in spec and k != "id"}, \
        (f'{spec["text"]} ({rule.split("] ", 1)[-1]})' if spec.get("text") else rule)


def add(s, items, version=None):
    """Replace the open advice with items [{"category", "title", "text", "action"?}]; returns the stored items."""
    st = load(s)
    if len(items) > MAX_ITEMS:
        raise ValueError(f"at most {MAX_ITEMS} pieces of advice: keep the ones that matter")
    keep = [x for x in st["items"] if x["status"] != "open"]
    used = {x["id"] for x in keep}
    v = s.state(version)["version"]
    new = []
    for it in items:
        cat = str(it.get("category", "")).strip().lower()
        if cat not in CATEGORIES:
            raise ValueError(f"category must be one of {list(CATEGORIES)}, got {cat!r}")
        title, text = str(it.get("title", "")).strip(), str(it.get("text", "")).strip()
        if not title or not text:
            raise ValueError("every piece of advice needs a title and a text")
        act = dict(it.get("action") or {})
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
            reg = insights.registration(s)
            fit = reg.get("affine_fit")
            if not fit:
                raise ValueError("no registration fit for this sample")
            out |= {"correction_px": fit["correction_px"], "max_offset_um": fit["max_offset_um"]}
        elif kind == "claude":
            if not str(act.get("request", "")).strip():
                raise ValueError("a claude action needs 'request' (what Claude Code should do when the user asks)")
            out["request"] = str(act["request"]).strip()
        k = 1
        while f"H{k}" in used:
            k += 1
        used.add(f"H{k}")
        new.append({"id": f"H{k}", "category": cat, "title": title, "text": text, "action": out, "status": "open",
                    "version": v, "time": _now()})
    st["items"] = keep + new
    st["runs"].append({"version": v, "labels_key": labels_key(s, version), "time": _now(), "n_items": len(new)})
    _save(s, st)
    return new


def items(s, status=None):
    return [x for x in load(s)["items"] if status is None or x["status"] == status]


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


def request(s, version=None):
    """Note that Claude was asked for advice on this version's labels."""
    st = load(s)
    st["requests"].append({"version": s.state(version)["version"], "labels_key": labels_key(s, version), "time": _now()})
    _save(s, st)


def status(s, version=None):
    """Whether advice for these labels exists, was asked for, or is being prepared."""
    st = load(s)
    key = labels_key(s, version)
    done = any(r["labels_key"] == key for r in st["runs"])
    asked = any(r["labels_key"] == key for r in st["requests"])
    return {"current": done, "preparing": asked and not done, "n_open": sum(x["status"] == "open" for x in st["items"]),
            "last_run": st["runs"][-1] if st["runs"] else None, "n_runs": len(st["runs"])}


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
