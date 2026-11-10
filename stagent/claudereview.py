"""Claude's check of the computed review candidates, so the Review panel shows only real problems.

When the user clicks Review annotation, the app asks the user's own Claude Code to look at every candidate area on
the H&E (show_area) and judge every marker check; Claude records a verdict per item with record_review (keep = a
real annotation problem, or a false alarm, each with a one-line reason). The Review panel lists only kept items;
items not checked yet are not listed while Claude is checking, and false alarms never appear (one collapsed line
lists them with Claude's reasons). Verdicts carry over to later versions: areas by lineage and tiles (their ids
change), marker checks by key. Nothing here changes the annotation.
"""
import json, time
from . import plausibility, markerqc
from .osutil import dump_json

OVERLAP = 0.5  # a candidate of a later version shares this share of its tiles with an area Claude judged


def _file(s):
    return s.dir / "claude_review.json"


def load(s):
    p = _file(s)
    return json.load(open(p)) if p.exists() else {"requests": [], "areas": [], "markers": {}}


def _save(s, st):
    dump_json(st, _file(s), indent=1)


def _area_verdict(a, rows):
    tiles = {tuple(t) for t in a["tile_set"]}
    hits = [r for r in rows if r["lineage"] == a["lineage"]
            and len(tiles & {tuple(t) for t in r["tile_set"]}) >= OVERLAP * len(tiles)]
    return hits[-1] if hits else None


def items(s, version=None):
    """The open candidates (not solved by the user) of a version, each with Claude's verdict or None."""
    st = load(s)
    areas = [(a, _area_verdict(a, st["areas"])) for a in plausibility.open_areas(s, version)]
    try:
        flags, unassigned = markerqc.open_checks(s, version)
    except ValueError:  # annotations made without a marker list have no marker checks
        flags, unassigned = [], None
    checks = [(f | {"key": markerqc.check_key(f)}, st["markers"].get(markerqc.check_key(f))) for f in flags]
    if unassigned:
        checks.append(({"key": "unassigned", "node": None, "marker": None, **unassigned}, st["markers"].get("unassigned")))
    return areas, checks


def status(s, version=None):
    """What the Review panel needs: kept / filtered-out / unchecked counts, and whether Claude was asked."""
    areas, checks = items(s, version)
    st = load(s)
    v = plausibility.check(s, version)["version"]
    asked = [r for r in st["requests"] if r["version"] == v]
    n_un = sum(1 for _, x in areas + checks if x is None)
    return {"asked": bool(asked), "asked_at": asked[-1]["time"] if asked else None, "unchecked": n_un,
            "checked": sum(1 for _, x in areas + checks if x is not None), "checking": bool(asked) and n_un > 0}


def request(s, version=None):
    """Note that Claude was asked to check this version's candidates."""
    st = load(s)
    v = plausibility.check(s, version)["version"]
    st["requests"].append({"version": v, "time": time.strftime("%Y-%m-%dT%H:%M:%S")})
    _save(s, st)
    return v


def record(s, areas=(), marker_checks=(), version=None):
    """Claude's verdicts: areas [{"id", "keep", "reason"}] (ids of this version), marker_checks [{"key", "keep",
    "reason"}]."""
    st = load(s)
    sc = plausibility.scan(s, version)
    by_id = {a["id"]: a for a in sc["areas"]}
    now = time.strftime("%Y-%m-%dT%H:%M:%S")
    done = []
    for x in areas:
        a = by_id.get(str(x.get("id", "")).strip().upper())
        if a is None:
            raise ValueError(f"no candidate area {x.get('id')!r} in version {sc['version']}")
        st["areas"].append({"lineage": a["lineage"], "tile_set": [list(t) for t in a["tile_set"]], "tiles": a["tiles"],
                            "area_id": a["id"], "keep": bool(x.get("keep")), "reason": str(x.get("reason", "")).strip(),
                            "version": sc["version"], "time": now})
        done.append(a["id"])
    for x in marker_checks:
        key = str(x.get("key", "")).strip()
        if not key:
            raise ValueError("a marker check verdict needs its key ('node|marker' or 'unassigned')")
        st["markers"][key] = {"keep": bool(x.get("keep")), "reason": str(x.get("reason", "")).strip(),
                              "version": sc["version"], "time": now}
        done.append(key)
    _save(s, st)
    return {"recorded": done, **status(s, version)}
