"""MCP server exposing the annotation engine as tools (stdio transport).

Run: python -m celltypeagent.mcp_server
The server keeps one active session; every tool that changes state creates a new version.
Tool results are compact JSON; figure tools also return the PNG so the client can display it.
"""
import json, os
from typing import Optional, Union
from mcp.server.fastmcp import FastMCP, Image
from . import vocab
from .session import Session, list_samples as _list_samples, list_sessions as _list_sessions

INSTRUCTIONS = """You help a biologist annotate cell types in imaging spatial transcriptomics (e.g. Xenium) on their own
computer. The tools compute every result; never invent cell types, regions, numbers or files.

Getting started
- list_samples. If the sample is not there, ask for: the h5ad file(s) (raw counts, cell coordinates in um), the
  registered H&E image, and the registration file (affine CSV or Xenium Explorer matrix); for an h5ad holding several
  sections also the obs column and value of the section. Then call import_sample with those paths.
- start_session(sample_id) runs a marker-based annotation (no labeled reference is used; the data is treated as
  unlabeled) with the default marker list or one the user gives (check it first with check_marker_list), or
  resume_session from list_sessions. The first annotation takes about a minute on a laptop.
- Call open_viewer and give the user the link: the viewer shows the cell types on the H&E; the user draws and names
  regions there and can review conflicts. Changes made through these tools appear in the viewer automatically.

Cell types and markers
- Cell types are the most specific entries of the marker tree; a type with subtypes is a lineage / group: it can be
  used in priors but is never a label. Use only names from get_vocabulary; ask when the user's words are ambiguous.
- list_clusters shows the marker-step clusters, their marker scores and the Leiden resolution chosen per level.
- For a mixed cluster, a type to split further, or a cluster the user keeps correcting: preview_recluster, report
  the marker-label changes and new clusters, and call recluster only after the user agrees. Report changes inside
  the re-clustered cells and elsewhere separately. New or removed cell types need a revised marker list and a new
  annotation (start_session with markers_text).

Priors and corrections
- A general prior ("epithelial cells are only in the mucosa", "Tfh cells sit next to B cells"): map it to a
  constraint (kind region / depth / neighborhood), call preview_constraint, explain your reading and the predicted
  changes in one or two sentences, and wait for an explicit go-ahead before add_constraint (text = the user's words).
  Region priors only act on cells inside drawn regions; say which regions exist (list_regions).
- A correction ("the smooth muscle in mucosa_2 is muscularis mucosae", "root/c27 is neurons"): label_cells
  (region=..., only_types=... for part of a region) or label_cluster; the instruction is the go-ahead. Report how many
  cells you labeled and how many others changed because the model was updated; suggest a region or depth prior
  when the types differ mainly by position. You cannot draw regions: ask the user to draw and name them in the viewer.
- Conflicts (conflict_summary, show_conflicts) are cells whose expression contradicts a prior; describe H&E crops
  cautiously and record verdicts only when the user states them.

Plausibility check ("does the spatial distribution make sense?", or the viewer's Check plausibility button)
- check_plausibility lists groups of isolated cells of cohesive types (types that form sheets, e.g. epithelium,
  smooth muscle), with the lineage around them, a suggested type, confidence, and lineages per drawn region.
  "likely mislabel" groups were in another lineage about half the time on test data; confident isolated cells
  much less often, so judge each group, never relabel everything.
- Judge with biology and the H&E (show_plausibility_groups, up to 4 groups per call): epithelium in submucosa or
  muscle is suspicious in normal tissue but can be real in IBD (deep crypts, fistulas, ulcer edges); smooth muscle
  cells around vessels (vascular smooth muscle) and the thin muscularis mucosae are real; B cells and naive T cells
  outside follicles are normal. Also use lineages_by_region with the user's region names.
- Report three lists: likely mislabels (with the suggested type), unusual but possibly real (why), and unsure (for
  the user or a pathologist to look at in the viewer). Relabel only after the user agrees:
  label_plausibility_groups (each group gets its suggested type unless you give cell_type).

Every change is a version (history, diff, revert). export writes a CSV with one row per cell. Reply briefly, in the
language of the user's message, keeping cell type and region names as they are."""


mcp = FastMCP("celltypeagent", instructions=INSTRUCTIONS)
_active: dict = {}


def _s() -> Session:
    if "session" not in _active:
        raise ValueError("no active session; call start_session or resume_session first")
    return _active["session"]


def _j(obj) -> str:
    return json.dumps(obj, default=str, ensure_ascii=False)


def _records(df, top=None):
    df = df if top is None else df.head(top)
    return df.to_dict("records")


def _spec(text, types, mode, weight, kind, regions, max_depth_um, neighbor_types, min_frac):
    return {"text": text, "types": types, "mode": mode, "weight": weight, "kind": kind, "regions": regions,
            "max_depth_um": max_depth_um, "neighbor_types": neighbor_types, "min_frac": min_frac}


# ---------------- data and viewer ----------------
@mcp.tool()
def import_sample(name: str, h5ad: list[str], he: str, registration: str, section_column: Optional[str] = None,
                  section: Optional[str] = None, region_column: Optional[str] = None, counts_layer: Optional[str] = None,
                  description: str = "") -> str:
    """Import one tissue section from files on this computer: h5ad file(s) with raw counts and cell coordinates in
    um (several files are merged), the registered H&E image (pyramidal TIFF / OME-TIFF), and the registration file
    (affine CSV a11,a12,a21,a22,tx,ty from cell um to H&E pixels, or a Xenium Explorer 3x3 matrix). For an h5ad with
    several sections give section_column and section. region_column optionally imports region annotations.
    Cell type annotations in the h5ad are ignored. Takes from seconds to a few minutes for large files."""
    from .importer import import_sample as imp
    log = []
    meta = imp(name, h5ad, he, registration, section_column, section, None, None, counts_layer, region_column,
               description=description, progress=log.append)
    return _j({"imported": meta["sample_id"], "n_cells": meta["n_cells"], "n_genes": meta["n_genes"],
               "he_um_per_px": round(meta["he_um_per_px"], 4), "log": log})


@mcp.tool()
def open_viewer(port: int = 8765) -> str:
    """Start the local web viewer if it is not running and return its link (opens the active session if any).
    The user draws regions, looks at the H&E and reviews conflicts there."""
    import subprocess, sys, time, urllib.request
    from . import config
    url = f"http://localhost:{port}"
    def up():
        try:
            urllib.request.urlopen(f"http://127.0.0.1:{port}/api/samples", timeout=2)
            return True
        except Exception:
            return False
    if not up():
        log = open(config.HOME / "viewer.log", "a")
        subprocess.Popen([sys.executable, "-m", "celltypeagent.webapp", "--port", str(port)], stdout=log, stderr=log,
                         stdin=subprocess.DEVNULL, start_new_session=True, env=dict(os.environ))
        for _ in range(60):
            time.sleep(0.5)
            if up():
                break
        else:
            raise RuntimeError(f"the viewer did not start; see {config.HOME / 'viewer.log'}")
    if "session" in _active:
        s = _active["session"]
        url += f"/?sample={s.sample_id}&session={s.id}"
    return _j({"url": url, "note": "open the link in a browser; changes made here show up there within a few seconds"})


# ---------------- sessions ----------------
@mcp.tool()
def list_samples() -> str:
    """List available samples (sections) with disease, location, cell count, and whether pathologist regions / H&E exist."""
    return _j(_records(_list_samples()))


@mcp.tool()
def list_sessions() -> str:
    """List saved annotation sessions that can be resumed."""
    return _j(_records(_list_sessions()))


@mcp.tool()
def start_session(sample_id: str, region_source: str = "none", markers_text: Optional[str] = None) -> str:
    """Start a new annotation session on a sample. The data is treated as unlabeled: cell types come from a
    hierarchical marker gene list (markers_text, format as in check_marker_list; default: the project's draft list).

    region_source: 'none' (the user draws regions on the map) or 'pathologist' (start from the dataset annotation).
    Version 1 has no constraints. The first call to any tool that needs labels runs the annotation (10-20 s).
    """
    _active["session"] = Session.create(sample_id, region_source, method="markers", markers_text=markers_text)
    return _j(_active["session"].summary())


@mcp.tool()
def resume_session(session_id: str) -> str:
    """Resume a saved session (latest version)."""
    _active["session"] = Session(session_id)
    return _j(_active["session"].summary())


@mcp.tool()
def get_vocabulary() -> str:
    """Exact cell type names, type groups, region names, constraint modes, strength presets and verdict options.
    Always map the user's words onto these names; ask the user when a mapping is ambiguous."""
    return _j(_s().vocabulary())


@mcp.tool()
def session_summary(version: Optional[int] = None) -> str:
    """Summary of a version (default latest): region composition, constraints, number of changed cells and conflicts."""
    return _j(_s().summary(version))


# ---------------- constraints ----------------
@mcp.tool()
def preview_constraint(text: str, types: list[str], mode: str, weight: Union[str, float] = "strong",
                       kind: str = "region", regions: Optional[list[str]] = None,
                       max_depth_um: Optional[float] = None, neighbor_types: Optional[list[str]] = None,
                       min_frac: Optional[float] = None) -> str:
    """Validate a constraint and show its effect WITHOUT committing it. Show the interpretation and the
    predicted label changes to the user and get explicit confirmation before calling add_constraint.

    text: the user's own words, verbatim. types: cell type or group names from get_vocabulary.
    kind: 'region' (needs regions), 'depth' (needs max_depth_um: distance to the epithelium / crypt core),
          'neighborhood' (needs neighbor_types and min_frac, the share of the 15 nearest neighbors).
    mode: 'only' (types unlikely where the condition is false), 'exclude' (unlikely where true),
          'enrich' (more likely where true).
    weight: 'strong' | 'moderate' | 'weak' or a number (only/exclude: 0-1, enrich: >1).
    """
    return _j(_s().preview_constraint(_spec(text, types, mode, weight, kind, regions, max_depth_um, neighbor_types, min_frac)))


@mcp.tool()
def add_constraint(text: str, types: list[str], mode: str, weight: Union[str, float] = "strong",
                   kind: str = "region", regions: Optional[list[str]] = None,
                   max_depth_um: Optional[float] = None, neighbor_types: Optional[list[str]] = None,
                   min_frac: Optional[float] = None, author: Optional[str] = None,
                   message: Optional[str] = None) -> str:
    """Commit a constraint the user has confirmed (same arguments as preview_constraint, plus the author's
    name if known). Creates a new version."""
    s = _s()
    spec = _spec(text, types, mode, weight, kind, regions, max_depth_um, neighbor_types, min_frac) | {"author": author}
    v = s.add_constraint(spec, message)
    return _j({"new_version": v, "diff": s.diff(v - 1, v), "n_conflicts": s.summary()["n_conflicts"]})


@mcp.tool()
def update_constraint(constraint_id: str, weight: Union[str, float, None] = None, types: Optional[list[str]] = None,
                      regions: Optional[list[str]] = None, mode: Optional[str] = None,
                      max_depth_um: Optional[float] = None, neighbor_types: Optional[list[str]] = None,
                      min_frac: Optional[float] = None, text: Optional[str] = None,
                      message: Optional[str] = None) -> str:
    """Change fields of an existing constraint (e.g. weaken it or change its threshold). Creates a new version."""
    changes = {k: v for k, v in dict(weight=weight, types=types, regions=regions, mode=mode, max_depth_um=max_depth_um,
                                     neighbor_types=neighbor_types, min_frac=min_frac, text=text).items() if v is not None}
    s = _s()
    v = s.update_constraint(constraint_id, changes, message)
    return _j({"new_version": v, "diff": s.diff(v - 1, v)})


@mcp.tool()
def remove_constraint(constraint_id: str, message: Optional[str] = None) -> str:
    """Remove a constraint. Creates a new version."""
    s = _s()
    v = s.remove_constraint(constraint_id, message)
    return _j({"new_version": v, "diff": s.diff(v - 1, v)})


# ---------------- regions ----------------
@mcp.tool()
def list_regions(version: Optional[int] = None) -> str:
    """Regions of the session. Regions are drawn and named by the user on the map (the agent cannot draw them);
    cells outside every region are 'unknown' and no region prior applies to them."""
    info = _s().regions_info(version)
    return _j({"source": info["source"], "regions": info["regions"], "n_unknown": info["n_unknown"],
               "polygons": [{k: p[k] for k in ("id", "name", "n_cells_drawn")} for p in info["polygons"]]})


@mcp.tool()
def label_cells(region: str, cell_type: str, only_types: Optional[list[str]] = None, text: str = "",
                author: Optional[str] = None) -> str:
    """Hand-label the cells of a user-drawn region as cell_type. only_types restricts it to cells currently labeled
    as those types (names or groups), e.g. relabel only the smooth muscle in that region. The labeled cells are
    locked and the expression model is updated from them, so similar cells elsewhere may change too.
    Creates a new version. text: the user's words, verbatim."""
    s = _s()
    v = s.label_cells(cell_type, region=region, only_types=only_types, text=text, author=author)
    corr = s.state(v)["corrections"][-1]
    d = s.diff(v - 1, v)
    return _j({"new_version": v, "correction_id": corr["id"], "n_cells_labeled": len(corr["cells"]),
               "n_label_changes_total": d["n_label_changes"],
               "n_other_cells_changed": max(0, d["n_label_changes"] - len(corr["cells"])),
               "top_transitions": d["top_transitions"][:8]})


@mcp.tool()
def set_region_source(source: str, message: Optional[str] = None) -> str:
    """Start regions from nothing ('none') or from the dataset's pathologist annotation ('pathologist');
    regions the user drew are kept on top."""
    s = _s()
    v = s.set_region_source(source, message)
    return _j({"new_version": v, "diff": s.diff(v - 1, v)})


@mcp.tool()
def show_region_map(conflict_group: Optional[str] = None, version: Optional[int] = None) -> list:
    """Spatial map colored by region. conflict_group: '' highlights all conflicts, a group_id highlights one group."""
    p = _s().fig_region_map(version, conflict_group)
    return [_j({"path": p}), Image(path=p)]


# ---------------- conflicts and review ----------------
@mcp.tool()
def conflict_summary(by: list[str] = ["constraint", "region"], top: int = 15, version: Optional[int] = None) -> str:
    """Conflicts (confident expression-based labels that violate a constraint) grouped by any of
    constraint, region, category, status, base_label. Each row has a group_id for show_conflicts."""
    s = _s()
    t = s.conflict_summary(tuple(by), version)
    return _j({"n_conflicts": int(t.n.sum()) if len(t) else 0, "n_groups": len(t), "groups": _records(t, top)})


@mcp.tool()
def show_conflicts(group_id: str, n: int = 6, version: Optional[int] = None) -> list:
    """H&E crops of sampled conflict cells in a group (raw H&E on top, lineage overlay below), with each cell's
    expression-based label, current label and reference label. Returns cell ids for record_verdicts."""
    png, q = _s().fig_conflicts(group_id, n, version)
    keep = ["cell", "constraint", "region", "category", "base_label", "base_conf", "label",
            "status", "already_reviewed"]
    return [_j({"group_id": group_id, "cells": [{k: r[k] for k in keep} for r in q], "path": png}), Image(path=png)]


@mcp.tool()
def check_plausibility(top: int = 15, version: Optional[int] = None) -> str:
    """Find cells whose label does not fit where they lie: isolated cells of cohesive cell types (types that form
    sheets), grouped in space and ranked (likely mislabels first). Each group lists its cell types, the lineage
    around it, a suggested type (best type outside its lineage by expression, weighted by the lineages around),
    confidence, regions and distance to the epithelium. Changes nothing."""
    from . import plausibility
    return _j(plausibility.summary(_s(), version, top))


@mcp.tool()
def show_plausibility_groups(group_ids: list[str], version: Optional[int] = None) -> list:
    """H&E crops of up to 4 plausibility groups (raw H&E on top, cells colored by lineage below; the group's cells
    ringed in cyan)."""
    png, gs = _s().fig_plausibility(group_ids[:4], version)
    keep = ("id", "n_cells", "lineage", "types", "surrounded_by", "suggested_type", "median_confidence", "category",
            "regions", "median_depth_um")
    return [_j({"groups": [{k: g[k] for k in keep} for g in gs], "path": png}), Image(path=png)]


@mcp.tool()
def label_plausibility_groups(group_ids: list[str], cell_type: Optional[str] = None, text: str = "") -> str:
    """Hand-label the cells of plausibility groups of the latest version: each group gets its suggested type, or
    cell_type for all of them. One new version (the cells are locked and the model is updated). Only after the user
    agreed. text: the user's words."""
    from . import plausibility
    s = _s()
    v = plausibility.label_groups(s, group_ids, cell_type, text=text, author="agent")
    d = s.diff(v - 1, v)
    n = sum(len(c["cells"]) for c in s.state(v)["corrections"][-len(group_ids):])
    return _j({"new_version": v, "n_cells_labeled": n, "n_label_changes_total": d["n_label_changes"],
               "n_other_cells_changed": max(0, d["n_label_changes"] - n), "top_transitions": d["top_transitions"][:8]})


@mcp.tool()
def record_verdicts(cells: list[int], verdict: str, note: str = "", reviewer: str = "user") -> str:
    """Record the reviewer's verdict for conflict cells. Only record what the user (or pathologist) decided.
    verdict: constraint_right | reference_right | region_wrong | artifact | new_state | unclear."""
    n = _s().record_verdicts(cells, verdict, note, reviewer)
    return _j({"recorded": n, "verdict": verdict, "meaning": vocab.VERDICTS[verdict]})


@mcp.tool()
def review_progress() -> str:
    """Latest verdict per (cell, reviewer), counted by verdict, reviewer and conflict group (from the agent
    and the review page). Use it to suggest weakening or narrowing constraints the reviewers contradicted."""
    v = _s().verdicts(latest=True)
    if v.empty:
        return _j({"n": 0})
    by_group = v.assign(group=v.group.fillna("unknown")).groupby(["group", "verdict"]).size().unstack(fill_value=0)
    return _j({"n": len(v), "by_verdict": v.verdict.value_counts().to_dict(),
               "by_reviewer": v.reviewer.value_counts().to_dict(),
               "by_group": by_group.reset_index().to_dict("records")})


@mcp.tool()
def remove_correction(correction_id: str) -> str:
    """Undo one hand correction (e.g. 'L2'); the expression model is recomputed without it. Creates a new version."""
    s = _s()
    v = s.remove_correction(correction_id)
    return _j({"new_version": v, "diff": s.diff(v - 1, v)})


@mcp.tool()
def get_panel_genes(sample_id: Optional[str] = None) -> str:
    """All genes of a sample's panel (default: the active session's sample). Marker lists may use only these genes."""
    from .session import panel_genes
    sid = sample_id or (_active["session"].sample_id if "session" in _active else None)
    genes = panel_genes(sid)
    return _j({"n_genes": len(genes), "genes": genes})


@mcp.tool()
def check_marker_list(text: str, sample_id: Optional[str] = None) -> str:
    """Validate a marker list (format: 'Name: GENE+ GENE- A|B+' per line, children indented two spaces) against a
    sample's panel (default: the active session's sample). Returns errors (with suggestions for genes not in the
    panel), warnings, and the parsed outline."""
    from .session import panel_genes
    from . import markers
    sid = sample_id or (_active["session"].sample_id if "session" in _active else None)
    p = markers.parse(text, panel_genes(sid))
    return _j({k: p[k] for k in ("ok", "errors", "warnings", "n_types", "depth", "canonical")})


@mcp.tool()
def list_clusters(parent: Optional[str] = None, only_unassigned: bool = False, top: int = 30) -> str:
    """Clusters of the marker-annotation step (marker sessions only): id, parent type, size, assigned type or
    None, marker score per candidate type, Leiden resolution, and current labels inside; plus the resolution chosen
    at every level. parent filters to one level (e.g. 'Epithelial'; 'root' = top level)."""
    rec = _s().clusters()
    if parent:
        rec = [r for r in rec if (r["parent"] or "root") == parent]
    if only_unassigned:
        rec = [r for r in rec if r["assigned"] is None]
    rec = sorted(rec, key=lambda r: -r["n_cells"])[:top]
    return _j({"n": len(rec), "clusters": rec,
               "levels": [{k: n[k] for k in ("node", "n_cells", "resolution", "tried", "assigned_frac")}
                          for n in _s().marker_levels()]})


@mcp.tool()
def label_cluster(cluster_id: str, cell_type: str, text: str = "") -> str:
    """Label all cells of a marker-step cluster (e.g. 'root/c27') as cell_type; works like label_cells (cells locked,
    expression model updated). Creates a new version."""
    s = _s()
    mres = s.marker_result(s.state())
    if cluster_id not in mres["members"]:
        raise ValueError(f"no cluster {cluster_id}; see list_clusters")
    cells = s.cells.cell.to_numpy()[mres["members"][cluster_id]].tolist()
    v = s.label_cells(cell_type, cells=cells, text=text or f"cluster {cluster_id}")
    d = s.diff(v - 1, v)
    return _j({"new_version": v, "n_cells_labeled": len(cells), "n_label_changes_total": d["n_label_changes"],
               "top_transitions": d["top_transitions"][:8]})


@mcp.tool()
def preview_recluster(target: str, resolution: Union[str, float] = "auto", marker_edits: Optional[dict] = None) -> str:
    """Preview re-clustering part of the marker step WITHOUT committing: target is 'root', a cell type (its cells are
    split again into its subtypes) or a cluster id from list_clusters (only that cluster is re-clustered).
    resolution: 'auto' (tries 0.5/1/2/4 and keeps the one the markers resolve best) or a number (higher = finer).
    marker_edits: optional {type: "GENE+ GENE- A|B+"} replacing the markers of existing subtypes for this re-run.
    Show the user the number of cells whose marker label changes and the new clusters, then ask for a go-ahead."""
    return _j(_s().preview_recluster(target, resolution, marker_edits))


@mcp.tool()
def recluster(target: str, resolution: Union[str, float] = "auto", marker_edits: Optional[dict] = None,
              text: str = "", author: Optional[str] = None) -> str:
    """Commit a re-clustering the user agreed to (same arguments as preview_recluster; text = the user's words).
    Regions, priors and corrections are kept; the expression model is retrained on the new marker labels.
    Reports changes inside the re-clustered cells separately from changes elsewhere (model update)."""
    import numpy as np
    s = _s()
    v = s.recluster(target, resolution, marker_edits, text, author)
    mres = s.marker_result(s.state(v))
    op = s.state(v)["marker_step"]["ops"][-1]
    rows = np.unique(np.concatenate([m for k, m in mres["members"].items() if k.startswith(op["id"] + ":")]))
    a, b = s.result(v - 1)["pred"], s.result(v)["pred"]
    ch = a != b
    inside = np.zeros(len(a), bool); inside[rows] = True
    d = s.diff(v - 1, v)
    return _j({"new_version": v, "n_cells_reclustered": int(len(rows)), "n_changed_inside": int((ch & inside).sum()),
               "n_changed_elsewhere": int((ch & ~inside).sum()), "top_transitions": d["top_transitions"][:8],
               "levels": [n for n in mres["nodes"] if n.get("prefix") == op["id"] + ":"]})


@mcp.tool()
def list_rule_suggestions(status: Optional[str] = "open") -> str:
    """Rules that reviewers wrote in plain words on the review page (status: open | accepted | rejected | None for all),
    with the cell and conflict group they were looking at. Turn each into a constraint with the user, then resolve it."""
    return _j(_s().suggestions(status))


@mcp.tool()
def resolve_rule_suggestion(suggestion_id: str, status: str, note: str = "") -> str:
    """Mark a reviewer's rule suggestion as accepted or rejected (or reopen it), with a short note on what was done."""
    _s().resolve_suggestion(suggestion_id, status, note)
    return _j({"id": suggestion_id, "status": status})


# ---------------- history ----------------
@mcp.tool()
def history() -> str:
    """All versions of the session with the action and message of each."""
    return _j(_s().history())


@mcp.tool()
def diff(from_version: int, to_version: int) -> str:
    """Label and region changes between two versions."""
    return _j(_s().diff(from_version, to_version))


@mcp.tool()
def revert(version: int, message: Optional[str] = None) -> str:
    """Restore the regions and constraints of an earlier version as a new version (history is never deleted)."""
    s = _s()
    v = s.revert(version, message)
    return _j({"new_version": v, "summary": s.summary()})


@mcp.tool()
def export(version: Optional[int] = None) -> str:
    """Write cells.csv (cell_id, region, cell_type, lineage, type levels, confidence, ...), conflicts.csv, state.json
    and CHANGELOG.md for a version; returns the folder."""
    return _j({"folder": _s().export(version)})


def main():
    mcp.run()


if __name__ == "__main__":
    main()
