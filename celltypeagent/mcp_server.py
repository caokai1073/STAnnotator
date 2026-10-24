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
- The user usually has the Spatial Annotator app (or the web viewer) open: call use_viewer_session to load what they
  are looking at. After starting a new annotation call open_viewer: a running app switches to it by itself;
  otherwise give the user the link. The viewer shows the cell types on the H&E (or DAPI); the user draws and names
  regions there and can review conflicts. Changes made through these tools appear in the viewer automatically.
- A new gene panel needs its own marker list: draft one from get_panel_genes, validate it with check_marker_list,
  and save it with save_marker_list; the user loads it in the viewer's "Markers..." dialog and clicks Annotate (or
  start_session(markers_text=...) when the user asks you to run it).

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

Reviewing an annotation ("does this make sense?", "check my annotation", or "Ask Claude" in the viewer's Review panel)
- Use your knowledge of tissue biology, marker genes and histology; the tools give the evidence. Cover:
  1. Markers (review_markers): for every lineage and cell type, how often its listed markers are detected compared
     with its siblings (flags: weak, not specific, negative marker detected), and the genes that actually separate
     those cells here (top_genes). Ask: do the top genes fit the label (e.g. a "dendritic cell" group whose top genes
     are C1QC, LYVE1 is macrophages)? Are markers weak in this tissue or panel, with a better panel gene available?
     Do Unassigned cells (or an unassigned cluster in list_clusters) express a recognizable type missing from the
     list (e.g. CD177, CSF3R, S100A12: neutrophils)? Is the composition plausible for the tissue and site?
  2. Space (show_tissue_overview, check_plausibility, show_area): identify compartments on the H&E (mucosa with
     crypts, submucosa, muscularis propria, lymphoid aggregates, ulcers) and areas where the annotation contradicts
     them (e.g. epithelial cells in the muscle layer). Candidate areas come from the annotation alone; vascular
     smooth muscle, the muscularis mucosae, and deep crypts or fistulas in IBD are real. Lineages per drawn region
     show e.g. epithelial cells outside the user's "mucosa".
- Then give a short, prioritized list of concrete suggestions, each with its evidence and the action:
  - a revised marker list: propose_marker_list(text, summary) validates it and puts it in the viewer's Review panel,
    where the user opens it in "Markers..." and annotates again keeping regions and priors (in Claude Desktop:
    start_session(markers_text=..., carry_from_current=True) after the user agrees);
  - a local re-clustering with sharper subtype markers (preview_recluster);
  - areas: suggest_areas (tiles, title, reason, action, region_name) draws them on the user's map; once the user
    adopted or drew the region (list_area_suggestions, list_regions), apply the action after a go-ahead: usually an
    exclude prior (preview_constraint mode="exclude", then add_constraint) so each cell gets its next most likely
    type, or label_cells(region=..., only_types=...) when the user names the target type;
  - priors or cluster labels.
  Change nothing without the user's agreement.

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
def use_viewer_session() -> str:
    """Load the sample and annotation the user is looking at in the Spatial Annotator app (or the web viewer) and
    make it the active session. Call it first when the user refers to what they see."""
    from . import claudelink
    st = claudelink.viewer_state()
    if not st.get("session"):
        return _j({"error": "the viewer shows no annotation; ask the user to open one (or start_session)",
                   "viewer_running": claudelink.viewer_alive(st)})
    _active["session"] = Session(st["session"])
    return _j({"session_id": st["session"], "sample": st.get("sample"), "summary": _active["session"].summary()})


@mcp.tool()
def open_viewer(port: int = 8765) -> str:
    """Show the active session to the user. If the Spatial Annotator app (or a web viewer) is running, it switches to
    the session by itself; otherwise a web viewer is started and its link returned."""
    import subprocess, sys, time
    from . import claudelink, config
    st = claudelink.viewer_state()
    s = _active.get("session")
    if claudelink.viewer_alive(st):
        if s is not None:
            claudelink.request_show(s.sample_id, s.id)
        where = "the Spatial Annotator app" if st.get("mode") == "desktop" else f"http://localhost:{st['port']}"
        return _j({"shown_in": where, "session": s.id if s else None,
                   "note": "the viewer switches to this session within a few seconds"})
    url = f"http://localhost:{port}"
    if not claudelink.viewer_alive({"port": port}):
        log = open(config.HOME / "viewer.log", "a")
        subprocess.Popen(config.self_command("viewer", "--port", str(port)), stdout=log, stderr=log,
                         stdin=subprocess.DEVNULL, start_new_session=True, env=dict(os.environ))
        for _ in range(60):
            time.sleep(0.5)
            if claudelink.viewer_alive({"port": port}):
                break
        else:
            raise RuntimeError(f"the viewer did not start; see {config.HOME / 'viewer.log'}")
    if s is not None:
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
def start_session(sample_id: str, region_source: str = "none", markers_text: Optional[str] = None,
                  carry_from_current: bool = False) -> str:
    """Start a new annotation session on a sample. The data is treated as unlabeled: cell types come from a
    hierarchical marker gene list (markers_text, format as in check_marker_list; default: the project's draft list).

    region_source: 'none' (the user draws regions on the map) or 'pathologist' (start from the dataset annotation).
    carry_from_current: copy the regions, priors and corrections of the active session of the same sample (use it
    when re-annotating with a revised marker list). The first call to any tool that needs labels runs the annotation
    (about a minute). Give the user a new viewer link afterwards (open_viewer).
    """
    carry = _active["session"].id if carry_from_current and _active.get("session") else None
    _active["session"] = Session.create(sample_id, region_source, method="markers", markers_text=markers_text,
                                        carry_from=carry)
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
def review_markers(node: Optional[str] = None, version: Optional[int] = None) -> str:
    """Marker check of the current annotation. For every node of the marker tree (lineage, group, cell type): its
    cell count, each listed marker's detection share in its cells vs its siblings' cells with a flag (weak / not
    specific / negative marker detected), and the panel genes that best separate it from its siblings here; plus the
    genes of Unassigned cells. node: only that node and its subtypes."""
    from . import markerqc
    q = markerqc.check(_s(), version)
    if node:
        names = {node} | set(_s().type_groups.get(node, []))
        q = q | {"nodes": [n for n in q["nodes"] if n["node"] in names or n["parent"] in names],
                 "flags": [f for f in q["flags"] if f["node"] in names]}
    return _j(q)


@mcp.tool()
def propose_marker_list(text: str, summary: str) -> str:
    """Validate a revised marker list (full text, same format as the current one) and offer it to the user: it
    appears in the viewer's Review panel, from where the user opens it in "Markers..." and annotates again.
    summary: what changed and why, in a few lines. Returns the changed lines."""
    from . import markerqc
    rec = markerqc.propose(_s(), text, summary)
    return _j({k: rec[k] for k in ("diff", "warnings", "n_types", "status")})


@mcp.tool()
def check_plausibility(top: int = 12, version: Optional[int] = None) -> str:
    """Candidate areas where the annotation may contradict the tissue: grid tiles with cells of a sheet-forming type
    (epithelium, smooth muscle) that sit apart from their own lineage, merged into areas, with composition,
    suggested types and drawn regions; plus lineage counts per drawn region. Changes nothing."""
    from . import plausibility
    return _j(plausibility.summary(_s(), version, top))


@mcp.tool()
def show_tissue_overview(highlight: Optional[list[str]] = None, version: Optional[int] = None) -> list:
    """Whole sample: H&E with the tile grid (left) and cells colored by lineage with candidate areas (right).
    highlight: lineages or types to color (others gray), e.g. ["Epithelial"]."""
    from . import plausibility
    s = _s()
    png = plausibility.fig_overview(s, version, highlight)
    g = plausibility.scan(s, version)["grid"]
    return [_j({"path": png, "tile_um": g["tile_um"], "rows": f"{g['rows'][0]}-{g['rows'][-1]}", "columns": f"1-{g['nx']}"}),
            Image(path=png)]


@mcp.tool()
def show_area(area: str, highlight: Optional[list[str]] = None, version: Optional[int] = None) -> list:
    """Zoom on an area: a candidate id ("A3") or tiles ("F3:H6", "F3, G4:G6"). H&E (left) and cells colored by
    lineage on a faded H&E (right; a candidate's own lineage is highlighted); out-of-place cells ringed."""
    from . import plausibility
    png, a = plausibility.fig_area(_s(), area, version, highlight)
    return [_j({k: a[k] for k in plausibility.AREA_KEYS if k in a} | {"path": png}), Image(path=png)]


@mcp.tool()
def suggest_areas(areas: list[dict]) -> str:
    """Show area suggestions on the user's map (replaces your earlier open suggestions). Each item:
    {"area": "A3"} or {"tiles": "F3:H6"}, plus "title" (short), "reason" (what the H&E and the annotation show),
    "action" (what to do once it is a region, e.g. "exclude Epithelial: cells get their next most likely type"),
    "region_name" (proposed name, e.g. "muscle_1"). The user adopts a suggestion as a region or draws it."""
    from . import plausibility
    new = plausibility.suggest(_s(), areas)
    return _j([{k: x[k] for k in ("id", "title", "tiles", "region_name")} for x in new])


@mcp.tool()
def list_area_suggestions() -> str:
    """Area suggestions and their status: open, done (adopted as a region; region_name says which) or dismissed."""
    from . import plausibility
    keep = ("id", "title", "reason", "action", "region_name", "tiles", "status")
    return _j([{k: x.get(k) for k in keep} for x in plausibility.suggestions(_s())])


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
def save_marker_list(name: str, text: str, sample_id: Optional[str] = None) -> str:
    """Save a marker list (validated against the sample's panel) in the user's library, where the viewer's
    "Markers..." dialog lists it ("load a list..."). Use it for a new panel: draft a list from get_panel_genes,
    check it with check_marker_list, save it, and tell the user to load it and click Annotate."""
    from . import config, markers
    from .session import panel_genes
    import re
    p = markers.parse(text, panel_genes(sample_id))
    if not p["ok"]:
        raise vocab.VocabError("marker list has errors: " + "; ".join(p["errors"][:8]))
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name.strip()).strip("_") or "marker_list"
    config.USER_MARKERS.mkdir(parents=True, exist_ok=True)
    path = config.USER_MARKERS / f"{safe}.txt"
    path.write_text(text.strip() + "\n")
    return _j({"saved": str(path), "name": safe, "n_types": p["n_types"], "warnings": p["warnings"]})


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
