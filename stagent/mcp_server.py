"""MCP server exposing the annotation engine as tools (stdio transport).

Run: python -m stagent.mcp_server
The server keeps one active session; every tool that changes state creates a new version.
Tool results are compact JSON; figure tools also return the PNG so the client can display it.
"""
import json, os, re
from typing import Optional, Union
import numpy as np
from mcp.server.fastmcp import FastMCP, Image
from . import vocab, render
from .session import Session, list_samples as _list_samples, list_sessions as _list_sessions

INSTRUCTIONS = """You help a biologist annotate cell types in imaging spatial transcriptomics (e.g. Xenium) on their own
computer. The tools compute every result; never invent cell types, regions, numbers or files.

Getting started
- list_samples. If the sample is not there: CosMx or MERFISH (MERSCOPE) output goes in with import_cosmx_merfish (its
  flat files or folders, optionally the user's H&E and registration); a Visium HD / Visium Space Ranger output with
  import_visium (cells from cell outlines when there are, the 2 um bins split by area; else bins at 8 um; no
  registration needed); a Stereo-seq SAW output with import_stereoseq (cells from its cell bins, else 8 um bins). Otherwise ask for: the h5ad file(s) or an R file (.rds / .rda with a Seurat or
  SpatialExperiment object; read without R) (raw counts, cell coordinates in um), the
  registered H&E image, and the registration file (affine CSV or Xenium Explorer matrix); for an h5ad holding several
  sections also the obs column and value of the section. Then call import_sample with those paths.
- start_session(sample_id, markers_text) runs a marker-based annotation (no labeled reference is used; the data is
  treated as unlabeled). There is no default marker list: use the one the user gives (text or CSV with columns
  cell_type, parent, positive_markers, negative_markers; check it first with check_marker_list) or draft one. Or,
  with a labeled single-cell reference h5ad of the tissue, start_session(sample_id, reference_h5ad=..., label_column=...)
  (inspect_reference first: shared genes, label columns). Or load an existing annotation (CSV, one row per cell):
  start_session(sample_id, annotation_csv=..., id_column=..., label_column=...) after inspect_annotation. Or
  resume_session from list_sessions. The review window of start_session shows how long the annotation should take
  on this computer (about a minute with the default approximate neighbor search; with knn="exact", under a minute
  with a GPU, minutes for large samples without one).
- The user usually has the STAnnotator app (or the web viewer) open: call use_viewer_session to load what they
  are looking at. After starting a new annotation call open_viewer: a running app switches to it by itself;
  otherwise give the user the link. The viewer shows the cell types on the H&E (or DAPI); the user draws and names
  regions there and can review conflicts. Changes made through these tools appear in the viewer automatically.
- A sample not annotated yet opens in the app as its placeholder (method "none": every cell "Not annotated"):
  regions, check_registration_map, the registration correction and fit_piece / move_piece work on it; tools about
  cell types do not. A new annotation of the sample (start_session) starts with the placeholder's regions.
- A sample of an H&E only (list_samples: image_only; import_he_only) has no cells: the user (e.g. a pathologist)
  draws and names regions on it. Asked to outline a compartment there (e.g. "draw the mucosa"), do it from the image
  and the tissue's priors: tissue_priors(organ) for how each compartment looks on the H&E, segment_he for the
  classes of the tissue (hematoxylin, eosin, texture, distance to the edge), then propose_he_region with the classes
  that match (within a box for one section), or with points read on show_he; say it is a draft from the image for the
  user to check and edit. The tools about cells do not apply. Its regions go to the data registered to the same
  image: the app's Import dialog ("Regions from"), or export_regions and the regions file of the import.
- A new gene panel needs its own marker list: start from organ_marker_draft (built-in lists of common organs pruned
  to the panel; search the literature for types it lacks), or draft one from get_panel_genes; validate it with check_marker_list,
  and save it with save_marker_list; the user loads it in the viewer's "Markers..." dialog and clicks Annotate (or
  start_session(markers_text=...) when the user asks you to run it).

Changes to the data are reviewed by the user
- Every tool that changes data (import_sample, import_cosmx_merfish, import_visium, import_stereoseq, import_he_only, start_session, add / update / remove_constraint, propose_he_region, label_cells,
  label_cluster, remove_correction, recluster, revert, apply_alignment, run_niches, name_niches, record_verdicts,
  resolve_rule_suggestion, assign_region_cells, propose_region, move_piece, undo_piece_move) only proposes the change: the app shows it in a review window with what it does, and it
  is applied only when the user clicks Accept there (Dismiss drops it). The tool waits a little for the decision and
  returns the outcome (status "accepted by the user and applied" with the result, "dismissed by the user", or
  "waiting for the user's review"). While waiting, tell the user in one or two sentences what you proposed and why;
  read the outcome later with pending_changes. The review window is the user's go-ahead: do not ask for it twice in
  the chat, and never propose a dismissed change again unchanged.

Cell types and markers
- Cell types are the most specific entries of the marker tree; a type with subtypes is a lineage / group: it can be
  used in priors but is never a label. Use only names from get_vocabulary; ask when the user's words are ambiguous.
- list_clusters shows the marker-step clusters, their marker scores and the Leiden resolution chosen per level.
- For a mixed cluster, a type to split further, or a cluster the user keeps correcting: preview_recluster, report
  the marker-label changes and new clusters, then recluster (the user accepts it in the app). Report changes inside
  the re-clustered cells and elsewhere separately. New or removed cell types need a revised marker list and a new
  annotation (start_session with markers_text).

Priors and corrections
- A general prior ("epithelial cells are only in the mucosa", "Tfh cells sit next to B cells"): map it to a
  constraint (kind region / depth / neighborhood), call preview_constraint, explain your reading and the predicted
  changes in one or two sentences, then add_constraint (text = the user's words): the user accepts it in the app.
  Region priors only act on cells inside drawn regions; say which regions exist (list_regions).
- A correction ("the smooth muscle in mucosa_2 is muscularis mucosae", "root/c27 is neurons"): label_cells
  (region=..., only_types=... for part of a region) or label_cluster (the user accepts it in the app). Report how many
  cells you labeled and how many others changed because the model was updated; suggest a region or depth prior
  when the types differ mainly by position.
- Regions are the user's. When the user asks you to draw or outline one ("draw the mucosa", "outline the follicles",
  "make a region of tiles F3:H6"): compute it with preview_compartment (tissue a lineage forms) or take tiles / a
  candidate area, check the figure against the H&E, then propose_region: the user sees the outline on the H&E in the
  app's review window and Accept adds the region; the user can then edit its points.
  Never add regions the user did not ask for, and never offer region or cell type columns stored with the data
  (they are kept only for offline evaluation).
- When the user asks to put the cells between their regions (gaps, ragged boundaries) into the region they belong
  to: preview_region_fill (by position only), look at the patches that matter with show_region_fill on the H&E,
  then assign_region_cells with the patches that hold up (only on request; the user reviews it, outlines stay).
  Do not move cells because of their cell type: a region prior would then confirm itself.
- A torn or folded piece of tissue on the H&E (cells beside a piece of tissue that holds no cells, which the local
  registration correction cannot follow): the user outlines the piece's cells as a region; fit_piece(region) searches
  its rotation, shift and mirroring (folded over) and shows it; when the verdict is "fits" and the yellow cells sit
  on the piece's nuclei in the figure, propose it with move_piece (the user reviews it with the figure); otherwise
  say why it does not fit (a folded piece out of focus has no nuclei to match: the cells can stay where they are).
- Conflicts (conflict_summary, show_conflicts) are cells whose expression contradicts a prior; describe H&E crops
  cautiously and record verdicts only when the user states them.

Reviewing an annotation ("does this make sense?", "check my annotation", "review my annotation")
- Use your knowledge of tissue biology, marker genes and histology; the tools give the evidence. Cover:
  1. Markers (review_markers): for every lineage and cell type, how often its listed markers are detected compared
     with its siblings (flags: weak, not specific, negative marker detected), and the genes that actually separate
     those cells here (top_genes). Ask: do the top genes fit the label (e.g. a "dendritic cell" group whose top genes
     are C1QC, LYVE1 is macrophages)? Are markers weak in this tissue or panel, with a better panel gene available?
     Do Unassigned cells (or an unassigned cluster in list_clusters) express a recognizable type missing from the
     list (e.g. CD177, CSF3R, S100A12: neutrophils)? Is the composition plausible for the tissue and site? Checks the
     user marked as solved (solved_by_user) are settled: mention them only if asked.
  2. Space (show_tissue_overview, check_plausibility, show_area): identify compartments on the H&E (mucosa with
     crypts, submucosa, muscularis propria, lymphoid aggregates, ulcers) and areas where the annotation contradicts
     them (e.g. epithelial cells in the muscle layer). Candidate areas come from the annotation alone; vascular
     smooth muscle, the muscularis mucosae, and deep crypts or fistulas in IBD are real. Lineages per drawn region
     show e.g. epithelial cells outside the user's "mucosa".
- When the app asks you to check the review candidates (or the user asks for a review): look at every candidate area
  on the H&E with show_area and judge every marker check, then record each verdict with record_review (keep only real
  annotation problems; one-line reasons). Out-of-place cells are computed from the labels of their neighbors, so
  cells at the edge of a sheet, next to mislabeled neighbors, or where few cells were segmented (e.g. mucus-rich
  crypts) are false alarms. Items already judged carry a claude_verdict; skip them. Never change the annotation
  while checking. The Review panel shows only what you kept.
- Then give a short, prioritized list of concrete suggestions, each with its evidence and the action:
  - a revised marker list: propose_marker_list(text, summary) validates it and puts it in the viewer's Review panel,
    where the user opens it in "Markers..." and annotates again keeping regions and priors (or
    start_session(markers_text=..., carry_from_current=True) when the user asks you to run it; the app follows);
  - a local re-clustering with sharper subtype markers (preview_recluster);
  - areas: suggest_areas (tiles, title, reason, action, region_name) draws them on the user's map; once the user
    adopted or drew the region (list_area_suggestions, list_regions), apply the action after a go-ahead: usually an
    exclude prior (preview_constraint mode="exclude", then add_constraint) so each cell gets its next most likely
    type, or label_cells(region=..., only_types=...) when the user names the target type;
  - priors or cluster labels.
  Change nothing without the user's agreement.

Reviewing (Agent Review: the user picks one part of the analysis; or when asked what you think of it)
- Review only the part asked for, as a pathologist and a single-cell analyst together. Gather the evidence:
  annotation: sample_context, type_confusion, type_structure (large or suspicious types), review_markers,
    qc_overview, check_registration; categories taxonomy (types to merge, split or add, citing sub-cluster sizes),
    context (site and disease and what they imply: inflamed IBD neutrophils, CXCL1+ colonocytes, plasma cell shifts;
    ileum vs colon Paneth cells), qc (registration, low-quality tiles, spillover, segmentation; only real issues);
  organization: show_tissue_overview, preview_compartment, spatial_profile, list_regions, composition(rows="region"),
    check_registration_map (crypts, vessels and similar structures first; judge each area on its close-up);
    compartments (outlines that fit the H&E, proposed as regions with the prior that would use them; check every
    outline on the H&E), organization (expected gradients and neighbors: stem cells at the crypt base, B cells in
    follicles, ICC between muscle layers, found or contradicted), registration (areas whose cells really sit off
    their nuclei: the local correction, action {"kind": "registration", "mode": "local"});
  embedding: embedding_facts, show_embedding (by cell type, leiden, slice); clusters (mixed clusters, types split
    over clusters, clusters that are another type), integration (clusters of one slice, batch), parameters
    (resolution, neighbors);
  niches: niche_summary, show_niches, group_markers(by="niche"), composition(rows="niche", cols="region"),
    niche_context; niches (what each is on the H&E, names, niches to merge or split), parameters (k, number);
  compare: recent_analyses(section="compare"), then rerun them (de_genes, group_markers, composition); findings
    (what the genes or shares mean), design (replicate slices, spillover, confounded groups), next (the comparison
    that would settle it);
  signaling: recent_analyses(section="signaling"), cell_communication / compare_communication, show_communication
    for striking pairs; interactions (which fit the tissue, which are artifacts), design (panel coverage, distances,
    replicate slices), next.
- Then review your own advice: draft_advice(scope, items), check every draft as a skeptical reviewer (rerun the tool
  that would confirm or refute it, look at the figure, weigh spillover, segmentation, batch, too few cells, panel
  limits), and record check_advice with keep / drop, a reason and the evidence. Only kept items reach the user; a
  few well-supported items beat many. If you can start a subagent, let it check the drafts without your reasoning.
  Do not repeat what the review candidates already show. Regions are the user's: propose outlines to add, never
  region columns of the dataset. Say "about" when a number is an estimate. Change nothing while reviewing.
- A report (the viewer's Report button, or when asked): gather annotation_facts, session_summary, sample_context,
  list_advice, check_plausibility and review_markers, and save it with write_report.

Aligning cells to an image (the app's "Align with Claude", or when asked): alignment_status, then
search_alignment; look at the figures of the best candidates (show_alignment), pick the one whose footprint covers
the right section in the right orientation, refine_alignment it, check the close-ups (centroids on nuclei) and the
refinement numbers, then apply_alignment. If no candidate fits, say so and ask for the registration file.

Niches (the viewer's Niche tab, or when asked): run_niches averages each cell's identity (NMF gene programs, or
its cell type) over its k nearest cells and clusters the averages. When asked to choose k: niche_context, then k
per drawn region or radius_um so neighborhoods span the tissue scale that matters (small in dense epithelium, larger
in sparse stroma). Afterwards niche_summary and show_niches; name the niches (name_niches) from their composition,
genes and place, cautiously; group_markers(by="niche") gives each niche's genes.

Comparing groups (when asked, e.g. "what differs between N3 and N5", "fibroblasts in the mucosa vs elsewhere",
"which niches does UC have more of"): analysis_groupings lists what cells can be grouped by (cell type, lineage,
region, niche, leiden, selection = the cells the user selected in the viewer, slice, and the slices' metadata such as
condition) and their values. de_genes compares a group with the rest or with another group, within a subset;
group_markers gives the top genes of every group; composition splits one grouping over another (cell types per niche,
niches per region, niches per slice). scope="project" pools the latest annotation of every slice of the project.
Rank DE genes by AUC and fold change rather than p-values (cells are not independent replicates); with several
slices trust what holds in every slice (slices_same_direction / slices_separate), and with one slice per side say a
slice effect cannot be told from a group effect. A gene of the neighboring cell types may be spillover. Interpret the
genes with care and say which claims come from the data and which from general knowledge.

Trajectories (when asked, or the app's Trajectory tab): trajectory orders cells on a spatial axis between two groups
(e.g. the crypt axis from stem cells to colonocytes, portal to central regions), or by diffusion pseudotime from a
root group within one lineage, optionally with tissue neighbors joined in (mixed, alpha). Read the rising, falling
and transient genes and the composition along it against the known organization; an axis or root imposes the order,
so a gradient along a spatial axis is position, not time, and spillover from neighboring cells follows the axis too.

Cell-cell communication (when asked, or the app's Signal tab): cell_communication scores ligand-receptor pairs
on neighboring cells between groups (cell types by default); compare_communication repeats it per condition, slice,
region or niche; show_communication shows on the H&E where one signal is received or sent. Co-expression in
neighbors suggests signaling, it does not prove it; mention how few pairs a small panel covers, broadly expressed
genes and spillover, and check striking results on the map.

Every change is a version (history, diff, revert). export writes a CSV with one row per cell. Reply briefly, in the
language of the user's message, keeping cell type and region names as they are."""


mcp = FastMCP("stagent", instructions=INSTRUCTIONS)
_active: dict = {}


def _log_calls(manager):
    """Every tool call goes to agent_calls.jsonl (report.py), so a problem report shows what the agent did."""
    import time
    call = manager.call_tool

    async def logged(name, arguments, *a, **kw):
        from . import report
        t0 = time.time()
        s = _active.get("session")
        try:
            out = await call(name, arguments, *a, **kw)
        except Exception as e:
            import traceback
            cause = e.__cause__ or e  # FastMCP wraps the tool's own exception in a ToolError
            report.log_call(name, arguments, time.time() - t0, error=f"{type(e).__name__}: {e}",
                            session=s.id if s else None, via=_active.get("via"),
                            trace="".join(traceback.format_exception(cause)))
            raise
        res = out[0] if isinstance(out, tuple) and len(out) == 2 and isinstance(out[0], (list, tuple)) else out
        if isinstance(res, (list, tuple)):  # content blocks: their text; images are left out
            res = [getattr(b, "text", f"<{getattr(b, 'type', 'content')}>") for b in res]
        s = _active.get("session") or s
        report.log_call(name, arguments, time.time() - t0, result=res, session=s.id if s else None,
                        via=_active.get("via"))
        return out
    manager.call_tool = logged


_log_calls(mcp._tool_manager)


def _heartbeat():
    """Note that Claude is connected (the app's About shows when Claude last used the tools)."""
    import time
    now = time.time()
    if now - _active.get("beat", 0) < 20:
        return
    _active["beat"] = now
    try:
        from . import config
        json.dump({"pid": os.getpid(), "time": now, "via": _active.get("via", "app")}, open(config.HOME / "mcp_seen.json", "w"))
    except OSError:
        pass


def _leave_placeholder():
    """The sample's placeholder (not annotated yet) gives way once the sample is annotated: to the annotation open in
    the viewer when it is of this sample, else the sample's latest annotation."""
    from . import claudelink
    from .session import list_sessions
    cur = _active["session"]
    t = list_sessions()
    real = t[(t.sample_id == cur.sample_id) & t.method.isin(["markers", "scref", "existing"])] if len(t) else t
    if not len(real):
        return
    shown = claudelink.viewer_state().get("session")
    _active["session"] = Session(shown if shown in set(real.session_id) else
                                 real.sort_values("created_at").session_id.iloc[-1])


def _s() -> Session:
    _heartbeat()
    if "session" not in _active:
        raise ValueError("no active session; call start_session or resume_session first")
    if _active["session"].scheme["kind"] == "none":
        _leave_placeholder()
    if _active["session"].stale():  # the cells moved (registration corrected in the app): reload
        _active["session"] = Session(_active["session"].id)
    _active["session"].actor = "agent"  # its changes for one request are undone together in the viewer
    return _active["session"]


def _j(obj) -> str:
    return json.dumps(obj, default=str, ensure_ascii=False)


def _records(df, top=None):
    df = df if top is None else df.head(top)
    return df.to_dict("records")


def _spec(text, types, mode, weight, kind, regions, max_depth_um, neighbor_types, min_frac):
    return {"text": text, "types": types, "mode": mode, "weight": weight, "kind": kind, "regions": regions,
            "max_depth_um": max_depth_um, "neighbor_types": neighbor_types, "min_frac": min_frac}


# ---------------- changes to the data wait for the user's Accept in the app ----------------
# A tool that changes data stores a proposal (proposals.py) with a description of what it would do; the app shows it
# in a review window and runs the tool only when the user clicks Accept (apply_proposal). The tool waits up to
# REVIEW_WAIT seconds for the decision and returns the outcome; tools served by the app itself (/mcp, its chat box) run
# on the app's event loop and cannot wait, so they return at once.
GATED = {}
SAME_PROCESS = False  # set by the webapp
REVIEW_WAIT = float(os.environ.get("STAGENT_REVIEW_WAIT", 45))


def _reviewed(describe):
    """describe(args) -> {"title", "details": [lines], "preview": {...}} (raises ValueError for invalid input)."""
    import functools, inspect

    def deco(fn):
        sig = inspect.signature(fn)

        @functools.wraps(fn)
        def tool(*a, **kw):
            from . import proposals, claudelink
            b = sig.bind(*a, **kw); b.apply_defaults()
            args = dict(b.arguments)
            d = describe(args)
            s = _active.get("session")
            rec = proposals.create(tool=fn.__name__, args=d.pop("args", args), title=d["title"],
                                   details=d.get("details", []), preview=d.get("preview"),
                                   session_id=d.get("session_id", s.id if s else None),
                                   sample_id=d.get("sample_id", s.sample_id if s else None),
                                   version=d.get("version"), via=_active.get("via", "app"))
            if not SAME_PROCESS and claudelink.viewer_alive():
                rec = proposals.wait(rec["id"], REVIEW_WAIT)
            return _proposal_outcome(rec)
        GATED[fn.__name__] = fn
        return tool
    return deco


def _proposal_outcome(rec):
    from . import claudelink
    st = rec["status"]
    if st == "done":
        if rec.get("session_after") and (not _active.get("session") or _active["session"].id != rec["session_after"]):
            _active["session"] = Session(rec["session_after"])  # e.g. the annotation start_session made
        try:
            result = json.loads(rec["result"])
        except (TypeError, ValueError):
            result = rec.get("result")
        return _j({"status": "accepted by the user and applied", "proposal": rec["id"], "result": result})
    if st == "dismissed":
        return _j({"status": "dismissed by the user", "proposal": rec["id"],
                   "note": "nothing was changed; ask the user what to do differently instead of proposing it again"})
    if st == "failed":
        return _j({"status": "accepted but failed", "proposal": rec["id"], "error": rec.get("error")})
    return _j({"status": "waiting for the user's review", "proposal": rec["id"],
               "app_open": claudelink.viewer_alive(),
               "next": "the change is shown in a review window in the STAnnotator app and applies only when the user "
                       "clicks Accept (Dismiss drops it). Tell the user what you proposed and why; later call "
                       "pending_changes to see the outcome." + ("" if claudelink.viewer_alive() else
                       " The app is not open: ask the user to open it to review the change.")})


def apply_proposal(rec, get_session=None):
    """Run a proposal's tool (the app does it when the user clicks Accept). Returns (result, active session id)."""
    fn = GATED[rec["tool"]]
    if rec.get("session_id"):
        _active["session"] = get_session(rec["session_id"]) if get_session else Session(rec["session_id"])
    out = fn(**rec["args"])
    return out, (_active["session"].id if _active.get("session") else None)


def _words(args, key="text"):
    return [f"Your words: \u201c{args[key]}\u201d"] if args.get(key) else []


def _names(xs, n=6):
    xs = list(xs or [])
    return ", ".join(map(str, xs[:n])) + (f" and {len(xs) - n} more" if len(xs) > n else "")


def _transitions(rows, n=8):
    return [{"from": r["from"], "to": r["to"], "region": None if r.get("region") in (None, "unknown", "none") else r["region"],
             "cells": int(r.get("count", r.get(0, 0)))}
            for r in rows[:n]]


def _d_import_sample(a):
    return {"title": f"Import a new sample {a['name']}", "session_id": None, "sample_id": a["name"],
            "details": [f"Expression and cell positions from {_names(a['h5ad'])}",
                        f"Image: {a['he'] or 'none'}; registration: {a['registration'] or 'none'}",
                        *([f"Section {a['section']} of column {a['section_column']}"] if a.get("section") else []),
                        "Adds a sample; nothing that exists now changes."]}


def _d_import_cosmx_merfish(a):
    from . import importer
    sl = importer.inspect_flat_many(a["kind"], a["paths"])["slices"]
    det = [f"{x['key']}: {x['n_cells']:,} cells, {x['n_genes']} genes ({x['n_controls']} control probes left out), "
           f"positions in {x['positions']}" for x in sl]
    if len(sl) > 1:
        raise ValueError(f"these paths hold {len(sl)} slides or regions: import one at a time")
    det += [f"Image: {a['he'] or 'none (cells on a blank background)'}; registration: {a['registration'] or 'none'}",
            "Adds a sample; nothing that exists now changes."]
    return {"title": f"Import {importer.FLAT_NAMES[a['kind']]} output as sample {a['name']}",
            "session_id": None, "sample_id": a["name"], "details": det}


def _d_start_session(a):
    from . import markers
    from .session import panel_genes, filtered_genes
    det = []
    if a.get("annotation_csv"):
        det.append(f"Loads the cell types of {a['annotation_csv']} (column {a.get('label_column')})")
    elif a.get("reference_h5ad"):
        det.append(f"Labels the cells with the single-cell reference {a['reference_h5ad']} (labels: {a.get('label_column')})")
    else:
        p = markers.parse(a.get("markers_text") or "", panel_genes(a["sample_id"]), filtered_genes(a["sample_id"]))
        if not p["ok"]:
            raise ValueError("marker list has errors: " + "; ".join(p["errors"][:8]))
        det.append(f"Marker list with {p['n_types']} cell types under {len(p['tree'])} lineages: "
                   + _names([n["name"] for n in p["tree"]], 8))
        det += [f"Warning: {w}" for w in p["warnings"][:4]]
        from . import config, timing
        from .session import knn_mode
        knn = knn_mode(a.get("knn"))
        n = json.load(open(config.SAMPLES / a["sample_id"] / "meta.json")).get("n_cells") or 0
        est = timing.estimate(n, len(panel_genes(a["sample_id"])), p["tree"], "auto", knn, len(p["leaves"]))
        det.append(f"Neighbor search: {est['backend']}"
                   + (" (the same result on every computer)" if knn == "exact" else "")
                   + f". Estimated time for {n:,} cells: {timing.words(est['seconds'])} ({est['basis']}).")
    cur = _active.get("session")
    if a.get("carry_from_current") and cur:
        det.append(f"Copies the regions, priors and corrections of the open annotation {cur.id}")
    det.append("A new annotation; earlier annotations stay and can be reopened under Annotation.")
    return {"title": f"Start a new annotation of {a['sample_id']}", "sample_id": a["sample_id"], "details": det}


def _d_constraint(a):
    s = _s()
    spec = _spec(a["text"], a["types"], a["mode"], a["weight"], a["kind"], a["regions"], a["max_depth_um"],
                 a["neighbor_types"], a["min_frac"])
    pv = s.preview_constraint(spec)
    return {"title": f"Add a prior: {pv['interpretation']}", "version": pv["based_on_version"],
            "details": [*_words(a), f"Types: {_names(pv['types_resolved'], 10)}",
                        f"{pv['n_cells_currently_of_these_types']:,} cells are of these types now; "
                        f"{pv['n_cells_of_these_types_penalized']:,} of them are where the prior speaks against them",
                        "Creates a new version (Undo takes it back)."],
            "preview": {"n_label_changes": pv["n_label_changes"], "top_transitions": _transitions(pv["top_transitions"])}}


def _constraint(s, cid):
    c = next((c for c in s.state()["constraints"] if c["id"] == cid), None)
    if c is None:
        raise ValueError(f"no constraint {cid}; constraints: {[x['id'] for x in s.state()['constraints']]}")
    return c


def _d_update_constraint(a):
    from . import constraints as cmod
    s = _s()
    c = _constraint(s, a["constraint_id"])
    ch = {k: a[k] for k in ("weight", "types", "regions", "mode", "max_depth_um", "neighbor_types", "min_frac", "text")
          if a.get(k) is not None}
    return {"title": f"Change prior {c['id']}: {cmod.describe(c)}", "version": s.state()["version"],
            "details": [f"{k}: {c.get(k)!r} \u2192 {v!r}" for k, v in ch.items()] + ["Creates a new version."]}


def _d_remove_constraint(a):
    from . import constraints as cmod
    s = _s()
    c = _constraint(s, a["constraint_id"])
    return {"title": f"Remove prior {c['id']}", "version": s.state()["version"],
            "details": [cmod.describe(c), *([f"Its words: \u201c{c['text']}\u201d"] if c.get("text") else []),
                        "The cells are labeled again without it; creates a new version."]}


def _label_breakdown(s, rows_or_cells, by_cell=True):
    import numpy as np, pandas as pd
    pred = s.result()["pred"]
    rows = [s.row_of(c) for c in rows_or_cells] if by_cell else rows_or_cells
    t = pd.Series(np.array(s.types)[pred[rows]]).value_counts()
    return _names([f"{k} {v:,}" for k, v in t.items()], 6)


def _d_label_cells(a):
    s = _s()
    cells = s._cells_for(a["region"], None, a.get("only_types"))
    return {"title": f"Label {len(cells):,} cells of region {a['region']} as {a['cell_type']}",
            "version": s.state()["version"],
            "details": [*_words(a), *([f"Only the cells now labeled {_names(a['only_types'])}"] if a.get("only_types") else []),
                        f"Labeled now: {_label_breakdown(s, cells)}",
                        "These cells are locked to the new type and the expression model learns from them, so similar "
                        "cells elsewhere may change too. Creates a new version (Undo takes it back)."]}


def _d_label_cluster(a):
    s = _s()
    mres = s.marker_result(s.state())
    if a["cluster_id"] not in mres["members"]:
        raise ValueError(f"no cluster {a['cluster_id']}; see list_clusters")
    rows = list(mres["members"][a["cluster_id"]])
    return {"title": f"Label cluster {a['cluster_id']} ({len(rows):,} cells) as {a['cell_type']}",
            "version": s.state()["version"],
            "details": [*_words(a), f"Labeled now: {_label_breakdown(s, rows, by_cell=False)}",
                        "These cells are locked to the new type and the expression model learns from them, so similar "
                        "cells elsewhere may change too. Creates a new version."]}


def _d_remove_correction(a):
    s = _s()
    c = next((c for c in s.state().get("corrections", []) if c["id"] == a["correction_id"]), None)
    if c is None:
        raise ValueError(f"no correction {a['correction_id']}")
    return {"title": f"Undo hand correction {c['id']}", "version": s.state()["version"],
            "details": [f"{len(c['cells']):,} cells labeled {c['label']}" + (f" in region {c['region']}" if c.get("region") else ""),
                        *([f"Its words: \u201c{c['text']}\u201d"] if c.get("text") else []),
                        "The expression model is computed again without it; creates a new version."]}


def _d_recluster(a):
    s = _s()
    det = [*_words(a), f"Resolution: {a['resolution']}"]
    if a.get("marker_edits"):
        det.append("Marker changes: " + "; ".join(f"{k}: {v}" for k, v in a["marker_edits"].items()))
    det.append("Regions, priors and corrections stay; the expression model learns from the new marker labels, so "
               "cells elsewhere may change too. Creates a new version.")
    return {"title": f"Re-cluster {a['target']} in the marker step", "version": s.state()["version"], "details": det}


def _d_revert(a):
    s = _s()
    cur = s.state()["version"]
    if a["version"] not in s.versions():
        raise ValueError(f"no version {a['version']}")
    d = s.diff(cur, a["version"])
    return {"title": f"Go back to version {a['version']} (now v{cur})", "version": cur,
            "details": [f"v{a['version']}: {s.state(a['version'])['message']}",
                        "Its regions, priors and corrections come back as a new version; the history stays."],
            "preview": {"n_label_changes": d["n_label_changes"], "top_transitions": _transitions(d["top_transitions"])}}


def _d_apply_alignment(a):
    sid = _sample_id(a.get("sample_id"))
    return {"title": f"Put the cells of {sid} on the image with alignment {a['candidate']}", "sample_id": sid,
            "args": a | {"sample_id": sid},
            "details": ["The cells move to where this alignment places them on the image; regions and outlines drawn "
                        "before move with them.", "Undo last alignment (Import dialog) takes it back."]}


def _d_run_niches(a):
    s = _s()
    det = [f"Identity: {a['identity']}" + (f" ({a['n_programs']} NMF programs)" if a["identity"] == "nmf" else ""),
           f"Neighborhood: within {a['radius_um']} um" if a.get("radius_um") else f"Neighborhood: {a['k']} nearest cells",
           f"{a['n_niches']} niches"]
    if a.get("k_by_region"):
        det.append("Per region: " + ", ".join(f"{k}: {v} cells" for k, v in a["k_by_region"].items()))
    det.append("Replaces the current niches of this annotation (their names too).")
    return {"title": "Find niches", "version": s.state()["version"], "details": det}


def _d_name_niches(a):
    _s()
    return {"title": f"Name {len(a['names'])} niches", "details": [f"{k} \u2192 {v}" for k, v in a["names"].items()]}


def _d_record_verdicts(a):
    _s()
    if a["verdict"] not in vocab.VERDICTS:
        raise ValueError(f"verdict must be one of {list(vocab.VERDICTS)}")
    return {"title": f"Record the verdict {a['verdict']} for {len(a['cells']):,} conflict cells",
            "details": [vocab.VERDICTS[a["verdict"]], f"Reviewer: {a['reviewer']}", *_words(a, "note")]}


def _d_resolve_rule_suggestion(a):
    s = _s()
    sg = next((x for x in s.suggestions() if x["id"] == a["suggestion_id"]), None)
    if sg is None:
        raise ValueError(f"no rule suggestion {a['suggestion_id']}")
    return {"title": f"Mark the reviewer's rule {a['suggestion_id']} as {a['status']}",
            "details": [f"The rule: \u201c{sg.get('text', '')}\u201d", *_words(a, "note")]}


def _d_import_visium(a):
    from . import importer
    x = importer.inspect_visium(a["path"], a.get("segmentation"))
    cells = a.get("unit", "auto") != "bins" and x.get("segmentation") and x.get("fine_bins")
    what = "Visium spots" if x["visium_spots"] else f"{x['bin_um']:g} um bins"
    how = (f"cells from the outlines {x['segmentation']} (about {x.get('n_outlines_about') or 0:,}), the 2 um bins split "
           "between them by area" if cells else
           f"pooled into {x['analysis_um']:g} um bins (about {x.get('n_analysis_bins_about') or 0:,})"
           if x["analysis_um"] != x["bin_um"] else "used as they are")
    return {"title": f"Import a new sample {a['name']} (Visium HD)", "session_id": None, "sample_id": a["name"],
            "details": [f"{x['folder']}: {x.get('n_bins', 0):,} {what} under the tissue, {x.get('n_genes', 0):,} genes; {how}",
                        f"Image: {a['he'] or 'the output tissue_hires_image.png'} (placed by Space Ranger's alignment)",
                        "Adds a sample; nothing that exists now changes."]}


@mcp.tool()
@_reviewed(_d_import_visium)
def import_visium(name: str, path: str, he: Optional[str] = None, description: str = "",
                  segmentation: Optional[str] = None, unit: str = "auto") -> str:
    """Import a Visium HD (or Visium) Space Ranger output: path = outs/, binned_outputs/ or one square_XXXum folder.
    With cell outlines (Space Ranger 4's segmented_outputs/cell_segmentations.geojson, found by itself, or
    segmentation = a GeoJSON of outlines in full-resolution pixels) the cells are the unit (unit "auto" or "cells"):
    the 2 um bins are split between the cells by area. Otherwise (or unit "bins") the bins are the unit at 8 um: 2 um
    bins pooled into 8 um, bins of 8 um or more (16 um, Visium spots) used as they are. he: the
    full-resolution image Space Ranger was given (no registration needed); without it, the output's
    tissue_hires_image.png. Whole-transcriptome data: the PCA and the expression model use the 3,000 most variable
    genes (plus the marker list's). The nuclei-based registration check does not apply to bins."""
    from . import importer
    log = []
    meta = importer.import_visium(name, path, he=he, description=description, progress=log.append,
                                  segmentation=segmentation, unit=unit)
    return _j({"imported": meta["sample_id"], "unit": meta["source"]["unit"], "n": meta["n_cells"],
               "n_genes": meta["n_genes"], "bin_um": meta["source"]["bin_um"], "log": log[-12:]})


def _d_import_he_only(a):
    from . import importer
    upp = a.get("um_per_px") or importer.image_pixel_size(a["image"])
    return {"title": f"Import a new sample {a['name']} (H&E only)", "session_id": None, "sample_id": a["name"],
            "details": [f"Image: {a['image']}, " + (f"{upp:.4f} um per pixel" if upp else "pixel size unknown"),
                        "No cells: regions are drawn and named on it, exported as GeoJSON, or given to the import of the "
                        "data registered to the same image (import_sample / import_xenium ... regions_from).",
                        "Adds a sample; nothing that exists now changes."]}


@mcp.tool()
@_reviewed(_d_import_he_only)
def import_he_only(name: str, image: str, um_per_px: Optional[float] = None, description: str = "") -> str:
    """Import an H&E (TIFF / OME-TIFF, SVS, NDPI, PNG, JPEG) without any cells, e.g. for a pathologist to draw and
    name regions before the expression data exist. um_per_px: the pixel size, read from the image when not given (else
    unknown: no micrometers are shown). The regions drawn on it go to the data later: in the app's Import dialog
    ("Regions from"), or exported as GeoJSON (export_regions) and given as the import's regions file. Tools about cells
    do not apply to such a sample."""
    from . import importer
    log = []
    m = importer.import_image(name, image, um_per_px=um_per_px, description=description, progress=log.append)
    return _j({"imported": m["sample_id"], "image_only": True, "um_per_px": m["he_um_per_px"] if m.get("um_per_px_known") else None,
               "log": log[-6:]})


def _d_import_stereoseq(a):
    from . import importer
    x = importer.inspect_stereoseq(a["paths"])
    cells = a.get("unit", "auto") != "bins" and (x["cellbin_gef"] or x["cell_gem"])
    return {"title": f"Import a new sample {a['name']} (Stereo-seq)", "session_id": None, "sample_id": a["name"],
            "details": [("Cells from " + (x["cellbin_gef"] or x["cell_gem"]) + (f" ({x['n_cells']:,} cells)" if x.get("n_cells") else ""))
                        if cells else f"8 um bins pooled from the DNBs of {x['gem']}",
                        f"Image: {a.get('image') or x['image'] or 'none (a blank background)'}" + (" (in the DNBs' coordinates)"
                        if a.get("image") or x["image_registered"] else " (not SAW's registered image: it waits for an alignment)"),
                        "Adds a sample; nothing that exists now changes."]}


@mcp.tool()
@_reviewed(_d_import_stereoseq)
def import_stereoseq(name: str, paths: list[str], image: Optional[str] = None, description: str = "",
                     unit: str = "auto") -> str:
    """Import a Stereo-seq (STOmics SAW) output: paths = its folder, or its files (<SN>.cellbin.gef, <SN>.tissue.gem.gz,
    <SN>_HE_regist.tif or the ssDNA _regist.tif). The cells are the unit when SAW segmented them (a cell-bin GEF, or a
    GEM with a CellID column; unit "auto" or "cells"), else (or unit "bins") the DNBs (0.5 um) are pooled into 8 um
    bins. The registered image is in the DNBs' coordinates (no registration needed); image: another such image."""
    from . import importer
    log = []
    meta = importer.import_stereoseq(name, paths, image=image, description=description, progress=log.append, unit=unit)
    return _j({"imported": meta["sample_id"], "unit": meta["source"]["unit"], "n": meta["n_cells"],
               "n_genes": meta["n_genes"], "image": meta["image_kind"], "log": log[-12:]})


# ---------------- data and viewer ----------------
@mcp.tool()
@_reviewed(_d_import_cosmx_merfish)
def import_cosmx_merfish(kind: str, name: str, paths: list[str], he: Optional[str] = None,
                         registration: Optional[str] = None, description: str = "") -> str:
    """Import one CosMx slide (kind "cosmx": its AtoMx flat files <slide>_exprMat_file.csv and
    <slide>_metadata_file.csv) or one MERFISH region (kind "merfish": a MERSCOPE region's cell_by_gene.csv and
    cell_metadata.csv): paths = those two files or their folder; several slides or regions are imported one after the
    other. The instruments' own images are not used: the cells are shown on a blank background, or on the H&E given
    with a registration file (affine CSV a11,a12,a21,a22,tx,ty from cell um to H&E pixels); an H&E without
    registration waits for Align with Claude."""
    from . import importer
    if kind not in importer.FLAT_FILES:
        raise ValueError("kind is 'cosmx' or 'merfish'")
    log = []
    src = importer.flat_sources(kind, paths)
    if len(src) > 1:
        raise ValueError(f"these paths hold {len(src)} slides or regions ({', '.join(x['key'] for x in src)}): import one at a time")
    meta = importer.import_flat(kind, name, src[0], he=he, registration=registration, description=description,
                                progress=log.append)
    return _j({"imported": meta["sample_id"], "n_cells": meta["n_cells"], "n_genes": meta["n_genes"],
               "image": meta["image_kind"], "log": log[-12:]})


@mcp.tool()
@_reviewed(_d_import_sample)
def import_sample(name: str, h5ad: list[str], he: str, registration: str, section_column: Optional[str] = None,
                  section: Optional[str] = None, region_column: Optional[str] = None, counts_layer: Optional[str] = None,
                  description: str = "") -> str:
    """Import one tissue section from files on this computer: h5ad file(s) with raw counts and cell coordinates in
    um (several files are merged; an .rds / .rda with a Seurat or SpatialExperiment object works too, read without R), the registered H&E image (pyramidal TIFF / OME-TIFF), and the registration file
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
    """Load the sample and annotation the user is looking at in the STAnnotator app (or the web viewer) and
    make it the active session. Call it first when the user refers to what they see."""
    _heartbeat()
    from . import claudelink
    st = claudelink.viewer_state()
    if not st.get("session"):
        return _j({"error": "the viewer shows no annotation; ask the user to open one (or start_session)",
                   "viewer_running": claudelink.viewer_alive(st)})
    _active["session"] = Session(st["session"])
    return _j({"session_id": st["session"], "sample": st.get("sample"), "summary": _active["session"].summary()})


@mcp.tool()
def open_viewer(port: int = 8765) -> str:
    """Show the active session to the user. If the STAnnotator app (or a web viewer) is running, it switches to
    the session by itself; otherwise a web viewer is started and its link returned."""
    import subprocess, time
    from . import claudelink, config
    st = claudelink.viewer_state()
    s = _active.get("session")
    if claudelink.viewer_alive(st):
        if s is not None:
            claudelink.request_show(s.sample_id, s.id)
        where = "the STAnnotator app" if st.get("mode") == "desktop" else f"http://localhost:{st['port']}"
        return _j({"shown_in": where, "session": s.id if s else None,
                   "note": "the viewer switches to this session within a few seconds"})
    url = f"http://localhost:{port}"
    if not claudelink.viewer_alive({"port": port}):
        log = open(config.HOME / "viewer.log", "a")
        from . import osutil
        subprocess.Popen(config.self_command("viewer", "--port", str(port)), stdout=log, stderr=log,
                         stdin=subprocess.DEVNULL, start_new_session=not osutil.WINDOWS, env=dict(os.environ),
                         creationflags=osutil.NO_WINDOW)
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
    """List available samples (sections) with description, cell count and whether an H&E exists."""
    return _j(_records(_list_samples().drop(columns=["has_pathologist_regions"], errors="ignore")))


@mcp.tool()
def list_sessions() -> str:
    """List saved annotation sessions that can be resumed."""
    return _j(_records(_list_sessions()))


@mcp.tool()
@_reviewed(_d_start_session)
def start_session(sample_id: str, markers_text: Optional[str] = None,
                  carry_from_current: bool = False, reference_h5ad: Optional[str] = None,
                  label_column: Optional[str] = None, group_column: Optional[str] = None,
                  majority_voting: bool = False, annotation_csv: Optional[str] = None, id_column: Optional[str] = None,
                  lineage_column: Optional[str] = None, knn: Optional[str] = None) -> str:
    """Start a new annotation session on a sample. The data is treated as unlabeled: cell types come from a
    hierarchical marker gene list (markers_text: indented text as in check_marker_list, or CSV with columns
    cell_type, parent, positive_markers, negative_markers; required).

    Single-cell reference instead (CellTypist-like classifier): reference_h5ad + label_column (the cell type column),
    optionally group_column (coarse classes, used as lineages) and majority_voting; check it with inspect_reference.
    An existing annotation instead: annotation_csv (one row per cell) with id_column, label_column and optionally
    lineage_column (check it with inspect_annotation); its labels stay the labels.
    carry_from_current: copy the regions, priors and corrections of the active session of the same sample (use it
    when re-annotating with a revised marker list). The first call to any tool that needs labels runs the annotation
    (the review window shows how long it should take here). Give the user a new viewer link afterwards (open_viewer).
    knn (marker list only): the neighbor search of the marker step, "approx" (default; pynndescent: fast on any
    computer, the same result only with the same pynndescent / numba versions) or "exact" (the same result on every
    computer, GPU or not; slow without a GPU for large samples). Use "exact" only when the user asks for it.
    """
    carry = _active["session"].id if carry_from_current and _active.get("session") else None
    ref = ex = None
    if annotation_csv:
        ex = {"path": annotation_csv, "id_column": id_column, "label_column": label_column,
              "lineage_column": lineage_column}
    elif reference_h5ad:
        ref = {"path": reference_h5ad, "label_column": label_column, "group_column": group_column,
               "majority_voting": majority_voting}
    method = "existing" if ex else "scref" if ref else "markers"
    _active["session"] = Session.create(sample_id, "none", method=method, markers_text=markers_text,
                                        carry_from=carry, reference=ref, existing=ex, knn=knn)
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
    predicted label changes to the user, then call add_constraint (the user accepts it in the app's review window).

    text: the user's own words, verbatim. types: cell type or group names from get_vocabulary.
    kind: 'region' (needs regions), 'depth' (needs max_depth_um: distance to the epithelium / crypt core),
          'neighborhood' (needs neighbor_types and min_frac, the share of the 15 nearest neighbors).
    mode: 'only' (types unlikely where the condition is false), 'exclude' (unlikely where true),
          'enrich' (more likely where true).
    weight: 'strong' | 'moderate' | 'weak' or a number (only/exclude: 0-1, enrich: >1).
    """
    return _j(_s().preview_constraint(_spec(text, types, mode, weight, kind, regions, max_depth_um, neighbor_types, min_frac)))


@mcp.tool()
@_reviewed(_d_constraint)
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
@_reviewed(_d_update_constraint)
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
@_reviewed(_d_remove_constraint)
def remove_constraint(constraint_id: str, message: Optional[str] = None) -> str:
    """Remove a constraint. Creates a new version."""
    s = _s()
    v = s.remove_constraint(constraint_id, message)
    return _j({"new_version": v, "diff": s.diff(v - 1, v)})


# ---------------- regions ----------------
@mcp.tool()
def list_regions(version: Optional[int] = None) -> str:
    """Regions of the session. Regions are drawn by the user on the map, or added from an outline you proposed
    (propose_region) when the user accepts it;
    cells outside every region are 'unknown' and no region prior applies to them."""
    info = _s().regions_info(version)
    return _j({"source": info["source"], "regions": info["regions"], "n_unknown": info["n_unknown"],
               "polygons": [{k: p[k] for k in ("id", "name", "n_cells_drawn")} for p in info["polygons"]]})


@mcp.tool()
@_reviewed(_d_label_cells)
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


def _d_assign_region_cells(a):
    from . import regionfill
    s = _s()
    moves = a.get("resolved") or regionfill.pick(s, a["patches"], a.get("to_region"))
    n = sum(len(m["cells"]) for m in moves)
    return {"title": f"Put {n:,} cells between regions into " + ", ".join(m["to"] for m in moves),
            "version": s.state()["version"],
            "args": dict(a, resolved=moves),  # applied exactly as reviewed, even if the candidates are computed again
            "details": [*_words(a), *[f"{len(m['cells']):,} cells into {m['to']} (patches {_names(m['patches'])}; now in "
                                      + ", ".join(f"{k if k != 'unknown' else 'no region'}: {v:,}" for k, v in m["from"].items()) + ")"
                                      for m in moves],
                        "Found by position (distance to the regions and the regions of each cell's neighbors), not by "
                        "cell type. The outlines stay as drawn; the cells are listed under their new region in the "
                        "Regions panel. Creates a new version (Undo takes it back).",
                        "Show on map highlights these cells before you decide."],
            "preview": {"kind": "cells", "n_cells": n}}


@mcp.tool()
def preview_region_fill(max_gap_um: float = 30.0, min_share: float = 0.7, regions: Optional[list[str]] = None,
                        include_boundary: bool = True, between_only: bool = True) -> list:
    """Cells between the user's regions that belong to one of them, by position only (nothing changes): gaps (cells
    in no region within max_gap_um of two regions or more; between_only=False: of any region, which also grows the
    regions into tissue the user did not outline, only when asked: to the nearest region, or the one most of their
    in-region neighbors are in) and boundary cells (in one region but with >= min_share of their 15 nearest neighbors in another: an
    outline drawn a little off). Grouped into patches G1, G2, ... (connected areas: a gap's cells go to either
    side) with how many cells go into each region, the regions nearby and the cell types in them (for your judgment only: moving cells because of their type would let a
    region prior confirm itself). regions: only these targets. Figure: the patches on the sample and four of them on
    the H&E. Look at the patches that matter with show_region_fill, then propose the ones that hold up with
    assign_region_cells; the user reviews it in the app."""
    from . import regionfill
    s = _s()
    info = regionfill.candidates(s, max_gap_um, min_share, regions, include_boundary, between_only)
    png = regionfill.figure(s) if info["n_patches"] else None
    return [_j(info | {"path": png})] + ([Image(path=png)] if png else [])


@mcp.tool()
def show_region_fill(patches: list[str]) -> list:
    """Up to four patches of preview_region_fill on the H&E (their cells ringed in the target region's color, the
    other cells as dots in their region's color): check that each patch belongs to its target on the tissue."""
    from . import regionfill
    png = regionfill.figure(_s(), patches)
    return [_j({"path": png}), Image(path=png)]


@mcp.tool()
@_reviewed(_d_assign_region_cells)
def assign_region_cells(patches: list[str], to_region: Optional[str] = None, text: str = "",
                        resolved: Optional[list[dict]] = None) -> str:
    """Propose putting the cells of patches (ids from preview_region_fill, checked on the H&E) into their target
    region (to_region: one region for all of them instead). The user sees it in the app's review window, can show
    the cells on the map, and accepts or dismisses it; the outlines are not changed. text: the user's words.
    resolved: filled in by the app (leave it out)."""
    from . import regionfill
    s = _s()
    moves = resolved or regionfill.pick(s, patches, to_region)
    v = s.assign_region_cells(moves, author="agent (accepted by the user)", text=text)
    return _j({"new_version": v, "message": s.state(v)["message"], "regions": s.regions_info(v)["regions"]})


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
    genes of Unassigned cells. node: only that node and its subtypes. Checks the user marked as solved in the
    viewer's Review panel are left out of flags and listed in solved_by_user: do not raise them again."""
    from . import markerqc
    s = _s()
    q = dict(markerqc.check(s, version))
    flags, unassigned = markerqc.open_checks(s, version)
    from . import claudereview
    judged = claudereview.load(s)["markers"]
    flags = [f | {"key": markerqc.check_key(f), "claude_verdict": judged.get(markerqc.check_key(f))} for f in flags]
    q.update(flags=flags, unassigned=unassigned, unassigned_verdict=judged.get("unassigned"),
             solved_by_user=markerqc.solved(s))
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
    suggested types and drawn regions; plus lineage counts per drawn region. Areas the user marked as solved in the
    viewer are left out and listed in solved_by_user: do not raise them again. Changes nothing."""
    from . import plausibility, claudereview
    s = _s()
    out = plausibility.summary(s, version, top)
    verdicts = {a["id"]: v for a, v in claudereview.items(s, version)[0]}
    for a in out["candidate_areas"]:
        v = verdicts.get(a["id"])
        a["claude_verdict"] = None if v is None else {"keep": v["keep"], "reason": v["reason"]}
    return _j(out)


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
def record_review(areas: Optional[list[dict]] = None, marker_checks: Optional[list[dict]] = None) -> str:
    """Record your verdicts on the computed review candidates, so the user's Review panel shows only real problems:
    areas = [{"id": "A3", "keep": true/false, "reason": "..."}] (ids from check_plausibility, after looking at each
    on the H&E with show_area), marker_checks = [{"key": "Cycling TAs|UBE2C+" or "unassigned", "keep": ..., "reason":
    ...}] (keys from review_markers). keep = a real annotation problem the user should see; false alarms are hidden.
    One-line reasons. Call it as you go (in batches); it never changes the annotation."""
    from . import claudereview
    return _j(claudereview.record(_s(), areas or [], marker_checks or []))


# ---------------- high-level review (advice and report) ----------------
@mcp.tool()
def sample_context() -> str:
    """What the sample is (description, panel size, image), the annotation method, and its lineage and type
    composition: the start of a high-level review (is this composition plausible for the tissue, site and disease?)."""
    from . import insights
    return _j(insights.sample_context(_s()))


@mcp.tool()
def preview_compartment(name: str, lineages: list[str], max_dist_um: float = 30.0, min_cells: int = 200) -> list:
    """Outline where the given lineages / types form tissue (cells within max_dist_um of where they are dense, holes
    filled), e.g. name="mucosa", lineages=["Epithelial"]; "muscularis" from smooth muscle; "follicles" from B cells
    (max_dist_um ~15). Returns the outline on the H&E and the stats (pieces, cells, composition). Check the figure
    against the H&E before proposing it (draft_advice scope "organization", action {"kind": "region", "compartment": name}). Changes
    nothing."""
    from . import insights, advice
    s = _s()
    comp = insights.compartment(s, lineages, max_dist_um=max_dist_um, min_cells=min_cells)
    spec = {"lineages": lineages, "max_dist_um": max_dist_um, "min_cells": min_cells}
    advice.store_compartment(s, name, comp, spec)
    png = insights.fig_compartment(s, comp, name)
    return [_j({k: v for k, v in comp.items() if k not in ("polygons_px", "cells")} | {"name": name, "path": png}),
            Image(path=png)]


def _region_outline(s, compartment, area):
    from . import advice, plausibility
    if bool(compartment) == bool(area):
        raise ValueError("give either compartment (from preview_compartment) or area (tiles or a candidate id)")
    if compartment:
        comp = advice.load(s)["compartments"].get(compartment)
        if comp is None:
            raise ValueError(f"no previewed compartment {compartment!r}: call preview_compartment first")
        return comp["polygons_px"], comp["n_cells"]
    info = plausibility.area_info(s, area)
    return [info["outline_px"]], info.get("n_cells", 0)


def _d_propose_region(a):
    from . import insights, regions as regmod
    s = _s()
    name = regmod.validate_name(re.sub(r"\s+", "_", str(a["name"]).strip()))
    polys, n = (a["resolved"], a.get("n_cells", 0)) if a.get("resolved") else _region_outline(s, a.get("compartment"), a.get("area"))
    exists = name in s.region_names()
    return {"title": f"Add region {name}" + (f" ({len(polys)} pieces)" if len(polys) > 1 else "") + f": {n:,} cells",
            "version": s.state()["version"],
            "args": dict(a, name=name, resolved=polys, n_cells=n),  # applied exactly as reviewed
            "details": [*_words(a), "Outlined by the agent from " + (f"the {a['compartment']} compartment (where its "
                        "lineages form tissue)" if a.get("compartment") else f"tiles {a['area']}") + "; the figure "
                        "shows the outline on the H&E.",
                        f"Adds {len(polys)} polygon(s) " + (f"to the existing region {name}" if exists else
                        f"as region {name}") + "; you can edit its points afterwards, and Undo takes it back."],
            "preview": {"kind": "image", "path": insights.fig_compartment(s, {"polygons_px": polys}, name)}}


@mcp.tool()
@_reviewed(_d_propose_region)
def propose_region(name: str, compartment: Optional[str] = None, area: Optional[str] = None, text: str = "",
                   resolved: Optional[list] = None, n_cells: int = 0) -> str:
    """Propose a region the user asked for: the outline of a compartment computed with preview_compartment
    (compartment=<its name>) or of tiles / a candidate area (area="F3:H6" or "A3"). The user sees the outline on the
    H&E in the app's review window and Accept adds it as region `name`. text: the user's words. A prior on the new
    region: propose it with add_constraint once the region exists. resolved, n_cells: filled in by the app (leave
    them out). Only when the user asked for the region."""
    s = _s()
    v = s.add_region_pieces(name, resolved, author="agent (accepted by the user)")
    return _j({"new_version": v, "message": s.state(v)["message"], "regions": s.regions_info(v)["regions"]})


# ---------------- regions from the H&E itself (also a sample of an H&E only) ----------------
@mcp.tool()
def tissue_priors(organ: Optional[str] = None) -> str:
    """The app's literature-based tissue priors (data/priors/tissue_priors.md) of an organ: its compartments to draw
    with how each looks on the H&E, and spatial priors (where each cell type lives, as rules of the constraint
    language, with exceptions and references). Without organ: the organs covered and how to read a rule. Use it to
    outline compartments on the H&E (segment_he, propose_he_region) and to suggest priors (preview_constraint),
    always with the user's agreement."""
    from . import hetissue
    return _j(hetissue.tissue_priors(organ))


@mcp.tool()
def segment_he(n_classes: int = 6) -> list:
    """Classes of the tissue read from the H&E itself (no cells needed: works on a sample of an H&E only): color
    deconvolution into hematoxylin and eosin and their texture, k-means into
    n_classes (1 = the most hematoxylin). Returns each class's area, hematoxylin, eosin, texture, distance from the
    edge and share at the edge, and a figure: the H&E with ticks in full-resolution pixels beside the class map.
    Match classes to the compartments of tissue_priors(organ) (e.g. colon: mucosa = dense purple crypts at the
    luminal edge, often 2-3 classes with high hematoxylin; lymphoid follicles = the darkest round blue spots;
    submucosa = pale, loose, low stain; muscularis propria = thick bands of the highest eosin and no texture), then
    propose_he_region. Changes nothing; the classes are kept for propose_he_region."""
    from . import hetissue
    s = _s()
    meta = hetissue.segment(s.sample_id, n_classes)
    png = hetissue.fig_classes(s.sample_id, s)
    return [_j({k: meta[k] for k in ("n_classes", "um_per_px", "classes", "pixel_size_known")} | {"path": png}),
            Image(path=png)]


@mcp.tool()
def show_he(box_px: Optional[list[float]] = None) -> list:
    """The H&E (whole, or box_px = [x0, y0, x1, y1] in full-resolution pixels, at about 1 um per pixel) with a grid
    labelled in full-resolution pixels and the regions drawn so far, to check an outline or read the points of one
    (propose_he_region polygon_px). Changes nothing."""
    from . import hetissue
    s = _s()
    ops = [o for o in s.state()["regions"].get("ops", []) if o.get("op") == "polygon"]
    regs = {}
    for o in ops:
        regs.setdefault(o["name"], []).append(o["polygon"])
    png = hetissue.fig_view(s.sample_id, box_px, s, regs)
    return [_j({"path": png, "box_px": box_px}), Image(path=png)]


def _he_outline(s, a):
    from . import hetissue
    if a.get("polygon_px"):
        outer = [[float(x), float(y)] for x, y in a["polygon_px"]]
        if len(outer) < 3:
            raise ValueError("polygon_px needs at least 3 points")
        return [{"outer": outer, "holes": []}], "points read on the H&E"
    if not a.get("classes"):
        raise ValueError("give classes (from segment_he) or polygon_px")
    polys = hetissue.outline(s.sample_id, a["classes"], a.get("box_px"), a.get("within_px"),
                             float(a.get("min_area_um2") or hetissue.MIN_AREA_UM2))
    return polys, f"H&E classes {', '.join(map(str, a['classes']))}" + (" within a box" if a.get("box_px") else "")


def _d_propose_he_region(a):
    from . import hetissue, regions as regmod
    s = _s()
    name = regmod.validate_name(re.sub(r"\s+", "_", str(a["name"]).strip()))
    if a.get("resolved"):
        polys, how = a["resolved"], a.get("how") or "the H&E"
    else:
        polys, how = _he_outline(s, a)
    exists = name in s.region_names()
    return {"title": f"Add region {name}" + (f" ({len(polys)} pieces)" if len(polys) > 1 else "") + " outlined on the H&E",
            "version": s.state()["version"],
            "args": dict(a, name=name, resolved=polys, how=how),  # applied exactly as reviewed
            "details": [*_words(a), f"Outlined by the agent from the image alone ({how}), guided by the tissue's "
                        "priors: check it against the H&E; it is a draft to correct.",
                        f"Adds {len(polys)} polygon(s) " + (f"to the existing region {name}" if exists else
                        f"as region {name}") + "; you can edit its points afterwards, and Undo takes it back."],
            "preview": {"kind": "image", "path": hetissue.fig_outline(s.sample_id, polys, name, s)}}


@mcp.tool()
@_reviewed(_d_propose_he_region)
def propose_he_region(name: str, classes: Optional[list[int]] = None, box_px: Optional[list[float]] = None,
                      within_px: Optional[list[list[float]]] = None, polygon_px: Optional[list[list[float]]] = None,
                      min_area_um2: float = 20000.0, text: str = "", resolved: Optional[list] = None,
                      how: Optional[str] = None) -> str:
    """Propose a region outlined on the H&E (works without cells, e.g. a sample of an H&E only): the pixels of
    classes (from segment_he; e.g. the classes of the mucosa), optionally only within box_px [x0, y0, x1, y1] or the
    polygon within_px (full-resolution pixels; one section of a slide, one side of the tissue), gaps closed, holes
    filled, pieces under min_area_um2 dropped; or polygon_px, points you read on show_he. The user sees the outline
    on the H&E in the review window and Accept adds it as region `name` (editable afterwards). text: the user's
    words. resolved, how: filled in by the app (leave them out). Only when the user asked for the region; say it is a
    draft from the image."""
    s = _s()
    from . import regions as regmod
    pieces = [dict(p, name=regmod.validate_name(re.sub(r"\s+", "_", str(name).strip()))) for p in resolved]
    r = s.import_regions(pieces, source="the H&E", author="agent from the H&E (accepted by the user)")
    return _j({"new_version": r["new_version"], "message": r["message"], "regions": s.regions_info(r["new_version"])["regions"]})


# ---------------- aligning cells to an image without a registration file ----------------
def _sample_id(sample_id):
    if sample_id:
        return sample_id
    if "session" in _active:
        return _s().sample_id
    from . import claudelink
    sid = claudelink.viewer_state().get("sample")
    if not sid:
        raise ValueError("which sample? give sample_id (list_samples)")
    return sid


@mcp.tool()
def alignment_status(sample_id: Optional[str] = None) -> str:
    """Whether the sample's cells sit on an image: the image shown, an image waiting for alignment, whether a
    registration file was given and who placed the cells (file / Claude), the cells' extent in um and the image pixel
    size stored in the file. Default sample: the active session's, else the one open in the app."""
    from . import align
    return _j(align.status(_sample_id(sample_id)))


@mcp.tool()
def search_alignment(sample_id: Optional[str] = None) -> list:
    """Find where the cells lie on the sample's image without a registration file: the tissue on the image is matched
    with the cells' footprint over 8 orientations (rotations, mirrored) and a range of pixel sizes. Returns the best
    candidates (K1, K2, ... with orientation, um per pixel, score) and the figure of K1. Scores near 0.7-0.9 with a
    clear lead are typical of the right one; a slide often holds several sections, so check the place too."""
    from . import align
    sid = _sample_id(sample_id)
    st = align.search(align.Sample(sid))
    png = align.show(align.Sample(sid), st["candidates"][0]["id"], refined=False) if st["candidates"] else None
    out = {k: v for k, v in st.items() if k != "candidates"} | {"candidates": [{k: c[k] for k in (
        "id", "orientation_name", "um_per_px", "score")} for c in st["candidates"]], "figure": png}
    return [_j(out)] + ([Image(path=png)] if png else [])


@mcp.tool()
def show_alignment(candidate: str, sample_id: Optional[str] = None) -> list:
    """Figure of a candidate (refined if refine_alignment ran): the cells' footprint on the whole image and two
    close-ups where the cell centroids should sit on nuclei."""
    from . import align
    png = align.show(align.Sample(_sample_id(sample_id)), candidate)
    return [_j({"figure": png}), Image(path=png)]


@mcp.tool()
def refine_alignment(candidate: str, sample_id: Optional[str] = None) -> str:
    """Bring a candidate to micrometer precision: nuclei on the image vs cell centroids in windows over the tissue,
    an affine fit of the offsets, three rounds (60, 25, 12 um search). A good result ends with median shifts of a
    few um and correlations clearly above 0.2; then look at show_alignment."""
    from . import align
    return _j(align.refine(align.Sample(_sample_id(sample_id)), candidate))


@mcp.tool()
@_reviewed(_d_apply_alignment)
def apply_alignment(candidate: str, sample_id: Optional[str] = None) -> str:
    """Put the cells on the image with this candidate (refined if it was): the app shows them there at once and
    outlines drawn before move along. The user can undo it in the app (Import, bottom part: Undo last alignment)."""
    from . import align
    sid = _sample_id(sample_id)
    r = align.apply(align.Sample(sid), candidate)
    if "session" in _active and _active["session"].sample_id == sid:
        _active["session"] = Session(_active["session"].id)
    return _j(r)


# ---------------- differential expression and composition ----------------
def _cells(scope, uses_selection=False):
    from . import analysis, claudelink, projects
    s = _s()
    sel = claudelink.viewer_selection() if uses_selection else None
    if scope == "project":
        pid = projects.project_of(s.sample_id)
        if not pid:
            raise ValueError("this sample is not part of a project: use scope='slice'")
        return analysis.cells(pid=pid, selection=sel)
    if scope != "slice":
        raise ValueError("scope must be 'slice' or 'project'")
    return analysis.cells(s, selection=(sel or {}).get(s.sample_id) if sel else None)


def _uses_selection(*names):
    return any(n == "selection" for n in names) or any(isinstance(n, dict) and "selection" in n for n in names)


@mcp.tool()
def analysis_groupings(scope: str = "slice") -> str:
    """What the cells can be grouped by for de_genes, group_markers and composition, with the values of each (largest
    first): cell type, lineage, region, niche, leiden, selection (the cells the user selected in the viewer when
    they last ran an analysis there), slice, and the slices' metadata columns that differ between slices.
    scope = "slice" (the active annotation) or "project" (the latest annotation of every slice of its project)."""
    c = _cells(scope, uses_selection=True)
    opts = c.options()
    return _j({"scope": scope, "n_cells": len(c.labels), "groupings": {k: v[:40] for k, v in opts.items()},
               "missing": c.missing})


@mcp.tool()
def de_genes(by: str, group: Union[str, list], reference: Union[str, list] = "rest", within: Optional[dict] = None,
             scope: str = "slice", top: int = 20) -> str:
    """Differentially expressed genes of a group of cells: group = a value (or list) of grouping by, reference =
    "rest" (every other cell with a value of by) or another value (or list); within = {grouping: value(s)} limits
    both sides (e.g. by="region", group="mucosa", within={"cell type": "Fibroblasts"}). Wilcoxon test on
    log-normalized expression; returns the top up and down genes ranked by AUC, with log2 fold change, detection
    rates, FDR and, with several slices, the per-slice check. group="all": every group of by against the other
    cells in one run (groups under 30 cells skipped; the test on up to 5,000 cells per side), the top genes per
    group. scope as in analysis_groupings."""
    from . import analysis
    c = _cells(scope, _uses_selection(by, within))
    return _j(analysis.de(c, by, group, reference, within, max(1, min(int(top), 100))))


@mcp.tool()
def trajectory(method: str = "spatial", within: Optional[dict] = None, start: Optional[dict] = None,
               end: Optional[dict] = None, root: Optional[dict] = None, alpha: float = 0.3, by: str = "cell type",
               scope: str = "slice") -> str:
    """Order cells along a trajectory and find the genes that change along it.
    method "spatial": each cell's place on an axis in the tissue between start and end = {"by": grouping, "values":
    [...]} (e.g. start cell type Stem cells, end Colonocytes for the crypt axis; regions; niches): t = d_start /
    (d_start + d_end), 0..1; without end, the distance to start in um. "expression": diffusion pseudotime from root
    ({"by", "values"}, among the cells ordered) on the PCA neighbor graph; "mixed": the same with tissue neighbors
    joined in, alpha = the weight of space (0..1). within = {grouping: value(s)} orders only those cells (e.g.
    {"lineage": "Epithelial"}). Returns the rising, falling and transient genes (r = correlation with the order,
    eta2 = share of variance between bins, peak), the composition of `by` along it and a heatmap of bin means. The
    result also becomes the map's "trajectory: <method>" color in the app (each method keeps its own last run). A spatial axis or a root forces the order: say
    which; pooled slices are not batch-corrected."""
    from . import trajectory as tr
    c = _cells(scope, _uses_selection(within))
    r, _, t, raw, masks = tr.run(c, method, within, start, end, root, alpha, by)
    import uuid as _u
    run_id = _u.uuid4().hex[:12]
    i = 0
    for s in c.sessions:
        n = len(s.cells)
        tr.save(s, t[i:i + n], raw[i:i + n], r, run_id, {k: v[i:i + n] for k, v in masks.items() if v is not None},
                {"scope": scope})
        i += n
    r["label"] = tr.label(r)
    r["heatmap"] = {"genes": r["heatmap"]["genes"]}  # the bin means stay in the app
    return _j(r)


@mcp.tool()
def group_markers(by: str, top: int = 8, within: Optional[dict] = None, scope: str = "slice") -> str:
    """The top genes (by AUC) of every group of grouping by against the other cells, e.g. by="niche" for what sets
    each niche apart, by="leiden" for the clusters; groups under 30 cells are skipped."""
    from . import analysis
    c = _cells(scope, _uses_selection(by, within))
    return _j(analysis.markers(c, by, max(1, min(int(top), 30)), within))


@mcp.tool()
def composition(rows: str, cols: str, within: Optional[dict] = None, scope: str = "slice") -> str:
    """How the cells of one grouping split over another: per value of rows its cell count, the share of each value of
    cols and the log2 enrichment over independence, e.g. rows="niche", cols="cell type"; rows="region",
    cols="niche"; scope="project", rows="slice" (or a metadata column such as condition, which also gives every
    slice's shares), cols="niche"."""
    from . import analysis
    c = _cells(scope, _uses_selection(rows, cols, within))
    return _j(analysis.crosstab(c, rows, cols, within))


@mcp.tool()
def cell_communication(by: str = "cell type", within: Optional[dict] = None, scope: str = "slice",
                       resource: str = "all", contact_um: float = 15, secreted_um: float = 50, min_pct: float = 0.1,
                       nonprotein: bool = False, top: int = 30) -> str:
    """Cell-cell communication between the groups of grouping by (cells inside within): ligand-receptor pairs whose
    genes are in the panel (resource "all" = CellChatDB v2 + CellPhoneDB v5 + OmniPath consensus, or one of them),
    scored on actual neighbor pairs: contact signals within contact_um, secreted and ECM signals within secreted_um
    (Gaussian weights). Returns the top significant sender -> receiver interactions (score, specificity over random
    labels, colocalization of expressing cells, detection, neighbor pairs, FDR), the sender x receiver network
    (count and summed score), the leading pathways, and how many pairs the panel covers."""
    from . import communication
    c = _cells(scope, _uses_selection(by, within))
    r = communication.run(c, by, within, resource, contact_um, secreted_um, min_pct, nonprotein, max(1, min(int(top), 100)))
    r.pop("table")
    return _j(r)


@mcp.tool()
def compare_communication(split: str, by: str = "cell type", within: Optional[dict] = None, scope: str = "project",
                          resource: str = "all", top: int = 30) -> str:
    """The communication analysis in every group of split (e.g. a metadata column such as condition, slice, region,
    niche): the interactions significant in at least one group with their score in each, largest differences first."""
    from . import communication
    c = _cells(scope, _uses_selection(split, by, within))
    r = communication.compare(c, split, by, within, max(1, min(int(top), 100)), resource=resource)
    r.pop("table")
    return _j(r)


@mcp.tool()
def show_communication(ligand: str, receptor: str, side: str = "both") -> list:
    """Whole sample: the H&E (left) and every cell colored by its part in one ligand-receptor signal. side="both":
    cells that send it (their ligand times the receptor around them), receive it (their receptor times the ligand
    around them), both, or only weakly (strong = above the median of the cells with that signal); side="receiver" or
    "sender": percentiles of that side. ligand / receptor are gene symbols (a complex's subunits joined by "_", as
    cell_communication lists them with "+")."""
    from . import communication, plausibility
    s = _s()
    lig, rec = ligand.replace("+", "_"), receptor.replace("+", "_")
    a = communication.activity(s, lig, rec, side)
    if side == "both":
        b, table = communication.role_bins(*a), communication.ROLE_BINS
    else:
        b, table = communication.bins(a), communication.BINS
    names = np.array([x[1] for x in table], dtype=object)
    pal = {x[1]: x[2] for x in table}
    out = s.dir / "figures" / f"ccc_{lig}_{rec}_{side}.jpg"
    out.parent.mkdir(exist_ok=True)
    png = render.tissue_overview(s.meta.get("he_path"), s.cells.px.to_numpy(), s.cells.py.to_numpy(), names[b], pal,
                                 plausibility.grid(s), [], out,
                                 title=f"{lig.replace('_', '+')} -> {rec.replace('_', '+')}",
                                 image_kind=s.meta.get("image_kind", "H&E"),
                                 colored_by="sender / receiver role" if side == "both" else f"{side} signal")
    return [_j({"path": png, "cells": {x[1]: int((b == i).sum()) for i, x in enumerate(table)}}), Image(path=png)]


# ---------------- niches ----------------
@mcp.tool()
def niche_context() -> str:
    """For choosing the niche neighborhood size k: per drawn region (and the cells in none), cell density, the
    median distance covered by k = 10, 15, 20, 30, 50 neighbors, and the main lineages."""
    from . import niches
    return _j(niches.context(_s()))


@mcp.tool()
@_reviewed(_d_run_niches)
def run_niches(k: int = 15, n_niches: int = 8, identity: str = "nmf", n_programs: int = 20,
               k_by_region: Optional[dict] = None, radius_um: Optional[float] = None) -> str:
    """Niche analysis: every cell's identity (identity="nmf": its weights on n_programs NMF gene programs, scaled to
    sum 1; or "cell types": its label) is averaged over the cell and its k nearest cells, and the averages are
    clustered (k-means) into n_niches niches. k_by_region = {"<drawn region>": k} sets k per region (cells in no
    region use k); radius_um instead takes every cell within that distance (5-80 cells), adapting to density. The
    result replaces the session's current niches and shows in the viewer (color by niche). Returns niche_summary."""
    from . import niches
    return _j(niches.run(_s(), k, n_niches, identity, n_programs, k_by_region, radius_um))


@mcp.tool()
def niche_summary() -> str:
    """The current niches: per niche its size, cell types, identity profile (NMF programs with their top genes),
    regions and distance to the epithelium, plus the parameters used."""
    from . import niches
    return _j(niches.summary(_s()) or {"error": "no niches yet: run_niches first"})


@mcp.tool()
def show_niches() -> list:
    """Whole sample: the H&E (left) and the cells colored by niche (right)."""
    from . import niches, plausibility
    s = _s()
    meta, lab = niches.current(s)
    if meta is None:
        raise ValueError("no niches yet: run_niches first")
    names = np.array(niches.niche_names(s), dtype=object)
    pal = {n: tuple(int(c[i:i + 2], 16) for i in (1, 3, 5)) for n, c in zip(names, (render.DISTINCT * 3)[:len(names)])}
    out = s.dir / "figures" / f"niches_{meta['run_id']}.jpg"
    png = render.tissue_overview(s.meta.get("he_path"), s.cells.px.to_numpy(), s.cells.py.to_numpy(), names[lab], pal,
                                 plausibility.grid(s), [], out, title=f"niches ({meta['identity']}, k={meta['k']})",
                                 image_kind=s.meta.get("image_kind", "H&E"), colored_by="niche")
    return [_j({"path": png}), Image(path=png)]


@mcp.tool()
@_reviewed(_d_name_niches)
def name_niches(names: dict) -> str:
    """Give niches names the user sees in the legend, e.g. {"N1": "lamina propria", "N3": "crypt base"}; base them
    on niche_summary and the H&E."""
    from . import niches
    return _j(niches.set_names(_s(), names))


@mcp.tool()
def type_confusion(top: int = 15) -> str:
    """Pairs of cell types the expression model hardly tells apart (mean probability of B in cells labelled A and the
    reverse): candidates to merge, or to separate with sharper markers."""
    from . import insights
    return _j(insights.type_confusion(_s(), top))


@mcp.tool()
def type_structure(cell_type: str, resolution: float = 0.5) -> str:
    """Sub-clusters of one cell type with the genes that set each apart, their distance to the epithelium, regions
    and spatial clustering: shows a type that should be split (distinct subtypes, or another type mislabelled, e.g.
    PDGFRA+ C7+ fibroblasts inside 'smooth muscle')."""
    from . import insights
    return _j(insights.type_structure(_s(), cell_type, resolution))


@mcp.tool()
def spatial_profile() -> str:
    """Per cell type: distance to the epithelium and to the smooth muscle (quartiles, um), regions, and the types it
    neighbors more or less than expected. Compare with the expected organization of the tissue (crypt base to top,
    lamina propria, follicles, muscle layers): contradictions point at mislabels, consistent patterns at priors."""
    from . import insights
    return _j(insights.spatial_profile(_s()))


@mcp.tool()
def qc_overview() -> list:
    """Technical quality per tile (transcripts per cell, cell density, Unassigned share, mixed-lineage share) with
    hotspot tiles, and epithelial-marker spillover into nearby non-epithelial cells. Returns the tile maps."""
    from . import insights
    s = _s()
    png = insights.fig_qc(s)
    return [_j(insights.qc_overview(s) | {"path": png}), Image(path=png)]


@mcp.tool()
def check_registration() -> str:
    """Offset between the nuclei on the H&E (or DAPI) and the cell centroids, per window and as an affine fit (scale,
    rotation, offset at the corners). A max offset above ~8 um with a small fit residual means the cells sit off their
    nuclei: propose the correction (draft_advice scope "annotation", category "qc", action {"kind": "registration"})."""
    from . import insights
    r = insights.registration(_s())
    if r.get("affine_fit"):
        r["affine_fit"] = {k: v for k, v in r["affine_fit"].items() if k != "correction_px"}
    return _j(r)


@mcp.tool()
def check_registration_map(focus: Optional[dict] = None) -> list:
    """The cells on their nuclei over all the tissue (Agent Review: Spatial organization): the nuclei on the H&E (or
    DAPI) against the cell centroids in 200 um windows everywhere, the part one registration fixes (affine: scale,
    rotation, shift) and areas off their nuclei (L1, L2, ...: tissue stretched, folded or torn between the run and
    the staining), with the structures and regions in them; offsets by structure and by region; what the local
    correction would do. Shifts beyond the 25 um the 200 um windows look (up to 100 um) are found in 600 um windows
    (areas with "large_shift"; dashed in the figure; a crypt-sized shift there fits only where the match is
    unique); the correction moves those first, then measures again and corrects the rest. Figure: windows off their nuclei on the section, close-ups of the largest areas (cells now
    cyan, corrected yellow: yellow should sit on the nuclei).
    focus: structures to check first, as {name: [types or lineages from get_vocabulary]}, e.g. {"crypts":
    ["Epithelial"], "vessels": ["Endothelial", "Pericytes"]}; by default epithelium, vessels and smooth muscle by
    their names. Crypts and vessels show a shift first (their nuclei form a clear pattern).
    Judge every area on its close-up: real when the cyan cells sit beside the nuclei (off the tissue edge, across a
    lumen) and the yellow ones on them; a false alarm when cyan already sits on nuclei. If real areas remain, propose
    the local correction (draft_advice scope "organization", category "registration", action {"kind":
    "registration", "mode": "local"}, naming the areas and structures); a verdict "global" alone is check_registration's
    correction ({"kind": "registration"}). The first call on a placement takes up to a minute or two."""
    from . import regcheck
    s = _s()
    info = regcheck.summary(s, focus)
    if info.get("error"):
        return [_j(info)]
    png = regcheck.figure(s, info)
    return [_j(info | {"path": png}), Image(path=png)]


def _d_move_piece(a):
    from . import piecefit
    s = _s()
    p = piecefit.get(s.sample_id, a["piece"])
    if not p["verdict"].startswith("fits"):
        raise ValueError(f"piece {p['id']} does not fit ({p['verdict']}): tell the user instead of proposing it")
    how = f"rotated {p['angle_deg']:.0f} degrees" + (", mirrored (folded over)" if p["mirrored"] else "") + \
        f", shifted {float(np.hypot(*p['shift_um'])):.0f} um"
    return {"title": f"Move the {p['n_cells']:,} cells of region {p['region']} onto their torn piece of tissue ({how})",
            "sample_id": s.sample_id,
            "details": [*_words(a), f"Fit {p['id']}: correlation with the image {p['corr_now']:.2f} where the cells are, "
                        f"{p['corr']:.2f} after the move" + ("" if "another placement" in p["verdict"] else
                        f"; the best other placement scores {p['second']:.0%} of it") + f" ({p['verdict']}).",
                        f"The cells move {p['moved_um_median']:.0f} um (median, up to {p['moved_um_max']:.0f} um) on "
                        "the image; their expression, cell types and regions stay (the outline stays where it was "
                        "drawn). The figure: yellow should sit on the nuclei of the piece.",
                        "Undo: ask the agent (undo_piece_move), or undo correction in Review > Registration (back "
                        "to the registration as imported)."],
            "preview": {"kind": "image", "path": piecefit.figure(s, p)}}


def _d_undo_piece_move(a):
    from . import importer
    s = _s()
    last = importer.last_piece_move(s.sample_id)
    return {"title": f"Put the {last['n_cells']:,} cells of region {last['region']} back where they were before their "
                     "piece was moved", "sample_id": s.sample_id,
            "details": [*_words(a), "Undoes the last move of a torn piece (the cells only move on the image)."]}


@mcp.tool()
def fit_piece(region: str, max_shift_um: float = 300.0, mirror: bool = True) -> list:
    """A torn or folded piece of tissue: the rigid move (rotation, shift; mirrored too when mirror: a piece folded over
    shows its back) that puts the cells of `region` (the user's outline around the piece's cells) onto the image within
    max_shift_um, scored by the nuclei and the tissue no other cell covers (nothing changes). Returns the fit (id P1,
    P2, ...): angle, shift, correlation now vs after, the best other placement as a share of it (crypts repeat), and
    a verdict ("fits", "no better than where the cells are", "weak", "not unique"; "at the edge of the search": try a
    larger max_shift_um). Figure: the piece now (cyan) and moved (yellow); right, the moved cells on the H&E: they
    should sit on the nuclei of the torn piece. Then propose it with move_piece, or explain why it does not fit."""
    from . import piecefit
    s = _s()
    p = piecefit.fit(s, region, max_shift_um, mirror)
    png = piecefit.figure(s, p)
    return [_j(piecefit.public(p) | {"path": png}), Image(path=png)]


@mcp.tool()
@_reviewed(_d_move_piece)
def move_piece(piece: str, text: str = "") -> str:
    """Propose moving the cells of a fitted piece (fit_piece id, e.g. "P1") as fitted. The user sees the fit and its
    figure in the app's review window and accepts or dismisses it. text: the user's words."""
    from . import importer, piecefit
    s = _s()
    r = importer.move_piece(s.sample_id, piecefit.get(s.sample_id, piece))
    _active["session"] = Session(s.id)
    return _j(r)


@mcp.tool()
@_reviewed(_d_undo_piece_move)
def undo_piece_move(text: str = "") -> str:
    """Propose putting back the cells of the last piece moved (move_piece); the user accepts it in the app."""
    from . import importer
    s = _s()
    r = importer.undo_piece_move(s.sample_id)
    _active["session"] = Session(s.id)
    return _j(r)


@mcp.tool()
def draft_advice(scope: str, items: list[dict]) -> str:
    """Your draft advice on one part of the analysis (replaces your earlier drafts on it; at most 12). Drafts are not
    shown to the user: check each with check_advice first. scope and its categories:
      annotation: taxonomy | context | qc;   organization: compartments | organization | registration;
      embedding: clusters | integration | parameters;   niches: niches | parameters;
      compare: findings | design | next;   signaling: interactions | design | next.
    Each item: {"category", "title": short, "text": the finding with its evidence and what to do (2-4 sentences),
    "action": optional, one of
      {"kind": "region", "compartment": <name from preview_compartment>, "name": region name, "prior": optional spec}
      {"kind": "prior", "prior": {"kind", "text", "types", "mode", "weight", "regions" | "max_depth_um" |
                                  "neighbor_types" + "min_frac"}}   (as for preview_constraint)
      {"kind": "marker_list"}   (after propose_marker_list)
      {"kind": "registration"}  (after check_registration)
      {"kind": "registration", "mode": "local"}  (after check_registration_map: also the local shifts)
      {"kind": "niche_names", "names": {"N1": "...", ...}}
      {"kind": "show", "mode": "gene:<gene>" | "lr:<ligand>~<receptor>~receiver" | "niche" | "leiden" | "type"}
      {"kind": "rerun", "tab": "embed" | "niches" | "compare" | "signal", "params": the tab's settings, e.g.
         embed {"resolution", "min_dist", "n_pcs", "n_neighbors"}, niches {"k", "n_niches", "identity"},
         compare {"kind": "de" | "markers" | "composition", "by", "group", "reference", "rows", "cols", "within"},
         signal {"kind": "run" | "compare", "by", "split", "within", "contact_um", "secreted_um"}}
      {"kind": "claude", "request": what you would do when the user asks, e.g. a re-clustering}}."""
    from . import advice
    new = advice.draft(_s(), scope, items)
    return _j({"drafts": [{k: x[k] for k in ("id", "category", "title")} for x in new],
               "next": "check every draft against the data (check_advice): keep only what holds up"})


@mcp.tool()
def check_advice(verdicts: list[dict]) -> str:
    """Your check of your drafts, as a skeptical reviewer: [{"id", "keep": true/false, "reason": one line,
    "evidence": what you checked (tool, numbers, figure) for a kept item}]. Before keeping a draft, test it: rerun
    the tool that would confirm or refute it, look at the figure on the H&E, and weigh other explanations (spillover,
    segmentation, slice or batch effects, too few cells, panel limits, broadly expressed genes). Keep only advice
    that holds up and changes what the user should do. Kept items replace your earlier open advice on that part and
    appear in the Review panel; dropped ones are folded away there with your reason. The user acts with one click;
    nothing changes before."""
    from . import advice
    return _j(advice.check(_s(), verdicts))


@mcp.tool()
def list_advice(scope: Optional[str] = None, status: Optional[str] = None) -> str:
    """Your advice and what the user did with it (status draft | open | filtered | done | dismissed), for one scope
    (annotation, organization, embedding, niches, compare, signaling) or all."""
    from . import advice
    return _j([{k: v for k, v in x.items() if k != "action"} | {"action": {k: v for k, v in (x["action"] or {}).items()
                                                                      if k not in ("polygons_px", "correction_px")}}
               for x in advice.items(_s(), status, scope)])


@mcp.tool()
def recent_analyses(section: Optional[str] = None) -> str:
    """The analyses the user ran in the app on this annotation, newest last: section "compare" (DE genes, top genes,
    composition), "signaling" (cell-cell communication) or "trajectory", each with its settings and headline result.
    Rerun one with de_genes, group_markers, composition, cell_communication, compare_communication or trajectory to
    see it in full."""
    from . import analysis
    return _j(analysis.recent(_s(), section))


def _embedding_cells(scope):
    from . import embedding, analysis
    s = _s()
    meta = embedding.current(s.sample_id)
    if meta is None:
        raise ValueError("no PCA + UMAP yet: the user runs it in the Embed tab")
    if scope == "project" or (scope == "auto" and meta.get("project")):
        c = _cells("project")
    else:
        c = analysis.cells(s)
    if "leiden" not in c.labels.columns:
        raise ValueError("no Leiden clusters in these cells")
    return c, meta


@mcp.tool()
def embedding_facts(scope: str = "auto") -> str:
    """The PCA + UMAP of the sample (or of all slices, for a joint embedding; scope "auto" follows how it was made):
    its parameters, every Leiden cluster with its size, purity (share of its main cell type) and types, every type
    with the clusters holding 10% or more of it, mixed clusters, types split over 3+ clusters, and with several slices
    the clusters made almost only of one slice."""
    from . import analysis
    c, meta = _embedding_cells(scope)
    return _j(analysis.embedding_facts(c, meta))


@mcp.tool()
def show_embedding(color: str = "cell type", scope: str = "auto") -> list:
    """The UMAP colored by a grouping (cell type, lineage, leiden, niche, slice, region, ...)."""
    from . import embedding
    c, meta = _embedding_cells(scope)
    if color not in c.labels.columns:
        raise ValueError(f"color must be one of {list(c.labels.columns)}")
    Y = np.vstack([embedding.coords(x.sample_id, "umap") for x in c.sessions])
    s = _s()
    out = s.dir / "figures" / f"umap_{color.replace(' ', '_')}_{'all' if c.pid else 'slice'}.jpg"
    out.parent.mkdir(exist_ok=True)
    png = render.embedding_plot(Y, c.labels[color].to_numpy(), out,
                                title=f"UMAP ({'all slices' if c.pid else s.sample_id}), colored by {color}")
    return [_j({"path": png, "n_cells": int(len(Y))}), Image(path=png)]


@mcp.tool()
def annotation_facts(version: Optional[int] = None) -> str:
    """How the annotation was made (method, marker list or reference, priors, corrections, regions, history): the
    facts for an annotation report."""
    from . import insights
    return _j(insights.annotation_facts(_s(), version))


@mcp.tool()
def write_report(markdown: str) -> str:
    """Save the annotation report (Markdown) the user sees under Report in the viewer and can save as a file. Write it
    for a reader who was not there: sample; method and its inputs (annotation_facts); composition; compartments and
    spatial organization; quality and limitations (QC, registration, open review items); what the user changed and
    why; open questions. Facts from the tools only; cite numbers."""
    from . import advice
    return _j(advice.write_report(_s(), markdown))


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
@_reviewed(_d_record_verdicts)
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
@_reviewed(_d_remove_correction)
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
def organ_marker_draft(organ: str, sample_id: Optional[str] = None) -> str:
    """The built-in marker list of an organ, pruned to a sample's panel (default: the active session's sample): a
    marker list CSV with the cell types the panel can resolve (two markers in the panel, or one decisive marker), and
    the types left out with their markers (in_panel: those the panel has). organ: one of the built-in organs (listed
    when the name is unknown), or "other" for immune and stromal cells only. A starting point: add the tissue's own
    cell types the panel can resolve, check_marker_list, then save_marker_list."""
    from . import organs
    from .session import panel_genes
    sid = sample_id or (_active["session"].sample_id if "session" in _active else None)
    try:
        d = organs.draft(organs.safe_name(organ), panel_genes(sid))
    except ValueError as e:
        return _j({"error": str(e), "organs": organs.organs()})
    return _j(d)


@mcp.tool()
def inspect_reference(path: str, sample_id: Optional[str] = None) -> str:
    """A single-cell reference h5ad for start_session(reference_h5ad=...): number of cells, how many of the sample's
    panel genes it has (and where the names come from), and the obs columns that could hold the cell types."""
    from . import scref
    from .session import panel_genes
    return _j(scref.inspect(path, panel_genes(sample_id)))


@mcp.tool()
def inspect_annotation(path: str, sample_id: str) -> str:
    """An existing annotation CSV for start_session(annotation_csv=...): its columns (with examples), the guessed
    cell id / cell type / lineage columns, and how many of the sample's cells it covers."""
    from . import existing
    from .session import sample_cell_ids
    return _j(existing.inspect(path, sample_cell_ids(sample_id)))


@mcp.tool()
def save_marker_list(name: str, text: str, sample_id: Optional[str] = None) -> str:
    """Save a marker list (validated against the sample's panel) in the user's library, where the viewer's
    "Markers..." dialog lists it ("load a list..."). Use it for a new panel: draft a list from get_panel_genes,
    check it with check_marker_list, save it, and tell the user to load it and click Annotate."""
    from . import config, markers
    from .session import panel_genes, filtered_genes
    import re
    p = markers.parse(text, panel_genes(sample_id), filtered_genes(sample_id))
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
    panel), warnings (also for genes the sample's filter removed: they are left out), and the parsed outline."""
    from .session import panel_genes, filtered_genes
    from . import markers
    sid = sample_id or (_active["session"].sample_id if "session" in _active else None)
    p = markers.parse(text, panel_genes(sid), filtered_genes(sid))
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
@_reviewed(_d_label_cluster)
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
    Show the user the number of cells whose marker label changes and the new clusters; recluster then proposes it
    (the user accepts it in the app's review window)."""
    return _j(_s().preview_recluster(target, resolution, marker_edits))


@mcp.tool()
@_reviewed(_d_recluster)
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
@_reviewed(_d_resolve_rule_suggestion)
def resolve_rule_suggestion(suggestion_id: str, status: str, note: str = "") -> str:
    """Mark a reviewer's rule suggestion as accepted or rejected (or reopen it), with a short note on what was done."""
    _s().resolve_suggestion(suggestion_id, status, note)
    return _j({"id": suggestion_id, "status": status})


@mcp.tool()
def pending_changes(proposal_id: Optional[str] = None) -> str:
    """The data changes you proposed and what the user did with them: open (waiting in the app's review window),
    applying, done (accepted and applied, with the result), failed, or dismissed. proposal_id: one of them; without
    it the last 10."""
    from . import proposals
    if proposal_id:
        return _proposal_outcome(proposals.get(proposal_id))
    return _j([{k: r.get(k) for k in ("id", "tool", "title", "status", "error")} for r in proposals.all(limit=10)])


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
@_reviewed(_d_revert)
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


@mcp.tool()
def export_regions(version: Optional[int] = None) -> str:
    """Write the region outlines of a version as GeoJSON (the image's full-resolution pixels, named by region, as QuPath
    writes annotations) into the session's exports; returns its path. It is what a later import of the data on the same
    image takes as its regions file."""
    from . import regions as regmod, importer
    s = _s()
    v = version or s.versions()[-1]
    m = importer._meta(s.sample_id)
    size = importer._he_size(m["he_path"]) if m.get("he_path") and m.get("image_kind") not in (None, "none") else None
    text, info = regmod.to_geojson(s.state(v)["regions"], s.region_colors(), size)
    p = s.dir / "exports" / f"{s.sample_id}_v{v}_regions.geojson"
    p.parent.mkdir(exist_ok=True)
    p.write_text(text, encoding="utf-8")
    return _j({"path": str(p)} | info)


@mcp.tool()
def export_notebook(version: Optional[int] = None) -> str:
    """Write the analysis as a Jupyter notebook: the code of every step with this annotation's parameters (data,
    placement on the image, annotation replayed with its regions, priors and corrections, PCA / UMAP / Leiden,
    niches) and its results; running it repeats the analysis and checks the labels. Returns its path."""
    from . import notebook
    return _j({"notebook": notebook.build(_s(), version)})


def main():
    _active["via"] = "Claude Code"
    _heartbeat()
    mcp.run()


if __name__ == "__main__":
    main()
