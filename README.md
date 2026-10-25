# celltypeagent

Interactive cell type annotation for spatial transcriptomics (Xenium) on the H&E image. Cells are annotated
from a hierarchical marker gene list (no labeled reference needed). You then refine the result in plain language
with Claude: spatial priors ("epithelial cells are only in the mucosa"), corrections ("these cells are plasma
cells"), or requests to split a mixed cluster. Claude turns each request into an explicit constraint, correction
or re-clustering step; the probabilistic model does the annotation, flags conflicts for review, and versions
every step. Everything runs on your own computer; images and expression data never leave it.

## 1. Install

Requirements: a Mac with Apple silicon (macOS 11+), about 8 GB of RAM, and [Claude Code](https://claude.com/claude-code)
logged in with your own Claude account. (Linux and Windows can run it from source.)

**The app (recommended).** Open `SpatialAnnotator.dmg` and drag **Spatial Annotator** to Applications. It
contains everything it needs; no Python installation. If the build is not notarized, the first time right-click
the app and choose Open (or allow it in System Settings -> Privacy & Security).

**From source** (developers, Linux/Windows, or before a .dmg exists):

1. Install [uv](https://docs.astral.sh/uv/) (it also installs Python for you), then load it into the current
   terminal:

   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   source $HOME/.local/bin/env
   ```

2. Unzip the downloaded `STCellTypeAgent-main.zip` and move the folder somewhere permanent (the app and Claude
   are pointed at it), then install it with the desktop window (`[app]`):

   ```bash
   mkdir -p ~/Applications && mv ~/Downloads/STCellTypeAgent-main ~/Applications/
   cd ~/Applications/STCellTypeAgent-main && uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -e ".[app]"
   ```

   On Windows use `.venv\Scripts\...` in place of `.venv/bin/...`. Add `gpu` (`".[app,gpu]"`) to install PyTorch,
   which only speeds up the marker step on an NVIDIA or Apple GPU.

3. macOS: `.venv/bin/celltypeagent make-app` puts a **Spatial Annotator** launcher in `~/Applications` (it runs this
   folder's environment). Elsewhere: `.venv/bin/celltypeagent app` (`--browser` for the browser instead of a window).

Your data and annotations are kept in one data folder, `~/celltypeagent` by default (`CELLTYPEAGENT_HOME` sets
another one).

## 2. Import a sample

Click **Import…** in the app (or use the commands below). Two kinds of input:

| input | what it is |
|---|---|
| Xenium output folder | the folder written by the Xenium instrument / Xenium Ranger: `cells.parquet` or `cells.csv.gz` and `cell_feature_matrix.h5` (or `cell_feature_matrix/`); control probes are dropped. The DAPI image (`morphology_focus/`) is the background unless you add an H&E with its Xenium Explorer alignment file (Image alignment -> export). |
| h5ad + H&E + registration | AnnData file(s) with raw counts (layer `counts`, or `X`) and cell centroids in micrometers (`obs` columns `x_centroid`/`y_centroid` or `x`/`y`, or `obsm["spatial"]`); several files (e.g. one per lineage) are merged and cell type labels in them are ignored. The H&E as a pyramidal TIFF / OME-TIFF. The registration either as an affine CSV with columns `a11,a12,a21,a22,tx,ty` (`px = a11*x + a12*y + tx`, `py = a21*x + a22*y + ty`; cell um to H&E full-resolution pixels; with several rows, the row whose other columns contain the section value is used), or as a Xenium Explorer alignment matrix (3x3 CSV, H&E pixels to morphology pixels). |

In the dialog, the h5ad's columns are listed so you can pick the section (when a file holds several) and an
optional region annotation (regions start empty; Claude can start them from this column on request). The import checks that at least 90% of the cells land on the image (a
wrong registration or section fails here) and takes 10-30 s. From the command line:

```bash
.venv/bin/celltypeagent import-xenium --name ln1 --xenium /path/to/output-XETG... [--he he.ome.tif --registration he_imagealignment.csv]
.venv/bin/celltypeagent import --name <sample_id> --h5ad cells.h5ad --he he.ome.tif --registration alignment.csv \
  [--section-column donor_id --section 147328] [--region-column histopathology]
```

Then give the sample a marker list (there is no default): **Markers…** -> **Import…** a CSV (or text) file, or
ask Claude Code "draft a marker list for this panel" and load it under "load a list…". Imported and saved lists
stay in that list for the next samples. The CSV has one row per cell type:

| cell_type | parent | positive_markers | negative_markers |
|---|---|---|---|
| Epithelial | | EPCAM KRT8 | |
| Absorptive | Epithelial | CA1\|SLC26A3\|FABP2\|AQP8 | SPINK4 |
| Goblet cells | Epithelial | SPINK4 RETNLB ITLN1 | |

`parent` is empty for lineages; markers are separated by spaces (or semicolons); `A|B` means any one of A or B.
**Export CSV** in the same dialog saves the list being edited.

### Example: the IBD Xenium atlas (SCP2927)

Data of Pokatayev et al., Nature 2026. Download the five `scp_xenium_human_*.h5ad` files from Single Cell Portal
study [SCP2927](https://singlecell.broadinstitute.org/single_cell/study/SCP2927), and the H&E image of the slide
from [Zenodo record 17518435](https://zenodo.org/records/17518435) (`examples/scp2927/he_files.tsv` lists the
eight images with sizes, MD5 checksums and download links; 1.4-2.1 GB each). Each slide holds two sections
(`human_slideid-<slide>-<donor>-<donor>`). The registration file is shipped here,
`examples/scp2927/registration.csv` (16 sections; columns `slide`, `donor`):

```bash
D=~/Downloads/SCP2927
.venv/bin/celltypeagent import --name 147328 \
  --h5ad $D/scp_xenium_human_{epi,fib,lymph,myeloid,stroma}.h5ad \
  --he $D/human_slideid-0050423-108731-147328.ome.tif \
  --registration examples/scp2927/registration.csv \
  --section-column donor_id --section 147328 --region-column histopathology
```

Marker list for the 480-gene IBD panel: `examples/scp2927/ibd_xenium480_markers.csv` (import it in **Markers…**).

## 3. Work with Claude

Claude runs in your own Claude Code, with your own Claude account; the app connects it and shows every change.
Click **Claude Code** in the app: Terminal opens with Claude Code connected to the app and to the sample you are
looking at (install Claude Code and log in once beforehand).

Claude loads what you are looking at (`use_viewer_session`), and when it starts a new annotation (for example with
a revised marker list) the app switches to it by itself. Examples:

- "Review my annotation." (markers and layout; suggestions appear in the app's Review panel)
- "Draft a marker list for this panel." / "Epithelial cells are only found in the mucosa." (priors are previewed;
  say go ahead to apply)
- "Cluster 7 of the T cells is mixed, split it." / "These cells in region crypt_1 are goblet cells."
- "Show me the biggest conflict group." / "Go back to version 3." / "Export the result."

Claude answers in the language you write in.

## 4. The viewer

The app's window (or the web viewer: `.venv/bin/celltypeagent viewer`, then http://localhost:8765). Changes made
by Claude show up within a few seconds.

1. Pick a sample; earlier annotations are under **Session**. Choose the **Method**: *Marker genes* (a marker
   list, **Markers…**), *Single-cell reference* (a labeled reference h5ad, **Reference…**) or *Load existing* (a
   CSV with one cell type per cell, **Annotation file…**), then click
   **Annotate** (about 1 min for 350k cells with markers, seconds with a reference). Regions start empty: you
   draw them.
2. The result is drawn on the H&E: pan and zoom down to single cells; color by cell type, lineage, region,
   changes, conflicts, confidence or clusters. Click a legend entry to show only that type; click a cell for
   its label, alternatives and region. Opacity, dot size, hide cells, hide H&E.
3. **Regions are drawn by you.** In the Regions panel click **Draw region**, click points on the map, close
   with a double-click, Enter or the first point (Backspace removes a point, Esc cancels), and name it. A
   name can have several polygons; regions can be renamed, deleted, zoomed to and outlined. Cells in no
   region are "unknown" and region priors never apply to them. Refer to regions by name when talking to Claude.
   **edit** next to a polygon shows its points: drag a point to move it, drag a small circle between two points to
   add one, double-click a point to remove it, then **Save** (Enter) or **Cancel** (Esc). The region keeps its name
   and priors; its cells are recomputed, and Undo brings the old outline back. Click a region's color square to
   choose one of 15 colors; outlines, cells colored by region and the legend follow, and the color stays with the
   region when it is renamed.
   You can also ask Claude in Claude Code to draw one ("draw the mucosa", "outline the lymphoid follicles", "make
   a region of tiles F3:H6"): Claude computes the outline from where that tissue's cells lie, checks it on the H&E
   and offers it in the Review panel, which opens by itself and shows the outline on the map (other advice stays).
   Click **add region** to keep it, then **edit** to adjust it.
4. Hand-labeled cells are locked and fine-tune the expression model (only the involved cell types can move),
   so similar cells elsewhere update too. Every step is a version; **Undo** steps back one change at a time (clicking it again goes further
   back; it never undoes the undo).
5. **Export CSV** saves one row per cell
   (`cell_id`, `region`, `cell_type`, `lineage`, `level_1`… the full type path, `confidence`,
   `expression_label`, `marker_cluster`, `hand_labeled`, `x_um`, `y_um`).
6. **Reset** removes every annotation and chat of the selected sample (moved to `<data folder>/trash/`).
   **Remove all…** (next to Import, desktop app) deletes every sample with all its annotations and chats for good,
   to import from scratch (it cannot be undone); marker lists stay.
7. **Markers…** edits, checks, saves and loads marker lists (saved lists go to `<data folder>/marker_lists/`).
8. **Review annotation** opens the Review panel: computed marker checks (markers that are weak or unspecific in
   their cells, negative markers that are detected, genes of Unassigned cells) and areas to check (where cells of
   types that form sheets, such as epithelium or smooth muscle, sit apart from their own lineage), outlined on the
   map. Clicking it also asks your Claude Code to check every candidate on the H&E and with marker biology
   (record_review); while Claude checks, unchecked candidates are not listed, and afterwards only the real problems
   are shown, each with Claude's reason (false alarms are folded under "Filtered out by Claude"). Verdicts carry over
   to later versions, so only new candidates are checked again. **solve** removes an item you have dealt with. Asking Claude Code "review my annotation" lets Claude review markers and layout with the
   H&E and its knowledge of the tissue. Its suggestions appear in the panel: a revised marker list (**open in
   Markers…**, then Annotate with "keep my regions, priors & corrections") and areas drawn in orange (**use as region**, then
   Claude applies the agreed prior or correction to that region).

   When the labels changed since Claude's last look, Review annotation also asks Claude for **high-level advice**,
   listed at the top of the panel in five groups, each with one button:
   - *Compartments*: outlines of the mucosa, muscularis or lymphoid follicles, computed from where a lineage forms
     tissue and checked by Claude on the H&E (drawn in purple; **add region**, then **add prior…** for the prior
     that uses it).
   - *Cell types*: types to merge (the model hardly tells them apart), to split (sub-clusters with their own genes
     and place), or missing (**open in Markers…** for a revised list, **ask Claude Code** for a re-clustering).
   - *Spatial organization*: expected gradients and neighbors, found or contradicted (**add prior…**).
   - *Sample context*: what the composition says about site and disease, and what to check.
   - *Technical*: low-quality tiles, marker spillover next to the epithelium, and the registration. The registration
     check correlates the nuclei on the H&E (hematoxylin) or DAPI with the cell centroids in windows across the
     section and fits an affine; **correct registration…** moves every cell onto its nucleus (regions keep their
     cells; **undo correction** puts them back). On SCP2927 section 147328 it found a 0.28% scale mismatch (cells up
     to 30 um off their nuclei at the corners); after the correction the windows agree within 0.5 um.
   Advice changes nothing until you click; **dismiss** hides it.
9. **Niches** (fourth row): places where cells of similar make-up live together. Every cell gets an identity, its
   weights on NMF gene programs (*Identity: NMF programs*, independent of the annotation) or its cell type; the
   identities are averaged over the cell and its *k* nearest cells and clustered (k-means) into the chosen number of
   niches. **Run niches** uses one k everywhere; **Claude picks k** asks Claude Code to choose k per drawn region
   from cell density and cell types (or a radius, so neighborhoods cover a similar stretch of tissue in dense and
   sparse areas), run it and name the niches. The result appears as *Color by: niche* (names in the legend) and as
   a `niche` column in Export CSV. On SCP2927 section 147328 (356,629 cells) a run takes about 6 s (1.5 s once the
   NMF is cached) and the 8 niches match the tissue layers (surface epithelium, crypt base, lamina propria, lymphoid
   follicles, muscularis mucosae, muscularis propria, submucosal vessels and stroma).
10. **Report** shows the annotation report Claude writes (**Ask Claude to write it**): sample, method and inputs,
   composition, compartments, quality and open questions, from the tools' facts. **Save…** writes it as Markdown.

The app has no chat box of its own: you talk to Claude in Claude Code. The web viewer on a
shared server shows a chat box that uses the server's own logged-in Claude Code CLI (for the server owner's use).

## Annotation methods

**Load existing.** A CSV (or TSV) with one row per cell: the cell id as in the sample (h5ad obs names or the
Xenium `cell_id`), the cell type, and optionally a lineage column; the dialog guesses the columns and shows how
many cells match. Its labels become the annotation (cells not in the file are Unassigned); an expression model
trained on them adds a confidence per cell and lets priors and corrections act. A CSV exported by the app can be
loaded back. Example: `share/celltype.csv` holds the published cell types of SCP2927 section 147328
(`cell_id, cell_type, lineage`).

**Single-cell reference** (like CellTypist). Give a labeled single-cell h5ad of the same tissue, its cell type
column, and optionally a coarser group column (used as lineages for colors and priors). A multinomial logistic
regression is trained on up to 3,000 reference cells per label, on the genes the reference shares with the panel,
and labels every cell; both sides are normalized on the shared genes only, so a whole-transcriptome reference
fits a panel of a few hundred genes. Raw counts are read from `layers["counts"]`, `raw.X` or `X` (log-normalized
values are converted back), and gene names from the var index or the var column that matches the panel best
(e.g. `feature_name` next to Ensembl ids). *Majority voting* over-clusters the sample and lets each cluster agree
on one label. Hand corrections fine-tune the classifier (replaying the reference cells). Example: a reference built
from three other SCP2927 sections labeled section 147328 with 80% agreement on 56 fine types and 97% on the five
lineages (98% with majority voting); a reference from another platform (scRNA-seq) will agree less.

**Marker genes.** A hierarchical marker list. Format in text (indent two spaces per level; the CSV above is equivalent):

```
# Name: GENE+ GENE- A|B+    (+ expressed, - not expressed, A|B+ = any of A or B)
Epithelial: EPCAM+ KRT8+
  Absorptive: CA1|SLC26A3+ SPINK4-
    BEST4+ colonocytes: BEST4+ OTOP2+
T and NK cells: CD3E|NKG7+ EPCAM-
  CD4 T cells: CD4+ CD8A-
```

Claude can draft a list from the sample's panel genes and checks it before use. Annotation is top-down: at
each level the parent's cells are re-clustered (PCA, kNN graph, Leiden; the resolution is chosen per level
from 0.5 / 1 / 2 / 4 as the one whose clusters the markers assign best), each cluster is scored against the
child types' markers, and assigned only if the best score and its margin are high enough. **Cell types are
the leaves of the tree** (plus "Unassigned"); a type with subtypes is a lineage, usable in priors but never a
label, so the subtypes of a type should cover all its cells (add a catch-all subtype where needed). The leaf
labels train an expression model (multinomial logistic regression) on the sample's own cells, which labels
the cells the marker step left at a parent level within their lineage. Priors reweight its probabilities:

| kind | condition | example |
|---|---|---|
| `region` | cell lies in one of the named regions | epithelial cells only in the mucosa |
| `depth` | within a distance of the epithelium (estimated from EPCAM/KRT8) | goblet cells within 30 um of the epithelium |
| `neighborhood` | a fraction of the 15 nearest cells are given types | Tfh cells next to B cells |

Each works as `only`, `exclude` or `enrich` with strength `strong` / `moderate` / `weak`. Local re-clustering
("this cluster is mixed") re-clusters one cluster or level, is previewed, and keeps regions, priors and
corrections. To add or remove types, revise the marker list and Annotate again with "keep my regions, priors & corrections".

**Review.** Claude reviews an annotation with three kinds of evidence. *Markers:* for every lineage and cell type,
how often each listed marker is detected in its cells compared with its siblings, flags (weak, not specific,
negative marker detected), and the genes that actually separate those cells in this sample; Unassigned cells are
described the same way. *Space:* a whole-sample overview (H&E with a lettered tile grid next to the lineage map) and
computed candidate areas, i.e. tiles where cells of a type that forms sheets sit apart from their own lineage, at
least 3x more than expected from cell density. *Biology:* Claude's knowledge of the tissue, cell types and
histology, which also tells real structures (vascular smooth muscle, muscularis mucosae, deep crypts in IBD) from
mislabels. Nothing changes until the user agrees. Example (SCP2927 section 146685, ulcerative colitis, scored
afterwards against the dataset's own labels): the review found that 15,514 Unassigned cells expressed CD177,
CSF3R and S100A12 (neutrophils missing the myeloid markers) and that the "dendritic cells" had macrophage genes
(C1QC, LYVE1, MS4A7); with the revised markers, Unassigned cells fell to 195, 97% of the dataset's neutrophils
were labeled neutrophils (15% before), and lineage agreement rose from 91.4% to 98.6%.

## Performance

Everything runs on the CPU; with PyTorch installed (`gpu` extra) the marker step's neighbor search uses an NVIDIA
or Apple GPU. On an idle CPU a 356k-cell section takes about a minute to annotate (3.3 GB of RAM); a correction or
a local re-clustering takes a few seconds. The viewer draws the cells with WebGL, so panning, zooming and changing
colors are immediate.

## Troubleshooting

- The app does not open: see `<data folder>/app.log`; run `.venv/bin/celltypeagent app` in a terminal to see
  errors, or `--browser` to use the browser. The web viewer logs to `<data folder>/viewer.log`.
- Claude Code does not open: it must be installed and logged in (`claude` works in Terminal). By hand:
  `cd ~/celltypeagent/claude_workspace && claude --mcp-config celltypeagent_mcp.json`.
- Import says too few cells land on the H&E: check the registration file, `--section`, and that the cell
  coordinates are in micrometers.

## Other clients

- Claude Code without the app: `claude mcp add celltypeagent -e CELLTYPEAGENT_HOME=$HOME/celltypeagent -- $PWD/.venv/bin/python -m celltypeagent.mcp_server`
  (run in this folder). `.claude/skills/spatial-annotation/` holds the agent's working rules for Claude Code.
- Shared server: the viewer binds to localhost only. Others reach it through an SSH tunnel
  (`ssh -N -L 8765:localhost:8765 <user>@<server>`), then open http://localhost:8765. The conflict review page for
  pathologists is at http://localhost:8765/review (verdict keys 1-6, rules form, CSV download).

## Building the app (maintainers)

On a Mac with Apple silicon, from this folder:

```bash
packaging/macos/build_app.sh            # -> dist/Spatial Annotator.app and dist/SpatialAnnotator.dmg (ad-hoc signed)
```

It builds with PyInstaller (`packaging/spatial_annotator.spec`) in its own environment, without PyTorch. The app's
executable also takes the command-line subcommands (`Spatial Annotator mcp` is what Claude Code starts). To give it a Developer ID signature and Apple's notarization, so that it opens without warnings on other
Macs (needs an Apple Developer account):

```bash
xcrun notarytool store-credentials spatial-annotator --apple-id YOU@EXAMPLE.COM --team-id TEAMID   # once
DEVELOPER_ID="Developer ID Application: Your Name (TEAMID)" packaging/macos/sign_and_notarize.sh
```

The hardened-runtime exceptions it uses (`packaging/macos/entitlements.plist`) allow numba's run-time compilation.

## Layout

```
celltypeagent/     Python package
  cli.py           celltypeagent app / make-app / import / import-xenium / samples / viewer / mcp
  desktop.py       the desktop app (pywebview window around the local server)
  claudelink.py    connection to the user's Claude Code; viewer <-> agent handoff files
  importer.py      Xenium output folder, or h5ad + H&E + registration -> sample bundle
  markers.py       marker list parsing, hierarchical Leiden annotation, re-clustering
  model.py         expression model (training, fine-tuning on corrections)
  constraints.py   priors and Bayesian reweighting
  session.py       sessions, append-only versions, corrections, export
  regions.py       user-drawn regions
  conflicts.py     conflict detection and grouping
  plausibility.py  review of the spatial layout: tile grid, candidate areas, overview figures, area suggestions
  markerqc.py      review of the markers: support per node, flags, top genes, marker list proposals
  mcp_server.py    MCP tools (Claude Code)
  webapp.py        viewer, review page, chat box
  data/            app icon
examples/scp2927/  registration of the SCP2927 sections
packaging/         PyInstaller spec, macOS build / signing scripts
tests/             smoke tests (python tests/smoke_engine.py, tests/smoke_mcp.py; they expect sample 147328 / 167799)
```
