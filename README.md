# STAnnotator

Interactive cell type annotation for spatial transcriptomics (Xenium) on the H&E image. Cells are annotated
from a hierarchical marker gene list (no labeled reference needed). You then refine the result in plain language
with Claude: spatial priors ("epithelial cells are only in the mucosa"), corrections ("these cells are plasma
cells"), or requests to split a mixed cluster. Claude turns each request into an explicit constraint, correction
or re-clustering step; the probabilistic model does the annotation, flags conflicts for review, and versions
every step. Everything runs on your own computer; images and expression data never leave it. (The Python package
and its command are called `stagent`; the older name `celltypeagent` still works.)

## 1. Install

Requirements: a Mac with Apple silicon (macOS 11+), about 8 GB of RAM, and [Claude Code](https://claude.com/claude-code)
logged in with your own Claude account. (Linux and Windows can run it from source.)

**The app (recommended).** Open `STAnnotator.dmg` and drag **STAnnotator** to Applications. It
contains everything it needs; no Python installation. If the build is not notarized, the first time right-click
the app and choose Open (or allow it in System Settings -> Privacy & Security).

**From source** (developers, Linux/Windows, or before a .dmg exists):

1. Install [uv](https://docs.astral.sh/uv/) (it also installs Python for you), then load it into the current
   terminal:

   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   source $HOME/.local/bin/env
   ```

2. Unzip the downloaded `STAnnotator-main.zip` and move the folder somewhere permanent (the app and Claude
   are pointed at it), then install it with the desktop window (`[app]`):

   ```bash
   mkdir -p ~/Applications && mv ~/Downloads/STAnnotator-main ~/Applications/
   cd ~/Applications/STAnnotator-main && uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -e ".[app]"
   ```

   On Windows use `.venv\Scripts\...` in place of `.venv/bin/...`. Add `gpu` (`".[app,gpu]"`) to install PyTorch,
   which only speeds up the marker step on an NVIDIA or Apple GPU.

3. macOS: `.venv/bin/stagent make-app` puts an **STAnnotator** launcher in `~/Applications` (it runs this
   folder's environment). Elsewhere: `.venv/bin/stagent app` (`--browser` for the browser instead of a window).

Your data and annotations are kept in one data folder, `~/stagent` by default (`STAGENT_HOME` sets
another one; a data folder of an older version, `~/celltypeagent`, is used as it is when there is no `~/stagent`).

## 2. Import a sample

Click **Import…** in the app (or use the commands below): the cells, then **No image** or an image file, then a
**Registration file** or **No registration**:

| part | what it takes |
|---|---|
| cells | a Xenium output folder (`cells.parquet` or `cells.csv.gz` and `cell_feature_matrix.h5` or `cell_feature_matrix/`; control probes are dropped; its DAPI image `morphology_focus/` is shown unless an H&E is aligned), or AnnData file(s) with raw counts (layer `counts`, `X`, or `raw.X` when `X` is normalized, as scanpy leaves it) and cell centroids in micrometers (`obs` columns `x_centroid`/`y_centroid` or `x`/`y`, or `obsm["spatial"]`); several files (e.g. one per lineage) are merged and cell type labels in them are ignored. The dialog lists the h5ad's columns to pick the section when a file holds several. An R file (`.rds`, `.rda` / `.RData`) works in place of an h5ad, read without R: a Seurat object (v3/v4 `Assay` counts, or v5 `Assay5` counts layers, also split per sample; coordinates from an FOV's centroids as `LoadXenium` / `CreateFOV` make them, SlideSeq coordinates, or `meta.data` columns) or a SpatialExperiment / SingleCellExperiment (assay `counts`, `spatialCoords` or `colData` columns); its `meta.data` / `colData` columns act as obs columns. It is converted once into `converted/` of the data folder; normalized and scaled data, graphs, reductions and commands are skipped while reading. |
| image | the H&E (or a fluorescence image) as a pyramidal TIFF / OME-TIFF, or **No image**. |
| registration | an affine CSV with columns `a11,a12,a21,a22,tx,ty` (`px = a11*x + a12*y + tx`, `py = a21*x + a22*y + ty`; cell um to image full-resolution pixels; with one row per section, the row of the chosen section) or a Xenium Explorer alignment matrix (3x3 CSV, image pixels to morphology pixels); or **No registration**. |

The folded part at the bottom of the dialog, **Already imported: change the image or registration**, is for the sample
already open: add or replace its image, give it a registration file, align it with the agent, or undo an alignment.
Until the image is registered the cells are shown on a blank background (or the DAPI) with a yellow bar over the
map. **No registration file?** **Align with Agent** (next to No registration, or the bar) asks the agent (your Claude Code) to align it:
the tissue on the image is matched with the cells' footprint over 4 rotations, mirroring and a range of pixel sizes
(around the size stored in an OME-TIFF), Claude picks the right section and orientation from the candidates' figures,
and a refinement moves the cells onto their nuclei (nuclei on the image vs cell centroids in windows, affine fit,
three rounds). On SCP2927 section 147328 (two sections on the slide, rotated 180 degrees) the search took 8 s and the
refinement 3 s; the result agreed with the nuclei-corrected registration to 0.5 um (median; 1.6 um for 95% of
cells), while the shipped registration file was 13 um off. **Undo last alignment** goes back. A registration needs
at least 90% of the cells on the image. Regions drawn before the cells move are moved along with them.

From the command line (`--he` and `--registration` are optional):

```bash
.venv/bin/stagent import-xenium --name ln1 --xenium /path/to/output-XETG... [--he he.ome.tif --registration he_imagealignment.csv]
.venv/bin/stagent import --name <sample_id> --h5ad cells.h5ad [--he he.ome.tif --registration alignment.csv] \
  [--section-column donor_id --section 147328]
```

### Several slices (projects)

Tick several sections of an h5ad in **Import…** (or give several Xenium folders, one per line): each slice becomes
a sample `<project>_<section>` of one project, with the metadata that is constant within the slice (condition, day,
mouse, …). Images: one file for all slices, one file per slice whose name holds the section (e.g.
`..._147328_HE.ome.tif`; choose them together), or a table (CSV/TSV: section, image path); a registration file has
one row per section. Without an image, the cells are drawn on a blank background wherever their coordinates are
(negative ones too). In the app, **◀ ▶** (or `[` `]`) step through the slices and keep what is shown (Color by, legend filters, hidden H&E, the UMAP
window and, on a UMAP of all slices, its zoom and the cells selected on it), **Slices…** shows them side by side and
**Compare…** opens the **Compare** tab on the composition of every slice.

With **all slices together** ticked, **Annotate** pools the slices (one marker step and one model; each slice gets its
own session), **Run PCA + UMAP** makes one embedding of all slices (with **Harmony** on the PCs by slice, optional) and
**Run niches** finds niches shared by all slices. On four slices of the Cadinu et al. MERFISH mouse colon (healthy,
DSS day 3, 9 and 21; 66,818 cells) the joint annotation took about 15-40 s, the joint embedding 15 s and the niches
9 s; lineage agreement with the authors' labels was 87-91% on the healthy, day-3 and day-9 slices (74% for another healthy slice
annotated alone) and 46% on the day-21 slice, whose cells have few transcripts (median 47 vs 67-93).

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
.venv/bin/stagent import --name 147328 \
  --h5ad $D/scp_xenium_human_{epi,fib,lymph,myeloid,stroma}.h5ad \
  --he $D/human_slideid-0050423-108731-147328.ome.tif \
  --registration examples/scp2927/registration.csv \
  --section-column donor_id --section 147328 --region-column histopathology
```

Marker list for the 480-gene IBD panel: `examples/scp2927/ibd_xenium480_markers.csv` (import it in **Markers…**).

## 3. Work with Claude

Claude runs in your own Claude Code, with your own Claude account; the app connects it and shows every change.
Click **Ask Agent** (the bar on the map) in the app: Terminal opens with Claude Code connected to the app and to the sample you are
looking at (install Claude Code and log in once beforehand).

Claude loads what you are looking at (`use_viewer_session`), and when it starts a new annotation (for example with
a revised marker list) the app switches to it by itself. Examples:

- "Review my annotation." (markers and layout; suggestions appear in the app's Review panel)
- "Draft a marker list for this panel." / "Epithelial cells are only found in the mucosa." (priors are previewed;
  say go ahead to apply)
- "Cluster 7 of the T cells is mixed, split it." / "These cells in region crypt_1 are goblet cells."
- "Show me the biggest conflict group." / "Go back to version 3." / "Export the result."

Claude answers in the language you write in.

### Working on a server

The app can work on a server over SSH: the data and every computation stay on the server, this computer shows the
viewer and runs Claude Code. Click **This computer** (in the toolbar on the right) -> **A server over SSH**, enter the server (a `Host`
of `~/.ssh/config` or `user@server`), the full path of STAnnotator's `stagent` command there and, optionally,
its data folder (default `~/stagent` on the server), then **Test** and **Connect**. The app remembers the
server and connects again when it starts; **This computer** -> **Use this computer** goes back.

On the server, install STAnnotator once (the same version as the app; no window needed):

```bash
cd ~ && unzip -oq STAnnotator-<version>.zip && cd STAnnotator
uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -e .
# the command to enter in the app: ~/STAnnotator/.venv/bin/stagent
```

The SSH login must work without a password prompt (a key: `ssh-keygen`, then `ssh-copy-id <server>`; or ssh-agent).
While connected, **Choose…** buttons browse the server's files, **Export…** and the other Save buttons download into
files on this computer, and **Ask Agent** opens Claude Code here with its tools running on the server (through `ssh`), on the
same data as the viewer. The viewer keeps running on the server between connections (`stagent remote-serve`
starts it, on the server's own localhost only; its process id is in `<data folder>/remote_viewer.json`).

## 4. The viewer

The toolbar is a panel on the right of the map: **Import** and **Export** at the top, the server choice below, then
the steps **Data**, **Annotate**, **Embed**, **Niches**, **Compare** and **Signaling** as tabs (a dot on a tab whose analysis is running), and at
the bottom **Reset** and **Remove all**. On the map, next to the **Color by** bar, a bar holds **Undo**, **Ask Agent**
(Claude Code), **Agent Review** and **Report**. **›** at its top (or `t`) folds it to the right edge for a bigger map,
**‹** brings it back.

The app's window (or the web viewer: `.venv/bin/stagent viewer`, then http://localhost:8765). Changes made
by Claude show up within a few seconds.

1. Pick a sample; earlier annotations are under **Annotation** (Data tab). Choose the **Method** (Annotate tab): *Marker genes* (a marker
   list, **Markers…**), *Single-cell reference* (a labeled reference h5ad, **Reference…**) or *Load existing* (a
   CSV with one cell type per cell, **Annotation file…**), then click
   **Annotate** (about 1 min for 350k cells with markers, seconds with a reference). Regions start empty: you
   draw them.
2. The result is drawn on the H&E: pan and zoom down to single cells; color by cell type, lineage, region,
   changes, conflicts, confidence or clusters, or by **gene expression…** (type any panel gene: cells are colored by
   their transcript count, 0 / 1 / 2–3 / 4–7 / 8+, and a clicked cell shows its count). Click a legend entry to show only that type; click a cell for
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
5. **Export…** saves the **cells table** (CSV, one row per cell: `cell_id`, `region`, `cell_type`, `lineage`,
   `level_1`… the full type path, `confidence`, `expression_label`, `marker_cluster`, `hand_labeled`, and when
   computed `niche`, `leiden`, `umap_1`, `umap_2`; then `x_um`, `y_um`), the **registration** (CSV: the affine
   a11..ty that places the cells on the image now, whether it came from the file, Claude's alignment or a
   correction; it can be imported again), or the **analysis notebook** (Jupyter). The notebook runs without this
   app: it needs the data files and `numpy pandas scipy scikit-learn python-igraph pynndescent umap-learn h5py
   matplotlib`. Standard steps are plain code (reading the cells, the registration affine, PCA / neighbors /
   Leiden / UMAP, NMF niches, figures, the table); the app's annotation algorithm (marker clustering, expression
   model, priors) is carried as source code in a folded cell and replayed with this annotation's marker list,
   regions, priors, corrections and re-clusterings. The results at export time are filled in, and a rerun checks
   that labels, Leiden clusters and niches come out the same: identical when run on the same kind of machine (with
   a GPU and PyTorch, neighbors are exact as in the app; without, large samples use approximate neighbors);
   annotations made with an earlier app version can differ. UMAP coordinates differ slightly (all cores).
6. **Reset** removes every annotation and chat of the selected sample (moved to `<data folder>/trash/`); in a project,
   of this slice or of all slices, and optionally the PCA / UMAP and Leiden clusters too (they are computed from the
   expression alone and stay otherwise). **Compare…** lists a slice without results as "no annotation".
   **Remove all…** (next to Import, desktop app) deletes every sample with all its annotations and chats for good,
   to import from scratch (it cannot be undone); marker lists stay.
7. **Markers…** edits, checks, saves and loads marker lists (saved lists go to `<data folder>/marker_lists/`).
8. **Agent Review** (the bar on the map, with **Undo**, **Ask Agent** and **Report**) asks which part of the analysis
   the agent should review: **Annotation**, **Spatial organization**, **Embedding**, **Niches**, **Compare** or
   **Signaling** (the last four once they have been run; Compare and Signaling review the analyses you ran in
   those tabs, which the app logs for the agent). The agent looks only at that part, drafts its advice
   (`draft_advice`), then checks every draft as a skeptical reviewer: it reruns the tool that would confirm or
   refute it, looks at the figure on the H&E and weighs other explanations (spillover, segmentation, batch, too few
   cells, panel limits), and records `check_advice`. Only the drafts that hold up reach the Review panel, each with
   what was checked ("✓ checked: ..."); the dropped ones are folded under "Dropped by the agent after checking" with
   the reason. The panel's menu switches between the parts (with the number of open items of each); **Show panel**
   in the dialog opens a part without asking the agent.
   - *Annotation*: computed marker checks (markers that are weak or unspecific in their cells, negative markers that
     are detected, genes of Unassigned cells) and areas to check (cells of sheet-forming types, such as epithelium
     or smooth muscle, apart from their own lineage), outlined on the map. The agent checks every candidate on the
     H&E and with marker biology (`record_review`); afterwards only the real problems are shown, each with its
     reason (false alarms are folded away). Verdicts carry over to later versions, so only new candidates are
     checked again; **solve** removes an item you have dealt with. When the labels changed since its last look it
     also gives advice on *Cell types* (types to merge, split or add: **open in Markers…** for a revised list, **ask
     the agent** for a re-clustering), *Sample context* and *Technical* problems: low-quality tiles, marker
     spillover next to the epithelium, and the registration. The registration check correlates the nuclei on the
     H&E (hematoxylin) or DAPI with the cell centroids in windows across the section and fits an affine;
     **correct registration…** moves every cell onto its nucleus (regions keep their cells; **undo correction**
     puts them back). On SCP2927 section 147328 it found a 0.28% scale mismatch (cells up to 30 um off their nuclei
     at the corners); after the correction the windows agree within 0.5 um.
   - *Spatial organization*: compartments (outlines of the mucosa, muscularis or lymphoid follicles, computed from
     where a lineage forms tissue and checked on the H&E; **add region**, then **add prior…** for the prior that
     uses it), and gradients and neighbors, found or contradicted (**add prior…**).
   - *Embedding*: Leiden clusters against the cell types (mixed clusters, types split over clusters) and, for a
     joint embedding, clusters made of one slice (`embedding_facts`, `show_embedding`); **put in Embed** fills in
     suggested settings.
   - *Niches*: what each niche is on the H&E, names (**name the niches**), k and the number of niches (**put in
     Niches**).
   - *Compare* and *Signaling*: what the DE genes, compositions or interactions you looked at mean, the caveats
     (replicate slices, spillover, panel coverage) and the analysis that would settle a question (**put in
     Compare / Signaling**); **show on map** colors the map by a gene or a signal.
   Advice changes nothing until you click; **dismiss** hides it.
   The **Priors** list under Regions shows the priors of the annotation (the user's words and how they were read);
   **×** removes one (Undo brings it back). **?** (top of the toolbar on the right) explains every control and, under **About**, shows the
   build, the data folder, whether Claude Code was found and when Claude last used the tools, with a button that
   copies these details for a bug report.
9. **Embedding** (its own block of the toolbar): **Run PCA + UMAP** computes a PCA, a UMAP and Leiden clusters of the sample from
   expression only (shared by all its annotations; about 40–50 s for 356,629 cells: the kNN graph comes from the
   fast GPU / pynndescent search and is handed to UMAP, which starts from PC1/PC2 and runs on all cores).
   **Show UMAP** opens a window on the right (UMAP or PC1/PC2, colored as the map, by cell type, Leiden cluster or
   niche). **Select** (or Shift-drag) draws around cells: they light up on the spatial map and the rest turn gray;
   **Show on map** zooms there. Leiden clusters are also a *Color by* option.
10. **Niches** (fourth row): places where cells of similar make-up live together. Every cell gets an identity, its
   weights on NMF gene programs (*Identity: NMF programs*, independent of the annotation) or its cell type; the
   identities are averaged over the cell and its *k* nearest cells and clustered (k-means) into the chosen number of
   niches. **Run niches** uses one k everywhere (**k: Set k**); with **k: Agent chooses** the same button asks the agent to choose k per drawn region
   from cell density and cell types (or a radius, so neighborhoods cover a similar stretch of tissue in dense and
   sparse areas), run it and name the niches. The result appears as *Color by: niche* (names in the legend) and as
   a `niche` column in Export CSV. On SCP2927 section 147328 (356,629 cells) a run takes about 6 s (1.5 s once the
   NMF is cached) and the 8 niches match the tissue layers (surface epithelium, crypt base, lamina propria, lymphoid
   follicles, muscularis mucosae, muscularis propria, submucosal vessels and stroma).
11. **Compare** (fifth tab), for this slice or all slices of a project pooled:
   - **DE genes**: the genes that differ between one group and the other cells (or another group). A group is a cell
     type, lineage, region, niche, Leiden cluster, the cells selected in the UMAP window, a slice, or a value of a
     slice metadata column such as condition; **Only within** restricts both sides (e.g. fibroblasts in the mucosa vs
     fibroblasts elsewhere). Wilcoxon rank-sum test on log-normalized expression (up to 20,000 random cells per side;
     means and detection on all cells), genes ranked by AUC with log2 fold change and FDR. With pooled slices the
     change is checked in every slice that has both groups, or, when the groups are different slices, every slice
     of the group against every slice of the reference; with one slice per side the app says that a slice effect
     cannot be told from a group effect. Click a gene to color the map by it.
   - **Top genes**: every group against the rest (what sets each niche or cluster apart).
   - **Composition**: how one grouping splits over another (cell types per niche, niches per region, niches per
     slice or per condition, with every slice's own bar), as shares or as log2 enrichment over independence.
   **Ask Agent** asks the agent to interpret the result (it calls the same analysis through `de_genes`,
   `group_markers` or `composition`); **Save CSV…** keeps the full table. On SCP2927 (559,539 cells, 480 genes) a DE
   run takes 1-2 s and the top genes of 8 niches about 10 s.
12. **Signaling** (sixth tab): cell-cell communication, for this slice or all slices of a project pooled.
   - Ligand-receptor pairs: CellChatDB v2 (pathways and the kind of signaling), CellPhoneDB v5 (human) and the
     OmniPath consensus that LIANA uses (CellPhoneDB, CellChatDB, ICELLNET, connectomeDB2020, CellTalkDB), merged into
     `stagent/data/lr_human.csv.gz` (6,521 pairs) and `lr_mouse.csv.gz` (5,792) by `python -m stagent.lr_resource`;
     or one of them, or your own CSV (agent tools). A pair counts when every gene of the ligand and receptor
     (complexes have several subunits; a complex's level is its lowest subunit) is in the panel: 109 pairs for the
     480-gene Xenium IBD panel, about 1,200 for the 990-gene MERFISH mouse colon panel. Small-molecule signals (an
     enzyme stands in for the molecule) are off unless switched on.
   - Scoring at single-cell resolution: neighbor pairs are actual cells, never across slices. Contact signals:
     cells within 15 µm; secreted and ECM signals: the 30 nearest cells within 50 µm, weighted by a Gaussian of the
     distance. The score of sender group A → receiver group B is the mean of ligand (sender) × receptor (receiver)
     over their neighbor pairs (log-normalized expression); it is tested against random group labels (the
     CellPhoneDB idea) with a z-score that counts distinct cells rather than pairs (on SCP2927 it is more
     conservative than 100 label permutations), FDR over all tests; tested only when the ligand is detected in 10% of
     A and the receptor in 10% of B. Colocalization above 1 means the expressing cells sit together more than the
     group averages imply.
   - Results: **Interactions** (score, specificity, colocalization, detection, neighbor pairs, FDR), a sender ×
     receiver **Network**, **Pathways**; **receivers** / **senders** colors the map by the signal each cell receives
     (its receptor × the ligand around it) or sends. On SCP2927 section 147328, CCL19 → CCR7 is received in the
     lymphoid follicles. **Compare groups** runs it per condition, slice, region or niche (largest differences
     first; with one slice per group the app says a slice effect cannot be told apart). Agent tools:
     `cell_communication`, `compare_communication`, `show_communication` (the signal on the H&E). A run on 356,629
     cells takes about 10 s.
13. **Report** shows the annotation report the agent writes (**Ask Agent to write it**): sample, method and inputs,
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

**Single-cell reference** (like CellTypist). Give a labeled single-cell h5ad (or a Seurat / SingleCellExperiment `.rds`) of the same tissue, its cell type
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

- The app does not open: see `<data folder>/app.log`; run `.venv/bin/stagent app` in a terminal to see
  errors, or `--browser` to use the browser. The web viewer logs to `<data folder>/viewer.log`.
- Claude Code does not open: it must be installed and logged in (`claude` works in Terminal). By hand:
  `cd ~/stagent/claude_workspace && claude --mcp-config stagent_mcp.json`.
- Import says too few cells land on the H&E: check the registration file, `--section`, and that the cell
  coordinates are in micrometers.

## Other clients

- Claude Code without the app: `claude mcp add stagent -e STAGENT_HOME=$HOME/stagent -- $PWD/.venv/bin/python -m stagent.mcp_server`
  (run in this folder). `.claude/skills/spatial-annotation/` holds the agent's working rules for Claude Code.
- Shared server: the viewer binds to localhost only. Others reach it through an SSH tunnel
  (`ssh -N -L 8765:localhost:8765 <user>@<server>`), then open http://localhost:8765. The conflict review page for
  pathologists is at http://localhost:8765/review (verdict keys 1-6, rules form, CSV download).

## Building the app (maintainers)

On a Mac with Apple silicon, from this folder:

```bash
packaging/macos/build_app.sh            # -> dist/STAnnotator.app and dist/STAnnotator.dmg (ad-hoc signed)
```

It builds with PyInstaller (`packaging/stannotator.spec`) in its own environment, without PyTorch. The app's
executable also takes the command-line subcommands (`STAnnotator mcp` is what Claude Code starts). To give it a Developer ID signature and Apple's notarization, so that it opens without warnings on other
Macs (needs an Apple Developer account):

```bash
xcrun notarytool store-credentials spatial-annotator --apple-id YOU@EXAMPLE.COM --team-id TEAMID   # once
DEVELOPER_ID="Developer ID Application: Your Name (TEAMID)" packaging/macos/sign_and_notarize.sh
```

The hardened-runtime exceptions it uses (`packaging/macos/entitlements.plist`) allow numba's run-time compilation.

## Layout

```
stagent/     Python package
  cli.py           stagent app / make-app / import / import-xenium / samples / viewer / mcp
  desktop.py       the desktop app (pywebview window around the local server)
  claudelink.py    connection to the user's Claude Code; viewer <-> agent handoff files
  importer.py      Xenium output folder, or h5ad / RDS + H&E + registration -> sample bundle
  rds.py           R data files (.rds / .rda) read without R: Seurat, SpatialExperiment -> cells, genes, counts, xy
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
tests/             smoke tests (python tests/smoke_engine.py, tests/smoke_mcp.py expect sample 147328 / 167799; tests/smoke_joint.py runs in a temporary folder)
```
