# STAnnotator user guide

The full description of the app: import, the viewer, the agent, the analyses, and how each method works. Installation
and the local / server setups are in the [README](../README.md).

**Contents**

1. [Import a sample](#1-import-a-sample): Xenium and AnnData / R files, Visium HD, Stereo-seq, CosMx and MERFISH,
   H&E only, aligning the image, regions drawn elsewhere, several slices, an example
2. [Marker lists](#2-marker-lists)
3. [Work with the agent](#3-work-with-the-agent): asking, reviewing Claude's changes, working on a server
4. [The viewer](#4-the-viewer): the window, samples, annotating, the map, regions, corrections and Undo, Export,
   Agent Review, Help and About, storage
5. [Analyses](#5-analyses): Embedding, Niche, Compare, Signal, Trajectory, Report
6. [Annotation methods](#6-annotation-methods)
7. [Performance](#7-performance)
8. [Troubleshooting](#8-troubleshooting): bug reports, usage sharing, common problems
9. [Other clients](#other-clients), [On Windows](#on-windows), [Layout](#layout)

## 1. Import a sample

Click **Import** in the app (or use the [commands](#from-the-command-line) below). The dialog asks for three parts:

1. **the cells** (the platform's output, or an h5ad / R file);
2. **No H&E** or an **H&E file**;
3. a **Registration file** or **No registration**.

The file fields of Import grow with what they hold, so long paths and several lines need no scrolling. To give a
sample another image or registration file, import it again with them under the same name.

### Xenium

- Only `cells.parquet` (or `cells.csv.gz`) and `cell_feature_matrix.h5` (or `matrix.mtx.gz`, `barcodes.tsv.gz`,
  `features.tsv.gz`) are needed; optionally `experiment.xenium` (the pixel size) and
  `morphology_focus/ch0000_dapi.ome.tif` (shown unless an H&E is aligned).
- Give each file in its field of **Import** (any names). For several slices, one file per line in each field, matched
  by a name prefix such as `GSM1234_`. The command line also takes the whole output folder.
- Control probes are dropped.

### AnnData (h5ad)

- Raw counts: layer `counts`, `X`, or `raw.X` when `X` is normalized (as scanpy leaves it).
- Cell centroids in micrometers: `obs` columns `x_centroid` / `y_centroid` or `x` / `y`, or `obsm["spatial"]`.
- Several files (e.g. one per lineage) are merged; cell type labels in them are ignored.
- When a file holds several sections, the dialog lists the h5ad's columns to pick the section.

### R files (.rds, .rda / .RData)

An R file works in place of an h5ad and is read without R:

- a **Seurat** object: v3/v4 `Assay` counts, or v5 `Assay5` counts layers (also split per sample); coordinates from an
  FOV's centroids as `LoadXenium` / `CreateFOV` make them, SlideSeq coordinates, or `meta.data` columns;
- a **SpatialExperiment / SingleCellExperiment**: assay `counts`, `spatialCoords` or `colData` columns.

Its `meta.data` / `colData` columns act as obs columns. It is converted once into `converted/` of the data folder;
normalized and scaled data, graphs, reductions and commands are skipped while reading.

### H&E image and registration file

- **H&E image**: a pyramidal TIFF / OME-TIFF, or **No H&E**. Without an H&E, Xenium shows the DAPI image of its own
  field; for h5ad / RDS a DAPI or fluorescence image can stand in for the H&E (**Image kind**).
- **Registration file**, one of:
  - an affine CSV with columns `a11,a12,a21,a22,tx,ty` (`px = a11*x + a12*y + tx`, `py = a21*x + a22*y + ty`; cell um
    to the image's full-resolution pixels; with one row per section, the row of the chosen section);
  - a Xenium Explorer alignment matrix (3x3 CSV, image pixels to morphology pixels);
  - or **No registration**.

### No registration file: Align with Agent

Until the image is registered, the cells are shown on a blank background (or the DAPI) with a yellow bar over the
map. **Align with Agent** (next to No registration, or on the bar) asks the agent (your Claude Code) to align it:

1. the tissue on the image is matched with the cells' footprint over 4 rotations, mirroring and a range of pixel
   sizes (around the size stored in an OME-TIFF);
2. Claude picks the right section and orientation from the candidates' figures;
3. a refinement moves the cells onto their nuclei (nuclei on the image vs cell centroids in windows, affine fit,
   three rounds).

A registration needs at least 90% of the cells on the image. Regions drawn before the cells move are moved along
with them.

On SCP2927 section 147328 (two sections on the slide, rotated 180 degrees) the search took 8 s and the refinement 3
s; the result agreed with the nuclei-corrected registration to 0.5 um (median; 1.6 um for 95% of cells), while the
shipped registration file was 13 um off.

### Visium HD and Visium

Choose Space Ranger's `outs/`, `binned_outputs/` or one `square_XXXum` folder.

**Cells from outlines** (the default when there are outlines):

- Space Ranger 4's `segmented_outputs/cell_segmentations.geojson` is found by itself, or **Cell outlines** takes a
  GeoJSON of outlines in the full-resolution image's pixels (QuPath, StarDist, Cellpose).
- Every 2 um bin is split between the cells it overlaps by area (16 points per bin against the outlines drawn at 0.5
  um): each of its counts goes to a cell with the probability of that cell's share, so counts stay whole and their
  expectation is area-proportional. The part of a bin outside every cell is left out.
- On the mouse small intestine example (Space Ranger 4's 282,346 outlines over its 5,479,660 bins of 2 um) this gave
  281,984 cells (median 52 um2, 335 counts) in 68 s, 92% of all counts in a cell. The cells' totals agree with Space
  Ranger's own whole-bin assignment (sum 0.999, per-cell correlation 0.95; the outlines of its 4.0.1 run sit about 0.6
  um off the 3.0 run's bins, which the area split absorbs). The marker annotation then finds the cells that 8 um bins
  blur: macrophages (16,163 against almost none), enteroendocrine cells (461 against 62), monocytes, neutrophils,
  dendritic cells.

**Bins** (without outlines, or with **Unit: 8 um bins**):

- 2 um bins are pooled into 8 um bins (counts summed by grid position; the pooled bins match Space Ranger's own 8 um
  bins, the same counts for 95% of them). 8 um bins and larger ones (16 um, Visium's 55 um spots) are used as they
  are; given the whole output, its 8 um bins are used.

**Image**: Space Ranger's alignment places everything on its image, so no registration is needed. **H&E file** is the
full-resolution image Space Ranger was given (a copy with zoom levels is made in the sample's folder when it has
none); without it the output's `tissue_hires_image.png` (about 1 um per pixel) is shown.

**Whole transcriptome** (some 19,000 genes): the PCA uses the 3,000 most variable genes, the expression model those
and the marker list's genes (a panel keeps all its genes, Xenium 5K included). The marker annotation of the 282k cells
took 118 s (5.3 GB of memory). The nuclei-based registration check of Agent Review does not apply.

Command line: `stagent import-visium --name si --path outs/ [--segmentation cells.geojson] [--unit bins] [--he image.btf]`.

### Stereo-seq

Give the STOmics SAW output: its folder, or its files. Its spots (DNBs) are 0.5 um apart.

- **Cells**: when SAW segmented the cells (`<SN>.cellbin.gef`, the adjusted one when there is one, or a GEM with a
  cell column, `CellID` or MOSTA's `label`). A DNB is a point, so it belongs to one cell; DNBs outside every cell are
  left out.
- **Bins**: otherwise, or with **Unit: 8 um bins**, the DNBs of the GEM (`<SN>.tissue.gem.gz`, bin1) are pooled into
  8 um bins (16 DNBs a side).
- Genes are named by `geneName` when the file has it (Ensembl ids sharing a name are summed).
- **Image**: the registered image SAW writes (`<SN>_HE_regist.tif`, or the ssDNA one) is in the DNBs' coordinates,
  so no registration is needed; a copy with zoom levels is made in the sample's folder. Another image found (not a
  `_regist.tif`, such as MOSTA's raw ssDNA scan) is not assumed to be aligned: it waits for **Align with Agent** (on
  MOSTA the search finds it mirrored and rotated 270 degrees, x and y swapped, and the refinement ends at 0 um).
- The bin GEF (`tissue.gef`) is not read: use the GEM SAW writes beside it.

Tested on the STOmics demo mouse brain C04042E3 (cell-bin GEF: 130,588 cells, the same median counts as SAW reports,
35 s; its GEM of 247,751,988 rows as 816,465 bins of 8 um, 229 s, 6.7 GB; the cells sit on the H&E's nuclei) and the
MOSTA adult mouse brain (cell-bin GEM: 56,731 cells, 62 s).

Command line: `stagent import-stereoseq --name brain --path saw_outs/ [--unit bins] [--image regist.tif]`.

### CosMx and MERFISH (Vizgen MERSCOPE)

Choose the kind in **Import** and give its two files, each in its field (Counts, Cell metadata); one slide or region
per import (import several one after the other).

- **CosMx**: a slide's AtoMx flat files `<slide>_exprMat_file.csv` and `<slide>_metadata_file.csv`. A cell is
  `fov_cellID` (cell_ID 0, transcripts outside cells, is left out); its position comes from `CenterX/Y_global_px` at
  0.12028 um per pixel (or `x/y_slide_mm`).
- **MERSCOPE**: a region's `cell_by_gene.csv` and `cell_metadata.csv` (`center_x`, `center_y` in um).
- Control probes (NegPrb, Negative, SystemControl, FalseCode, Blank) are left out.
- The instruments' own images (CosMx FOV images, MERSCOPE mosaics) are not read yet: the cells are drawn on a blank
  background, or on your H&E with an affine registration CSV (cell um to H&E pixels), or aligned by Claude (**Align
  with Agent**).
- Claude can import them too (`import_cosmx_merfish`, through the review window).

### H&E only

An image without any cells, e.g. for a pathologist to draw and name regions before the expression data exist.

- Formats: TIFF / OME-TIFF, SVS, NDPI, PNG or JPEG. Other whole-slide formats such as MRXS or CZI must be converted to
  OME-TIFF first. An image without zoom levels is copied into the sample folder as a pyramid.
- The pixel size is read from the image (OME-TIFF, the SVS `MPP`, or the TIFF resolution tags) or typed in; without
  one the app shows no micrometers and gives areas in pixels.
- Such a sample has only the Data tab and the regions (the panel shows each outline's area, mm² with a pixel size).
  **Export** writes the outlines as GeoJSON.
- **Its regions go to the data** in two ways:
  - choose the sample under **Regions from** in **Import** when importing the data registered to the same H&E (the
    same scan: the image size is checked); every annotation of the data starts with them;
  - or give the exported GeoJSON as the import's **Regions** (QuPath opens it too: every outline is an annotation
    classified by its region, with its color).
- **The agent can outline compartments** on such an image when you ask (e.g. "draw the mucosa"):
  1. it reads the organ's priors (`tissue_priors`: literature-based priors of 13 organs, with how each compartment
     looks on the H&E);
  2. it splits the tissue into classes by its stains (`segment_he`: hematoxylin and eosin by color deconvolution and
     their texture, k-means; each class with its distance to the tissue edge);
  3. it proposes the classes that match as a region (`propose_he_region`, or points read on `show_he`). You see the
     outline in the review window, and **Accept** adds it as a draft to edit.

Agent tools: `import_he_only`, `export_regions`, `tissue_priors`, `segment_he`, `show_he`, `propose_he_region`.

### Regions drawn elsewhere (QuPath)

E.g. QuPath annotations on the H&E, exported as GeoJSON in the image's full-resolution pixels.

- **At import**: give the file under **Regions** in **Import** (several slices: one per slice, in their order), and
  every new annotation of the sample starts with them.
- **Later**: click **Import** in the map's Regions panel to add them to the open annotation.
- A region's name is the object's classification (else its name), its color QuPath's. Several outlines of one class
  form one region; holes stay out of it (drawn dashed); points, lines and cell detections are skipped.
- Outlines traced along the image's pixels (QuPath's brush and wand) are a staircase when zoomed in; they are made
  smooth curves at import (99% within about 2 pixels of the original, well under 1 um; a few cells in ten thousand
  change side). Outlines with free points are kept point for point.
- The H&E must be the image the outlines were drawn on. Drawing before the registration is fine (the outlines are in
  the image's own pixels): they apply once the cells sit on it. An H&E still waiting for its alignment keeps the file,
  and new annotations start with the regions after the alignment (for an annotation made before: Regions panel,
  **Import**).

### From the command line

`--he` and `--registration` are optional:

```bash
.venv/bin/stagent import-xenium --name ln1 --xenium /path/to/output-XETG... [--he he.ome.tif --registration he_imagealignment.csv]
.venv/bin/stagent import --name <sample_id> --h5ad cells.h5ad [--he he.ome.tif --registration alignment.csv] \
  [--section-column donor_id --section 147328]
```

### Several slices (projects)

**Importing.** Tick several sections of an h5ad in **Import** (or give the Xenium files of several slices, one per
line in each field). Each slice becomes a sample `<project>_<section>` of one project, with the metadata that is
constant within the slice (condition, day, mouse, …).

- Images: one file for all slices, one image per slice in the order of the slices (one per line), or a table
  (CSV/TSV: section, image path).
- Registration: one file for all slices (a row per slice, matched by its name) or one file per slice in their order.
- The dialog lists which image and registration file each slice gets.
- Without an image, the cells are drawn on a blank background wherever their coordinates are (negative ones too).

**Organizing.** **Projects** (Data tab) puts any samples into any project, moves a slice to another project or takes
it out ("no project"). The slices of a project need one gene panel; a project left without slices goes away.

**Browsing.** **◀ ▶** (or `[` `]`) step through the slices and keep what is shown (Color by, legend filters, hidden
H&E, the UMAP window and, on a UMAP of all slices, its zoom and the cells selected on it). **Overview of all samples**
shows them side by side, and the **Compare** tab compares them (all slices, by slice).

**All slices together.** With this ticked:

- **Annotate** pools the slices (one marker step and one model; each slice gets its own session);
- **Run PCA + UMAP** makes one embedding of all slices (optionally with **Harmony** on the PCs by slice);
- **Run niches** finds niches shared by all slices.

On four slices of the Cadinu et al. MERFISH mouse colon (healthy, DSS day 3, 9 and 21; 66,818 cells) the joint
annotation took about 15-40 s, the joint embedding 15 s and the niches 9 s. Lineage agreement with the authors' labels
was 87-91% on the healthy, day-3 and day-9 slices (74% for another healthy slice annotated alone) and 46% on the day-21
slice, whose cells have few transcripts (median 47 vs 67-93).

### Example: the IBD Xenium atlas (SCP2927)

Data of Pokatayev et al., Nature 2026.

- Download the five `scp_xenium_human_*.h5ad` files from Single Cell Portal study
  [SCP2927](https://singlecell.broadinstitute.org/single_cell/study/SCP2927).
- Download the H&E image of the slide from [Zenodo record 17518435](https://zenodo.org/records/17518435)
  (`examples/scp2927/he_files.tsv` lists the eight images with sizes, MD5 checksums and download links; 1.4-2.1 GB
  each). Each slide holds two sections (`human_slideid-<slide>-<donor>-<donor>`).
- The registration file is shipped here: `examples/scp2927/registration.csv` (16 sections; columns `slide`, `donor`).

```bash
D=~/Downloads/SCP2927
.venv/bin/stagent import --name 147328 \
  --h5ad $D/scp_xenium_human_{epi,fib,lymph,myeloid,stroma}.h5ad \
  --he $D/human_slideid-0050423-108731-147328.ome.tif \
  --registration examples/scp2927/registration.csv \
  --section-column donor_id --section 147328 --region-column histopathology
```

Marker list for the 480-gene IBD panel: `examples/scp2927/ibd_xenium480_markers.csv` (import it in **Markers…**).

## 2. Marker lists

A sample annotated with marker genes needs a marker list. **Markers…** edits, checks, saves and loads marker lists
(saved lists go to `<data folder>/marker_lists/`).

- **Import…** a CSV (or text) file.
- **Start from an organ**: choose the organ (brain, breast, colon, heart, kidney, liver, lung, lymph node / tonsil,
  pancreas, skin, small intestine, stomach). **Use built-in list** keeps the cell types your panel can resolve (two
  of their markers in the panel, or one decisive marker such as INS) with their first markers, in the editor's
  indented format (`Name: GENE+ GENE-`), and lists what it left out. The built-in lists are in
  `stagent/data/organs/` (sources in its README).
- **Refine with the agent** asks the agent to complete it for the tissue and panel (searching the literature for
  types the built-in list lacks) and save it under "load a list…". For an organ that is not built in (**Other
  organ…**) the agent starts from the immune and stromal cells and adds the organ's own.
- You can also ask Claude Code "draft a marker list for this panel" and load it under "load a list…".
- Imported and saved lists stay in that list for the next samples; **×** after a name removes a saved list (into the
  trash of the data folder).
- **Export CSV** saves the list being edited.

The CSV has one row per cell type:

| cell_type | parent | positive_markers | negative_markers |
|---|---|---|---|
| Epithelial | | EPCAM KRT8 | |
| Absorptive | Epithelial | CA1\|SLC26A3\|FABP2\|AQP8 | SPINK4 |
| Goblet cells | Epithelial | SPINK4 RETNLB ITLN1 | |

`parent` is empty for lineages; markers are separated by spaces (or semicolons); `A|B` means any one of A or B. The
text format and how the list is used are in [Annotation methods](#6-annotation-methods).

## 3. Work with the agent

### Asking

Claude runs in your own Claude Code, with your own Claude account; the app connects it and shows every change. Click
**Ask Agent** (the bar on the map): Terminal opens with Claude Code connected to the app and to the sample you are
looking at. Install Claude Code and log in once beforehand.

Claude loads what you are looking at (`use_viewer_session`), and when it starts a new annotation (for example with a
revised marker list) the app switches to it by itself. Changes made by Claude show up within a few seconds. Examples:

- "Review my annotation." (markers and layout; suggestions appear in the app's Review panel)
- "Draft a marker list for this panel." / "Epithelial cells are only found in the mucosa." (priors are previewed; say
  go ahead to apply)
- "Cluster 7 of the T cells is mixed, split it." / "These cells in region crypt_1 are goblet cells."
- "Show me the biggest conflict group." / "Go back to version 3." / "Export the result."

Claude answers in the language you write in. **?** → **Agent** in the app lists what to ask and what comes back.

The app has no chat box of its own: you talk to Claude in Claude Code. The web viewer on a shared server shows a chat
box that uses the server's own logged-in Claude Code CLI (for the server owner's use).

### Reviewing Claude's changes

Claude never changes the data on its own. Every change it makes from the conversation (a prior, labeling cells or a
cluster, re-clustering, going back to a version, aligning the cells to the image, niches and their names, verdicts, a
new annotation or sample) first opens a **Review Claude's change** window in the app with what it does: the cells it
touches, their current labels, the predicted label changes where they can be computed, your words.

- Nothing changes until you click **Accept**.
- **Dismiss** drops it and tells Claude.
- Several waiting changes come one after the other.
- Claude learns the outcome by itself when you decide within about 45 seconds, otherwise when it next checks.

### Working on a server

The app can work on a server over SSH: the data and every computation stay on the server, this computer shows the
viewer and runs Claude Code.

**Connecting**

1. Click **Local** (in the toolbar on the right, beside **Export**) → **A server over SSH**.
2. Enter the server (a `Host` of `~/.ssh/config` or `user@server`), the full path of STAnnotator's `stagent` command
   there and, optionally, its data folder (default `~/stagent` on the server).
3. **Test**, then **Connect**.

The app remembers the server and connects again when it starts; **Remote** → **Use this computer** goes back.

**On the server**, install STAnnotator once (the same version as the app; no window needed):

```bash
git clone --branch v<version> https://github.com/caokai1073/STAnnotator.git ~/STAnnotator && cd ~/STAnnotator
uv venv .venv --python 3.12 && uv pip install --python .venv/bin/python -e .
# the command to enter in the app: ~/STAnnotator/.venv/bin/stagent
```

The SSH login must work without a password prompt (a key: `ssh-keygen`, then `ssh-copy-id <server>`; or ssh-agent).

**While connected**

- **Choose…** buttons browse the server's files.
- **Export** and the other Save buttons download into files on this computer.
- **Ask Agent** opens Claude Code here with its tools running on the server (through `ssh`), on the same data as the
  viewer.
- The viewer keeps running on the server between connections (`stagent remote-serve` starts it, on the server's own
  localhost only; its process id is in `<data folder>/remote_viewer.json`).

## 4. The viewer

The app's window, or the web viewer (`.venv/bin/stagent viewer`, then http://localhost:8765).

### The window

- **The toolbar** is a panel on the right of the map:
  - at the top: **?** (help), the bug icon (report a bug) and **›** (folds the panel to the right edge for a bigger
    map; also `t`; **‹** brings it back);
  - **Import**, **Export** and **Local** / **Remote** (where the data is: this computer or a server);
  - two tabs, **Data** and **Analysis**. Analysis has the steps **Annotate**, **Embed**, **Niche**, **Compare**,
    **Signal** and **Trajectory** as a second row and reopens the one last used (a dot on a tab whose analysis is
    running);
  - at the bottom, **Reset**, **Remove sample** and **Remove all**.
- **On the map**: the **Color by** bar (with opacity, dot size, hide cells, hide H&E), and a bar with **Undo**,
  **Refresh**, **Ask Agent** (Claude Code), **Agent Review** and **Report**.
- **Panels on the map**: Regions (with the Priors list) at the top left, the cell type legend at the bottom left, the
  overview of the whole sample at the bottom right.

### Samples (Data tab)

- **Pick a sample**; earlier annotations are under **Annotation**.
- **Rename** changes the name the app shows (a slice: its name in the project, and the project's name); the id, files
  and annotations stay.
- **Overview of all samples** shows every sample side by side, each over its image (the **image** box) and colored the
  same way: cell type, lineage, Leiden clusters, niche or region (from the annotation open, else each sample's latest
  one), with one legend for all.
  - Tick legend entries to show only them (e.g. only ILCs; the other cells turn light gray on every sample); **Show
    all** clears the ticks.
  - Click a sample to open it.
  - **Save figure** writes the overview as shown (every sample, its name, the legend) as a PDF.
- **Unload H&E** frees the memory of the sample's image: the cells stay where they are on a blank background of the
  image's size, and the image path and registration are kept, so **Load H&E** puts it back as it was (no new
  registration or alignment). While unloaded, Claude's figures show the cells without the image.
- **Memory**: the image itself takes little memory (tens of MB of tiles); an open annotation takes much more (its
  counts, model and results: a few hundred MB for 350k cells). The app keeps only the last 3 annotations and 4
  samples' cell positions in memory and reads older ones again from disk when they are opened
  (`STAGENT_OPEN_ANNOTATIONS` sets the number of annotations).
- **Remove sample** (bottom of the toolbar) removes only the open sample: its cells, annotations with their chats,
  regions, niches and PCA / UMAP go to the data folder's `trash/`. A project keeps its other slices, and the image and
  source files stay where they are.

### Saving and loading an annotation (Data tab)

- **Save…** under Annotation writes the open annotation into one file (`.stannot`): every version with its regions,
  priors and corrections, the marker step and model, niches and their names, review verdicts and advice, the PCA /
  UMAP / Leiden clusters, and the last result of each Compare and Signaling analysis (these are kept with the
  annotation, so **Last result** also works after a restart). Chats are left out.
- **Load…** opens such a file as a new annotation of the sample with the same cells and genes (here or in another data
  folder); it is marked "loaded copy" in the list. A different PCA / UMAP of that sample goes to the trash first.

### Filtering cells and genes (Data tab)

**Filter cells and genes** deletes:

- cells with fewer or more transcripts than you set, with fewer genes, or with a mitochondrial share (MT- genes, when
  the panel has them) above a percentage;
- genes detected in fewer cells than you set.

Histograms show where each threshold cuts and the dialog counts what goes; empty fields are not applied.

- The data as imported stays in the sample folder (`original/`): new thresholds always start from it, and **Restore
  all** brings every cell and gene back.
- In a project one filter covers all slices (they keep one gene panel).
- Annotations and the PCA / UMAP of the filtered samples go to the trash (they number the cells).
- Marker genes the filter removed are left out of marker lists with a warning.

### Annotating (Analysis → Annotate)

1. Choose the **Method**:
   - *Marker genes*: a marker list (**Markers…**, see [Marker lists](#2-marker-lists));
   - *Single-cell reference*: a labeled reference h5ad (**Reference…**);
   - *Load existing*: a CSV with one cell type per cell (**Annotation file…**).
2. With marker genes, choose **Neighbors**, the neighbor search of the marker step:
   - *Approximate* (pynndescent, the default) is fast on any computer and gives the same result only with the same
     pynndescent and numba versions;
   - *Exact* gives the same neighbors, hence the same clusters and labels, on every computer, GPU or not, but is slow
     without a GPU for large samples.
3. Click **Annotate**.

**Time estimate.** Under Neighbors the app estimates how long the whole annotation will take here, for both choices
(hover for the steps), and the progress bar shows the time so far against that estimate. The estimate starts from a
short speed test of the computer and learns from every annotation run on it. For 357k cells: under a minute with a GPU
either way; without one, about 5 min exact and about 1 min approximate on a 12-core CPU with PyTorch (numpy alone
makes the exact search about 4 times slower).

**Stop** under Annotate ends a run at once and removes what it made: the annotation open before stays as it was. The
same Stop button appears under Run PCA + UMAP, Run niches, the Signaling run and Import (a stopped run puts back the
embedding or niches it would have replaced; a stopped import removes the samples it wrote).

The methods themselves are described in [Annotation methods](#6-annotation-methods).

### The map

- The result is drawn on the H&E: pan and zoom down to single cells.
- **The overview** at the bottom right shows the whole sample as the map does (the image, the cells in their colors
  and the region outlines), with the box of the current view.
- **Color by**: cell type, lineage, region, changes, conflicts, confidence, clusters, niche, trajectory, or **gene
  expression…** (type any panel gene: cells are colored by their transcript count, 0 / 1 / 2–3 / 4–7 / 8+, and a
  clicked cell shows its count).
- Click a legend entry to show only that type; click a cell for its label, alternatives and region.
- Opacity, dot size, hide cells, hide H&E.

### Regions

**Regions are drawn by you**; they start empty.

**Drawing** (the tools and keys of QuPath)

The Regions panel has a row of tools:

| Tool | Key | What it does |
|---|---|---|
| Move | M | Drag moves the map; a click on a region selects it (Shift-click selects more) |
| Rectangle | R | Drag a rectangle (Shift: a square) |
| Ellipse | O | Drag an ellipse (Shift: a circle) |
| Polygon | P | Click points; or turn on the **Pen** (P again) and drag freehand |
| Brush | B | Paint the area |
| Wand | W | Paint, and the area follows the colors of the H&E under the brush (± sets how different a color may be) |

- A Brush or Wand stroke that starts on a region (inside it, or with the brush over its edge) changes **that
  region**, as in QuPath, imported outlines included; elsewhere it starts a **new shape**. Rectangle, Ellipse and
  Polygon add to the selected region, or start a new shape. Hold **Alt** (or turn on **−**) and the tools erase
  instead: brush with Alt over a region to take a piece out of it.
- On a large outline (thousands of points, e.g. from QuPath) a brush stroke shows as a band while you paint and is
  applied when you let go.
- Nothing is saved until **Save** (Enter); a new shape is named then (an existing name adds it to that region).
  **Undo** (Backspace) takes back the last stroke, **Cancel** (Esc) drops the changes. Hold Space and drag to move
  the map while drawing; the wheel zooms. The slider sets the brush size (in screen pixels, so zooming in makes it
  finer).
- With the Polygon tool, letting go of the Pen only pauses: the next stroke goes on from the end, until **Finish**,
  Enter or a click on the first point closes it. The outline is kept while you draw, so **Pause** (or another sample,
  or a reload) keeps it.

**Changing regions**

- A selected region (click it on the map, or its row in the list) offers **Expand…** (grow, or shrink with a negative
  distance, in µm), **Fill holes**, **Merge** (the selected polygons of one region into one) and **Delete**; then
  **Save**.
- A region keeps its name and priors when its outline changes; its cells are recomputed, and Undo (the toolbar's)
  brings the old outline back. A shape cut in two keeps its id on the larger piece.
- A name can have several polygons; regions can be renamed, deleted, zoomed to and outlined.
- Click a region's color square to choose one of 15 colors; outlines, cells colored by region and the legend follow,
  and the color stays with the region when it is renamed.
- The selected region is outlined in bright red, so no region is red: a red color (e.g. QuPath's "Tumor" class in an
  imported file) shows as orange or magenta.

**Which cells a region holds**

- Where outlines overlap the smaller one takes the cells, so a region drawn inside another (e.g. a follicle inside
  the mucosa) keeps its cells whichever was drawn first.
- Cells in no region are "unknown", and region priors never apply to them.
- Refer to regions by name when talking to Claude.

**Drawn by the agent.** Ask Claude in Claude Code to draw one ("draw the mucosa", "outline the lymphoid follicles",
"make a region of tiles F3:H6"). Claude computes the outline from where that tissue's cells lie, checks it on the H&E
and proposes it: the review window opens with the outline on the H&E, **Accept** adds the region (Undo takes it back),
then **edit** adjusts it. This needs the sample annotated (the outline follows where a lineage's cells form tissue);
before that, Claude can offer tiles of the grid.

**Before the first annotation.** A sample not annotated yet opens with its cells gray ("Not annotated"). Regions can
already be drawn, edited and imported there, and **Agent Review** offers the registration (Spatial organization: the
cells on their nuclei, the local correction, torn pieces). The first annotation of the sample starts with these
regions.

**Priors.** The **Priors** list under Regions shows the priors of the annotation (the user's words and how they were
read); **×** removes one (Undo brings it back).

### Corrections and Undo

- Hand-labeled cells are locked and fine-tune the expression model (only the involved cell types can move), so similar
  cells elsewhere update too.
- Every step is a version. **Undo** takes back your last action; clicking it again goes further back (it never undoes
  the undo).
- One action is one version made in the app (an outline of the agent added as a region is one version, however many
  pieces it has), or all the versions the agent made for one request (its changes less than 90 s apart): they go back
  together.

### Export

**Export** saves one of:

- **Cells table** (CSV, one row per cell): `cell_id`, `region`, `cell_type`, `lineage`, `level_1`… (the full type
  path), `confidence`, `expression_label`, `marker_cluster`, `hand_labeled`, and when computed `niche`, `leiden`,
  `umap_1`, `umap_2`; then `x_um`, `y_um`.
- **Registration** (CSV): the affine a11..ty that places the cells on the image now, whether it came from the file,
  Claude's alignment or a correction; it can be imported again.
- **Cell positions** (CSV: `cell_id`, `x_um`, `y_um`, `px`, `py`): every cell's place on the image as the app shows
  it. After a local correction or a piece move this is the complete placement, since local shifts and moved pieces
  are not part of the affine.
- **Regions** (GeoJSON): the outlines of the version in the image's full-resolution pixels, named by region with their
  colors, as QuPath writes annotations. Cells moved into a region without an outline are not in it.
- **Analysis notebook** (Jupyter), described below.

**The analysis notebook** runs without this app.

- It needs the data files and `numpy pandas scipy scikit-learn python-igraph pynndescent umap-learn h5py matplotlib`
  (and `pillow` for Visium HD cells built from outlines).
- It reads the cells from the original files of every platform as the import did (Xenium, h5ad, CosMx, MERSCOPE,
  Visium HD cells or bins, Stereo-seq cells or bins) and applies the sample's filter when it has one.
- Standard steps are plain code: reading the cells, the registration affine and the corrections made in the app
  (replayed in order), PCA / neighbors / Leiden / UMAP, NMF niches, figures, the table. It also replays the
  annotation's last trajectory from its cells.
- The app's annotation algorithm (marker clustering, expression model, priors) is carried as source code in a folded
  cell and replayed with this annotation's marker list, regions, priors, corrections and re-clusterings.
- The results at export time are filled in, and a rerun checks that labels, Leiden clusters and niches come out the
  same: identical on any computer for annotations made with exact neighbors; with approximate neighbors, identical
  with the pynndescent and numba versions named in the notebook; annotations made before v0.2.22 (no choice of
  neighbors) only on the same kind of machine. UMAP coordinates differ slightly (all cores).

### Agent Review

**Agent Review** (the bar on the map) asks which part of the analysis the agent should review: **Annotation**,
**Spatial organization**, **Embedding**, **Niches**, **Compare** or **Signaling**. The last four appear once they have
been run; Compare and Signaling review the analyses you ran in those tabs, which the app logs for the agent.

**How the agent reviews**

1. It looks only at that part and drafts its advice (`draft_advice`).
2. It checks every draft as a skeptical reviewer: it reruns the tool that would confirm or refute it, looks at the
   figure on the H&E and weighs other explanations (spillover, segmentation, batch, too few cells, panel limits), and
   records `check_advice`.
3. Only the drafts that hold up reach the Review panel, each with what was checked ("✓ checked: ..."). The dropped
   ones are folded under "Dropped by the agent after checking" with the reason.

The panel's menu switches between the parts (with the number of open items of each); **Show panel** in the dialog
opens a part without asking the agent. Advice changes nothing until you click; **dismiss** hides it.

**What each part covers**

- *Annotation*:
  - computed marker checks (markers that are weak or unspecific in their cells, negative markers that are detected,
    genes of Unassigned cells) and areas to check (cells of sheet-forming types, such as epithelium or smooth muscle,
    apart from their own lineage), outlined on the map;
  - the agent checks every candidate on the H&E and with marker biology (`record_review`); afterwards only the real
    problems are shown, each with its reason (false alarms are folded away). Verdicts carry over to later versions,
    so only new candidates are checked again; **solve** removes an item you have dealt with;
  - when the labels changed since its last look, it also gives advice on *Cell types* (types to merge, split or add:
    **open in Markers…** for a revised list, **ask the agent** for a re-clustering), *Sample context* and *Technical*
    problems (low-quality tiles, marker spillover next to the epithelium, and the registration, below).
- *Spatial organization*: compartments (outlines of the mucosa, muscularis or lymphoid follicles, computed from where a
  lineage forms tissue and checked on the H&E; **add region**, then **add prior…** for the prior that uses it), and
  gradients and neighbors, found or contradicted (**add prior…**).
- *Embedding*: Leiden clusters against the cell types (mixed clusters, types split over clusters) and, for a joint
  embedding, clusters made of one slice (`embedding_facts`, `show_embedding`); **put in Embed** fills in suggested
  settings.
- *Niches*: what each niche is on the H&E, names (**name the niches**), k and the number of niches (**put in Niche**).
- *Compare* and *Signaling*: what the DE genes, compositions or interactions you looked at mean, the caveats (replicate
  slices, spillover, panel coverage) and the analysis that would settle a question (**put in Compare / Signal**);
  **show on map** colors the map by a gene or a signal.

### Registration check and corrections

Part of the Annotation and Spatial organization reviews.

- **The check** correlates the nuclei on the H&E (hematoxylin) or DAPI with the cell centroids in windows across the
  section and fits an affine.
- **correct registration…** moves every cell onto its nucleus (regions keep their cells; **undo correction** puts them
  back). On SCP2927 section 147328 it found a 0.28% scale mismatch (cells up to 30 um off their nuclei at the corners);
  after the correction the windows agree within 0.5 um.
- **Local correction**: where the H&E is warped, a smooth shift field is added on top of the affine. The 200 um
  windows look up to 25 um away, so a second pass of 600 um windows (up to 100 um) finds larger shifts: those areas
  are moved first and measured again (they show as dashed *large shift* areas in the agent's figure). In a test, a 70
  um bump (Gaussian, sigma 250 um) was corrected to 2 um; a narrower one (sigma 150 um) only to 31 um (correcting
  twice does not help there). A torn or folded piece of tissue cannot be followed by a smooth field.
- Every correction is measured anew on the cells as they are, so correcting again refines the last one.
- **Torn, turned or folded pieces** (cells beside a piece of tissue that holds no cells): outline the piece's cells
  as a region and ask the agent to move it. It searches the rotation, shift (up to 300 um) and mirroring (a folded
  piece shows its back) that put the cells on the nuclei and on the tissue no other cell covers, and proposes the move
  with its figure in the review window (cyan: the cells now, yellow: moved).
  - In tests, pieces 160-600 um across torn off by up to 300 um and turned by up to 180 degrees, mirrored or not, were
    put back within 1-6 um.
  - A fit that is no better than where the cells are, or as good elsewhere, is not proposed: a folded piece out of
    focus shows no nuclei to match.
  - The cells keep their region (the outline stays where it was drawn); the agent can undo the last piece move, and
    **undo correction** undoes all.

### Help and About

**?** (top of the toolbar on the right) has four pages:

- **How to use**: every control.
- **Agent**: what to ask the agent and what comes back (start, data in, alignment, marker lists, priors and
  corrections, regions, analyses, review, reports), with example requests.
- **Methods**: the algorithm of every analysis (data, registration, regions, annotation, embedding, niches, compare,
  signal, trajectory, review).
- **About**: the build, the data folder, whether Claude Code was found and when Claude last used the tools, with a
  button that copies these details. It also holds the settings below.

**Updates** (desktop app). The app asks GitHub for a newer release when it starts and every 6 hours while it stays
open, and offers one by itself at most once a day (**Later** waits until the next day). **Check for updates** in About
asks now and always shows it. Without an internet connection nothing is shown.

- The window shows the notes of every release since the installed one, newest first; **Download and install** fetches the installer in the background.
- **Install and restart** closes the app, puts the new version in its place (macOS: the STAnnotator.app it runs from;
  Windows: the installer, silently, into the same folder) and opens it again. The data folder and the settings stay.
- **Skip this version** stops offering it.
- A copy run from source is pointed to the release page instead.

**Data folder** (desktop app on this computer).

- **Open** beside the data folder shows it in Explorer / Finder.
- **Change…** picks another folder (another drive) for the data from the next start: **Copy my data there** copies the
  samples and annotations (in the background, with Stop; the old folder stays as it was), or **Use it as it is**
  (empty, or the data folder of an earlier installation).
- The choice is kept in `%APPDATA%\STAnnotator\settings.json` (macOS: `~/Library/Application Support/STAnnotator/`,
  Linux: `~/.config/stagent/`); the `STAGENT_HOME` environment variable, when set, takes precedence.

**Usage sharing** and **Daily count** (desktop app): see [Troubleshooting](#usage-sharing).

### Storage, Reset and Remove all

- **Caches** of the data folder (figures, exports, the expression model of earlier versions, saved files waiting for
  their download, converted R files, alignment figures: all made again when needed) stay under a limit, 5 GB by
  default: past it the oldest go first (checked every 10 minutes). About shows their size, sets the limit and clears
  them.
- **Trash** (`trash/` of the data folder: what Reset, Remove sample and the cell filter moved away, one folder per
  removal) has a limit of its own, 10 GB by default: past it the oldest removals are deleted for good. About shows its
  size, sets the limit and empties it.
- **Reset** removes every annotation and chat of the selected sample (moved to `<data folder>/trash/`); in a project,
  of this slice or of all slices, and optionally the PCA / UMAP and Leiden clusters too (they are computed from the
  expression alone and stay otherwise).
- **Remove all** (bottom of the toolbar, desktop app) deletes every sample with all its annotations and chats for
  good, to import from scratch (it cannot be undone); marker lists stay.

## 5. Analyses

The steps of the **Analysis** tab after Annotate. Each has **Ask Agent** to interpret its result.

### Embedding (Embed)

- **Run PCA + UMAP** computes a PCA, a UMAP and Leiden clusters of the sample from expression only (shared by all its
  annotations). About 40–50 s for 356,629 cells: the kNN graph comes from the fast GPU / pynndescent search and is
  handed to UMAP, which starts from PC1/PC2 and runs on all cores.
- **Show UMAP** opens a window on the right (UMAP or PC1/PC2, colored as the map, by cell type, Leiden cluster or
  niche).
- **Select** (or Shift-drag) draws around cells: they light up on the spatial map and the rest turn gray; **Show on
  map** zooms there.
- Leiden clusters are also a *Color by* option.

### Niche

Places where cells of similar make-up live together.

- **How**: every cell gets an identity, its weights on NMF gene programs (*Identity: NMF programs*, independent of the
  annotation) or its cell type. The identities are averaged over the cell and its *k* nearest cells and clustered
  (k-means) into the chosen number of niches.
- **Run niches** uses one k everywhere (**k: Set k**). With **k: Agent chooses** the same button asks the agent to
  choose k per drawn region from cell density and cell types (or a radius, so neighborhoods cover a similar stretch of
  tissue in dense and sparse areas), run it and name the niches.
- **Result**: *Color by: niche* (names in the legend) and a `niche` column in Export CSV.
- On SCP2927 section 147328 (356,629 cells) a run takes about 6 s (1.5 s once the NMF is cached), and the 8 niches
  match the tissue layers (surface epithelium, crypt base, lamina propria, lymphoid follicles, muscularis mucosae,
  muscularis propria, submucosal vessels and stroma).

### Compare

For this slice, or all slices of a project pooled.

**DE genes**: the genes that differ between one group and the other cells (or another group).

- A group is a cell type, lineage, region, niche, Leiden cluster, the cells selected in the UMAP window, a slice, or a
  value of a slice metadata column such as condition. **Only within** restricts both sides (e.g. fibroblasts in the
  mucosa vs fibroblasts elsewhere).
- Test: Wilcoxon rank-sum on log-normalized expression (up to 20,000 random cells per side; means and detection on all
  cells); genes ranked by AUC with log2 fold change and FDR.
- With pooled slices, the change is checked in every slice that has both groups, or, when the groups are different
  slices, every slice of the group against every slice of the reference. With one slice per side the app says that a
  slice effect cannot be told from a group effect.
- Group **All** runs this for every group of the grouping at once, each against the other cells (groups under 30 cells
  skipped; up to 5,000 cells per side in the test, so that many groups take seconds to a minute or two): the top up
  and down genes of each group, and every group's genes in Save CSV (column `group`).

**Top genes**: every group against the rest (what sets each niche or cluster apart).

**Composition**: how one grouping splits over another (cell types per niche, niches per region, niches per slice or per
condition, with every slice's own bar), as shares or as log2 enrichment over independence.

**Dot plot** (DE genes and Top genes): one row per group, dot size = share of cells detecting the gene, color = its
mean expression scaled per gene. For one DE comparison only its two sides (the group and the reference or all other
cells). The number in the window sets how many genes the table and the plot show (DE: 5–60, the plot the top N up and
down; DE of All: 1–30 per group; Top genes: 1–15 per group).

**Buttons**

- **Save figure…** writes the plot as PDF; **Save CSV…** keeps the full table.
- A click on a gene colors the map by it.
- **Last result** reopens the last result of the analysis chosen in the tab (DE genes, Top genes or Composition).
- **Ask Agent** asks the agent to interpret the result (it calls the same analysis through `de_genes`,
  `group_markers` or `composition`).

On SCP2927 (559,539 cells, 480 genes) a DE run takes 1-2 s and the top genes of 8 niches about 10 s.

### Signal

Cell-cell communication, for this slice or all slices of a project pooled.

**Ligand-receptor pairs**

- CellChatDB v2 (pathways and the kind of signaling), CellPhoneDB v5 (human) and the OmniPath consensus that LIANA uses
  (CellPhoneDB, CellChatDB, ICELLNET, connectomeDB2020, CellTalkDB), merged into `stagent/data/lr_human.csv.gz`
  (6,521 pairs) and `lr_mouse.csv.gz` (5,792) by `python -m stagent.lr_resource`; or one of them, or your own CSV
  (agent tools).
- A pair counts when every gene of the ligand and receptor is in the panel (complexes have several subunits; a
  complex's level is its lowest subunit): 109 pairs for the 480-gene Xenium IBD panel, about 1,200 for the 990-gene
  MERFISH mouse colon panel.
- Small-molecule signals (an enzyme stands in for the molecule) are off unless switched on.

**Scoring** at single-cell resolution: neighbor pairs are actual cells, never across slices.

- Contact signals: cells within 15 µm. Secreted and ECM signals: the 30 nearest cells within 50 µm, weighted by a
  Gaussian of the distance.
- The score of sender group A → receiver group B is the mean of ligand (sender) × receptor (receiver) over their
  neighbor pairs (log-normalized expression).
- It is tested against random group labels (the CellPhoneDB idea) with a z-score that counts distinct cells rather
  than pairs (on SCP2927 it is more conservative than 100 label permutations), FDR over all tests. Tested only when
  the ligand is detected in 10% of A and the receptor in 10% of B.
- Colocalization above 1 means the expressing cells sit together more than the group averages imply.

**Results**

- **Interactions**: score, specificity, colocalization, detection, neighbor pairs, FDR.
- **Bubble plot**: the strongest interactions over the sender → receiver pairs (color: score, size: FDR, hollow: not
  significant).
- **Network**: sender × receiver.
- **Circle plot**: arrows from sender to receiver groups as wide as their number of significant interactions (the
  number sets how many links).
- **Pathways**.
- **Save figure…** writes a plot as PDF.
- **show on map** colors the cells by their part in one signal: sending it (its ligand × the receptor around it;
  blue), receiving it (its receptor × the ligand around it; red), both (purple), or only weakly (above / below the
  median of the cells with that signal); a legend entry shows one role alone. On SCP2927 section 147328, CCL19 → CCR7
  is received in the lymphoid follicles.
- **Compare groups** runs it per condition, slice, region or niche (largest differences first, as a table or a bubble
  plot with one column per group; with one slice per group the app says a slice effect cannot be told apart).
- **Last result** reopens the last Interactions or Compare groups result.

Agent tools: `cell_communication`, `compare_communication`, `show_communication` (the signal on the H&E). A run on
356,629 cells takes about 10 s.

### Trajectory

The order of cells and the genes that change along it, for this slice or all slices of a project pooled, optionally
only within a subset (**Only within**, e.g. one lineage).

**Methods**

- **Spatial axis**: each cell's place between two groups of cells of any grouping (cell types, regions, niches), t =
  d_From / (d_From + d_To) with d the distance to the nearest cell of the group in the same slice: 0 at From, 1 at To.
  Without To, the distance to From in um. E.g. the crypt axis from stem cells to colonocytes, or portal to central
  regions in the liver. Under a second for 100,000 cells.
- **Pseudotime**: diffusion pseudotime (Haghverdi et al. 2016): a diffusion map of the neighbor graph of the cells'
  PCA (15 neighbors, adaptive Gaussian kernel, density-normalized, 10 components), and t = the diffusion distance from
  the root group's mean, scaled to 0..1.
  - Up to 1,000,000 cells the graph holds every cell; above, a random sample (the same for every root), and the other
    cells take the mean t of their neighbors that have one. Cells the root does not reach get none.
  - About 3 s for 97,000 epithelial cells, 17 s for all 356,000 cells of a section (with or without space). The PCA,
    graphs and diffusion maps are kept, so another root, weight of space or composition on the same cells takes about
    a second (5 s for 356,000 cells).
- **Expression + space**: the same on W = (1 - w) W_expression + w W_space, the spatial graph joining each cell's 15
  nearest cells in the tissue; **Weight of space** w from 0 (expression only) to 1 (the tissue alone).

**Genes along it**: each gene's correlation with the rank of t (rising, falling) and its mean in 20 bins of equal cell
numbers; the share of its variance between bins (eta squared) finds genes that peak in the middle and fall again on
both sides (transient).

**Results**

- **Gene heatmap**: z-scores of the bin means, genes in the order of their peak.
- **Genes**: r, eta squared, peak t, % cells; a click colors the map by the gene.
- **Composition** along it: shares of a grouping per bin.
- **Save CSV…** (every gene), **Save figure…**, **Ask Agent**.
- **Show on the map** colors the cells in ten steps of t (Color by: trajectory: spatial axis, pseudotime or
  expression + space).
- The three methods are kept apart: each has its own last result (**Last result** opens the one of the method chosen),
  its own color on the map, and its own section in the exported notebook, which replays it from its cells with the
  same code. A saved annotation carries them.

On SCP2927 section 147328 the epithelial cells ordered from stem cells rise in APOB, FABP2, CYP3A4, fall in OLFM4, GPX2,
LGR5, and MKI67, TOP2A, RRM2 peak in between (the transit-amplifying zone). An axis or a root sets the order: a
gradient along a spatial axis is position, not time, and spillover from neighboring cells follows the axis too.

Agent tool: `trajectory`.

### Report

**Report** (the bar on the map) shows the annotation report the agent writes (**Ask Agent to write it**): sample,
method and inputs, composition, compartments, quality and open questions, from the tools' facts. **Save…** writes it
as Markdown.

## 6. Annotation methods

### Load existing

A CSV (or TSV) with one row per cell: the cell id as in the sample (h5ad obs names or the Xenium `cell_id`), the cell
type, and optionally a lineage column; the dialog guesses the columns and shows how many cells match.

- Its labels become the annotation (cells not in the file are Unassigned).
- An expression model trained on them adds a confidence per cell and lets priors and corrections act.
- A CSV exported by the app can be loaded back.

Example: `share/celltype.csv` holds the published cell types of SCP2927 section 147328 (`cell_id, cell_type,
lineage`).

### Single-cell reference (like CellTypist)

Give a labeled single-cell h5ad (or a Seurat / SingleCellExperiment `.rds`) of the same tissue, its cell type column,
and optionally a coarser group column (used as lineages for colors and priors).

- A multinomial logistic regression is trained on up to 3,000 reference cells per label, on the genes the reference
  shares with the panel, and labels every cell. Both sides are normalized on the shared genes only, so a
  whole-transcriptome reference fits a panel of a few hundred genes.
- Raw counts are read from `layers["counts"]`, `raw.X` or `X` (log-normalized values are converted back), and gene
  names from the var index or the var column that matches the panel best (e.g. `feature_name` next to Ensembl ids).
- *Majority voting* over-clusters the sample and lets each cluster agree on one label.
- Hand corrections fine-tune the classifier (replaying the reference cells).

Example: a reference built from three other SCP2927 sections labeled section 147328 with 80% agreement on 56 fine types
and 97% on the five lineages (98% with majority voting); a reference from another platform (scRNA-seq) will agree
less.

### Marker genes

A hierarchical marker list. The text format (indent two spaces per level; the CSV in [Marker lists](#2-marker-lists)
is equivalent):

```
# Name: GENE+ GENE- A|B+    (+ expressed, - not expressed, A|B+ = any of A or B)
Epithelial: EPCAM+ KRT8+
  Absorptive: CA1|SLC26A3+ SPINK4-
    BEST4+ colonocytes: BEST4+ OTOP2+
T and NK cells: CD3E|NKG7+ EPCAM-
  CD4 T cells: CD4+ CD8A-
```

Claude can draft a list from the sample's panel genes and checks it before use.

**Top-down annotation.** At each level the parent's cells are re-clustered (PCA, kNN graph, Leiden; the resolution is
chosen per level from 0.5 / 1 / 2 / 4 as the one whose clusters the markers assign best). Each cluster is scored
against the child types' markers and assigned only if the best score and its margin are high enough.

**Cell types are the leaves of the tree** (plus "Unassigned"). A type with subtypes is a lineage, usable in priors but
never a label, so the subtypes of a type should cover all its cells (add a catch-all subtype where needed).

**Expression model.** The leaf labels train an expression model (multinomial logistic regression) on the sample's own
cells, which labels the cells the marker step left at a parent level within their lineage.

**Priors** reweight its probabilities:

| kind | condition | example |
|---|---|---|
| `region` | cell lies in one of the named regions | epithelial cells only in the mucosa |
| `depth` | within a distance of the epithelium (estimated from EPCAM/KRT8) | goblet cells within 30 um of the epithelium |
| `neighborhood` | a fraction of the 15 nearest cells are given types | Tfh cells next to B cells |

Each works as `only`, `exclude` or `enrich` with strength `strong` / `moderate` / `weak`.

**Changing the result**

- Local re-clustering ("this cluster is mixed") re-clusters one cluster or level, is previewed, and keeps regions,
  priors and corrections.
- To add or remove types, revise the marker list and Annotate again with "keep my regions, priors & corrections".

### Review

Claude reviews an annotation with three kinds of evidence:

- *Markers*: for every lineage and cell type, how often each listed marker is detected in its cells compared with its
  siblings, flags (weak, not specific, negative marker detected), and the genes that actually separate those cells in
  this sample; Unassigned cells are described the same way.
- *Space*: a whole-sample overview (H&E with a lettered tile grid next to the lineage map) and computed candidate
  areas, i.e. tiles where cells of a type that forms sheets sit apart from their own lineage, at least 3x more than
  expected from cell density.
- *Biology*: Claude's knowledge of the tissue, cell types and histology, which also tells real structures (vascular
  smooth muscle, muscularis mucosae, deep crypts in IBD) from mislabels.

Nothing changes until the user agrees.

Example (SCP2927 section 146685, ulcerative colitis, scored afterwards against the dataset's own labels): the review
found that 15,514 Unassigned cells expressed CD177, CSF3R and S100A12 (neutrophils missing the myeloid markers) and
that the "dendritic cells" had macrophage genes (C1QC, LYVE1, MS4A7). With the revised markers, Unassigned cells fell
to 195, 97% of the dataset's neutrophils were labeled neutrophils (15% before), and lineage agreement rose from 91.4%
to 98.6%.

## 7. Performance

- Everything runs on the CPU. With PyTorch installed (`gpu` extra) the marker step's exact neighbor search uses an
  NVIDIA or Apple GPU (same neighbors, only faster); without a GPU, PyTorch still makes it about 4 times faster than
  numpy.
- A 357k-cell section takes under a minute to annotate with a GPU, about 5 min with exact neighbors on a 12-core CPU
  and about 1 min with approximate ones (3.3 GB of RAM). The Annotate tab shows the estimate for the sample and
  computer at hand.
- A correction or a local re-clustering takes a few seconds.
- The viewer draws the cells with WebGL, so panning, zooming and changing colors are immediate.

## 8. Troubleshooting

### Reporting a bug

Something went wrong with the agent or the app: click the bug icon (top right, between **?** and **›**).

1. Describe the bug (what you did, what happened, what you expected) and optionally give your name or email.
   Screenshots help with anything you can see: add up to 3 with **Add screenshot…**, or drop or paste them into the
   dialog (any size).
2. **What goes in the report** lists each part and lets you leave it out. The report adds, from the last 1 to 30 days
   (3 by default):
   - your conversations with the agent (Claude Code's own records of them, images left out);
   - every call of the app's tools by the agent (arguments, results, or errors with where in the code they happened);
   - the changes it proposed and what became of them;
   - the open annotation's history (versions, priors, corrections, analyses run; only the latest 10 versions) and its
     sample's description (platform, cells, genes, image, registration);
   - errors in the app's window (scripts, failed requests) and the end of the app's logs.

   It holds no expression data or cell tables, and no images but your screenshots; conversations can include file
   paths, gene names and numbers the agent read.
3. Send it:
   - **Upload to the developers** sends it to the developers' Google Drive folder (the installed app knows where; run
     from source, it asks once for the upload link the developers gave you). Each screenshot goes as a file of its own
     next to the report; one file can be up to 35 MB;
   - **Save to file…** keeps it on this computer to send yourself (privately: an issue on GitHub is public).

With the app working on a server, the report is made there and this computer's conversations and log are added to it.

### Usage sharing

During testing, the desktop app asks once whether it may send the developers a report of its use every 3 days.

- The contents are those of a bug report, for the time since the last one; nothing is sent when the app was not used.
- The reports carry a random id of the installation, not the user's name.
- **?** → **About** → **Usage sharing** shows when the last one was sent, turns it off or on, and sends one now.

**Daily count.** Once a day the desktop app sends the developers its random installation id, its version and the
system (macOS or Windows), nothing else, so they know how many installations are in use. It does not depend on usage
sharing; **?** → **About** → **Daily count** turns it off or on.

### Common problems

- **The app does not open**: see its log, `~/Library/Caches/STAnnotator/app.log` on macOS and
  `%LOCALAPPDATA%\STAnnotator\app.log` on Windows (`<data folder>/app.log` when run from source). Run
  `.venv/bin/stagent app` in a terminal to see errors, or `--browser` to use the browser. The web viewer logs to
  `<data folder>/viewer.log`.
- **Claude Code does not open**: it must be installed and logged in (`claude` works in Terminal). By hand:
  `cd ~/stagent/claude_workspace && claude --mcp-config stagent_mcp.json`.
- **Import says too few cells land on the H&E**: check the registration file, `--section`, and that the cell
  coordinates are in micrometers.

## Other clients

- Claude Code without the app: `claude mcp add stagent -e STAGENT_HOME=$HOME/stagent -- $PWD/.venv/bin/python -m stagent.mcp_server`
  (run in this folder). `.claude/skills/spatial-annotation/` holds the agent's working rules for Claude Code.
- Shared server: the viewer binds to localhost only. Others reach it through an SSH tunnel
  (`ssh -N -L 8765:localhost:8765 <user>@<server>`), then open http://localhost:8765. The conflict review page for
  pathologists is at http://localhost:8765/review (verdict keys 1-6, rules form, CSV download).

## On Windows

- The installer and the app are not signed yet: Windows asks once (**More info → Run anyway**).
- **Ask Agent** opens Claude Code in its own console window. A request sent while it is open reopens it with the
  request (`claude --continue`, the same conversation), since text cannot be typed into another console.
- Sample names are limited to 100 characters (Windows' 260-character paths).
- Harmony in a project's **Run PCA + UMAP** uses `stagent/harmony.py`, a NumPy version of harmonypy 2.1.0 (which has
  no Windows build): on two slices of 559,539 cells about 20 s, and the same neighbors as harmonypy up to its
  randomness.
- The app writes its start-up log to `%LOCALAPPDATA%\STAnnotator\app.log`; the data folder chosen in the installer is
  kept in `%APPDATA%\STAnnotator\settings.json`.

## Layout

```
stagent/     Python package
  cli.py           stagent app / make-app / import / import-xenium / samples / viewer / mcp
  desktop.py       the desktop app (pywebview window around the local server)
  claudelink.py    connection to the user's Claude Code; viewer <-> agent handoff files
  importer.py      Xenium files (or output folder), or h5ad / RDS + H&E + registration -> sample bundle
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
  autoshare.py     usage sharing during testing (with the user's consent); the daily count
  report.py        bug reports: contents, upload
  updater.py       the desktop app's updates
  data/            app icon, built-in marker lists (organs/), ligand-receptor pairs
examples/scp2927/  registration of the SCP2927 sections
packaging/         packaging of the macOS and Windows apps
tests/             smoke tests (python tests/smoke_engine.py, tests/smoke_mcp.py expect sample 147328 / 167799; tests/smoke_joint.py runs in a temporary folder)
```
