# celltypeagent

Interactive cell type annotation for spatial transcriptomics (Xenium) on the H&E image. Cells are annotated
from a hierarchical marker gene list (no labeled reference needed). You then refine the result in plain language
with Claude: spatial priors ("epithelial cells are only in the mucosa"), corrections ("these cells are plasma
cells"), or requests to split a mixed cluster. Claude turns each request into an explicit constraint, correction
or re-clustering step; the probabilistic model does the annotation, flags conflicts for review, and versions
every step. Everything runs on your own computer; images and expression data never leave it.

## 1. Install

Requirements: macOS (Apple silicon or Intel), Linux or Windows; about 8 GB of RAM; [Claude Desktop](https://claude.ai/download).

1. Install [uv](https://docs.astral.sh/uv/) (it also installs Python for you):

   ```bash
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. Unzip the downloaded `STCellTypeAgent-main.zip` somewhere permanent (for example into `~/Applications`; Claude Desktop is
   pointed at this folder, so do not move it afterwards) and install it:

   ```bash
   cd ~/Applications/STCellTypeAgent-main
   uv venv .venv --python 3.12
   uv pip install --python .venv/bin/python -e .
   ```

   On Windows use `.venv\Scripts\python` and `.venv\Scripts\celltypeagent` in place of `.venv/bin/...`.

Your data and annotations are kept in one data folder, `~/celltypeagent` by default (set `CELLTYPEAGENT_HOME`
to use another one, before importing and before step 3).

## 2. Import a sample

Each sample (tissue section) needs three inputs:

| input | what it is |
|---|---|
| h5ad | AnnData file(s) with raw counts (layer `counts`, or `X`) and cell centroids in micrometers (`obs` columns `x_centroid`/`y_centroid` or `x`/`y`, or `obsm["spatial"]`). Several files (e.g. one per lineage) are merged. Cell type labels in them are ignored. |
| H&E | the H&E image of the section as a pyramidal TIFF / OME-TIFF |
| registration | how cells map onto the H&E: either an affine CSV with columns `a11,a12,a21,a22,tx,ty` (`px = a11*x + a12*y + tx`, `py = a21*x + a22*y + ty`; cell um to H&E full-resolution pixels; with several rows, the row whose other columns contain the `--section` value is used), or a Xenium Explorer alignment matrix (3x3 CSV, H&E pixels to morphology pixels of 0.2125 um) |

```bash
.venv/bin/celltypeagent import --name <sample_id> --h5ad cells.h5ad --he he.ome.tif --registration alignment.csv
```

Options: `--section-column donor_id --section 147328` when the h5ad files hold several sections;
`--region-column histopathology` to keep a region annotation from the h5ad as an optional starting point;
`--counts-layer`, `--x-column`/`--y-column` for other layouts; `--description "..."`. The import checks that at
least 90% of the cells land on the H&E (a wrong registration or section fails here) and takes about 15 s.
`.venv/bin/celltypeagent samples` lists what is imported.

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

`ibd_xenium480_draft.txt`, a marker list for the 480-gene IBD panel, is the default list.

## 3. Connect Claude Desktop

```bash
.venv/bin/celltypeagent setup-desktop
```

This adds a `celltypeagent` entry to Claude Desktop's `claude_desktop_config.json` (other entries are kept and a
backup is written next to it). Quit and reopen Claude Desktop; `celltypeagent` then shows in the tools
menu of a new chat. Then talk to it, for example:

- "Which samples do I have?" / "Annotate sample 147328 and open the viewer."
- "Epithelial cells are only found in the mucosa." (priors are previewed; say go ahead to apply)
- "Cluster 7 of the T cells is mixed, split it." / "These cells in region crypt_1 are goblet cells."
- "Show me the biggest conflict group." / "Go back to version 3." / "Export the result."

Claude answers in the language you write in. Allow the tools when Claude Desktop asks.

## 4. The viewer

Claude opens it with the `open_viewer` tool, or start it yourself with `.venv/bin/celltypeagent viewer` and
open http://localhost:8765. Changes made from Claude Desktop show up in the viewer within a few seconds.

1. Pick a sample; earlier annotations are under **Session**. Click **Annotate** for a new one (about 1 min
   for 350k cells). Regions start empty (or from the h5ad's region column, if imported).
2. The result is drawn on the H&E: pan and zoom down to single cells; color by cell type, lineage, region,
   changes, conflicts, confidence or clusters. Click a legend entry to show only that type; click a cell for
   its label, alternatives and region. Opacity, dot size, hide cells, hide H&E.
3. **Regions are drawn by you.** In the Regions panel click **Draw region**, click points on the map, close
   with a double-click, Enter or the first point (Backspace removes a point, Esc cancels), and name it. A
   name can have several polygons; regions can be renamed, deleted, zoomed to and outlined. Cells in no
   region are "unknown" and region priors never apply to them. Refer to regions by name when talking to Claude.
4. Hand-labeled cells are locked and fine-tune the expression model (only the involved cell types can move),
   so similar cells elsewhere update too. Every step is a version; **Undo** reverts the last one.
5. **Review conflicts** opens the conflict review page. **Export CSV** downloads one row per cell
   (`cell_id`, `region`, `cell_type`, `lineage`, `level_1`… the full type path, `confidence`,
   `expression_label`, `marker_cluster`, `hand_labeled`, `x_um`, `y_um`).
6. **Reset** removes every annotation and chat of the selected sample (moved to `<data folder>/trash/`).
7. **Markers…** edits, checks, saves and loads marker lists (saved lists go to `<data folder>/marker_lists/`).

The chat box in the viewer appears only when the [Claude Code](https://claude.com/claude-code) CLI (`claude`)
is installed and logged in; without it, talk to Claude Desktop instead.

## Annotation method

Marker list format (indent two spaces per level):

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
corrections. To add or remove types, revise the marker list and Annotate again with "keep my regions & priors".

## Performance

The expression model uses the GPU when there is one (NVIDIA CUDA or Apple silicon MPS); clustering runs on the
CPU. On a CPU alone a 356k-cell section takes about 55 s to annotate (3.3 GB of RAM), a correction about 3 s,
and a re-clustering about 5 s.

## Troubleshooting

- Claude Desktop does not list `celltypeagent`: run `setup-desktop` again, fully quit Claude Desktop and
  reopen it. Its server log is `~/Library/Logs/Claude/mcp-server-celltypeagent.log` on macOS
  (`%APPDATA%\Claude\logs` on Windows).
- The viewer does not open: see `<data folder>/viewer.log`, or start it by hand (`celltypeagent viewer`); use
  `--port` if 8765 is taken.
- Import says too few cells land on the H&E: check the registration file, `--section`, and that the cell
  coordinates are in micrometers.

## Other clients

- Claude Code: `claude mcp add celltypeagent -e CELLTYPEAGENT_HOME=$HOME/celltypeagent -- $PWD/.venv/bin/python -m celltypeagent.mcp_server`
  (run in this folder). `.claude/skills/spatial-annotation/` holds the agent's working rules for Claude Code.
- Shared server: the viewer binds to localhost only. Others reach it through an SSH tunnel
  (`ssh -N -L 8765:localhost:8765 <user>@<server>`), then open http://localhost:8765. The review page for
  pathologists is at http://localhost:8765/review (verdict keys 1-6, rules form, CSV download).

## Layout

```
celltypeagent/     Python package
  cli.py           celltypeagent import / samples / setup-desktop / viewer / mcp
  importer.py      h5ad + H&E + registration -> sample bundle
  markers.py       marker list parsing, hierarchical Leiden annotation, re-clustering
  model.py         expression model (training, fine-tuning on corrections)
  constraints.py   priors and Bayesian reweighting
  session.py       sessions, append-only versions, corrections, export
  regions.py       user-drawn regions
  conflicts.py     conflict detection and grouping
  mcp_server.py    MCP tools (Claude Desktop / Claude Code)
  webapp.py        viewer, review page, chat box
  data/            default marker list
examples/scp2927/  registration of the SCP2927 sections
tests/             smoke tests (python tests/smoke_engine.py, tests/smoke_mcp.py; they expect sample 147328 / 167799)
```
