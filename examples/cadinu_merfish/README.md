# Example: MERFISH mouse colon (Cadinu et al., Cell 2024)

MERFISH of the mouse colon during DSS colitis (990-gene panel focused on signaling: receptors, ligands, cytokines).
`cadinu_merfish_D0_m3a_1_slice_2.h5ad` (8 MB) is one healthy slice of the authors' `adata.h5ad` (18 GB, 52 slices):
32,002 cells x 990 genes, raw counts in `X` (from the source's `raw.X`; its `X` holds normalized log values), the
source's cell metadata in `obs` (`x`, `y` in um, also in `obsm["spatial"]`; the authors' labels `Tier1`-`Tier3`,
`Leiden_neigh`, kept to compare with, ignored by the import). No H&E: the cells are shown on a blank background.

In the app: **Import... > New sample > h5ad**, choose the file, name it. Or:

```bash
.venv/bin/stagent import --name cadinu_D0_m3a_s2 --h5ad examples/cadinu_merfish/cadinu_merfish_D0_m3a_1_slice_2.h5ad
```

(From the full file: `--h5ad /path/to/Cadinu_MERFISH_DSS/adata.h5ad --section-column Slice_ID --section
062921_D0_m3a_1_slice_2`; the importer reads `raw.X` there.) Then load `mouse_colon_merfish_markers.txt` in
**Markers... > Import...** and click **Annotate**.

`mouse_colon_merfish_markers.txt` is a marker list for this panel. The panel lacks most classic markers (Epcam,
Myh11, Acta2, Col1a1, Muc2, ...), so it uses what the panel has: Vil1 for epithelium, Flnc / Tnnt2 for smooth
muscle, Col1a2 / Dpt / Pdgfra for fibroblasts, Itgb2 / Rac2 / Ptprc for immune cells (Ptprc alone is detected in
few immune cells here). It was written from gene knowledge and the top genes of unsupervised Leiden clusters,
without the dataset's labels. Compared offline with the authors' Tier1 classes on this slice: 73.7% of cells in the
same lineage (epithelium 80%, immune 79%, smooth muscle 72%, fibroblasts 71%, endothelium 61%), 16% Unassigned;
stem and goblet cells are not separated from transit-amplifying cells (Lgr5, Spdef are rarely detected), and
T cells and neutrophils are not found. A review by Claude (Review annotation) is the next step.
