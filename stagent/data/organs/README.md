# Built-in marker lists of common organs

One file per organ (and shared blocks of immune and stromal cells, included with `@name`), in the app's marker list
format: one type per line, children indented by two spaces. Each line lists more markers than a panel will hold,
most specific first; `markers.organ_draft()` keeps, for one gene panel, the types it can resolve and their first
markers. `GENE!` is a marker that identifies the type on its own; `-GENE` a negative marker. Genes are human
symbols; for a mouse panel they are matched case-insensitively, and mouse-only symbols are listed too (e.g. LY6G,
CD8B1, ADGRE1). Sources: canonical markers of the Human Protein Atlas single-cell types, CellMarker 2.0, PanglaoDB
and organ atlases (Human Cell Atlas, Tabula Sapiens, the Gut Cell Atlas, the Human Lung Cell Atlas, the Human Kidney
atlas, Allen Brain Cell atlas).

`@immune` inserts the lines of `immune.txt`; `Myeloid > Kupffer cells: ...` adds a type under one of an included file.

The Xenium Prime 5K panel leaves out most very abundant genes (ALB, INS, MYH11, KRT5, COL1A1, …); lower-abundance
markers it has follow the classic ones on each line (e.g. smooth muscle MYLK SMTN KCNMB1 SYNM, fibroblasts PI16 MFAP5,
hepatocytes AHSG ARG1 CYP2C8, beta cells PDX1 SLC30A8), so the lists hold up on 5K panels too.

Rules the lists follow (checked 2026-11-07, again 2026-11-14):
- Current HGNC symbols, a previous symbol after the current one where older panels may use it (`LORICRIN LOR`,
  `IRAG2 LRMP`, `CIMAP3 PIFO`, `CBLIF GIF`). Every symbol is a current or previous HGNC symbol, or one of the
  mouse-only symbols listed on purpose.
- Only the first markers in a panel are used (four for a type, five for a group), so a gene shared with a sibling or
  with a neighbouring cell type comes late on the line (e.g. CD14 and S100A8/9 late for monocytes since macrophages
  and neutrophils have them, DES last for smooth muscle since cardiomyocytes have it, PROX1 not decisive in the liver
  since hepatocytes have it).
- A negative marker is one the type does not express in any state (mast cells are CD45+: no -PTPRC), and not a gene
  whose probes also catch a relative it has (FCGR3B: Xenium's probes see FCGR3A too, which CD16+ monocytes express;
  as a negative it emptied the monocytes of a colon section); a type sits
  under a parent whose markers it has (podocytes are PAX8-negative: not under the kidney's Epithelial).
- A marker is judged by its detection rate (cells with a transcript) in the type against its siblings in
  single-cell atlases of healthy and diseased tissue (colitis, Crohn's disease, interstitial lung disease, lung
  adenocarcinoma, follicular lymphoma), since most samples are diseased; a gene detected as often in a sibling as in
  the type goes late or out (e.g. PECAM1 in myeloid cells, TYROBP in NK and mast cells, IL7R in CD8 T cells and ILCs,
  CD14 in colon epithelium, CD27 in germinal center B cells), and a gene of a cell state (VCAM1 in injured proximal
  tubule) is not a type's first marker.
- Group lines end with alternatives for panels without the classic genes (epithelium CDH17 ELF3 LGALS4 KRT20
  CLDN4 …; smooth muscle FLNC), so a type is not dropped when its first markers are missing.

Checked 2026-11-14 against public single-cell atlases (CELLxGENE: Tabula Sapiens lymph node, large and small
intestine, lung, heart, kidney, liver, pancreas, skin, breast, stomach; colitis / PSC colon; Crohn's colon; tonsil B
cells; follicular lymphoma; the LNCS lymph node atlas; interstitial lung disease; lung adenocarcinoma; Tabula Muris
Senis), the Human Protein Atlas single cell types, and drafted for the Xenium Prime 5K, Xenium 480 (colon), CosMx
960, CosMx 6K and mouse MERFISH 990 panels. Each change rests on the atlases; three were found by looking at errors on
the data below (mature DCs of the lymph node, -VIL1 for immune cells and COL1A2 before PI16 for fibroblasts on the
mouse panel). Measured against the authors' labels and the protein stains (before -> after): SCP2927 Xenium colon,
three sections, cell-type agreement 0.852 / 0.763 / 0.867 -> 0.854 / 0.764 / 0.867, macro-F1 0.736 / 0.651 / 0.689 ->
0.740 / 0.664 / 0.689; Xenium 5K lymph node agreement 0.603 -> 0.610, macro-F1 0.503 -> 0.524 (mature DCs typed as
NK cells 51% -> as mregDC 70%), lineage 0.688 -> 0.695; mouse colon MERFISH, four slices, lineage 0.727 / 0.887 /
0.596 / 0.468 -> 0.883 / 0.892 / 0.922 / 0.510 (Cd14, detected in mouse colon epithelium, no longer makes epithelium
monocytes); CosMx lung and colitis colon, PanCK / CD45 protein AUROC of the epithelial / immune calls unchanged
(0.976 / 0.837 -> 0.973 / 0.836; 0.974 / 0.775 -> 0.974 / 0.774).
Brain (2026-11-19): rewritten from cortex classes only (excitatory / inhibitory neurons, four interneuron types,
glia) to the cortical layer types (L2/3 IT, L4/5 IT, L5 ET, L5/6 NP, L6 CT, L6b), hippocampal CA and dentate gyrus
neurons, thalamic neurons, cerebellar granule cells, Purkinje cells and molecular layer interneurons, striatal D1 / D2
medium spiny neurons, SNCG interneurons, dopaminergic, serotonergic and noradrenergic neurons, the oligodendrocyte
lineage (OPCs, newly formed oligodendrocytes), ependymal cells, choroid plexus epithelium and VSMCs; with alternatives
for panels without the classic genes (e.g. C1QL1 for OPCs, CLDN11 ERMN for oligodendrocytes, SALL1 for microglia,
FGD5 for endothelium, PPP1R17 SKOR2 for Purkinje cells). Checked against the Allen Brain Cell atlas MERFISH
Zhuang-ABCA-4 (three mouse brain sections, 162,578 cells, 1,122 genes, the atlas's subclasses), scored over the cells
of subclasses the list has: agreement 0.147 -> 0.947, macro-F1 0.044 -> 0.648 (the old list had no neuron subclass,
and on this panel lost oligodendrocytes, OPCs, microglia and endothelium for want of their classic genes). On the
MOSTA Stereo-seq adult mouse brain (cell bins, no labels) the layer types lie in their order across the cortex and the
dentate gyrus, CA and thalamic types in their places. Known limits: L5 IT and L6 IT neurons have no type (they go to
L4/5 IT, L6 CT or L6b), SST and VIP need SST / VIP in the panel, and Stereo-seq's high background of Pvalb and Sst
makes too many PVALB / SST interneurons there (as with the old list).
Earlier (2026-11-07), against the authors' labels of SCP2927 (Xenium, 480 genes): the small intestine list on an ileum
section and the colon list on a colon section; and the colon list on a mouse colon MERFISH slice (990 signaling genes,
at the lineage level). Agreement with the authors, before -> after that check: ileum 88.0% -> 88.0%, colon 89.0% -> 89.1%
(plasma cells 66% -> 76%), mouse colon 15% -> 35% (smooth muscle found again).
