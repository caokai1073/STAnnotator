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

Rules the lists follow (checked 2026-11-07):
- Current HGNC symbols, a previous symbol after the current one where older panels may use it (`LORICRIN LOR`,
  `IRAG2 LRMP`, `CIMAP3 PIFO`, `CBLIF GIF`). Every symbol is a current or previous HGNC symbol, or one of the
  mouse-only symbols listed on purpose.
- Only the first markers in a panel are used (four for a type, five for a group), so a gene shared with a sibling or
  with a neighbouring cell type comes late on the line (e.g. VCAN after CCR2 for monocytes, DES last for smooth
  muscle since cardiomyocytes have it, PROX1 not decisive in the liver since hepatocytes have it).
- A negative marker is one the type does not express in any state (mast cells are CD45+: no -PTPRC), and not a gene
  whose probes also catch a relative it has (FCGR3B: Xenium's probes see FCGR3A too, which CD16+ monocytes express;
  as a negative it emptied the monocytes of a colon section); a type sits
  under a parent whose markers it has (podocytes are PAX8-negative: not under the kidney's Epithelial).
- Group lines end with alternatives for panels without the classic genes (epithelium CDH17 ELF3 LGALS4 KRT20
  CLDN4 …; smooth muscle FLNC), so a type is not dropped when its first markers are missing.

Checked against the authors' labels of SCP2927 (Xenium, 480 genes): the small intestine list on an ileum section and
the colon list on a colon section; and the colon list on a mouse colon MERFISH slice (990 signaling genes, at the
lineage level). Agreement with the authors, before -> after this check: ileum 88.0% -> 88.0%, colon 89.0% -> 89.1%
(plasma cells 66% -> 76%), mouse colon 15% -> 35% (smooth muscle found again). The other organs were checked against the
Human Protein Atlas single cell types and organ atlases.
