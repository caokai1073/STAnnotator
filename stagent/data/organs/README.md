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
