# Tissue spatial priors for cell-type annotation

A reference collection of spatial priors per organ, gathered from the literature for imaging spatial
transcriptomics (Xenium, CosMx, MERFISH, Visium HD) on an H&E. The agent reads it with the tool `tissue_priors`:
each organ's "Compartments to draw" (how each compartment looks on the H&E) guide outlining regions on an image
(also an H&E alone), and each prior can be turned into one rule of the app's constraint language
(`stagent/constraints.py`) with `preview_constraint` / `add_constraint`, after the user agrees.

Compiled 2026-11-06 by Claude from web-verified sources. Each organ has its own numbered reference list. Claims
that could not be checked online are marked "(textbook histology, not verified online)". Review a prior before
using it, especially its exceptions.

## How to read a rule

- **Cell-type names** are those of the built-in marker lists (`stagent/data/organs/*.txt`). A group name (e.g.
  `Epithelial`) covers its subtypes. Types marked "(not in the built-in list)" or written `<like this>` need a
  marker list that has them.
- **kind**
  - `region`: the cell lies in a region the user drew and named. The names under "Compartments to draw" are
    suggestions and must match the names the user gives.
  - `depth, max_depth_um=X`: the cell is within X um of the epithelial core (cells dense in EPCAM / KRT8). This is
    useful in the gut and breast, weak in lung, and not meaningful in liver, pancreas, kidney, heart, skin, brain
    and lymph node (see each organ's notes).
  - `neighborhood, neighbor_types=[...], min_frac=f`: at least fraction f of the 15 nearest neighbors are of
    those types (typed by expression only).
- **mode**
  - `only`: the types are unlikely where the condition is false.
  - `exclude`: the types are unlikely where the condition is true.
  - `enrich`: the types are more likely where the condition is true.
- **strength**
  - `only` and `exclude`: strong x0.01, moderate x0.1, weak x0.5.
  - `enrich`: strong x10, moderate x3, weak x1.5.
- **Depth values** are starting points: calibrate them on each sample, for example from the depth of cells in a
  drawn muscularis region.

## Using priors with care

- Priors describe normal tissue. Tumors, inflammation, metaplasia, fibrosis and injury break many of them; each
  prior lists the cases it knows.
- A prior that fires on many confident cells is a signal, not only an error to suppress. It can mean
  misregistration, segmentation spill-over (e.g. plasma-cell transcripts in epithelium), a wrong region outline, or
  real biology. Check the conflicts before raising a prior's strength.
- Prefer moderate or weak strengths unless the anatomy is unambiguous (e.g. podocytes only in glomeruli).

## Contents

- Gastrointestinal tract: colon, small intestine, stomach
- Liver, pancreas, kidney
- Lung, heart, breast, skin
- Brain (cortex), lymph node and tonsil
- General cross-organ priors for immune and stromal cells


---

## Gastrointestinal tract: shared notes

Cell-type names are taken from `colon.txt`, `small_intestine.txt` and `stomach.txt`. Each of these includes
`@immune`, `@stroma` and `@enteric`, which add Glia, Neurons and ICC.

**How to read the depth suggestions (applies to all three organs)**
- `depth_um` is the distance from a cell to the nearest "epithelial core" cell. A cell is a core cell when at least
  30% of its 30 nearest neighbours are EPCAM/KRT8-positive. As a result:
  - Epithelium has depth ~0.
  - The narrow lamina propria between crypts, inside villus cores and between gastric glands also has depth ~0 to
    a few tens of um.
  - Depth increases below the crypt or gland bases: muscularis mucosae, then submucosa, then muscularis propria,
    then serosa.
- Depth cannot separate positions along the crypt-villus or pit-gland axis, because all of those cells are at
  depth ~0. Zonation priors therefore need drawn regions.
- Layer thicknesses:
  - Human colon: crypts are about 370 um (proximal) to 510 um (distal) long [COL 5]. The whole wall between the
    teniae is about 0.5 to 3.5 mm thick (mean about 1.7 mm), and the muscularis propria is its thickest layer
    [COL 6].
  - Mouse: the whole wall is only a few hundred um thick. No mouse per-layer value was verified online.
- The `max_depth_um` values below are deliberately conservative: they still hold in mouse. In human tissue, larger
  values can be used. Calibrate them on the sample by checking the depth of cells inside a muscularis region the
  user has drawn.
- Lymphoid aggregates and Peyer's patches create large epithelium-free zones inside the mucosa, so their cells get
  a large depth. Never use a depth `only` rule for lymphocytes.

---

## Colon

**Compartments to draw**
- **mucosa**: the crypt layer, from the lumen down to the bottom of the crypts. Straight, parallel, test-tube
  crypts full of pale goblet cells, with no villi. The lamina propria between the crypts is pink-purple and
  cellular.
- **muscularis mucosae** (optional; otherwise include it in the mucosa): a thin pink band of smooth muscle
  directly below the crypt bases.
- **submucosa**: loose, pale pink connective tissue with large vessels and small ganglia, between the
  muscularis mucosae and the thick muscle.
- **muscularis propria**: thick, dense eosinophilic smooth muscle. It has an inner circular layer and an outer
  longitudinal layer, and the longitudinal layer is gathered into the teniae.
- **myenteric plexus** (optional): pale, small clusters of large neurons with glia, lying between the circular
  and longitudinal muscle layers.
- **serosa**: the outer adventitia or subserosa, with fat and a mesothelial surface.
- **lymphoid follicle**: a dense, basophilic (blue) nodule that often straddles the muscularis mucosae. It
  sometimes has a pale germinal center, which can be drawn separately as **germinal center**.
- **crypt base** and **surface epithelium** (optional, at high zoom): the lower third of the crypts, and the top
  of the crypts together with the luminal surface.

**Priors**

- **COL-01. Epithelial cells (all subtypes) lie in the mucosa, not in the submucosa, muscularis propria or serosa.**
  - Rule: `kind=region, types=[Epithelial], regions=[mucosa], mode=only, strength=strong`
  - Without drawn regions: `kind=depth, types=[Epithelial], max_depth_um=30, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Invasive adenocarcinoma, and tumour cells inside lymphatics or vessels.
    - Pseudoinvasion (misplaced epithelium) in pedunculated polyps.
    - Lymphoglandular complexes: normal crypts herniating through the muscularis mucosae into a lymphoid
      follicle.
    - Colitis cystica profunda, and endometriosis.
    - Tangential sections.
    - Spill-over of EPCAM/KRT8 transcripts into the neighbours of crypt cells, which causes false-positive
      Epithelial calls in the lamina propria. This rule does not remove those calls.
    - Do not use this rule on tumour samples.
  - Confidence: high in normal tissue (strong). Low in neoplasia.
  - Evidence: [1, 3, 4]; layer anatomy is textbook histology, not verified online.

- **COL-02. Stem cells (LGR5+/OLFM4+) sit at the crypt base, and Cycling TAs sit in the lower half of the crypt.**
  - Rule:
    - `kind=region, types=[Stem cells], regions=[crypt base], mode=enrich, strength=moderate`
    - `kind=region, types=[Stem cells, Cycling TAs], regions=[surface epithelium], mode=exclude, strength=moderate`
  - Without drawn regions: none. Depth is ~0 along the whole crypt.
  - Exceptions and caveats:
    - Crypts cut transversely show only one level of the crypt.
    - Regeneration, hyperplastic and serrated polyps, and adenomas expand or shift the proliferative zone. In
      adenomas, proliferation can reach the surface ("top-down" morphology).
    - LGR5 is detected poorly by most panels.
    - Some MKI67+ cells are germinal-center B cells and are not Cycling TAs (see COL-03).
  - Confidence: high biologically, medium in practice because crypt-base regions are laborious to draw.
    Strength: moderate.
  - Evidence: [1, 3, 19]

- **COL-03. Differentiated, mature absorptive cells are at the crypt top and the luminal surface. This covers
  Colonocytes (CA1, SLC26A3, AQP8, CEACAM7) and especially BEST4+ colonocytes.**
  - Rule:
    - `kind=region, types=[BEST4+ colonocytes, Colonocytes], regions=[surface epithelium], mode=enrich, strength=moderate`
    - `kind=region, types=[BEST4+ colonocytes], regions=[crypt base], mode=exclude, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - BEST4/OTOP2 cells are lost in IBD and in colorectal cancer [2], so absence is a real finding.
    - SPIB also marks M cells (not in the colon list).
    - Immature colonocytes (progenitors) do occur in the mid crypt, so do not use an `only` rule.
  - Confidence: medium-high (moderate).
  - Evidence: [2, 3, 19]

- **COL-04. Goblet cells, Enteroendocrine cells, Tuft cells and Cycling TAs are scattered inside the epithelium.
  The same call in a cell surrounded by stroma is usually spill-over or a mis-assigned immune cell.**
  - Rule (no regions needed): `kind=neighborhood, types=[Goblet cells, Enteroendocrine cells, Tuft cells, Cycling TAs], neighbor_types=[Epithelial], min_frac=0.3, mode=only, strength=moderate`
  - Without drawn regions: the rule above. A depth version is `kind=depth, max_depth_um=30, mode=only`.
  - Exceptions and caveats:
    - In human colon, goblet cells are about 23-25% of crypt cells [5]. They line whole crypts, so neighbourhood
      fractions are high.
    - Mucinous carcinoma and signet-ring cells float in mucin pools far from intact epithelium.
    - MKI67+ germinal-center B cells and other proliferating immune cells must not become Cycling TAs: this rule
      helps here.
    - TFF3 is shared with Lymphatic EC.
  - Confidence: high (moderate).
  - Evidence: [3, 4, 5]

- **COL-05. Paneth cells (not in the built-in colon list) are normally present only in the proximal colon (caecum
  and ascending colon) of humans, at crypt bases. They are absent from the normal distal colon and from the
  mouse colon.**
  - Rule (only if the user adds "Paneth cells" to the colon list):
    - `kind=region, types=[Paneth cells], regions=[crypt base], mode=only, strength=moderate`
    - In a sample known to be distal colon or rectum, `kind=region, types=[Paneth cells], regions=[mucosa], mode=exclude, strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[Paneth cells], neighbor_types=[Epithelial], min_frac=0.4, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Paneth cell metaplasia in the distal colon is a real and frequent finding in IBD. It was found in 11 of 13
      paediatric UC and 14 of 15 paediatric Crohn's cases, and it also occurs in other chronic colitis [7]. Never
      use a strong exclude.
    - In mouse colon, REG4+ deep crypt secretory cells take the Paneth niche role [8]. With this list they will
      be called Goblet cells.
    - LYZ is not specific to Paneth cells (macrophages express it too) [19].
  - Confidence: high for the anatomy (weak exclude).
  - Evidence: [7, 8, 19]

- **COL-06. Plasma cells (mostly IgA) and most effector immune cells are concentrated in the lamina propria of the
  mucosa. They are sparse in the submucosa and the muscularis propria.**
  - Rule: `kind=region, types=[Plasma cells], regions=[mucosa], mode=enrich, strength=moderate`
  - Without drawn regions: `kind=depth, types=[Plasma cells], max_depth_um=60, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - Basal plasmacytosis (plasma cells between the crypt bases and the muscularis mucosae) is a hallmark of
      chronic IBD. It stays within the mucosa but at a slightly larger depth.
    - Transmural Crohn's disease, diverticulitis and tumours bring plasma cells deep into the wall.
    - Do not use `only`.
  - Confidence: high (moderate as region enrich; weak as depth enrich).
  - Evidence: [4, 9, 18]

- **COL-07. B cells (naive, MS4A1+) are concentrated in lymphoid follicles. In the colon these are isolated
  lymphoid follicles and colonic patches, and they often straddle the muscularis mucosae. Outside follicles the
  lamina propria has few B cells but many Plasma cells.**
  - Rule: `kind=region, types=[B cells], regions=[lymphoid follicle], mode=enrich, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[B cells], neighbor_types=[B and plasma cells, T and NK cells], min_frac=0.6, mode=enrich, strength=moderate`
  - Exceptions and caveats:
    - In IBD the lymphoid network expands strongly, with tertiary-like follicles [10].
    - MALT lymphoma, and follicles cut at their edge.
    - Do not use a depth rule: follicle centres are far from the epithelium.
    - Naive CD4 T cells also gather in the follicle T-zones (see COL-08).
  - Confidence: high (strong when the follicles are drawn).
  - Evidence: [4, 9, 10]

- **COL-08. Naive CD4 T cells (CCR7+, SELL+) sit in the T-cell zones of lymphoid follicles and aggregates. The
  diffuse lamina propria holds mainly Effector/memory CD4 T cells and Tregs.**
  - Rule:
    - `kind=region, types=[Naive CD4 T cells], regions=[lymphoid follicle], mode=enrich, strength=moderate`
    - `kind=region, types=[Effector/memory CD4 T cells], regions=[lymphoid follicle], mode=enrich, strength=weak`
      (not needed if not wanted)
  - Without drawn regions: `kind=neighborhood, types=[Naive CD4 T cells], neighbor_types=[B cells, T and NK cells], min_frac=0.6, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - CCR7 also marks mature (LAMP3+) dendritic cells.
    - Naive T cells pass through vessels in every layer.
    - Inflamed mucosa recruits many T cells of all states.
  - Confidence: medium (weak-moderate).
  - Evidence: [9]; inductive versus effector sites is textbook mucosal immunology, not verified online beyond the
    abstract of [9].

- **COL-09. Intraepithelial lymphocytes are almost exclusively T cells, mainly CD8 T cells and gd T cells. They sit
  between the epithelial cells. The colon has fewer of them than the small intestine.**
  - Rule (no regions needed): `kind=neighborhood, types=[CD8 T cells, gd T cells], neighbor_types=[Epithelial], min_frac=0.6, mode=enrich, strength=weak`
  - Without drawn regions: as above.
  - Exceptions and caveats:
    - Lymphocytic and microscopic colitis raise IEL counts.
    - Segmentation spill-over gives IELs epithelial transcripts, and the list uses -EPCAM in "T and NK cells".
      Check that the IELs are not lost to Epithelial calls.
    - Narrow lamina propria between crypts can give lamina propria T cells high epithelial fractions. Keep the
      strength weak.
  - Confidence: medium (weak).
  - Evidence: [9, 11]

- **COL-10. Smooth muscle cells form the muscularis propria (and the thin muscularis mucosae). Within the mucosa
  proper, only scattered muscularis mucosae fibres and vascular smooth muscle occur.**
  - Rule:
    - `kind=region, types=[Smooth muscle cells], regions=[muscularis propria, muscularis mucosae], mode=enrich, strength=moderate`
    - Optional: `kind=region, types=[Smooth muscle cells], regions=[mucosa], mode=exclude, strength=weak`, only
      when the muscularis mucosae was drawn as its own region.
  - Without drawn regions: none recommended. The muscularis mucosae sits at depth ~0-40 um, just under the crypt
    bases.
  - Exceptions and caveats:
    - Prolapse-type changes and diverticular disease cause muscularis mucosae fibres to extend between crypts.
    - Pericytes and myofibroblasts share ACTA2 and TAGLN.
  - Confidence: high (moderate).
  - Evidence: [6, 16]; layer anatomy is textbook histology, not verified online.

- **COL-11. Neurons (enteric neuronal cell bodies) are in the ganglia of the submucosal plexus and the myenteric
  plexus. The myenteric plexus lies between the circular and longitudinal muscle. The mucosa contains nerve fibres
  but no neuronal cell bodies. Glia occur in all layers, including the mucosa.**
  - Rule:
    - `kind=region, types=[Neurons], regions=[submucosa, muscularis propria, myenteric plexus], mode=only, strength=moderate`
    - `kind=region, types=[Neurons, Glia], regions=[myenteric plexus], mode=enrich, strength=strong`
  - Without drawn regions: `kind=depth, types=[Neurons], max_depth_um=40, mode=exclude, strength=moderate`. In
    human tissue, 100 um is safe.
  - Exceptions and caveats:
    - Do NOT exclude Glia from the mucosa: mucosal enteric glia sit beneath the epithelium and support the
      barrier [12].
    - Hirschsprung disease leaves aganglionic segments with no neurons. Mucosal ganglioneuromas occur in
      ganglioneuromatosis.
    - PRPH, UCHL1 and SNAP25 transcripts in nerve fibres can be assigned to nearby non-neuronal cells.
  - Confidence: high (moderate).
  - Evidence: [12, 13]

- **COL-12. ICC (interstitial cells of Cajal, KIT+ ANO1+) are part of the smooth-muscle syncytium of the
  muscularis propria: around the myenteric plexus, within the muscle layers, and, in the colon, at the submucosal
  border of the circular muscle. They are not in the mucosa.**
  - Rule:
    - `kind=region, types=[ICC], regions=[muscularis propria, myenteric plexus], mode=only, strength=moderate`
  - Without drawn regions: `kind=depth, types=[ICC], max_depth_um=40, mode=exclude, strength=moderate`. In human
    tissue, 150 um is safe.
  - Exceptions and caveats:
    - KIT also marks Mast cells and ILCs, which sit in the mucosa and submucosa. The list already uses -CPA3 and
      -PTPRC, and this prior helps further.
    - ANO1 is expressed by some epithelial and smooth muscle cells.
    - GIST arises from ICC lineage cells and forms masses.
  - Confidence: high (moderate).
  - Evidence: [14]; ICC sub-locations are textbook physiology, not verified online beyond the abstract of [14].

- **COL-13. Lymphatic EC are essentially absent from the normal colonic lamina propria above the muscularis
  mucosae. Colonic lymphatics start at or below the muscularis mucosae. This differs from the small intestine,
  which has lacteals.**
  - Rule:
    - `kind=region, types=[Lymphatic EC], regions=[mucosa], mode=exclude, strength=weak`, drawing the mucosa down
      to the crypt bases only (not including the muscularis mucosae).
  - Without drawn regions: `kind=depth, types=[Lymphatic EC], max_depth_um=15, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Lymphatics appear in the lamina propria in IBD, polyps and neoplasia [15].
    - The listed marker TFF3 is also a goblet-cell gene, so spill-over from goblet cells creates false Lymphatic
      EC calls. Partly for this reason the prior is useful, but keep it weak.
  - Confidence: medium (weak).
  - Evidence: [15]

- **COL-14. Types deliberately left unconstrained across layers: Macrophages, Fibroblasts, Blood vessel EC,
  Pericytes and Mast cells.**
  - Rule: none.
  - Without drawn regions: none.
  - Why no rule:
    - Macrophages sit in the lamina propria, the submucosa and the muscularis. Muscularis and submucosal
      macrophages sit next to the plexuses and vessels [17].
    - Fibroblast states differ by layer (for example, crypt-niche SOX6+/F3+ fibroblasts) [16], but these
      subtypes are not in the built-in list.
  - Confidence: high (no rule).
  - Evidence: [16, 17]

**References**
- [1] Barker N et al., Identification of stem cells in small intestine and colon by marker gene Lgr5, Nature, 2007. doi:10.1038/nature06196, PMID 17934449.
- [2] Parikh K et al., Colonic epithelial cell diversity in health and inflammatory bowel disease, Nature, 2019. doi:10.1038/s41586-019-0992-y, PMID 30814735.
- [3] Elmentaite R et al., Cells of the human intestinal tract mapped across space and time, Nature, 2021. doi:10.1038/s41586-021-03852-1, PMID 34497389.
- [4] Hickey JW et al., Organization of the human intestine at single-cell resolution, Nature, 2023. doi:10.1038/s41586-023-05915-x, PMID 37468586.
- [5] Arai T, Kino I, Morphometrical and cell kinetic studies of normal human colorectal mucosa. Comparison between the proximal and the distal large intestine, Acta Pathol Jpn, 1989. doi:10.1111/j.1440-1827.1989.tb02421.x, PMID 2618659.
- [6] Konovalov DIu, Morphometric characteristics of human colon wall as anatomic basis for microsurgical anastomoses (in Russian; English abstract), Morfologiia, 2007. PMID 17526266.
- [7] Simmonds N et al., Paneth cell metaplasia in newly diagnosed inflammatory bowel disease in children, BMC Gastroenterol, 2014. doi:10.1186/1471-230X-14-93, PMID 24885054.
- [8] Sasaki N et al., Reg4+ deep crypt secretory cells function as epithelial niche for Lgr5+ stem cells in colon, PNAS, 2016. doi:10.1073/pnas.1607327113, PMID 27573849.
- [9] Mowat AM, Agace WW, Regional specialization within the intestinal immune system, Nat Rev Immunol, 2014. doi:10.1038/nri3738, PMID 25234148.
- [10] Buettner M, Lochner M, Development and function of secondary and tertiary lymphoid organs in the small intestine and the colon, Front Immunol, 2016. doi:10.3389/fimmu.2016.00342, PMID 27656182.
- [11] Cheroutre H et al., The light and dark sides of intestinal intraepithelial lymphocytes, Nat Rev Immunol, 2011. doi:10.1038/nri3007, PMID 21681197.
- [12] Furness JB, The enteric nervous system and neurogastroenterology, Nat Rev Gastroenterol Hepatol, 2012. doi:10.1038/nrgastro.2012.32, PMID 22392290.
- [13] Drokhlyansky E et al., The human and mouse enteric nervous system at single-cell resolution, Cell, 2020. doi:10.1016/j.cell.2020.08.003, PMID 32888429.
- [14] Sanders KM et al., Interstitial cells: regulators of smooth muscle function, Physiol Rev, 2014. doi:10.1152/physrev.00037.2013, PMID 24987007.
- [15] Kenney BC, Jain D, Identification of lymphatics within the colonic lamina propria in inflammation and neoplasia using the monoclonal antibody D2-40, Yale J Biol Med, 2008. PMID 18827885.
- [16] Kinchen J et al., Structural remodeling of the human colonic mesenchyme in inflammatory bowel disease, Cell, 2018. doi:10.1016/j.cell.2018.08.067, PMID 30270042.
- [17] De Schepper S et al., Self-maintaining gut macrophages are essential for intestinal homeostasis, Cell, 2018. doi:10.1016/j.cell.2018.07.048, PMID 30173915.
- [18] Smillie CS et al., Intra- and inter-cellular rewiring of the human colon during ulcerative colitis, Cell, 2019. doi:10.1016/j.cell.2019.06.029, PMID 31348891.
- [19] Burclaff J et al., A proximal-to-distal survey of healthy adult human small intestine and colon epithelium by single-cell transcriptomics, Cell Mol Gastroenterol Hepatol, 2022. doi:10.1016/j.jcmgh.2022.02.007, PMID 35176508.

---

## Small intestine

**Compartments to draw**
- **villi**: finger-like projections into the lumen. Each has a columnar epithelium with goblet cells and a
  lamina propria core containing capillaries, a central lacteal and smooth muscle strands.
- **crypts**: the short tubular glands at the base of the mucosa, between the villus bases. Paneth cells at the
  crypt bottom have bright eosinophilic apical granules.
- **mucosa**: villi and crypts together, down to the muscularis mucosae.
- **Peyer's patch** / **lymphoid follicle**: large basophilic follicles, mainly in the ileum. They occupy the
  mucosa and often extend into the submucosa, and the mucosa above them is flattened into a dome without villi.
- **follicle-associated epithelium**: the epithelium covering the dome of a Peyer's patch or isolated follicle.
  It has fewer goblet cells and is infiltrated by lymphocytes.
- **Brunner's glands** (duodenum only): pale, mucous gland lobules in the submucosa.
- **submucosa**: loose connective tissue with large vessels and submucosal ganglia.
- **muscularis propria**: inner circular and outer longitudinal smooth muscle, with the myenteric plexus between
  them.
- **myenteric plexus** (optional): ganglia between the two muscle layers.
- **serosa**: mesothelium-covered outer surface and mesenteric fat.

**Priors**

- **SI-01. Epithelial cells lie in the mucosa (villi and crypts). The one normal exception is the duodenal
  Brunner's glands in the submucosa.**
  - Rule: `kind=region, types=[Epithelial], regions=[mucosa, Brunner's glands], mode=only, strength=strong`
  - Without drawn regions: `kind=depth, types=[Epithelial], max_depth_um=30, mode=only, strength=moderate`.
    Brunner's glands form their own epithelial core, so depth stays ~0 inside them.
  - Exceptions and caveats:
    - Brunner's glands are MUC6+. The SI list has no Brunner's gland type, so they will be called Goblet cells
      or Epithelial.
    - Adenocarcinoma, neuroendocrine tumours, endometriosis and pancreatic heterotopia.
    - Epithelium dragged into lymphoid follicles.
  - Confidence: high in normal tissue (strong).
  - Evidence: [3, 16, 17]; Brunner's glands are textbook histology, not verified online.

- **SI-02. Paneth cells sit only at the crypt base, intercalated with LGR5+ Stem cells. They are not on the villi.**
  - Rule:
    - `kind=region, types=[Paneth cells], regions=[crypts], mode=only, strength=strong`
    - `kind=region, types=[Paneth cells], regions=[villi], mode=exclude, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Paneth cells], neighbor_types=[Stem cells, Paneth cells, Cycling TAs], min_frac=0.2, mode=only, strength=moderate`
  - Exceptions and caveats:
    - In Crohn's ileitis and other chronic injury, Paneth cell hyperplasia and abnormal Paneth cells occur; in
      some disease states Paneth cells are mislocalised higher up the crypt.
    - LYZ (and in mouse Lyz1) also marks Macrophages and Monocytes [18], so keep DEFA5/DEFA6 as the required
      markers.
    - Gene names differ by species: human DEFA5/DEFA6; mouse Defa/Lyz1/Mptx2.
  - Confidence: high (strong).
  - Evidence: [1, 2, 3, 18]

- **SI-03. Stem cells (LGR5+/OLFM4+) are restricted to the crypt base. Cycling TAs fill the lower and middle
  crypt. Neither is on the villi.**
  - Rule: `kind=region, types=[Stem cells, Cycling TAs], regions=[villi], mode=exclude, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Stem cells, Cycling TAs], neighbor_types=[Epithelial], min_frac=0.3, mode=only, strength=moderate`
  - Exceptions and caveats:
    - In coeliac disease and other villous atrophy, crypt hyperplasia lengthens the proliferative zone.
    - OLFM4 is also expressed by a subset of neutrophils.
    - MKI67+ cells in germinal centres are B cells (germinal-center B cells are not in the list).
  - Confidence: high (strong for the region rule).
  - Evidence: [1, 3, 18]

- **SI-04. Enterocytes are the villus cell type, and their programmes are zoned along the villus axis.**
  - Rule:
    - `kind=region, types=[Enterocytes], regions=[villi], mode=enrich, strength=moderate`
    - `kind=region, types=[Enterocytes], regions=[crypts], mode=exclude, strength=weak`
  - Without drawn regions: none.
  - Zonation, mostly measured in mouse [4]:
    - The villus base has an antimicrobial programme (REG3).
    - The mid-villus has sequential carbohydrate, peptide and fat absorption programmes (FABP2, APOA4, APOB),
      with APOA4 and APOB toward the top.
    - The villus tips have an immune-modulatory programme (NT5E/CD73).
    - Single-subtype splits are not in the list.
  - Exceptions and caveats:
    - Mature enterocyte markers vary by segment: FABP6 is ileum-specific, and APOA4/APOB are highest in the
      jejunum [18].
    - Villous atrophy (coeliac disease) and short villi in biopsies.
    - Tangentially cut villi look like islands of epithelium.
  - Confidence: high for the villus localisation (moderate). The tip/base split was shown mainly in mouse.
  - Evidence: [3, 4, 18]

- **SI-05. Goblet cells, Enteroendocrine cells, Tuft cells and BEST4+ enterocytes are scattered single cells
  inside the epithelium. Enteroendocrine hormone content changes along the crypt-villus axis.**
  - Rule (no regions needed): `kind=neighborhood, types=[Goblet cells, Enteroendocrine cells, Tuft cells, BEST4+ enterocytes], neighbor_types=[Epithelial], min_frac=0.3, mode=only, strength=moderate`
  - Without drawn regions: as above.
  - Hormone switching along the crypt-villus axis [5]: GLP1 (GCG) is enriched in crypts and secretin on the villi,
    driven by the BMP gradient. Hormone subtypes are not in the list.
  - Exceptions and caveats:
    - Neuroendocrine tumours.
    - Spill-over of goblet-cell mucin transcripts (TFF3, SPINK4) into the lamina propria.
    - Duodenal BEST4+ cells (CFTR-high) are described in [19].
  - Confidence: high (moderate).
  - Evidence: [3, 5, 19]

- **SI-06. M cells occur only in the follicle-associated epithelium over Peyer's patches and isolated lymphoid
  follicles. They make up about 10% of follicle-associated epithelium cells. They are not in the ordinary villus
  epithelium.**
  - Rule:
    - `kind=region, types=[M cells], regions=[follicle-associated epithelium, Peyer's patch, lymphoid follicle], mode=only, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[M cells], neighbor_types=[B and plasma cells, T and NK cells, Myeloid], min_frac=0.2, mode=only, strength=moderate`
  - Exceptions and caveats:
    - A few villous M cells exist, and inflammation and infection can induce M cells outside follicles [6].
    - SPIB is shared with BEST4+ enterocytes, and CCL20 is also expressed by inflamed epithelium. A strong rule
      therefore mainly protects against false M-cell calls.
  - Confidence: high (strong with drawn regions).
  - Evidence: [6, 7, 8]

- **SI-07. B cells are concentrated in Peyer's patches and lymphoid follicles. Peyer's patches occur mainly in the
  ileum, and the follicles span the mucosa and the submucosa.**
  - Rule: `kind=region, types=[B cells], regions=[Peyer's patch, lymphoid follicle], mode=enrich, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[B cells], neighbor_types=[B and plasma cells, T and NK cells], min_frac=0.6, mode=enrich, strength=moderate`
  - Exceptions and caveats:
    - Never use a depth `only` rule for B cells or T cells: follicle centres can lie hundreds of um from the
      nearest epithelial core and in the submucosa.
    - Lymphoma, and reactive follicular hyperplasia (for example terminal ileum in children).
    - Isolated lymphoid follicles occur along the whole intestine, not only in the ileum [8].
  - Confidence: high (strong).
  - Evidence: [7, 8, 11, 17]

- **SI-08. The small intestine has many intraepithelial lymphocytes. Normal duodenum has about 11 (up to about 25)
  per 100 villus epithelial cells. Nearly all are T cells, mostly CD8 T cells and gd T cells, and there are more of
  them than in the colon.**
  - Rule (no regions needed):
    - `kind=neighborhood, types=[CD8 T cells, gd T cells], neighbor_types=[Epithelial], min_frac=0.6, mode=enrich, strength=moderate`
    - `kind=neighborhood, types=[Naive CD4 T cells, B cells], neighbor_types=[Epithelial], min_frac=0.6, mode=exclude, strength=weak`
  - Without drawn regions: as above.
  - Exceptions and caveats:
    - Coeliac disease, other enteropathies and drugs raise IEL counts above 25 per 100 [10]. That is a real
      finding, so the enrich rule is fine but must not be turned into an exclusion of lamina propria T cells.
    - Mouse IELs are proportionally richer in gd T cells and CD8aa cells.
    - Epithelial spill-over into IEL segments.
  - Confidence: high (moderate).
  - Evidence: [9, 10, 11]

- **SI-09. Plasma cells, Macrophages and most Dendritic cells fill the lamina propria of the villus cores and
  between the crypts.**
  - Rule: `kind=region, types=[Plasma cells], regions=[mucosa], mode=enrich, strength=moderate`
  - Without drawn regions: `kind=depth, types=[Plasma cells], max_depth_um=50, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - Plasma cells also gather around Peyer's patches.
    - In Crohn's disease, plasma cells and other immune cells extend transmurally.
    - Macrophages also live in the submucosa and around the myenteric plexus (muscularis macrophages) [21]. Do
      not constrain macrophages by layer.
  - Confidence: high (moderate).
  - Evidence: [11, 17, 21]

- **SI-10. Lymphatic EC include the central lacteal in every villus core, so lymphatics ARE expected in the
  small-intestinal mucosa, unlike the colon. Smooth muscle strands also run up the villus core along the lacteal.**
  - Rule:
    - `kind=region, types=[Lymphatic EC], regions=[villi, submucosa], mode=enrich, strength=weak`
    - Do not port COL-13 to the small intestine.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Lymphangiectasia and lacteal dilatation.
    - Villus-core smooth muscle means that a Smooth muscle cells exclusion from the mucosa must stay weak or be
      omitted.
  - Confidence: high (weak, mainly as a "do not exclude" note).
  - Evidence: [12]; villus smooth muscle is textbook histology, not verified online.

- **SI-11. Smooth muscle cells mainly form the muscularis propria (inner circular and outer longitudinal) and the
  muscularis mucosae.**
  - Rule: `kind=region, types=[Smooth muscle cells], regions=[muscularis propria], mode=enrich, strength=moderate`
  - Without drawn regions: none (see SI-10).
  - Exceptions and caveats:
    - Vascular smooth muscle in every layer.
    - Pericytes and myofibroblasts share ACTA2 and TAGLN.
    - Subepithelial telocytes and myofibroblasts are fibroblast-lineage cells under the epithelium [20]. They are
      not smooth muscle.
  - Confidence: high (moderate).
  - Evidence: [20]; layer anatomy is textbook histology, not verified online.

- **SI-12. Neurons are only in the submucosal and myenteric ganglia. Glia occur in all layers, including the
  mucosa.**
  - Rule:
    - `kind=region, types=[Neurons], regions=[submucosa, muscularis propria, myenteric plexus], mode=only, strength=moderate`
    - `kind=region, types=[Neurons, Glia], regions=[myenteric plexus], mode=enrich, strength=strong`
  - Without drawn regions: `kind=depth, types=[Neurons], max_depth_um=40, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Glia must not be excluded from the mucosa [13].
    - Small mammals have a single submucosal ganglion layer [13].
    - Hirschsprung disease and other aganglionosis, and ganglioneuroma.
  - Confidence: high (moderate).
  - Evidence: [13, 14]

- **SI-13. ICC are confined to the muscularis propria. In the small intestine this includes the myenteric plexus
  level and the deep muscular plexus inside the circular muscle.**
  - Rule: `kind=region, types=[ICC], regions=[muscularis propria, myenteric plexus], mode=only, strength=moderate`
  - Without drawn regions: `kind=depth, types=[ICC], max_depth_um=40, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Mast cells and ILCs share KIT.
    - GIST.
  - Confidence: high (moderate).
  - Evidence: [15]; ICC sub-locations are textbook physiology, not verified online beyond the abstract of [15].

**References**
- [1] Barker N et al., Identification of stem cells in small intestine and colon by marker gene Lgr5, Nature, 2007. doi:10.1038/nature06196, PMID 17934449.
- [2] Clevers HC, Bevins CL, Paneth cells: maestros of the small intestinal crypts, Annu Rev Physiol, 2013. doi:10.1146/annurev-physiol-030212-183744, PMID 23398152.
- [3] Haber AL et al., A single-cell survey of the small intestinal epithelium, Nature, 2017. doi:10.1038/nature24489, PMID 29144463.
- [4] Moor AE et al., Spatial reconstruction of single enterocytes uncovers broad zonation along the intestinal villus axis, Cell, 2018. doi:10.1016/j.cell.2018.08.063, PMID 30270040.
- [5] Beumer J et al., Enteroendocrine cells switch hormone expression along the crypt-to-villus BMP signalling gradient, Nat Cell Biol, 2018. doi:10.1038/s41556-018-0143-y, PMID 30038251.
- [6] Mabbott NA et al., Microfold (M) cells: important immunosurveillance posts in the intestinal epithelium, Mucosal Immunol, 2013. doi:10.1038/mi.2013.30, PMID 23695511.
- [7] Jung C, Hugot JP, Barreau F, Peyer's patches: the immune sensors of the intestine, Int J Inflam, 2010. doi:10.4061/2010/823710, PMID 21188221.
- [8] Buettner M, Lochner M, Development and function of secondary and tertiary lymphoid organs in the small intestine and the colon, Front Immunol, 2016. doi:10.3389/fimmu.2016.00342, PMID 27656182.
- [9] Cheroutre H et al., The light and dark sides of intestinal intraepithelial lymphocytes, Nat Rev Immunol, 2011. doi:10.1038/nri3007, PMID 21681197.
- [10] Hayat M et al., Quantitation of intraepithelial lymphocytes in human duodenum: what is normal?, J Clin Pathol, 2002. PMID 11986350.
- [11] Mowat AM, Agace WW, Regional specialization within the intestinal immune system, Nat Rev Immunol, 2014. doi:10.1038/nri3738, PMID 25234148.
- [12] Bernier-Latmani J, Petrova TV, Intestinal lymphatic vasculature: structure, mechanisms and functions, Nat Rev Gastroenterol Hepatol, 2017. doi:10.1038/nrgastro.2017.79, PMID 28655884.
- [13] Furness JB, The enteric nervous system and neurogastroenterology, Nat Rev Gastroenterol Hepatol, 2012. doi:10.1038/nrgastro.2012.32, PMID 22392290.
- [14] Drokhlyansky E et al., The human and mouse enteric nervous system at single-cell resolution, Cell, 2020. doi:10.1016/j.cell.2020.08.003, PMID 32888429.
- [15] Sanders KM et al., Interstitial cells: regulators of smooth muscle function, Physiol Rev, 2014. doi:10.1152/physrev.00037.2013, PMID 24987007.
- [16] Elmentaite R et al., Cells of the human intestinal tract mapped across space and time, Nature, 2021. doi:10.1038/s41586-021-03852-1, PMID 34497389.
- [17] Hickey JW et al., Organization of the human intestine at single-cell resolution, Nature, 2023. doi:10.1038/s41586-023-05915-x, PMID 37468586.
- [18] Burclaff J et al., A proximal-to-distal survey of healthy adult human small intestine and colon epithelium by single-cell transcriptomics, Cell Mol Gastroenterol Hepatol, 2022. doi:10.1016/j.jcmgh.2022.02.007, PMID 35176508.
- [19] Busslinger GA et al., Human gastrointestinal epithelia of the esophagus, stomach, and duodenum resolved at single-cell resolution, Cell Rep, 2021. doi:10.1016/j.celrep.2021.108819, PMID 33691112.
- [20] Shoshkes-Carmel M et al., Subepithelial telocytes are an important source of Wnts that supports intestinal crypts, Nature, 2018. doi:10.1038/s41586-018-0084-4, PMID 29720649.
- [21] De Schepper S et al., Self-maintaining gut macrophages are essential for intestinal homeostasis, Cell, 2018. doi:10.1016/j.cell.2018.07.048, PMID 30173915.

---

## Stomach

**Compartments to draw**
- **surface/pit**: the luminal surface epithelium and the gastric pits (foveolae). These are tall columnar cells
  with pale apical mucin, and the pits make up the upper part of the mucosa.
- **isthmus/neck**: the narrow zone just below the pits. Small mucous neck cells mix with parietal cells here.
- **gland base**: the lower part of the oxyntic glands. Chief cells there have basophilic (blue-purple) granular
  cytoplasm, mixed with parietal cells, which are large, round and pink-eosinophilic.
- **oxyntic mucosa** (fundus/body/corpus): long, straight, tightly packed glands rich in parietal cells. The
  mucosa is thick.
- **antral mucosa** (antrum/pylorus): longer pits and shorter, coiled, pale mucous glands, with few or no
  parietal cells. The "transitional" zone between the two has mixed glands.
- **muscularis mucosae**: a thin smooth muscle band under the gland bases.
- **submucosa**, **muscularis propria** (thick, three layers: inner oblique, middle circular, outer longitudinal),
  **serosa**.
- **lymphoid follicle**: basophilic nodules in the deep mucosa. In the stomach they mostly mean gastritis (see
  STO-09).
- **intestinal metaplasia** (if present): glands with goblet cells and a brush border in gastric mucosa.
- Mouse only: **forestomach**, a keratinised squamous epithelium sharply demarcated from the glandular stomach.
  Squamous cell types are not in the built-in list.

A note on depth in the stomach: the oxyntic and antral glands fill the whole mucosal thickness, so the entire
mucosa (pits, glands and the lamina propria between them) has depth ~0. Depth therefore separates mucosa from
submucosa and muscularis well, but it carries no information on pit-versus-base zonation.

**Priors**

- **STO-01. Gastric epithelial cells (all subtypes) lie in the mucosa only.**
  - Rule: `kind=region, types=[Epithelial], regions=[oxyntic mucosa, antral mucosa, surface/pit, isthmus/neck, gland base], mode=only, strength=strong`
  - Without drawn regions: `kind=depth, types=[Epithelial], max_depth_um=30, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Diffuse-type (signet-ring) gastric carcinoma infiltrates the wall as single cells, far from any epithelial
      core. Depth and region rules would suppress exactly the tumour cells, so do not use either rule in cancer
      samples.
    - Also gastritis cystica profunda and pancreatic heterotopia (acinar and ductal cells in the submucosa or
      muscularis).
  - Confidence: high in normal tissue (strong). Not for tumours.
  - Evidence: [1, 3]

- **STO-02. Pit mucous cells (MUC5AC, TFF1, GKN1/2) line the surface and the pits. They are not in the gland
  bases.**
  - Rule:
    - `kind=region, types=[Pit mucous cells], regions=[surface/pit], mode=enrich, strength=strong`
    - `kind=region, types=[Pit mucous cells], regions=[gland base], mode=exclude, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Foveolar hyperplasia (reactive gastropathy, near ulcers) lengthens the pit zone.
    - Antral pits are deeper than oxyntic pits, roughly half the mucosa. Draw the surface/pit region per area.
  - Confidence: high (strong for enrich, moderate for exclude).
  - Evidence: [1, 3, 4]

- **STO-03. Neck mucous cells (MUC6, TFF2) sit in the neck of oxyntic glands. In the antrum, MUC6+ mucous cells
  make up most of the gland below the pits.**
  - Rule:
    - Oxyntic mucosa: `kind=region, types=[Neck mucous cells], regions=[isthmus/neck], mode=enrich, strength=moderate`
    - Antrum: `kind=region, types=[Neck mucous cells], regions=[antral mucosa], mode=enrich, strength=moderate`.
      The antral gland mucous cells have no separate type in the list.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - SPEM (spasmolytic polypeptide-expressing metaplasia: TFF2+, MUC6+ cells at the base of corpus glands after
      parietal cell loss) puts this type at the gland base. This is a real and important finding [6], so never
      use `exclude` in the gland base.
    - Brunner's-like glands at the pylorus.
  - Confidence: medium-high (moderate).
  - Evidence: [1, 2, 6]

- **STO-04. Chief cells (PGC, PGA3, LIPF) sit at the base of the oxyntic glands (fundus/body). They are rare in the
  antrum.**
  - Rule:
    - `kind=region, types=[Chief cells], regions=[gland base], mode=enrich, strength=strong`
    - `kind=region, types=[Chief cells], regions=[surface/pit], mode=exclude, strength=strong`
    - `kind=region, types=[Chief cells], regions=[antral mucosa], mode=exclude, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - In humans, individual antral glands contain MIST1+ chief cells [2], so keep the antrum exclusion weak.
    - PGC is also made by antral and neck mucous cells.
    - In atrophic gastritis and SPEM, chief cells are lost or transdifferentiate [6].
  - Confidence: high for base versus surface (strong). Low for the antrum (weak).
  - Evidence: [1, 2, 6]

- **STO-05. Parietal cells (ATP4A/ATP4B) occupy the upper and middle parts of the oxyntic glands (isthmus, neck
  and upper base) of the fundus/body. They are not in the surface or pit epithelium. Human antral glands are
  frequently mixed.**
  - Rule:
    - `kind=region, types=[Parietal cells], regions=[surface/pit], mode=exclude, strength=strong`
    - `kind=region, types=[Parietal cells], regions=[oxyntic mucosa], mode=enrich, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - In humans, parietal cells occur in over 50% of antral glands [2]. Do NOT exclude Parietal cells from the
      antrum in human samples.
    - Rodents have a sharper body/antrum boundary [2]. In mouse, a weak antrum exclusion is acceptable.
    - Autoimmune and atrophic gastritis lose parietal cells.
    - PPI therapy causes parietal cell protrusions but does not move them.
  - Confidence: high (strong for the pit exclusion).
  - Evidence: [1, 2, 3]

- **STO-06. Proliferating cells (MKI67+) of the gastric epithelium sit in the isthmus, which is the stem/progenitor
  zone of corpus glands.**
  - Rule:
    - `kind=region, types=[Proliferating cells], regions=[isthmus/neck], mode=enrich, strength=moderate`
    - `kind=neighborhood, types=[Proliferating cells], neighbor_types=[Epithelial], min_frac=0.3, mode=only, strength=moderate`
  - Without drawn regions: the neighbourhood rule above. Proliferating lymphocytes in follicles must not become
    epithelial Proliferating cells.
  - Exceptions and caveats:
    - In the antrum the proliferative zone lies lower, in the lower pit to the gland base.
    - Human and mouse corpus glands also have a slower basal (chief-cell-associated) stem population.
    - Injury, SPEM and dysplasia expand the proliferative zone.
  - Confidence: medium-high (moderate).
  - Evidence: [1, 3, 5]

- **STO-07. Enteroendocrine cells are scattered single cells in the glands, and their subtypes are regional. Gastrin
  (G) cells define the antrum, ghrelin cells are in the body, and somatostatin cells are concentrated
  proximally.**
  - Rule:
    - All Enteroendocrine cells: `kind=neighborhood, types=[Enteroendocrine cells], neighbor_types=[Epithelial], min_frac=0.3, mode=only, strength=moderate`
    - If the user splits the type, for example into "G cells" (not in the built-in list):
      `kind=region, types=[G cells], regions=[antral mucosa], mode=only, strength=moderate`
  - Without drawn regions: the neighbourhood rule above.
  - Exceptions and caveats:
    - ECL-cell hyperplasia and neuroendocrine tumours occur in atrophic gastritis.
    - Antral and duodenal enterochromaffin cells differ in expression [4].
    - Human-specific LH+ EEC [4].
  - Confidence: high (moderate).
  - Evidence: [2, 4]

- **STO-08. Intestinal metaplasia (CDX2+, MUC2+, TFF3+) is absent from normal gastric mucosa. When present, it
  forms whole glands or gland patches, not single scattered cells.**
  - Rule: `kind=region, types=[Intestinal metaplasia], regions=[intestinal metaplasia], mode=enrich, strength=moderate`.
    Use only when the user drew such areas. Optionally, `kind=neighborhood, types=[Intestinal metaplasia], neighbor_types=[Intestinal metaplasia], min_frac=0.2, mode=only, strength=weak`
  - Without drawn regions: the weak neighbourhood rule above.
  - Exceptions and caveats:
    - Intestinal metaplasia is a real, common finding in H. pylori and atrophic gastritis and in the context of
      cancer [3, 6]. Never use a global exclude.
    - TFF3 is shared with Lymphatic EC, and MUC2 with goblet-like spill-over.
  - Confidence: medium (moderate as enrich).
  - Evidence: [3, 6]

- **STO-09. Lymphoid follicles are essentially absent from normal gastric mucosa. In biopsies, lymphoid aggregates
  are virtually always associated with H. pylori-related chronic active gastritis, except near ulcers. B cells are
  therefore sparse in normal stomach and concentrated in follicles when follicles exist.**
  - Rule: `kind=region, types=[B cells], regions=[lymphoid follicle], mode=enrich, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[B cells], neighbor_types=[B and plasma cells, T and NK cells], min_frac=0.6, mode=enrich, strength=moderate`
  - Exceptions and caveats:
    - Do not exclude B cells from the mucosa: the follicles are the finding.
    - MALT lymphoma.
    - Plasma cells are normal in the superficial lamina propria between the pits.
  - Confidence: high (strong for enrich within drawn follicles).
  - Evidence: [7]

- **STO-10. Neurons are in the myenteric plexus of the muscularis propria. Unlike the intestine, the stomach has no
  ganglionated submucosal plexus.**
  - Rule:
    - `kind=region, types=[Neurons], regions=[muscularis propria], mode=only, strength=moderate`
    - `kind=region, types=[Neurons], regions=[submucosa], mode=exclude, strength=weak`
  - Without drawn regions: `kind=depth, types=[Neurons], max_depth_um=40, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Scattered submucosal neurons may still occur, so keep the submucosa exclusion weak.
    - Glia occur along nerve fibres in all layers; do not exclude them.
    - In the pyloric region, the duodenal submucosal plexus begins.
  - Confidence: medium-high (moderate).
  - Evidence: [8]

- **STO-11. ICC lie in the muscularis propria. The gastric corpus and antrum pacemaker ICC are around the
  myenteric plexus and within the muscle layers.**
  - Rule: `kind=region, types=[ICC], regions=[muscularis propria], mode=only, strength=moderate`
  - Without drawn regions: `kind=depth, types=[ICC], max_depth_um=40, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Mast cells and ILCs (KIT+) in the mucosa and submucosa.
    - The stomach is the most common site of GIST, an ICC-lineage tumour that forms a muscularis-based mass.
  - Confidence: high (moderate).
  - Evidence: [9]; ICC sub-locations are textbook physiology, not verified online beyond the abstract of [9].

- **STO-12. Smooth muscle cells mainly form the three-layered muscularis propria and the muscularis mucosae.
  Strands from the muscularis mucosae also ascend between the gastric glands.**
  - Rule: `kind=region, types=[Smooth muscle cells], regions=[muscularis propria, muscularis mucosae], mode=enrich, strength=moderate`.
    Do not exclude smooth muscle from the mucosa.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Leiomyoma.
    - Pericytes and myofibroblasts share markers.
  - Confidence: high (moderate).
  - Evidence: textbook histology, not verified online.

**References**
- [1] Mills JC, Shivdasani RA, Gastric epithelial stem cells, Gastroenterology, 2011. doi:10.1053/j.gastro.2010.12.001, PMID 21144849.
- [2] Choi E et al., Cell lineage distribution atlas of the human stomach reveals heterogeneous gland populations in the gastric antrum, Gut, 2014. doi:10.1136/gutjnl-2013-305964, PMID 24488499.
- [3] Tsubosaka A et al., Stomach encyclopedia: Combined single-cell and spatial transcriptomics reveal cell diversity and homeostatic regulation of human stomach, Cell Rep, 2023. doi:10.1016/j.celrep.2023.113236, PMID 37819756.
- [4] Busslinger GA et al., Human gastrointestinal epithelia of the esophagus, stomach, and duodenum resolved at single-cell resolution, Cell Rep, 2021. doi:10.1016/j.celrep.2021.108819, PMID 33691112.
- [5] Dong J et al., Spatially resolved expression landscape and gene-regulatory network of human gastric corpus epithelium, Protein Cell, 2023. doi:10.1093/procel/pwac059, PMID 37402315.
- [6] Goldenring JR, Mills JC, Cellular plasticity, reprogramming, and regeneration: metaplasia in the stomach and beyond, Gastroenterology, 2022. PMID 34728185.
- [7] Genta RM, Hamner HW, The significance of lymphoid follicles in the interpretation of gastric biopsy specimens, Arch Pathol Lab Med, 1994. PMID 7517659.
- [8] Furness JB, The enteric nervous system and neurogastroenterology, Nat Rev Gastroenterol Hepatol, 2012. doi:10.1038/nrgastro.2012.32, PMID 22392290.
- [9] Sanders KM et al., Interstitial cells: regulators of smooth muscle function, Physiol Rev, 2014. doi:10.1152/physrev.00037.2013, PMID 24987007.

**Verification note**
- Every reference above was checked against PubMed or publisher records during this task. For most, only the
  abstract or bibliographic record was read.
- The full text was read only for one point: Furness 2012 [COL 12 / STO 8] gives the plexus locations and the
  absence of a ganglionated submucosal plexus in the stomach.
- Statements marked "textbook histology, not verified online" have no online source behind them.
- No mouse per-layer thickness was verified, so the mouse-safe depth values (30-40 um) are estimates to calibrate
  on the sample.

---

## Liver, pancreas, kidney: shared notes

Cell-type names are those of the built-in lists `liver.txt`, `pancreas.txt`, `kidney.txt` (with `@immune`, and
`@stroma` for pancreas and kidney). Note that `liver.txt` does not include `@stroma`: it defines its own
Endothelial (LSEC, Vascular EC, Lymphatic EC), Hepatic stellate cells, Fibroblasts and Smooth muscle cells, and has
no Pericytes. In `kidney.txt` the endothelium comes from `@stroma` (Endothelial > Blood vessel EC, Lymphatic EC);
there is no separate glomerular endothelial type.

**The depth rule is not meaningful in these three organs.** `depth_um` is the distance to cells dense in EPCAM+KRT8
(`stagent/importer.py`, `DEPTH_MARKERS`). Hepatocytes, acinar cells and every renal tubule segment express KRT8, so
the "core" covers almost the whole parenchyma and depth is about 0 nearly everywhere; it would at most pick out
large vessels, capsule and fat. All priors below therefore use drawn regions, with neighborhood rules as the
fallback. Neighborhood rules use the expression-only typing, so they only help when the neighbor types are
themselves called well, and they cannot fix a whole cluster that is wrong.

---

## Liver

**Compartments to draw**
- `portal tract`: triangular or round connective-tissue island with a bile duct (cuboidal epithelium), a hepatic
  artery branch (thick wall) and a portal vein branch (thin wall, wide lumen); draw it up to the limiting plate of
  hepatocytes.
- `central vein`: round thin-walled vein alone in the parenchyma, sinusoids radiating into it, no bile duct; draw
  the vein wall and lumen only.
- `periportal zone` (zone 1): the 2-4 hepatocyte plates around each portal tract.
- `pericentral zone` (zone 3): the 2-4 hepatocyte plates around each central vein.
- `parenchyma` (or `lobule`): hepatocyte plates and sinusoids, i.e. everything except portal tracts, central veins,
  large vessels and capsule (zone 2 is the part not covered by the two zones above).
- `capsule`: thin fibrous Glisson's capsule at the organ surface, with mesothelium.
- `hilum / large vessels`: large portal/hepatic vein branches and large ducts with thick fibrous sheaths.

Drawing every portal tract and central vein is laborious; draw them in a representative field and rely on the
neighborhood fallbacks elsewhere.

**Priors**

- **LIV-01. Hepatocytes are not inside portal tract stroma or inside central vein walls.**
  - Rule: `kind=region, types=[Hepatocytes], regions=[portal tract, central vein], mode=exclude, strength=moderate`
  - Without drawn regions: none (hepatocytes are the majority type everywhere else).
  - Exceptions and caveats: periportal hepatocytes of the limiting plate abut the tract, so a loosely drawn outline
    will include real hepatocytes; ALB/APO transcripts are abundant and spill into adjacent stromal and immune
    cells after segmentation, which is exactly what this rule corrects. Breaks in interface hepatitis, ductular
    reaction, cirrhosis (nodules surrounded by septa) and in tumors (HCC, metastases).
  - Confidence: high; moderate (strong only for tightly drawn outlines).
  - Evidence: [2], [12] (textbook lobule anatomy).

- **LIV-02. Periportal hepatocytes (CPS1, HAL, SDS, ASS1, PCK1) sit around portal tracts, not around central veins.**
  - Rule: `kind=region, types=[Periportal hepatocytes], regions=[periportal zone], mode=enrich, strength=moderate`;
    and `kind=region, types=[Periportal hepatocytes], regions=[pericentral zone], mode=exclude, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Periportal hepatocytes], neighbor_types=[Cholangiocytes,
    Fibroblasts, Smooth muscle cells, Vascular EC], min_frac=0.1, mode=enrich, strength=weak` (cells next to portal
    structures).
  - Exceptions and caveats: zonation is a continuous gradient and about half of liver genes are zonated, many
    peaking mid-lobule [1]; the zone outlines are approximate. Zonation is blunted or lost in cirrhosis, in
    regenerating nodules, in HCC and after pericentral injury (e.g. acetaminophen); feeding state and circadian
    phase shift the gradient. Zonation genes are largely conserved between mouse and human [3], but mouse lobules
    are smaller (fewer plates per zone).
  - Confidence: high; moderate.
  - Evidence: [1], [2], [3], [4], [15]

- **LIV-03. Pericentral hepatocytes (GLUL, CYP2E1, CYP1A2, OAT, RHBG) sit around central veins; GLUL+ hepatocytes
  form only a thin ring of 1-3 cell layers around the central vein.**
  - Rule: `kind=region, types=[Pericentral hepatocytes], regions=[pericentral zone], mode=enrich,
    strength=moderate`; and `kind=region, types=[Pericentral hepatocytes], regions=[periportal zone, portal tract],
    mode=exclude, strength=moderate`
  - Without drawn regions: none (central-vein endothelium is typed as Vascular EC, like portal vessels, so there is
    no specific neighbor type).
  - Exceptions and caveats: CYP2E1 extends further into zone 2 than GLUL; GLUL is also expressed by some
    non-hepatocytes. Pericentral loss in centrilobular necrosis (acetaminophen, ischemia) and in congestion;
    beta-catenin-mutated HCC and adenomas express GLUL diffusely, so the rule would wrongly suppress a real
    pericentral-like tumor program.
  - Confidence: high; moderate.
  - Evidence: [1], [2], [3], [4]

- **LIV-04. Cholangiocytes form bile ducts and ductules in portal tracts (canals of Hering at the portal interface);
  they are not in the pericentral parenchyma.**
  - Rule: `kind=region, types=[Cholangiocytes], regions=[portal tract, periportal zone, hilum / large vessels],
    mode=only, strength=moderate`; and `kind=region, types=[Cholangiocytes], regions=[pericentral zone, central
    vein], mode=exclude, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Cholangiocytes], neighbor_types=[Cholangiocytes, Fibroblasts,
    Smooth muscle cells, Vascular EC, Lymphatic EC], min_frac=0.3, mode=enrich, strength=moderate`
  - Exceptions and caveats: ductular reaction (cholestasis, PSC/PBC, advanced fibrosis, regeneration) puts KRT19+/
    KRT7+ cells deep in septa and parenchyma; hepatocytes can acquire biliary markers; cholangiocarcinoma and
    metastatic adenocarcinoma are EPCAM/KRT19+ anywhere. Use `exclude ... strong` only in healthy tissue.
  - Confidence: high; moderate (only) / strong (pericentral exclude, normal liver).
  - Evidence: [4], [12]

- **LIV-05. Fibroblasts (portal fibroblasts) and Smooth muscle cells (hepatic artery and portal vein walls) are in
  portal tracts, around central veins, in the capsule and in hilar vessels, not between hepatocyte plates.**
  - Rule: `kind=region, types=[Fibroblasts, Smooth muscle cells], regions=[portal tract, central vein, capsule,
    hilum / large vessels], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Fibroblasts, Smooth muscle cells], neighbor_types=[Hepatocytes],
    min_frac=0.8, mode=exclude, strength=weak`
  - Exceptions and caveats: in fibrosis, activated stellate cells become COL1A1+ myofibroblasts in septa and
    pericentral scars and may be typed Fibroblasts; the rule would then hide real fibrogenesis, so do not use it
    in fibrotic or cirrhotic liver. Central-vein adventitia is thin (few cells). ACTA2 is also on activated
    stellate cells.
  - Confidence: medium; moderate.
  - Evidence: [10], [12]

- **LIV-06. Hepatic stellate cells sit in the space of Disse, between sinusoids and hepatocytes, throughout the
  lobule (with portal-side and central-side subsets), not in portal tract stroma.**
  - Rule: `kind=region, types=[Hepatic stellate cells], regions=[parenchyma, periportal zone, pericentral zone],
    mode=enrich, strength=weak`; and `kind=region, types=[Hepatic stellate cells], regions=[portal tract,
    capsule], mode=exclude, strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[Hepatic stellate cells], neighbor_types=[Hepatocytes, LSEC],
    min_frac=0.5, mode=enrich, strength=weak`
  - Exceptions and caveats: stellate cells are thin and their transcripts mix with hepatocytes and LSEC after
    segmentation; RGS5/PDGFRB also mark vascular mural cells in portal tracts. In injury they migrate and become
    myofibroblasts; central-vein-associated stellate cells drive centrilobular fibrosis in mouse [10].
  - Confidence: medium; weak.
  - Evidence: [10], [11], [8]

- **LIV-07. LSEC (CLEC4G, STAB2, CRHBP, FCN2, OIT3) line the sinusoids inside the lobule; the endothelium of portal
  vessels, central veins and large vessels is Vascular EC.**
  - Rule: `kind=region, types=[LSEC], regions=[portal tract, hilum / large vessels, capsule], mode=exclude,
    strength=moderate`; and `kind=region, types=[Vascular EC], regions=[portal tract, central vein, hilum / large
    vessels], mode=enrich, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[LSEC], neighbor_types=[Hepatocytes, Hepatic stellate cells,
    Kupffer cells], min_frac=0.5, mode=enrich, strength=weak`
  - Exceptions and caveats: LSEC are zonated; periportal sinusoids are less LSEC-like and pericentral LSEC grade
    into central-vein endothelium [4], [6], so cells at the tract border and the central vein are mixed. Sinusoids
    capillarize (lose LSEC markers) in fibrosis, cirrhosis and HCC, where Vascular EC-like cells are real in the
    parenchyma; do not exclude Vascular EC from the parenchyma in diseased liver.
  - Confidence: high; moderate.
  - Evidence: [4], [5], [6]

- **LIV-08. Lymphatic EC are confined to portal tracts (plus capsule and around large veins); there are no
  lymphatics in the sinusoids.**
  - Rule: `kind=region, types=[Lymphatic EC], regions=[portal tract, capsule, hilum / large vessels], mode=only,
    strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Lymphatic EC], neighbor_types=[Hepatocytes, LSEC],
    min_frac=0.7, mode=exclude, strength=moderate`
  - Exceptions and caveats: important for marker artifacts: LYVE1 is expressed by normal LSEC and PROX1 by
    hepatocytes, so parenchymal "Lymphatic EC" calls are usually false [13], [14]; the rule removes them. PDPN is
    the more specific liver lymphatic marker if on the panel. Lymphatics expand in cirrhosis and at tumor margins
    [14].
  - Confidence: high; strong.
  - Evidence: [13], [14]

- **LIV-09. Kupffer cells (CD5L, VSIG4, MARCO, TIMD4, CLEC4F) reside in the sinusoids of the lobule, enriched
  periportally; they are not part of the portal tract stroma or capsule.**
  - Rule: `kind=region, types=[Kupffer cells], regions=[portal tract, capsule], mode=exclude, strength=moderate`;
    and `kind=region, types=[Kupffer cells], regions=[periportal zone], mode=enrich, strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[Kupffer cells], neighbor_types=[Hepatocytes, LSEC,
    Hepatic stellate cells], min_frac=0.5, mode=enrich, strength=weak`
  - Exceptions and caveats: the Kupffer-cell niche is formed by LSEC, stellate cells and hepatocytes [8]; other
    macrophage types (capsule macrophages, bile-duct/lipid-associated macrophages, monocyte-derived macrophages)
    live in portal tracts, capsule and steatotic areas [7] and should be typed Macrophages, not Kupffer cells.
    Periportal enrichment is shown in mouse and depends on commensal signals [9]. In NASH/steatosis and after
    Kupffer-cell loss, monocyte-derived cells fill the sinusoids.
  - Confidence: high (sinusoidal location); medium (periportal bias); moderate / weak.
  - Evidence: [5], [7], [8], [9]

- **LIV-10. In normal liver, B cells, plasma cells and most dendritic cells sit in portal tracts; lymphocyte and
  myeloid density is higher periportally.**
  - Rule: `kind=region, types=[B cells, Plasma cells, Dendritic cells], regions=[portal tract], mode=enrich,
    strength=moderate`; `kind=region, types=[T and NK cells], regions=[portal tract, periportal zone], mode=enrich,
    strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[B cells, Plasma cells], neighbor_types=[Cholangiocytes,
    Fibroblasts, Vascular EC, Lymphatic EC], min_frac=0.2, mode=enrich, strength=weak`
  - Exceptions and caveats: NK cells and many T cells (liver-resident NK, MAIT) patrol the sinusoids, so do not
    exclude T and NK cells from the parenchyma. Portal lymphoid aggregates and follicles are common in chronic
    hepatitis C, PBC, PSC and autoimmune hepatitis; lobular inflammation occurs in hepatitis and NASH. Use
    `enrich` only, never `only`.
  - Confidence: medium; moderate / weak.
  - Evidence: [9], [5] (portal lymphoid aggregates in disease: textbook histology, not verified online)

- **LIV-11. Mast cells are found in portal tracts and capsule, not in the sinusoids.**
  - Rule: `kind=region, types=[Mast cells], regions=[portal tract, capsule, hilum / large vessels], mode=enrich,
    strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats: numbers rise in fibrosis and cholestatic disease; very rare in normal parenchyma.
  - Confidence: low; weak.
  - Evidence: (textbook histology, not verified online)

**References**
[1] Halpern KB et al., Single-cell spatial reconstruction reveals global division of labour in the mammalian liver,
Nature, 2017. DOI 10.1038/nature21065. PMID 28166538.
[2] Ben-Moshe S, Itzkovitz S, Spatial heterogeneity in the mammalian liver, Nature Reviews Gastroenterology &
Hepatology, 2019. DOI 10.1038/s41575-019-0134-x.
[3] Brosch M et al., Epigenomic map of human liver reveals principles of zonated morphogenic and metabolic control,
Nature Communications, 2018. DOI 10.1038/s41467-018-06611-5.
[4] Aizarani N et al., A human liver cell atlas reveals heterogeneity and epithelial progenitors, Nature, 2019.
DOI 10.1038/s41586-019-1373-2.
[5] MacParland SA et al., Single cell RNA sequencing of human liver reveals distinct intrahepatic macrophage
populations, Nature Communications, 2018. DOI 10.1038/s41467-018-06318-7.
[6] Halpern KB et al., Paired-cell sequencing enables spatial gene expression mapping of liver endothelial cells,
Nature Biotechnology, 2018. DOI 10.1038/nbt.4231.
[7] Guilliams M et al., Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage
niches, Cell, 2022. DOI 10.1016/j.cell.2021.12.018.
[8] Bonnardel J et al., Stellate cells, hepatocytes, and endothelial cells imprint the Kupffer cell identity on
monocytes colonizing the liver macrophage niche, Immunity, 2019. DOI 10.1016/j.immuni.2019.08.017. PMID 31561945.
[9] Gola A et al., Commensal-driven immune zonation of the liver promotes host defence, Nature, 2021.
DOI 10.1038/s41586-020-2977-2. PMID 33239787.
[10] Dobie R et al., Single-cell transcriptomics uncovers zonation of function in the mesenchyme during liver
fibrosis, Cell Reports, 2019. DOI 10.1016/j.celrep.2019.10.024.
[11] Friedman SL, Hepatic stellate cells: protean, multifunctional, and enigmatic cells of the liver, Physiological
Reviews, 2008 (88:125-172). PMC2888531.
[12] Andrews TS et al., Single-cell, single-nucleus, and spatial RNA sequencing of the human liver identifies
cholangiocyte and mesenchymal heterogeneity, Hepatology Communications, 2022. DOI 10.1002/hep4.1854.
[13] Tanaka M, Iwakiri Y, The hepatic lymphatic vascular system: structure, function, markers, and
lymphangiogenesis, Cellular and Molecular Gastroenterology and Hepatology, 2016. DOI 10.1016/j.jcmgh.2016.09.002.
[14] Mouta Carreira C et al., LYVE-1 is not restricted to the lymph vessels: expression in normal liver blood
sinusoids and down-regulation in human liver cancer and cirrhosis, Cancer Research, 2001 (61:8079-8084).
[15] Yu S et al., Spatial transcriptome profiling of normal human liver, Scientific Data, 2022.
DOI 10.1038/s41597-022-01676-w.

---

## Pancreas

**Compartments to draw**
- `islet`: round or oval pale-staining cluster (about 50-300 um) of polygonal cells with pale cytoplasm and a dense
  capillary network, sitting among acini; draw each islet tightly.
- `acinar parenchyma`: the bulk lobules of basophilic (dark purple base, pink apical granules) acini.
- `duct`: interlobular and main ducts (cuboidal to columnar epithelium) with their periductal fibrous collar.
- `septum`: interlobular connective tissue with vessels, nerves and fat.
- `PP-rich region`: the posterior head / uncinate process (ventral pancreas), only if the section is known to come
  from there.
- `lymph node`: peripancreatic node if included.

Small islets are many; draw a representative set and rely on the neighborhood fallbacks for the rest.

**Priors**

- **PAN-01. Endocrine cells (Alpha, Beta, Delta, PP) are in islets; single extra-islet endocrine cells exist but are
  rare.**
  - Rule: `kind=region, types=[Endocrine], regions=[islet], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Endocrine], neighbor_types=[Endocrine], min_frac=0.3,
    mode=enrich, strength=moderate`
  - Exceptions and caveats: scattered endocrine cells and small clusters in duct walls and between acini are real
    (more in neonates [5]); neuroendocrine tumors (PanNETs) are endocrine outside islets. Keep `moderate` so single
    cells with a strong hormone signal (GCG!, INS!, SST!, PPY!) still win.
  - Confidence: high; moderate.
  - Evidence: [1], [2], [5]

- **PAN-02. Acinar cells are not inside islets.**
  - Rule: `kind=region, types=[Acinar cells], regions=[islet], mode=exclude, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Acinar cells], neighbor_types=[Endocrine], min_frac=0.6,
    mode=exclude, strength=moderate`
  - Exceptions and caveats: PRSS1/CPA1/CELA3A/CTRB transcripts are extremely abundant and spill into islet cells at
    the islet edge and through segmentation errors; this rule mainly fixes that. Acinar cells directly at the islet
    border are real (the "peri-insular" acini), so draw islets tightly.
  - Confidence: high; moderate.
  - Evidence: [1], [2] (textbook histology for islet-acinar boundary)

- **PAN-03. Acinar cells form the lobules; they are not in ducts, septa or lymph nodes.**
  - Rule: `kind=region, types=[Acinar cells], regions=[duct, septum, lymph node], mode=exclude, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats: in pancreatitis and around PDAC, acinar-to-ductal metaplasia gives cells with mixed
    acinar/ductal programs, and acinar atrophy with fibrosis leaves isolated acini in fibrous tissue [5].
  - Confidence: high; moderate.
  - Evidence: [5] (textbook histology)

- **PAN-04. Ductal cells line the duct tree, from interlobular and main ducts down to intercalated ducts and
  centroacinar cells inside acini; they are not inside islets.**
  - Rule: `kind=region, types=[Ductal cells], regions=[duct], mode=enrich, strength=moderate`; and
    `kind=region, types=[Ductal cells], regions=[islet], mode=exclude, strength=weak`
  - Without drawn regions: none (centroacinar and intercalated duct cells are surrounded by acinar cells, so a
    neighborhood rule would penalize real cells).
  - Exceptions and caveats: do not use `only` in duct regions: centroacinar cells sit in acinar parenchyma. Islets
    often lie next to ducts. Ductal expansion in chronic pancreatitis, PanIN and PDAC.
  - Confidence: high; moderate / weak.
  - Evidence: [4], [5]

- **PAN-05. In human islets, alpha, beta and delta cells are intermingled, with many cells along intra-islet vessels;
  there is no beta-cell core and non-beta mantle. Do not use core/mantle priors on human tissue.**
  - Rule: none for human (deliberately no rule).
  - Without drawn regions: none.
  - Exceptions and caveats: human islets also contain proportionally more alpha cells and fewer beta cells than mouse
    islets, and composition varies strongly between islets [1], [2]; a mouse-style rule would wrongly push
    central human alpha cells towards Beta cells. Non-human primate islets resemble human ones [2].
  - Confidence: high; (no rule).
  - Evidence: [1], [2], [3]

- **PAN-06. In mouse (and rat) islets, Beta cells form the core and Alpha and Delta cells sit in the mantle at the
  islet periphery.**
  - Rule (mouse only): `kind=neighborhood, types=[Alpha cells, Delta cells], neighbor_types=[Exocrine, Fibroblasts,
    Endothelial], min_frac=0.15, mode=enrich, strength=weak` (non-beta cells touch the exocrine border); and
    `kind=neighborhood, types=[Beta cells], neighbor_types=[Beta cells], min_frac=0.6, mode=enrich, strength=weak`
  - Without drawn regions: the rule above is already region-free.
  - Exceptions and caveats: species-specific: rodent only [3]. Even in mouse the mantle is incomplete in large
    islets and alpha cells can line intra-islet vessels; diabetes and obesity models disrupt the architecture
    (alpha cells move into the core) [3]. Weak strength only.
  - Confidence: medium (mouse); weak.
  - Evidence: [1], [2], [3]

- **PAN-07. PP cells are concentrated in islets of the posterior head / uncinate process (ventral pancreas), where
  alpha cells are scarce; elsewhere PP cells are rare.**
  - Rule: `kind=region, types=[PP cells], regions=[PP-rich region], mode=enrich, strength=moderate`; only if the
    section origin is known.
  - Without drawn regions: none (PP cells cluster with each other: `kind=neighborhood, types=[PP cells],
    neighbor_types=[PP cells], min_frac=0.2, mode=enrich, strength=weak` is possible but low value).
  - Exceptions and caveats: sections from body/tail will have very few PP cells; a body/tail sample with many PP
    calls suggests a marker problem (PPY is a single-gene type). Alpha-cell mass is lower in the PP-rich region,
    so do not expect the usual alpha/beta mix there.
  - Confidence: high (human); moderate.
  - Evidence: [6]

- **PAN-08. Islets are densely vascularized: endothelial cells are expected inside islets and must not be
  suppressed there.**
  - Rule: `kind=region, types=[Endothelial], regions=[islet], mode=enrich, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats: endothelial transcripts mix with endocrine cells in dense islets; use only enrich.
  - Confidence: medium; weak.
  - Evidence: [1] (cells aligned along islet vessels); (textbook histology for vascular density, not verified
    online)

- **PAN-09. Fibroblasts and Smooth muscle cells are in septa, around ducts and vessels; stellate cells sit
  periacinar but are scarce in normal tissue.**
  - Rule: `kind=region, types=[Smooth muscle cells], regions=[septum, duct], mode=only, strength=moderate`;
    `kind=region, types=[Fibroblasts], regions=[septum, duct], mode=enrich, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats: pancreatic stellate cells (not in the built-in list) are periacinar, around small ducts
    and vessels, about 4-7% of cells [7], [8], and will be typed Fibroblasts or Pericytes inside acinar lobules, so
    do not exclude Fibroblasts from acinar parenchyma. Activated stellate cells/myofibroblasts expand massively in
    chronic pancreatitis and PDAC desmoplasia [5], [7].
  - Confidence: medium; moderate.
  - Evidence: [5], [7], [8]

- **PAN-10. Lymphocytes are sparse in normal pancreas; B cells and lymphoid aggregates belong to septa and lymph
  nodes, not to islets or acini.**
  - Rule: `kind=region, types=[B cells], regions=[lymph node, septum], mode=enrich, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats: insulitis in type 1 diabetes puts CD8 T cells and macrophages in and around islets;
    pancreatitis and PDAC bring dense infiltrates and tertiary lymphoid structures. Resident islet macrophages are
    normal. Never use `only`.
  - Confidence: low; weak.
  - Evidence: [5] (macrophage-stellate interplay in pancreatitis); (textbook histology for the rest, not verified
    online)

**References**
[1] Cabrera O et al., The unique cytoarchitecture of human pancreatic islets has implications for islet cell
function, PNAS, 2006. DOI 10.1073/pnas.0510790103. PMID 16461897.
[2] Brissova M et al., Assessment of human pancreatic islet architecture and composition by laser scanning confocal
microscopy, Journal of Histochemistry & Cytochemistry, 2005. DOI 10.1369/jhc.5C6684.2005.
[3] Steiner DJ, Kim A, Miller K, Hara M, Pancreatic islet plasticity: interspecies comparison of islet architecture
and composition, Islets, 2010. DOI 10.4161/isl.2.3.11815. PMID 20657742.
[4] Baron M et al., A single-cell transcriptomic map of the human and mouse pancreas reveals inter- and intra-cell
population structure, Cell Systems, 2016. DOI 10.1016/j.cels.2016.08.011. PMID 27667365.
[5] Tosti L et al., Single-nucleus and in situ RNA-sequencing reveal cell topographies in the human pancreas,
Gastroenterology, 2021. DOI 10.1053/j.gastro.2020.11.010. PMID 33212097.
[6] Wang X et al., Quantitative analysis of pancreatic polypeptide cell distribution in the human pancreas, PLoS
ONE, 2013. DOI 10.1371/journal.pone.0055501.
[7] Apte MV, Pirola RC, Wilson JS, Pancreatic stellate cells: a starring role in normal and diseased pancreas,
Frontiers in Physiology, 2012. DOI 10.3389/fphys.2012.00344.
[8] Apte MV et al., Periacinar stellate shaped cells in rat pancreas: identification, isolation, and culture, Gut,
1998 (43:128).

---

## Kidney

**Compartments to draw**
- `cortex`: outer band with glomeruli and convoluted tubules (labyrinth), crossed by straight medullary rays.
- `outer stripe` (outer medulla, outer stripe): just below the cortex, no glomeruli, straight proximal tubules
  (pink, brush border) mixed with thick limbs and collecting ducts.
- `inner stripe` (outer medulla, inner stripe): no proximal tubules; thick ascending limbs, thin descending limbs,
  collecting ducts and vascular bundles.
- `inner medulla` (papilla): pale, tapering zone of thin limbs and large collecting ducts in abundant interstitium.
- `glomerulus`: round capillary tuft with Bowman's capsule; draw it including the capsule (or draw `glomerular tuft`
  and `Bowman's capsule` separately if resolution allows).
- `large vessels`: interlobar/arcuate arteries and veins at the corticomedullary junction, with their adventitia.
- `capsule`, `pelvis` (urothelium-lined cavity at the hilum).

Human needle biopsies usually contain only cortex (sometimes outer medulla); mouse sections usually contain all
zones. If only `cortex` and `medulla` are drawn, use `medulla` where the rules below say outer/inner stripe or inner
medulla, and drop the rules that need the finer split.

**Priors**

- **KID-01. Podocytes are only on the glomerular tuft.**
  - Rule: `kind=region, types=[Podocytes], regions=[glomerulus], mode=only, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Podocytes], neighbor_types=[Podocytes, Mesangial cells,
    Endothelial, Parietal epithelial cells], min_frac=0.5, mode=enrich, strength=moderate`
  - Exceptions and caveats: podocytes detach and appear in Bowman's space and in crescents; glomeruli cut at the
    edge of a section are small and easily missed when drawing. NPHS1/NPHS2 are specific, so false calls outside
    glomeruli are usually spill-over.
  - Confidence: high; strong.
  - Evidence: [1], [6], [7]

- **KID-02. Glomeruli (and therefore Podocytes, Mesangial cells and Parietal epithelial cells) are only in the
  cortex, including the juxtamedullary cortex; never in the outer or inner medulla.**
  - Rule: `kind=region, types=[Podocytes, Mesangial cells, Parietal epithelial cells], regions=[outer stripe,
    inner stripe, inner medulla], mode=exclude, strength=strong`
  - Without drawn regions: none (use KID-01).
  - Exceptions and caveats: the corticomedullary border is irregular; juxtamedullary glomeruli sit right at it.
  - Confidence: high; strong.
  - Evidence: [1], [2]

- **KID-03. Mesangial cells are inside the glomerular tuft (plus the extraglomerular mesangium at the vascular
  pole).**
  - Rule: `kind=region, types=[Mesangial cells], regions=[glomerulus], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Mesangial cells], neighbor_types=[Podocytes, Endothelial],
    min_frac=0.3, mode=enrich, strength=moderate`
  - Exceptions and caveats: extraglomerular mesangial (lacis) cells and renin cells (not in the built-in list) sit
    just outside the tuft at the vascular pole, so draw glomeruli generously there. The list's Mesangial markers
    overlap with interstitial cells (PDGFRB with pericytes/fibroblasts, GATA3 with principal and connecting tubule
    cells, PIEZO2 with renin cells); the rule fixes interstitial Pericytes/Fibroblasts called Mesangial.
    Mesangial expansion in diabetic nephropathy and IgA nephropathy stays within glomeruli.
  - Confidence: high; moderate.
  - Evidence: [9], [6]

- **KID-04. Parietal epithelial cells line Bowman's capsule; outside glomeruli they are not expected.**
  - Rule: `kind=region, types=[Parietal epithelial cells], regions=[glomerulus], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Parietal epithelial cells], neighbor_types=[Podocytes,
    Mesangial cells, Proximal tubule], min_frac=0.2, mode=enrich, strength=weak`
  - Exceptions and caveats: at the urinary pole PECs merge into proximal tubule. The PEC markers CLDN1 and VCAM1 are
    also expressed by injured/adaptive (failed-repair) proximal tubule cells [6], so in injured kidney many
    "Parietal epithelial cells" outside glomeruli are actually injured PT, and this rule helps; but PEC
    activation in crescents and FSGS is real and stays inside glomeruli.
  - Confidence: medium; moderate.
  - Evidence: [8], [6]

- **KID-05. Tubular epithelium is not on the glomerular tuft.**
  - Rule: `kind=region, types=[Proximal tubule, Thin limb, Thick ascending limb, Distal convoluted tubule,
    Connecting tubule, Principal cells, Intercalated cells], regions=[glomerulus], mode=exclude, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats: proximal tubule starts at the urinary pole and the macula densa (TAL segment, not in
    the built-in list) touches the vascular pole, so cells on the outline are mixed; draw the tuft tightly or keep
    moderate.
  - Confidence: high; moderate.
  - Evidence: [1], [2]

- **KID-06. Proximal tubule is in the cortex (convoluted part, S1-S2) and the outer stripe (straight part, S3);
  it is absent from the inner stripe and inner medulla.**
  - Rule: `kind=region, types=[Proximal tubule], regions=[cortex, outer stripe], mode=only, strength=strong`
  - Without drawn regions: none.
  - Exceptions and caveats: the S3 boundary defines the outer/inner stripe border, so cells at that border are
    mixed. Adaptive/injured PT loses LRP2/CUBN and gains VCAM1/HAVCR1 [6], which can push them to other labels.
  - Confidence: high; strong.
  - Evidence: [1], [2], [3], [5]

- **KID-07. Thin limbs are medullary: descending thin limbs begin at the outer/inner stripe border and run into the
  inner medulla; ascending thin limbs exist only in the inner medulla; no thin limbs in the cortex.**
  - Rule: `kind=region, types=[Thin limb], regions=[cortex], mode=exclude, strength=moderate`; and
    `kind=region, types=[Thin limb], regions=[inner stripe, inner medulla], mode=enrich, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats: the list's thin-limb markers (CRYAB, AKR1B1, SH3GL3) are weak and not specific, so
    cortical "Thin limb" calls are common artifacts. Short-loop nephrons turn in the inner stripe (their thin
    descending limb is short).
  - Confidence: high (anatomy), medium (calls); moderate.
  - Evidence: [1], [2], [5]

- **KID-08. Thick ascending limbs run in the outer medulla (medullary TAL, both stripes) and cortex (cortical TAL in
  medullary rays and labyrinth); they are absent from the inner medulla.**
  - Rule: `kind=region, types=[Thick ascending limb], regions=[inner medulla, glomerulus], mode=exclude,
    strength=strong`
  - Without drawn regions: none.
  - Exceptions and caveats: UMOD is secreted and can contaminate casts in tubule lumens; cells at the
    inner/outer medulla border are mixed. The macula densa (end of TAL) touches its glomerulus.
  - Confidence: high; strong.
  - Evidence: [1], [2], [3]

- **KID-09. Distal convoluted tubule and connecting tubule are in the cortex only (cortical labyrinth).**
  - Rule: `kind=region, types=[Distal convoluted tubule, Connecting tubule], regions=[cortex], mode=only,
    strength=strong` (use `moderate` for Connecting tubule, whose markers overlap principal cells)
  - Without drawn regions: none.
  - Exceptions and caveats: connecting tubule arcades of juxtamedullary nephrons reach the corticomedullary
    junction. CALB1 and SCNN1G are shared with neighboring segments.
  - Confidence: high; strong (DCT) / moderate (CNT).
  - Evidence: [1], [2], [3]

- **KID-10. Principal cells occur along the whole collecting duct (cortical collecting ducts in medullary rays,
  outer and inner medullary collecting ducts) and dominate the inner medulla; Intercalated cells are in connecting
  tubule, cortical and outer medullary collecting ducts, become sparse in the initial inner medullary collecting
  duct and are absent from the deeper inner medulla.**
  - Rule: `kind=region, types=[Intercalated cells], regions=[inner medulla], mode=exclude, strength=moderate`;
    `kind=region, types=[Principal cells, Intercalated cells], regions=[glomerulus], mode=exclude, strength=strong`
  - Without drawn regions: none.
  - Exceptions and caveats: evidence for IC loss in the inner medulla is mostly rat/mouse histology [11] (title
    verified; finer details from secondary sources); human data are sparse, so keep moderate and draw `inner
    medulla` only for the papilla. Principal-to-intercalated ratios are about 3:1 in cortical and 2:1 in outer
    medullary collecting duct [10]. Transitional IC-PC cells exist in mouse [4]. Urothelium of the pelvis (not in
    the built-in list) shares GATA3/KRT markers with principal cells; add a urothelial type if the pelvis is in
    the section.
  - Confidence: medium; moderate.
  - Evidence: [2], [4], [10], [11]

- **KID-11. Lymphatic EC are in the cortex along interlobular/arcuate vessels and around the hilum; the medulla has
  very few or none.**
  - Rule: `kind=region, types=[Lymphatic EC], regions=[outer stripe, inner stripe, inner medulla, glomerulus],
    mode=exclude, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats: lymphangiogenesis in chronic kidney disease, transplant rejection and fibrosis adds
    lymphatics in the cortical interstitium (not usually in the medulla). LYVE1 is also on some macrophages.
  - Confidence: high; moderate.
  - Evidence: [13]

- **KID-12. Smooth muscle cells belong to arteries and arterioles (interlobar, arcuate, interlobular, afferent and
  efferent arterioles), the capsule and the pelvis wall, not to tubules or glomerular tufts.**
  - Rule: `kind=region, types=[Smooth muscle cells], regions=[large vessels, capsule, pelvis], mode=enrich,
    strength=moderate`; `kind=region, types=[Smooth muscle cells], regions=[glomerulus], mode=exclude,
    strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[Smooth muscle cells], neighbor_types=[Endothelial,
    Fibroblasts, Smooth muscle cells], min_frac=0.3, mode=enrich, strength=weak`
  - Exceptions and caveats: afferent arterioles (and renin cells) enter at the vascular pole and touch the
    glomerular outline; descending vasa recta in the medulla have pericytes rather than smooth muscle. Mesangial
    cells share contractile genes (ACTA2 in injury).
  - Confidence: medium; moderate / weak.
  - Evidence: (textbook histology, not verified online)

- **KID-13. Interstitial Fibroblasts and Pericytes are sparse in the cortex and more abundant in the inner medulla
  (medullary interstitial cells); they are not on glomerular tufts.**
  - Rule: `kind=region, types=[Fibroblasts, Pericytes], regions=[inner medulla], mode=enrich, strength=weak`;
    `kind=region, types=[Fibroblasts], regions=[glomerulus], mode=exclude, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats: interstitial fibrosis (CKD, transplant) expands fibroblasts/myofibroblasts in the
    cortex [6]; periglomerular fibrosis rings glomeruli. PDGFRB is shared with Mesangial cells (see KID-03).
  - Confidence: medium; weak.
  - Evidence: [6], [7] (medullary interstitial cells: textbook histology, not verified online)

- **KID-14. Resident macrophages and dendritic cells are spread through the interstitium; antibacterial
  macrophages and neutrophils are positioned toward the medulla and pelvis, the regions most exposed to
  ascending infection.**
  - Rule: `kind=region, types=[Neutrophils], regions=[inner medulla, pelvis], mode=enrich, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats: neutrophils are rare in normal tissue; pyelonephritis, acute tubular injury,
    glomerulonephritis and rejection bring leukocytes everywhere, including glomeruli. Do not use exclude rules
    for immune cells in kidney.
  - Confidence: low; weak.
  - Evidence: [12]

**References**
[1] Kriz W, Bankir L (Renal Commission of the IUPS), A standard nomenclature for structures of the kidney, Kidney
International, 1988 (33:1-7).
[2] Chen L, Clark JZ, Nelson JW, Kaissling B, Ellison DH, Knepper MA, Renal-tubule epithelial cell nomenclature for
single-cell RNA-sequencing studies, Journal of the American Society of Nephrology, 2019 (30:1358-1364).
[3] Lee JW, Chou CL, Knepper MA, Deep sequencing in microdissected renal tubules identifies nephron
segment-specific transcriptomes, Journal of the American Society of Nephrology, 2015. PMID 25817355.
[4] Park J et al., Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney
disease, Science, 2018. DOI 10.1126/science.aar2131.
[5] Ransick A et al., Single-cell profiling reveals sex, lineage, and regional diversity in the mouse kidney,
Developmental Cell, 2019. DOI 10.1016/j.devcel.2019.10.005.
[6] Lake BB et al., An atlas of healthy and injured cell states and niches in the human kidney, Nature, 2023.
DOI 10.1038/s41586-023-05769-3. PMID 37468583.
[7] Hansen J et al., A reference tissue atlas for the human kidney, Science Advances, 2022.
DOI 10.1126/sciadv.abn4965. PMID 35675394.
[8] Shankland SJ, Smeets B, Pippin JW, Moeller MJ, The emergence of the glomerular parietal epithelial cell, Nature
Reviews Nephrology, 2014 (10:158-173).
[9] Schlondorff D, Banas B, The mesangial cell revisited: no cell is an island, Journal of the American Society of
Nephrology, 2009. DOI 10.1681/ASN.2008050549. PMID 19470685.
[10] Roy A, Al-bataineh MM, Pastor-Soler NM, Collecting duct intercalated cell function and regulation, Clinical
Journal of the American Society of Nephrology, 2015. PMID 25632105.
[11] Madsen KM, Clapp WL, Verlander JW, Structure and function of the inner medullary collecting duct, Kidney
International, 1988. DOI 10.1038/ki.1988.201.
[12] Stewart BJ et al., Spatiotemporal immune zonation of the human kidney, Science, 2019.
DOI 10.1126/science.aat5031. PMID 31604275.
[13] Russell PS et al., Renal lymphatics: anatomy, physiology, and clinical implications, Frontiers in Physiology,
2019. DOI 10.3389/fphys.2019.00251.

---

## Lung

Type names are from `lung.txt` (+ `@immune`, `@stroma`): Epithelial (AT1 cells, AT2 cells, Club cells, Ciliated
cells, Basal cells, Goblet cells), Mesothelial cells, Alveolar macrophages (under Myeloid), and the shared immune and
stromal types. Not in the built-in list: submucosal gland serous / mucous / duct cells, gland myoepithelial cells,
neuroendocrine cells, ionocytes, tuft cells, respiratory airway secretory (RAS) cells, interstitial macrophages,
chondrocytes, alveolar vs adventitial fibroblast subtypes, capillary subtypes.

**Compartments to draw**
- `airway epithelium`: the lining of bronchi and bronchioles. It is tall, pseudostratified and ciliated in bronchi
  (cilia and a dark basal nuclear row are visible) and becomes a single cuboidal layer in bronchioles.
- `airway wall`: the tissue under the airway lining, with lamina propria, a smooth-muscle band, cartilage plates and
  submucosal glands. Cartilage and glands occur only in bronchi.
- `submucosal gland`: acini of pale, foamy mucous cells and granular, eosinophilic serous cells in the bronchial wall,
  between the epithelium and the cartilage. They are found only in cartilaginous airways.
- `cartilage`: hyaline cartilage plates (glassy blue-grey matrix with chondrocytes in lacunae) in bronchial walls.
- `alveoli`: lace-like thin septa around empty airspaces, which is most of the section. AT2 cells are the cuboidal
  cells at septal corners. Alveolar macrophages are rounded cells lying free in the airspaces.
- `bronchovascular bundle`: a thick connective-tissue sheath around an airway and its paired pulmonary artery, with
  lymphatics and lymphoid aggregates.
- `pleura`: the thin fibrous layer on the outer lung surface, lined by a single row of flat mesothelial cells.
- `interlobular septum`: thin fibrous lines that run in from the pleura and carry veins and lymphatics.

**Depth rule.** The app's depth core is built from EPCAM and KRT8 density. In lung this core covers both the alveolar
parenchyma (AT2 cells are strongly EPCAM/KRT8-positive) and the airway lining, so almost every alveolar-region cell has
a depth near 0. Depth therefore separates only the thick epithelium-free structures (bronchovascular sheaths,
large-vessel walls, cartilage, thick pleura or fibrotic scars) from everything else. It is **weakly meaningful**: it
can stand in for "not in a thick connective-tissue sheath", but it cannot tell airways from alveoli. Use drawn regions
for that.

**Priors**

- **LUNG-01. AT1 and AT2 cells line the alveoli. They are absent from bronchial walls, cartilage, bronchovascular
  sheaths and the pleural surface.**
  - Rule: `kind=region, types=[AT1 cells, AT2 cells], regions=[alveoli], mode=only, strength=moderate`
  - Without drawn regions: `kind=depth, types=[AT1 cells, AT2 cells], max_depth_um=30, mode=only, strength=weak`.
    This only catches calls deep inside connective-tissue sheaths.
  - Exceptions and caveats:
    - At the bronchioalveolar duct junction, and in human respiratory bronchioles, AT2-like and airway secretory
      cells mix; RAS cells co-express SCGB3A2 with some AT2 genes [8].
    - In fibrosis ("bronchiolization" and honeycomb cysts), AT2 hyperplasia and lung adenocarcinoma (SFTPC/NAPSA can
      stay positive), cells keep their identity but leave the normal compartment.
    - A loosely drawn `alveoli` region will include small bronchioles.
  - Confidence: high. Strength: moderate (strong only for normal tissue with carefully drawn regions).
  - Evidence: [1][4][7]

- **LUNG-02. Basal cells sit on the basement membrane of the conducting airways. They are absent from alveoli. In
  human they extend from the trachea down into small airways; in mouse they are restricted to the trachea, the main
  bronchi and the first intralobar branches, and mouse bronchioles have essentially none.**
  - Rule: `kind=region, types=[Basal cells], regions=[airway epithelium, submucosal gland], mode=only,
    strength=moderate`. For mouse intralobar bronchioles, also
    `kind=region, types=[Basal cells], regions=[bronchiole], mode=exclude, strength=moderate`.
  - Without drawn regions: `kind=neighborhood, types=[Basal cells], neighbor_types=[Ciliated cells, Club cells,
    Goblet cells, Basal cells], min_frac=0.3, mode=only, strength=weak`
  - Exceptions and caveats:
    - After severe injury (influenza, bleomycin in mouse), KRT5+ "pods" appear in the alveolar region [4].
    - IPF contains KRT17+ aberrant basaloid cells and basal-cell-lined honeycomb cysts; squamous metaplasia and
      squamous carcinoma also break this rule.
    - Submucosal gland myoepithelial and duct basal cells express KRT5/KRT14 and TP63 and will be typed as Basal cells
      (not in the built-in list), which is why the rule includes the gland region [2][4].
  - Confidence: high. Strength: moderate.
  - Evidence: [1][2][4][5]

- **LUNG-03. Ciliated cells line bronchi and bronchioles, never alveolar septa.**
  - Rule: `kind=region, types=[Ciliated cells], regions=[airway epithelium, submucosal gland], mode=only,
    strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Ciliated cells], neighbor_types=[AT1 cells, AT2 cells],
    min_frac=0.6, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Bronchiolization of alveoli in fibrosis.
    - Ciliated cells shed into the airspaces in damaged tissue.
    - Gland ducts contain ciliated cells, which is why the gland region is included.
    - The atlas distinguishes proximal and distal ciliated states [1], but the built-in list does not.
  - Confidence: high. Strength: strong.
  - Evidence: [1][2][4]

- **LUNG-04. Goblet cells are a feature of the larger (cartilaginous) airways and become rare in small bronchioles.
  They are absent from alveoli.**
  - Rule: `kind=region, types=[Goblet cells], regions=[airway epithelium, submucosal gland], mode=only,
    strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Goblet cells], neighbor_types=[AT1 cells, AT2 cells],
    min_frac=0.6, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Goblet-cell hyperplasia and metaplasia in asthma, COPD and smokers extends goblet cells distally [4].
    - Mucinous adenocarcinoma expresses MUC5AC/MUC5B in the parenchyma.
    - Normal mouse airways have very few goblet cells.
    - Submucosal gland mucous cells (MUC5B, BPIFB1) will be typed as Goblet cells (gland cells are not in the
      built-in list), which is why the gland region is included [2].
  - Confidence: high. Strength: moderate.
  - Evidence: [1][2][4]

- **LUNG-05. Club cells are airway secretory cells. In human they are most abundant in small airways and depleted in
  the trachea; in mouse they are the main secretory cell all along the airways, including the trachea. In both
  species they are absent from the alveolar septa.**
  - Rule: `kind=region, types=[Club cells], regions=[airway epithelium, submucosal gland], mode=only,
    strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Club cells], neighbor_types=[AT1 cells, AT2 cells],
    min_frac=0.7, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Human terminal and respiratory bronchioles contain RAS cells (SCGB3A2+ SCGB1A1+) that sit right next to
      alveoli and can make AT2 cells [8]. Mice have no respiratory bronchioles [8].
    - Mouse Scgb1a1+ cells rarely contribute to alveoli at steady state [6], but after injury their descendants
      appear in the alveoli [4].
    - SCGB1A1/SCGB3A1 transcripts spill over from the airway lumen and mucus.
  - Confidence: high (normal), medium at the bronchioalveolar junction. Strength: moderate.
  - Evidence: [1][2][4][6][8]

- **LUNG-06. Airway submucosal glands are found only in cartilaginous bronchi (in mouse, only in the proximal
  trachea). IgA plasma cells are concentrated around these glands, especially their serous acini.**
  - Rule: `kind=region, types=[Plasma cells], regions=[submucosal gland], mode=enrich, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Plasma cells], neighbor_types=[Goblet cells], min_frac=0.3,
    mode=enrich, strength=weak`. Gland mucous cells type as Goblet cells.
  - Exceptions and caveats:
    - IgG plasma cells sit in immune infiltrates rather than in glands [2].
    - Peripheral lung samples contain no glands, so this rule applies only to proximal sections.
  - Confidence: high (human). Strength: moderate.
  - Evidence: [2]

- **LUNG-07. Alveolar macrophages lie free in the alveolar airspaces and on the alveolar lining. They are not
  resident in airway walls, the pleura or bronchovascular sheaths, where interstitial macrophages are found instead.**
  - Rule: `kind=region, types=[Alveolar macrophages], regions=[alveoli], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Alveolar macrophages], neighbor_types=[AT1 cells, AT2 cells,
    Alveolar macrophages, Endothelial], min_frac=0.4, mode=only, strength=weak`
  - Exceptions and caveats:
    - Alveolar macrophages also travel up the small airways and appear in airway lumens; draw `alveoli` to include
      bronchiolar lumens if needed.
    - In pneumonia, COVID-19 and fibrosis, monocyte-derived macrophages (e.g. SPP1+ profibrotic) appear in the
      interstitium and can resemble alveolar macrophages [3].
    - Interstitial macrophages are not in the built-in list and will be typed as Macrophages.
    - Airspace location is textbook histology; the lung-restriction of alveolar macrophages is from [1].
  - Confidence: high. Strength: moderate.
  - Evidence: [1][3]; airspace location: (textbook histology, not verified online)

- **LUNG-08. Mesothelial cells exist only on the pleural surface, as a single flat layer.**
  - Rule: `kind=region, types=[Mesothelial cells], regions=[pleura], mode=only, strength=strong`
  - Without drawn regions: none. There is no surface marker, and in this panel context the depth core says nothing
    about the pleura.
  - Exceptions and caveats:
    - Mesothelioma, pleural plaques and reactive mesothelial hyperplasia.
    - WT1/MSLN/CALB2 are also expressed at low levels by some fibroblasts and tumors, so check the marker evidence
      before applying a strong rule to calls inside the parenchyma.
    - Sections without a pleural edge should have essentially zero mesothelial cells.
  - Confidence: high. Strength: strong.
  - Evidence: [4] (the lung is encased in mesothelium)

- **LUNG-09. Smooth muscle cells form bands in airway walls and in the walls of arteries and veins. Alveolar septa
  contain no true smooth muscle, only myofibroblasts and pericytes.**
  - Rule: `kind=region, types=[Smooth muscle cells], regions=[alveoli], mode=exclude, strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[Smooth muscle cells], neighbor_types=[AT1 cells, AT2 cells],
    min_frac=0.6, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Small arterioles and alveolar-duct smooth-muscle knobs lie within any drawn "alveoli" region; hence the weak
      strength.
    - Pulmonary hypertension (muscularized arterioles), LAM and fibrosis (myofibroblast foci expressing ACTA2/TAGLN).
  - Confidence: medium. Strength: weak.
  - Evidence: [4] (smooth muscle with intralobar airways down to small bronchioles)

- **LUNG-10. Lymphatic EC are in bronchovascular bundles, interlobular septa and the pleura. They are rare in
  alveolar septa.**
  - Rule: `kind=region, types=[Lymphatic EC], regions=[alveoli], mode=exclude, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Small lymphatics accompany alveolar ducts and venules.
    - Lymphangiogenesis in fibrosis and LAM.
    - TFF3 is in the Lymphatic EC list and is also a goblet/secretory marker, so watch for airway spill-over.
  - Confidence: medium. Strength: weak.
  - Evidence: (textbook histology, not verified online)

- **LUNG-11. Organized lymphoid aggregates (B cells with T cells) sit along bronchovascular bundles and the pleura,
  not inside alveolar septa.**
  - Rule: `kind=region, types=[B cells], regions=[bronchovascular bundle, pleura, airway wall], mode=enrich,
    strength=weak`
  - Without drawn regions: `kind=neighborhood, types=[B cells], neighbor_types=[B cells, T and NK cells],
    min_frac=0.5, mode=enrich, strength=weak` (use only after a first pass).
  - Exceptions and caveats:
    - Lymphoid aggregates are inconspicuous in healthy adults and expand in COPD, rheumatoid lung and infection.
    - Intravascular and circulating B cells occur anywhere.
  - Confidence: medium. Strength: weak.
  - Evidence: [2] (immune infiltrates and gland-associated B cells); localization otherwise (textbook histology, not
    verified online)

- **LUNG-12. Cartilage plates of the bronchi contain chondrocytes (not in the built-in list). Cells inside a
  `cartilage` region should not be typed as epithelial or immune.**
  - Rule: `kind=region, types=[Epithelial, Myeloid, T and NK cells, B and plasma cells], regions=[cartilage],
    mode=exclude, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Chondrocytes will be forced into the nearest stromal type (usually Fibroblasts). Flag them for a custom marker
      list (COL2A1, ACAN, SOX9) rather than trusting the label.
  - Confidence: high. Strength: moderate.
  - Evidence: [4] (cartilage in intralobar airway walls)

**References**
[1] Travaglini KJ et al., A molecular cell atlas of the human lung from single-cell RNA sequencing, Nature, 2020.
doi:10.1038/s41586-020-2922-4 (PMID 33208946)
[2] Madissoon E et al., A spatially resolved atlas of the human lung characterizes a gland-associated immune niche,
Nature Genetics, 2023 (online Dec 2022). doi:10.1038/s41588-022-01243-4 (PMID 36543915)
[3] Sikkema L et al., An integrated cell atlas of the lung in health and disease, Nature Medicine, 2023.
doi:10.1038/s41591-023-02327-2 (PMID 37291214)
[4] Hogan BLM et al., Repair and regeneration of the respiratory system: complexity, plasticity, and mechanisms of
lung stem cell function, Cell Stem Cell, 2014. doi:10.1016/j.stem.2014.07.012 (PMID 25105578)
[5] Rock JR et al., Basal cells as stem cells of the mouse trachea and human airway epithelium, PNAS, 2009.
doi:10.1073/pnas.0906850106 (PMID 19625615)
[6] Rawlins EL et al., The role of Scgb1a1+ Clara cells in the long-term maintenance and repair of lung airway, but
not alveolar, epithelium, Cell Stem Cell, 2009. doi:10.1016/j.stem.2009.04.002 (PMID 19497281)
[7] Barkauskas CE et al., Type 2 alveolar cells are stem cells in adult lung, Journal of Clinical Investigation, 2013.
doi:10.1172/JCI68782 (PMID 23921127)
[8] Basil MC et al., Human distal airways contain a multipotent secretory cell that can regenerate alveoli, Nature,
2022. doi:10.1038/s41586-022-04552-0 (PMID 35355013)


## Heart

Type names are from `heart.txt` (+ `@immune`, `@stroma`): Cardiomyocytes (Ventricular cardiomyocytes, Atrial
cardiomyocytes), Mesothelial cells, Adipocytes, and the shared immune and stromal types. Not in the built-in list:
endocardial cells, pacemaker / conduction-system cells (SAN P cells, AVN, Purkinje), valve interstitial cells,
neuronal and Schwann cells of cardiac nerves and ganglia.

**Compartments to draw**
- `atrium`: thin-walled myocardium with pectinate muscles. Atrial cardiomyocytes are smaller and more loosely packed
  than ventricular ones.
- `ventricle`: thick, compact myocardium of large striated cardiomyocytes with central nuclei. It includes the
  interventricular septum and the apex.
- `epicardium`: the outer surface, with a single flat mesothelial layer over a subepicardial layer of fat, coronary
  vessels, nerves and lymphatics.
- `epicardial fat`: large clear adipocytes on the outer surface, around the coronary arteries and in the
  atrioventricular groove.
- `endocardium`: the thin inner lining of the chambers, with one layer of endothelium over a fibroelastic layer
  (thicker in atria). Subendocardial Purkinje fibres in the ventricles are large and pale.
- `valve`: avascular fibrous leaflets attached at the annulus.
- `large vessel`: the aorta or pulmonary trunk wall (elastic lamellae and smooth muscle), if present.

**Depth rule.** **Not meaningful.** Heart has no EPCAM/KRT8-dense epithelium. The only KRT8+ cells are the epicardial
mesothelium, a monolayer that rarely meets the core density threshold, so depth is usually undefined. Where it does
appear it marks only the outer surface. Use drawn regions and neighborhood rules.

**Priors**

- **HEART-01. Atrial cardiomyocytes (NPPA, MYL7, MYL4) occur only in the atria.**
  - Rule: `kind=region, types=[Atrial cardiomyocytes], regions=[atrium], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Atrial cardiomyocytes], neighbor_types=[Ventricular
    cardiomyocytes], min_frac=0.6, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Stressed or failing ventricles re-express NPPA and NPPB. Kanemaru et al. describe a ventricular "stress" niche
      with NPPB/XIRP2/ACTA1 [3]. NPPA is in the Atrial cardiomyocytes list, so in heart failure, hypertrophy or
      after infarction, NPPA-high ventricular cells are a real finding and this rule should not be strong.
    - Conduction-system cells (AVN, His) express atrial-like myosins.
    - Myocardial sleeves extend into the pulmonary veins and venae cavae.
  - Confidence: high (normal), medium (diseased). Strength: moderate.
  - Evidence: [1][2][3]

- **HEART-02. Ventricular cardiomyocytes (MYL2, MYL3, IRX4) occur only in the ventricles and septum.**
  - Rule: `kind=region, types=[Ventricular cardiomyocytes], regions=[ventricle], mode=only, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Ventricular cardiomyocytes], neighbor_types=[Atrial
    cardiomyocytes], min_frac=0.6, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - MYL2 has low-level atrial expression in some diseased atria.
    - Sections cut through the AV junction mix both types along the annulus; draw the border carefully.
  - Confidence: high. Strength: strong.
  - Evidence: [1] (ventricular CMs enriched in MYH7 and MYL2; markedly distinct atrial vs ventricular signatures)

- **HEART-03. Mesothelial cells exist only in the epicardium (the visceral pericardium on the heart surface).**
  - Rule: `kind=region, types=[Mesothelial cells], regions=[epicardium], mode=only, strength=strong`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - In adult heart the epicardium is essentially mesothelial and immune-responsive; fetal epicardium has
      additional migratory fibroblast-like states [4].
    - Pericarditis and adhesions thicken the layer.
    - Sections without the outer surface should contain no mesothelial cells.
    - WT1 is also re-expressed in some injured-heart fibroblasts.
  - Confidence: high. Strength: strong.
  - Evidence: [1] (BNC1+/CDH5- mesothelial cells localized to the epicardium by smFISH), [3], [4]

- **HEART-04. Adipocytes sit in the epicardial and subepicardial fat (around coronary vessels, in the AV groove).
  They are absent from compact myocardium and endocardium.**
  - Rule: `kind=region, types=[Adipocytes], regions=[epicardium, epicardial fat], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Adipocytes], neighbor_types=[Cardiomyocytes], min_frac=0.7,
    mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Intramyocardial fat in the right ventricle is common with age and obesity.
    - Arrhythmogenic cardiomyopathy causes fibrofatty replacement of the myocardium, and chronic infarcts show fatty
      metaplasia. In these, adipocytes inside the myocardium are the real finding.
    - Adipocytes are large and poorly captured by single-cell methods, so their transcripts often spill into
      neighbors.
  - Confidence: high (normal LV), medium (RV). Strength: moderate.
  - Evidence: [3] (adipocytes enriched in the epicardial niche of all four chambers); [1] (cardiac adipocytes, with
    large epicardial fat deposits intentionally excluded from sampling)

- **HEART-05. The epicardium is an immune niche. Plasma cells (IgG and IgA) are enriched there, together with
  fibroblasts and lymphatic endothelium.**
  - Rule: `kind=region, types=[Plasma cells, Lymphatic EC], regions=[epicardium, epicardial fat], mode=enrich,
    strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Plasma cells], neighbor_types=[Adipocytes, Mesothelial cells],
    min_frac=0.2, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - Plasma cells and lymphatics also occur sparsely in the myocardium, so use `enrich`, never `only`.
  - Confidence: medium-high. Strength: moderate.
  - Evidence: [3] (enrichment of IgG+ and IgA+ plasma cells, mesothelial cells, fibroblasts, lymphatic EC and
    adipocytes in the epicardium)

- **HEART-06. Endocardial cells (not in the built-in list) line the chamber cavities. In atria they are an NPR3+,
  SMOC1+ endothelial state, and they will be typed as Endothelial / Blood vessel EC.**
  - Rule: `kind=region, types=[Endothelial], regions=[endocardium], mode=enrich, strength=weak`. This is only useful
    for flagging the region; it should not change other types.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Endothelium is everywhere in myocardium (capillaries), so never use `only` for Endothelial.
    - Add an "Endocardial cells" subtype (NPR3, SMOC1, NPR3+ PECAM1+) to a custom list if needed.
  - Confidence: medium. Strength: weak.
  - Evidence: [1] (mainly atrial EC7_atria expressing SMOC1 and NPR3, suggestive of endocardial cells)

- **HEART-07. Cardiomyocytes are absent from valve leaflets, the epicardial fat pad, the pericardium and the walls of
  the great arteries.**
  - Rule: `kind=region, types=[Cardiomyocytes], regions=[valve, epicardial fat, large vessel], mode=exclude,
    strength=strong`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Myocardium extends a short way into the base of the AV valves and into the pulmonary-vein and caval sleeves.
    - Large cardiomyocyte transcript clouds (TTN, MYH7) spill into adjacent fat and endothelium; this rule cleans that
      up but can hide tissue-folding artifacts.
  - Confidence: high. Strength: strong.
  - Evidence: (textbook histology, not verified online)

- **HEART-08. Smooth muscle cells in the heart belong to arterial and arteriolar walls (coronary arteries run in the
  epicardial fat). They are not scattered among cardiomyocytes. Pericytes, by contrast, surround the capillaries
  throughout the myocardium.**
  - Rule: `kind=neighborhood, types=[Smooth muscle cells], neighbor_types=[Endothelial, Smooth muscle cells,
    Pericytes, Fibroblasts], min_frac=0.3, mode=only, strength=weak`
  - Without drawn regions: same neighborhood rule.
  - Exceptions and caveats:
    - Cardiomyocytes and SMC share ACTA2/TAGLN in stressed myocardium.
    - Myofibroblasts in infarct scars express ACTA2 and can type as smooth muscle.
  - Confidence: medium. Strength: weak.
  - Evidence: [3] (arterial SMC2_art in vascular niches); (textbook histology, not verified online)

- **HEART-09. Chamber identity extends to non-myocytes. Fibroblast states differ between atria and ventricles, and
  atria contain the endocardial-like EC state. This does not translate into a built-in-type rule.**
  - Rule: none. This is informational: do not suppress Fibroblasts in either chamber.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Atrial fibroblasts show stronger profibrotic responses [1]; atrial fibrosis is common in AF.
  - Confidence: high. Strength: n/a.
  - Evidence: [1][2]

- **HEART-10. The sinoatrial node is a discrete subepicardial structure at the junction of the crista terminalis and
  the right-atrial posterior wall. Its pacemaker cells (not in the built-in list) are small, rounder myocytes embedded
  in dense fibrous tissue.**
  - Rule: `kind=region, types=[Ventricular cardiomyocytes], regions=[SAN], mode=exclude, strength=strong`. Pacemaker
    cells will be typed Atrial cardiomyocytes or Cardiomyocytes.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Only relevant for sections that sample the SAN or AVN. Pacemaker cells need a custom subtype (HCN1, HCN4,
      SHOX2, TBX3).
  - Confidence: high. Strength: strong.
  - Evidence: [3]

- **HEART-11. Resident macrophages are distributed throughout the myocardium interstitium. Neutrophils and dense
  monocyte infiltrates are not part of normal myocardium.**
  - Rule: none for macrophages. For a normal-control sample only:
    `kind=region, types=[Neutrophils], regions=[ventricle, atrium], mode=exclude, strength=weak`.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Acute infarction, myocarditis and transplant rejection bring neutrophils and monocytes into the myocardium; do
      not use this rule on diseased tissue.
    - Intravascular leukocytes are present everywhere.
  - Confidence: medium. Strength: weak.
  - Evidence: [2] (two groups of resident macrophages); (textbook histology, not verified online)

**References**
[1] Litviňuková M et al., Cells of the adult human heart, Nature, 2020. doi:10.1038/s41586-020-2797-4 (PMID 32971526)
[2] Tucker NR et al., Transcriptional and cellular diversity of the human heart, Circulation, 2020.
doi:10.1161/CIRCULATIONAHA.119.045401 (PMID 32403949)
[3] Kanemaru K et al., Spatially resolved multiomics of human cardiac niches, Nature, 2023.
doi:10.1038/s41586-023-06311-1 (PMID 37438528)
[4] Knight-Schrijver VR et al., A single-cell comparison of adult and fetal human epicardium defines the
age-associated changes in epicardial activity, Nature Cardiovascular Research, 2022. doi:10.1038/s44161-022-00183-w
(PMID 36938497)


## Breast

Type names are from `breast.txt` (+ `@immune`, `@stroma`): Epithelial (Luminal hormone-sensing, Luminal progenitors,
Basal / myoepithelial cells), Adipocytes, and the shared immune and stromal types. "Luminal progenitors" corresponds
to the atlas's LumSec (luminal secretory / alveolar) lineage; "Luminal hormone-sensing" corresponds to LumHR. Not in
the built-in list: basal-luminal hybrid cells, intralobular vs interlobular fibroblast states, and nipple
keratinocytes (the skin over the nipple).

**Compartments to draw**
- `duct`: tubular, bilayered epithelium (a cuboidal or columnar luminal layer with a flatter myoepithelial layer
  under it) in dense collagenous stroma. Lumens are often open.
- `lobule`: a TDLU, a tight cluster of small acini (alveoli) in loose, paler, cellular intralobular stroma. It is the
  most cell-dense region.
- `intralobular stroma`: the loose, pale, cellular connective tissue inside a lobule, between acini.
- `connective tissue`: dense pink interlobular collagen with sparse spindle nuclei and vessels.
- `adipose`: large clear empty-looking adipocytes with thin rims.
- `nipple / skin`: squamous epidermis and lactiferous sinuses, if the section includes them.

**Depth rule.** **Meaningful.** Luminal cells are strongly EPCAM/KRT8-positive and densely packed, so the depth core
marks ducts and lobules. Depth then measures distance from the ductal-lobular epithelium: myoepithelial cells sit at
about one cell (under ~15 um), intralobular stroma within tens of um, and interlobular connective tissue and fat are
far away. The thresholds below are starting values to calibrate on the sample, not literature numbers.

**Priors**

- **BR-01. All epithelial cells (luminal and basal / myoepithelial) are confined to ducts and lobules. Normal
  adipose and interlobular connective tissue contain no epithelium.**
  - Rule: `kind=region, types=[Epithelial], regions=[duct, lobule], mode=only, strength=strong`
  - Without drawn regions: `kind=depth, types=[Epithelial], max_depth_um=20, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Invasive carcinoma (ductal or lobular) puts luminal-like cells in the stroma and fat. On a tumor section this
      rule would hide the invasion, so do not use it there, or only outside a drawn `tumor` region.
    - Lobules small enough to be missed when drawing.
  - Confidence: high (normal tissue). Strength: strong (normal), off (tumor).
  - Evidence: [1]

- **BR-02. The breast epithelium is bilayered. Basal / myoepithelial cells form a continuous outer layer against the
  basement membrane, and luminal cells form the inner layer facing the lumen.**
  - Rule: `kind=neighborhood, types=[Basal / myoepithelial cells], neighbor_types=[Luminal hormone-sensing, Luminal
    progenitors], min_frac=0.2, mode=only, strength=moderate`
  - Without drawn regions: same, or `kind=depth, types=[Basal / myoepithelial cells], max_depth_um=15, mode=only,
    strength=moderate`
  - Exceptions and caveats:
    - The Basal / myoepithelial markers include ACTA2, MYLK and KRT17, which are shared with vascular smooth muscle,
      pericytes and myofibroblasts. This prior is mainly what stops vessels being called myoepithelium.
    - In DCIS the myoepithelial layer persists around expanded ducts; in invasive carcinoma it is lost. Tumor-boundary
      myoepithelial states were described by Xenium [4].
    - Basal-luminal hybrid cells accumulate with age [3].
  - Confidence: high. Strength: moderate.
  - Evidence: [1][3][4]

- **BR-03. Luminal hormone-sensing (LumHR: ESR1, PGR, AR, ANKRD30A) cells are relatively enriched in lobules. Both
  luminal lineages exist in both ducts and lobules.**
  - Rule: `kind=region, types=[Luminal hormone-sensing], regions=[lobule], mode=enrich, strength=weak`
  - Without drawn regions: none (ducts and lobules both have depth 0).
  - Exceptions and caveats:
    - The difference is one of proportion, not presence, so never use `only` or `exclude`.
    - ER+ tumors expand LumHR-like cells anywhere.
    - Parity, age and menstrual phase shift the proportions [2][3].
  - Confidence: medium. Strength: weak.
  - Evidence: [1] (lobular regions enriched for LumHR genes ANKRD30A, AR, ESR1, PGR; LumHR-SCGB clusters in lobules)

- **BR-04. Luminal progenitor / secretory states are mixed in location. LTF+ SLPI+ KRT15+ LumSec-basal and HLA-high
  states sit mainly in ducts, while KIT+ and ELF5+ LumSec states sit mainly in lobules.**
  - Rule: none recommended for the built-in "Luminal progenitors" type, which merges both. If the list is split, use
    `kind=region, types=[<LTF+ subtype>], regions=[duct], mode=enrich, strength=weak` and
    `kind=region, types=[<KIT+/ELF5+ subtype>], regions=[lobule], mode=enrich, strength=weak`.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Lactation and pregnancy expand the ELF5+ secretory (alveolar) cells in lobules.
    - BRCA1 carriers show expanded luminal progenitors.
  - Confidence: medium. Strength: weak.
  - Evidence: [1][2]

- **BR-05. Adipocytes occupy the adipose regions (and the fat between lobules). They are absent from ducts, lobules
  and intralobular stroma. Breast adipocytes are white adipocytes.**
  - Rule: `kind=region, types=[Adipocytes], regions=[duct, lobule, intralobular stroma], mode=exclude,
    strength=strong`
  - Without drawn regions: `kind=depth, types=[Adipocytes], max_depth_um=30, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Adipocytes are >50 um, lipid-filled and fragile, and are missed by scRNA-seq [1]. In imaging data their thin rim
      gives few transcripts, segmentation often splits them or merges them with neighbors, and their transcripts
      (ADIPOQ, PLIN1, LPL) spill into nearby fibroblasts and endothelium.
    - Fat encroaching on involuted lobules in older women is real; keep the strength moderate for post-menopausal
      tissue.
  - Confidence: high. Strength: strong (premenopausal), moderate (involuted).
  - Evidence: [1]

- **BR-06. Cell density and composition differ by region. Connective (interlobular) tissue is sparse and dominated by
  fibroblasts and vascular endothelium. Intralobular fibroblasts are a distinct, more cellular population
  (FBLN1/SERPINF1-high) close to the acini.**
  - Rule: `kind=region, types=[Fibroblasts], regions=[connective tissue, intralobular stroma], mode=enrich,
    strength=weak`
  - Without drawn regions: none (fibroblasts occur at all depths).
  - Exceptions and caveats:
    - Desmoplastic stroma in tumors.
    - Fibroadenoma and phyllodes tumors have stromal overgrowth inside lobules.
  - Confidence: medium. Strength: weak.
  - Evidence: [1]

- **BR-07. Immune cells are denser in ductal and lobular regions than in the interlobular connective tissue. Tissue-
  resident CD8 T cells (ITGA1+, RUNX3+) are part of the epithelial niche. B cells sit mainly in the stroma around
  ducts and lobules (naive and memory B cells mostly lobular). Macrophages are found in both regions.**
  - Rule: `kind=depth, types=[CD8 T cells], max_depth_um=50, mode=enrich, strength=weak`
  - Without drawn regions: same depth rule; it is region-free.
  - Exceptions and caveats:
    - Tumor-infiltrating lymphocytes and tertiary lymphoid structures (in tumor stroma) change this completely.
    - BRCA carriers show an exhausted immune phenotype [2].
    - Intraductal macrophages and lymphocytes occur in normal ducts.
  - Confidence: medium. Strength: weak.
  - Evidence: [1][2]

- **BR-08. Smooth muscle cells in breast belong to vessel walls (and nipple-areola smooth muscle). They are not
  interleaved with luminal epithelium; a contractile cell touching luminal cells is myoepithelial.**
  - Rule: `kind=neighborhood, types=[Smooth muscle cells, Pericytes], neighbor_types=[Luminal hormone-sensing,
    Luminal progenitors], min_frac=0.3, mode=exclude, strength=moderate`
  - Without drawn regions: same.
  - Exceptions and caveats:
    - Periductal capillaries and pericytes lie within 1-2 cells of the myoepithelium, so do not use `strong`.
    - In the nipple, smooth-muscle bundles surround the lactiferous ducts.
  - Confidence: medium-high. Strength: moderate.
  - Evidence: [1] (bilayer organization); (textbook histology, not verified online)

- **BR-09. Pericytes and vascular endothelium belong to the stroma, both intralobular capillaries and the larger
  vessels of interlobular tissue. They are never inside the epithelial bilayer.**
  - Rule: `kind=region, types=[Endothelial, Pericytes], regions=[duct], mode=exclude, strength=weak`. Use this only if
    the `duct` region is drawn tightly around the epithelium.
  - Without drawn regions: `kind=neighborhood, types=[Endothelial, Pericytes], neighbor_types=[Luminal
    hormone-sensing, Luminal progenitors], min_frac=0.5, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Drawn duct and lobule regions always include periductal capillaries.
  - Confidence: medium. Strength: weak.
  - Evidence: [1] (connective region dominated by fibroblasts and vascular EC)

- **BR-10. The nipple and the skin over it contain squamous keratinocytes, which are not in the breast list. Cells in
  the epidermis will be forced into Basal / myoepithelial cells (KRT5/KRT14/KRT17/TP63).**
  - Rule: `kind=region, types=[Basal / myoepithelial cells, Luminal hormone-sensing, Luminal progenitors],
    regions=[nipple / skin], mode=exclude, strength=strong`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Better: add the skin list (Keratinocytes) to the marker list when the section includes nipple skin.
    - Paget disease puts luminal tumor cells into the nipple epidermis.
  - Confidence: high. Strength: strong.
  - Evidence: (textbook histology, not verified online)

- **BR-11. Physiological state changes composition. Age, parity and germline BRCA status shift epithelial and immune
  states, and postmenopausal involution replaces lobules with fat and collagen.**
  - Rule: none. Use this to lower the strengths of BR-03 to BR-07 for older, parous or lactating donors.
  - Confidence: high. Strength: n/a.
  - Evidence: [2][3]

**References**
[1] Kumar T et al., A spatially resolved single-cell genomic atlas of the adult human breast, Nature, 2023.
doi:10.1038/s41586-023-06252-9 (PMID 37380767)
[2] Reed AD et al., A single-cell atlas enables mapping of homeostatic cellular shifts in the adult human breast,
Nature Genetics, 2024. doi:10.1038/s41588-024-01688-9 (PMID 38548988)
[3] Gray GK et al., A human breast atlas integrating single-cell proteomics and transcriptomics, Developmental Cell,
2022. doi:10.1016/j.devcel.2022.05.003 (PMID 35617956)
[4] Janesick A et al., High resolution mapping of the tumor microenvironment using integrated single-cell, spatial
and in situ analysis, Nature Communications, 2023. doi:10.1038/s41467-023-43458-x


## Skin

Type names are from `skin.txt` (+ `@immune`, `@stroma`): Keratinocytes (Basal keratinocytes, Suprabasal
keratinocytes, Granular keratinocytes), Melanocytes, Sweat gland cells, Langerhans cells (under Myeloid), Adipocytes,
and the shared immune and stromal types. Not in the built-in list: hair-follicle keratinocytes (outer and inner root
sheath, matrix, bulge), dermal papilla and dermal sheath fibroblasts, sebocytes, Merkel cells, papillary vs reticular
fibroblast states, arrector pili muscle (typed as Smooth muscle cells).

**Compartments to draw**
- `epidermis`: the stratified squamous layer on top. From the bottom up it has a basal row of small dark cuboidal
  cells on the basement membrane, pink polygonal spinous cells with "prickles", a granular layer with purple
  keratohyalin granules, and anucleate stratum corneum.
- `papillary dermis`: the thin layer of loose, fine collagen just under the epidermis, in the dermal papillae, with
  capillary loops.
- `reticular dermis`: thick, coarse pink collagen bundles below the papillary dermis.
- `hair follicle`: an epithelial tube or bulb running down into the dermis or subcutis, with a dermal papilla at the
  bulb and sebaceous glands (foamy lobules) and the arrector pili attached.
- `sweat gland`: eccrine coils (clusters of small tubules in cross-section, with pale secretory cells and darker
  two-layered ducts) in the deep dermis or at the dermal-subcutaneous junction.
- `subcutis`: subcutaneous fat lobules of large clear adipocytes separated by fibrous septa.

**Depth rule.** **Not meaningful for the epidermis in its current form.** The app's depth core uses EPCAM and KRT8.
Interfollicular keratinocytes express KRT5/KRT14 (basal) and KRT1/KRT10 (suprabasal) and are essentially KRT8-negative
and EPCAM-low (textbook, not verified online). The core will therefore mostly mark sweat glands (KRT8/KRT18/KRT19+),
plus Merkel cells or a basal cell carcinoma if present, and "depth" becomes distance to sweat coils, not to the
epidermis. A distance-to-epidermis rule would need a keratinocyte-based core (e.g. KRT5, KRT14, KRT10), which is a
software change and is not proposed here. All skin priors below therefore use drawn regions or neighborhood rules.

**Priors**

- **SKIN-01. Keratinocytes are confined to the epidermis and the epithelium of skin appendages (hair follicle,
  sebaceous gland, sweat duct). Dermis and subcutis contain none.**
  - Rule: `kind=region, types=[Keratinocytes], regions=[epidermis, hair follicle, sweat gland], mode=only,
    strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Keratinocytes], neighbor_types=[Keratinocytes, Melanocytes,
    Langerhans cells, Sweat gland cells], min_frac=0.4, mode=only, strength=weak`
  - Exceptions and caveats:
    - Squamous cell carcinoma, BCC and epidermal inclusion cysts put keratinocytes in the dermis.
    - Tangential sections show islands of epidermis that look dermal. Draw them as epidermis.
    - Follicular keratinocytes are not in the list and will type as Basal or Suprabasal keratinocytes.
  - Confidence: high. Strength: strong (normal).
  - Evidence: [1][3]

- **SKIN-02. Basal keratinocytes (KRT15, COL17A1, TP63, KRT5/KRT14) form the single bottom layer of the epidermis on
  the basement membrane, which is the layer in contact with the dermis. In hair follicles the KRT5/KRT14/KRT15+ outer
  root sheath and bulge are basal-like.**
  - Rule: `kind=region, types=[Basal keratinocytes], regions=[epidermis, hair follicle, sweat gland], mode=only,
    strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Basal keratinocytes], neighbor_types=[Fibroblasts, Endothelial,
    Pericytes], min_frac=0.6, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Psoriasis and wound edges expand the proliferative basal compartment upward (acanthosis).
    - Sweat-duct and myoepithelial cells express KRT5/KRT14 and will type as Basal keratinocytes; this is why the
      sweat-gland region is included.
  - Confidence: high. Strength: strong.
  - Evidence: [1][3]

- **SKIN-03. Suprabasal (spinous) keratinocytes (KRT1, KRT10, DSG1) sit above the basal layer. They do not touch the
  dermis directly.**
  - Rule: `kind=region, types=[Suprabasal keratinocytes], regions=[epidermis, hair follicle], mode=only,
    strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Suprabasal keratinocytes], neighbor_types=[Fibroblasts,
    Endothelial, Pericytes], min_frac=0.3, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - At 15 nearest neighbors, the first suprabasal row is often within reach of papillary-dermis cells in thin skin
      (face, eyelid), so keep the neighborhood rule weak.
    - Follicular infundibulum keratinocytes are suprabasal-like.
  - Confidence: high. Strength: strong (region), weak (neighborhood).
  - Evidence: [1][3]

- **SKIN-04. Granular keratinocytes (LOR, FLG, IVL) are the uppermost living layer, just under the stratum corneum.
  They are the farthest keratinocytes from the dermis and are a small population.**
  - Rule: `kind=region, types=[Granular keratinocytes], regions=[epidermis, hair follicle], mode=only, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Granular keratinocytes], neighbor_types=[Fibroblasts,
    Endothelial, Pericytes, Smooth muscle cells], min_frac=0.15, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Psoriasis shows parakeratosis with a reduced or absent granular layer; ichthyoses and atopic dermatitis
      (FLG-low) change marker levels.
    - Mucosal (non-keratinized) sites have no granular layer.
    - Granular cells are few in dissociated data [1], but they form a continuous layer on sections.
  - Confidence: high. Strength: strong (region), moderate (neighborhood).
  - Evidence: [1][2]

- **SKIN-05. Melanocytes sit in the basal layer of the epidermis, about one per 10 basal keratinocytes, and in the
  hair-bulb matrix of pigmented hair. They are not normally in the dermis.**
  - Rule: `kind=region, types=[Melanocytes], regions=[epidermis, hair follicle], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Melanocytes], neighbor_types=[Basal keratinocytes,
    Keratinocytes], min_frac=0.4, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Intradermal and compound nevi, blue nevi, dermal melanocytosis and melanoma place real melanocytes in the
      dermis; a strong rule would erase them.
    - Melanophages (dermal macrophages full of melanin) can carry PMEL/TYRP1 transcript spill-over.
  - Confidence: high (normal). Strength: moderate.
  - Evidence: [7] (basal-layer position, epidermal melanin unit); [1][8] (melanocyte populations in epidermal and
    follicular single-cell data)

- **SKIN-06. Langerhans cells (CD207, CD1A) are the dendritic cells resident in the epidermis, mainly suprabasal.**
  - Rule: `kind=region, types=[Langerhans cells], regions=[epidermis, hair follicle], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Langerhans cells], neighbor_types=[Keratinocytes], min_frac=0.4,
    mode=only, strength=weak`
  - Exceptions and caveats:
    - Langerin (CD207) is also expressed by a dermal DC subset [6], and migrating LCs pass through dermal lymphatics,
      so some CD207+ dermal cells are real.
    - Langerhans cell histiocytosis.
    - Psoriatic epidermis has fewer LCs [1].
  - Confidence: high. Strength: moderate.
  - Evidence: [1][6]

- **SKIN-07. Sweat gland cells (DCD, SCGB2A2, MUCL1, KRT19) form the eccrine secretory coils deep in the dermis or at
  the dermal-subcutaneous junction. They are absent from the epidermis, except the intraepidermal duct (acrosyringium),
  which is keratinocyte-like.**
  - Rule: `kind=region, types=[Sweat gland cells], regions=[sweat gland], mode=only, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Sweat gland cells], neighbor_types=[Suprabasal keratinocytes,
    Granular keratinocytes], min_frac=0.3, mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Apocrine glands (axilla, groin) and adnexal tumors.
    - DCD protein is secreted onto the skin surface, but its RNA is gland-restricted.
    - PIP/SCGB2A2 are also breast markers; watch for spill-over only within glands.
  - Confidence: high. Strength: strong.
  - Evidence: [3] (KRT19 and DCD mark sweat glands; PTGDS+ fibroblasts near eccrine coils)

- **SKIN-08. Fibroblasts and vascular cells (Endothelial, Pericytes, Lymphatic EC) live in the dermis and subcutis.
  The epidermis is avascular and contains no fibroblasts.**
  - Rule: `kind=region, types=[Fibroblasts, Endothelial, Pericytes, Smooth muscle cells], regions=[epidermis],
    mode=exclude, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Fibroblasts, Endothelial], neighbor_types=[Suprabasal
    keratinocytes, Granular keratinocytes], min_frac=0.5, mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Dermal papillae push capillary loops and fibroblasts close to the epidermis. A hand-drawn epidermis border
      captures some of them, so use moderate, not strong.
  - Confidence: high. Strength: moderate.
  - Evidence: [3][4]; avascular epidermis: (textbook histology, not verified online)

- **SKIN-09. Dermal fibroblasts are zonated: PTGDS+ (papillary) fibroblasts lie in the papillary dermis and around
  follicles and eccrine coils, SFRP2+ fibroblasts lie between collagen bundles concentrated in the reticular dermis,
  and POSTN+ fibroblasts surround the lower hair follicle.**
  - Rule: none for the built-in Fibroblasts type, which merges these. If the list is split:
    `kind=region, types=[<PTGDS+ fibroblasts>], regions=[papillary dermis], mode=enrich, strength=weak` and
    `kind=region, types=[<SFRP2+ fibroblasts>], regions=[reticular dermis], mode=enrich, strength=weak`.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Aging blurs fibroblast identity ("loss of priming") [5].
    - BCC expands the POSTN+ population [3].
    - Scleroderma and keloid alter the states.
  - Confidence: medium. Strength: weak.
  - Evidence: [3][4][5]

- **SKIN-10. Smooth muscle cells in skin are restricted to larger vessel walls and the arrector pili muscles attached
  to hair follicles (plus the dartos and areolar muscle at special sites).**
  - Rule: `kind=region, types=[Smooth muscle cells], regions=[papillary dermis, epidermis], mode=exclude,
    strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Smooth muscle cells], neighbor_types=[Endothelial, Pericytes,
    Smooth muscle cells, Basal keratinocytes], min_frac=0.3, mode=only, strength=weak` (arrector pili sit next to
    the follicle's basal-like outer root sheath).
  - Exceptions and caveats:
    - Leiomyomas.
    - Myoepithelial cells of sweat coils express ACTA2/MYH11 and can be typed as Smooth muscle cells inside the
      `sweat gland` region.
  - Confidence: high. Strength: moderate.
  - Evidence: [3] (SMC restricted to larger vessels and the arrector pili muscle)

- **SKIN-11. Adipocytes are in the subcutis and in the deep dermis around large (terminal or anagen) hair follicles.
  They are absent from the epidermis and papillary dermis.**
  - Rule: `kind=region, types=[Adipocytes], regions=[epidermis, papillary dermis], mode=exclude, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Adipocytes], neighbor_types=[Keratinocytes], min_frac=0.3,
    mode=exclude, strength=moderate`
  - Exceptions and caveats:
    - Nevus lipomatosus.
    - Most biopsies (shave and punch) contain little or no subcutis, so few adipocytes are expected at all.
    - Adipocyte transcripts are sparse and spill into dermal fibroblasts.
  - Confidence: high. Strength: strong.
  - Evidence: (textbook histology, not verified online)

- **SKIN-12. Mast cells and most T cells, macrophages and dendritic cells of normal skin are dermal, clustered around
  the superficial vascular plexus. The epidermis normally holds only Langerhans cells and a few resident T cells.**
  - Rule: `kind=region, types=[Mast cells, Macrophages, Plasma cells, B cells], regions=[epidermis], mode=exclude,
    strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Mast cells], neighbor_types=[Endothelial, Pericytes,
    Fibroblasts], min_frac=0.3, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - Do not exclude T cells from the epidermis: CD8 tissue-resident memory T cells are intraepidermal.
    - Psoriasis, atopic dermatitis, lichenoid dermatitis and cutaneous T-cell lymphoma (epidermotropism) put many
      immune cells into the epidermis [2]. Neutrophil microabscesses in psoriasis are epidermal.
  - Confidence: medium-high (normal). Strength: moderate (normal), off (inflamed).
  - Evidence: [2]; perivascular dermal location: (textbook histology, not verified online)

- **SKIN-13. Hair follicles carry their own melanocytes (bulb), basal-like keratinocytes (outer root sheath, bulge)
  and specialized mesenchyme (dermal papilla, dermal sheath). Follicle-associated cells must not be read as invasion
  of the dermis.**
  - Rule: `kind=region, types=[Basal keratinocytes, Melanocytes], regions=[hair follicle], mode=enrich,
    strength=weak`. Mainly this is a reminder to draw follicles as their own region so that SKIN-01 and SKIN-05 do
    not suppress them.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Follicle density differs widely by site: abundant in scalp, sparse in trunk, nearly absent in foreskin [1].
    - BCC arises from follicle-like cells [3].
  - Confidence: high. Strength: weak.
  - Evidence: [1][3][8]

**References**
[1] Cheng JB et al., Transcriptional programming of normal and inflamed human epidermis at single-cell resolution,
Cell Reports, 2018. doi:10.1016/j.celrep.2018.09.006 (PMID 30355494)
[2] Reynolds G et al., Developmental cell programs are co-opted in inflammatory skin disease, Science, 2021.
doi:10.1126/science.aba6500 (PMID 33479125)
[3] Ganier C et al., Multiscale spatial mapping of cell populations across anatomical sites in healthy human skin and
basal cell carcinoma, PNAS, 2024. doi:10.1073/pnas.2313326120 (PMID 38165934)
[4] Tabib T et al., SFRP2/DPP4 and FMO1/LSP1 define major fibroblast populations in human skin, Journal of
Investigative Dermatology, 2018. doi:10.1016/j.jid.2017.09.045 (PMID 29080679)
[5] Solé-Boldo L et al., Single-cell transcriptomes of the human skin reveal age-related loss of fibroblast priming,
Communications Biology, 2020. doi:10.1038/s42003-020-0922-4 (PMID 32327715)
[6] Merad M et al., Origin, homeostasis and function of Langerhans cells and other langerin-expressing dendritic
cells, Nature Reviews Immunology, 2008. doi:10.1038/nri2455 (PMID 19029989)
[7] Cichorek M et al., Skin melanocytes: biology and development, Postepy Dermatologii i Alergologii, 2013.
doi:10.5114/pdia.2013.33376 (PMID 24278043)
[8] Takahashi R et al., Defining transcriptional signatures of human hair follicle cell states, Journal of
Investigative Dermatology, 2020. doi:10.1016/j.jid.2019.07.726

---

## Brain (cortex)

Type names are from `brain.txt`. That file does not include `@immune` or `@stroma`, so the brain types are Neurons
(Excitatory neurons; Inhibitory neurons: PVALB / SST / VIP / LAMP5 interneurons), Astrocytes, Oligodendrocytes, OPCs,
Microglia, Border macrophages, T cells, Endothelial cells, Pericytes and VLMCs. The built-in list does not split
excitatory neurons by layer.

**Compartments to draw**
- `meninges`: the thin fibrous pia-arachnoid layer on the outer cortical surface, with its large vessels. It often
  tears or folds away in sections.
- `L1`: the cell-sparse band, about 100-250 um thick in human and about 100 um in mouse, just below the pia. It has
  pale neuropil and few small nuclei.
- `L2/3`: densely packed small-to-medium pyramidal somata below L1.
- `L4`: a granular band of small round somata, prominent in sensory cortex (e.g. V1, S1). It is absent or very thin in
  motor and frontal agranular cortex, so do not draw it there.
- `L5`: large pyramidal somata, with Betz cells in human M1, and lower cell density than L2/3.
- `L6`: the deepest grey-matter layer, with dense small fusiform and pyramidal somata. Its lower part (L6b/subplate)
  borders the white matter.
- `white matter`: pale, eosinophilic fibre tracts with rows of small round dark nuclei (interfascicular
  oligodendrocytes) and almost no large somata. In human sections it is subcortical white matter; in mouse it is the
  corpus callosum or external capsule.
- On H&E, layer borders are gradual. Draw them generously and prefer moderate or weak strength for rules that use
  layer regions.

The `kind=depth` rule (epithelial-marker density) is not meaningful in brain cortex because there is no epithelium.
The only epithelium is the choroid plexus, which is not cortex. Use regions or neighborhood rules instead.

**Priors**

- **BRN-01. Neuronal somata are rare in white matter. A few interstitial neurons (mostly GABAergic, densest just
  below the cortex) are normal, especially in human.**
  - Rule: `kind=region, types=[Neurons], regions=[white matter], mode=exclude, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Neurons], neighbor_types=[Oligodendrocytes], min_frac=0.6,
    mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Interstitial white-matter neurons are far more numerous in primates than in rodents [9], so do not use
      "strong" on human tissue.
    - Neuronal transcripts (SNAP25, SYT1) carried in axons can spill into segmented white-matter cells.
    - Heterotopia, cortical dysplasia and gliomas break this rule.
  - Confidence: high. Strength: moderate.
  - Evidence: [3][9]

- **BRN-02. Oligodendrocytes are densest in white matter. In grey matter they are fewer and increase with depth
  (L5-L6 more than L2/3, and fewest in L1).**
  - Rule: `kind=region, types=[Oligodendrocytes], regions=[white matter], mode=enrich, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Oligodendrocytes], neighbor_types=[Oligodendrocytes],
    min_frac=0.5, mode=enrich, strength=weak` (self-reinforcing: use only after a first pass)
  - Exceptions and caveats:
    - Satellite oligodendrocytes beside neuronal somata are normal in grey matter, so never use `only`.
    - Demyelinating lesions (multiple sclerosis) remove oligodendrocytes from white-matter plaques.
  - Confidence: high. Strength: moderate.
  - Evidence: [3][4]

- **BRN-03. L1 contains no pyramidal (excitatory) somata. Its few neurons are almost all GABAergic interneurons.**
  - Rule: `kind=region, types=[Excitatory neurons], regions=[L1], mode=exclude, strength=moderate`
  - Without drawn regions: none (no surface marker; VLMC density at the pia is too sparse to stand in for "L1").
  - Exceptions and caveats:
    - Cajal-Retzius cells are glutamatergic, sit in L1 and are sparse in adult cortex. This rule would type them as
      Excitatory neurons and suppress them, but they are a real finding, mostly in developing or young tissue.
    - Upper-L2 pyramidal somata fall into a loosely drawn L1, which is why the strength is moderate rather than strong.
  - Confidence: high. Strength: moderate.
  - Evidence: [1][5][6]

- **BRN-04. L1 interneurons are mainly LAMP5 types (neurogliaform and canopy cells) plus some VIP cells. PVALB and
  SST interneuron somata are scarce in L1.**
  - Rules:
    - `kind=region, types=[LAMP5 interneurons], regions=[L1], mode=enrich, strength=moderate`
    - `kind=region, types=[PVALB interneurons, SST interneurons], regions=[L1], mode=exclude, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - LAMP5 interneurons also occur in L2/3, and some in deeper layers, so never use `only`.
    - Human L1 is more diverse than mouse L1. For example, the human rosehip cell has a marker signature not seen in
      mouse [7]. Hodge et al. report shifts in laminar distribution between human and mouse homologous types [2].
  - Confidence: medium-high (mouse); medium (human). Strength: moderate for enrich, weak for exclude.
  - Evidence: [1][2][5][6][7]

- **BRN-05. Each excitatory neuron class has a layer of residence:**
  - L2/3 IT (CUX2, LAMP5-low) in L2/3
  - L4/5 IT (RORB) in L4 and upper L5
  - L5 ET/PT (FEZF2, BCL11B, POU3F1) in L5
  - L5/6 NP in L5-L6
  - L6 CT (FOXP2, SYT6) in L6
  - L6b (CCN2/CTGF, NXPH4) at the L6-white-matter border

  None of these subtypes is in the built-in list.
  - Rule (if the user adds the subtypes): `kind=region, types=[L5 ET neurons (not in the built-in list)],
    regions=[L5], mode=only, strength=moderate`, and similarly for the other subtypes, each with its own layer and
    the adjacent layer(s) listed in `regions`.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - In mouse, IT neurons form a continuous gradient with cortical depth [4].
    - In human, excitatory types are less layer-restricted than expected and can span adjacent layers [2].
    - Always include the neighbouring layer(s) in `regions`, and do not use "strong".
    - Agranular cortex (motor, frontal) has no distinct L4.
  - Confidence: medium. Strength: moderate (weak for human).
  - Evidence: [1][2][3][4]

- **BRN-06. VIP interneurons are concentrated in superficial layers (L1-L2/3). SST and PVALB interneurons are spread
  over L2-L6 and are densest in L4-L5.**
  - Rule: `kind=region, types=[VIP interneurons], regions=[L1, L2/3], mode=enrich, strength=weak`
  - Without drawn regions: none.
  - Exceptions and caveats: These are density gradients, not boundaries. Every interneuron class is found in most
    layers.
  - Confidence: medium. Strength: weak.
  - Evidence: [1][4]

- **BRN-07. Astrocytes are present in every layer and in white matter. Glia limitans and fibrous astrocytes make the
  pial surface and white matter rich in GFAP-high astrocytes, and astrocytes themselves form superficial, middle and
  deep layers that do not match the neuronal layers.**
  - Rule: `kind=region, types=[Astrocytes], regions=[L1], mode=enrich, strength=weak` (optional)
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Do not write any `only` or `exclude` rule for astrocytes in cortex.
    - Interlaminar astrocytes (somata in L1, processes into L2-L4) are primate-specific [8].
    - Reactive gliosis near lesions, tumors and injection tracks increases astrocyte density anywhere.
  - Confidence: high (that the cells are ubiquitous); low for the enrichment. Strength: weak.
  - Evidence: [8][16]

- **BRN-08. OPCs and microglia are tiled roughly evenly through grey and white matter. Do not use region rules for
  them.**
  - Rule: none (warning prior).
  - Exceptions and caveats:
    - Microglia accumulate around plaques, infarcts and tumors (disease-associated microglia), so do not use density
      rules there either.
    - A `Microglia exclude in meninges` rule is reasonable. See BRN-10.
  - Confidence: high. Strength: none.
  - Evidence: [11][12]

- **BRN-09. Pericytes sit on capillaries. Every pericyte soma touches an endothelial tube.**
  - Rule: `kind=neighborhood, types=[Pericytes], neighbor_types=[Endothelial cells], min_frac=0.07 (>=1 of 15),
    mode=only, strength=moderate`
  - Without drawn regions: this rule is already region-free.
  - Exceptions and caveats:
    - Brain has the highest pericyte coverage of any organ, so this rule is safest here.
    - Arteriolar smooth-muscle-like mural cells (ACTA2-high) are also adjacent to endothelium and are typed as
      Pericytes in this list.
    - If the panel types endothelium poorly (few CLDN5/FLT1 counts), the rule wrongly suppresses real pericytes.
      Check the endothelial call rate first.
  - Confidence: high. Strength: moderate.
  - Evidence: [10][17]

- **BRN-10. Border-associated macrophages (MRC1/LYVE1/CD163) live at the borders: meninges, perivascular spaces of
  penetrating vessels, and choroid plexus. Microglia are in the parenchyma and not in the meninges.**
  - Rules:
    - `kind=region, types=[Border macrophages], regions=[meninges], mode=enrich, strength=strong`
    - `kind=neighborhood, types=[Border macrophages], neighbor_types=[Endothelial cells, Pericytes, VLMCs],
      min_frac=0.13 (>=2 of 15), mode=only, strength=moderate`
    - `kind=region, types=[Microglia], regions=[meninges], mode=exclude, strength=weak`
  - Exceptions and caveats:
    - Blood-borne monocyte-derived macrophages enter the parenchyma in stroke, tumors and multiple sclerosis lesions.
      These cells express MRC1/CD163 away from vessels.
    - Activated microglia can up-regulate CD163.
  - Confidence: high. Strength: as given.
  - Evidence: [11]

- **BRN-11. VLMCs (meningeal and perivascular fibroblasts) are in the pia-arachnoid and around larger penetrating
  arterioles and venules, never free in the neuropil. Perivascular fibroblasts are absent from capillaries.**
  - Rules:
    - `kind=region, types=[VLMCs], regions=[meninges], mode=enrich, strength=strong`
    - `kind=neighborhood, types=[VLMCs], neighbor_types=[Endothelial cells, Pericytes, VLMCs, Border macrophages],
      min_frac=0.2, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Fibrotic scars after injury and the stroma of meningioma or metastases break this rule.
    - VLMCs are PDGFRA+, so OPC/VLMC confusion is possible. Check DCN/COL1A2.
  - Confidence: high. Strength: as given.
  - Evidence: [10][12]

- **BRN-12. T cells are rare in normal brain. The few present are tissue-resident memory cells in perivascular
  (Virchow-Robin) spaces and meninges, not scattered in the neuropil.**
  - Rule: `kind=neighborhood, types=[T cells], neighbor_types=[Endothelial cells, Pericytes, VLMCs, Border
    macrophages], min_frac=0.13, mode=only, strength=weak`
  - Exceptions and caveats:
    - Do not use this rule in encephalitis, multiple sclerosis, gliomas or metastases, or after stroke, where
      parenchymal T cells are the finding.
    - A strong rule here would hide real pathology.
  - Confidence: medium (normal tissue). Strength: weak.
  - Evidence: [13]

- **BRN-13. Human vs mouse layout (drawing guidance, no rule).**
  - Human cortex is 2-4 mm thick with an expanded L2/3, and a section usually includes a gyrus and subcortical white
    matter. Mouse cortex is about 1-2 mm thick and its white matter is the corpus callosum.
  - Cell-type proportions and laminar distributions differ between species for homologous types [2][14].
  - Do not reuse mouse layer thickness, or mouse rules about which type is "absent from layer X", for human without
    checking.
  - Confidence: high. Evidence: [2][14][15]

**References**
- [1] Tasic B et al., Shared and distinct transcriptomic cell types across neocortical areas, Nature, 2018.
  doi:10.1038/s41586-018-0654-5. PMID 30382198.
- [2] Hodge RD et al., Conserved cell types with divergent features in human versus mouse cortex, Nature, 2019.
  doi:10.1038/s41586-019-1506-7. PMID 31435019.
- [3] Maynard KR et al., Transcriptome-scale spatial gene expression in the human dorsolateral prefrontal cortex, Nature
  Neuroscience, 2021. doi:10.1038/s41593-020-00787-0. PMID 33558695.
- [4] Zhang M et al., Spatially resolved cell atlas of the mouse primary motor cortex by MERFISH, Nature, 2021.
  doi:10.1038/s41586-021-03705-x.
- [5] Schuman B et al., Four unique interneuron populations reside in neocortical layer 1, Journal of Neuroscience,
  2019. doi:10.1523/JNEUROSCI.1613-18.2018. PMID 30413647.
- [6] Yao Z et al., A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain, Nature,
  2023. doi:10.1038/s41586-023-06812-z.
- [7] Boldog E et al., Transcriptomic and morphophysiological evidence for a specialized human cortical GABAergic cell
  type, Nature Neuroscience, 2018. doi:10.1038/s41593-018-0205-2.
- [8] Oberheim NA et al., Uniquely hominid features of adult human astrocytes, Journal of Neuroscience, 2009.
  doi:10.1523/JNEUROSCI.4707-08.2009.
- [9] Suarez-Sola ML et al., Neurons in the white matter of the adult human neocortex, Frontiers in Neuroanatomy, 2009.
  doi:10.3389/neuro.05.007.2009.
- [10] Vanlandewijck M et al., A molecular atlas of cell types and zonation in the brain vasculature, Nature, 2018.
  doi:10.1038/nature25739. PMID 29443965.
- [11] Van Hove H et al., A single-cell atlas of mouse brain macrophages reveals unique transcriptional identities
  shaped by ontogeny and tissue environment, Nature Neuroscience, 2019. doi:10.1038/s41593-019-0393-4. PMID 31061494.
- [12] Marques S et al., Oligodendrocyte heterogeneity in the mouse juvenile and adult central nervous system, Science,
  2016. doi:10.1126/science.aaf6463.
- [13] Smolders J et al., Tissue-resident memory T cells populate the human brain, Nature Communications, 2018.
  doi:10.1038/s41467-018-07053-9.
- [14] Jorstad NL et al., Comparative transcriptomics reveals human-specific cortical features, Science, 2023.
  doi:10.1126/science.ade9516.
- [15] Cortical thickness and agranular L4 (textbook histology, not verified online).
- [16] Bayraktar OA et al., Astrocyte layers in the mammalian cerebral cortex revealed by a single-cell in situ
  transcriptomic map, Nature Neuroscience, 2020. doi:10.1038/s41593-020-0602-1.
- [17] Armulik A, Genove G, Betsholtz C, Pericytes: developmental, physiological, and pathological perspectives,
  problems, and promises, Developmental Cell, 2011. doi:10.1016/j.devcel.2011.07.001.

---

## Lymph node and tonsil

Type names are from `lymph_node.txt`. It has its own immune and stromal lines and does not include
`@immune`/`@stroma`. The types are:
- T and NK cells: CD4 T cells (Naive CD4, Tfh, Tregs, Effector/memory CD4), CD8 T cells, NK cells
- B and plasma cells: Naive B, Germinal center B, Memory B, Plasma cells
- Myeloid: Macrophages, Dendritic cells, pDC
- Stromal: FDCs, FRCs
- HEV / endothelial
- Lymphatic EC

There is no epithelial type, so tonsil epithelium needs an added line: "Epithelial cells (not in the built-in list)".

**Compartments to draw**
- `capsule/subcapsular sinus` (LN): the thin fibrous capsule plus the pale, loosely cellular sinus just under it, which
  holds lymph, macrophages and lining cells.
- `follicle`: round, dark-blue aggregates of small lymphocytes in the outer cortex. A primary follicle, or the mantle
  zone around a germinal centre, is uniformly dark.
- `germinal center`: the pale centre of a secondary follicle, with larger cells, mitoses and tingible-body macrophages
  (apoptotic debris). Optionally split it into:
  - `dark zone`: centroblast-rich, many mitoses. In LN it faces the paracortex/medulla; in tonsil it faces away from
    the crypt.
  - `light zone`: paler, with centrocytes and the FDC meshwork. It faces the capsule, or the crypt epithelium in
    tonsil.
- `paracortex`: the diffuse, non-nodular T-cell zone between follicles and medulla, with cuboidal-endothelium high
  endothelial venules (HEVs). In tonsil this compartment is the `interfollicular area`.
- `medulla` (LN): cords of lymphocytes and plasma cells separated by pale medullary sinuses, next to the hilum.
- Tonsil only:
  - `surface epithelium`: non-keratinized stratified squamous epithelium.
  - `crypt epithelium`: the invaginated lining of the crypts, reticulated and patchy, and heavily infiltrated by
    lymphocytes.

The `kind=depth` rule (epithelial-marker density) is meaningless in lymph node. In tonsil it measures distance to
crypt or surface epithelium only if the panel has keratin markers (e.g. KRT5, KRT13), and even then it does not
separate follicles from the interfollicular zone. Prefer regions.

**Priors**

- **LN-01. Germinal center B cells are confined to germinal centres of secondary follicles.**
  - Rule: `kind=region, types=[Germinal center B cells], regions=[germinal center], mode=only, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Germinal center B cells], neighbor_types=[Germinal center B
    cells, FDCs, Tfh cells], min_frac=0.4, mode=only, strength=moderate`
  - Exceptions and caveats:
    - Follicular lymphoma, Burkitt lymphoma and GCB-type DLBCL spread GC-like cells everywhere. Turn the rule off for
      lymphoma.
    - Small or tangentially cut GCs are easy to miss on H&E. Check for MKI67/AICDA hot spots before drawing.
    - Mouse lymph nodes from unimmunized specific-pathogen-free mice have few or no GCs.
  - Confidence: high. Strength: strong (normal or reactive tissue).
  - Evidence: [1][2][11]

- **LN-02. FDCs form the meshwork of B-cell follicles and germinal centres, densest in the light zone, and are absent
  from the paracortex and medulla.**
  - Rules:
    - `kind=region, types=[FDCs], regions=[follicle, germinal center], mode=only, strength=strong`
    - `kind=region, types=[FDCs], regions=[dark zone], mode=exclude, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[FDCs], neighbor_types=[B and plasma cells, Tfh cells],
    min_frac=0.5, mode=only, strength=moderate`
  - Exceptions and caveats:
    - FDC dendrites are long. In segmentation, FDCSP/CR2 transcripts spill into neighbouring B cells, so B cells may
      be called FDCs, not the other way round.
    - FDC sarcoma and angioimmunoblastic T-cell lymphoma (an expanded FDC meshwork outside follicles) break the rule.
  - Confidence: high. Strength: as given.
  - Evidence: [1][3][4]

- **LN-03. GC B cells are polarized within the germinal centre. CXCR4-high proliferating centroblasts sit in the dark
  zone. Centrocytes, FDCs and Tfh cells sit in the light zone.**
  - Rules:
    - `kind=region, types=[Tfh cells], regions=[light zone], mode=enrich, strength=weak`
    - FDCs in the dark zone: see LN-02.
  - Without drawn regions: none.
  - Exceptions and caveats:
    - The dark and light zones are hard to separate on H&E. Draw them only when the polarity is obvious.
    - GC B cells are in both zones, so do not restrict GC B cells by zone.
  - Confidence: high (biology); medium (drawability). Strength: weak.
  - Evidence: [2][3]

- **LN-04. Naive B cells fill primary follicles and the mantle zone around germinal centres. They are not a
  paracortex or medulla population.**
  - Rule: `kind=region, types=[Naive B cells], regions=[follicle, germinal center], mode=only, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[Naive B cells], neighbor_types=[B and plasma cells, FDCs],
    min_frac=0.4, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - Recirculating B cells cross the paracortex after entering through HEVs, and they sit at the T-B border.
    - B-cell lymphomas (CLL/SLL, mantle cell) efface the architecture.
  - Confidence: high. Strength: moderate.
  - Evidence: [1][4]

- **LN-05. Tfh cells are in follicles and germinal centres, and at the T-B border.**
  - Rule: `kind=region, types=[Tfh cells], regions=[follicle, germinal center], mode=enrich, strength=strong`
  - Without drawn regions: `kind=neighborhood, types=[Tfh cells], neighbor_types=[B and plasma cells, FDCs],
    min_frac=0.4, mode=enrich, strength=moderate`
  - Exceptions and caveats:
    - Pre-Tfh and extrafollicular helper cells sit at the T-B border and interfollicular area, so use `enrich`, not
      `only`.
    - Angioimmunoblastic T-cell lymphoma is a Tfh-derived lymphoma that spreads into the paracortex.
  - Confidence: high. Strength: strong (enrich).
  - Evidence: [5][2]

- **LN-06. Naive CD4 T cells, CD8 T cells, CCR7+/LAMP3+ dendritic cells, FRCs (CCL19/CCL21 T-zone reticular cells)
  and HEVs make up the paracortex.**
  - Rules:
    - `kind=region, types=[Naive CD4 T cells, CD8 T cells, FRCs, Dendritic cells], regions=[paracortex,
      interfollicular area], mode=enrich, strength=moderate`
    - `kind=region, types=[Naive CD4 T cells], regions=[germinal center], mode=exclude, strength=moderate`
  - Without drawn regions: `kind=neighborhood, types=[FRCs], neighbor_types=[T and NK cells, Dendritic cells],
    min_frac=0.5, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - CD8 T cells also occur in the medulla and in germinal centres in small numbers.
    - Marginal reticular cells and the perivascular stroma of the medulla share CCL19/PDPN, so FRC calls outside
      the paracortex can be real stroma of another subset.
    - In tonsil the T zone is the interfollicular area.
  - Confidence: high. Strength: moderate.
  - Evidence: [1][4][6]

- **LN-07. HEVs are post-capillary venules of the paracortex and interfollicular zone. They are absent from follicle
  centres and from the medulla, where ordinary venules sit.**
  - Rule: not expressible as written, because `HEV / endothelial` covers all blood endothelium. If the user splits off
    an `HEV` subtype (CHST4, ACKR1, MECA79 epitope), use `kind=region, types=[HEV (not in the built-in list)],
    regions=[paracortex, interfollicular area], mode=only, strength=moderate`.
  - Exceptions and caveats:
    - HEV-like vessels appear in chronic inflammation and in tertiary lymphoid structures (see GEN-08).
    - GLYCAM1 is a pseudogene in human [12], so it is a marker in mouse only.
  - Confidence: high (biology). Strength: moderate.
  - Evidence: [6][12]

- **LN-08. Plasma cells collect in the medullary cords (LN) and in interfollicular and subepithelial areas (tonsil).
  They are rare in germinal centres.**
  - Rules:
    - `kind=region, types=[Plasma cells], regions=[medulla], mode=enrich, strength=strong`
    - `kind=region, types=[Plasma cells], regions=[germinal center], mode=exclude, strength=weak`
  - Without drawn regions: none (the medullary cords have no simple neighbor signature).
  - Exceptions and caveats:
    - Extrafollicular plasmablast foci in the paracortex are normal in active responses.
    - Plasma-cell neoplasms and Castleman disease break this rule.
    - JCHAIN/MZB1 transcripts spill into neighbours, giving false plasma calls at follicle edges.
  - Confidence: medium-high. Strength: as given.
  - Evidence: [7][8]

- **LN-09. Lymph-node macrophages are in three niches: the subcapsular sinus (a sheet over follicles), the medullary
  sinuses and the medullary cords. Germinal centres also contain tingible-body macrophages.**
  - Rule: `kind=region, types=[Macrophages], regions=[capsule/subcapsular sinus, medulla], mode=enrich,
    strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - Do not exclude macrophages from germinal centres, because tingible-body macrophages are a normal GC feature.
    - In human, the subcapsular sinus macrophage layer is thinner than in rodents, where it is dense.
    - Sinus histiocytosis and granulomas (e.g. sarcoid, TB) expand macrophages anywhere.
  - Confidence: high. Strength: moderate.
  - Evidence: [8]

- **LN-10. Lymphatic EC line the subcapsular, cortical and medullary sinuses and are absent from follicles, germinal
  centres and the T-zone parenchyma.**
  - Rule: `kind=region, types=[Lymphatic EC], regions=[germinal center, follicle], mode=exclude, strength=moderate`
  - Without drawn regions: none.
  - Exceptions and caveats:
    - The subcapsular sinus floor lies directly on follicles, so a loosely drawn follicle catches sinus LECs. Keep
      `moderate`.
    - Cortical and transverse sinuses cross the paracortex, so do not exclude LECs there.
    - LEC subsets differ by sinus [9].
  - Confidence: high. Strength: moderate.
  - Evidence: [9][6]

- **LN-11. Tonsil: the crypt epithelium is reticulated and heavily infiltrated by B cells, T cells, dendritic cells
  and macrophages. Lymphocytes inside the crypt epithelium are normal, not a segmentation error.**
  - Rule: no `exclude` for immune types in `crypt epithelium`. If the user adds Epithelial cells, use
    `kind=region, types=[Epithelial cells (not in the built-in list)], regions=[surface epithelium, crypt
    epithelium], mode=only, strength=strong`.
  - Without drawn regions: `kind=depth` (if the panel has keratins) is useful only for the epithelial type.
  - Exceptions and caveats:
    - Memory B cells are enriched in and beneath the crypt epithelium. That observation is not verified online here,
      so no rule is given.
    - Tonsil GCs face the crypt with their light zone.
  - Confidence: high. Strength: strong (epithelial `only`).
  - Evidence: [10][11]

- **LN-12. Lymphoid-tissue species and sample notes (drawing guidance).**
  - Human tonsils and reactive nodes from children have many secondary follicles. Adult and elderly nodes have fewer,
    plus fat and fibrosis in the hilum.
  - Mouse LN (no tonsil) is small, with a thin cortex.
  - Spleen (also in this list) has white pulp (follicles, PALS) and red pulp instead of cortex and medulla. Draw
    `white pulp`/`red pulp` and do not reuse the LN medulla rules.
  - Confidence: medium. Evidence: [11] (tonsil age effects); spleen layout: textbook histology, not verified online.

**References**
- [1] Rodda LB et al., Single-cell RNA sequencing of lymph node stromal cells reveals niche-associated heterogeneity,
  Immunity, 2018. doi:10.1016/j.immuni.2018.04.006. PMID 29752062.
- [2] Victora GD, Nussenzweig MC, Germinal centers, Annual Review of Immunology, 2022.
  doi:10.1146/annurev-immunol-120419-022408.
- [3] Allen CD et al., Germinal center dark and light zone organization is mediated by CXCR4 and CXCR5, Nature
  Immunology, 2004. doi:10.1038/ni1100. PMID 15300245.
- [4] Krishnamurty AT, Turley SJ, Lymph node stromal cells: cartographers of the immune system, Nature Immunology,
  2020. doi:10.1038/s41590-020-0635-3.
- [5] Crotty S, T follicular helper cell biology: a decade of discovery and diseases, Immunity, 2019.
  doi:10.1016/j.immuni.2019.04.011. PMID 31117010.
- [6] Girard JP, Moussion C, Forster R, HEVs, lymphatics and homeostatic immune cell trafficking in lymph nodes, Nature
  Reviews Immunology, 2012. doi:10.1038/nri3298.
- [7] Kleshchevnikov V et al., Cell2location maps fine-grained cell types in spatial transcriptomics, Nature
  Biotechnology, 2022. doi:10.1038/s41587-021-01139-4. (Human LN Visium compartment mapping.)
- [8] Gray EE, Cyster JG, Lymph node macrophages, Journal of Innate Immunity, 2012. doi:10.1159/000337007. PMID
  22488251.
- [9] Takeda A et al., Single-cell survey of human lymphatics unveils marked endothelial cell heterogeneity and
  mechanisms of homing for neutrophils, Immunity, 2019. doi:10.1016/j.immuni.2019.06.027. PMID 31402260.
- [10] Perry ME, The specialised structure of crypt epithelium in the human palatine tonsil and its functional
  significance, Journal of Anatomy, 1994. PMC1166820. (See also Perry M, Whyte A, Immunology of the tonsils,
  Immunology Today, 1998, PMID 9745205.)
- [11] Massoni-Badosa R et al., An atlas of cells in the human tonsil, Immunity, 2024.
  doi:10.1016/j.immuni.2024.01.006.
- [12] GeneCards, GLYCAM1 (Glycosylation Dependent Cell Adhesion Molecule 1, pseudogene), gene entry
  (genecards.org/card/GLYCAM1).

---

## General cross-organ priors (types of `immune.txt` and `stroma.txt`)

These apply in any organ whose list includes `@immune` and/or `@stroma`. Most are neighborhood rules, so they work
without drawn regions. Neighbor types are themselves assigned by expression only, so mis-typed neighbours propagate.
Keep these rules at moderate or weak strength and re-check after a first pass. With K = 15, `min_frac=0.07` means at
least 1 neighbour and `0.13` means at least 2.

**Compartments to draw**: these priors are organ-independent. Optionally draw `tertiary lymphoid structure` (a dense,
dark-blue nodular lymphoid aggregate, sometimes with a pale germinal centre, inside a non-lymphoid organ or tumor) and
`large vessel wall` (concentric smooth-muscle media around an arterial or venous lumen).

**Priors**

- **GEN-01. Pericytes are always adjacent to capillary or post-capillary endothelium.**
  - Rule: `kind=neighborhood, types=[Pericytes], neighbor_types=[Endothelial], min_frac=0.07, mode=only,
    strength=moderate`
  - Exceptions and caveats:
    - Pericyte coverage varies by organ: highest in brain and retina, low in skeletal muscle and some glands.
    - In tumors, pericytes detach and become loosely associated with vessels. Use `weak` there.
    - Hepatic stellate cells and kidney mesangial cells share RGS5/PDGFRB. Stellate cells are perisinusoidal, which
      is fine. Mesangial cells sit inside glomeruli with endothelium, which is also fine.
    - Fibroblasts are PDGFRA+ while mural cells are not [4]. Use PDGFRA to separate pericytes from fibroblasts, not
      location.
  - Confidence: high. Strength: moderate.
  - Evidence: [1][2][4]

- **GEN-02. Vascular smooth muscle forms the media of arteries, arterioles and veins, so vascular SMC are next to
  endothelium. Visceral smooth muscle (gut muscularis, airway, bladder, uterus, myometrium) forms thick layers with
  little endothelium.**
  - Rule: none as a hard rule. Optional: `kind=neighborhood, types=[Smooth muscle cells],
    neighbor_types=[Endothelial, Smooth muscle cells, Pericytes], min_frac=0.3, mode=enrich, strength=weak`
  - Warning:
    - Do not write `Smooth muscle cells only near Endothelial`, because it would wipe out the muscularis propria.
    - Do not exclude SMC from epithelial regions in breast, salivary, sweat and prostate glands. Myoepithelial and
      basal cells there express MYH11/ACTA2/CNN1 and will be typed as SMC if the list lacks a myoepithelial type.
    - Myofibroblasts (ACTA2+) in wounds and tumors also mimic SMC.
  - Confidence: high (warning). Strength: weak or none.
  - Evidence: [4] (textbook histology for the visceral layers)

- **GEN-03. Lymphatic EC vs blood EC.**
  - Initial lymphatics have no pericyte or SMC coat, so a cell typed Lymphatic EC whose neighbours are mostly
    pericytes is more likely a blood capillary EC.
  - The CNS parenchyma has no lymphatics. Meningeal lymphatics are in the dura.
  - Lymphatics in mucosae and skin lie in the lamina propria/submucosa or dermis, not in the epithelium.
  - Rule: `kind=neighborhood, types=[Lymphatic EC], neighbor_types=[Pericytes], min_frac=0.2, mode=exclude,
    strength=weak`
  - Warning (marker overlap, not location):
    - LYVE1 is expressed by perivascular/tissue macrophages [5], liver sinusoidal EC and splenic sinusoids.
    - PROX1 is expressed by hepatocytes (textbook, not verified online).
    - CCL21 and PDPN are expressed by fibroblastic reticular cells and other stroma.
    - Do not treat any of these alone as lymphatic identity. In liver, apply no Lymphatic EC prior outside portal
      tracts.
  - Confidence: medium. Strength: weak.
  - Evidence: [5][6][7]

- **GEN-04. ACKR1 (DARC) marks post-capillary and small collecting venules, the sites where leukocytes leave the
  blood. Arteries, arterioles, capillaries and most lymphatics lack it.**
  - Rule: none for cell typing. Use it for interpretation: neutrophils and monocytes in normal tissue cluster at
    ACKR1+ venules. Optional: `kind=neighborhood, types=[Neutrophils, Monocytes], neighbor_types=[Blood vessel EC],
    min_frac=0.07, mode=enrich, strength=weak`
  - Exceptions and caveats:
    - Under inflammation, neutrophils fill the tissue. Never `exclude` neutrophils in parenchyma, because they are
      the finding.
    - Lung and spleen hold marginated or resident neutrophils.
  - Confidence: medium. Strength: weak.
  - Evidence: [8][3]

- **GEN-05. Plasma cells sit in connective tissue next to, but not inside, epithelium: lamina propria of the gut and
  airways, around exocrine glands, in LN medullary cords, bone marrow and tertiary lymphoid structures, and around
  chronic inflammation.**
  - Rules:
    - `kind=neighborhood, types=[Plasma cells], neighbor_types=[Epithelial], min_frac=0.7, mode=exclude,
      strength=weak`
    - Gut: use a drawn `mucosa` region with `enrich` (see the gut sections).
  - Exceptions and caveats:
    - JCHAIN/IGHA1/IGKC transcripts are abundant and spill into adjacent epithelium. This rule partly protects
      against that error, which is why it is `exclude` rather than `enrich`.
    - Plasma cells inside epithelium are not expected; transmigration is of IgA, not cells.
  - Confidence: high. Strength: weak to moderate.
  - Evidence: [9]

- **GEN-06. Mast cells are connective-tissue residents, often perivascular and beneath epithelial surfaces: dermis,
  submucosa, lamina propria, airway wall and serosa. They are not inside epithelium or in lymphoid follicles in
  normal tissue.**
  - Rule: `kind=neighborhood, types=[Mast cells], neighbor_types=[Epithelial], min_frac=0.7, mode=exclude,
    strength=weak`
  - Exceptions and caveats:
    - Intraepithelial mast cells occur in allergic airway and oesophageal mucosa (asthma, eosinophilic oesophagitis).
      Turn the rule off for those conditions.
    - Mast cells vary considerably between tissues [10], so tissue-specific marker sets (TPSAB1 vs CMA1) do not
      change their niche.
  - Confidence: medium. Strength: weak.
  - Evidence: [10][11]

- **GEN-07. Intraepithelial lymphocytes are normal in gut, airway and skin. CD8 T cells, gd T cells and some NK/ILC
  cells sit within the epithelium.**
  - Rule (warning): do not write `T and NK cells exclude in [epithelium]` or `exclude where Epithelial >= f`. For
    the gut, `gd T cells` and `CD8 T cells` may be enriched at depth 0 (`kind=depth, max_depth_um=10, mode=enrich,
    strength=weak`).
  - Confidence: high (warning). Evidence: [9] (mucosal immune organization); textbook histology, not verified online,
    for intraepithelial lymphocyte frequencies.

- **GEN-08. Tertiary lymphoid structures (TLS) are organized B-cell aggregates with T-cell zones, FDCs and HEV-like
  vessels. They form in chronic inflammation, autoimmunity, transplants and many tumors.**
  - Rule:
    - Draw a `tertiary lymphoid structure` region and use `kind=region, types=[B cells], regions=[tertiary lymphoid
      structure], mode=enrich, strength=moderate`.
    - Without regions: `kind=neighborhood, types=[B cells], neighbor_types=[B cells, T and NK cells], min_frac=0.6,
      mode=enrich, strength=weak`.
  - Warning:
    - Never write organ-level rules such as "B cells only in lymph node / mucosa follicles" or "exclude B cells in
      tumor stroma". TLS are prognostically important and are a key finding the annotation must not hide.
    - GC B cells and FDCs are not in `immune.txt`. For mature TLS with germinal centres, add them from
      `lymph_node.txt` (mark them "not in the built-in list" for that organ).
  - Confidence: high. Strength: as given.
  - Evidence: [12]

- **GEN-09. Fibroblasts are present in the connective tissue of every organ and are absent from the inside of
  epithelial sheets and from glomerular tufts. Specialized fibroblast subtypes occupy discrete anatomical positions,
  e.g. PI16+ adventitial fibroblasts around larger vessels.**
  - Rule: `kind=neighborhood, types=[Fibroblasts], neighbor_types=[Epithelial], min_frac=0.8, mode=exclude,
    strength=weak`
  - Exceptions and caveats:
    - Epithelial-mesenchymal transition in carcinoma and sarcomatoid tumors breaks the rule.
    - Thin epithelia (alveoli, glomerular capsule) put fibroblasts within a few um of epithelium, so keep `min_frac`
      high.
    - PDGFRA is also expressed by gut telocytes and brain OPCs.
  - Confidence: medium-high. Strength: weak.
  - Evidence: [4][13]

- **GEN-10. Macrophages are tissue residents in every organ: Kupffer cells, alveolar macrophages, microglia,
  perivascular and serosal macrophages.**
  - Rule (warning): no `only` or `exclude` rule for Macrophages or Monocytes at organ level. At most use organ-specific
    `enrich` rules, such as alveolar lumen in lung or the LN sinuses.
  - LYVE1+ macrophages sit along vessels invested with smooth muscle [5]. This is a weak enrich candidate:
    `kind=neighborhood, types=[Macrophages], neighbor_types=[Smooth muscle cells, Endothelial], min_frac=0.2,
    mode=enrich, strength=weak`.
  - Confidence: high (warning); low (enrich). Evidence: [5][3]

- **GEN-11. Endothelial cells are everywhere at roughly capillary spacing, but their markers are tissue-specific
  [6]. Tissue type, more than vessel type, drives endothelial heterogeneity.**
  - Rule (warning): no location rule for Endothelial.
  - Check that the endothelial markers chosen exist on the panel and fit the organ:
    - CLDN5 is high in brain and lung and low elsewhere.
    - Liver sinusoidal EC are PECAM1-low and VWF-low.
    - Glomerular EC are fenestrated.
  - Confidence: high. Evidence: [6]

- **GEN-12. Lymphocyte aggregates in normal non-lymphoid organs are confined to defined sites.**
  - The sites are Peyer's patches and isolated lymphoid follicles in gut, BALT in airways (mostly in disease), and
    portal tracts in liver.
  - Rule: organ-specific, so use the organ sections.
  - Warning: in disease, ectopic aggregates (GEN-08) are expected. Pair any organ-level `exclude` of B cells with
    `weak` strength.
  - Confidence: medium. Evidence: [12]; organ specifics in the organ sections.

- **GEN-13. Method caveat for all neighborhood rules.**
  - Neighbour types come from expression only, so:
    - A rule can reinforce the first-pass error ("pericytes need endothelium", but endothelium was called pericyte).
    - Cells at tissue edges, in folds or next to tears have neighbours from the wrong compartment.
    - Small populations (pericytes, mast cells) fail `min_frac` thresholds if K = 15 neighbours span a large area in
      sparse tissue.
  - Prefer `min_frac` values of 0.07 to 0.2 for "touches a cell of type X" rules, and 0.5 or more for "embedded in
    tissue of type X" rules.
  - Confidence: high (method). No reference.

**References**
- [1] Armulik A, Genove G, Betsholtz C, Pericytes: developmental, physiological, and pathological perspectives,
  problems, and promises, Developmental Cell, 2011. doi:10.1016/j.devcel.2011.07.001.
- [2] Vanlandewijck M et al., A molecular atlas of cell types and zonation in the brain vasculature, Nature, 2018.
  doi:10.1038/nature25739.
- [3] Tabula Sapiens Consortium, The Tabula Sapiens: a multiple-organ, single-cell transcriptomic atlas of humans,
  Science, 2022. doi:10.1126/science.abl4896. PMID 35549404.
- [4] Muhl L et al., Single-cell analysis uncovers fibroblast heterogeneity and criteria for fibroblast and mural cell
  identification and discrimination, Nature Communications, 2020. doi:10.1038/s41467-020-17740-1.
- [5] Lim HY et al., Hyaluronan receptor LYVE-1-expressing macrophages maintain arterial tone through
  hyaluronan-mediated regulation of smooth muscle cell collagen, Immunity, 2018. doi:10.1016/j.immuni.2018.06.008.
- [6] Kalucka J et al., Single-cell transcriptome atlas of murine endothelial cells, Cell, 2020.
  doi:10.1016/j.cell.2020.01.015. PMID 32059779.
- [7] Oliver G, Kipnis J, Randolph GJ, Harvey NL, The lymphatic vasculature in the 21st century: novel functional roles
  in homeostasis and disease, Cell, 2020. doi:10.1016/j.cell.2020.06.039. PMID 32707093.
- [8] Thiriot A et al., Differential DARC/ACKR1 expression distinguishes venular from non-venular endothelial cells in
  murine tissues, BMC Biology, 2017. doi:10.1186/s12915-017-0381-7.
- [9] Brandtzaeg P, Johansen FE, Mucosal B cells: phenotypic characteristics, transcriptional regulation, and homing
  properties, Immunological Reviews, 2005. doi:10.1111/j.0105-2896.2005.00283.x.
- [10] Dwyer DF et al., Expression profiling of constitutive mast cells reveals a unique identity within the immune
  system, Nature Immunology, 2016. doi:10.1038/ni.3445.
- [11] Krystel-Whittemore M, Dileepan KN, Wood JG, Mast cell: a multi-functional master cell, Frontiers in Immunology,
  2016. doi:10.3389/fimmu.2015.00620.
- [12] Schumacher TN, Thommen DS, Tertiary lymphoid structures in cancer, Science, 2022. doi:10.1126/science.abf9419.
  PMID 34990248.
- [13] Buechler MB et al., Cross-tissue organization of the fibroblast lineage, Nature, 2021.
  doi:10.1038/s41586-021-03549-5.

Not verified online (stated as textbook): cortical thickness and agranular L4 (BRN-13), spleen layout (LN-12), PROX1
in hepatocytes (GEN-03), visceral smooth-muscle layers (GEN-02), intraepithelial lymphocyte frequencies (GEN-07), and
memory B cells beneath tonsil crypt epithelium (LN-11, no rule given). The perivascular location of mast cells
(GEN-06) rests on the Krystel-Whittemore abstract (vascular roles). The "beneath epithelial surfaces" and "near
nerves" parts were not seen in text.
