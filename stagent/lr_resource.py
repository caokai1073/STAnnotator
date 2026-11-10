"""Build the ligand-receptor resource of the cell-cell communication analysis: stagent/data/lr_human.csv.gz and
lr_mouse.csv.gz, one row per ligand-receptor pair with the genes of each side (complexes joined by "_"), a pathway,
the kind of signaling and the resources that list it.

Sources (downloaded):
  CellChatDB v2 (jinworks/CellChat, GPL-3.0): human and mouse; pathways and the kind of signaling (Secreted
    Signaling, ECM-Receptor, Cell-Cell Contact, Non-protein Signaling).  Reading its .rda needs `pip install rdata`.
  CellPhoneDB v5 (ventolab/cellphonedb-data, CellPhoneDB is MIT): human; classification used as the pathway;
    secreted partners give Secreted Signaling, others Cell-Cell Contact; small molecules made by an enzyme
    ("<molecule>_by<ENZYME>") are Non-protein Signaling.
  OmniPath consensus as used by LIANA (saezlab/liana-py, BSD-3-Clause): "consensus" (human) and "mouseconsensus",
    built from CellPhoneDB, CellChatDB, ICELLNET, connectomeDB2020 and CellTalkDB.
Pairs listed by several resources are merged (same genes on each side, in any order); pathway and kind come from
CellChatDB, else CellPhoneDB, else the most common annotation of the same ligand in CellChatDB.

  python -m stagent.lr_resource [download folder]
"""
import sys, urllib.request
from collections import Counter
from pathlib import Path
import pandas as pd

OUT = Path(__file__).resolve().parent / "data"
URLS = {
    "cellchat_human.rda": "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.human.rda",
    "cellchat_mouse.rda": "https://raw.githubusercontent.com/jinworks/CellChat/main/data/CellChatDB.mouse.rda",
    "omni_resource.csv": "https://raw.githubusercontent.com/saezlab/liana-py/main/src/liana/resource/omni_resource.csv",
    **{f"cpdb_{f}.csv": f"https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/{f}.csv"
       for f in ("interaction_input", "complex_input", "gene_input", "protein_input")},
}


def fetch(d):
    d.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        if not (d / name).exists():
            print("downloading", url)
            urllib.request.urlretrieve(url, d / name)


def cellchat(d, species):
    import rdata
    db = rdata.read_rda(str(d / f"cellchat_{species}.rda"))[f"CellChatDB.{species}"]
    it, cx = db["interaction"], db["complex"]
    sub = {n: [g for g in r if isinstance(g, str) and g] for n, r in cx.iterrows()}
    side = lambda x: "_".join(sub.get(x, [x]))
    return pd.DataFrame({"ligand": it.ligand.map(side), "receptor": it.receptor.map(side), "pathway": it.pathway_name,
                         "category": it.annotation, "source": "CellChatDB"})


def cellphonedb(d):
    g = pd.read_csv(d / "cpdb_gene_input.csv")
    p = pd.read_csv(d / "cpdb_protein_input.csv")
    c = pd.read_csv(d / "cpdb_complex_input.csv")
    it = pd.read_csv(d / "cpdb_interaction_input.csv")
    u2g = g.drop_duplicates("uniprot").set_index("uniprot").gene_name.to_dict()
    ucols = [x for x in c.columns if x.startswith("uniprot_")]
    cx = {r.complex_name: [u2g.get(r[k]) for k in ucols if isinstance(r[k], str)] for _, r in c.iterrows()}
    secreted = {**p.set_index("uniprot").secreted.to_dict(), **c.set_index("complex_name").secreted.to_dict()}

    def side(x):
        if x in cx:
            return "_".join(cx[x]) if all(cx[x]) else None
        return u2g.get(x)
    rows = []
    for r in it.itertuples():
        a, b = side(r.partner_a), side(r.partner_b)
        if not a or not b:
            continue
        kind = ("Non-protein Signaling" if "_by" in r.partner_a else
                "Secreted Signaling" if secreted.get(r.partner_a) in (True, "True") else "Cell-Cell Contact")
        rows.append({"ligand": a, "receptor": b, "pathway": ("" if pd.isna(r.classification) else str(r.classification)).replace("Signaling by ", ""),
                     "category": kind, "source": "CellPhoneDB"})
    return pd.DataFrame(rows)


def consensus(d, species):
    o = pd.read_csv(d / "omni_resource.csv", index_col=0)
    o = o[o.resource == ("consensus" if species == "human" else "mouseconsensus")]
    return pd.DataFrame({"ligand": o.source_genesymbol, "receptor": o.target_genesymbol, "pathway": "", "category": "",
                         "source": "consensus"})


def merge(parts):
    t = pd.concat(parts, ignore_index=True).fillna("")
    t["key"] = [("_".join(sorted(a.split("_"))), "_".join(sorted(b.split("_")))) for a, b in zip(t.ligand, t.receptor)]
    by_lig = {}
    for r in t[t.source == "CellChatDB"].itertuples():  # the usual pathway / kind of each ligand
        for g in r.ligand.split("_"):
            by_lig.setdefault(g.upper(), Counter())[(r.pathway, r.category)] += 1
    rows = []
    for key, grp in t.groupby("key", sort=False):
        first = grp.iloc[0]
        ann = grp[grp.category != ""]
        if len(ann):
            pathway, category = ann.iloc[0].pathway, ann.iloc[0].category
        else:
            guess = Counter()
            for g in first.ligand.split("_"):
                guess.update(by_lig.get(g.upper(), Counter()))
            pathway, category = guess.most_common(1)[0][0] if guess else ("", "")
        rows.append({"ligand": first.ligand, "receptor": first.receptor, "pathway": pathway, "category": category,
                     "sources": ";".join(dict.fromkeys(grp.source))})
    out = pd.DataFrame(rows)
    out.insert(0, "interaction", out.ligand.str.replace("_", "+") + " → " + out.receptor.str.replace("_", "+"))
    return out


def main():
    d = Path(sys.argv[1] if len(sys.argv) > 1 else "lr_sources")
    fetch(d)
    for species in ("human", "mouse"):
        parts = [cellchat(d, species)] + ([cellphonedb(d)] if species == "human" else []) + [consensus(d, species)]
        t = merge(parts)
        t.to_csv(OUT / f"lr_{species}.csv.gz", index=False, compression="gzip")
        print(species, len(t), "pairs;", t.category.replace("", "unknown").value_counts().to_dict())


if __name__ == "__main__":
    main()
