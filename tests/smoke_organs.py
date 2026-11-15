"""The built-in organ marker lists (stagent/data/organs): every list parses with its @includes, type names are unique,
a type does not list a gene both as a positive and a negative marker, every type has two or more positive markers or a
decisive one, and on a panel holding every gene of the list the draft keeps all of its types; on a panel of a few
genes the draft drops what it cannot resolve and says why.
Run: python tests/smoke_organs.py"""
import os, tempfile
os.environ.setdefault("STAGENT_HOME", tempfile.mkdtemp(prefix="stagent_organs_"))
os.environ["STAGENT_SETTINGS_DIR"] = tempfile.mkdtemp(prefix="stagent_cfg_")  # (never the user's settings)
from stagent import organs, markers  # noqa: E402

names = [o["organ"] for o in organs.organs()]
assert {"colon", "lymph_node", "lung", "other"} <= set(names), names
for organ in names:
    nodes = organs.tree(organ)
    seen = [n["name"] for n in nodes]
    assert len(seen) == len(set(seen)), f"{organ}: duplicate type names"
    genes = set()
    for n in nodes:
        pos = [g for g, _ in n["markers"]]
        assert not set(pos) & set(n["negative"]), f"{organ}/{n['name']}: a gene is positive and negative"
        assert len(pos) >= 2 or any(d for _, d in n["markers"]), f"{organ}/{n['name']}: fewer than two markers"
        genes |= set(pos) | set(n["negative"])
    d = organs.draft(organ, sorted(genes))
    assert not d["dropped"], f"{organ}: dropped on its own genes: {d['dropped']}"
    leaves = {n["name"] for n in nodes} - {n["parent"] for n in nodes}
    assert d["n_types"] == len(leaves), (organ, d["n_types"], len(leaves))
    text = markers.csv_to_text(d["csv"])
    r = markers.parse(text, sorted(genes))
    assert r["ok"], (organ, r["errors"])
small = organs.draft("lymph_node", ["CD3E", "CD2", "MS4A1", "CD79A", "CD4", "CD40LG", "CD8A", "CD8B"])
assert "CD4 T cells" in small["kept"] and "CD8 T cells" in small["kept"], small["kept"]
assert any(x["type"] == "Macrophages" for x in small["dropped"]), small["dropped"]
print(f"{len(names)} organ lists ok")
