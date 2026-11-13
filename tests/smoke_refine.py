"""Refining a marker annotation round by round (stagent/iterate.py), on the bundled mouse colon example:
1. every cell gets a cell type (markers.ASSIGN_ALL): no Unassigned;
2. a round's report: organ, types, weak evidence, marker flags, types not found;
3. the agent's proposal (edits checked against the list), the user's message, some edits accepted: the next marker
   list has exactly those edits; a bad edit is refused;
4. annotating with that list and carrying the rounds: the next round compares with the one before.
Run: python tests/smoke_refine.py  (a temporary data folder removed at the end)."""
import os, shutil, tempfile
from pathlib import Path

HOME = Path(tempfile.mkdtemp(prefix="stagent_refine_"))
os.environ["STAGENT_HOME"] = str(HOME)
import numpy as np  # noqa: E402
from stagent import importer, iterate, markers  # noqa: E402
from stagent.session import Session  # noqa: E402

EX = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish"
try:
    importer.import_sample("M", str(EX / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    text = (EX / "mouse_colon_merfish_markers.txt").read_text()
    s = Session.create("M", method="markers", markers_text=text, knn="exact")
    lab = np.asarray(s.labels())
    assert not (lab == markers.UNASSIGNED).any(), "cells left Unassigned"
    st = iterate.start_round(s, organ="Colon")
    rep = st["rounds"][-1]["report"]
    assert rep["n_cells"] == len(lab) and rep["types"] and "weak_types" in rep and "not_found" in rep
    print(f"round 1: {len(rep['types'])} types, {len(rep['weak_types'])} weak, {len(rep['not_found'])} not found, "
          f"{rep['untyped_by_markers']['n_cells']} cells the marker step could not type")
    first = next(n for n in markers.parse(text, None)["tree"] if n["children"])
    kid = first["children"][0]["name"]
    gene = first["children"][0]["markers"][0]["genes"][0]
    iterate.propose(s, "test", [
        {"kind": "remove_marker", "type": kid, "marker": gene + "+", "reason": "r", "evidence": "e"},
        {"kind": "rename_type", "type": kid, "new_name": kid + " (renamed)", "reason": "r", "evidence": "e"},
        {"kind": "add_type", "type": "Test type", "parent": first["name"], "markers": "Epcam+", "reason": "r", "evidence": "e"}],
        questions=["Healthy or inflamed?"])
    rd = iterate.add_message(s, "user", "Keep the name")
    assert [m["from"] for m in rd["messages"]] == ["agent", "user"]
    try:
        iterate.propose(s, "bad", [{"kind": "explode_type", "type": kid}])
        raise AssertionError("an unknown edit was accepted")
    except ValueError:
        pass
    iterate.propose(s, "again", [
        {"kind": "remove_marker", "type": kid, "marker": gene + "+", "reason": "r", "evidence": "e"},
        {"kind": "rename_type", "type": kid, "new_name": kid + " (renamed)", "reason": "r", "evidence": "e"}])
    out = iterate.accept(s, ["E1"])
    new = markers.parse(out["marker_text"], s.panel)
    node = markers._find(new["tree"], kid)
    assert node is not None and gene not in {g for m in node["markers"] for g in m["genes"]}, "E1 not applied"
    assert markers._find(new["tree"], kid + " (renamed)") is None, "E2 applied although not chosen"
    assert iterate.load(s)["rounds"][-1]["decision"] == "accepted"
    s2 = Session.create("M", method="markers", markers_text=out["marker_text"], knn="exact", carry_from=s.id)
    iterate.carry(s, s2)
    st2 = iterate.start_round(s2)
    assert len(st2["rounds"]) == 2 and st2["organ"] == "Colon" and "changes" in st2["rounds"][-1]["report"]
    print(f"round 2 after accepting E1: {len(st2['rounds'][-1]['report']['changes'])} types changed")
    print("ok")
finally:
    shutil.rmtree(HOME, ignore_errors=True)
