"""Refining a marker annotation round by round (stagent/iterate.py), on the bundled mouse colon example:
1. every cell the marker step could not type gets a type from its marker genes (markers.ASSIGN_ALL, no model), or
   Low quality when its marker genes do not say what it is (none detected, or a Bayes factor below 3);
2. a round's report: organ, types, weak evidence, marker flags, types not found;
3. the agent's proposal (edits checked against the list), the user's message, some edits accepted: the next marker
   list has exactly those edits; a bad edit is refused;
4. annotating with that list and carrying the rounds: the next round compares with the one before;
6. the round after Accept is on a new annotation: refine_report / propose_refinement follow the session_id of the
   request, or without it the annotation open in the app when the active one's round is decided;
5. a whole revised list the agent sends with propose_marker_list while the round waits (Agent review) goes to the
   Refine window as edits, which applied give that list.
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
    assert not (lab == markers.UNASSIGNED).any(), "a new annotation uses Low quality, not Unassigned"
    left = np.where(lab == markers.LOW_QUALITY)[0]
    res = s.marker_result()
    untyped = np.where(np.asarray(res["labels"]) == markers.UNASSIGNED)[0]
    assert set(left) <= set(untyped), "a cell typed by the marker step became Low quality"
    found = {s.types[j] for j in np.unique(markers.training_labels(res, s.types)) if j >= 0} - set(markers.NO_TYPE)
    assert set(np.unique(lab[untyped])) - set(markers.NO_TYPE) <= found, "a type the marker step did not find got cells"
    # Low quality: exactly the untyped cells with no marker gene detected or too little evidence (Bayes factor < 3)
    xy = s.cells[["x", "y"]].to_numpy(float)
    _, best, ev = markers.marker_assign(s.counts, s.panel, s._marker_tree([]), s.types, untyped,
                                        markers.training_labels(res, s.types),
                                        spill_xy=xy if s.scheme.get("spillover") else None, evidence=True)
    want_lq = (best < 0) | (ev < np.log(markers.LOW_QUALITY_BF))
    assert set(untyped[want_lq]) == set(left), "Low quality is not the cells without marker evidence"
    print(f"{len(untyped):,} cells the marker step could not type: {len(untyped) - len(left):,} got a type from "
          f"their marker genes, {len(left):,} are Low quality ({int((best < 0).sum()):,} with no marker detected)")
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

    from stagent import mcp_server, markerqc, claudelink
    import json
    mcp_server._active["session"] = s  # the agent still on the annotation of round 1 (decided)
    r = json.loads(mcp_server.refine_report(session_id=s2.id))
    assert r["session_id"] == s2.id and r["round"] == 2 and not r["decided"], r["session_id"]
    mcp_server._active["session"] = s
    real = claudelink.viewer_state
    claudelink.viewer_state = lambda: {"session": s2.id}
    try:
        r = json.loads(mcp_server.refine_report())
        assert r["session_id"] == s2.id and r["round"] == 2, "did not follow the app to the new round"
    finally:
        claudelink.viewer_state = real
    mcp_server._active["session"] = s
    try:
        mcp_server.propose_refinement("x", [], session_id=s.id)
        raise AssertionError("a proposal went to a decided round")
    except ValueError:
        pass
    print("the agent follows the round to the new annotation")
    old = st2["rounds"][-1]["marker_text"]
    tree = markers.parse(old, None)["tree"]
    epi = tree[0]["name"]
    kids = tree[0]["children"]
    g1, g2, g3 = [g for g in s2.panel if g not in old][:3]
    new = old.replace(f"  {kids[0]['name']}:", f"  {kids[0]['name']}: {g1}+", 1)  # a marker added
    new = "\n".join(l for l in new.splitlines() if not l.strip().startswith(kids[-1]["name"] + ":"))  # a type removed
    new += f"\nTest top type: {g2}+ {g3}+"  # a type added at the top
    mcp_server._active["session"] = s2
    assert not iterate.awaiting(s2)
    iterate.mark_asked(s2)
    assert iterate.awaiting(s2)
    r = json.loads(mcp_server.propose_marker_list(new, "a whole list"))
    assert "Refine window" in r["sent_to"], r
    assert not iterate.awaiting(s2) and markerqc.proposal(s2) is None, "the list went to the Review panel"
    edits = iterate.current(s2)[1]["proposal"]["edits"]
    got, _ = iterate.apply_edits(old, edits)
    flat = lambda t: sorted(l.strip() for l in markers.to_text(markers.parse(t, None)["tree"]).splitlines())
    assert flat(got) == flat(new), (flat(got), flat(new))
    print("whole list during Agent review ->", ", ".join(f"{e['kind']} {e['type']}" for e in edits))
    r = json.loads(mcp_server.propose_marker_list(new, "not asked"))  # not waiting: the Review panel as before
    assert "sent_to" not in r and markerqc.proposal(s2) is not None

    # 7. the round's decisions: Dismiss, then asking again reopens it; accept only the proposal seen, and once
    rd = iterate.current(s2)[1]
    iterate.decide(s2, "dismissed")
    assert not iterate.awaiting(s2)
    iterate.mark_asked(s2)
    assert iterate.awaiting(s2) and not iterate.current(s2)[1]["decision"], "asking again did not reopen the round"
    p1 = iterate.propose(s2, "again", [{"kind": "remove_marker", "type": kid, "marker": gene + "+", "reason": "r", "evidence": "e"}])
    try:
        iterate.accept(s2, ["E1"], proposal_time="2000-01-01T00:00:00")
        raise AssertionError("an edit of an older proposal was accepted")
    except ValueError:
        pass
    iterate.accept(s2, ["E1"], proposal_time=p1["time"])
    for bad in (lambda: iterate.accept(s2, ["E1"]), lambda: iterate.decide(s2, "dismissed"),
                lambda: iterate.propose(s2, "late", [])):
        try:
            bad()
            raise AssertionError("an accepted round was changed")
        except ValueError:
            pass
    # rounds are not carried to another sample's annotation
    importer.import_sample("N", str(EX / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    n = Session.create("N", method="markers", markers_text=text, knn="exact")
    iterate.carry(s2, n)
    assert not iterate.load(n)["rounds"], "rounds carried to another sample"
    # remove_marker keeps the other sign; a malformed revised list is refused, not read as removals
    neg = next((t, m["genes"][0]) for t in markers.descendants(tree) for m in markers._find(tree, t)["markers"] if m["sign"] < 0)
    out, _ = iterate.apply_edits(old, [{"kind": "remove_marker", "type": neg[0], "marker": neg[1] + "+"}])
    assert any(neg[1] in m["genes"] and m["sign"] < 0 for m in markers._find(markers.parse(out, None)["tree"], neg[0])["markers"]), \
        "removing GENE+ removed GENE-"
    try:
        iterate.edits_from_text(old, old.replace(f"{kids[0]['name']}:", kids[0]["name"], 1))
        raise AssertionError("a malformed list was read")
    except ValueError:
        pass
    print("round decisions, carry and edits checked")
    print("ok")
finally:
    shutil.rmtree(HOME, ignore_errors=True)
