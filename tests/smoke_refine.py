"""Refining a marker annotation round by round (stagent/iterate.py), on the bundled mouse colon example:
1. every cell the marker step could not type gets a type from its marker genes (markers.ASSIGN_ALL, no model), or
   stays Unassigned when its marker genes do not say what it is (none detected, or a Bayes factor below 3; 0.2.57-0.2.58
   labeled these Low quality);
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
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(os.environ["STAGENT_HOME"], "cfg")  # (never the user's settings)
import numpy as np  # noqa: E402
from stagent import importer, iterate, markers  # noqa: E402
from stagent.session import Session  # noqa: E402

EX = Path(__file__).resolve().parents[1] / "examples" / "cadinu_merfish"
try:
    importer.import_sample("M", str(EX / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    text = (EX / "mouse_colon_merfish_markers.txt").read_text()
    s = Session.create("M", method="markers", markers_text=text, knn="exact")
    lab = np.asarray(s.labels())
    assert markers.LOW_QUALITY not in s.types and not (lab == markers.LOW_QUALITY).any(), "a new annotation uses Low quality"
    left = np.where(lab == markers.UNASSIGNED)[0]
    res = s.marker_result()
    untyped = np.where(np.asarray(res["labels"]) == markers.UNASSIGNED)[0]
    assert set(left) <= set(untyped), "a cell typed by the marker step became Unassigned"
    found = {s.types[j] for j in np.unique(markers.training_labels(res, s.types)) if j >= 0} - set(markers.NO_TYPE)
    assert set(np.unique(lab[untyped])) - set(markers.NO_TYPE) <= found, "a type the marker step did not find got cells"
    # Unassigned: exactly the untyped cells with no marker gene detected or too little evidence (Bayes factor < 3)
    xy = s.cells[["x", "y"]].to_numpy(float)
    _, best, ev = markers.marker_assign(s.counts, s.panel, s._marker_tree([]), s.types, untyped,
                                        markers.training_labels(res, s.types),
                                        spill_xy=xy if s.scheme.get("spillover") else None, evidence=True)
    want_lq = (best < 0) | (ev < np.log(markers.LOW_QUALITY_BF))
    assert set(untyped[want_lq]) == set(left), "Unassigned is not the cells without marker evidence"
    print(f"{len(untyped):,} cells the marker step could not type: {len(untyped) - len(left):,} got a type from "
          f"their marker genes, {len(left):,} stay Unassigned ({int((best < 0).sum()):,} with no marker detected)")
    # a root cluster of the marker step whose cells all went elsewhere (rare-type search) has no members: skipped
    real = s.marker_result
    def emptied(*a, **k):
        ml = dict(real(*a, **k)); ml["members"] = dict(ml["members"])
        r = next(x for x in ml["records"] if x["parent"] is None and not x["assigned"])
        ml["members"][r["id"]] = np.zeros(0, np.int32)
        return ml
    s.marker_result = emptied
    assert iterate.report(s)["untyped_by_markers"]["n_cells"] >= 0
    del s.marker_result
    st = iterate.start_round(s, organ="Colon")
    rep = st["rounds"][-1]["report"]
    assert rep["n_cells"] == len(lab) and rep["types"] and "weak_types" in rep and "not_found" in rep
    print(f"round 1: {len(rep['types'])} types, {len(rep['weak_types'])} weak, {len(rep['not_found'])} not found, "
          f"{rep['untyped_by_markers']['n_cells']} cells the marker step could not type")
    first = next(n for n in markers.parse(text, None)["tree"] if n["children"])
    kid = first["children"][0]["name"]
    gene = first["children"][0]["markers"][0]["genes"][0]
    # the edits' groups: the problem the agent named, else the lineage the edit changes (a renamed type's too)
    other = next(n["name"] for n in markers.parse(text, None)["tree"] if n["name"] != first["name"])
    eds = iterate.propose(s, "groups", [
        {"kind": "rename_type", "type": kid, "new_name": "Renamed X"},
        {"kind": "add_marker", "type": "Renamed X", "marker": gene + "-"},
        {"kind": "add_type", "type": "Top Y", "parent": None, "markers": "Epcam+"},
        {"kind": "remove_marker", "type": other, "marker": "X+", "problem": "  One   problem "},
        {"kind": "add_marker", "type": first["name"], "marker": gene + "+", "problem": "One problem"}])["edits"]
    assert [e["group"] for e in eds] == [first["name"], first["name"], "Top Y", "One problem", "One problem"], eds
    assert [e["id"] for e in eds] == ["E1", "E2", "E3", "E4", "E5"]
    assert [e["optional"] for e in eds] == [False] * 5
    eo = iterate.propose(s, "optional", [{"kind": "add_type", "type": "Split A", "parent": first["name"], "markers": "Epcam+",
                                          "optional": "yes", "question": "  Split   it? "}])["edits"]
    assert eo[0]["optional"] is True and eo[0]["question"] == "Split it?", eo
    ef = iterate.propose(s, "optional", [{"kind": "add_type", "type": "Split A", "parent": first["name"], "markers": "Epcam+",
                                          "optional": "false"}])["edits"]
    assert ef[0]["optional"] is False, ef
    # where the types lie (computed), the user's rules gone against (none here)
    lie = rep["where_types_lie"]
    assert lie and all(set(x) >= {"type", "lineage", "n_cells", "among", "mostly_among"} for x in lie), lie[:2]
    assert all(x["mostly_among"] is None or x["among"][x["mostly_among"]] > 0.5 for x in lie)
    assert rep["rule_conflicts"] == []
    # a finding's map: its cells are those of the type lying among the lineage (as the report counts them)
    from starlette.testclient import TestClient
    from stagent import webapp
    c1 = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    lt = next(x for x in lie if x["among"])
    la = next(iter(lt["among"]))
    pf = iterate.propose(s, "figs", [], findings=[{"cells": "x", "type": lt["type"], "among": la},
                                                  {"cells": f"some {lt['type']} cells"}, {"cells": "nothing known"}])
    rows, of_type, title = iterate.finding_cells(s, pf["findings"][0])
    assert set(rows) <= set(of_type) and abs(len(rows) / len(of_type) - lt["among"][la]) < 0.002, (len(rows), len(of_type), lt)
    assert len(iterate.finding_cells(s, pf["findings"][1])[0]) == lt["n_cells"]  # (the type found in the text)
    r1 = c1.get(f"/api/s/{s.id}/refine/finding_fig", params={"key": pf["findings"][0]["key"]})
    assert r1.status_code == 200 and r1.content[:4] == b"\x89PNG", r1.status_code
    r2 = c1.get(f"/api/s/{s.id}/refine/finding_fig", params={"key": pf["findings"][0]["key"], "size": "large"})
    from PIL import Image as _I
    import io as _io
    assert r2.status_code == 200 and _I.open(_io.BytesIO(r2.content)).width > 2 * _I.open(_io.BytesIO(r1.content)).width
    assert c1.get(f"/api/s/{s.id}/refine/finding_fig", params={"key": pf["findings"][2]["key"]}).status_code == 404
    assert c1.get(f"/api/s/{s.id}/refine/finding_fig", params={"key": "nope"}).status_code == 404
    print(f"finding map: {title}")
    # findings against the priors; the user's verdicts carry over to later rounds; a bad finding is refused
    p = iterate.propose(s, "findings", [], findings=[
        {"cells": "Neurons", "where": "71% among epithelial cells", "rule": "COLON-03 neurons in the plexus",
         "judgement": "Likely mislabel", "evidence": "Nos1 in 23%"},
        {"cells": "Epithelial", "where": "in the submucosa", "rule": "COLON-01", "judgement": "whatever"}])
    assert [f["judgement"] for f in p["findings"]] == ["likely mislabel", "unsure"]
    k0 = p["findings"][0]["key"]
    assert iterate.finding_key({"cells": " neurons", "rule": "COLON-03  neurons in the PLEXUS"}) == k0
    iterate.set_verdict(s, k0, "biology")
    assert iterate.load(s)["prior_verdicts"][k0]["verdict"] == "biology"
    for bad in ((k0, "maybe"), ("no such", "biology")):
        try:
            iterate.set_verdict(s, *bad); raise AssertionError(f"accepted {bad}")
        except ValueError:
            pass
    try:
        iterate.propose(s, "bad", [], findings=[{"where": "x"}]); raise AssertionError("a finding without cells")
    except ValueError:
        pass
    iterate.set_verdict(s, k0, None)
    assert k0 not in iterate.load(s)["prior_verdicts"]
    iterate.set_verdict(s, k0, "mislabel")
    iterate.propose(s, "test", [
        {"kind": "remove_marker", "type": kid, "marker": gene + "+", "reason": "r", "evidence": "e"},
        {"kind": "rename_type", "type": kid, "new_name": kid + " (renamed)", "reason": "r", "evidence": "e"},
        {"kind": "add_type", "type": "Test type", "parent": first["name"], "markers": "Epcam+", "reason": "r", "evidence": "e"}],
        questions=["Healthy or inflamed?"])
    rd = iterate.add_message(s, "user", "Keep the name")
    assert [m["from"] for m in rd["messages"]] == ["agent"] * 6 + ["user"], [m["from"] for m in rd["messages"]]  # (the groups, maps and findings checks proposed too)
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

    from stagent import mcp_server, markerqc, claudelink, advice
    import json
    mcp_server._active["session"] = s  # the agent still on the annotation of round 1 (decided)
    r = json.loads(mcp_server.refine_report(session_id=s2.id))
    assert r["session_id"] == s2.id and r["round"] == 2 and not r["decided"], r["session_id"]
    # the verdict of round 1 is in round 2's report for the agent
    assert [v["verdict"] for v in r["user_verdicts_on_prior_findings"]] == ["mislabel"], r["user_verdicts_on_prior_findings"]
    from starlette.testclient import TestClient
    from stagent import webapp
    c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
    for body in (None, [], {"key": 3}, {"key": "x", "verdict": "maybe"}):
        assert c.post(f"/api/s/{s2.id}/refine/finding", json=body).status_code == 400, body
    assert c.post(f"/api/s/{s2.id}/refine/finding", json={"key": "nope", "verdict": "biology"}).status_code == 400
    k = r["user_verdicts_on_prior_findings"] and iterate.load(s2)["prior_verdicts"]
    k0 = next(iter(k))
    out = c.post(f"/api/s/{s2.id}/refine/finding", json={"key": k0, "verdict": "biology"}).json()
    assert out.get("prior_verdicts", {}).get(k0, {}).get("verdict") == "biology", out
    # offline (the desktop app; the computer has no network at all): nothing is sent, the round is not waiting for the
    # agent, the user's words are not stored; any other request to the agent: the reason (503). A network that is up
    # is never refused (Claude Code may reach Anthropic through a proxy the app does not know of).
    os.environ["STAGENT_DESKTOP"] = "1"
    route = claudelink._has_route
    claudelink._has_route = lambda: False
    try:
        n_msgs = len(iterate.current(s2)[1]["messages"])
        r = c.post(f"/api/s/{s2.id}/refine/ask", json={"message": "look at the B cells"})
        assert r.status_code == 400 and r.json()["error"].startswith("No internet"), r.text
        r = c.post(f"/api/s/{s2.id}/refine/message", json={"text": "and the T cells"})
        assert r.status_code == 400 and r.json()["error"].startswith("No internet"), r.text
        assert not iterate.awaiting(s2) and len(iterate.current(s2)[1]["messages"]) == n_msgs, "recorded as sent"
        r = c.post("/api/claude/report", json={"session": s2.id})
        assert r.status_code == 503 and r.json()["offline"], r.text
        before = advice.status(s2, "annotation")
        r = c.post("/api/claude/review", json={"session": s2.id, "scope": "annotation"})
        assert r.status_code == 503 and advice.status(s2, "annotation") == before, "the review recorded as asked"
        os.environ["STAGENT_SKIP_ONLINE_CHECK"] = "1"
        claudelink.check_online()
    finally:
        claudelink._has_route = route
        os.environ.pop("STAGENT_DESKTOP", None); os.environ.pop("STAGENT_SKIP_ONLINE_CHECK", None)
    assert claudelink._has_route() in (True, False)
    # the time the round has waited for the agent (the window points to Claude Code's window after 10 minutes)
    iterate.mark_asked(s2)
    st_, rd_ = iterate.current(s2); rd_["asked_t"] -= 700; iterate.save(s2, st_)
    w = c.get(f"/api/s/{s2.id}/refine").json()
    assert w["awaiting"] and w["waited_s"] >= 700, (w["awaiting"], w["waited_s"])
    st_, rd_ = iterate.current(s2); rd_.pop("asked"); rd_.pop("asked_t"); iterate.save(s2, st_)
    assert c.get(f"/api/s/{s2.id}/refine").json()["waited_s"] is None
    print("offline: nothing sent; the time waited")
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
