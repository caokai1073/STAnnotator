"""A project of two slices (the example MERFISH slice cut in halves): import, joint annotation, joint embedding with
Harmony, joint niches, composition, DE genes and cell-cell communication. Runs in a temporary data folder."""
import os, shutil, tempfile, time
import h5py
import numpy as np
EX = os.path.join(os.path.dirname(__file__), "..", "examples", "cadinu_merfish")
home = tempfile.mkdtemp(prefix="stagent_joint_")
os.environ["STAGENT_HOME"] = home
os.environ["STAGENT_SETTINGS_DIR"] = os.path.join(os.environ["STAGENT_HOME"], "cfg")  # (never the user's settings)
from stagent import analysis, communication, importer, joint, projects
try:
    h5 = os.path.join(home, "two.h5ad")
    shutil.copy(os.path.join(EX, "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), h5)
    with h5py.File(h5, "r+") as h:  # a categorical obs column "half": left / right of the middle
        x = h["obs"]["x"][:]
        g = h["obs"].create_group("half")
        g.attrs.update({"encoding-type": "categorical", "encoding-version": "0.2.0", "ordered": False})
        g["codes"] = (x >= (x.min() + x.max()) / 2).astype(np.int8)
        g.create_dataset("categories", data=np.array(["left", "right"], dtype=object), dtype=h5py.string_dtype())
    t = time.time()
    importer.import_project("demo", "h5ad", [h5], sections=["left", "right"], section_column="half", progress=lambda *a: None)
    p = projects.list_projects()[0]
    print("project:", p["project_id"], p["samples"], [p["slices"][s]["metadata"] for s in p["samples"]], f"{time.time() - t:.0f}s")
    # niches of NMF programs over all slices before any annotation: kept in each slice's placeholder, then taken over
    # by the slices' first annotation (cell types as the identity still needs every slice annotated)
    from stagent import niches as nmod, session as smod
    try:
        joint.niches("demo", n_niches=4, identity="cell types")
        raise AssertionError("niches of cell types ran on slices not annotated")
    except ValueError as e:
        assert "annotated" in str(e), e
    joint.niches("demo", n_niches=4, n_programs=8)
    ph = {sl: smod.placeholder(sl) for sl in p["samples"]}
    assert all(ph.values()) and all(nmod.summary(smod.Session(x))["params"]["identity"] == "nmf" for x in ph.values()), ph
    # the agent's tools too, on a slice not annotated yet: each slice's density, then NMF niches of all slices
    import json as _json
    from stagent import mcp_server
    mcp_server._active["session"] = smod.Session(ph[p["samples"][0]])
    ctx = _json.loads(mcp_server.niche_context(all_slices=True))
    assert set(ctx["slices"]) == set(p["samples"]) and all(v["regions"] for v in ctx["slices"].values()), ctx
    assert mcp_server._d_run_project_niches({"k": 20, "n_niches": 5, "identity": "nmf", "n_programs": 8})["title"]
    r = _json.loads(mcp_server.GATED["run_project_niches"](k=20, n_niches=5, identity="nmf", n_programs=8))
    assert r["params"]["n_niches"] == 5 and set(r["slices"]) == set(p["samples"]), r
    try:
        mcp_server._d_run_project_niches({"k": 20, "n_niches": 4, "identity": "cell types", "n_programs": 8})
        raise AssertionError("the agent's niches of cell types ran on slices not annotated")
    except ValueError as e:
        assert "annotated" in str(e)
    mcp_server._active.pop("session", None)
    joint.niches("demo", n_niches=4, n_programs=8)  # (the run the annotation below takes over)
    t = time.time()
    joint.annotate("demo", open(os.path.join(EX, "mouse_colon_merfish_markers.txt")).read())
    for sl, sid in joint.latest_sessions("demo").items():
        sm = nmod.summary(smod.Session(sid))
        assert sm and sm["params"]["identity"] == "nmf" and sm["params"]["n_niches"] == 4, (sl, sm and sm["params"])
    print("joint annotation:", {s: sid for s, sid in joint.latest_sessions("demo").items()}, f"{time.time() - t:.0f}s")
    t = time.time()
    e = joint.embed("demo", harmony=True)
    print("joint embedding:", {k: e[k] for k in e if k not in ("samples",)}, f"{time.time() - t:.0f}s")
    t = time.time()
    joint.niches("demo", n_niches=5)
    print("joint niches:", f"{time.time() - t:.0f}s")
    # niches of one slice after a run of all slices use the slice's own NMF programs, not the project's (an exported
    # notebook fits the slice again): the project's are kept apart; a file without the counts key (as runs of all
    # slices wrote before) is not taken for the slice's own
    import numpy as np
    from stagent import config
    s0 = smod.Session(list(joint.latest_sessions("demo").values())[0])
    nd = s0.dir / "niches"
    Hp = np.load(nd / "nmf_20_project.npz")["H"]
    nmod.run(s0, n_niches=5)
    W1, H1 = nmod.programs(s0)
    assert H1.shape == Hp.shape and not np.allclose(H1, Hp), "the slice took the project's programs"
    lab1 = nmod.current(s0)[1].copy()
    np.savez(nd / "nmf_20.npz", W=np.load(nd / "nmf_20_project.npz")["W"], H=Hp)  # as before 0.2.63: no key
    for f in (config.SAMPLES / s0.sample_id / "niches").glob("nmf_*.npz"):
        f.unlink()  # (refitted, not read from the sample's copy)
    nmod.run(s0, n_niches=5)
    assert np.allclose(nmod.programs(s0)[1], H1) and (nmod.current(s0)[1] == lab1).all(), "an old project file was used"
    print("one slice's niches after all slices: its own programs, reproduced")
    joint.niches("demo", n_niches=5)  # (all slices' niches again, for the analyses below)
    c = analysis.cells(pid="demo")
    for mode in ("lineage", "niche", "leiden"):
        x = analysis.crosstab(c, "slice", mode)
        print(mode, {r["row"]: (r["n_cells"], {k: r["shares"].get(k, 0) for k in x["keys"][:3]}) for r in x["table"]})
    t = time.time()
    d = analysis.de(c, "niche", c.options()["niche"][0], top=5)
    print("DE of", d["group"], [g["gene"] for g in d["up"]], d["per_slice_check"], f"{time.time() - t:.0f}s")
    m = analysis.markers(c, "leiden", top=3)
    print("top genes per Leiden cluster:", {g["group"]: [x["gene"] for x in g["genes"]] for g in m["groups"][:4]})
    t = time.time()
    r = communication.run(c, "lineage", top=5)
    print("communication:", r["n_pairs_panel"], "pairs in the panel,", r["n_significant"], "significant;",
          [(x["interaction"], x["sender"], x["receiver"]) for x in r["interactions"][:3]], f"{time.time() - t:.0f}s")
    r = communication.compare(c, "slice", "lineage", top=3)
    print("communication by slice:", r["n_significant"], [x["interaction"] for x in r["interactions"]])
finally:
    shutil.rmtree(home)
