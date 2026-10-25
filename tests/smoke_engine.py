"""End-to-end smoke test of the engine on one sample with user-drawn regions."""
import json, shutil, time
from pathlib import Path
import numpy as np
from celltypeagent.session import Session, list_samples

t0 = time.time()
print(list_samples().head(3).to_string())
MARKERS = open(Path(__file__).parent.parent / "examples" / "scp2927" / "ibd_xenium480_markers.csv").read()
s = Session.create("147328", region_source="none", markers_text=MARKERS)  # marker annotation (IBD 480-gene list)
print("created", s.id, f"{time.time() - t0:.1f}s")
print(json.dumps(s.summary(), indent=1, default=str)[:600])

# draw two regions around the median cell position (H&E level-0 pixels), like the viewer does
c = s.cells
cx, cy = float(np.median(c.px)), float(np.median(c.py))
s.add_region("patch A", [[cx - 900, cy - 900], [cx + 900, cy - 900], [cx + 900, cy + 900], [cx - 900, cy + 900]])
s.add_region("patch B", [[cx + 1000, cy - 900], [cx + 2600, cy - 900], [cx + 1800, cy + 900]])
print("regions:", [(r["name"], r["n_cells"]) for r in s.regions_info()["regions"]], "| unknown:", s.regions_info()["n_unknown"])

c1 = {"text": "Epithelial cells are only found in patch A", "types": ["epithelial"], "regions": ["patch A"],
      "mode": "only", "weight": "strong"}
print(json.dumps(s.preview_constraint(c1), indent=1, default=str)[:700])
v = s.add_constraint(c1)
for spec in [{"kind": "depth", "text": "Goblet cells stay close to the epithelium", "types": ["Goblet cells"],
              "max_depth_um": 30, "mode": "only", "weight": "moderate"},
             {"kind": "neighborhood", "text": "CD4 T cells sit next to B cells", "types": ["CD4 T cells"],
              "neighbor_types": ["B cells"], "min_frac": 0.2, "mode": "only", "weight": "weak"}]:
    p = s.preview_constraint(spec)
    print(p["interpretation"], "| changes:", p["n_label_changes"])
print(s.conflict_summary().head(5).to_string())
try:
    s.preview_constraint({"text": "typo", "types": ["Gobelt cells"], "regions": ["patch A"], "mode": "only"})
except Exception as e:
    print("expected error:", e)
v_before = s.versions()[-1]
before = s.result()["pred"].copy()
v = s.delete_region("R1")  # patch A is used by C1: the prior keeps applying to its cells
print("delete region:", s.state(v)["message"], "| labels kept:", bool((s.result(v)["pred"] == before).all()))
s.revert(v_before)
s.rename_region("patch B", "patch B2")
cs = s.conflict_summary()
if len(cs):
    png, q = s.fig_conflicts(cs.group_id.iloc[0], n=3)
    print("conflict sheet", png, len(q))
print("map", s.fig_region_map(conflict_group=""))
v = s.label_cells(s.types[int(np.bincount(s.result()["pred"][s.result()["region"] == "patch A"]).argmax())],
                  region="patch A", text="smoke test")
print("label_cells in region:", s.diff(v - 1, v)["n_label_changes"], "changes")
s.record_verdicts([int(c.cell.iloc[0])], "reference_right", note="smoke test")
from celltypeagent import plausibility, markerqc
ps = plausibility.summary(s, top=3)
print("plausibility:", ps["n_out_of_place_cells"], "out-of-place cells,", ps["n_candidate_areas"], "candidate areas; grid", ps["grid"])
print("overview", plausibility.fig_overview(s, highlight=["Epithelial"]))
if ps["candidate_areas"]:
    print("area", plausibility.fig_area(s, ps["candidate_areas"][0]["id"])[0])
    sg = plausibility.suggest(s, [{"area": ps["candidate_areas"][0]["id"], "title": "smoke test", "region_name": "smoke"}])
    print("suggestion", sg[0]["id"], sg[0]["tiles"])
print("area by tiles", plausibility.area_info(s, "B2:C3")["tiles"])
q = markerqc.check(s)
print("marker check:", len(q["nodes"]), "nodes,", len(q["flags"]), "flags; unassigned", (q["unassigned"] or {}).get("n_cells"))
print("proposal diff:", markerqc.propose(s, s.scheme["text"].replace("Goblet cells: ", "Goblet cells: MUC2+ ", 1)
                                         if "MUC2" in s.panel else s.scheme["text"], "smoke test")["diff"])
from celltypeagent import insights, advice
comp = insights.compartment(s, ["Epithelial"])
print("compartment:", comp["n_pieces"], "pieces,", comp["n_cells"], "cells; fig", insights.fig_compartment(s, comp, "mucosa"))
advice.store_compartment(s, "mucosa", comp, {"lineages": ["Epithelial"]})
print("confusion:", insights.type_confusion(s, 3)["pairs"][0]["types"])
print("structure:", len(insights.type_structure(s, "Smooth muscle cells")["subclusters"]), "sub-clusters")
print("profile:", len(insights.spatial_profile(s)["types"]), "types; qc hotspots", insights.qc_overview(s)["hotspot_tiles"]["low_density"][:3])
reg = insights.registration(s)
print("registration:", reg["offset_um"], (reg["affine_fit"] or {}).get("max_offset_um"))
new = advice.add(s, [{"category": "compartments", "title": "mucosa", "text": "smoke test",
                      "action": {"kind": "region", "compartment": "mucosa", "name": "mucosa_smoke",
                                 "prior": {"kind": "region", "text": "Epithelial only in mucosa_smoke", "types": ["Epithelial"],
                                           "regions": ["mucosa_smoke"], "mode": "only", "weight": "strong"}}},
                     {"category": "qc", "title": "registration", "text": "smoke test", "action": {"kind": "registration"}}])
print("advice:", [(x["id"], x["action"]["kind"]) for x in new], advice.status(s))
print("report:", advice.write_report(s, "# Report\n\n" + json.dumps(insights.annotation_facts(s), default=str)[:400])["path"])
v = s.revert(2)
print([(h["version"], h["action"]) for h in s.history()])
print("export", s.export())
print(f"total {time.time() - t0:.1f}s")
shutil.rmtree(s.dir)  # keep work/sessions free of test sessions
