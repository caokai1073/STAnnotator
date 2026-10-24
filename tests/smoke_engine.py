"""End-to-end smoke test of the engine on one sample with user-drawn regions."""
import json, shutil, time
import numpy as np
from celltypeagent.session import Session, list_samples

t0 = time.time()
print(list_samples().head(3).to_string())
s = Session.create("147328", region_source="none")  # marker annotation with the project's draft list
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
try:
    s.delete_region("R1")
except Exception as e:
    print("expected refusal:", e)
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
from celltypeagent import plausibility
ps = plausibility.summary(s, top=3)
print("plausibility:", ps["n_flagged_cells"], "flagged cells,", ps["n_groups"], "groups; cohesive:", ps["cohesive_types"])
if ps["groups"]:
    print("plausibility sheet", s.fig_plausibility([g["id"] for g in ps["groups"]])[0])
    v = plausibility.label_groups(s, [ps["groups"][0]["id"]])
    print("plausibility label:", s.state(v)["message"])
v = s.revert(2)
print([(h["version"], h["action"]) for h in s.history()])
print("export", s.export())
print(f"total {time.time() - t0:.1f}s")
shutil.rmtree(s.dir)  # keep work/sessions free of test sessions
