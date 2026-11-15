"""Projects renamed (the Change dialog of the Data tab): the name changes, the id and slices stay; empty, too long
and duplicate names are refused, nothing is written then."""
import json, os, sys, tempfile
home = tempfile.mkdtemp()
os.environ["STAGENT_HOME"] = home; os.environ["XDG_CONFIG_HOME"] = home + "/cfg"
from stagent import config, projects

for sid in ("a1", "a2", "b1"):
    (config.SAMPLES / sid).mkdir(parents=True)
    json.dump({"sample_id": sid, "n_cells": 10}, open(config.SAMPLES / sid / "meta.json", "w"))
projects.save("pa", "Project A", [("a1", "s1", {}), ("a2", "s2", {})])
projects.save("pb", "Project B", [("b1", "s1", {})])
r = {p["project_id"]: p for p in projects.rename_projects({"pa": "  DSS time course "})}
assert r["pa"]["name"] == "DSS time course" and r["pa"]["samples"] == ["a1", "a2"] and r["pb"]["name"] == "Project B"
for bad in ({"pa": ""}, {"pa": "x" * 81}, {"pa": "Project B"}, {"zz": "new"}):
    try:
        projects.rename_projects(bad); raise SystemExit(f"accepted {bad}")
    except ValueError:
        pass
assert projects.load("pa")["name"] == "DSS time course"
r = {p["project_id"]: p for p in projects.rename_projects({"pa": "Project B", "pb": "Project A"})}  # a swap is fine
assert (r["pa"]["name"], r["pb"]["name"]) == ("Project B", "Project A")
print("ok")
