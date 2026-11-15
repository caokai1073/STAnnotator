"""Projects renamed (the Change dialog of the Data tab): the name changes, the id and slices stay; empty, too long
and duplicate names are refused, nothing is written then. Assigning by name: a name before an id; a project added in
the dialog (new) is always a new project, also when its name is the id of another one. The endpoints refuse bodies
of the wrong shape (400, not 500)."""
import json, os, sys, tempfile
home = tempfile.mkdtemp()
os.environ["STAGENT_HOME"] = home; os.environ["XDG_CONFIG_HOME"] = home + "/cfg"
from stagent import config, projects

for sid in ("a1", "a2", "b1", "c1", "c2"):
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

# a project added in the dialog under the name "pa" (another project's id): a new project, not pa
r = {p["project_id"]: p for p in projects.assign({"c1": "pa", "c2": "pa"}, new=["pa"])}
newp = [pid for pid, p in r.items() if p["name"] == "pa"]
assert len(newp) == 1 and newp[0] != "pa" and r[newp[0]]["samples"] == ["c1", "c2"], r
assert r["pa"]["samples"] == ["a1", "a2"]
# by name: "pa" is now the name of the new project, so it wins over the id pa
r = {p["project_id"]: p for p in projects.assign({"b1": "pa"})}
assert "b1" in r[newp[0]]["samples"] and "b1" not in r["pa"]["samples"]
try:  # a new project under the name of an existing one: refused
    projects.assign({"b1": "Project B"}, new=["Project B"]); raise SystemExit("accepted a second Project B")
except ValueError:
    pass

from starlette.testclient import TestClient
from stagent import webapp
c = TestClient(webapp.build_app(8799), base_url="http://127.0.0.1:8799")
for url, body in (("/api/projects/rename", {"names": ["x"]}), ("/api/projects/rename", ["x"]),
                  ("/api/projects/assign", {"groups": ["c1"]}), ("/api/projects/assign", {"groups": {}, "new": "x"})):
    assert c.post(url, json=body).status_code == 400, (url, body)
assert c.post("/api/sample/c1/remove", json={"check": True}).json() == {"ok": True}
assert (config.SAMPLES / "c1" / "meta.json").exists()  # (a check removes nothing)
print("ok")
