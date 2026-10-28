"""Projects: several slices (samples) analyzed together.

Every slice stays a sample of its own (its cells, image and registration); a project lists them with what is known
of each (the section value and the metadata that is constant within it, e.g. condition, day, mouse), so the viewer
can step through the slices, show them side by side, and later annotate and compare them jointly.

$STAGENT_HOME/projects/<project_id>.json:
  {"project_id", "name", "samples": [sample ids in order], "slices": {sample_id: {"section", "metadata"}}, "created_at"}
"""
import json, re, time
from . import config

PROJECTS = config.HOME / "projects"


def _file(pid):
    return PROJECTS / f"{pid}.json"


def list_projects():
    if not PROJECTS.exists():
        return []
    out = []
    for f in sorted(PROJECTS.glob("*.json")):
        p = json.load(open(f))
        p["samples"] = [s for s in p["samples"] if (config.SAMPLES / s / "meta.json").exists()]  # removed slices drop out
        if p["samples"]:
            out.append(p)
    return out


def load(pid):
    if not _file(pid).exists():
        raise ValueError(f"no project {pid!r}")
    return json.load(open(_file(pid)))


def project_of(sample_id):
    m = config.SAMPLES / sample_id / "meta.json"
    return json.load(open(m)).get("project") if m.exists() else None


def save(pid, name, slices):
    """Create or extend a project; slices = [(sample_id, section, metadata)] in order."""
    PROJECTS.mkdir(parents=True, exist_ok=True)
    p = json.load(open(_file(pid))) if _file(pid).exists() else {
        "project_id": pid, "name": name, "samples": [], "slices": {}, "created_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    for sid, section, meta in slices:
        if sid not in p["samples"]:
            p["samples"].append(sid)
        p["slices"][sid] = {"section": section, "metadata": meta}
        mf = config.SAMPLES / sid / "meta.json"
        m = json.load(open(mf))
        m.update(project=pid, project_section=section, metadata=meta)
        json.dump(m, open(mf, "w"), indent=1)
    json.dump(p, open(_file(pid), "w"), indent=1)
    return p


def drop_sample(sample_id, trash):
    """Take a removed sample out of its project; a project left without slices goes to trash with its folder."""
    pid = project_of(sample_id)
    if not pid or not _file(pid).exists():
        return None
    import shutil
    p = load(pid)
    p["samples"] = [s for s in p["samples"] if s != sample_id]
    p["slices"].pop(sample_id, None)
    if p["samples"]:
        json.dump(p, open(_file(pid), "w"), indent=1)
        return {"project": pid, "left": len(p["samples"])}
    trash.mkdir(parents=True, exist_ok=True)
    shutil.move(str(_file(pid)), str(trash / f"project_{pid}.json"))
    if (PROJECTS / pid).exists():
        shutil.move(str(PROJECTS / pid), str(trash / f"project_{pid}"))
    return {"project": pid, "left": 0}


def slice_name(project, value):
    """Sample id of a section of a project: <project>_<section>, letters, digits, dot, dash, underscore only."""
    return f"{project}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value)).strip('_')}"[:80]


def rename(sample_id, name, project_name=None):
    """Rename a sample as the app shows it; its id (folder, annotations) stays. A slice of a project gets a new section
    name (its name in the project, on the map and in the analyses); a sample on its own a display name ("" for its
    id). project_name renames the slice's project."""
    mf = config.SAMPLES / sample_id / "meta.json"
    if not mf.exists():
        raise ValueError(f"no sample {sample_id!r}")
    m = json.load(open(mf))
    name = str(name or "").strip()
    if len(name) > 80:
        raise ValueError("the name is longer than 80 characters")
    pid = m.get("project")
    p = load(pid) if pid and _file(pid).exists() else None
    if p:
        if not name:
            raise ValueError("a slice needs a name")
        if any(v["section"] == name for k, v in p["slices"].items() if k != sample_id):
            raise ValueError(f"another slice of project {p['name']} is called {name!r}")
        p["slices"].setdefault(sample_id, {"metadata": {}})["section"] = name
        m["project_section"] = name
        if project_name is not None:
            if not str(project_name).strip():
                raise ValueError("the project needs a name")
            p["name"] = str(project_name).strip()[:80]
        json.dump(p, open(_file(pid), "w"), indent=1)
    else:
        taken = set()
        for d in config.SAMPLES.iterdir():
            if d.name != sample_id and (d / "meta.json").exists():
                o = json.load(open(d / "meta.json"))
                taken |= {d.name, o.get("name") or d.name}
        if name and name in taken:
            raise ValueError(f"another sample is called {name!r}")
        if name and name != sample_id:
            m["name"] = name
        else:
            m.pop("name", None)
    json.dump(m, open(mf, "w"), indent=1)
    return {"sample_id": sample_id, "name": m.get("name") or m.get("project_section") or sample_id,
            "project_name": p["name"] if p else None}
