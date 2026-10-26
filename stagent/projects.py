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


def slice_name(project, value):
    """Sample id of a section of a project: <project>_<section>, letters, digits, dot, dash, underscore only."""
    return f"{project}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value)).strip('_')}"[:80]
