"""Projects: several slices (samples) analyzed together.

Every slice stays a sample of its own (its cells, image and registration); a project lists them with what is known
of each (the section value and the metadata that is constant within it, e.g. condition, day, mouse), so the viewer
can step through the slices, show them side by side, and later annotate and compare them jointly.

$STAGENT_HOME/projects/<project_id>.json:
  {"project_id", "name", "samples": [sample ids in order], "slices": {sample_id: {"section", "metadata"}}, "created_at"}
"""
import json, re, time
from . import config
from .osutil import dump_json

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
        dump_json(m, mf, indent=1)
    dump_json(p, _file(pid), indent=1)
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
        dump_json(p, _file(pid), indent=1)
        return {"project": pid, "left": len(p["samples"])}
    trash.mkdir(parents=True, exist_ok=True)
    shutil.move(str(_file(pid)), str(trash / f"project_{pid}.json"))
    if (PROJECTS / pid).exists():
        shutil.move(str(PROJECTS / pid), str(trash / f"project_{pid}"))
    return {"project": pid, "left": 0}


def assign(groups, reordered=None):
    """Put samples into projects: groups = {sample_id: project name, or "" for none}; samples not named keep their
    project. A project is matched by its name (or id), a new name makes a new project; the slices of a project must
    share one gene panel. A slice that leaves a project keeps its name as the sample's name; a project left without
    slices goes away (its joint results to trash/). A slice with the project's genes in another order is put in the
    project's order (reorder_genes; its id appended to the list reordered). Returns the projects."""
    import shutil
    metas = {}
    for sid in groups:
        mf = config.SAMPLES / sid / "meta.json"
        if not mf.exists():
            raise ValueError(f"no sample {sid!r}")
        metas[sid] = json.load(open(mf))
    projs = {p["project_id"]: p for p in (json.load(open(f)) for f in sorted(PROJECTS.glob("*.json")))} if PROJECTS.exists() else {}
    by_name = {p["name"]: pid for pid, p in projs.items()} | {pid: pid for pid in projs}
    target = {}  # sample -> project id (None: on its own)
    for sid, name in groups.items():
        name = str(name or "").strip()[:80]
        if not name:
            target[sid] = None
            continue
        pid = by_name.get(name)
        if pid is None:
            pid = re.sub(r"[^A-Za-z0-9_.-]+", "_", name).strip("_") or "project"
            while pid in projs or (config.SAMPLES / pid).exists():
                pid += "_p"
            projs[pid] = {"project_id": pid, "name": name, "samples": [], "slices": {},
                          "created_at": time.strftime("%Y-%m-%dT%H:%M:%S")}
            by_name[name] = by_name[pid] = pid
        target[sid] = pid
    for sid, pid in target.items():  # leave the old project
        old = metas[sid].get("project")
        if old and old != pid and old in projs:
            sec = projs[old]["slices"].get(sid, {}).get("section")
            projs[old]["samples"] = [x for x in projs[old]["samples"] if x != sid]
            projs[old]["slices"].pop(sid, None)
            if pid is None and sec and sec != sid and not metas[sid].get("name"):
                metas[sid]["name"] = sec  # keeps the name it had in the project
    for sid, pid in target.items():  # join the new one
        if pid and sid not in projs[pid]["samples"]:
            sec = metas[sid].get("project_section") if metas[sid].get("project") == pid else None
            sec = sec or metas[sid].get("name") or metas[sid].get("project_section") or sid
            taken = {v["section"] for v in projs[pid]["slices"].values()}
            base, i = sec, 2
            while sec in taken:
                sec, i = f"{base}_{i}", i + 1
            projs[pid]["samples"].append(sid)
            projs[pid]["slices"][sid] = {"section": sec, "metadata": metas[sid].get("metadata") or {}}
    reorder = {}  # one gene panel per project: the same genes in another order are put in the project's order
    for pid in set(x for x in target.values() if x):
        reorder.update(panel_plan(projs[pid]["name"], projs[pid]["samples"]))
    for sid, ref in reorder.items():
        reorder_genes(sid, ref)
    if reordered is not None:
        reordered.extend(reorder)
    PROJECTS.mkdir(parents=True, exist_ok=True)
    for pid, p in projs.items():
        if p["samples"]:
            dump_json(p, _file(pid), indent=1)
        elif _file(pid).exists():  # empty now: gone, its joint results kept in the trash
            trash = config.TRASH / f"{time.strftime('%Y%m%d-%H%M%S')}_project_{pid}"
            trash.mkdir(parents=True, exist_ok=True)
            shutil.move(str(_file(pid)), str(trash / f"{pid}.json"))
            if (PROJECTS / pid).exists():
                shutil.move(str(PROJECTS / pid), str(trash / pid))
    for sid, pid in target.items():
        m = metas[sid]
        if pid:
            m.update(project=pid, project_section=projs[pid]["slices"][sid]["section"])
        else:
            m.pop("project", None); m.pop("project_section", None)
        dump_json(m, config.SAMPLES / sid / "meta.json", indent=1)
    return list_projects()


# ---------------- one gene panel per project ----------------
LIST_DIFFERENT = 20  # differing gene names listed in the error below this many


def _annotated(sample_id):
    """The sample has annotations other than its placeholder (their models and analyses use its gene columns)."""
    from .session import list_sessions
    ss = list_sessions()
    return bool(len(ss)) and bool(((ss.sample_id == sample_id) & (ss.method != "none")).any())


def panel_plan(name, sample_ids):
    """The slices of project `name` must share one gene panel, gene for gene in one order (the joint analyses stack
    their counts). Returns {sample_id: the slice whose gene order to put it in} for the slices that have the same genes in another
    order (as CosMx exports can list them); raises ValueError when gene names differ (how many, and which when fewer
    than LIST_DIFFERENT), or when a slice to reorder is annotated in its own order. The order kept is that of the
    annotated slices, else of the project's first slice."""
    panels = {}
    for sid in sample_ids:
        g = config.SAMPLES / sid / "genes.txt"
        if g.exists():
            panels.setdefault(tuple(config.read_genes(g)), []).append(sid)
    if len(panels) < 2:
        return {}
    annotated = {sid for sids in panels.values() for sid in sids if _annotated(sid)}
    keep = [k for k, sids in panels.items() if annotated & set(sids)]
    ref = keep[0] if keep else next(iter(panels))
    rs, ref_set = ", ".join(panels[ref]), set(ref)
    problems, plan = [], {}
    for genes, sids in panels.items():
        if genes == ref:
            continue
        here = ", ".join(sids)
        only_ref = [g for g in ref if g not in set(genes)]
        only_here = [g for g in genes if g not in ref_set]
        if not only_ref and not only_here and len(set(genes)) == len(genes) == len(ref):
            if genes in keep:
                problems.append(f"{here} and {rs} have the same genes in another order, and both are annotated in "
                                "their own order (an annotation's model uses the gene columns): join them without "
                                "these annotations, or annotate the slices again after joining")
            else:
                plan.update({sid: panels[ref][0] for sid in sids})
            continue
        n = len(only_ref) + len(only_here)
        msg = (f"{n} gene names differ between {here} ({len(genes)} genes) and {rs} ({len(ref)} genes)")
        if n < LIST_DIFFERENT:
            parts = ([f"only in {rs}: {', '.join(only_ref)}"] if only_ref else []) + \
                    ([f"only in {here}: {', '.join(only_here)}"] if only_here else [])
            msg += ": " + "; ".join(parts)
        else:
            msg += f" ({len(only_ref)} only in {rs}, {len(only_here)} only in {here})"
        problems.append(msg)
    if problems:
        raise ValueError(f"the slices of project {name} must share one gene panel: " + ". ".join(problems))
    return plan


def reorder_genes(sample_id, like):
    """Put the sample's gene columns in the order of sample `like` (the same genes): counts.npz and genes.txt, and the
    data as imported (original/, kept by a filter) in the order of like's (genes like does not have stay after the
    ones they followed). The counts are the same numbers in other columns."""
    import numpy as np
    import scipy.sparse as sp
    from .osutil import save_counts
    d, L = config.SAMPLES / sample_id, config.SAMPLES / like
    for sub in (d / "original", d):  # (the data as imported first: the current files tell what still needs it)
        if not (sub / "counts.npz").exists():
            continue
        genes = config.read_genes(sub / "genes.txt")
        lf = L / "original" / "genes.txt" if sub.name == "original" and (L / "original" / "genes.txt").exists() \
            else L / "genes.txt"
        pos = {g: i for i, g in enumerate(config.read_genes(lf))}
        # genes like does not have (its filter removed them, this slice's did not) after the ones before them
        order = sorted(range(len(genes)), key=lambda j: (pos[genes[j]], 0) if genes[j] in pos else
                       (max([pos[g] for g in genes[:j] if g in pos], default=-1), 1 + j))
        if order == list(range(len(genes))):
            continue
        inv = np.empty(len(order), np.int64)
        inv[order] = np.arange(len(order))
        X = sp.load_npz(sub / "counts.npz").tocsr()  # (its own copy: the open annotations' shared one stays as it is)
        X.indices = inv[X.indices].astype(X.indices.dtype)
        X.has_sorted_indices = False
        X.sort_indices()
        save_counts(sub / "counts.npz", X)
        tmp = sub / "genes.txt.part"
        config.write_genes(tmp, [genes[j] for j in order])
        tmp.replace(sub / "genes.txt")


def slice_name(project, value):
    """Sample id of a section of a project: <project>_<section>, letters, digits, dot, dash, underscore only."""
    return f"{project}_{re.sub(r'[^A-Za-z0-9_.-]+', '_', str(value)).strip('_')}"  # never shortened: that made slices with a common prefix clash


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
        dump_json(p, _file(pid), indent=1)
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
    dump_json(m, mf, indent=1)
    return {"sample_id": sample_id, "name": m.get("name") or m.get("project_section") or sample_id,
            "project_name": p["name"] if p else None}
