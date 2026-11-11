"""Removing one sample to the data folder's trash (recoverable), putting it back, and replacing it with a new import.

remove(): the sample's folder (cells, counts, placement, PCA / UMAP, ...), its annotations and their Claude Code
transcripts go to trash/<time>_<sample>_removed/ (or _replaced/), with removed.json saying what went where; a project
keeps its other slices (keep_project: the sample stays a slice of its project, for a replacement; the project's joint
results, made with the old cells, go to the trash). restore(): puts it all back (a replacement that failed or was
stopped). finish_replace(): the new import takes the old sample's place in its project and its display name.
The viewer drops what it keeps in memory itself (webapp); the agent asks a running viewer to do it (mcp_server).
"""
import json, os, shutil, time
from pathlib import Path
from . import config, projects
from .osutil import dump_json

def claude_project():
    """Claude Code's folder of the chat box's conversations (its working folder is the app's root): named as Claude Code
    names it, every character but letters and digits made '-' (report.project_dir), also on Windows and in an app
    bundle."""
    from .report import project_dir
    return project_dir(config.ROOT)
CARRIED = ("project", "project_section", "name", "organ")  # what a replacement keeps of the old sample's meta.json


def sessions_of(sample):
    out = []
    for d in sorted(config.SESSIONS.iterdir()) if config.SESSIONS.exists() else []:
        try:
            if (d / "session.json").exists() and json.load(open(d / "session.json", encoding="utf-8")).get("sample_id") == sample:
                out.append(d)
        except (OSError, ValueError):
            continue
    return out


def trash_sessions(sess, trash):
    """Move sessions and their Claude Code transcripts into trash (None: delete them for good). Returns the number of
    transcripts and the transcript files moved."""
    n_chats, moved = 0, []
    for d in sess:
        ids = set()
        for f in [d / "chat_meta.json", *(d / "chat_archive").glob("*_chat_meta.json")]:
            if f.exists():
                try:
                    ids.add(json.load(open(f, encoding="utf-8")).get("claude_session"))
                except (OSError, ValueError):
                    pass
        for cid in filter(None, ids):
            for t in (claude_project() / f"{cid}.jsonl", claude_project() / cid):
                if t.exists():
                    n_chats += t.suffix == ".jsonl"
                    if trash is None:
                        shutil.rmtree(t) if t.is_dir() else t.unlink()
                        continue
                    (trash / "claude_transcripts").mkdir(exist_ok=True)
                    shutil.move(str(t), str(trash / "claude_transcripts" / t.name))
                    moved.append(t.name)
        shutil.rmtree(d) if trash is None else shutil.move(str(d), str(trash / d.name))
    return n_chats, moved


def remove(sample, keep_project=False, why="removed"):
    """The sample, its annotations and their chats to the trash. Returns what went where (removed.json). removed.json
    is written before anything moves; when a move fails (e.g. Windows: a file another program has open), what moved
    already is put back and the error raised."""
    from . import render
    d = config.SAMPLES / sample
    if not (d / "meta.json").exists():
        raise ValueError(f"no sample {sample!r}")
    meta = json.load(open(d / "meta.json", encoding="utf-8"))
    trash = config.TRASH / f"{time.strftime('%Y%m%d-%H%M%S')}_{sample}_{why}"
    trash.mkdir(parents=True)
    sess = sessions_of(sample)
    pid = projects.project_of(sample)
    member = None  # its place in its project, to put it back there
    if pid and not keep_project:
        try:
            p = projects.load(pid)
            member = {"index": p["samples"].index(sample), "slice": p["slices"].get(sample)}
        except (OSError, ValueError, KeyError):
            pass
    joint = f"project_{pid}" if keep_project and pid and (projects.PROJECTS / pid).exists() else None
    info = {"sample": sample, "trash": str(trash), "sessions": [x.name for x in sess], "transcripts": [],
            "n_chats": 0, "project": {"project": pid, "kept": True} if keep_project and pid else None, "joint": joint,
            "pid": pid, "member": member, "keep_project": keep_project,
            "carried": {k: meta[k] for k in CARRIED if k in meta}, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    dump_json(info, trash / "removed.json", indent=1)
    try:
        info["n_chats"], info["transcripts"] = trash_sessions(sess, trash)
        if joint:  # joint results made with the old cells
            shutil.move(str(projects.PROJECTS / pid), str(trash / joint))
        if not keep_project:
            info["project"] = projects.drop_sample(sample, trash)
        render.forget(d)  # its image (Windows cannot move an open file)
        shutil.move(str(d), str(trash / "sample"))
    except BaseException:
        try:
            restore(info)
        except BaseException as e2:
            print(f"remove {sample}: could not put back what moved ({e2}); it is in {trash}")
        raise
    dump_json(info, trash / "removed.json", indent=1)
    return info


def _plain(name):
    """A file or folder name of its own (no path), or None."""
    if name is None:
        return None
    name = str(name)
    if not name or name in (".", "..") or "/" in name or "\\" in name or Path(name).name != name:
        raise ValueError(f"not a plain name: {name!r}")
    return name


def _in_trash(info):
    t = Path(info["trash"]).resolve()
    if config.TRASH.resolve() not in t.parents or not (t / "removed.json").exists():
        raise ValueError("not a removed sample in the data folder's trash")
    return t


def restore(info):
    """Put a removed sample back as it was (its folder, annotations, chats, its project, a project's joint results).
    What goes where is read from the trash folder's removed.json (info only names that folder), also when remove()
    failed half way and only part of it moved."""
    from . import render
    t = _in_trash(info)
    info = json.load(open(t / "removed.json", encoding="utf-8"))
    sample, pid = _plain(info["sample"]), _plain(info.get("pid"))
    d = config.SAMPLES / sample
    if (t / "sample").exists():
        if d.exists():  # what a failed replacement wrote (normally its import took it away already)
            render.forget(d)
            shutil.move(str(d), str(t / f"sample_new_{os.getpid()}"))
        shutil.move(str(t / "sample"), str(d))
    config.SESSIONS.mkdir(parents=True, exist_ok=True)
    for s in map(_plain, info["sessions"]):
        if (t / s).exists():
            shutil.move(str(t / s), str(config.SESSIONS / s))
    ct = t / "claude_transcripts"
    if ct.exists():  # every transcript moved (also those of a remove() that failed before listing them)
        cp = claude_project()
        cp.mkdir(parents=True, exist_ok=True)
        for f in list(ct.iterdir()):
            shutil.move(str(f), str(cp / f.name))
        ct.rmdir()
    if info.get("joint") and pid and (t / _plain(info["joint"])).exists():
        shutil.move(str(t / info["joint"]), str(projects.PROJECTS / pid))
    if pid and not info.get("keep_project"):  # back into its project
        if (t / f"project_{pid}.json").exists():  # it was the project's last slice
            shutil.move(str(t / f"project_{pid}.json"), str(projects.PROJECTS / f"{pid}.json"))
            if (t / f"project_{pid}").exists():
                shutil.move(str(t / f"project_{pid}"), str(projects.PROJECTS / pid))
        elif info.get("member") and (projects.PROJECTS / f"{pid}.json").exists():
            p, m = projects.load(pid), info["member"]
            if sample not in p["samples"]:
                p["samples"].insert(min(int(m.get("index", len(p["samples"]))), len(p["samples"])), sample)
                p["slices"].setdefault(sample, m.get("slice") or {})
                dump_json(p, projects.PROJECTS / f"{pid}.json", indent=1)
    (t / "removed.json").rename(t / "restored.json")
    left = [x for x in t.iterdir() if x.name != "restored.json"]
    if not left:  # all of it went back
        shutil.rmtree(t)
    return {"restored": sample}


def finish_replace(info):
    """The new import of a replaced sample keeps the old one's project, slice name, display name and organ."""
    f = config.SAMPLES / info["sample"] / "meta.json"
    if not f.exists() or not info.get("carried"):
        return
    m = json.load(open(f, encoding="utf-8"))
    m.update(info["carried"])
    dump_json(m, f, indent=1)  # (atomic, with a temporary name of its own)
