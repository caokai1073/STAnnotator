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

CLAUDE_PROJECT = Path.home() / ".claude" / "projects" / str(config.ROOT).replace("/", "-")
CARRIED = ("project", "project_section", "name", "organ")  # what a replacement keeps of the old sample's meta.json


def sessions_of(sample):
    out = []
    for d in sorted(config.SESSIONS.iterdir()) if config.SESSIONS.exists() else []:
        try:
            if (d / "session.json").exists() and json.load(open(d / "session.json")).get("sample_id") == sample:
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
                    ids.add(json.load(open(f)).get("claude_session"))
                except (OSError, ValueError):
                    pass
        for cid in filter(None, ids):
            for t in (CLAUDE_PROJECT / f"{cid}.jsonl", CLAUDE_PROJECT / cid):
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
    """The sample, its annotations and their chats to the trash. Returns what went where (removed.json)."""
    from . import render
    d = config.SAMPLES / sample
    if not (d / "meta.json").exists():
        raise ValueError(f"no sample {sample!r}")
    meta = json.load(open(d / "meta.json"))
    trash = config.TRASH / f"{time.strftime('%Y%m%d-%H%M%S')}_{sample}_{why}"
    trash.mkdir(parents=True)
    sess = sessions_of(sample)
    n_chats, moved = trash_sessions(sess, trash)
    pid = projects.project_of(sample)
    joint = None
    if keep_project:
        if pid and (projects.PROJECTS / pid).exists():  # joint results made with the old cells
            joint = f"project_{pid}"
            shutil.move(str(projects.PROJECTS / pid), str(trash / joint))
        proj = {"project": pid, "kept": True} if pid else None
    else:
        proj = projects.drop_sample(sample, trash)
    render.forget(d)  # its image (Windows cannot move an open file)
    shutil.move(str(d), str(trash / "sample"))
    info = {"sample": sample, "trash": str(trash), "sessions": [x.name for x in sess], "transcripts": moved,
            "n_chats": n_chats, "project": proj, "joint": joint, "pid": pid, "keep_project": keep_project,
            "carried": {k: meta[k] for k in CARRIED if k in meta}, "time": time.strftime("%Y-%m-%dT%H:%M:%S")}
    json.dump(info, open(trash / "removed.json", "w"), indent=1)
    return info


def _in_trash(info):
    t = Path(info["trash"]).resolve()
    if config.TRASH.resolve() not in t.parents or not (t / "removed.json").exists():
        raise ValueError("not a removed sample in the data folder's trash")
    return t


def restore(info):
    """Put a removed sample back as it was (its folder, annotations, chats, a project's joint results)."""
    from . import render
    t = _in_trash(info)
    sample = info["sample"]
    d = config.SAMPLES / sample
    if d.exists():  # what a failed replacement wrote (normally its import took it away already)
        render.forget(d)
        shutil.move(str(d), str(t / f"sample_new_{os.getpid()}"))
    shutil.move(str(t / "sample"), str(d))
    config.SESSIONS.mkdir(parents=True, exist_ok=True)
    for s in info["sessions"]:
        if (t / s).exists():
            shutil.move(str(t / s), str(config.SESSIONS / s))
    for name in info.get("transcripts", []):
        if (t / "claude_transcripts" / name).exists():
            CLAUDE_PROJECT.mkdir(parents=True, exist_ok=True)
            shutil.move(str(t / "claude_transcripts" / name), str(CLAUDE_PROJECT / name))
    if info.get("joint") and (t / info["joint"]).exists() and info.get("pid"):
        shutil.move(str(t / info["joint"]), str(projects.PROJECTS / info["pid"]))
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
    m = json.load(open(f))
    m.update(info["carried"])
    tmp = f.with_name(f".meta.{os.getpid()}.tmp")
    json.dump(m, open(tmp, "w"), indent=1)
    os.replace(tmp, f)
