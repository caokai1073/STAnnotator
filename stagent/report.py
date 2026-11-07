"""A problem report: what the user did with the agent and what the app logged, as one zip the user chooses to send.

Report bug... (right of Export in the toolbar) saves it. It holds text only, no expression data and no images:
- README.txt: what the report holds, the version and the system; note.txt: what the user wrote;
- about.json: the About information (version, build, data folder, Claude Code found / connected, libraries);
- agent/<id>.jsonl: the Claude Code conversations of the agent (Ask Agent's workspace, the app's chat box) touched
  in the last days, with images left out;
- agent_calls.jsonl: every call of the app's agent tools in those days (tool, arguments, time, result or error),
  written by log_call;
- proposals/: the changes the agent proposed in those days and what became of them (Accept, Dismiss, error);
- annotation/: the open annotation's history (versions, priors, corrections, analyses run, chat), no cell tables;
- app.log: the end of the app's log.
The user sees what goes in before saving; nothing is sent by the app.
"""
import json
import os
import platform
import re
import time
import zipfile
from pathlib import Path

from . import config

CALLS = config.HOME / "agent_calls.jsonl"
CALLS_MAX = 8 * 2 ** 20      # agent_calls.jsonl keeps its newer half past this size
ARG_MAX = 4000               # characters of one argument / result kept in the call log
LOG_TAIL = 2 * 2 ** 20       # the end of app.log that goes into a report
TRANSCRIPT_MAX = 40 * 2 ** 20  # all conversations together; the newest first
SESSION_FILE_MAX = 2 * 2 ** 20
SESSION_TYPES = {".json", ".jsonl", ".md", ".txt"}


# ---------------- the agent's tool calls ----------------
def _short(v):
    s = v if isinstance(v, str) else json.dumps(v, default=str, ensure_ascii=False)
    return s if len(s) <= ARG_MAX else s[:ARG_MAX] + f"… ({len(s)} characters)"


TRACE_MAX = 12000


def log_call(tool, args, seconds, result=None, error=None, session=None, via=None, trace=None):
    """One line per call of an agent tool (mcp_server wraps every tool with it), with the traceback when it failed:
    the tools run in the process Claude Code starts, whose own output reaches no log. Never raises."""
    try:
        rec = {"time": time.strftime("%Y-%m-%dT%H:%M:%S"), "tool": tool, "seconds": round(seconds, 2),
               "session": session, "via": via, "args": {k: _short(v) for k, v in (args or {}).items()}}
        if error is not None:
            rec["error"] = _short(error)
            if trace:
                rec["traceback"] = trace[-TRACE_MAX:]
        else:
            rec["result"] = _short(result)
        CALLS.parent.mkdir(parents=True, exist_ok=True)
        with open(CALLS, "a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
        if CALLS.stat().st_size > CALLS_MAX:
            lines = CALLS.read_text(encoding="utf-8").splitlines(keepends=True)
            CALLS.write_text("".join(lines[len(lines) // 2:]), encoding="utf-8")
    except Exception:
        pass


def _recent_calls(since):
    if not CALLS.exists():
        return []
    stamp = time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(since))
    out = []
    for line in CALLS.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            if json.loads(line).get("time", "") >= stamp:
                out.append(line)
        except ValueError:
            continue
    return out


# ---------------- the Claude Code conversations ----------------
def _claude_dir():
    return Path(os.environ.get("CLAUDE_CONFIG_DIR") or Path.home() / ".claude")


def project_dir(folder):
    """Where Claude Code keeps the conversations it had in this working folder: the folder's path with every
    character other than a letter or digit made '-' (long names are cut at 200 characters and given a hash)."""
    name = re.sub(r"[^A-Za-z0-9]", "-", str(Path(folder).resolve()))
    base = _claude_dir() / "projects"
    if len(name) <= 200:
        return base / name
    hits = sorted(base.glob(name[:200] + "*")) if base.exists() else []
    return hits[0] if hits else base / name


def _agent_folders():
    """The working folders the app starts Claude Code in: Ask Agent's workspace (one per server) and the chat box."""
    return [*sorted(config.HOME.glob("claude_workspace*")), config.ROOT]


def _chat_box_ids():
    """The conversations of the app's chat box (their ids are kept with each annotation)."""
    ids = set()
    if config.SESSIONS.is_dir():
        for f in [*config.SESSIONS.glob("*/chat_meta.json"), *config.SESSIONS.glob("*/chat_archive/*_chat_meta.json")]:
            try:
                ids.add(json.loads(f.read_text(encoding="utf-8")).get("claude_session"))
            except (OSError, ValueError):
                continue
    return {i for i in ids if i}


def transcripts(since):
    """Conversation files touched since then, newest first: (path, name in the report). Ask Agent's workspaces hold
    only the app's conversations, so all of them count; the chat box runs in the app's own folder, which (run from
    source) may hold other Claude Code work, so only its own conversations count there."""
    out = []
    for folder in _agent_folders():
        d = project_dir(folder)
        if not d.is_dir():
            continue
        box = Path(folder) == config.ROOT
        ids = _chat_box_ids() if box else None
        tag = "chat_box" if box else Path(folder).name
        for f in d.rglob("*.jsonl"):  # sub-agents' conversations sit in <id>/subagents/
            rel = f.relative_to(d)
            if box and rel.parts[0].removesuffix(".jsonl") not in ids:
                continue
            try:
                if f.stat().st_mtime >= since:
                    out.append((f, f"agent/{tag}/{rel.as_posix()}"))
            except OSError:
                continue
    return sorted(out, key=lambda x: -x[0].stat().st_mtime)


def _drop_images(o):
    """Images the tools showed the agent (figures, H&E views) are large and not needed to follow a conversation."""
    if isinstance(o, dict):
        if o.get("type") == "image" and isinstance(o.get("source"), dict):
            return {"type": "image", "source": {"type": "omitted from the report"}}
        return {k: _drop_images(v) for k, v in o.items()}
    if isinstance(o, list):
        return [_drop_images(v) for v in o]
    return o


def _clean_transcript(path):
    lines = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        try:
            lines.append(json.dumps(_drop_images(json.loads(line)), ensure_ascii=False))
        except ValueError:
            lines.append(line)
    return "\n".join(lines) + "\n"


# ---------------- the rest ----------------
def app_logs():
    """The logs and their names in the report: app.log (the data folder when run from source, else the built app's
    own folder), viewer.log (the viewer in a browser, or started by the agent) and remote_viewer.log (the viewer this
    computer serves to the desktop app as a server)."""
    out = {}
    for p in (config.HOME / "app.log", Path.home() / "Library" / "Caches" / "STAnnotator" / "app.log",
              Path(os.environ.get("LOCALAPPDATA") or Path.home() / "AppData" / "Local") / "STAnnotator" / "app.log"):
        if p.is_file() and "app.log" not in out.values():
            out[p] = "app.log"
    for n in ("viewer.log", "remote_viewer.log"):
        if (config.HOME / n).is_file():
            out[config.HOME / n] = n
    return list(out.items())


def _proposals(since):
    d = config.HOME / "proposals"
    return sorted((p for p in d.glob("*.json") if p.stat().st_mtime >= since), key=lambda p: p.stat().st_mtime) \
        if d.is_dir() else []


def _session_files(sid):
    if not sid or not re.fullmatch(r"[A-Za-z0-9_.-]+", sid):
        return []
    d = config.SESSIONS / sid
    if not d.is_dir():
        return []
    fs = [f for f in d.rglob("*") if f.is_file() and f.suffix in SESSION_TYPES and f.stat().st_size <= SESSION_FILE_MAX
          and "cache" not in f.relative_to(d).parts]
    return sorted(fs)


def contents(days=7, sid=None):
    """What a report with these settings would hold: counts and sizes per part, for the dialog."""
    since = time.time() - days * 86400
    tr = transcripts(since)
    props = _proposals(since)
    sess = _session_files(sid)
    logs = [p for p, _ in app_logs()]
    return {"days": days, "session": sid if sess else None,
            "conversations": {"files": len([1 for _, n in tr if "/subagents/" not in n]), "bytes": sum(p.stat().st_size for p, _ in tr),
                              "folders": [str(project_dir(f)) for f in _agent_folders() if project_dir(f).is_dir()]},
            "calls": {"n": len(_recent_calls(since))},
            "proposals": {"n": len(props)},
            "annotation": {"files": len(sess)},
            "app_log": {"files": [str(p) for p in logs], "bytes": sum(min(p.stat().st_size, LOG_TAIL) for p in logs)}}


def _readme(about, days, parts, sid):
    lines = ["STAnnotator problem report", "",
             f"Made {time.strftime('%Y-%m-%d %H:%M')} with STAnnotator {about.get('version')} (build {about.get('build')})",
             f"System: {platform.system()} {platform.release()} ({platform.machine()}), Python {platform.python_version()}",
             f"Covers the last {days} day(s){f'; annotation {sid}' if sid else ''}", "", "Holds:"]
    what = {"note": "note.txt: what the user wrote", "about": "about.json: the About information",
            "conversations": "agent/: the agent's Claude Code conversations (JSON lines; images left out)",
            "calls": "agent_calls.jsonl: the app's agent tools called (arguments, results or errors with tracebacks, seconds)",
            "proposals": "proposals/: changes the agent proposed and what became of them",
            "annotation": "annotation/: the annotation's history (versions, priors, corrections, analyses, chat); "
                          "sample/meta.json: its sample (platform, cells, genes, image, registration)",
            "page_errors": "page_errors.json: errors in the app's window (scripts, failed requests)",
            "app_log": "logs/: the end of the app's logs (app.log, viewer.log, remote_viewer.log)"}
    lines += [f"- {what[p]}" for p in what if p in parts]
    lines += ["", "No expression data, cell tables or images are included."]
    return "\n".join(lines) + "\n"


def _sample_meta(sid):
    """The open annotation's sample: platform, cells, genes, image and its scale, how it was registered (no data)."""
    try:
        sample = json.loads((config.SESSIONS / sid / "session.json").read_text(encoding="utf-8"))["sample_id"]
    except (OSError, ValueError, KeyError, TypeError):
        return None
    f = config.SAMPLES / str(sample) / "meta.json"
    return f if re.fullmatch(r"[A-Za-z0-9_.-]+", str(sample)) and f.is_file() else None


def _page_errors(errors):
    """Errors the app's window noted (scripts, failed requests): a list of small records, kept as they came."""
    if not isinstance(errors, list):
        return []
    out = []
    for e in errors[-200:]:
        if isinstance(e, dict):
            out.append({str(k)[:40]: _short(v) for k, v in list(e.items())[:12]})
    return out


def build(dest, about, days=7, parts=("conversations", "calls", "proposals", "annotation", "app_log"), sid=None,
          note="", page_errors=None):
    """Write the report to dest (a .zip); returns what went in."""
    since = time.time() - days * 86400
    parts = set(parts) | {"about"} | ({"note"} if note.strip() else set())
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    n = {"conversations": 0, "calls": 0, "proposals": 0, "annotation": 0, "app_log": 0, "page_errors": 0}
    skipped = []
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
        if note.strip():
            z.writestr("note.txt", note.strip() + "\n")
        z.writestr("about.json", json.dumps(about, indent=1, default=str, ensure_ascii=False))
        if "conversations" in parts:
            total = 0
            for p, name in transcripts(since):
                try:
                    text = _clean_transcript(p)
                except OSError:
                    continue
                if total + len(text) > TRANSCRIPT_MAX:
                    skipped.append(name)
                    continue
                total += len(text)
                z.writestr(name, text)
                n["conversations"] += "/subagents/" not in name
        if "calls" in parts:
            calls = _recent_calls(since)
            if calls:
                z.writestr("agent_calls.jsonl", "\n".join(calls) + "\n")
            n["calls"] = len(calls)
        if "proposals" in parts:
            for p in _proposals(since):
                z.write(p, f"proposals/{p.name}")
                n["proposals"] += 1
        if "annotation" in parts:
            d = config.SESSIONS / (sid or "")
            for f in _session_files(sid):
                z.write(f, f"annotation/{f.relative_to(d).as_posix()}")
                n["annotation"] += 1
            meta = _sample_meta(sid) if n["annotation"] else None
            if meta:
                z.write(meta, "sample/meta.json")
        if "page_errors" in parts:
            errs = _page_errors(page_errors)
            if errs:
                z.writestr("page_errors.json", json.dumps(errs, indent=1, ensure_ascii=False))
            n["page_errors"] = len(errs)
        if "app_log" in parts:
            for p, name in app_logs():
                with open(p, "rb") as f:
                    f.seek(max(0, p.stat().st_size - LOG_TAIL))
                    z.writestr(f"logs/{name}", f.read())
                n["app_log"] += 1
        z.writestr("README.txt", _readme(about, days, parts, sid if n["annotation"] else None))
    os.replace(tmp, dest)
    return {"saved": str(dest), "bytes": dest.stat().st_size, "included": n, "left_out": skipped}


def add_local(zip_path, days=7, parts=("conversations", "app_log")):
    """With the app working on a server, the report is made there; the conversations (Claude Code runs on this
    computer) and this computer's app.log are added to it here."""
    since = time.time() - days * 86400
    n = {"conversations": 0, "app_log": 0}
    with zipfile.ZipFile(zip_path, "a", zipfile.ZIP_DEFLATED) as z:
        have = set(z.namelist())
        if "conversations" in parts:
            total = 0
            for p, name in transcripts(since):
                name = "local/" + name
                if name in have:
                    continue
                text = _clean_transcript(p)
                if total + len(text) > TRANSCRIPT_MAX:
                    continue
                total += len(text)
                z.writestr(name, text)
                n["conversations"] += "/subagents/" not in name
        if "app_log" in parts:
            for p, name in app_logs():
                with open(p, "rb") as f:
                    f.seek(max(0, p.stat().st_size - LOG_TAIL))
                    z.writestr(f"local/logs/{name}", f.read())
                n["app_log"] += 1
    return n


# ---------------- upload to the developers (packaging/report_receiver.gs, a Google Apps Script web app) ----------------
UPLOAD_MAX = 40 * 2 ** 20
_URL = re.compile(r"https://script\.google(usercontent)?\.com/macros/s/[A-Za-z0-9_-]+/exec")


BUILT_IN = config.PACKAGE / "data" / "report_upload_url.txt"  # put there before building the app; not in git


def _built_in():
    try:
        u = BUILT_IN.read_text(encoding="utf-8").strip()
        return u if _URL.fullmatch(u) else ""
    except OSError:
        return ""


def upload_url():
    """The link the tester saved, else the one the app was built with."""
    return config.read_settings().get("report_upload_url", "") or _built_in()


def link_status():
    """Whether uploading works without asking for a link (the link itself stays out of the page)."""
    if config.read_settings().get("report_upload_url"):
        return {"set": True, "source": "settings"}
    return {"set": bool(_built_in()), "source": "built-in" if _built_in() else None}


def _check_url(url):
    url = (url or "").strip()
    if not _URL.fullmatch(url):
        raise ValueError("this is not an upload link from the developers (it looks like "
                         "https://script.google.com/macros/s/…/exec)")
    return url


def _request(url, payload=None, timeout=120):
    """GET (payload None) or POST JSON; Apps Script answers through a redirect, which urllib follows."""
    import urllib.request
    data = None if payload is None else json.dumps(payload).encode()
    req = urllib.request.Request(url, data=data, headers={"Content-Type": "application/json"} if data else {})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        text = r.read().decode("utf-8", errors="replace")
    try:
        return json.loads(text)
    except ValueError:
        raise ValueError("the upload link did not answer as expected: check it, or save the report and send it")


def set_upload_url(url):
    """Check the link (the receiver answers a GET) and keep it in the app's settings."""
    url = _check_url(url)
    r = _request(url, timeout=30)
    if not r.get("ok"):
        raise ValueError(r.get("error") or "the upload link did not answer as expected")
    config.write_settings(report_upload_url=url)
    return url


def upload(path, name, about=None, note="", sender="", url=None):
    """Send a saved report to the developers' Drive folder."""
    import base64
    url = _check_url(url or upload_url())
    path = Path(path)
    size = path.stat().st_size
    if size > UPLOAD_MAX:
        raise ValueError(f"the report is {size / 2 ** 20:.0f} MB, more than the upload takes (40 MB): "
                         "choose fewer days or leave the conversations out, or save it and send it")
    about = about or {}
    r = _request(url, {"name": name, "data": base64.b64encode(path.read_bytes()).decode(),
                       "version": about.get("version"), "system": about.get("system"),
                       "from": sender[:200], "note": note[:2000]})
    if not r.get("ok"):
        raise ValueError(r.get("error") or "the upload failed")
    return {"uploaded": r.get("name", name), "bytes": size}
