"""Changes Claude proposes to the data, waiting for the user's decision in the app.

A tool of the agent that changes data (labels, priors, regions' cells, clusters, niches, alignment, new samples or
annotations) does not change anything itself: it stores a proposal here, with a description of what it would do.
The app shows each open proposal in a review window; Accept applies it (the app runs the tool with the stored
arguments, see webapp.api_proposal_decide), Dismiss drops it. The agent's tool waits a little for the decision and
returns the outcome, or tells the agent that the proposal is still waiting (pending_changes reads it later).

One JSON file per proposal in proposals/ of the data folder; the agent (MCP server) and the app may be separate
processes. status: open -> applying -> done | failed, or open -> dismissed.
"""
import json, os, time, uuid
from . import config

DIR = config.HOME / "proposals"
FINAL = ("done", "failed", "dismissed")


def _path(pid):
    if not pid or any(ch in pid for ch in "/\\."):
        raise ValueError(f"no proposal {pid!r}")
    return DIR / f"{pid}.json"


def _write(rec):
    DIR.mkdir(parents=True, exist_ok=True)
    tmp = DIR / f".{rec['id']}.{os.getpid()}.tmp"
    tmp.write_text(json.dumps(rec, default=str, ensure_ascii=False, indent=1))
    os.replace(tmp, _path(rec["id"]))  # readers never see half a file


def create(**rec):
    rec = {"id": f"P{time.strftime('%m%d%H%M%S')}{uuid.uuid4().hex[:4]}", "status": "open", "created": time.time()} | rec
    _write(rec)
    return rec


def get(pid):
    p = _path(pid)
    if not p.exists():
        raise ValueError(f"no proposal {pid!r}")
    return json.loads(p.read_text())


def update(pid, **kw):
    rec = get(pid) | kw
    _write(rec)
    return rec


def all(status=None, limit=50):
    """Proposals, oldest first (status: one status or a tuple of them)."""
    out = []
    for p in sorted(DIR.glob("P*.json")) if DIR.exists() else []:
        try:
            r = json.loads(p.read_text())
        except (OSError, ValueError):
            continue
        if status is None or r["status"] == status or (isinstance(status, tuple) and r["status"] in status):
            out.append(r)
    out.sort(key=lambda r: r["created"])
    return out[-limit:]


def wait(pid, seconds):
    """The proposal once the user decided and it was applied (or dismissed), or as it is after `seconds`."""
    end = time.time() + seconds
    while True:
        rec = get(pid)
        if rec["status"] in FINAL or time.time() >= end:
            return rec
        time.sleep(0.5)
