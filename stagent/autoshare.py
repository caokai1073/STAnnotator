"""Usage sharing during testing: with the user's consent, the desktop app sends the developers a report of its use every
few days, the same report as Report a bug (report.py: conversations with the agent without images, the agent's tool
calls, the changes it proposed, the annotation in use, the app's log; no expression data, cell tables or images).

The user is asked once (the app's first start with this version); About shows the state, turns it off or on, and sends
one now. The receiver files these reports under autosave/<installation id>/ in the developers' Drive folder.
Settings: share_usage ("yes" / "no"; absent: not asked yet), share_last (the last look, sent or with nothing to send),
share_sent (the last report sent).
"""
import tempfile
import threading
import time
from pathlib import Path

from . import config, report

EVERY_DAYS = 3
MAX_DAYS = 7          # a report covers the time since the last one, at most this
CHECK_SECONDS = 6 * 3600
_lock = threading.Lock()


def status():
    st = config.read_settings()
    return {"consent": st.get("share_usage"), "last_sent": st.get("share_sent"), "every_days": EVERY_DAYS,
            "link": report.link_status().get("set", False), "install_id": report.install_id()}


def set_consent(yes):
    config.write_settings(share_usage="yes" if yes else "no")
    return status()


def _active_session(since):
    """The annotation in use: the one the app shows, else the one changed last in this time."""
    from . import claudelink
    sid = claudelink.viewer_state().get("session")
    if sid and (config.SESSIONS / sid).is_dir():
        return sid
    recent = [d for d in config.SESSIONS.glob("*") if d.is_dir() and d.stat().st_mtime >= since] \
        if config.SESSIONS.is_dir() else []
    return max(recent, key=lambda d: d.stat().st_mtime).name if recent else None


def _used_since(since):
    """Whether the app was used in this time: a conversation, an agent tool call, or an annotation changed."""
    if report.transcripts(since) or report._recent_calls(since):
        return True
    return config.SESSIONS.is_dir() and any(d.stat().st_mtime >= since for d in config.SESSIONS.glob("*"))


def send(about, force=False):
    """Make and upload one report if due (or force). Returns what happened; never raises."""
    if not _lock.acquire(blocking=False):
        return {"skipped": "already sending"}
    try:
        st = config.read_settings()
        if st.get("share_usage") != "yes":
            return {"skipped": "not allowed"}
        if not report.upload_url():
            return {"skipped": "no upload link"}
        now = time.time()
        last = float(st.get("share_last") or 0)
        if not force and now - last < EVERY_DAYS * 86400:
            return {"skipped": "not due", "next": last + EVERY_DAYS * 86400}
        # since the last report (at most MAX_DAYS); Send now: the last EVERY_DAYS days, used or not
        since = now - EVERY_DAYS * 86400 if force else max(last, now - MAX_DAYS * 86400)
        days = max(1, -(-(now - since) // 86400))
        if not force and not _used_since(since):
            config.write_settings(share_last=now)   # nothing to send: the next look is in a few days
            return {"skipped": "not used"}
        sid = _active_session(since)
        stamp = time.strftime("%Y%m%d-%H%M")
        name = f"STAnnotator-usage-{stamp}-{report.install_id()}.zip"
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / name
            parts = ["conversations", "calls", "proposals", "annotation", "app_log"]
            report.build(path, about, int(days), parts, sid,
                         f"Automatic usage report ({int(days)} day(s); the user allowed usage sharing).")
            if path.stat().st_size > report.UPLOAD_MAX:   # too large: without the conversations
                report.build(path, about, int(days), [p for p in parts if p != "conversations"], sid,
                             "Automatic usage report (conversations left out: too large).")
            up = report.upload(path, name, about, extra={"auto": True, "install": report.install_id()})
        config.write_settings(share_last=now, share_sent=now)
        return {"sent": up["uploaded"], "bytes": up["bytes"], "days": int(days), "session": sid}
    except Exception as e:  # offline, the receiver away: try again at the next look
        return {"error": str(e)}
    finally:
        _lock.release()


def start(about_fn):
    """The desktop app's background loop: a first look a minute after start, then every few hours."""
    def loop():
        time.sleep(60)
        while True:
            try:
                r = send(about_fn())
                if "sent" in r or "error" in r:
                    print(f"usage sharing: {r}", flush=True)
            except Exception as e:
                print(f"usage sharing: {e}", flush=True)
            time.sleep(CHECK_SECONDS)
    threading.Thread(target=loop, daemon=True, name="usage-sharing").start()
