"""The viewer for a recording: the app as is, except that a request to the agent (Align, Markers' Refine with the agent,
Refine's Agent review and its next rounds, Review, a prior in your words) is written to requests.jsonl in the data
folder instead of opening Claude Code in Terminal. The app answers as if Claude Code had opened, so these buttons can be
clicked while recording; the agent then runs from the saved prompt (README: Record a session).

    STAGENT_HOME=work/tut/home STAGENT_DESKTOP=1 .venv/bin/python tools/tutorial/record_viewer.py 8775
"""
import json, os, sys, time
from stagent import claudelink, config


def asked(sample=None, session=None, request=None, send_if_open=False, remote=None):
    claudelink.workspace()
    prompt = None  # what open_claude_code would give Claude Code
    if session:
        prompt = (f"I am looking at sample {sample} (session {session}) in the STAnnotator app. "
                  "Call use_viewer_session to load it, then " + (f"do this: {request}" if request else "wait for my request."))
    elif request:
        prompt = request
    with open(config.HOME / "requests.jsonl", "a", encoding="utf-8") as f:
        f.write(json.dumps({"t": time.strftime("%H:%M:%S"), "sample": sample, "session": session, "prompt": prompt}) + "\n")
    return {"opened": True, "reused": True}


claudelink.open_claude_code = asked
claudelink.workspace()
sys.argv = ["stagent-web", "--port", sys.argv[1] if len(sys.argv) > 1 else "8775"]
from stagent.webapp import main  # noqa: E402

main()
