"""Recording proxy for the interactive tutorial (developer tool, not part of the app).

Sits between a browser and a running STAnnotator viewer and records every request the page makes, grouped by
"scene" (a named checkpoint of the walkthrough), so the tutorial player can replay the real app offline:

    STAGENT_HOME=work/tut/home .venv/bin/stagent viewer --port 8765
    .venv/bin/python tools/tutorial/record_proxy.py --upstream http://127.0.0.1:8765 --port 8766 --out work/tut/rec

Control (from the shell or the browser):  GET /__rec/scene?name=<scene>   start recording into a scene
                                          GET /__rec/off | /__rec/on      pause / resume recording
Output: <out>/rec.jsonl (one line per request) and <out>/blobs/<sha1> (bodies, stored once).
"""
import argparse, hashlib, json, time
from pathlib import Path

import httpx, uvicorn
from starlette.applications import Starlette
from starlette.requests import Request
from starlette.responses import JSONResponse, Response
from starlette.routing import Route

HOP = {"origin", "referer", "connection", "keep-alive", "transfer-encoding", "content-encoding", "content-length", "host"}


def make_app(upstream, out):
    out = Path(out)
    (out / "blobs").mkdir(parents=True, exist_ok=True)
    log = open(out / "rec.jsonl", "a")
    state = {"scene": "start", "on": True, "seq": 0}
    client = httpx.AsyncClient(base_url=upstream, timeout=None)

    def blob(data):
        h = hashlib.sha1(data).hexdigest()
        p = out / "blobs" / h
        if not p.exists():
            p.write_bytes(data)
        return h

    async def control(req: Request):
        cmd = req.path_params["cmd"]
        if cmd == "scene":
            state["scene"] = req.query_params["name"]
        elif cmd in ("on", "off"):
            state["on"] = cmd == "on"
        return JSONResponse(state)

    async def proxy(req: Request):
        body = await req.body()
        headers = {k: v for k, v in req.headers.items() if k.lower() not in HOP}
        r = await client.request(req.method, req.url.path, params=req.url.query, content=body, headers=headers)
        data = r.content
        ctype = r.headers.get("content-type", "")
        if state["on"] and (req.url.path.startswith("/api/") or req.url.path.startswith("/review")):
            state["seq"] += 1
            log.write(json.dumps({
                "seq": state["seq"], "t": round(time.time(), 3), "scene": state["scene"], "method": req.method,
                "path": req.url.path, "query": req.url.query, "req": blob(body) if body else None,
                "status": r.status_code, "ctype": ctype, "body": blob(data)}) + "\n")
            log.flush()
        hdr = {k: v for k, v in r.headers.items() if k.lower() not in HOP}
        return Response(data, status_code=r.status_code, headers=hdr)

    return Starlette(routes=[Route("/__rec/{cmd}", control),
                             Route("/{path:path}", proxy, methods=["GET", "POST", "PUT", "DELETE", "PATCH"])])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--upstream", default="http://127.0.0.1:8765")
    ap.add_argument("--port", type=int, default=8766)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    uvicorn.run(make_app(a.upstream, a.out), host="127.0.0.1", port=a.port, log_level="warning")
