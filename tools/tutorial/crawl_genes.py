"""Color every recorded annotation version by the usual modes (cell type, lineage, niche, region, Leiden) and a batch of genes, into the recording (scene "_genes"), so viewers of
the tutorial can pause and color the map by any of these genes (developer tool).

    .venv/bin/python tools/tutorial/crawl_genes.py --server http://127.0.0.1:8765 --rec work/tut/rec --genes EPCAM,MYH11,...
"""
import argparse, hashlib, json, re, time
from pathlib import Path
from urllib.parse import parse_qsl

import httpx

ap = argparse.ArgumentParser()
ap.add_argument("--server", default="http://127.0.0.1:8765")
ap.add_argument("--rec", required=True)
ap.add_argument("--genes", required=True)
a = ap.parse_args()
rec = Path(a.rec)
c = httpx.Client(base_url=a.server, timeout=120)
pairs = set()  # (session, version) the recording showed
for line in open(rec / "rec.jsonl"):
    r = json.loads(line)
    m = re.match(r"/api/s/([^/]+)/colors\.bin$", r["path"])
    if m:
        q = dict(parse_qsl(r["query"], keep_blank_values=True))
        pairs.add((m.group(1), q.get("version")))
rows = []
alive = {s["session_id"] for s in c.get("/api/sessions").json()} if c.get("/api/sessions").status_code == 200 else None
for sid, v in sorted(pairs):
    if alive is not None and sid not in alive:
        continue
    sample = sid.split("_")[0]
    panel = set(c.get(f"/api/sample/{sample}/genes").json().get("genes", [])) if True else set()
    modes = ["type", "lineage", "niche", "region", "leiden"] + [f"gene:{g}" for g in a.genes.split(",") if not panel or g in panel]
    for mode in modes:
        for path, q in ((f"/api/s/{sid}/colors.bin", {"version": v, "mode": mode, "only": ""}),
                        (f"/api/s/{sid}/legend", {"version": v, "mode": mode})):
            try:
                r = c.get(path, params=q)
            except httpx.HTTPError:
                continue
            if r.status_code != 200:
                continue
            sha = hashlib.sha1(r.content).hexdigest()
            if not (rec / "blobs" / sha).exists():
                (rec / "blobs" / sha).write_bytes(r.content)
            rows.append({"seq": 0, "t": time.time(), "scene": "_genes", "method": "GET", "path": path,
                         "query": "&".join(f"{k}={val}" for k, val in q.items()), "req": None, "status": 200,
                         "ctype": r.headers.get("content-type", ""), "body": sha})
    print(sid, v, len(rows))
with open(rec / "rec_genes.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
