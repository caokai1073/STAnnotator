"""Fetch every H&E tile of every sample from a running viewer into the recording (scene "_tiles"), so the tutorial
can zoom anywhere on the image, not only where the recorded session looked (developer tool).

    .venv/bin/python tools/tutorial/crawl_tiles.py --server http://127.0.0.1:8765 --out work/tut/rec
"""
import argparse, hashlib, json, math, time
from pathlib import Path

import httpx

ap = argparse.ArgumentParser()
ap.add_argument("--server", default="http://127.0.0.1:8765")
ap.add_argument("--out", required=True)
a = ap.parse_args()
out = Path(a.out)
(out / "blobs").mkdir(parents=True, exist_ok=True)
c = httpx.Client(base_url=a.server, timeout=60)
rows = []
for s in c.get("/api/samples").json():
    sid = s["sample_id"]
    info = c.get(f"/api/sample/{sid}/he_info").json()
    w, h, ts, rev = info["width"], info["height"], info["tile_size"], info["rev"]
    top = math.ceil(math.log2(max(w, h)))
    n = 0
    for L in range(0, top + 1):
        k = 2 ** (L - top)
        for x in range(math.ceil(w * k / ts)):
            for y in range(math.ceil(h * k / ts)):
                path, q = f"/api/sample/{sid}/he/{L}/{x}_{y}.jpg", f"r={rev}"
                r = c.get(path, params={"r": rev})
                if r.status_code != 200:
                    continue
                sha = hashlib.sha1(r.content).hexdigest()
                p = out / "blobs" / sha
                if not p.exists():
                    p.write_bytes(r.content)
                rows.append({"seq": 0, "t": time.time(), "scene": "_tiles", "method": "GET", "path": path, "query": q,
                             "req": None, "status": 200, "ctype": r.headers.get("content-type", "image/jpeg"), "body": sha})
                n += 1
    print(sid, w, h, "levels 0-%d" % top, n, "tiles")
with open(out / "rec_tiles.jsonl", "w") as f:
    for r in rows:
        f.write(json.dumps(r) + "\n")
