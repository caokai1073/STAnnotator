"""Self-check of the interactive tutorial in a headless browser (developer tool; needs playwright + chromium).

    .venv/bin/python tools/tutorial/check.py http://localhost:8790/tutorial/ [--speed 3] [--shots work/tut/shots]

Plays every chapter fast (?check), prints the steps whose target was not found and the requests that were not
recorded, and saves a screenshot of every step for a look.
"""
import argparse, json, time
from pathlib import Path

from playwright.sync_api import sync_playwright

ap = argparse.ArgumentParser()
ap.add_argument("url")
ap.add_argument("--speed", default="3")
ap.add_argument("--chapter", default="")
ap.add_argument("--shots", default="")
ap.add_argument("--timeout", type=int, default=3600)
a = ap.parse_args()
shots = Path(a.shots) if a.shots else None
if shots:
    shots.mkdir(parents=True, exist_ok=True)
with sync_playwright() as p:
    b = p.chromium.launch()
    pg = b.new_page(viewport={"width": 1440, "height": 900})
    errors = []
    pg.on("pageerror", lambda e: errors.append(str(e)[:200]))
    q = f"?check&speed={a.speed}" + (f"&chapter={a.chapter}" if a.chapter else "")
    pg.goto(a.url + q)
    t0, last = time.time(), None
    while time.time() - t0 < a.timeout:
        time.sleep(1.5)
        s = pg.evaluate("() => { const S = window.__player && window.__player.S; return S ? {ci: S.ci, si: S.si, done: !!S.done, ch: window.TUTORIAL.chapters[S.ci].id} : null }")
        if not s:
            continue
        if (s["ci"], s["si"]) != last:
            last = (s["ci"], s["si"])
            print(f"{time.time() - t0:6.0f}s  {s['ch']} step {s['si'] + 1}", flush=True)
            if shots:
                time.sleep(2.5)
                pg.screenshot(path=str(shots / f"{s['ci'] + 1:02d}_{s['ch']}_{s['si'] + 1:02d}.jpg"), quality=70, type="jpeg")
        if s["done"]:
            break
    res = pg.evaluate("() => ({ fails: __player.S.fails || [], misses: __player.S.misses || [] })")
    print(json.dumps(res, indent=1))
    print("page errors:", json.dumps(errors[:20], indent=1))
    b.close()
