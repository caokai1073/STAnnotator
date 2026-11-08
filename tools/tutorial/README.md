# Interactive tutorial (developers)

`docs/tutorial/` is the interactive tutorial on the website: the real app page (`stagent/static/app.html`, unchanged)
in the browser, answered from a **recorded session** instead of a server, with the agent's real conversations replayed
in a Claude Code pane, a virtual cursor, captions, notes and chapters. Pause at any step and the app is the viewer's.

Nothing here is part of the app or its package.

## How it works

| File | Role |
|---|---|
| `record_proxy.py` | sits between a browser and a running viewer; records every request and answer, grouped by **scene** |
| `crawl_tiles.py`, `crawl_genes.py` | add every H&E tile, and the map colored by a batch of genes, for the viewer's own exploring |
| `scenes.json` | the scenes in order; per scene: `from` / `until` (time window), `source` (a slice of another scene), `skip`, `inject` (answers the recording could not make); `rewrite` (local paths shown as generic ones), `synthetic`, `harvest` |
| `build.py` | builds `docs/tutorial/` (the app page with `player/shim.js` first, the player) and the data (the recorded answers, each body once, and the agent's conversations) |
| `player/shim.js` | in the app page: answers `fetch` and `<img src>` for `/api/...` from the recording; keeps the app's settings in memory |
| `player/engine.js` | the player: scenes, cursor, spotlight, terminal replay, timeline, pause-and-try, the self-check |
| `player/chapters.js` | **the script**: chapters, steps, captions, what the cursor does |

A step's `scene` switches what the app's requests get back, e.g. to "after the alignment was accepted". In a scene,
`GET /api/job/...` answers in the recorded order (progress); every other request gets what it was at the end of the
scene. The agent's own work is not replayed against the app: its conversation plays in the left pane, and at a tool
call (or at its end) the step switches the scene, so the app shows the agent's result as it did.

## Record a session

```bash
STAGENT_HOME=work/tut/home STAGENT_DESKTOP=1 .venv/bin/stagent viewer --port 8765
.venv/bin/python tools/tutorial/record_proxy.py --upstream http://127.0.0.1:8765 --port 8766 --out work/tut/rec
# open http://localhost:8766 in a browser; before each action:
curl "localhost:8766/__rec/scene?name=c08-mucosa"     # what the app shows after this action
curl "localhost:8766/__rec/off"                       # while the agent works (then a new scene for its result)
```

The agent runs as the desktop app would run it, in `$STAGENT_HOME/claude_workspace`, with the request the app sends
(e.g. `ALIGN_REQUEST` in `webapp.py`), saving its stream:

```bash
cd work/tut/home/claude_workspace && claude -p "<request>" --continue --mcp-config stagent_mcp.json \
  --allowedTools mcp__stagent WebSearch WebFetch --output-format stream-json --verbose > work/tut/transcripts/<name>.jsonl
```

Buttons that would open a Terminal (`/api/claude/...`) are not clicked while recording: their answer goes in
`scenes.json` (`synthetic`, or `inject` for one scene).

Demo data: two sections of SCP2927 (147328 terminal ileum, 167799 sigmoid colon), cropped to a window with mucosa,
submucosa and muscle (cells and H&E; the crop keeps the H&E pixel size in its OME metadata): `demo_crop.py`, and
`demo_geojson.py` for the QuPath-style outlines the Regions chapter imports.

## Build and check

```bash
.venv/bin/python tools/tutorial/crawl_tiles.py --out work/tut/rec
.venv/bin/python tools/tutorial/crawl_genes.py --rec work/tut/rec --genes EPCAM,MYH11,...
# the page in docs/tutorial (0.8 MB); the recorded session (55 MB) in the media repository, cloned in media/
.venv/bin/python tools/tutorial/build.py --rec work/tut/rec --agent work/tut/transcripts --out docs/tutorial \
  --data-out media/tutorial/data --data-url https://www.caokai.site/STAnnotator-media/tutorial/data/
# the website's layout, locally: http://localhost:8790/STAnnotator/tutorial/
mkdir -p work/tut/www && ln -sfn "$PWD/docs" work/tut/www/STAnnotator && ln -sfn "$PWD/media" work/tut/www/STAnnotator-media
cd work/tut/www && python3 -m http.server 8790
```

The data lives in [STAnnotator-media](https://github.com/caokai1073/STAnnotator-media) (`tutorial/data/`), so the
source repository stays small. The page reads it at its full address, so the website works wherever it is served
(www.caokai.site, or a fork's github.io: GitHub Pages lets any site read it). Without `--data-out` the data goes next to the page (`docs/tutorial/data`, `--data-url data/`).

`?check&speed=3` plays every chapter fast and lists the steps whose target was not found and the requests that were
not recorded (`?debug` logs every request). `check.py` does it in a headless browser and saves a screenshot of every
step (`uv pip install playwright && python -m playwright install chromium` once):

```bash
.venv/bin/python tools/tutorial/check.py http://localhost:8790/STAnnotator/tutorial/ --shots work/tut/shots
```
 Run it after changing the app's page: a renamed button or a moved control
shows up there, and the fix is in `chapters.js` (a target) or a new recording of that chapter.
