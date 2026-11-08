"""Build the interactive tutorial (docs/tutorial/) from a recorded session (developer tool).

    .venv/bin/python tools/tutorial/build.py --rec work/tut/rec --agent work/tut/transcripts --out docs/tutorial

- the app page: stagent/static/app.html as it is, with the replay shim (player/shim.js) loaded first
- the data: every recorded answer, grouped by scene (tools/tutorial/scenes.json), each body stored once
- the agent: the Claude Code conversations (stream-json) as short event lists for the terminal pane
- the player: tools/tutorial/player/*
"""
import argparse, hashlib, json, re, shutil
from collections import OrderedDict
from pathlib import Path
from urllib.parse import parse_qsl

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
EXT = {"application/json": "json", "image/jpeg": "jpg", "image/png": "png", "application/octet-stream": "bin",
       "text/plain": "txt", "text/markdown": "md", "text/csv": "csv", "application/pdf": "pdf", "image/svg+xml": "svg",
       "application/geo+json": "json", "text/html": "html"}
VOLATILE = {"t", "_"}  # cache busters: not part of the key


def norm(path, query):
    """The key of a request: path and decoded query pairs in order (the player builds the same key)."""
    q = sorted((k, v) for k, v in parse_qsl(query, keep_blank_values=True) if k not in VOLATILE)
    return path + ("?" + "&".join(f"{k}={v}" for k, v in q) if q else "")


def load_rewrites(cfg):
    return [(a, b) for a, b in cfg.get("rewrite", [])]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rec", required=True)
    ap.add_argument("--agent", required=True)
    ap.add_argument("--out", default=str(ROOT / "docs" / "tutorial"))
    ap.add_argument("--data-out", help="where the recorded data goes (default: <out>/data); on the website it lives in "
                                       "the media repository, e.g. media/tutorial/data")
    ap.add_argument("--data-url", default="data/", help="where the page finds it, e.g. /STAnnotator-media/tutorial/data/")
    a = ap.parse_args()
    rec, out = Path(a.rec), Path(a.out)
    dout = Path(a.data_out) if a.data_out else out / "data"
    cfg = json.loads((HERE / "scenes.json").read_text())
    rewrites = load_rewrites(cfg)
    for d in (out, dout):
        if d.exists():
            shutil.rmtree(d)
    out.mkdir(parents=True)
    (dout / "b").mkdir(parents=True)

    # ---------- recorded answers ----------
    rows = [json.loads(l) for l in open(rec / "rec.jsonl")]
    for extra in ("rec_tiles.jsonl", "rec_genes.jsonl"):  # crawl_tiles.py, crawl_genes.py: for the viewer's own looking
        if (rec / extra).exists():
            rows += [json.loads(l) for l in open(rec / extra)]
    rules = {s["name"]: s for s in cfg["scenes"]}
    rules["_tiles"] = {"name": "_tiles"}
    rules["_genes"] = {"name": "_genes"}
    order = [s["name"] for s in cfg["scenes"]]
    stored, entries = {}, OrderedDict((n, OrderedDict()) for n in order + ["_tiles", "_genes"])

    def store(sha, ctype):
        if (sha, ctype) in stored:
            return stored[(sha, ctype)]
        data = (rec / "blobs" / sha).read_bytes()
        ct = ctype.split(";")[0].strip()
        if ct in ("application/json", "text/plain", "text/markdown", "text/csv", "application/geo+json", "text/html"):
            txt = data.decode("utf-8", "replace")
            for x, y in rewrites:
                txt = txt.replace(x, y)
            data = txt.encode()
        h = hashlib.sha1(data).hexdigest()[:16]
        name = f"{h}.{EXT.get(ct, 'bin')}"
        p = dout / "b" / name
        if not p.exists():
            p.write_bytes(data)
        stored[(sha, ctype)] = name
        return name

    entries["_static"] = OrderedDict()  # lists that do not change (organs, options): from any scene, dropped ones too
    for r in rows:
        key = r["method"] + " " + norm(r["path"], r["query"])
        if r["status"] == 200 and any(key.startswith(k) for k in cfg.get("harvest", [])):
            entries["_static"][key] = [[r["status"], r["ctype"].split(";")[0], store(r["body"], r["ctype"]), ""]]
    by_source = {}  # a scene can be part of a recorded one (source + from/until), e.g. two imports in one recording
    for rule in cfg["scenes"]:
        by_source.setdefault(rule.get("source", rule["name"]), []).append(rule)
    by_source["_tiles"], by_source["_genes"] = [rules["_tiles"]], [rules["_genes"]]
    for r in rows:
        for rule in by_source.get(r["scene"], []):
            if rule.get("drop"):
                continue
            if rule.get("until") and r["t"] > rule["until"]:
                continue
            if rule.get("from") and r["t"] < rule["from"]:
                continue
            key = r["method"] + " " + (r["path"] if rule["name"] == "_tiles" else norm(r["path"], r["query"]))
            if any(key.startswith(k) for k in rule.get("skip", [])):
                continue
            f = store(r["body"], r["ctype"])
            ans = [r["status"], r["ctype"].split(";")[0], f]
            if r["method"] != "GET":  # the request body picks among the answers of the same path
                body = (rec / "blobs" / r["req"]).read_text("utf-8", "replace") if r.get("req") else ""
                for x, y in rewrites:
                    body = body.replace(x, y)
                ans.append(body)
            lst = entries[rule["name"]].setdefault(key, [])
            if not lst or lst[-1] != ans:
                lst.append(ans)
    # answers the recording could not make (they would open a real Terminal): per scene (inject), or in every scene
    for rule in cfg["scenes"]:
        for key, val in rule.get("inject", {}).items():
            data = json.dumps(val).encode()
            name = hashlib.sha1(data).hexdigest()[:16] + ".json"
            (dout / "b" / name).write_bytes(data)
            entries[rule["name"]][key] = [[200, "application/json", name, ""]]
    entries["_synthetic"] = OrderedDict()
    for key, val in cfg.get("synthetic", {}).items():
        data = json.dumps(val).encode()
        name = hashlib.sha1(data).hexdigest()[:16] + ".json"
        (dout / "b" / name).write_bytes(data)
        entries["_synthetic"][key] = [[200, "application/json", name, ""]]
    manifest = {"scenes": order, "global": ["_synthetic", "_tiles", "_genes", "_static"], "entries": entries}
    (dout / "manifest.json").write_text(json.dumps(manifest, separators=(",", ":")))

    # ---------- the agent's conversations ----------
    agent = {}
    for p in sorted(Path(a.agent).glob("*.jsonl")):
        if p.name.startswith("_"):
            continue
        prompt = (p.with_suffix(".prompt.txt")).read_text() if p.with_suffix(".prompt.txt").exists() else ""
        agent[p.stem] = conversation(p, prompt, rewrites)
    (dout / "agent.json").write_text(json.dumps(agent, separators=(",", ":")))

    # ---------- the app page, unchanged but for the shim and the vendor path ----------
    (out / "app").mkdir()
    html = (ROOT / "stagent" / "static" / "app.html").read_text()
    html = html.replace("<head>", '<head>\n<script src="../shim.js"></script>', 1)
    html = html.replace('src="/static/vendor/', 'src="vendor/')
    (out / "app" / "index.html").write_text(html)
    shutil.copytree(ROOT / "stagent" / "static" / "vendor", out / "app" / "vendor")
    for f in (HERE / "player").iterdir():
        if f.is_file():
            shutil.copy(f, out / f.name)
    (out / "config.js").write_text(f"window.TUTORIAL_DATA = {json.dumps(a.data_url)};  // where the recorded session is (build.py --data-url)\n")
    n = sum(1 for _ in (dout / "b").iterdir())
    size = lambda d: sum(f.stat().st_size for f in d.rglob("*") if f.is_file()) / 1e6
    print(f"{len(rows)} requests, {sum(len(e) for e in entries.values())} keys, {n} files: page {size(out):.1f} MB -> {out}, "
          f"data {size(dout):.1f} MB -> {dout} (served at {a.data_url})")


def summarize(content):
    if isinstance(content, list):
        parts = []
        for c in content:
            if c.get("type") == "text":
                parts.append(c["text"])
            elif c.get("type") == "image":
                parts.append("[Image]")
        content = "\n".join(parts)
    return str(content)


def conversation(path, prompt, rewrites):
    """stream-json -> [{k:'user'|'text'|'tool'|'result', ...}] for the terminal pane."""
    ev, by_id = [{"k": "user", "text": prompt}], {}
    for line in open(path):
        e = json.loads(line)
        if e["type"] == "assistant":
            for c in e["message"]["content"]:
                if c["type"] == "text" and c["text"].strip():
                    ev.append({"k": "text", "text": c["text"].strip()})
                elif c["type"] == "tool_use" and c["name"] != "ToolSearch":
                    t = {"k": "tool", "name": c["name"], "input": c["input"], "out": None}
                    by_id[c["id"]] = t
                    ev.append(t)
        elif e["type"] == "user":
            for c in e["message"]["content"]:
                if isinstance(c, dict) and c.get("type") == "tool_result" and c.get("tool_use_id") in by_id:
                    by_id[c["tool_use_id"]]["out"] = summarize(c.get("content"))[:4000]
                    by_id[c["tool_use_id"]]["error"] = bool(c.get("is_error"))
        elif e["type"] == "result":
            ev.append({"k": "done", "ms": e.get("duration_ms"), "turns": e.get("num_turns")})
    txt = json.dumps(ev)
    for x, y in rewrites:
        txt = txt.replace(json.dumps(x)[1:-1], json.dumps(y)[1:-1])
    return json.loads(txt)


if __name__ == "__main__":
    main()
