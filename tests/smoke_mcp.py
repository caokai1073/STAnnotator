"""Smoke test of the MCP server over stdio, as a client (e.g. Claude Code) would call it. Tools that change data only
propose the change; accept() plays the user clicking Accept in the app (the app's own apply_proposal), and the agent
then reloads the session.

Run: python tests/smoke_mcp.py: the repo's mouse colon example in a temporary data folder. Or on a data folder that
holds an SCP2927 IBD Xenium sample: STAGENT_HOME=<folder> SMOKE_SAMPLE=<its id> python tests/smoke_mcp.py"""
import asyncio, json, os, shutil, sys
import tempfile
from pathlib import Path
os.environ.setdefault("STAGENT_SETTINGS_DIR", tempfile.mkdtemp(prefix="stagent_cfg_"))  # (never the user's settings)
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

PY = sys.executable
EX = Path(__file__).parent.parent / "examples"
TMP_HOME = None
if "STAGENT_HOME" not in os.environ:  # on its own: the example, imported into a temporary data folder
    TMP_HOME = tempfile.mkdtemp(prefix="stagent_mcp_")
    os.environ["STAGENT_HOME"] = TMP_HOME
    from stagent import importer
    importer.import_sample("smoke", str(EX / "cadinu_merfish" / "cadinu_merfish_D0_m3a_1_slice_2.h5ad"), progress=lambda *a: None)
    SAMPLE, MARKERS = "smoke", open(EX / "cadinu_merfish" / "mouse_colon_merfish_markers.txt").read()
else:
    SAMPLE = os.environ.get("SMOKE_SAMPLE", "167799")
    MARKERS = open(EX / "scp2927" / "ibd_xenium480_markers.csv").read()


def show(res, n=400):
    for c in res.content:
        print("  ", c.type, (c.text[:n] if c.type == "text" else f"{c.mimeType} {len(c.data)} b64 chars"))
    if res.isError:
        print("   (isError)")


async def accept(s, res):
    """The user clicks Accept in the app's review window: the proposal's tool runs, then the agent reloads."""
    from stagent import mcp_server, proposals
    rec = json.loads(res.content[0].text)
    assert rec.get("status") == "waiting for the user's review", rec
    out, sid = mcp_server.apply_proposal(proposals.get(rec["proposal"]))
    proposals.update(rec["proposal"], status="done", result=out, session_after=sid)
    r = await s.call_tool("resume_session", {"session_id": sid})
    assert not r.isError, r
    print("   accepted", rec["proposal"], "->", sid)
    return sid


async def main():
    params = StdioServerParameters(command=PY, args=["-m", "stagent.mcp_server"], env=dict(os.environ))
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            tools = await s.list_tools()
            print("tools:", [t.name for t in tools.tools])
            calls = [
                ("start_session", {"sample_id": SAMPLE, "markers_text": MARKERS}),
                ("list_regions", {}),
                ("preview_constraint", {"text": "Epithelial cells are only found in the mucosa", "types": ["epithelial"],
                                        "regions": ["mucosa"], "mode": "only"}),
                ("add_constraint", {"text": "Epithelial cells are only found in the mucosa", "types": ["epithelial"],
                                    "regions": ["mucosa"], "mode": "only"}),
                ("conflict_summary", {"by": ["constraint", "region"], "top": 5}),
                ("preview_constraint", {"text": "CD4 T cells next to B cells", "types": ["CD4 T cells"], "mode": "only",
                                        "weight": "weak", "kind": "neighborhood", "neighbor_types": ["B cells"], "min_frac": 0.2}),
                ("list_rule_suggestions", {}),
                ("preview_constraint", {"text": "bad", "types": ["Gobelt cells"], "regions": ["mucosa"], "mode": "only"}),
            ]
            gid = None
            for name, args in calls:
                print(f"\n>> {name}")
                res = await s.call_tool(name, args)
                show(res)
                assert not res.isError or "Gobelt" in json.dumps(args), f"{name} failed"
                if name == "start_session":
                    sid = await accept(s, res)
                    from stagent.session import Session  # the user draws a region in the viewer
                    import numpy as np
                    c = Session(sid).cells
                    cx, cy, h = float(np.median(c.px)), float(np.median(c.py)), 2500.0
                    Session(sid).add_region("mucosa", [[cx - h, cy - h], [cx + h, cy - h], [cx + h, cy + h], [cx - h, cy + h]])
                    await s.call_tool("resume_session", {"session_id": sid})  # the agent sees the new region
                if name == "add_constraint":
                    await accept(s, res)
                if name == "conflict_summary" and not res.isError:
                    groups = json.loads(res.content[0].text)["groups"]
                    gid = groups[0]["group_id"] if groups else None
            if gid:
                print(f"\n>> show_conflicts {gid}")
                show(await s.call_tool("show_conflicts", {"group_id": gid, "n": 3}))
            print("\n>> history")
            show(await s.call_tool("history", {}), 600)
            print("\n>> pending_changes")
            show(await s.call_tool("pending_changes", {}), 600)
            print("ok")
    from stagent import config
    shutil.rmtree(config.SESSIONS / sid)  # keep the sessions folder free of test sessions


try:
    asyncio.run(main())
finally:
    if TMP_HOME:
        shutil.rmtree(TMP_HOME, ignore_errors=True)
