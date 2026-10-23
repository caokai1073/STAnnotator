"""Smoke test of the MCP server over stdio, as a client (e.g. Claude Code) would call it."""
import asyncio, json, os, shutil, sys
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

PY = sys.executable


def show(res, n=400):
    for c in res.content:
        print("  ", c.type, (c.text[:n] if c.type == "text" else f"{c.mimeType} {len(c.data)} b64 chars"))
    if res.isError:
        print("   (isError)")


async def main():
    params = StdioServerParameters(command=PY, args=["-m", "celltypeagent.mcp_server"], env=dict(os.environ))
    async with stdio_client(params) as (r, w):
        async with ClientSession(r, w) as s:
            await s.initialize()
            tools = await s.list_tools()
            print("tools:", [t.name for t in tools.tools])
            calls = [
                ("start_session", {"sample_id": "167799", "region_source": "pathologist"}),
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
                if name == "start_session":
                    sid = json.loads(res.content[0].text)["session_id"]
                if name == "conflict_summary" and not res.isError:
                    groups = json.loads(res.content[0].text)["groups"]
                    gid = groups[0]["group_id"] if groups else None
            if gid:
                print(f"\n>> show_conflicts {gid}")
                show(await s.call_tool("show_conflicts", {"group_id": gid, "n": 3}))
            print("\n>> history")
            show(await s.call_tool("history", {}), 600)
    from celltypeagent import config
    shutil.rmtree(config.SESSIONS / sid)  # keep the sessions folder free of test sessions


asyncio.run(main())
