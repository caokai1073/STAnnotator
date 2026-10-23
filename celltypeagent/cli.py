"""Command line: celltypeagent {import, samples, setup-desktop, viewer, mcp}."""
import argparse, json, os, platform, shutil, sys, time
from pathlib import Path
from . import config


def desktop_config_path():
    system = platform.system()
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "Claude" / "claude_desktop_config.json"
    if system == "Windows":
        return Path(os.environ.get("APPDATA", Path.home())) / "Claude" / "claude_desktop_config.json"
    return Path.home() / ".config" / "Claude" / "claude_desktop_config.json"


def setup_desktop(config_path=None, name="celltypeagent"):
    """Register this installation's MCP server in Claude Desktop (existing entries are kept; a backup is written)."""
    path = Path(config_path) if config_path else desktop_config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    cfg = json.loads(path.read_text()) if path.exists() and path.read_text().strip() else {}
    if path.exists():
        shutil.copy(path, path.with_suffix(f".backup-{time.strftime('%Y%m%d-%H%M%S')}.json"))
    cfg.setdefault("mcpServers", {})[name] = {
        "command": sys.executable, "args": ["-m", "celltypeagent.mcp_server"],
        "env": {"CELLTYPEAGENT_HOME": str(config.HOME)}}
    path.write_text(json.dumps(cfg, indent=2))
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(prog="celltypeagent", description="Spatial cell type annotation with Claude.")
    sub = ap.add_subparsers(dest="cmd", required=True)
    im = sub.add_parser("import", help="import one sample from h5ad + registered H&E + registration file")
    im.add_argument("--name", required=True, help="sample id used in the viewer and in Claude")
    im.add_argument("--h5ad", required=True, nargs="+", help="one or more h5ad files (merged)")
    im.add_argument("--he", required=True, help="registered H&E image (pyramidal TIFF / OME-TIFF)")
    im.add_argument("--registration", required=True, help="affine CSV (a11,a12,a21,a22,tx,ty) or Xenium Explorer 3x3 matrix")
    im.add_argument("--section-column", help="obs column that identifies the section (if the h5ad holds several)")
    im.add_argument("--section", help="value of --section-column to import")
    im.add_argument("--x-column"); im.add_argument("--y-column")
    im.add_argument("--counts-layer", help="layer with raw counts (default: 'counts' if present, else X)")
    im.add_argument("--region-column", help="optional obs column with region annotations to start from")
    im.add_argument("--xenium-px-um", type=float, default=0.2125, help="morphology pixel size for Xenium Explorer matrices")
    im.add_argument("--description", default="")
    im.add_argument("--overwrite", action="store_true")
    sub.add_parser("samples", help="list imported samples")
    sd = sub.add_parser("setup-desktop", help="register the MCP server in Claude Desktop")
    sd.add_argument("--config", help="path of claude_desktop_config.json (default: the platform's location)")
    vw = sub.add_parser("viewer", help="start the web viewer")
    vw.add_argument("--port", type=int, default=8765)
    sub.add_parser("mcp", help="run the MCP server on stdio (Claude Desktop starts this itself)")
    a = ap.parse_args(argv)
    config.ensure_dirs()
    if a.cmd == "import":
        from .importer import import_sample
        import_sample(a.name, a.h5ad, a.he, a.registration, a.section_column, a.section, a.x_column, a.y_column,
                      a.counts_layer, a.region_column, a.xenium_px_um, a.description, a.overwrite)
    elif a.cmd == "samples":
        from .session import list_samples
        t = list_samples()
        print(t.to_string(index=False) if len(t) else f"no samples in {config.SAMPLES}")
    elif a.cmd == "setup-desktop":
        p = setup_desktop(a.config)
        print(f"registered 'celltypeagent' in {p}\ndata folder: {config.HOME}\nRestart Claude Desktop to load it.")
    elif a.cmd == "viewer":
        sys.argv = ["celltypeagent-web", "--port", str(a.port)]
        from .webapp import main as web
        web()
    elif a.cmd == "mcp":
        from .mcp_server import main as mcp
        mcp()


if __name__ == "__main__":
    main()
